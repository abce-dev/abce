# Agent decision model

##########################################################################
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################

println("\n-----------------------------------------------------------")
println("Julia agent choice algorithm: starting")
println("Loading packages...")
using JuMP, GLPK, LinearAlgebra, DataFrames, CSV, Printf, YAML, SQLite

# Load settings and file locations from the settings file
settings_file = ARGS[1]
settings = YAML.load_file(settings_file)

# Include local ABCE functions module
julia_ABCE_module = joinpath(settings["ABCE_abs_path"], "ABCEfunctions.jl")
include(julia_ABCE_module)
using .ABCEfunctions
println("Packages loaded successfully.")

###### Set up inputs
println("Initializing data...")

# File names
db_file = joinpath(settings["ABCE_abs_path"], settings["db_file"])
# Constants
hours_per_year = settings["hours_per_year"]
consider_future_projects = settings["consider_future_projects"]
if consider_future_projects
    num_lags = settings["num_future_periods_considered"]
else
    num_lags = 0
end
MW2kW = 1000   # Converts MW to kW
MMBTU2BTU = 1000   # Converts MMBTU to BTU

# Load the inputs
db = load_db(db_file)
pd = get_current_period()
agent_id = get_agent_id()

# Set up agent-specific data
# Get a list of all ongoing construction projects for the current agent
agent_projects = get_WIP_projects_list(db, pd, agent_id)

# Get agent financial parameters
agent_params = get_agent_params(db, agent_id)
d = agent_params[1, :discount_rate]

# System parameters
# Read unit operational data (unit_data) and number of unit types (num_types)
unit_data, num_types = get_unit_specs(db)
num_alternatives = num_types * (num_lags + 1)

# Ensure that forecast horizon is long enough to accommodate the end of life
#   for the most long-lived possible unit
fc_pd = set_forecast_period(unit_data, num_lags)

# Load the demand data
available_demand = get_demand_forecast(db, pd, agent_id, fc_pd, settings)

# Extend the unserved demand data to match the total forecast period (constant projection)
available_demand = get_net_demand(db, pd, agent_id, fc_pd, available_demand)
println(DataFrame(net_demand = available_demand)[1:10, :])

# Load the price data
price_curve = DBInterface.execute(db, "SELECT * FROM price_curve") |> DataFrame

# Add empty column for project NPVs in unit_data
unit_data[!, :FCF_NPV] = zeros(Float64, num_types)

# NPV is the only decision criterion, so create a dataframe to hold the results
#    for each alternative
alternative_names, NPV_results = create_NPV_results_df(unit_data, num_lags)

# Create per-unit financial statement tables
println("Creating and populating unit financial statements for NPV calculation")
unit_FS_dict = create_unit_FS_dict(unit_data, fc_pd, num_lags)

# Populate financial statements with top-line data
# This data is deterministic, and is on a per-unit basis (not per-kW or per-kWh)
# Therefore, you can multiply each of these dataframes by its corresponding
#    u[] value to determine the actual impact of the units chosen on the GC's
#    final financial statement
for i = 1:num_types
    for j = 0:num_lags
        # Set up parameters for this alternative
        unit_type = unit_data[i, :unit_type]
        name = string(unit_type, "_lag-", j)
        fs = unit_FS_dict[name]
        unit_type_data = filter(row -> row.unit_type == unit_type, unit_data)

        # Generate the alternative's construction expenditure profile and save
        #   it to the FS
        fs[!, :xtr_exp] = generate_xtr_exp_profile(unit_type_data, j, fc_pd)

        # Set up the time-series of outstanding debt principal based on this
        #   expenditure profile: sets unit_fs[!, :remaining_debt_principal]
        #   for all construction periods
        set_initial_debt_principal_series(fs, unit_type_data, j, agent_params)

        # Generate prime-mover events from the end of construction to the end of the unit's life
        generate_prime_movers(unit_type_data, fs, j, agent_params[1, :cost_of_debt])

        # Forecast unit revenue ($/period) and generation (kWh/period)
        forecast_unit_revenue_and_gen(unit_type_data, fs, price_curve, db, pd, j)


        # Apply reactive functions to the rest of the dataframe
        # Compute total cost of fuel used
        transform!(fs, [:gen] => ((gen) -> unit_data[i, :FC_per_MMBTU] .* unit_data[i, :heat_rate] ./ (MW2kW * MMBTU2BTU) .* gen) => :Fuel_Cost)
        # Compute total VOM cost incurred during generation
        transform!(fs, [:gen] => ((gen) -> unit_data[i, :VOM] .* gen) => :VOM_Cost)
        # Compute total FOM cost for the year (non-reactive)
        fs[!, :FOM_Cost] = zeros(size(fs)[1])
        fs[(j + unit_data[i, :d_x]+1):(j + unit_data[i, :d_x]+unit_data[i, :unit_life]), :FOM_Cost] .= unit_data[i, :FOM] * unit_data[i, :capacity] * MW2kW
        # Compute EBITDA
        transform!(fs, [:Revenue, :Fuel_Cost, :VOM_Cost, :FOM_Cost] => ((rev, fc, VOM, FOM) -> rev - fc - VOM - FOM) => :EBITDA)
        # Compute EBIT
        transform!(fs, [:EBITDA, :depreciation] => ((EBITDA, dep) -> EBITDA - dep) => :EBIT)
        # Compute EBT
        transform!(fs, [:EBIT, :interest_payment] => ((EBIT, interest) -> EBIT - interest) => :EBT)
        # Compute taxes owed
        transform!(fs, [:EBT] => ((EBT) -> EBT .* agent_params[1, :tax_rate]) => :tax_owed)
        # Compute net income
        transform!(fs, [:EBT, :tax_owed] => ((EBT, tax) -> EBT - tax) => :Net_Income)
        # Compute FCF
        transform!(fs, [:Net_Income, :interest_payment, :xtr_exp] => ((NI, interest, xtr_exp) -> NI + interest - xtr_exp) => :FCF)

        # Add column of compounded discount factors
        transform!(fs, [:year] => ((year) -> (1+d) .^ (-1 .* (year .- 1))) => :d_factor)
        # Discount unit FCF NPV values and save to the NPV_results dataframe
        NPV_results[findall(NPV_results.name .== name)[1], :NPV] = transpose(fs[!, :FCF]) * fs[!, :d_factor]
        unit_data[i, :FCF_NPV] = transpose(fs[!, :FCF]) * fs[!, :d_factor]
    end
end

println(NPV_results)

if pd == 0
    println(unit_data)
end

println("Data initialized.")

###### Set up the model
# Create the model
println("Setting up model...")
m = Model(GLPK.Optimizer)
# Turn on higher model verbosity for debugging
#set_optimizer_attribute(m, "msg_lev", GLPK.GLP_MSG_ALL)
@variable(m, u[1:num_alternatives] >= 0, Int)
@variable(m, z[1:fc_pd])

# Restrict total construction to be less than maximum available demand (subject to capacity factor)
# To prevent unwanted infeasibility, convert nonpositive available_demand values to 0
for i = 1:size(available_demand)[1]
    if available_demand[i] < 0
        available_demand[i] = 0
    end
end

for i = 1:size(unit_FS_dict[alternative_names[1]])[1]
    @constraint(m, sum(u[j] * unit_FS_dict[alternative_names[j]][i, :gen] for j=1:num_alternatives) / (hours_per_year*MW2kW) <= available_demand[i] * 2)
end
# Constraint on max amount of interest payable per year
#@constraint(m, transpose(u) * (unit_data[!, :uc_x] .* unit_data[!, :capacity] .* agent_params[1, :debt_fraction] .* d ./ (1 .- (1+d) .^ (-1 .* unit_data[!, :unit_life]))) <= agent_params[1, :interest_cap])

#@objective(m, Max, transpose(u) * unit_data[!, :FCF_NPV] - 0.25 * sum((available_demand[i] - sum(u))^2 for i=1:size(available_demand)[0]))
@objective(m, Max, transpose(u) * NPV_results[!, :NPV])
println("Model set up.")


###### Solve the model
println("Solving problem...")
optimize!(m)
status = termination_status.(m)
unit_qty = value.(u)


###### Display the results
println(status)
println("Units to build:")
println(hcat(alternative_names, DataFrame(units = unit_qty)))


###### Save the new units into the `assets` and `WIP_projects` DB tables
for i = 1:num_alternatives
    unit_type = join(deleteat!(split(alternative_names[i], "_"), size(split(alternative_names[i], "_"))[1]), "_")
    unit_index = findall(unit_data.unit_type .== unit_type)[1]
    if occursin("0", alternative_names[i])
        # Only record projects starting this period
        for j = 1:unit_qty[i]
            next_id = get_next_asset_id(db)
            # Update `WIP_projects` table
            rcec = unit_data[unit_index, :uc_x] * unit_data[unit_index, :capacity] * MW2kW
            rtec = unit_data[unit_index, :d_x]
            WIP_projects_vals = (next_id, agent_id, pd, rcec, rtec, rcec / 10)
            DBInterface.execute(db, "INSERT INTO WIP_projects VALUES (?, ?, ?, ?, ?, ?)", WIP_projects_vals)

            # Update `assets` table
            revealed = "false"
            completion_pd = pd + unit_data[unit_index, :d_x]
            cancellation_pd = 9999
            retirement_pd = pd + unit_data[unit_index, :d_x] + unit_data[unit_index, :unit_life]
            total_capex = 0    # Only updated once project is complete
            cap_pmt = 0
            assets_vals = (next_id, agent_id, unit_type, revealed, completion_pd, cancellation_pd, retirement_pd, total_capex, cap_pmt)
            DBInterface.execute(db, "INSERT INTO assets VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", assets_vals)
        end
    end
end

##### Authorize ANPE for all current WIP projects
# Retrieve all WIP projects
WIP_projects = get_WIP_projects_list(db, pd, agent_id)
# Authorize ANPE for the upcoming period (default: $1B/year)
authorize_anpe(db, agent_id, pd, WIP_projects, unit_data)

# End
println("\n Julia: finishing")
println("\n-----------------------------------------------------------")

