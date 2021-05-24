# Agent decision model

println("\n-----------------------------------------------------------")
println("Julia agent choice algorithm: starting")
println("Loading packages...")
using JuMP, GLPK, LinearAlgebra, DataFrames, CSV, Printf, YAML, SQLite
# Include localy module of ABCE functions
include("./ABCEfunctions.jl")
using .ABCEfunctions
println("Packages loaded successfully.")

###### Set up inputs
println("Initializing data...")

# Load settings and file locations from the settings file
settings_file = ARGS[1]
settings = YAML.load_file(settings_file)
# File names
unit_specs_file = settings["unit_specs_file"]
fuel_cost_file = settings["fuel_data_file"]
demand_data_file = settings["demand_data_file"]
price_curve_data_file = settings["price_curve_data_file"]
db_file = settings["db_file"]
# Constants
hours_per_year = settings["hours_per_year"]
demand_vis_horizon = settings["demand_visibility_horizon"]
hist_demand_growth = settings["historical_demand_growth_rate"]
hist_demand_weight = settings["historical_demand_weight"]
consider_future_projects = settings["consider_future_projects"]
if consider_future_projects
    num_lags = settings["num_future_periods_considered"]
else
    num_lags = 0
end

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

# Set an average e- price
avg_e_price = 0.08    # $/kWh, avg electricity price

# System parameters
# Read unit operational data (unit_data) and number of unit types (num_types)
unit_data, num_types = get_unit_specs(db)
num_alternatives = num_types * (num_lags + 1)

# Ensure that forecast horizon is long enough to accommodate the end of life
#   for the most long-lived possible unit
fc_pd = set_forecast_period(unit_data, num_lags)

# Load the demand data
available_demand = get_demand_forecast(db, pd, demand_vis_horizon, agent_id, fc_pd, "exponential")

# Extend the unserved demand data to match the total forecast period (constant projection)
available_demand = get_net_demand(db, pd, agent_id, fc_pd, available_demand)
println(DataFrame(net_demand = available_demand)[1:10, :])

# Load the price data
price_curve = DBInterface.execute(db, "SELECT * FROM price_curve") |> DataFrame

# Add empty column for project NPVs in unit_data
unit_data[!, :FCF_NPV] = zeros(Float64, num_types)

# NPV is the only decision criterion, so create a dataframe to hold the results
#    for each alternative
alternative_names = Vector{String}()
for i = 1:num_types
    for j = 0:num_lags
        name = string(unit_data[i, :unit_type], "_lag-", j)
        push!(alternative_names, name)
    end
end
NPV_results = DataFrame(name = alternative_names, NPV = zeros(num_alternatives))

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
        name = string(unit_data[i, :unit_type], "_lag-", j)
        fs = unit_FS_dict[name]

        # Generate events during the construction period
        head_zeros_series = zeros(j)
        xtr_exp_per_pd = unit_data[i, :uc_x] * unit_data[i, :capacity] * 1000 / unit_data[i, :d_x]
        xtr_exp_series = ones(unit_data[i, :d_x]) .* xtr_exp_per_pd
        tail_zeros_series = zeros(fc_pd - j - unit_data[i, :d_x])
        xtr_exp = vcat(head_zeros_series, xtr_exp_series, tail_zeros_series)
        fs[!, :xtr_exp] .= xtr_exp
        for k = j+1:j+unit_data[i, :d_x]
            # Uniformly distribute construction costs over the construction duration
            fs[k, :remaining_debt_principal] = sum(fs[1:k, :xtr_exp]) * agent_params[1, :debt_fraction]
        end

        # Generate prime-mover events from the end of construction to the end of the unit's life
        for k = (j + unit_data[i, :d_x] + 1):(j + unit_data[i, :d_x] + unit_data[i, :unit_life])
            # Apply constant debt payment (sinking fund at cost of debt)
            # This amount is always calculated based on the amount of debt outstanding at
            #    the project's completion
            fs[k, :debt_payment] = fs[j + unit_data[i, :d_x], :remaining_debt_principal] .* d ./ (1 - (1+d) .^ (-1*unit_data[i, :unit_life]))
            # Determine portion of payment which pays down interest (instead of principal)
            fs[k, :interest_due] = fs[k-1, :remaining_debt_principal] * d
            # Update the amount of principal remaining at the end of the year
            fs[k, :remaining_debt_principal] = fs[k-1, :remaining_debt_principal] - (fs[k, :debt_payment] - fs[k, :interest_due])

            # Apply straight-line depreciation, based on debt outstanding at
            #    the project's completion
            fs[k, :depreciation] = fs[j + unit_data[i, :d_x], :xtr_exp] ./ unit_data[i, :unit_life]
        end

        # Compute unit revenue, based on the price duration curve loaded from file
        submarginal_hours = filter(row -> row.lamda > unit_data[i, :VOM] * 1000 + (unit_data[i, :FC_per_MMBTU] * unit_data[i, :heat_rate]i / 1000), price_curve)
        marginal_hours = filter(row -> row.lamda == unit_data[i, :VOM] * 1000 + (unit_data[i, :FC_per_MMBTU] * unit_data[i, :heat_rate] / 1000), price_curve)
        if size(marginal_hours)[1] != 0
            marginal_hours_revenue = sum(marginal_hours[!, :lamda]) * unit_data[i, :capacity] / (unit_data[i, :capacity] + size(marginal_hours)[1])
        else
            # There were no hours where this unit (or an equivalently-priced one) was marginal
            marginal_hours_revenue = 0
        end
        submarginal_hours_revenue = sum(submarginal_hours[!, :lamda]) * unit_data[i, :capacity]
        fs[!, :Revenue] .= 0.0
        fs[(j + unit_data[i, :d_x] + 1):(j + unit_data[i, :d_x] + unit_data[i, :unit_life]), :Revenue] .= (submarginal_hours_revenue + marginal_hours_revenue) * hours_per_year / size(price_curve)[1]

        # Unit generates during all marginal and sub-marginal hours
        num_active_hours = (size(submarginal_hours)[1] + size(marginal_hours)[1] * unit_data[i, :capacity] / (unit_data[i, :capacity] + size(marginal_hours)[1])) * hours_per_year / size(price_curve)[1]
        gen = num_active_hours * unit_data[i, :capacity] * 1000   # kWh
        fs[(j + unit_data[i, :d_x] + 1):(j + unit_data[i, :d_x] + unit_data[i, :unit_life]), :gen] .= gen


        # Apply reactive functions to the rest of the dataframe
        # Compute total cost of fuel used
        transform!(fs, [:gen] => ((gen) -> unit_data[i, :FC_per_MMBTU] .* unit_data[i, :heat_rate] ./ 1000000 .* gen) => :Fuel_Cost)
        # Compute total VOM cost incurred during generation
        transform!(fs, [:gen] => ((gen) -> unit_data[i, :VOM] .* gen) => :VOM_Cost)
        # Compute total FOM cost for the year (non-reactive)
        fs[!, :FOM_Cost] = zeros(size(fs)[1])
        fs[(j + unit_data[i, :d_x]+1):(j + unit_data[i, :d_x]+unit_data[i, :unit_life]), :FOM_Cost] .= unit_data[i, :FOM] * unit_data[i, :capacity] * 1000
        # Compute EBITDA
        transform!(fs, [:Revenue, :Fuel_Cost, :VOM_Cost, :FOM_Cost] => ((rev, fc, VOM, FOM) -> rev - fc - VOM - FOM) => :EBITDA)
        # Compute EBIT
        transform!(fs, [:EBITDA, :depreciation] => ((EBITDA, dep) -> EBITDA - dep) => :EBIT)
        # Compute EBT
        transform!(fs, [:EBIT, :interest_due] => ((EBIT, interest) -> EBIT - interest) => :EBT)
        # Compute taxes owed
        transform!(fs, [:EBT] => ((EBT) -> EBT .* agent_params[1, :tax_rate]) => :tax_owed)
        # Compute net income
        transform!(fs, [:EBT, :tax_owed] => ((EBT, tax) -> EBT - tax) => :Net_Income)
        # Compute FCF
        transform!(fs, [:Net_Income, :interest_due, :xtr_exp] => ((NI, interest, xtr_exp) -> NI + interest - xtr_exp) => :FCF)

        # Add column of compounded discount factors
        transform!(fs, [:year] => ((year) -> (1+d) .^ (-1 .* (year .- 1))) => :d_factor)
        # Discount unit FCF NPV values and save to the NPV_results dataframe
        NPV_results[findall(NPV_results.name .== name)[1], :NPV] = transpose(fs[!, :FCF]) * fs[!, :d_factor]
        unit_data[i, :FCF_NPV] = transpose(fs[!, :FCF]) * fs[!, :d_factor]
    end
end

#for i = 1:num_alternatives
#    if occursin("smr", alternative_names[i])
#        print(unit_FS_dict[alternative_names[i]])
#    end
#end
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
    @constraint(m, sum(u[j] * unit_FS_dict[alternative_names[j]][i, :gen] for j=1:num_alternatives) / (hours_per_year*1000) <= available_demand[i] * 2)
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
            rcec = unit_data[unit_index, :uc_x] * unit_data[unit_index, :capacity] * 1000
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

# Show final assets and WIP_projects tables
#show_table(db, "assets")
#show_table(db, "WIP_projects")

# End
println("\n Julia: finishing")
println("\n-----------------------------------------------------------")

