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

using Logging

@info "-----------------------------------------------------------"
@info "Julia agent choice algorithm: starting"
@info "Loading packages..."
using JuMP, GLPK, LinearAlgebra, DataFrames, CSV, YAML, SQLite, ArgParse

# Set up command-line parser
s = ArgParseSettings()
@add_arg_table s begin
    "--settings_file"
        help = "absolute path to the settings file"
        required = true
    "--agent_id"
        help = "unique ID number of the agent"
        arg_type = Int
        required = true
    "--current_pd"
        help = "current ABCE time period"
        arg_type = Int
        required = true
end

# Retrieve parsed arguments from command line
CLI_args = parse_args(s)

# Load settings and file locations from the settings file
settings_file = CLI_args["settings_file"]
settings = YAML.load_file(settings_file)

# Include local ABCE functions module
julia_ABCE_module = joinpath(settings["ABCE_abs_path"], "ABCEfunctions.jl")
include(julia_ABCE_module)
using .ABCEfunctions

@info "Packages loaded successfully."

###### Set up inputs
@info "Initializing data..."

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
#pd = get_current_period()
pd = CLI_args["current_pd"]
agent_id = CLI_args["agent_id"]
#agent_id = get_agent_id()

# Set up agent-specific data
# Get a list of all ongoing construction projects for the current agent
agent_projects = get_WIP_projects_list(db, pd, agent_id)
# Get a list of all operating assets owned by the current agent
agent_assets, asset_counts = get_current_assets_list(db, pd, agent_id)

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
@info "Available demand:"
@info DataFrame(net_demand = available_demand)[1:10, :]

# Load the price data
price_curve = DBInterface.execute(db, "SELECT * FROM price_curve") |> DataFrame

# Add empty column for project NPVs in unit_data
unit_data[!, :FCF_NPV] = zeros(Float64, num_types)

# NPV is the only decision criterion, so create a dataframe to hold the results
#    for each alternative
alternative_names, NPV_results = create_NPV_results_df(unit_data, num_lags)

# Create per-unit financial statement tables
@info "Creating and populating unit financial statements for NPV calculation"
unit_FS_dict = create_FS_dict(unit_data, fc_pd, num_lags)

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

        # Generate "prime movers" (debt payments and depreciation)
        generate_prime_movers(unit_type_data, fs, j, agent_params[1, :cost_of_debt])

        # Forecast unit revenue ($/period) and generation (kWh/period)
        forecast_unit_revenue_and_gen(unit_type_data, fs, price_curve, db, pd, j)

        # Forecast unit costs: fuel cost, VOM, and FOM
        forecast_unit_op_costs(unit_type_data, fs, j)

        # Propagate the accounting logic (EBITDA --> FCF)
        propagate_accounting_line_items(fs, db)

        # Compute this unit alternative's FCF NPV
        FCF_NPV = compute_alternative_NPV(fs, agent_params)

        # Save the NPV result
        NPV_results[findall(NPV_results.name .== name)[1], :NPV] = FCF_NPV
        unit_data[i, :FCF_NPV] = FCF_NPV
    end
end

@info "xtr NPV results:"
@info NPV_results


# Create a dataframe to hold NPV results for each retirement alternative
ret_alt_names, ret_NPV_results = create_NPV_results_df(asset_counts, num_lags; mode="retire")

# Create a dataframe to store results for retirement NPV calculations
ret_FS_dict = create_FS_dict(asset_counts, fc_pd, num_lags; mode="retire")

# Compute dataframes for retiring existing assets
for i = 1:size(asset_counts)[1]
    for j = 0:num_lags
        name = string(asset_counts[i, :unit_type], "_", asset_counts[i, :retirement_pd], "_lag-", j)
        original_ret_pd = asset_counts[i, :retirement_pd]
        fs = ret_FS_dict[name]
        unit_type_data = filter(row -> row.unit_type == asset_counts[i, :unit_type], unit_data)

        # Implies any retiring unit is 100% paid off; need to implement tracking of debt repayments

        # Forecast unit revenue ($/period) and generation (kWh/period)
        forecast_unit_revenue_and_gen(unit_type_data, fs, price_curve, db, pd, j; mode="retire", orig_ret_pd=original_ret_pd)

        # Forecast unit costs: fuel cost, VOM, and FOM
        forecast_unit_op_costs(unit_type_data, fs, j; mode="retire", orig_ret_pd=original_ret_pd)

        # Convert to marginal deltas
        convert_to_marginal_delta_FS(fs, j)

        # Propagate the accounting logic (EBITDA --> FCF)
        propagate_accounting_line_items(fs, db)

        # Compute this unit alternative's FCF NPV
        FCF_NPV = compute_alternative_NPV(fs, agent_params)

        # Save the NPV result
        ret_NPV_results[findall(ret_NPV_results.name .== name)[1], :NPV] = FCF_NPV
    end
end



@info "ret NPV results:"
@info ret_NPV_results

if pd == 0
    @info "Unit data loaded:"
    @info unit_data
end

@info "Data initialized."

###### Set up the model
m = set_up_model(unit_FS_dict, ret_FS_dict, available_demand, NPV_results, ret_NPV_results, asset_counts)

###### Solve the model
@info "Solving optimization problem..."
optimize!(m)
status = termination_status.(m)
unit_qty = value.(m[:u])


###### Display the results
all_alternatives = vcat(NPV_results, ret_NPV_results)[!, :name]
all_results = hcat(all_alternatives, DataFrame(units = unit_qty))
@info status
@info "Units to build:"
@info all_results


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
@info "Julia: finishing"
@info "-----------------------------------------------------------"



