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
using JuMP, GLPK, LinearAlgebra, DataFrames, CSV, YAML, SQLite, ArgParse, CPLEX

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
dispatch_module = joinpath(settings["ABCE_abs_path"], "dispatch.jl")
include(dispatch_module)
using .ABCEfunctions, .Dispatch

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
pd = CLI_args["current_pd"]
agent_id = CLI_args["agent_id"]

# Set up agent-specific data
# Get a list of all ongoing construction projects for the current agent
agent_projects = get_WIP_projects_list(db, pd, agent_id)
# Get a list of all operating assets owned by the current agent
agent_assets, asset_counts = get_current_assets_list(db, pd, agent_id)

# Get agent financial parameters
agent_params = get_agent_params(db, agent_id)

# System parameters
# Read unit operational data (unit_specs) and number of unit types (num_types)
unit_specs, num_types = get_unit_specs(db)
if pd == 0
    @info unit_specs
end
num_alternatives = num_types * (num_lags + 1)

# Ensure that forecast horizon is long enough to accommodate the end of life
#   for the most long-lived possible unit
fc_pd = set_forecast_period(unit_specs, num_lags)

# Load the demand data
total_demand = get_demand_forecast(db, pd, agent_id, fc_pd, settings)

# Extend the unserved demand data to match the total forecast period (constant projection)
available_demand = get_net_demand(db, pd, agent_id, fc_pd, total_demand)
@info "Available demand:"
@info DataFrame(net_demand = available_demand)[1:10, :]

# Load the price data
price_curve = DBInterface.execute(db, "SELECT * FROM price_curve") |> DataFrame

# Add empty column for project NPVs in unit_specs
unit_specs[!, :FCF_NPV] = zeros(Float64, num_types)

@info "Data initialized."

# Set up portfolio projections
@info "Setting up dispatch portfolios..."
system_portfolios = Dict()
for y = pd:settings["num_dispatch_years"]+pd
    # Retrieve a list of all units expected to be operational during this year,
    #   grouped by unit type
    year_portfolio = DBInterface.execute(db, "SELECT unit_type, COUNT(unit_type) FROM assets WHERE completion_pd <= $y AND retirement_pd > $y AND cancellation_pd > $y GROUP BY unit_type") |> DataFrame

    # Rename the COUNT(unit_type) column to num_units
    rename!(year_portfolio, Symbol("COUNT(unit_type)") => :num_units)

    # Add the year-i portfolio to the dictionary
    system_portfolios[y] = year_portfolio
    @debug y
    @debug system_portfolios[y]
end

@info "Dispatch portfolios set up."

# Run dispatch for each forecast year
@info string("Running the dispatch simulation for ", settings["num_dispatch_years"], " years...")
num_repdays = settings["num_repdays"]
ts_data = Dispatch.load_ts_data(
              joinpath(settings["ABCE_abs_path"],
                       "inputs",
                       "ALEAF_inputs"
              ),
              num_repdays
          )

# Set up dataframes to record all results
all_prices, all_gc_results = Dispatch.set_up_results_dfs()

run_next_year = true
y = pd
while (run_next_year) && (y <=settings["num_dispatch_years"])
    @info "Start of year $y dispatch simulation..."
    # Retrieve the current year's expected portfolio
    year_portfolio = system_portfolios[y]

    # Set up and run the dispatch simulation for this year
    # total_demand is 0-indexed
    global run_next_year = Dispatch.run_annual_dispatch(y, year_portfolio, total_demand[y+1, :demand], ts_data, unit_specs, all_gc_results, all_prices)
    @info "Year $y dispatch complete."

    # Manually increment the year counter
    global y = y + 1
end

if size(all_gc_results)[1] == 0
    throw(ErrorException("No dispatch simulations were able to be run. Re-check your inputs."))
end

# Propagate the results dataframes out to the end of the projection horizon
# Assume no change after the last modeled year
all_gc_results, all_prices = Dispatch.propagate_all_results(fc_pd-1, all_gc_results, all_prices, pd)

Dispatch.save_raw_results(all_prices, all_gc_results)

@info "Postprocessing dispatch simulation..."
long_econ_results = Dispatch.postprocess_results(system_portfolios, all_prices, all_gc_results, ts_data, unit_specs, fc_pd, pd)

CSV.write("long_econ_results_$pd.csv", long_econ_results)

@info "Setting up project alternatives..."
PA_uids, PA_fs_dict = set_up_project_alternatives(unit_specs, asset_counts, num_lags, fc_pd, agent_params, price_curve, db, pd, long_econ_results, settings["allowed_xtr_types"])

@info "Project alternatives:"
@info PA_uids

###### Set up the model
system_portfolios = Dict()
agent_system_portfolios = Dict()
for i = pd:pd+fc_pd-1
    # Retrieve a list of all units expected to be operational during this year,
    #   grouped by unit type
    year_portfolio = DBInterface.execute(db, "SELECT unit_type, COUNT(unit_type) FROM assets WHERE completion_pd <= $i AND retirement_pd > $i AND cancellation_pd > $i GROUP BY unit_type") |> DataFrame
    agent_year_portfolio = DBInterface.execute(db, "SELECT unit_type, COUNT(unit_type) FROM assets WHERE agent_id = $agent_id AND completion_pd <= $i AND retirement_pd > $i AND cancellation_pd > $i GROUP BY unit_type") |> DataFrame

    # Rename the COUNT(unit_type) column to num_units
    rename!(year_portfolio, Symbol("COUNT(unit_type)") => :num_units)
    rename!(agent_year_portfolio, Symbol("COUNT(unit_type)") => :num_units)

    year_portfolio = innerjoin(unit_specs, year_portfolio, on = :unit_type)

    year_portfolio = select(year_portfolio, [:unit_type, :num_units, :CF, :capacity])

    transform!(year_portfolio, [:num_units, :CF, :capacity] => ((num_units, CF, cap) -> num_units .* CF .* cap) => :effective_capacity)

    # Add the year-i portfolio to the dictionary
    system_portfolios[i] = year_portfolio
    agent_system_portfolios[i] = agent_year_portfolio

end

agent_all_year_portfolios = Dispatch.create_all_year_portfolios(agent_system_portfolios, fc_pd, pd)

agent_fs = update_agent_financial_statement(agent_id, db, unit_specs, pd, fc_pd, long_econ_results, agent_all_year_portfolios)

m = set_up_model(settings, PA_uids, PA_fs_dict, available_demand, asset_counts, agent_params, unit_specs, pd, total_demand[1, :demand], system_portfolios, db, agent_id, agent_fs, fc_pd)

###### Solve the model
@info "Solving optimization problem..."
optimize!(m)
status = string(termination_status.(m))
if status == "OPTIMAL"
    # This MILP should always return integral solutions; convert the float values
    #   to integers to avoid some future TypeErrors
    unit_qty = Int.(round.(value.(m[:u])))
else
    # If the agent has no valid options, do nothing
    unit_qty = zeros(Int64, size(PA_uids)[1])
end

###### Display the results
all_results = hcat(PA_uids, DataFrame(units_to_execute = unit_qty))
@info status
@info "Alternatives to execute:"
@info all_results


###### Save the new units into the `assets` and `WIP_projects` DB tables
postprocess_agent_decisions(all_results, unit_specs, db, pd, agent_id)

##### Authorize ANPE for all current WIP projects
# Retrieve all WIP projects
#WIP_projects = get_WIP_projects_list(db, pd, agent_id)
# Authorize ANPE for the upcoming period (default: $1B/year)
#authorize_anpe(db, agent_id, pd, WIP_projects, unit_specs)

# End
@info "Julia: finishing"
@info "-----------------------------------------------------------"



