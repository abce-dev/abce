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
pd = CLI_args["current_pd"]
agent_id = CLI_args["agent_id"]

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
@info unit_data
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

PA_uids, PA_fs_dict = set_up_project_alternatives(unit_data, asset_counts, num_lags, fc_pd, agent_params, price_curve, db, pd)

@info "Project alternatives:"
@info PA_uids

if pd == 0
    @info "Unit data loaded:"
    @info unit_data
end

@info "Data initialized."

###### Set up the model
m = set_up_model(settings, PA_uids, PA_fs_dict, available_demand, asset_counts)

###### Solve the model
@info "Solving optimization problem..."
optimize!(m)
status = termination_status.(m)
# This MILP should always return integral solutions; convert the float values
#   to integers to avoid some future TypeErrors
unit_qty = Int.(round.(value.(m[:u])))


###### Display the results
all_results = hcat(PA_uids, DataFrame(units_to_execute = unit_qty))
@info status
@info "Units to build:"
@info all_results


###### Save the new units into the `assets` and `WIP_projects` DB tables
postprocess_agent_decisions(all_results, unit_data, db, pd, agent_id)

##### Authorize ANPE for all current WIP projects
# Retrieve all WIP projects
WIP_projects = get_WIP_projects_list(db, pd, agent_id)
# Authorize ANPE for the upcoming period (default: $1B/year)
authorize_anpe(db, agent_id, pd, WIP_projects, unit_data)

# End
@info "Julia: finishing"
@info "-----------------------------------------------------------"



