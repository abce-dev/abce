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
C2N_module = joinpath(settings["ABCE_abs_path"], "C2N_projects.jl")
include(C2N_module)
using .ABCEfunctions, .Dispatch, .C2N

@info "Packages loaded successfully."

###### Set up inputs
@info "Initializing data..."

# File names
db_file = joinpath(settings["ABCE_abs_path"], settings["db_file"])
C2N_specs_file = joinpath(
                     settings["ABCE_abs_path"],
                     "inputs",
                     "C2N_project_definitions.yml"
                 )
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

# Load C2N specs data
C2N_specs = YAML.load_file(C2N_specs_file)

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

# Load the price data
price_pd = pd - 1
price_curve = DBInterface.execute(db, "SELECT * FROM price_curve WHERE base_pd = $price_pd") |> DataFrame

# Add empty column for project NPVs in unit_specs
unit_specs[!, :FCF_NPV] = zeros(Float64, num_types)

@info "Data initialized."

all_year_system_portfolios, all_year_agent_portfolios = Dispatch.set_up_dispatch_portfolios(db, pd, fc_pd, agent_id, unit_specs)

# Load the demand data
total_demand = get_demand_forecast(db, pd, agent_id, fc_pd, settings)

# Extend the unserved demand data to match the total forecast period (constant projection)
total_demand = get_net_demand(db, pd, agent_id, fc_pd, total_demand, all_year_system_portfolios, unit_specs)
@info "Demand data:"
@info total_demand[1:10, :]

long_econ_results = Dispatch.execute_dispatch_economic_projection(db, settings, pd, fc_pd, total_demand, unit_specs, all_year_system_portfolios)

@info "Setting up project alternatives..."
PA_uids, PA_fs_dict = set_up_project_alternatives(settings, unit_specs, asset_counts, num_lags, fc_pd, agent_params, price_curve, db, pd, long_econ_results, settings["allowed_xtr_types"], C2N_specs)

@info "Project alternatives:"
@info PA_uids

###### Set up the model
unified_agent_portfolios = Dispatch.create_all_year_portfolios(all_year_agent_portfolios, fc_pd, pd)

agent_fs = update_agent_financial_statement(agent_id, db, unit_specs, pd, fc_pd, long_econ_results, unified_agent_portfolios, price_curve, settings)

m = set_up_model(settings, PA_uids, PA_fs_dict, total_demand, asset_counts, agent_params, unit_specs, pd, all_year_system_portfolios, db, agent_id, agent_fs, fc_pd)

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



