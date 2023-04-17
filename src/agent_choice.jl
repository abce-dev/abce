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

@debug "-----------------------------------------------------------"
@debug "Julia agent choice algorithm: starting"
@debug "Loading packages..."
using JuMP, LinearAlgebra, DataFrames, CSV, YAML, SQLite, ArgParse

# Include local ABCE functions modules
include("ABCEfunctions.jl")
include("dispatch.jl")
include("C2N_projects.jl")
using .ABCEfunctions, .Dispatch, .C2N

CLI_args = ABCEfunctions.get_CL_args()

ABCEfunctions.set_verbosity(CLI_args["verbosity"])

# Load settings and file locations from the settings file
settings = YAML.load_file(CLI_args["settings_file"])

@debug "Julia modules loaded successfully."

###### Set up inputs
@info "Initializing data..."

settings = ABCEfunctions.set_up_local_paths(settings, CLI_args["abce_abs_path"])

# File names
db_file = joinpath(pwd(), "outputs", settings["simulation"]["scenario_name"], settings["file_paths"]["db_file"])
C2N_specs_file = joinpath(
                     settings["file_paths"]["ABCE_abs_path"],
                     "inputs",
                     "C2N_project_definitions.yml"
                 )
# Constants
hours_per_year = settings["constants"]["hours_per_year"]
num_lags = settings["agent_opt"]["num_future_periods_considered"]

# Load the inputs
db = ABCEfunctions.load_db(db_file)
pd = CLI_args["current_pd"]
agent_id = CLI_args["agent_id"]

# Load C2N specs data
C2N_specs = YAML.load_file(C2N_specs_file)

# Set up agent-specific data
# Get a list of all ongoing construction projects for the current agent
agent_projects = ABCEfunctions.get_WIP_projects_list(db, pd, agent_id)
# Get a list of all operating assets owned by the current agent
agent_assets, asset_counts = ABCEfunctions.get_current_assets_list(
                                 db,
                                 pd,
                                 agent_id
                             )

# Get agent financial parameters
agent_params = ABCEfunctions.get_agent_params(db, agent_id)

# System parameters
# Read unit operational data (unit_specs)
unit_specs = ABCEfunctions.get_unit_specs(db)

# Ensure that forecast horizon is long enough to accommodate the end of life
#   for the most long-lived possible unit
fc_pd = ABCEfunctions.set_forecast_period(unit_specs, num_lags)

@info "Data initialized."
@info "Setting up dispatch simulation..."

all_year_system_portfolios, all_year_agent_portfolios = Dispatch.set_up_dispatch_portfolios(db, pd, fc_pd, agent_id, unit_specs)

# Load the demand data
total_demand = ABCEfunctions.get_demand_forecast(db, pd, agent_id, fc_pd, settings)

# Extend the unserved demand data to match the total forecast period (constant projection)
total_demand = ABCEfunctions.get_net_demand(db, pd, agent_id, fc_pd, total_demand, all_year_system_portfolios, unit_specs)
@debug "Demand data:"
@debug total_demand[1:10, :]

@info "Running dispatch simulation..."
long_econ_results = Dispatch.execute_dispatch_economic_projection(db, settings, pd, fc_pd, total_demand, unit_specs, all_year_system_portfolios)
@info "Dispatch projections complete."

@info "Setting up project alternatives..."
PA_uids, PA_fs_dict = ABCEfunctions.set_up_project_alternatives(settings, unit_specs, asset_counts, num_lags, fc_pd, agent_params, db, pd, long_econ_results, C2N_specs)

@info "Project alternatives set up."

@debug "Project alternatives:"
@debug PA_uids

###### Set up the model
@info "Setting up the agent's decision optimization model..."
unified_agent_portfolios = Dispatch.combine_and_extend_year_portfolios(all_year_agent_portfolios, pd+fc_pd)

agent_fs = ABCEfunctions.update_agent_financial_statement(agent_id, db, unit_specs, pd, fc_pd, long_econ_results, unified_agent_portfolios)

m = ABCEfunctions.set_up_model(settings, PA_uids, PA_fs_dict, total_demand, asset_counts, agent_params, unit_specs, pd, all_year_system_portfolios, db, agent_id, agent_fs, fc_pd)

###### Solve the model
@info "Solving agent's decision optimization problem..."
optimize!(m)

all_results = ABCEfunctions.finalize_results_dataframe(m, PA_uids)

###### Display the results
ABCEfunctions.display_agent_choice_results(CLI_args, all_results)

###### Save the new units into the `assets` and `WIP_projects` DB tables
ABCEfunctions.postprocess_agent_decisions(settings, all_results, unit_specs, db, pd, agent_id)

