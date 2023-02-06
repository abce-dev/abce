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

# Set up command-line parser
s = ArgParseSettings()
@add_arg_table s begin
    "--settings_file"
        help = "absolute path to the settings file"
        required = false
        default = joinpath(pwd(), "settings.yml")
    "--agent_id"
        help = "unique ID number of the agent"
        arg_type = Int
        required = true
    "--current_pd"
        help = "current ABCE time period"
        arg_type = Int
        required = true
    "--verbosity"
        help = "level of output logged to the console"
        arg_type = Int
        required = false
        default = 1
        range_tester = x -> x in [0, 1, 2, 3]
end

# Retrieve parsed arguments from command line
CLI_args = parse_args(s)

lvl = 0
if CLI_args["verbosity"] == 0
    # Only show Logging messages of severity Error (level 2000) and above
    lvl = 2000
elseif CLI_args["verbosity"] == 1
    # Only show Logging messages of severity Warning (level 1000) and above
    lvl = 1000
elseif CLI_args["verbosity"] == 3
    # Show Logging messages of severity Debug (level -1000) and above
    lvl = -1000
end
global_logger(ConsoleLogger(lvl))

# Load settings and file locations from the settings file
settings_file = CLI_args["settings_file"]
settings = YAML.load_file(settings_file)

settings_file = CLI_args["settings_file"]
settings = YAML.load_file(settings_file)

# Include local ABCE functions module
julia_ABCE_module = "ABCEfunctions.jl"
include(julia_ABCE_module)
dispatch_module = "dispatch.jl"
include(dispatch_module)
C2N_module = "C2N_projects.jl"
include(C2N_module)
using .ABCEfunctions, .Dispatch, .C2N

@debug "Julia modules loaded successfully."

###### Set up inputs
@info "Initializing data..."

settings = ABCEfunctions.set_up_local_paths(settings)

solver = lowercase(settings["simulation"]["solver"])
@debug string("Solver is `$solver`")
if solver == "cplex"
    try
        using CPLEX
    catch LoadError
        throw(error("CPLEX is not available!"))
    end
elseif solver == "glpk"
    using GLPK
elseif solver == "cbc"
    using Cbc
elseif solver == "highs"
    using HiGHS
else
    throw(error("Solver `$solver` not supported. Try `cplex` instead."))
end

# File names
db_file = joinpath(pwd(), settings["file_paths"]["db_file"])
C2N_specs_file = joinpath(
                     @__DIR__,
                     "inputs",
                     "C2N_project_definitions.yml"
                 )
# Constants
hours_per_year = settings["constants"]["hours_per_year"]
consider_future_projects = settings["agent_opt"]["consider_future_projects"]
if consider_future_projects
    num_lags = settings["agent_opt"]["num_future_periods_considered"]
else
    num_lags = 0
end

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
agent_assets, asset_counts = ABCEfunctions.get_current_assets_list(db, pd, agent_id)

# Get agent financial parameters
agent_params = ABCEfunctions.get_agent_params(db, agent_id)

# System parameters
# Read unit operational data (unit_specs) and number of unit types (num_types)
unit_specs, num_types = ABCEfunctions.get_unit_specs(db)
num_alternatives = num_types * (num_lags + 1)

# Ensure that forecast horizon is long enough to accommodate the end of life
#   for the most long-lived possible unit
fc_pd = ABCEfunctions.set_forecast_period(unit_specs, num_lags)

# Add empty column for project NPVs in unit_specs
unit_specs[!, :FCF_NPV] = zeros(Float64, num_types)

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
long_econ_results = Dispatch.execute_dispatch_economic_projection(db, settings, pd, fc_pd, total_demand, unit_specs, all_year_system_portfolios, solver)
@info "Dispatch projections complete."

@info "Setting up project alternatives..."
PA_uids, PA_fs_dict = ABCEfunctions.set_up_project_alternatives(settings, unit_specs, asset_counts, num_lags, fc_pd, agent_params, db, pd, long_econ_results, C2N_specs)

@info "Project alternatives set up."

@debug "Project alternatives:"
@debug PA_uids

###### Set up the model
@info "Setting up the agent's decision optimization model..."
unified_agent_portfolios = Dispatch.create_all_year_portfolios(all_year_agent_portfolios, fc_pd, pd)

agent_fs = ABCEfunctions.update_agent_financial_statement(agent_id, db, unit_specs, pd, fc_pd, long_econ_results, unified_agent_portfolios)

m = ABCEfunctions.set_up_model(settings, PA_uids, PA_fs_dict, total_demand, asset_counts, agent_params, unit_specs, pd, all_year_system_portfolios, db, agent_id, agent_fs, fc_pd)

###### Solve the model
@info "Solving agent's decision optimization problem..."
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
short_results = filter(:units_to_execute => u -> u > 0, all_results)
@info status
if CLI_args["verbosity"] == 2
    @info "Alternatives to execute:"
    @info short_results
elseif CLI_args["verbosity"] == 3
    @debug "Alternatives to execute:"
    @debug all_results
end

###### Save the new units into the `assets` and `WIP_projects` DB tables
ABCEfunctions.postprocess_agent_decisions(all_results, unit_specs, db, pd, agent_id)

##### Authorize ANPE for all current WIP projects
# Retrieve all WIP projects
#WIP_projects = get_WIP_projects_list(db, pd, agent_id)
# Authorize ANPE for the upcoming period (default: $1B/year)
#authorize_anpe(db, agent_id, pd, WIP_projects, unit_specs)


