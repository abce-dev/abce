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

using Logging, JuMP, LinearAlgebra, DataFrames, CSV, YAML, SQLite, ArgParse

# Include local ABCE functions modules
include("ABCEfunctions.jl")
include("dispatch.jl")
include("C2N_projects.jl")
using .ABCEfunctions, .Dispatch, .C2N



function set_up_run(CLI_args)
    ABCEfunctions.set_verbosity(CLI_args["verbosity"])

    # Load settings and file locations from the settings file
    settings = YAML.load_file(CLI_args["settings_file"])

    settings = ABCEfunctions.set_up_local_paths(settings, CLI_args["abce_abs_path"])

    # File names
    db_file = joinpath(pwd(), "outputs", settings["simulation"]["scenario_name"], settings["file_paths"]["db_file"])
    C2N_specs_file = joinpath(
                         settings["file_paths"]["ABCE_abs_path"],
                         "inputs",
                         "C2N_project_definitions.yml"
                     )

    # Load the database
    db = ABCEfunctions.load_db(db_file)

    # Load C2N specs data
    C2N_specs = YAML.load_file(C2N_specs_file)

    return settings, db, C2N_specs
end


function get_raw_db_data(db, CLI_args)
    # Get agent financial parameters
    agent_params = ABCEfunctions.get_agent_params(db, CLI_args["agent_id"])

    # System parameters
    # Read unit operational data (unit_specs)
    unit_specs = ABCEfunctions.get_unit_specs(db)

    return agent_params, unit_specs
end


function process_results(settings, CLI_args, m, db, PA_uids, unit_specs)
    # Ensure model results data is valid and of correct type
    all_results = ABCEfunctions.finalize_results_dataframe(m, PA_uids)

    # Display the results
    ABCEfunctions.display_agent_choice_results(CLI_args, all_results)

    # Save newly-selected project alternatives happening in the current period
    #   to the database
    ABCEfunctions.postprocess_agent_decisions(
        settings,
        all_results,
        unit_specs,
        db,
        CLI_args["agent_id"],
        CLI_args["current_pd"]
    )
end


function run_agent_choice()
    # Read in the command-line arguments
    CLI_args = ABCEfunctions.get_CL_args()

    # Read in data and the database from file
    settings, db, C2N_specs = set_up_run(CLI_args)

    # Read in some raw data from the database
    agent_params, unit_specs = get_raw_db_data(db, CLI_args)

    # Retrieve a list of the agent's currently-operating assets, grouped by
    #   type and mandatory retirement date
    grouped_agent_assets = ABCEfunctions.get_grouped_current_assets(
                               db,
                               CLI_args["current_pd"],
                               CLI_args["agent_id"]
                           )

    # Ensure that forecast horizon is long enough to accommodate the end of life
    #   for the most long-lived possible unit
    fc_pd = ABCEfunctions.set_forecast_period(
                unit_specs,
                settings["agent_opt"]["num_future_periods_considered"]
            )

    # Retrieve the year-by-year system generation portfolio based on currently
    #   available data
    system_portfolios = Dispatch.get_system_portfolios(
                            db,
                            CLI_args["current_pd"],
                            fc_pd,
                            unit_specs
                        )

    # Load the demand data
    total_demand = ABCEfunctions.get_demand_forecast(
                       db,
                       CLI_args["current_pd"],
                       fc_pd,
                       settings
                   )

    # Extend the unserved demand data to match the total forecast period (constant projection)
    total_demand = ABCEfunctions.get_net_demand(
                       db,
                       CLI_args["current_pd"],
                       fc_pd,
                       total_demand,
                       system_portfolios,
                       unit_specs
                   )

    # Use the agent's internal dispatch forecast generator to project dispatch
    #   results in the system over the forecast horizon
    long_econ_results = Dispatch.execute_dispatch_economic_projection(
                            db,
                            settings,
                            CLI_args["current_pd"],
                            fc_pd,
                            total_demand,
                            unit_specs,
                            system_portfolios
                        )

    # Set up all available project alternatives, including computing marginal
    #   NPV for all potential projects (new construction and retirements)
    PA_uids, PA_fs_dict = ABCEfunctions.set_up_project_alternatives(
                              settings,
                              unit_specs,
                              grouped_agent_assets,
                              fc_pd,
                              agent_params,
                              db,
                              CLI_args["current_pd"],
                              long_econ_results,
                              C2N_specs
                          )

    # Update the agent's baseline projected financial statements, to use in
    #   the decision optimization model
    agent_fs = ABCEfunctions.update_agent_financial_statement(
                   CLI_args["agent_id"],
                   db,
                   unit_specs,
                   CLI_args["current_pd"],
                   fc_pd,
                   long_econ_results
               )

    # Set up the agent's decision optimization model
    m = ABCEfunctions.set_up_model(
            settings,
            PA_uids,
            PA_fs_dict,
            total_demand,
            grouped_agent_assets,
            agent_params,
            unit_specs,
            CLI_args["current_pd"],
            system_portfolios,
            db,
            CLI_args["agent_id"],
            agent_fs,
            fc_pd
        )

    # Solve the model
    optimize!(m)

    # Process the model outputs
    process_results(settings, CLI_args, m, db, PA_uids, unit_specs)
end


run_agent_choice()


