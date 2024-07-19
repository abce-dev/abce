##########################################################################
# Copyright 2023 Argonne National Laboratory
#
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

# Agent decision model

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

    settings =
        ABCEfunctions.set_up_local_paths(settings, CLI_args["abce_abs_path"])

    # File names
    db_file = joinpath(
        pwd(),
        "outputs",
        settings["simulation"]["scenario_name"],
        settings["file_paths"]["db_file"],
    )
    C2N_specs_file = joinpath(
        CLI_args["inputs_path"],
        "C2N_project_definitions.yml",
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


function process_results(settings, CLI_args, final_model, final_mode, db, PA_uids, unit_specs)
    # Ensure model results data is valid and of correct type
    all_results = ABCEfunctions.finalize_results_dataframe(final_model, final_mode, PA_uids)

    # Display the results
    ABCEfunctions.display_agent_choice_results(CLI_args, final_model, all_results)

    # Save newly-selected project alternatives happening in the current period
    #   to the database
    ABCEfunctions.postprocess_agent_decisions(
        settings,
        all_results,
        unit_specs,
        db,
        CLI_args["current_pd"],
        CLI_args["agent_id"],
    )
end


function run_agent_choice()
    @info "Setting up data..."

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
        CLI_args["agent_id"],
    )

    # Ensure that forecast horizon is long enough to accommodate the end of life
    #   for the most long-lived possible unit
    fc_pd = ABCEfunctions.set_forecast_period(
        unit_specs,
        settings["agent_opt"]["num_future_periods_considered"],
    )

    # Retrieve the year-by-year system generation portfolio based on currently
    #   available data
    system_portfolios = ABCEfunctions.get_portfolio_forecast(
        db,
        settings,
        CLI_args["current_pd"],
        unit_specs,
    )

    adj_system_portfolios = ABCEfunctions.fill_portfolios_missing_units(
        CLI_args["current_pd"],
        deepcopy(system_portfolios),
        unit_specs,
    )

    # Retrieve the year-by-year projected portfolio for the current agent
    agent_portfolios = ABCEfunctions.get_portfolio_forecast(
        db,
        settings,
        CLI_args["current_pd"],
        unit_specs,
        agent_id=CLI_args["agent_id"],
    )

    # Load the demand data
    demand_forecast = ABCEfunctions.get_demand_forecast(
        db,
        CLI_args["current_pd"],
        fc_pd,
        settings,
    )

    adj_system_portfolios = ABCEfunctions.forecast_balance_of_market_investment(
        db,
        adj_system_portfolios,
        agent_portfolios,
        agent_params,
        CLI_args["current_pd"],
        settings,
        demand_forecast,
    )

    # Save all system portfolio forecasts to the cnerg groupspace
    pfs = deepcopy(adj_system_portfolios[CLI_args["current_pd"]])
    for i=CLI_args["current_pd"]+1:maximum(keys(adj_system_portfolios))
        append!(pfs, adj_system_portfolios[i])
    end

    za = CLI_args["current_pd"]
    zb = CLI_args["agent_id"]
    filename = joinpath(
        settings["file_paths"]["output_logging_dir"],
        settings["simulation"]["scenario_name"],
        string("agent_", zb, "_pd_", za, "_pf_forecast.csv"),
    )
    CSV.write(filename, pfs)


    # Use the agent's internal dispatch forecast generator to project dispatch
    #   results in the system over the forecast horizon
    @info "Simulating future market dispatch..."
    long_econ_results, dispatch_results = Dispatch.execute_dispatch_economic_projection(
        CLI_args,
        db,
        settings,
        fc_pd,
        demand_forecast,
        unit_specs,
        adj_system_portfolios;
        run_mode="forecast",
        downselection_mode=settings["dispatch"]["downselection"]
    )

    # Save all agents' long econ results to the cnerg groupspace
    za = CLI_args["current_pd"]
    zb = CLI_args["agent_id"]
    CSV.write(
        joinpath(
            settings["file_paths"]["output_logging_dir"],
            settings["simulation"]["scenario_name"],
            string("agent_", zb, "_pd_", za, "_long_econ_results.csv"),
        ),
        long_econ_results,
    )


    # Set up all available project alternatives, including computing marginal
    #   NPV for all potential projects (new construction and retirements)
    @info "Setting up all project alternatives..."
    PA_uids, PA_fs_dict = ABCEfunctions.set_up_project_alternatives(
        settings,
        unit_specs,
        grouped_agent_assets,
        fc_pd,
        CLI_args["agent_id"],
        agent_params,
        db,
        CLI_args["current_pd"],
        C2N_specs,
        dispatch_results,
        CLI_args["verbosity"],
    )

    # Update the agent's baseline projected financial statements, to use in
    #   the decision optimization model
    @info "Generating the agent's projected financial statements..."
    agent_fs = ABCEfunctions.forecast_agent_financial_statement(
        settings,
        CLI_args["agent_id"],
        db,
        unit_specs,
        CLI_args["current_pd"],
        fc_pd,
        dispatch_results,
        agent_params,
    )

    # Set up the agent's decision optimization model
    @info "Setting up the agent's decision optimization problem..."
    m = ABCEfunctions.set_up_model(
        settings,
        CLI_args,
        PA_uids,
        PA_fs_dict,
        demand_forecast,
        grouped_agent_assets,
        agent_params,
        unit_specs,
        CLI_args["current_pd"],
        adj_system_portfolios,
        db,
        CLI_args["agent_id"],
        agent_fs,
        fc_pd,
    )

    # Solve the model
    @info "Solving optimization problem..."
    m = ABCEfunctions.solve_model(m, CLI_args["verbosity"])

    status = string(termination_status.(m))
    if status == "OPTIMAL"
        final_model = m
        final_mode = "normal"
    else
        m_ret = ABCEfunctions.set_up_model(
            settings,
            CLI_args,
            PA_uids,
            PA_fs_dict,
            demand_forecast,
            grouped_agent_assets,
            agent_params,
            unit_specs,
            CLI_args["current_pd"],
            adj_system_portfolios,
            db,
            CLI_args["agent_id"],
            agent_fs,
            fc_pd;
            mode="ret_only"
        )

        m_ret = ABCEfunctions.solve_model(m_ret, CLI_args["verbosity"])

        final_model = m_ret
        final_mode = "ret_only"
    end

    # Process the model outputs
    @debug "Postprocessing model results..."
    process_results(settings, CLI_args, final_model, final_mode, db, PA_uids, unit_specs)
end


run_agent_choice()
