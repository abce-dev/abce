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
using .ABCEfunctions, .Dispatch


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

    # Load the database
    db = ABCEfunctions.load_db(db_file)

    return settings, db
end


function get_raw_db_data(db, CLI_args)
    # Get agent financial parameters
    agent_params = ABCEfunctions.get_agent_params(db, CLI_args["agent_id"])

    # System parameters
    # Read unit operational data (unit_specs)
    unit_specs = ABCEfunctions.get_unit_specs(db)

    return agent_params, unit_specs
end


function compute_last_year_results(db, settings, CLI_args, agent_params)
    @info "Computing realized financial results for last year..."
    id = CLI_args["agent_id"]
    y = CLI_args["current_pd"] - 1
    unit_fin_results = DBInterface.execute(db, "SELECT * FROM annual_dispatch_unit_summary WHERE period = $y") |> DataFrame
    pf = DBInterface.execute(db, "SELECT unit_type, COUNT(unit_type) FROM assets WHERE agent_id = $id AND completion_pd <= $y AND retirement_pd > $y AND cancellation_pd > $y GROUP BY unit_type") |> DataFrame
    pf = rename(pf, "COUNT(unit_type)" => "num_units")

    pf_res = innerjoin(pf, unit_fin_results, on = "unit_type")

    ufs = select(pf_res, ["unit_type", "num_units", "revenue", "VOM", "fuel_cost", "FOM", "carbon_tax", "tax_credits"])
    for col in filter!(e -> !(e in ["unit_type", "num_units"]), names(ufs))
        transform!(
            ufs, 
            [col, "num_units"] =>
            ((d, num_units) -> d .* num_units)
            => col
        )
    end

    # Ensure all costs are negative numbers
    for col in ["VOM", "FOM", "fuel_cost", "carbon_tax"]
        transform!(
            ufs,
            [col] =>
            (c -> -1 .* abs.(c))
            => col
        )
    end

    # Sum over all unit types to determine total operational results
    revtot = sum(ufs[!, "revenue"])
    VOMtot = sum(ufs[!, "VOM"])
    fctot = sum(ufs[!, "fuel_cost"])
    FOMtot = sum(ufs[!, "FOM"])
    ctaxtot = sum(ufs[!, "carbon_tax"])
    tctot = sum(ufs[!, "tax_credits"])

    push!(ufs, ["total" 9999 revtot VOMtot fctot FOMtot ctaxtot tctot])

    fs = filter("unit_type" => u -> u == "total", ufs)

    fs[!, "base_pd"] .= CLI_args["current_pd"] - 1
    fs[!, "projected_pd"] .= CLI_args["current_pd"] - 1

    select!(fs, ["base_pd", "projected_pd", "revenue", "VOM", "fuel_cost", "FOM", "carbon_tax", "tax_credits"])

    for col in names(fs)
        rename!(fs, col => Symbol(col))
    end

    # Add in scheduled financing factors: depreciation, interest payments, and capex
    fs = ABCEfunctions.compute_scheduled_financing_factors(db, fs, CLI_args["agent_id"], y)

    # Compute accounting line items
    fs = ABCEfunctions.compute_accounting_line_items(db, y, settings, fs, agent_params, CLI_args["agent_id"])
    fs = ABCEfunctions.compute_credit_indicator_scores(settings, fs)

    ABCEfunctions.save_agent_fs!(fs, CLI_args["agent_id"], db, "realized")


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


function save_intermediate_outputs(settings, CLI_args, adj_system_portfolios, long_econ_results)
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

end


function run_agent_choice()
    @info "Setting up data..."

    # Read in the command-line arguments
    CLI_args = ABCEfunctions.get_CL_args()

    # Read in data and the database from file
    settings, db = set_up_run(CLI_args)

    # Read in some raw data from the database
    agent_params, unit_specs = get_raw_db_data(db, CLI_args)

    # Set last year's realized financial results
    if CLI_args["current_pd"] != 0
        compute_last_year_results(db, settings, CLI_args, agent_params)
    end

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
        system_portfolios,
        agent_portfolios,
        agent_params,
        CLI_args["current_pd"],
        settings,
        demand_forecast,
    )


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

    if settings["simulation"]["file_logging_level"] > 0
        save_intermediate_outputs(settings, CLI_args, adj_system_portfolios, long_econ_results)
    end

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
        dispatch_results,
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
