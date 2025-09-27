using ArgParse, Logging, YAML, DataFrames, CPLEX, CSV, SQLite

include(joinpath(ENV["ABCE_DIR"], "src", "ABCEfunctions.jl"))
using .ABCEfunctions

include(joinpath(ENV["ABCE_DIR"], "src", "dispatch.jl"))
using .Dispatch

function get_CL_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--current_pd"
        help = "current simulation period"
        required = true
        arg_type = Int

        "--ABCE_dir"
        help = "absolute path to the top-level ABCE directory"
        required = false
        default = ENV["ABCE_DIR"]

        "--settings_file"
        help = "absolute path to the settings file"
        required = false
        default = joinpath(ENV["ABCE_DIR"], "settings.yml")

        "--inputs_path"
        help = "absolute path to the input files"
        required = false
        default = joinpath(ENV["ABCE_DIR"], "inputs")
    end

    return parse_args(s)
end


function show_header()
    @info "==============================="
    @info "ABCE dispatch simulation module"
    @info "==============================="
end


function get_settings(CL_args)
    # Retrieve settings
    settings = YAML.load_file(CL_args["settings_file"])

    return settings
end

function get_db(CL_args, settings)
    # Load database
    db_file = joinpath(
        CL_args["ABCE_dir"], 
        "outputs",
        settings["simulation"]["scenario_name"], 
        settings["file_paths"]["db_file"],
    )
    db = ABCEfunctions.load_db(db_file)

    return db
end

function get_unit_specs(db)
    sql_cmd = "SELECT * FROM unit_specs"
    unit_specs = DBInterface.execute(db, sql_cmd) |> DataFrame

    return unit_specs
end


function get_year_portfolio(db, current_pd, unit_specs)
    brief_unit_specs = unit_specs[!, [:unit_type, :capacity, :capacity_factor]]

    sql_cmd = string(
        "SELECT unit_type, COUNT(unit_type) FROM assets ",
        "WHERE completion_pd <= $current_pd AND retirement_pd > $current_pd ",
        "AND cancellation_pd > $current_pd GROUP BY unit_type",
    )

    system_portfolio = DBInterface.execute(db, sql_cmd) |> DataFrame

    # Tidy up column names
    rename!(system_portfolio, Symbol("COUNT(unit_type)") => :num_units)

    # Add a column to indicate that all units are real
    system_portfolio[!, :real] = ones(size(system_portfolio)[1])
    system_portfolio[!, :esc_num_units] = system_portfolio[!, :num_units]

    # Join in the unit specs data
    system_portfolio = innerjoin(system_portfolio, brief_unit_specs, on = :unit_type)

    # Compute total capacity
    transform!(system_portfolio, [:num_units, :capacity] => ((num_units, cap) -> num_units .* cap) => :total_capacity)

    system_portfolio_dict = Dict()

    system_portfolio_dict[current_pd] = system_portfolio

    return system_portfolio_dict
end

function save_intermediate_outputs(settings, CL_args, long_econ_results, dispatch_results)
    @info "Saving annual dispatch results to the database..."

    CSV.write(
        joinpath(
            settings["file_paths"]["output_logging_dir"],
            settings["simulation"]["scenario_name"],
            string(
                "annual_dispatch_LER_pd_",
                CL_args["current_pd"],
                ".csv",
            ),
        ),
        long_econ_results,
    )

    CSV.write(
        joinpath(
            settings["file_paths"]["output_logging_dir"],
            settings["simulation"]["scenario_name"],
            string(
                "annual_dispatch_dispres_pd_",
                CL_args["current_pd"],
                ".csv",
            )
        ),
        dispatch_results,
    )
end


function run_true_annual_dispatch()
    show_header()

    @info "Setting up data..."
    CL_args = get_CL_args()
    settings = get_settings(CL_args)
    db = get_db(CL_args, settings)
    unit_specs = get_unit_specs(db)
    total_demand = ABCEfunctions.get_demand_forecast(db, CL_args["current_pd"], 1, settings)
    system_portfolio_dict = get_year_portfolio(db, CL_args["current_pd"], unit_specs)

    @info "Running the year-long dispatch simulation (this may take a little while)..."
    # Run the year's UC/ED problem
    long_econ_results, dispatch_results = Dispatch.execute_dispatch_economic_projection(
        CL_args,
        db,
        settings,
        1,   # forecast period is always 1 for current-year runs
        total_demand,
        unit_specs,
        system_portfolio_dict;
        run_mode="current",
        downselection_mode="exact"
    )

    if settings["simulation"]["file_logging_level"] > 0
        save_intermediate_outputs(settings, CL_args, long_econ_results, dispatch_results)
    end

    # Adjust formatting and save to the database
    Dispatch.finalize_annual_dispatch_results(db, CL_args["current_pd"], long_econ_results, dispatch_results)

    @info "Done!"
end

run_true_annual_dispatch()   
