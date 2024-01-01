#########################################################################
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

module ABCEfunctions

using ArgParse, CPLEX,
    Requires, SQLite, DataFrames, CSV, JuMP, GLPK, Cbc, Logging, Tables, HiGHS, Statistics

# Use CPLEX if available
#function __init__()
#    @require CPLEX = "a076750e-1247-5638-91d2-ce28b192dca0" @eval using CPLEX
#end

include("./dispatch.jl")
using .Dispatch
include("./C2N_projects.jl")
using .C2N


#####
# Agent turn setup functions
#####
function get_CL_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--settings_file"
        help = "absolute path to the settings file"
        required = false
        default = joinpath(pwd(), "settings.yml")

        "--inputs_path"
        help = "relative path to the input files (aside from settings.yml)"
        required = false
        default = joinpath(pwd(), "inputs")

        "--agent_id"
        help = "unique ID of the agent"
        required = true
        arg_type = Int

        "--current_pd"
        help = "current ABCE time period"
        required = true
        arg_type = Int

        "--verbosity"
        help = "level of output logged to the console"
        required = false
        arg_type = Int
        range_tester = x -> x in [0, 1, 2, 3]
        default = 1

        "--abce_abs_path"
        help = "absolute path to the top-level ABCE directory"
        required = true
        arg_type = String

    end

    return parse_args(s)

end


function set_verbosity(vlevel)
    # The default logging level is Info (level 0)
    lvl = 0

    if vlevel == 0
        # Only show Logging messages of severity Error (level 2000) and above
        lvl = 2000
    elseif vlevel == 1
        # Only show Logging messages of severity Warning (level 1000) and above
        lvl = 1000
    elseif vlevel == 3
        # Show Logging messages of severity Debug (level -1000) and above
        lvl = -1000
    end

    # Initialize the Logger within the Julia scope
    global_logger(ConsoleLogger(lvl))
end


function set_up_local_paths(settings, abce_abs_path)
    settings["file_paths"]["ABCE_abs_path"] = abce_abs_path
    if settings["simulation"]["annual_dispatch_engine"] == "ALEAF"
        try
            settings["file_paths"]["ALEAF_abs_path"] = ENV["ALEAF_DIR"]
        catch LoadError
            @error string(
                "The environment variable ALEAF_abs_path does not ",
                "appear to be set. Please make sure it points to ",
                "the correct directory.",
            )
        end
    else
        settings["file_paths"]["ALEAF_abs_path"] = "NULL_PATH"
    end

    return settings

end


function validate_project_data(db, settings, unit_specs, C2N_specs)
    for unit_type in unit_specs[!, :unit_type]
        if occursin("C2N", unit_type)
            unit_type_data =
                filter(:unit_type => x -> x == unit_type, unit_specs)[1, :]

            # Dummy data for schedule retrieval
            lag = 0
            fc_pd = 70
            current_pd = 0

            capex_tl, activity_schedule = project_C2N_capex(
                db,
                settings,
                unit_type_data,
                lag,
                fc_pd,
                current_pd,
                C2N_specs,
            )

            unit_specs[
                (unit_specs.unit_type .== unit_type),
                :construction_duration,
            ] .= size(capex_tl)[1]
        end
    end

    return unit_specs

end


function set_forecast_period(unit_specs, num_lags)
    # Compute forecast period as the maximum possible project horizon, based
    #   on the sum of maximum lead time, maximum construction duration, 
    #   and maximum unit life

    transform!(
        unit_specs,
        [:construction_duration, :unit_life] =>
            ((lead_time, unit_life) -> lead_time + unit_life) => :full_life,
    )
    max_horizon =
        convert(Int64, ceil(maximum(unit_specs[!, :full_life]) + num_lags))

    return max_horizon
end


#####
# Database interaction functions
#####
function load_db(db_file)
    try
        db = SQLite.DB(db_file)
        return db
    catch e
        @error "Couldn't load the database:"
        @error e
        exit()
    end
end


function get_agent_params(db, agent_id)
    try
        command = "SELECT * FROM agent_params WHERE agent_id = $agent_id"
        df = DBInterface.execute(db, command) |> DataFrame
    catch e
        @error "Could not get agent parameters from file:"
        @error e
        exit()
    end
end


function get_grouped_current_assets(db, pd, agent_id)
    # Get a list of all of the agent's currently-operating assets, plus their
    #   unit type and mandatory retirement date
    # Time-period convention:
    #   <status>_pd == "the first period where this asset has status <status>"
    cmd = string(
        "SELECT asset_id, unit_type, retirement_pd, C2N_reserved ",
        "FROM assets WHERE agent_id = ",
        agent_id,
        " AND completion_pd <= ",
        pd,
        " AND cancellation_pd > ",
        pd,
        " AND retirement_pd > ",
        pd,
    )

    asset_list = DBInterface.execute(db, cmd) |> DataFrame

    # Count the number of assets by type
    grouped_agent_assets = combine(
        groupby(asset_list, [:unit_type, :retirement_pd, :C2N_reserved]),
        nrow => :count,
    )

    return grouped_agent_assets
end


#####
# Preprocessing
#####
function get_portfolio_forecast(db, settings, current_pd, unit_specs; agent_id=nothing)
    # Retrieve the year-by-year projected portfolio for the system or
    #   current agent

    # If an agent id has been supplied, use "agent mode"; otherwise,
    #   "system mode"
    mode = "system"
    sql_filter = ""

    if agent_id != nothing
        mode = "agent"
        sql_filter = " AND agent_id = $agent_id "
    end

    @debug "Retrieving $mode portfolio forecast..."
    portfolios_dict = Dict()

    end_year = (current_pd + convert(Int64, settings["dispatch"]["num_dispatch_years"]) - 1)

    # Retrieve only the necessary unit_specs columns
    brief_unit_specs = unit_specs[!, [:unit_type, :capacity, :capacity_factor]]

    for y = current_pd:end_year
        # Retrieve annual portfolios
        sql_command = string(
            "SELECT unit_type, COUNT(unit_type) FROM assets ",
            "WHERE completion_pd <= $y AND retirement_PD > $y ",
            "AND cancellation_pd > $y ",
            sql_filter,
            "GROUP BY unit_type"
       )

        portfolios_dict[y] = DBInterface.execute(db, sql_command) |> DataFrame

        # Tidy up column names
        rename!(portfolios_dict[y], Symbol("COUNT(unit_type)") => :num_units)

        # Add a column to indicate real vs fictitious units
        portfolios_dict[y][!, :real] = ones(size(portfolios_dict[y])[1])

        # Join in the unit specs data for existing units
        portfolios_dict[y] = innerjoin(portfolios_dict[y], brief_unit_specs, on = :unit_type)

        # Compute total capacity
        transform!(portfolios_dict[y], [:num_units, :capacity] => ((num_units, cap) -> num_units .* cap) => :total_capacity)
    end

    return portfolios_dict

end


function extrapolate_demand(visible_demand, db, pd, fc_pd, settings)
    mode = settings["demand"]["demand_projection_mode"]
    if mode == "flat"
        demand = project_demand_flat(visible_demand, fc_pd)
    elseif mode == "exp_termrate"
        term_demand_gr = settings["demand"]["terminal_demand_growth_rate"]
        demand =
            project_demand_exp_termrate(visible_demand, fc_pd, term_demand_gr)
    elseif mode == "exp_fitted"
        demand =
            project_demand_exp_fitted(visible_demand, db, pd, fc_pd, settings)
    else
        @error string(
            "The specified demand extrapolation mode, $mode, is not ",
            "valid. Use 'flat', 'exp_termrate', or 'exp_fitted'.",
        )
        exit()
    end

    return demand
end


function project_demand_flat(visible_demand, fc_pd)
    known_pd = size(visible_demand)[1]
    demand = DataFrame(demand = zeros(Float64, convert(Int64, fc_pd)))
    demand[1:known_pd, :total_demand] .= visible_demand[!, :total_demand]
    demand[(known_pd + 1):fc_pd, :total_demand] .=
        demand[known_pd, :total_demand]
    return demand
end


function project_demand_exp_termrate(visible_demand, fc_pd, term_demand_gr)
    total_demand = deepcopy(visible_demand)
    prev_pd = visible_demand[size(visible_demand)[1], :period]
    prev_real_demand = last(visible_demand[!, :real_demand])
    for i = (prev_pd + 1):fc_pd
        next_real_demand =
            (prev_real_demand * (1 + term_demand_gr)^(i - prev_pd - 1))
        push!(total_demand, [i, next_real_demand])
    end
    return total_demand
end


function project_demand_exp_fitted(visible_demand, db, pd, fc_pd, settings)
    # Retrieve historical data, if available
    demand_projection_window = settings["demand"]["demand_projection_window"]
    demand_history_start =
        max(0, pd - settings["demand"]["demand_visibility_horizon"])
    start_and_end = (demand_history_start, pd)
    demand_history =
        DBInterface.execute(
            db,
            string(
                "SELECT demand FROM demand ",
                "WHERE period >= ? AND period < ? ",
                "ORDER BY period ASC",
            ),
            start_and_end,
        ) |> DataFrame

    # Create suitable arrays for x (including intercept) and y
    num_obs = size(visible_demand)[1] + size(demand_history)[1]
    x = hcat(ones(num_obs), [i for i = 0:(num_obs - 1)])
    if size(demand_history)[1] == 0
        y = visible_demand[:, :demand]
    else
        y = vcat(demand_history[:, "demand"], visible_demand[:, :demand])
    end

    # Take the log of y to make the x-y relationship linear
    y_log = log.(y)

    # Compute the estimated regression coefficient between x and log(y)
    beta = inv(transpose(x) * x) * transpose(x) * y_log

    # Compute the estimated residual (for display purposes only)
    error_est = transpose(transpose(beta) * transpose(x)) - y_log

    # If there isn't enough demand history to fulfill the desired projection
    #   window, take a weighted average of the computed beta with the known
    #   historical beta
    beta[2] = (
        (
            beta[2] * size(y)[1] +
            settings["demand"]["historical_demand_growth_rate"] *
            (demand_projection_window - size(y)[1])
        ) / demand_projection_window
    )

    # Project future demand
    proj_horiz = fc_pd - size(visible_demand)[1]
    x_proj = hcat(
        ones(proj_horiz),
        [
            i for i =
                (size(visible_demand)[1] + size(demand_history)[1] + 1):(fc_pd + size(
                    demand_history,
                )[1])
        ],
    )
    y_log_proj = x_proj[:, 1] .* beta[1] + x_proj[:, 2] .* beta[2]
    y_proj = exp.(x_proj[:, 1] .* beta[1] + x_proj[:, 2] .* beta[2])

    # Concatenate input and projected demand data into a single DataFrame
    demand = DataFrame(demand = vcat(visible_demand[!, :demand], y_proj))
    return demand
end


function get_demand_forecast(db, pd, fc_pd, settings)
    # Retrieve 
    vals = (pd, pd + settings["demand"]["demand_visibility_horizon"])
    visible_demand =
        DBInterface.execute(
            db,
            string(
                "SELECT period, demand FROM demand ",
                "WHERE period >= ? AND period < ?",
            ),
            vals,
        ) |> DataFrame

    rename!(visible_demand, :demand => :real_demand)

    demand_forecast =
        extrapolate_demand(visible_demand, db, pd, fc_pd, settings)
    prm =
        DBInterface.execute(
            db,
            "SELECT * FROM model_params WHERE parameter = 'PRM'",
        ) |> DataFrame

    transform!(
        demand_forecast,
        :real_demand =>
            (
                (dem) -> (
                    dem .* (1 + prm[1, :value]) .+
                    settings["system"]["peak_initial_reserves"]
                )
            ) => :total_demand,
    )

    return demand_forecast
end


function fill_portfolios_missing_units(current_pd, system_portfolios, unit_specs)
    # Ensure that at least 1 unit of every available type in unit_specs is
    #   represented in every year of the system portfolio, by adding 1 instance
    #   of each missing unit type.
    # Units are only added starting at the earliest year in which construction
    #   of such a unit could have been completed.

    # Retrieve only the necessary unit_specs columns
    brief_unit_specs = unit_specs[!, [:unit_type, :capacity, :capacity_factor]]

    for y = minimum(keys(system_portfolios)):maximum(keys(system_portfolios))
        for unit_type_specs in eachrow(unit_specs)
            if !in(unit_type_specs.unit_type, system_portfolios[y][!, :unit_type])
                if y - current_pd >= unit_type_specs.construction_duration
                    # Target dataframe columns:
                    # unit_type, num_units, real, capacity, capacity_factor, total_capacity
                    push!(system_portfolios[y], (unit_type_specs.unit_type, 1, 0, unit_type_specs.capacity, unit_type_specs.capacity_factor, unit_type_specs.capacity))
                end
            end
        end
    end

    return system_portfolios
end


function forecast_balance_of_market_investment(adj_system_portfolios, agent_portfolios, agent_params, current_pd, settings, demand_forecast)
    end_year = (current_pd + convert(Int64, settings["dispatch"]["num_dispatch_years"]) - 1)

    for y = current_pd:end_year
        apf = select(agent_portfolios[y], [:unit_type, :total_capacity])
        rename!(apf, :total_capacity => :agent_total_capacity)

        adj_system_portfolios[y] = coalesce.(outerjoin(adj_system_portfolios[y], apf, on = :unit_type), 0)

        transform!(adj_system_portfolios[y], [:total_capacity, :capacity_factor] => ((cap, cf) -> cap .* cf) => :total_derated_capacity)
        transform!(adj_system_portfolios[y], [:total_capacity, :agent_total_capacity] => ((total_cap, agent_cap) -> (total_cap - agent_cap) ./ total_cap) => :my)

        d_y = filter(:period => x -> x == y, demand_forecast)[1, :total_demand]
        c_y = sum(adj_system_portfolios[y][!, :total_derated_capacity])

        d_0 = filter(:period => x -> x == current_pd, demand_forecast)[1, :total_demand]
        c_0 = sum(adj_system_portfolios[current_pd][!, :total_derated_capacity])


        if c_y / d_y < c_0 / d_0
            transform!(adj_system_portfolios[y], [:total_derated_capacity, :my, :capacity_factor, :real] => ((c_iy, my, cf, real) -> c_iy .+ real .* (my .* (c_iy ./ c_y) .* (d_y .* c_0 ./ d_0 .- c_y))) => :total_esc_der_capacity)
        else
            adj_system_portfolios[y][!, :total_esc_der_capacity] .= adj_system_portfolios[y][!, :total_derated_capacity]
        end

        transform!(adj_system_portfolios[y], [:total_esc_der_capacity, :capacity_factor] => ((c_iy, cf) -> c_iy ./ cf) => :total_esc_capacity)

        transform!(adj_system_portfolios[y], [:total_esc_capacity, :capacity] => ((total_esc_cap, cap) -> ceil.(total_esc_cap ./cap)) => :esc_num_units)
    end

    return adj_system_portfolios

end


function get_next_asset_id(db)
    tables_to_check = [
        "assets",
        "WIP_projects",
        "asset_updates",
        "WIP_updates",
        "depreciation_projections",
    ]

    id_vals = Vector{Any}()

    for table in tables_to_check
        id_val =
            DBInterface.execute(db, "SELECT MAX(asset_id) FROM $table") |>
            DataFrame
        id_val = id_val[1, Symbol("MAX(asset_id)")]
        push!(id_vals, id_val)
    end

    # Convert id_vals into a skipmissing object
    id_vals = skipmissing(id_vals)

    # Return the next available asset ID (one greater than the current 
    #   largest ID)
    next_id = maximum(id_vals) + 1

    return next_id
end


function get_unit_specs(db)
    # Retrieve the table of unit specifications from the DB
    unit_specs =
        DBInterface.execute(db, "SELECT * FROM unit_specs") |> DataFrame

    # Convert the Int-type columns to Int64
    for col in names(unit_specs)
        if eltype(unit_specs[!, col]) <: Integer
            unit_specs[!, col] = convert.(Int64, unit_specs[:, col])
        end
    end

    return unit_specs
end


#####
# NPV functions
#####
function set_up_project_alternatives(
    settings,
    unit_specs,
    asset_counts,
    fc_pd,
    agent_id,
    agent_params,
    db,
    current_pd,
    C2N_specs,
    dispatch_results,
    verbosity,
)
    PA_summaries = create_PA_summaries(settings, unit_specs, asset_counts)
    PA_fs_dict = Dict()

    # Create the 
    PA_subprojects = Dict()

    for PA in eachrow(PA_summaries)
        # Create a vector of subprojects for this project alternative
        PA_subprojects[PA.uid] = create_PA_subprojects(
            settings,
            db,
            unit_specs,
            PA,
            fc_pd,
            current_pd,
            C2N_specs,
            agent_id,
            agent_params,
            dispatch_results,
        )

        # Create an aggregated financial statement for this project alternative
        #   based on its subprojects
        PA_fs_dict[PA.uid] = create_PA_aggregated_fs(PA_subprojects[PA.uid])

        # Save raw project fs results, if verbose outputs are enabled
        if verbosity > 2
            CSV.write(
                joinpath(
                    "tmp",
                    string(
                        PA["unit_type"],
                        "_",
                        PA["project_type"],
                        "_",
                        PA["lag"],
                        "_basepd_",
                        current_pd,
                        ".csv"
                    ),
                ),
                PA_fs_dict[PA.uid]
            )
        end

        # Compute the project alternative's overall NPV based on its
        #   subprojects' financial statements
        PA.NPV = compute_PA_NPV(PA_fs_dict[PA.uid])
    end

    return PA_summaries, PA_fs_dict

end

function create_PA_summaries(settings, unit_specs, asset_counts)
    PA_summaries = DataFrame(
        unit_type = String[],
        project_type = String[],
        lag = Int64[],
        ret_pd = Union{Int64,Nothing}[],
        uid = Int64[],
        NPV = Float64[],
        allowed = Bool[],
    )

    # First project ID is 1
    uid = 1

    # Always initialize NPV as zero
    NPV = 0.0

    # Set up new unique IDs for project alternatives
    # New construction
    project_type = "new_xtr"
    for unit_type in unit_specs[!, :unit_type]
        # Assume alternative is alterable unless set up otherwise
        allowed = true
        if !(unit_type in settings["scenario"]["allowed_xtr_types"])
            allowed = false
        end

        for lag = 0:settings["agent_opt"]["num_future_periods_considered"]
            push!(
                PA_summaries,
                [unit_type project_type lag nothing uid NPV allowed],
            )
            uid += 1
        end
    end

    # Retirement
    project_type = "retirement"
    allowed = true   # by default, retirements are always allowed
    for unit_type in unique(asset_counts[!, :unit_type])
        ret_pds = filter(
            [:unit_type, :C2N_reserved] =>
                ((x, reserved) -> (x == unit_type) && (reserved == 0)),
            asset_counts,
        )[
            !,
            :retirement_pd,
        ]

        for ret_pd in ret_pds
            for lag = 0:settings["agent_opt"]["num_future_periods_considered"]
                row = [unit_type project_type lag ret_pd uid NPV allowed]
                push!(PA_summaries, row)
                uid += 1
            end
        end
    end

    # Filter out all disallowed project alternatives
    PA_summaries = deepcopy(filter(:allowed => a -> a == true, PA_summaries))

    return PA_summaries

end


function create_PA_subprojects(
    settings,
    db,
    unit_specs,
    PA,
    fc_pd,
    current_pd,
    C2N_specs,
    agent_id,
    agent_params,
    dispatch_results,
)
    # Initialize the metadata and empty financial statements for all
    #   subprojects for this project alternative
    subprojects = initialize_subprojects(unit_specs, PA, fc_pd)

    # Retrieve historical dispatch results data
    historical_dispatch_results = average_historical_dispatch_results(settings, db)

    for subproject in subprojects
        # Retrieve unit type data for convenience
        unit_type_data =
            filter(:unit_type => x -> x == subproject["unit_type"], unit_specs)[1, :]

        # Forecast financial results for each subproject within this project
        #   alternative
        subproject["financial_statement"] = forecast_subproject_financials(
            settings,
            db,
            unit_type_data,
            subproject,
            fc_pd,
            current_pd,
            C2N_specs,
            agent_id,
            agent_params,
            dispatch_results,
            historical_dispatch_results,
        )
    end

    return subprojects
end


function initialize_subprojects(unit_specs, PA, fc_pd)
    # Create a vector to store the subprojects for this project
    #   alternative
    subprojects = []

    # Create the principal subproject and add it to the vector
    subproject = Dict(
        "unit_type" => PA.unit_type,
        "project_type" => PA.project_type,
        "lag" => PA.lag,
        "ret_pd" => PA.ret_pd,
    )
    push!(subprojects, subproject)

    # If the current project alternative is a new coal-to-nuclear
    #   conversion project, add a second subproject to represent the coal
    #   unit retirement
    if occursin("C2N", PA.unit_type)
        # Retrieve coal retirement scheduling information from the 
        #   unit_specs dataframe
        unit_type_data =
            filter(:unit_type => x -> x == PA.unit_type, unit_specs)[1, :]
        construction_duration = unit_type_data["construction_duration"]

        # Set the coal retirement period, depending on whether the project is
        #   greenfield or brownfield
        if occursin("C2N0", PA.unit_type)
            lag = PA.lag + construction_duration - 1
            ret_pd = PA.lag + 12
        else
            lag = PA.lag + unit_type_data["cpp_ret_lead"]
            ret_pd = PA.lag + 12
        end

        coal_retirement = Dict(
            "unit_type" => "coal",
            "project_type" => "retirement",
            "lag" => lag,
            "ret_pd" => ret_pd,
        )

        # Add the subproject to the subprojects vector
        push!(subprojects, coal_retirement)
    end

    return subprojects
end


function forecast_subproject_financials(
    settings,
    db,
    unit_type_data,
    subproject,
    fc_pd,
    current_pd,
    C2N_specs,
    agent_id,
    agent_params,
    dispatch_results,
    historical_dispatch_results,
)
    # Create a blank DataFrame for the subproject's financial statement
    subproject_fs = DataFrame(
        year = 1:fc_pd,
        weight = zeros(fc_pd),
        capex = zeros(fc_pd),
        remaining_debt_principal = zeros(fc_pd),
        debt_payment = zeros(fc_pd),
        interest_payment = zeros(fc_pd),
        depreciation = zeros(fc_pd),
        generation = zeros(fc_pd),
        revenue = zeros(fc_pd),
        VOM = zeros(fc_pd),
        fuel_cost = zeros(fc_pd),
        FOM = zeros(fc_pd),
        policy_adj = zeros(fc_pd),
        tax_credits = zeros(fc_pd),
    )

    # Compute the series of DCF weights (discounting by cost of equity, as DCF
    #   is calculated on post-interest net cash flow)
    coe = agent_params[1, :cost_of_equity]

    transform!(
        subproject_fs,
        [:year] => ((yr) -> 1 ./ (1 .+ coe) .^ (yr .- 1)) => :weight,
    )

    if subproject["project_type"] == "new_xtr"
        # Add capital expenditures projection
        subproject_fs = forecast_capex(
            settings,
            db,
            unit_type_data,
            subproject,
            current_pd,
            C2N_specs,
            deepcopy(subproject_fs),
        )

        # Set up the time-series of outstanding debt based on this
        #   construction expenditure profile
        subproject_fs = forecast_construction_debt_principal(
            unit_type_data,
            subproject,
            agent_params,
            deepcopy(subproject_fs),
        )

        # Forecast the debt repayment schedule
        subproject_fs = forecast_debt_schedule(
            agent_params,
            unit_type_data,
            subproject,
            deepcopy(subproject_fs),
        )

        # Forecast the depreciation schedule
        subproject_fs = forecast_depreciation(settings, deepcopy(subproject_fs))

    end

    # Forecast all operational results for this subproject: generation,
    #   revenue, all cost types, policy adjustments
    subproject_fs = forecast_subproject_operations(
        settings,
        current_pd,
        subproject,
        unit_type_data,
        dispatch_results,
        historical_dispatch_results,
        deepcopy(subproject_fs),
    )

    # Forecast all post-facto policy adjustments
    subproject_fs = forecast_subproject_pf_policy_adj(
        settings,
        current_pd,
        subproject,
        unit_type_data,
        deepcopy(subproject_fs)
    )

    subproject_fs =
        compute_accounting_line_items(db, deepcopy(subproject_fs), agent_id, agent_params, mode="subproject")

    return subproject_fs
end


function forecast_capex(
    settings,
    db,
    unit_type_data,
    subproject,
    current_pd,
    C2N_specs,
    fs_copy,
)
    capex_per_pd = (
        unit_type_data[:overnight_capital_cost] *
        unit_type_data[:capacity] *
        settings["constants"]["MW2kW"] /
        unit_type_data[:construction_duration]
    )
    capex_timeline =
        ones(convert(Int64, unit_type_data[:construction_duration])) * capex_per_pd

    head_zeros = zeros(subproject["lag"])
    tail_zeros =
        zeros(size(fs_copy)[1] - subproject["lag"] - size(capex_timeline)[1])

    fs_copy[!, :capex] = vcat(head_zeros, capex_timeline, tail_zeros)

    return fs_copy
end


function project_C2N_capex(
    db,
    settings,
    unit_type_data,
    lag,
    fc_pd,
    current_pd,
    C2N_specs,
)
    assumption = settings["simulation"]["C2N_assumption"]
    # Set the project parameters
    if occursin("C2N0", unit_type_data[:unit_type])
        conversion_type = "greenfield"
        if occursin("PWR", unit_type_data[:unit_type])
            rxtr_type = "PWR"
        elseif occursin("HTGR", unit_type_data[:unit_type])
            rxtr_type = "HTGR"
        else
            rxtr_type = "SFR"
        end

        data = deepcopy(C2N_specs[conversion_type][rxtr_type])

    else
        if occursin("C2N1", unit_type_data[:unit_type])
            conversion_type = "electrical"
            rxtr_type = "PWR"
        elseif occursin("C2N2", unit_type_data[:unit_type])
            conversion_type = "steam_noTES"
            rxtr_type = "HTGR"
        else
            conversion_type = "steam_TES"
            rxtr_type = "SFR"
        end
        data = deepcopy(C2N_specs[conversion_type][assumption])
    end

    # Scale cost component data to match unit type specification data
    for key in keys(data)
        data[key]["cost_rem"] = data[key]["cost_rem"] * unit_type_data[:overnight_capital_cost]
    end

    # Develop the C2N capex profile
    capex_tl, activity_schedule = C2N.create_C2N_capex_timeline(
        db,
        conversion_type,
        rxtr_type,
        current_pd,
        lag,
        fc_pd,
        data,
        unit_type_data,
    )

    return capex_tl, activity_schedule
end


function get_capex_end(fs_copy)
    # Find the end of the capex accumulation period
    capex_end = size(fs_copy)[1]
    for i = 1:size(fs_copy)[1]
        if (fs_copy[i, "capex"]) != 0 && (fs_copy[i + 1, "capex"] == 0)
            capex_end = i
        end
    end

    return capex_end
end


function forecast_construction_debt_principal(
    unit_type_data,
    subproject,
    agent_params,
    fs_copy,
)
    # Find the end of the capex accumulation period
    capex_end = get_capex_end(fs_copy)

    # Compute a rolling sum of accumulated debt principal
    debt_timeline = zeros(capex_end)
    for i = 1:capex_end
        debt_timeline[i] =
            sum(fs_copy[1:i, :capex]) * agent_params[1, :debt_fraction]
    end

    # Pad the rest of the series with zeros (will be populated later)
    tail_zeros = zeros(size(fs_copy)[1] - capex_end)

    # Add an appropriate column to the dataframe
    fs_copy[!, :remaining_debt_principal] = vcat(debt_timeline, tail_zeros)

    return fs_copy
end


function forecast_debt_schedule(
    agent_params,
    unit_type_data,
    subproject,
    fs_copy,
)
    # Retrieve maximum debt level as the principal
    principal = maximum(fs_copy[!, :remaining_debt_principal])

    # Find the end of the capex accumulation
    capex_end = get_capex_end(fs_copy)

    # Compute the constant recurring payment into the sinking fund
    pmt =
        agent_params[1, :cost_of_debt] * principal /
        (1 - (1 + agent_params[1, :cost_of_debt])^(-unit_type_data[:unit_life]))

    # The end of the schedule is the sooner of the end of the financial
    #   statement or the end of the unit's life
    fin_end = convert(
        Int64,
        min(
            size(fs_copy)[1],
            subproject["lag"] +
            unit_type_data[:construction_duration] +
            unit_type_data[:unit_life],
        ),
    )

    # Seed post-construction financing series
    fs_copy[capex_end + 1, :remaining_debt_principal] =
        fs_copy[capex_end, :remaining_debt_principal]

    for i = (capex_end + 1):fin_end
        fs_copy[i, :debt_payment] = pmt
        fs_copy[i, :interest_payment] =
            fs_copy[i, :remaining_debt_principal] *
            agent_params[1, :cost_of_debt]
        if i != fin_end
            fs_copy[i + 1, :remaining_debt_principal] =
                fs_copy[i, :remaining_debt_principal] -
                (fs_copy[i, :debt_payment] - fs_copy[i, :interest_payment])
        end
    end

    return fs_copy
end


function forecast_depreciation(settings, fs_copy)
    capex_end = get_capex_end(fs_copy)
    total_capex = sum(fs_copy[!, :capex])

    dep_end = min(
        size(fs_copy)[1],
        capex_end + settings["financing"]["depreciation_horizon"],
    )

    for i = (capex_end + 1):dep_end
        fs_copy[i, :depreciation] =
            total_capex / settings["financing"]["depreciation_horizon"]
    end

    return fs_copy
end


function forecast_subproject_operations(
    settings,
    current_pd,
    subproject,
    unit_type_data,
    dispatch_results,
    historical_dispatch_results,
    fs_copy,
)
    mode = subproject["project_type"]
    hist_wt = settings["dispatch"]["hist_wt"]
    data_to_get =
        ["generation", "revenue", "VOM", "fuel_cost", "FOM", "policy_adj", "tax_credits"]

    # Get historical dispatch results for this unit type
    historical_dispatch_results = filter(
        :unit_type => unit_type -> unit_type == subproject["unit_type"],
        historical_dispatch_results,
    )

    # Get projected dispatch results for this unit type
    ABCE_dispatch_results = filter(
        :unit_type => unit_type -> unit_type == subproject["unit_type"],
        dispatch_results,
    )

    # Set up timeline start/end and value sign based on project type
    if subproject["project_type"] == "new_xtr"
        # Record marginal additional generation
        # get_capex_end gives the relative index, not the absolute year
        series_start = get_capex_end(fs_copy) + 1
        series_end = convert(
            Int64,
            round(
                min(
                    size(fs_copy)[1],
                    series_start + unit_type_data[:unit_life],
                ),
                digits = 3,
            ),
        )
        sign = 1
    elseif subproject["project_type"] == "retirement"
        # Record foregone generation as the negative of the projection
        series_start = convert(Int64, subproject["lag"] + 1)
        series_end = convert(Int64, min(size(fs_copy)[1], subproject["ret_pd"]))
        sign = -1
    end

    # Update the operations results data in the subproject's financial statement
    # series_start and series_end are relative indices, not absolute years
    #   --they need to be converted with current_pd
    for i = series_start:series_end
        # Convert relative period to absolute
        yr = current_pd + i - 1

        for data_type in data_to_get
            # Get the corresponding data value for the year yr from the ABCE
            #   dispatch projection
            # If year yr was actually simulated by dispatch.jl, retrieve
            #   the corresponding result
            if yr in ABCE_dispatch_results[!, :y]
                ABCE_data_value = filter(
                    [:y, :dispatch_result] =>
                        (y, disp_res) ->
                            (y == yr) && (disp_res == data_type),
                    ABCE_dispatch_results,
                )
                ABCE_data_value = ABCE_data_value[1, :qty]

            # If year yr was not simulated by dispatch.jl, use the dispatch
            #   results from the last simulated year
            else
                try
                    last_dispatch_year = maximum(
                        filter(
                            :unit_type => unit_type -> unit_type == subproject["unit_type"],
                            ABCE_dispatch_results
                        )[!, :y]
                    )

                    ABCE_data_value = filter(
                        [:y, :dispatch_result] =>
                            (y, disp_res) ->
                                (y == last_dispatch_year) &&
                                    (disp_res == data_type),
                        ABCE_dispatch_results,
                    )
                    ABCE_data_value = ABCE_data_value[1, :qty]
                catch
                    ABCE_data_value = 0.0
                end
            end

            # Get the corresponding cumulative historical estimate from the 
            #   A-LEAF aggregated dispatch histories
            if size(historical_dispatch_results)[1] > 0
                historical_data_value = sum(
                    historical_dispatch_results[
                        !,
                        Symbol(string("wtd_", data_type)),
                    ],
                )
                hist_wt = settings["dispatch"]["hist_wt"]
            else
                # If this unit type does not appear in the A-LEAF historical 
                #   results, rely only on the ABCE dispatch forecast
                historical_data_value = 0
                hist_wt = 0
            end

            # Save the signed value into the financial statement
            fs_copy[i, Symbol(data_type)] =
                sign * ABCE_data_value * (1 - hist_wt) +
                sign * historical_data_value * (hist_wt)
        end

        # Zero out production tax credits if the subproject would not actually be
        #   eligible for them during this year
        if "policies" in keys(settings["scenario"])
            if "PTC" in keys(settings["scenario"]["policies"])
                if subproject["unit_type"] in settings["scenario"]["policies"]["PTC"]["eligibility"]["unit_type"]
                    PTC_start = settings["scenario"]["policies"]["PTC"]["start_year"]
                    PTC_end = settings["scenario"]["policies"]["PTC"]["end_year"]
                    PTC_duration = settings["scenario"]["policies"]["PTC"]["duration"]
                    unit_xtr_duration = unit_type_data[:construction_duration]

                    # Check if this unit is eligible for the PTC by its
                    #   construction start date
                    if (PTC_start > current_pd + subproject["lag"]) || (PTC_end < current_pd + subproject["lag"]) || (yr >= current_pd + subproject["lag"] + unit_xtr_duration + PTC_duration)
                        fs_copy[i, "tax_credits"] = 0
                    end
                end
            end
        end

    end

    return fs_copy

end


function forecast_subproject_pf_policy_adj(settings, current_pd, subproject, unit_type_data, fs_copy)
    try
        ITC_data = settings["scenario"]["policies"]["ITC"]
        ITC_qty = ITC_data["qty"]

        if ((ITC_data["enabled"]) && (subproject["unit_type"] in ITC_data["eligibility"]["unit_type"]) && (ITC["start_year"] <= current_pd + subproject["lag"]) &&(ITC["end_year"] >= current_pd + subproject["lag"]))
            capex_end = get_capex_end(fs_copy)
            total_capex = sum(fs_copy[!, :capex])

            if (capex_end != nothing) && (capex_end <= size(fs_copy)[1] - 1)
                fs_copy[capex_end+1, :tax_credits] = ITC_qty * total_capex
            end
        end
    catch
        # do nothing
    end

    return fs_copy
end


function create_PA_aggregated_fs(subprojects)
    # Initialize the project alternative's financial statement with the 
    #   data from the first subproject
    PA_aggregated_fs = deepcopy(subprojects[1]["financial_statement"])

    # If other subprojects are present, add them in
    cols_to_collect = names(PA_aggregated_fs, Not([:year, :weight]))
    for subproject in deleteat!(subprojects, 1)
        for col in cols_to_collect
            PA_aggregated_fs[!, col] .+=
                subproject["financial_statement"][!, col]
        end
    end

    return PA_aggregated_fs
end


function compute_PA_NPV(fs_copy)
    NPV = 0

    transform!(
        fs_copy,
        [:FCF, :weight] => ((fcf, wt) -> fcf .* wt) => :wtd_net_cash,
    )
    NPV += sum(fs_copy[!, :wtd_net_cash])

    return NPV
end


function average_historical_dispatch_results(settings, db)
    historical_dispatch_results =
        DBInterface.execute(db, "SELECT * FROM annual_dispatch_unit_summary") |>
        DataFrame

    if size(historical_dispatch_results)[1] != 0
        # Get a list of all unique years in the dataframe
        years = unique(historical_dispatch_results, :period)[!, :period]

        hist_decay = settings["dispatch"]["hist_decay"]

        # Compute scaling factor to normalize to 1
        c = 1 / sum(1 / (1 + hist_decay)^k for k in years)

        # Add a column with the diminishing historical weighting factor
        transform!(
            historical_dispatch_results,
            [:period] =>
                ((pd) -> c ./ (1 + hist_decay) .^ pd) => :hist_wt_coeffs,
        )

        # Weight the data columns
        data_to_weight =
            ["generation", "revenue", "VOM", "fuel_cost", "FOM", "policy_adj", "tax_credits"]

        for data_type in data_to_weight
            transform!(
                historical_dispatch_results,
                [Symbol(data_type), :hist_wt_coeffs] =>
                    ((rev, wt) -> rev .* wt) =>
                        Symbol(string("wtd_", data_type)),
            )
        end
    end

    return historical_dispatch_results

end


function compute_moodys_score(icr, fcf_debt, re_debt)
    icr_wt = 0.1
    fcf_debt_wt = 0.2
    re_debt_wt = 0.1

    icr_score = 0
    if icr >= 18
        icr_score = 1
    elseif icr >= 13
        icr_score = 3
    elseif icr >= 8
        icr_score = 6
    elseif icr >= 4.2
        icr_score = 9
    elseif icr >= 2.8
        icr_score = 12
    elseif icr >= 1
        icr_score = 15
    else
        icr_score = 18
    end

    fcf_debt_score = 0
    if fcf_debt >= 0.9
        fcf_debt_score = 1
    elseif fcf_debt >= 0.6
        fcf_debt_score = 3
    elseif fcf_debt >= 0.35
        fcf_debt_score = 6
    elseif fcf_debt >= 0.2
        fcf_debt_score = 9
    elseif fcf_debt >= 0.12
        fcf_debt_score = 12
    elseif fcf_debt >= 0.05
        fcf_debt_score = 15
    else
        fcf_debt_score = 18
    end

    re_debt_score = 0
    if re_debt >= 0.6
        re_debt_score = 1
    elseif re_debt >= 0.45
        re_debt_score = 3
    elseif re_debt >= 0.25
        re_debt_score = 6
    elseif re_debt >= 0.15
        re_debt_score = 9
    elseif re_debt >= 0.08
        re_debt_score = 12
    elseif re_debt >= 0.03
        re_debt_score = 15
    else
        re_debt_score = 18
    end

    score = (icr_wt * icr_score + fcf_debt_wt * fcf_debt_score + re_debt_wt * re_debt_score) / (icr_wt + fcf_debt_wt + re_debt_wt)

    return score
end


#####
# JuMP optimization model initialization
#####
function create_model_with_optimizer(settings)
    # Determine which solver to use, based on the settings file
    solver = lowercase(settings["simulation"]["solver"])

    if solver == "cplex"
        m = Model(CPLEX.Optimizer)
    elseif solver == "glpk"
        m = Model(GLPK.Optimizer)
    elseif solver == "cbc"
        m = Model(Cbc.Optimizer)
    elseif solver == "highs"
        m = Model(HiGHS.Optimizer)
    else
        throw(error("Solver `$solver` not supported. Try `cplex` instead."))
    end

    set_silent(m)

    return m

end


"""
    set_up_model(unit_FS_dict, ret_fs_dict, fc_pd, available_demand, NPV_results, ret_NPV_results)

Set up the JuMP optimization model, including variables, constraints, and the
objective function.

Returns:
  m (JuMP model object)
"""
function set_up_model(
    settings,
    PA_summaries,
    PA_fs_dict,
    total_demand,
    asset_counts,
    agent_params,
    unit_specs,
    current_pd,
    adj_system_portfolios,
    db,
    agent_id,
    agent_fs,
    fc_pd;
    mode="normal"
)
    # Create the model object
    @debug "Setting up model..."

    # Initialize a model with the settings-specified optimizer
    m = create_model_with_optimizer(settings)

    # If running in "ret_only" mode, change the "allowed" value for all
    #   "new_xtr" projects to 0 to forbid them
    if mode == "ret_only"
        for i=1:size(PA_summaries)[1]
            if PA_summaries[i, :project_type] == "new_xtr"
                PA_summaries[i, :allowed] = false
            end
        end
    end

    # Parameter names
    num_alternatives = size(PA_summaries)[1]
    num_time_periods = size(PA_fs_dict[PA_summaries[1, :uid]])[1]

    # Set up variables
    # Number of units of each type to build: must be Integer
    @variable(m, u[1:num_alternatives] >= 0, Int)

    # Compute expected marginal generation and effective nameplate capacity
    #   contribution per alternative type
    marg_gen = zeros(num_alternatives, num_time_periods)
    marg_eff_cap = zeros(num_alternatives, num_time_periods)
    for i = 1:size(PA_summaries)[1]
        for j = 1:num_time_periods
            # Marginal generation in kWh
            marg_gen[i, j] = (
                PA_fs_dict[PA_summaries[i, :uid]][j, :generation] /
                settings["constants"]["MW2kW"]
            )
            # Convert total anticipated marginal generation to an effective
            #   nameplate capacity and save to the appropriate entry in the
            #   marginal generation array
            unit_type_data = filter(
                :unit_type => x -> x == PA_summaries[i, :unit_type],
                unit_specs,
            )

            if PA_summaries[i, :project_type] == "new_xtr"
                op_start =
                    PA_summaries[i, :lag] +
                    unit_type_data[1, :construction_duration]
                if j > op_start
                    marg_eff_cap[i, j] = (
                        unit_type_data[1, :capacity] *
                        unit_type_data[1, :capacity_factor]
                    )
                end
            elseif PA_summaries[i, :project_type] == "retirement"
                if (j > PA_summaries[i, :lag]) &&
                   (j <= PA_summaries[i, :ret_pd])
                    marg_eff_cap[i, j] = (
                        (-1) *
                        unit_type_data[1, :capacity] *
                        unit_type_data[1, :capacity_factor]
                    )
                end
            end
        end
    end

    # Record which elements take effect immediately
    PA_summaries[!, :current] .= 0
    for i = 1:size(PA_summaries)[1]
        if PA_summaries[i, :project_type] == "retirement"
            PA_summaries[i, :current] = 1
        elseif (
            (PA_summaries[i, :project_type] == "new_xtr") &&
            (PA_summaries[i, :lag] == 0)
        )
            PA_summaries[i, :current] = 1
        end
    end

    # Prevent the agent from intentionally causing foreseeable energy shortages
    for i = 1:settings["agent_opt"]["shortage_protection_period"]
        # Get extant capacity sufficiency measures (projected demand and
        #   capacity)
        sufficiency =
            filter(:period => x -> x == current_pd + i - 1, total_demand)
        pd_total_demand = sufficiency[1, :total_demand]
        total_eff_cap = sum(adj_system_portfolios[current_pd + i][!, :total_esc_der_capacity])

        # Enforce different capacity assurance requirements based on extant
        #   reserve margin
        if (
            total_eff_cap / pd_total_demand >
            settings["agent_opt"]["cap_decrease_threshold"]
        )
            margin = settings["agent_opt"]["cap_decrease_margin"]
        else
            margin = settings["agent_opt"]["cap_maintain_margin"]
        end

        # Get total agent portfolio size for this year
        stmt = "SELECT unit_type, COUNT(unit_type) FROM assets WHERE agent_id = $agent_id AND completion_pd <= $i AND retirement_pd > $i GROUP BY unit_type"
        agent_pf = DBInterface.execute(db, stmt) |> DataFrame
        if size(agent_pf)[1] != 0
            rename!(agent_pf, Symbol("COUNT(unit_type)") => :num_units)
            agent_pf = leftjoin(agent_pf, select(unit_specs, [:unit_type, :capacity, :capacity_factor]), on=:unit_type)
            transform!(agent_pf, [:num_units, :capacity, :capacity_factor] => ((num, cap, cf) -> num .* cap .* cf) => :total_derated_cap)
            agent_year_derated_cap = sum(agent_pf[!, :total_derated_cap])

            @constraint(
                m,
                transpose(u .* PA_summaries[:, :current]) * marg_eff_cap[:, i] >=
                agent_year_derated_cap * margin
            )
        end
    end

    
    # Create arrays of expected marginal debt, interest, dividends, and FCF
    #   per unit type
    marg_debt = zeros(num_alternatives, num_time_periods)
    marg_int = zeros(num_alternatives, num_time_periods)
    marg_FCF = zeros(num_alternatives, num_time_periods)
    marg_retained_earnings = zeros(num_alternatives, num_time_periods)
    for i = 1:size(PA_summaries)[1]
        project = PA_fs_dict[PA_summaries[i, :uid]]
        for j = 1:num_time_periods
            # Retrieve the marginal value of interest
            # Scale to units of $B
            marg_debt[i, j] = project[j, :remaining_debt_principal] / 1e9
            marg_int[i, j] = project[j, :interest_payment] / 1e9
            marg_FCF[i, j] = project[j, :FCF] / 1e9
            marg_retained_earnings[i, j] = project[j, :FCF] / 1e9 * (1 - agent_params[1, :cost_of_equity])
        end
    end


    # Prevent the agent from reducing its aggregated financial metrics score
    #   below the weighted sum of the Moody's Ba rating thresholds (from the 
    #   Unregulated Power Companies ratings grid)
    # Getting some values for conciseness
    ICR_floor = settings["agent_opt"]["icr_floor"]
    CDR_floor = settings["agent_opt"]["fcf_debt_floor"]
    RCDR_floor = settings["agent_opt"]["re_debt_floor"]

    ICR_solo_floor = settings["agent_opt"]["icr_solo_floor"]
    CDR_solo_floor = settings["agent_opt"]["fcf_debt_solo_floor"]
    RCDR_solo_floor = settings["agent_opt"]["re_debt_solo_floor"]


    if mode == "normal"
        for i = 1:settings["agent_opt"]["fin_metric_horizon"]
            # Limit aggregated score
            @constraint(
                m,
                0.1 / abs(ICR_floor) * (
                    agent_fs[i, :FCF] / 1e9 + sum(u .* marg_FCF[:, i])
                    + (1 - ICR_floor) * (agent_fs[i, :interest_payment] / 1e9 + sum(u .* marg_int[:, i]))
                )

                + 0.2 / abs(CDR_floor) * (
                    (agent_fs[i, :FCF] / 1e9 + sum(u .* marg_FCF[:, i])) 
                    - CDR_floor * (agent_fs[i, :remaining_debt_principal] / 1e9 + sum(u .* marg_debt[:, i])) 
                )

                + 0.1 / abs(RCDR_floor) * (
                    (agent_fs[i, :retained_earnings] / 1e9 + sum(u .* marg_retained_earnings[:, i])) 
                    - RCDR_floor * (agent_fs[i, :remaining_debt_principal] / 1e9 + sum(u .* marg_debt[:, i]))
                )

                >= 0
            )

            # Limit individual financial metrics
            @constraint(
                m,
                agent_fs[i, :FCF] / 1e9 + sum(u .* marg_FCF[:, i])
                    + (1 - ICR_solo_floor) * (agent_fs[i, :interest_payment] / 1e9 + sum(u .* marg_int[:, i]))
                    >= 0
            )

            @constraint(
                m,
                (agent_fs[i, :FCF] / 1e9 + sum(u .* marg_FCF[:, i])) 
                    - CDR_solo_floor * (agent_fs[i, :remaining_debt_principal] / 1e9 + sum(u .* marg_debt[:, i])) 
                    >= 0
            )

            @constraint(
                m,
                (agent_fs[i, :retained_earnings] / 1e9 + sum(u .* marg_retained_earnings[:, i])) 
                    - RCDR_solo_floor * (agent_fs[i, :remaining_debt_principal] / 1e9 + sum(u .* marg_debt[:, i]))
                    >= 0
            )
        end
    end

    # Enforce the user-specified maximum number of new construction/retirement
    #   projects by type per period, and the :allowed field in PA_summaries
    exp_PA_summaries = leftjoin(PA_summaries, unit_specs, on=:unit_type)

    for i = 1:num_alternatives
        if PA_summaries[i, :project_type] == "new_xtr"
            @constraint(
                m,
                u[i] .* exp_PA_summaries[i, :capacity] .<=
                convert(Int64, exp_PA_summaries[i, :allowed]) .*
                settings["agent_opt"]["max_type_newcap_per_pd"]
            )
        elseif PA_summaries[i, :project_type] == "retirement"
            unit_type = PA_summaries[i, :unit_type]
            ret_pd = PA_summaries[i, :ret_pd]
            asset_count = filter(
                [:unit_type, :retirement_pd] =>
                    (x, y) -> x == unit_type && y == ret_pd,
                asset_counts,
            )[
                1,
                :count,
            ]
            max_retirement = (
                convert(Int64, PA_summaries[i, :allowed]) .* min(
                    asset_count,
                    settings["agent_opt"]["max_type_rets_per_pd"],
                )
            )
            @constraint(m, u[i] .<= max_retirement)
        end
    end

    # The total number of assets of each type to be retired cannot exceed
    #   the total number of assets of that type owned by the agent

    # Setup
    # Convenient variable for the number of distinct retirement alternatives,
    #   excluding assets reserved for C2N conversion
    retireable_asset_counts =
        filter([:C2N_reserved] => ((reserved) -> reserved == 0), asset_counts)

    # Shortened name for the number of lag periods to consider
    #   1 is added, as the user-set value only specifies future periods,
    #   with the "lag = 0" instance being implied
    num_lags = settings["agent_opt"]["num_future_periods_considered"] + 1

    # Create the matrix to collapse lagged options
    # This matrix has one long row for each retiring-asset category,
    #   with 1's in each element where the corresponding element of u[] is
    #   one of these units
    ret_summation_matrix =
        zeros(size(retireable_asset_counts)[1], size(PA_summaries)[1])
    for i = 1:size(retireable_asset_counts)[1]
        for j = 1:size(PA_summaries)[1]
            if (
                (PA_summaries[j, :project_type] == "retirement") &
                (
                    PA_summaries[j, :unit_type] ==
                    retireable_asset_counts[i, :unit_type]
                ) &
                (
                    PA_summaries[j, :ret_pd] ==
                    retireable_asset_counts[i, :retirement_pd]
                )
            )
                ret_summation_matrix[i, j] = 1
            end

        end
    end

    # Specify constraint: the agent cannot plan to retire more units (during
    #   all lag periods) than exist of that unit type
    for i = 1:size(retireable_asset_counts)[1]
        @constraint(
            m,
            sum(ret_summation_matrix[i, :] .* u) <=
            retireable_asset_counts[i, :count]
        )
    end

    # Set the maximum lookahead for coal retirements to be after the latest
    #   possible coal unit retirement from any project alternative
    coal_horizon = convert(Int64, ceil(round(maximum(PA_summaries[!, :lag]) + maximum(unit_specs[!, :construction_duration]) + 2, digits=3)))

    # Set up the coal retirements matrix, to ensure that the total "pool"
    #   of coal plants available for retirement is respected across the entire
    #   visible horizon
    coal_retirements = zeros(size(PA_summaries)[1], coal_horizon)
    planned_coal_units_operating = DataFrame(pd = Int64[], num_units = Int64[])
    planned_coal_retirements = DataFrame(pd = Int64[], num_units = Int64[])

    for i = 1:size(PA_summaries)[1]
        # Mark all of the direct coal retirements in the matrix
        if (PA_summaries[i, :project_type] == "retirement") &&
           (PA_summaries[i, :unit_type] == "coal")
            coal_retirements[i, PA_summaries[i, :lag] + 1] = 1
        end

        # Mark all of the C2N-forced coal retirements in the matrix
        if occursin("C2N", PA_summaries[i, :unit_type])
            unit_type_data = filter(
                :unit_type => ((ut) -> ut == PA_summaries[i, :unit_type]),
                unit_specs,
            )[1, :]

            if occursin("C2N0", PA_summaries[i, :unit_type])
                target_coal_ret_pd = convert(Int64, PA_summaries[i, :lag] + ceil(round(unit_type_data[:construction_duration], digits=3))) + 1
            else
                target_coal_ret_pd =
                    convert(Int64, PA_summaries[i, :lag] + ceil(round(unit_type_data[:cpp_ret_lead], digits=3))) + 1
            end

            if target_coal_ret_pd <= size(coal_retirements)[2]
                coal_retirements[i, target_coal_ret_pd] = 1
            end
        end

    end

    CSV.write("./coal_rets_indicators.csv", DataFrame(coal_retirements, :auto))

    for i = 1:size(coal_retirements)[2]
        y = current_pd + i - 1  # -1 converts i from julia dataframe index back to true relative year
        num_units = 0
        num_rets = 0

        # Get number of planned coal units operating during this period
        coal_ops =
            DBInterface.execute(
                db,
                "SELECT unit_type, COUNT(unit_type) FROM assets WHERE unit_type = 'coal' AND C2N_reserved = 0 AND agent_id = $agent_id AND completion_pd <= $y AND retirement_pd >= $y AND cancellation_pd >= $y GROUP BY unit_type",
            ) |> DataFrame

        coal_rets = DBInterface.execute(
            db,
            "SELECT unit_type, COUNT(unit_type) FROM assets WHERE unit_type = 'coal' AND agent_id = $agent_id AND C2N_reserved = 0 AND retirement_pd = $y GROUP BY unit_type",
        ) |> DataFrame

        if size(coal_ops)[1] > 0
            num_units = coal_ops[1, "COUNT(unit_type)"]
        end

        if size(coal_rets)[1] > 0
            num_rets = coal_rets[1, "COUNT(unit_type)"]
        end

        # Record the number of operating coal units in this period
        push!(planned_coal_units_operating, (y, num_units))
        push!(planned_coal_retirements, (y, num_rets))
    end

    # Introduce a slack variable for limiting coal retirements
    @variable(m, z[1:size(coal_retirements)[2]] >= 0, Int)

    # Constrain rolling total of coal retirements
    for i = 1:size(coal_retirements)[2]
        @constraint(
            m,
            sum(transpose(u) * coal_retirements[:, i]) <=
            planned_coal_units_operating[i, :num_units]
        )

        # The following three constraints are equivalent to:
        # O(t=0) - sum{t=0 -> i} (max(R(t), P(t)) >= 0
        #   or, more verbosely:
        # @constraint(m, max(planned_coal_retirements[i, :num_units], sum(transpose(u) * coal_retirements[:, i]
        #                <= planned_coal_units_operating[1, :num_units])
        @constraint(
            m,
            z[i] >= planned_coal_retirements[i, :num_units]
        )

        @constraint(
            m,
            z[i] >= sum(transpose(u) * coal_retirements[:, i])
        )

        @constraint(
            m,
            sum(z[j] for j=1:i) <= planned_coal_units_operating[1, :num_units]
        )
    end

    # Create the objective function 
    fin_metric_horizon = settings["agent_opt"]["fin_metric_horizon"]
    int_bound = settings["agent_opt"]["int_bound"]

    # Average credit rating over the horizon
    avg_cr = mean(agent_fs[1:fin_metric_horizon, :moodys_score])
    cr_adj = 10.5 / avg_cr

    # Set relative valuation-vs-credit metrics priority, depending on current
    #   credit grade
    profit_lamda = settings["agent_opt"]["profit_lamda"] / 1e9 * cr_adj
    credit_rating_lamda = settings["agent_opt"]["credit_rating_lamda"] / cr_adj

    @objective(
        m,
        Max,
        (
            profit_lamda * (transpose(u) * PA_summaries[!, :NPV]) +
            credit_rating_lamda * (
                sum(agent_fs[1:fin_metric_horizon, :FCF]) / 1e9 +
                sum(transpose(u) * marg_FCF[:, 1:fin_metric_horizon]) +
                sum(agent_fs[1:fin_metric_horizon, :interest_payment]) / 1e9 +
                sum(transpose(u) * marg_int[:, 1:fin_metric_horizon]) -
                (int_bound) * (
                    sum(agent_fs[1:fin_metric_horizon, :interest_payment]) / 1e9 +
                    sum(transpose(u) * marg_int[:, 1:fin_metric_horizon])
                )
            )
        )
    )


    @debug "Optimization model set up."

    return m
end





### Postprocessing
function finalize_results_dataframe(m, mode, PA_summaries)
    # Check solve status of model
    status = string(termination_status.(m))

    # If the model solved to optimality, convert the results to type Int64
    if status == "OPTIMAL"
        unit_qty = Int64.(round.(value.(m[:u])))
    else
        # If the model did not solve to optimality, the agent does nothing. Return
        #   a vector of all zeroes instead.
        unit_qty = zeros(Int64, size(PA_summaries)[1])
    end

    all_results = hcat(PA_summaries, DataFrame(units_to_execute = unit_qty))
    all_results[!, :mode] .= mode

    return all_results
end


function postprocess_agent_decisions(
    settings,
    all_results,
    unit_specs,
    db,
    current_pd,
    agent_id,
)
    # Break the results dataframe into non-C2N (processed first) and
    #   C2N (processed second)
    regular_results = filter(:unit_type => unit_type -> !occursin("C2N", unit_type), all_results)
    C2N_results = filter(:unit_type => unit_type -> occursin("C2N", unit_type), all_results)

    for i = 1:size(regular_results)[1]
        # Retrieve the individual result row for convenience
        result = regular_results[i, :]

        if result[:project_type] == "new_xtr"
            # New construction decisions are binding for this period only
            if (result[:lag] == 0) && (result[:units_to_execute] != 0)
                record_new_construction_projects(
                    settings,
                    result,
                    unit_specs,
                    db,
                    current_pd,
                    agent_id,
                )

                # If the project is a C2N project, retire the necessary number
                #   of coal units at the appropriate time
                if occursin("C2N", result[:unit_type])
                    record_asset_retirements(
                        result,
                        db,
                        agent_id,
                        unit_specs,
                        current_pd,
                        mode = "C2N_newbuild",
                    )
                end
            end
        elseif result[:project_type] == "retirement"
            # Retirements are recorded as binding, whether the retirement date is
            #   this period or in the future
            if result[:units_to_execute] != 0
                record_asset_retirements(
                    result,
                    db,
                    agent_id,
                    unit_specs,
                    current_pd,
                )
            end
        else
            @warn "I'm not sure what to do with the following project type:"
            @warn result
            @warn "This decision entry will be skipped."
        end

    end

    # Process all C2N results
    # Join in some coal retirement info from unit_specs
    C2N_results = leftjoin(C2N_results, unit_specs[:, [:unit_type, :cpp_ret_lead, :num_cpp_rets]], on = :unit_type)

    # Sort from earliest CPP retirement date to latest
    sort!(C2N_results, :cpp_ret_lead)

    # Record updates from all selected C2N alternatives
    for i = 1:size(C2N_results)[1]
        result = C2N_results[i, :]

        if (result[:project_type] == "new_xtr") && (result[:lag] == 0) && (result[:units_to_execute] != 0)
            # Record start of new construction project
                record_new_construction_projects(
                    settings,
                    result,
                    unit_specs,
                    db,
                    current_pd,
                    agent_id,
                )

            # Record CPP retirement(s)
            for j = 1:result[:num_cpp_rets]
                record_asset_retirements(
                    result,
                    db,
                    agent_id,
                    unit_specs,
                    current_pd,
                    mode = "C2N_newbuild",
                )

            end

        elseif result[:project_type] == "retirement"
            # Retire the nuclear unit
            if result[:units_to_execute] != 0
                record_asset_retirements(
                    result,
                    db,
                    agent_id,
                    unit_specs,
                    current_pd,
                )
            end
        end 

    end

    # Save the agent's decisions to the database
    save_agent_decisions(db, current_pd, agent_id, all_results)


end


function record_new_construction_projects(
    settings,
    result,
    unit_data,
    db,
    current_pd,
    agent_id,
)
    # Retrieve unit_specs data for this unit type
    unit_type_specs =
        filter(:unit_type => x -> x == result[:unit_type], unit_data)

    # Set default initial values
    cum_occ = (
        unit_type_specs[1, :overnight_capital_cost] *
        unit_type_specs[1, :capacity] *
        settings["constants"]["MW2kW"]
    )
    rcec = cum_occ
    cum_construction_exp = 0
    cum_construction_duration = unit_type_specs[1, :construction_duration]
    rtec = cum_construction_duration
    start_pd = current_pd
    completion_pd = current_pd + unit_type_specs[1, :construction_duration]
    cancellation_pd = settings["constants"]["distant_time"]
    retirement_pd = (
        current_pd +
        unit_type_specs[1, :construction_duration] +
        unit_type_specs[1, :unit_life]
    )
    total_capex = 0
    cap_pmt = 0
    anpe = (
        unit_type_specs[1, :overnight_capital_cost] *
        unit_type_specs[1, :capacity] *
        settings["constants"]["MW2kW"] /
        unit_type_specs[1, :construction_duration]
    )
    C2N_reserved = 0

    # Add a number of project instances equal to the 'units_to_execute'
    #   value from the decision solution
    for j = 1:result[:units_to_execute]
        next_id = get_next_asset_id(db)

        # Update the `WIP_updates` table
        WIP_projects_vals = (
            next_id,
            agent_id,
            current_pd,
            cum_occ,
            rcec,
            cum_construction_duration,
            rtec,
            cum_construction_exp,
            anpe,
        )
        DBInterface.execute(
            db,
            "INSERT INTO WIP_updates VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            WIP_projects_vals,
        )

        # Update the `asset_updates` table
        assets_vals = (
            next_id,
            agent_id,
            result[:unit_type],
            start_pd,
            completion_pd,
            cancellation_pd,
            retirement_pd,
            total_capex,
            cap_pmt,
            C2N_reserved,
        )
        DBInterface.execute(
            db,
            "INSERT INTO asset_updates VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            assets_vals,
        )
    end

end


function record_asset_retirements(
    result,
    db,
    agent_id,
    unit_specs,
    current_pd;
    mode = "standalone",
)
    if mode == "standalone"
        # Generate a list of assets which match 'result' on agent owner,
        #   unit type, and mandatory retirement date
        match_vals = (result[:unit_type], result[:ret_pd], agent_id)
        ret_candidates =
            DBInterface.execute(
                db,
                string(
                    "SELECT asset_id FROM assets WHERE unit_type = ? ",
                    "AND retirement_pd = ? AND agent_id = ? AND C2N_reserved = 0",
                ),
                match_vals,
            ) |> DataFrame

        # Set the number of units to execute
        units_to_execute = result[:units_to_execute]

        # Set the new retirement pd
        new_ret_pd = result[:lag]

    elseif mode == "C2N_newbuild"
        # Determine the period in which the coal units must retire
        unit_type_specs =
            filter(:unit_type => x -> x == result[:unit_type], unit_specs)[1, :]

        if occursin("C2N0", result[:unit_type])
            new_ret_pd = convert(
                Int64,
                (
                    current_pd +
                    result[:lag] + 
                    unit_type_specs[:construction_duration]
                )
            )
        else
            new_ret_pd = convert(
                Int64, 
                (
                    current_pd +
                    result[:lag] +
                    ceil(round(unit_type_specs[:cpp_ret_lead], digits = 3))
                )
            )
        end

        C2N_reserved = 0

        # Generate a list of coal units which will still be operational by
        #   the necessary period
        match_vals = (new_ret_pd, new_ret_pd, agent_id, C2N_reserved)
        ret_candidates =
            DBInterface.execute(
                db,
                string(
                    "SELECT asset_id, retirement_pd FROM assets WHERE unit_type = 'coal' ",
                    "AND completion_pd <= ? AND retirement_pd >= ? ",
                    "AND agent_id = ? AND C2N_reserved = ?",
                ),
                match_vals,
            ) |> DataFrame

        # Check the assets_updates table
        # If any units have been marked as C2N_reserved by another operation
        #   this round, ensure those units are not double-allocated
        claimed_coal_assets = DBInterface.execute(
            db,
            string("SELECT asset_id FROM asset_updates WHERE unit_type = 'coal' AND C2N_reserved = 1")
        ) |> DataFrame

        for k in 1:size(claimed_coal_assets)[1]
            ret_candidates = filter(:asset_id => x -> x != claimed_coal_assets[k, :asset_id], ret_candidates)
        end

        # Sort the final list of eligible coal retirement candidates from
        #   earliest initial retirement period to latest
        sort!(ret_candidates, :retirement_pd)

        # Set the number of units to execute
        units_to_execute =
            unit_type_specs["num_cpp_rets"] * result[:units_to_execute]
    end

    # Retire as many of these matching assets as is indicated by the agent
    #   optimization result
    for j = 1:units_to_execute
        asset_to_retire = ret_candidates[convert(Int64, j), :asset_id]
        asset_data =
            DBInterface.execute(
                db,
                "SELECT * FROM assets WHERE asset_id = $asset_to_retire",
            ) |> DataFrame

        # Overwrite the original record's retirement period with the new
        #   retirement period based on the C2N project's scheduling needs
        asset_data[1, :retirement_pd] = new_ret_pd
        # If this is a C2N-related coal retirement, mark it as such
        if mode == "C2N_newbuild"
            asset_data[1, :C2N_reserved] = 1
        end
        replacement_data = [item for item in asset_data[1, :]]

        # Save this new record to the asset_updates table
        DBInterface.execute(
            db,
            "INSERT INTO asset_updates VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            replacement_data,
        )
    end

end


function get_agent_portfolio_forecast(agent_id, db, current_pd, fc_pd)
    agent_portfolios = DataFrame()
    # Retrieve the agent's projected portfolios
    for y = current_pd:(current_pd + fc_pd)
        agent_portfolio =
            DBInterface.execute(
                db,
                string(
                    "SELECT unit_type, COUNT(unit_type) FROM assets ",
                    "WHERE agent_id = $agent_id AND completion_pd <= $y ",
                    "AND retirement_pd > $y AND cancellation_pd > $y ",
                    "GROUP BY unit_type",
                ),
            ) |> DataFrame

        rename!(agent_portfolio, Symbol("COUNT(unit_type)") => :num_units)
        agent_portfolio[!, :y] .= y

        append!(agent_portfolios, agent_portfolio)
    end

    return agent_portfolios
end


function update_agent_financial_statement(
    settings,
    agent_id,
    db,
    unit_specs,
    current_pd,
    fc_pd,
    dispatch_results,
    agent_params,
)
    # Filter out any dispatch results extending beyond the forecast period
    dispatch_results =
        filter(:y => y -> y <= current_pd + fc_pd, dispatch_results)

    # Retrieve horizontally-abbreviated dataframes
    short_unit_specs = select(unit_specs, [:unit_type, :FOM])

    # Retrieve the agent's complete portfolio forecast
    agent_portfolio_forecast =
        get_agent_portfolio_forecast(agent_id, db, current_pd, fc_pd)

    # Inner join the year's portfolio with financial pivot
    fin_results = coalesce.(outerjoin(
        dispatch_results,
        agent_portfolio_forecast,
        on = [:y, :unit_type],
    ), 0)

    # Compute total quantities based on agent's number of owned units
    transform!(
        fin_results,
        [:qty, :num_units] =>
            ((qty, num_units) -> qty .* num_units) => :total_qty,
    )
    fin_results =
        select(fin_results, [:y, :unit_type, :dispatch_result, :total_qty])

    # Add in FOM data for all unit types
    fin_results = add_FOM(
        settings,
        fin_results,
        agent_portfolio_forecast,
        short_unit_specs,
    )

    # Create the financial statement with unit types aggregated
    agent_fs = unstack(fin_results, :y, :dispatch_result, :total_qty)
    rename!(agent_fs, :y => :projected_pd)
    agent_fs[!, :base_pd] .= current_pd

    # TODO: in post-facto policies, recompute all PTC awards based on units' actual operational dates
    agent_fs = compute_policy_contributions(settings, db, deepcopy(agent_fs), agent_id, current_pd, unit_specs, agent_portfolio_forecast, dispatch_results)

    # Add in scheduled financing factors: depreciation, interest payments, and capex
    agent_fs =
        compute_scheduled_financing_factors(db, agent_fs, agent_id, current_pd)

    # Replace missing values with 0s
    for col in names(agent_fs)
        agent_fs[!, col] = coalesce.(agent_fs[!, col], 0)
    end

    # Propagate out accounting line items (EBITDA through FCF)
    agent_fs = compute_accounting_line_items(db, agent_fs, agent_id, agent_params, mode="agent")

    agent_fs = compute_credit_indicator_scores(agent_fs)

    # Save the dataframe to the database
    save_agent_fs!(agent_fs, agent_id, db)

    return agent_fs


end


function add_FOM(
    settings,
    fin_results,
    agent_portfolio_forecast,
    short_unit_specs,
)
    # Fill in total FOM
    # Inner join the FOM table with the agent portfolios
    FOM_df =
        innerjoin(agent_portfolio_forecast, short_unit_specs, on = :unit_type)
    transform!(
        FOM_df,
        [:FOM, :num_units] =>
            (
                (FOM, num_units) ->
                    FOM .* num_units .* settings["constants"]["MW2kW"]
            ) => :total_FOM,
    )

    # Restructure the FOM dataframe to match the other results
    FOM_df = combine(
        groupby(FOM_df, [:y, :unit_type]),
        [:num_units, :total_FOM] .=> sum;
        renamecols = false,
    )
    rename!(FOM_df, :total_FOM => :total_qty)
    FOM_df[!, :dispatch_result] .= "FOM"
    FOM_df = select(FOM_df, [:y, :unit_type, :dispatch_result, :total_qty])

    # Combine the FOM dataframe into the other financial results from dispatch
    fin_results = vcat(fin_results, FOM_df)
    fin_results = combine(
        groupby(fin_results, [:y, :dispatch_result]),
        :total_qty => sum;
        renamecols = false,
    )

    return fin_results
end


function compute_policy_contributions(settings, db, fs_copy, agent_id, current_pd, unit_specs, agent_portfolio_forecast, dispatch_results)
    # Check whether any policies are specified in the settings file
    if "policies" in keys(settings["scenario"])
        policies = settings["scenario"]["policies"]

        # Initialize the tax credits column with all zeroes
        fs_copy[!, :tax_credits] = zeros(size(fs_copy)[1])

        # Get a list of all WIP and operating assets owned by this agent
        stmt = string(
            "SELECT * FROM assets WHERE agent_id = $agent_id ",
            "AND cancellation_pd > $current_pd ",
            "AND retirement_pd > $current_pd",
        )
        assets = DBInterface.execute(db, stmt) |> DataFrame

        stmt = "SELECT * FROM WIP_projects WHERE agent_id = $agent_id AND period = $current_pd"
        WIP_projects = DBInterface.execute(db, stmt) |> DataFrame

        # Iterate through all assets and individually compute policy incentives
        #   for which it qualifies
        for asset in eachrow(assets)
            # If PTC exists
            if haskey(policies, "PTC") && policies["PTC"]["enabled"]
                # If unit is listed as PTC-eligible
                if asset[:unit_type] in policies["PTC"]["eligibility"]["unit_type"]
                    # Determine asset's completion date
                    if asset["completion_pd"] == 0
                        # If the asset starts the simulation already operational,
                        #   assume that its completion period was its
                        #   retirement date minus the unit life
                        completion_pd = asset["retirement_pd"] - filter(:unit_type => unit_type -> unit_type == asset["unit_type"], unit_specs)[1, :unit_life]
                    else
                        completion_pd = asset["completion_pd"]
                    end
                    # If unit's completion date is within the PTC eligibility period
                    if (policies["PTC"]["start_year"] <= completion_pd) && (policies["PTC"]["end_year"] >= completion_pd)
                        # If current period is still within the award period for the PTC
                        if current_pd <= completion_pd + policies["PTC"]["duration"]
                            # Calculate the start and end of the award period,
                            #   in relative years
                            award_start = max(completion_pd - current_pd + 1, policies["PTC"]["start_year"] - current_pd + 1, 1)
                            award_end = min(completion_pd + policies["PTC"]["duration"] - current_pd , size(fs_copy)[1])

                            for y = award_start:award_end
                                # Get asset type generation for the year
                                asset_generation = filter([:unit_type, :y, :dispatch_result] => ((unit_type, y, disp_res) -> ((unit_type == asset["unit_type"]) && (y == y) && (disp_res == "generation"))), dispatch_results)[1, :qty]

                                # Add the asset's total PTC earnings 
                                #   contribution to the year's total
                                fs_copy[y, :tax_credits] += asset_generation * policies["PTC"]["qty"]
                            end
                        end
                    end
                end
            end

            # If ITC exists
            if haskey(policies, "ITC") && policies["ITC"]["enabled"]
                # If unit is listed as ITC-eligible
                if asset["unit_type"] in policies["ITC"]["eligibility"]["unit_type"]
                    # If unit's completion date is within the ITC eligibility period
                    if (policies["ITC"]["start_year"] <= asset["completion_pd"]) && (policices["ITC"]["end_year"] >= asset["completion_pd"])
                        # If unit is operational, use its total_capex value
                        if asset["completion_pd"] <= current_pd
                            total_capex = asset["total_capex"]
                        else
                        # If unit is not operational, use its total estimated OCC
                        # TODO: estimate fully financed cost
                            WIP_data = filter(:asset_id => asset_id -> asset_id == asset["asset_id"], WIP_projects)[1, :]
                            total_capex = WIP_data["cum_occ"]
                        end

                        # Add the appropriate ITC award to the tax_credits column
                        fs_copy[asset["completion_pd"] + 1, :tax_credits] += total_capex * policies["ITC"]["qty"]
                    end
                end
            end

            if (occursin("C2N", asset["unit_type"])) && asset["completion_pd"] >= current_pd + 1
                unit_type_data = filter(:unit_type => unit_type -> unit_type == asset["unit_type"], unit_specs)[1, :]
                total_capex = unit_type_data[:overnight_capital_cost] * unit_type_data[:capacity] * 1000
                


                fs_copy[asset["completion_pd"] + 1, :tax_credits] += 0.4 * total_capex
            end
        end
    end

    return fs_copy
end


function compute_scheduled_financing_factors(db, agent_fs, agent_id, current_pd)
    # Fill in total capex
    capex_ts =
        DBInterface.execute(
            db,
            string(
                "SELECT projected_pd, SUM(capex) ",
                "FROM capex_projections ",
                "WHERE agent_id = $agent_id ",
                "AND base_pd = $current_pd GROUP BY projected_pd",
            ),
        ) |> DataFrame

    remaining_debt_principal_ts = compute_debt_principal_ts(db, agent_id, current_pd)

    # Fill in total depreciation
    depreciation_ts =
        DBInterface.execute(
            db,
            string(
                "SELECT projected_pd, SUM(depreciation) ",
                "FROM depreciation_projections ",
                "WHERE agent_id = $agent_id ",
                "AND base_pd = $current_pd ",
                "GROUP BY projected_pd",
            ),
        ) |> DataFrame

    # Fill in total interest payments
    interest_payment_ts =
        DBInterface.execute(
            db,
            string(
                "SELECT projected_pd, SUM(interest_payment) ",
                "FROM financing_schedule ",
                "WHERE agent_id = $agent_id ",
                "AND base_pd = $current_pd ",
                "GROUP BY projected_pd",
            ),
        ) |> DataFrame

    # Join the scheduled columns to the financial statement
    agent_fs = leftjoin(agent_fs, capex_ts, on = :projected_pd)
    agent_fs = leftjoin(agent_fs, remaining_debt_principal_ts, on = :projected_pd)
    agent_fs = leftjoin(agent_fs, depreciation_ts, on = :projected_pd)
    agent_fs = leftjoin(agent_fs, interest_payment_ts, on = :projected_pd)

    # Standardize column names
    rename!(
        agent_fs,
        Symbol("SUM(capex)") => :capex,
        Symbol("SUM(depreciation)") => :depreciation,
        Symbol("SUM(interest_payment)") => :interest_payment,
    )

    return agent_fs
end


function compute_debt_principal_ts(db, agent_id, current_pd)
    # Get debt principal repayments by year
    principal_payments = DBInterface.execute(db, string("SELECT projected_pd, SUM(principal_payment) FROM financing_schedule WHERE agent_id = $agent_id AND base_pd = $current_pd GROUP BY projected_pd")) |> DataFrame
    rename!(principal_payments, Symbol("SUM(principal_payment)") => :principal_payment)

    # Get cumulative debt issued by year
    # Get the list of all debt instruments issued by this agent
    inst_manifest = DBInterface.execute(db, string("SELECT * FROM financial_instrument_manifest WHERE agent_id = $agent_id AND instrument_type='debt'")) |> DataFrame

    debt_principal_ts = DataFrame(projected_pd = Int[], remaining_debt_principal = Float64[])

    if size(principal_payments)[1] > 0
        for y=current_pd:maximum(principal_payments[!, :projected_pd])
            debt_by_year = filter(:pd_issued => pd -> pd <= y, inst_manifest)
            cum_debt = sum(debt_by_year[!, :initial_principal])
            push!(debt_principal_ts, [y, cum_debt])
        end
    end

    debt_principal_ts = outerjoin(debt_principal_ts, principal_payments, on = :projected_pd)
    transform!(debt_principal_ts, [:remaining_debt_principal, :principal_payment] => ((rem_prin, prin_pmt) -> rem_prin - prin_pmt) => :remaining_debt_principal)

    return debt_principal_ts
end


function compute_accounting_line_items(db, agent_fs, agent_id, agent_params; mode=nothing)
    ### Computed FS quantities
    # EBITDA
    transform!(
        agent_fs,
        [:revenue, :VOM, :FOM, :fuel_cost, :policy_adj, :tax_credits] =>
            ((rev, VOM, FOM, FC, pol, tc) -> (rev - VOM - FOM - FC + pol + tc)) =>
                :EBITDA,
    )

    # EBIT
    transform!(
        agent_fs,
        [:EBITDA, :depreciation] => ((EBITDA, dep) -> (EBITDA - dep)) => :EBIT,
    )

    # EBT
    transform!(
        agent_fs,
        [:EBIT, :interest_payment] =>
            ((EBIT, interest) -> EBIT - interest) => :EBT,
    )

    # Retrieve some system parameters from the database
    model_params = DBInterface.execute(db, "SELECT * FROM model_params") |> DataFrame
    tax_rate = filter(:parameter => x -> x == "tax_rate", model_params)[1, :value]

    tax_credits_discount = filter(:parameter => x -> x == "tax_credits_discount", model_params)[1, :value]

    # Compute nominal tax owed
    # If this is a subproject, allow negative tax owed, and assume no tax credits are resold
    if mode == "subproject"
        transform!(agent_fs, :EBT => ((EBT) -> EBT * tax_rate) => :tax_owed)
        transform!(agent_fs, [:tax_owed, :tax_credits] => ((T, C) -> T - C) => :realized_tax_owed)
        agent_fs[!, :tax_credit_sale_revenue] .= 0.0
    # If this is an agent, require tax to be nonnegative, and allow resale of tax credits if applicable
    else
        transform!(agent_fs, :EBT => ByRow(EBT -> ifelse.(EBT >= 0, EBT * tax_rate, 0)) => :tax_owed)
        transform!(agent_fs, [:tax_owed, :tax_credits] => ByRow((T, C) -> ifelse.(T >= C, T - C, 0)) => :realized_tax_owed)
        transform!(agent_fs, [:tax_owed, :tax_credits] => ByRow((T, C) -> ifelse.(C > T, (C - T) * (1 - tax_credits_discount) * (1 - tax_rate), 0)) => :tax_credit_sale_revenue)
    end

    # Net Income
    transform!(agent_fs, [:EBT, :realized_tax_owed, :tax_credit_sale_revenue] => ((EBT, real_tax, C_salerev) -> (EBT - real_tax + C_salerev)) => :net_income)

    # Free Cash Flow
    transform!(
        agent_fs,
        [:net_income, :depreciation, :capex] =>
            ((NI, dep, capex) -> (NI + dep - 0.5 * capex)) => :nominal_FCF,
    )

    # Dividends and Retained Earnings
    # Set up data
    agent_id = agent_params[1, "agent_id"]
    cost_of_equity =
        DBInterface.execute(
            db,
            "SELECT cost_of_equity FROM agent_params WHERE agent_id = $agent_id",
        ) |> DataFrame
    cost_of_equity = cost_of_equity[1, :cost_of_equity]

    # Compute dividends
    # Allow negative dividends if this is a project (negative FCF accruing to
    #   the project is assumed to reduce total dividend across the agent's
    #   finances)
    if mode == "subproject"
        transform!(agent_fs, :nominal_FCF => ((nfcf) -> nfcf .* cost_of_equity) => :dividends)
    else
        # If this is the agent's financial statement, the minimum actual
        #   dividend is zero
        agent_fs.dividends = ifelse.(agent_fs.nominal_FCF .> 0, agent_fs.nominal_FCF .* cost_of_equity, 0)
    end

    # Compute annual contribution to retained earnings and finalized FCF
    if mode == "subproject"
        # If this is a subproject, nominal FCF = FCF because retained earnings
        #   are not used to tide over cash flow
        agent_fs[!, :FCF] = agent_fs[!, :nominal_FCF]

        # Record each year's delta-RE as the RE value
        transform!(
            agent_fs,
            [:FCF, :dividends] => ((fcf, dividends) -> fcf .- dividends) => :delta_RE,
        )
        agent_fs[!, :retained_earnings] = agent_fs[!, :delta_RE]
    else
        # If this is the agent, track the full net change in RE starting 
        #   from the current year's RE value
        transform!(
            agent_fs,
            [:nominal_FCF, :dividends] => ((nfcf, div) -> nfcf .- div) => :delta_RE,
        )

        # Get the preceding period's RE
        hist_agent_fss = DBInterface.execute(db, "SELECT base_pd, projected_pd, retained_earnings FROM agent_financial_statements WHERE agent_id == $agent_id") |> DataFrame

        if size(hist_agent_fss)[1] == 0
            # If there are no financial statements recorded for this agent yet,
            #   this must be the first simulation period; use data from
            #   agent_params
            init_RE = agent_params[1, :starting_RE]
        else
            # If some financial statements were found, choose the most recent
            #   immediate-year prediction
            re_vals = filter([:base_pd, :projected_pd] => ((base_pd, pr_pd) -> base_pd == pr_pd), hist_agent_fss)
            init_RE = filter(:base_pd => base_pd -> base_pd == maximum(hist_agent_fss[!, :base_pd]), re_vals)[1, :retained_earnings]
        end

        # Set up a blank column for retained earnings
        agent_fs[!, :retained_earnings] .= 0.0
        agent_fs[!, :FCF] .= 0.0

        # Propagate the running retained earnings total through the agent
        #   financial statement, while using extra RE to bolster negative FCF
        running_RE = init_RE
        for i=1:size(agent_fs)[1]
            if (running_RE > 0) && (agent_fs[i, :nominal_FCF] < 0)
                # If the nominal FCF is negative, and there is RE cash
                #   available, use it to make up the shortfall (as far as possible)
                RE_transfer = min(running_RE, (-1) * agent_fs[i, :nominal_FCF])

                FCF = agent_fs[i, :nominal_FCF] + RE_transfer
                running_RE = running_RE - RE_transfer
            else
                # If the nominal FCF is positive, there is no transfer
                #   and RE and FCF take on their normal values
                FCF = agent_fs[i, :nominal_FCF]
                running_RE = running_RE + agent_fs[i, :delta_RE]
            end

            agent_fs[i, :retained_earnings] = running_RE
            agent_fs[i, :FCF] = FCF
        end
    end

    return agent_fs
end

function compute_credit_indicator_scores(agent_fs)
    # Credit ratings indicators
    # Compute ICR metric column
    transform!(agent_fs, [:interest_payment, :FCF] => ((int, fcf) -> (fcf .+ int) ./ int) => :ICR)

    # Compute FCF-to-debt metric column
    transform!(agent_fs, [:remaining_debt_principal, :FCF] => ((debt, fcf) -> fcf ./ debt) => :FCF_debt_ratio)

    # Compute retained earnings-to-debt metric column
    transform!(agent_fs, [:remaining_debt_principal, :retained_earnings] => ((debt, re) -> re ./ debt) => :RE_debt_ratio)

    # Compute scaled Moody's indicator score from metrics
    transform!(agent_fs, [:ICR, :FCF_debt_ratio, :RE_debt_ratio] => ((icr, fcf_debt, re_debt) -> compute_moodys_score.(icr, fcf_debt, re_debt)) => :moodys_score)

    return agent_fs
end



function save_agent_fs!(fs, agent_id, db)
    # Add the agent id to the dataframe
    fs[!, :agent_id] .= agent_id

    # Get the column order from the agent_financial_statements database table
    fs_col_order =
        DBInterface.execute(
            db,
            string(
                "SELECT name FROM ",
                "PRAGMA_TABLE_INFO('agent_financial_statements')",
            ),
        ) |> DataFrame
    fs_col_order = collect(fs_col_order[!, "name"]) # convert DF col to vector

    # Reorder the agent fs to match the DB table
    fs = select(fs, fs_col_order)

    # Set up the SQL "(?, ?, ... , ?)" string of correct length
    fill_tuple = string("(", repeat("?, ", size(fs_col_order)[1] - 1), "?)")

    for row in Tuple.(eachrow(fs))
        DBInterface.execute(
            db,
            string(
                "INSERT INTO agent_financial_statements VALUES $fill_tuple"
            ),
            row,
        )
    end

end

function save_agent_decisions(db, current_pd, agent_id, decision_df)
    decision_df[!, :agent_id] .= agent_id
    decision_df[!, :base_pd] .= current_pd
    cols_to_ignore = [:uid, :current]
    select!(decision_df, [:agent_id, :base_pd], Not(vcat([:agent_id], cols_to_ignore)))
    for row in Tuple.(eachrow(decision_df))
        DBInterface.execute(
            db,
            string(
                "INSERT INTO agent_decisions ",
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ),
            row,
        )
    end
end


function display_agent_choice_results(CLI_args, m, all_results)
    status = string(termination_status.(m))
    @debug "Model termination status: $status"
    println("=== Results: ===")

    agent_id = CLI_args["agent_id"]

    if CLI_args["verbosity"] == 2
        msg = nothing
        units_to_execute = nothing

        if status == "OPTIMAL"
            results = filter(:units_to_execute => u -> u > 0, all_results)

            if size(results)[1] == 0
                msg = "Agent $agent_id's optimal decision is to take no actions this turn."
            else
                msg = "Project alternatives to execute:"
                units_to_execute = results
            end
        else
            msg = "No feasible solution found for the decision optimization problem.\nAgent $agent_id will take no actions this turn."
        end

        if msg != nothing
            @info msg
        end
        if units_to_execute != nothing
            @info units_to_execute
        end

    elseif CLI_args["verbosity"] == 3
        @debug "Alternatives to execute:"
        @debug all_results
    end

end

end
