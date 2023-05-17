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

module ABCEfunctions

using ArgParse,
    Requires, SQLite, DataFrames, CSV, JuMP, GLPK, Cbc, Logging, Tables, HiGHS

# Use CPLEX if available
function __init__()
    @require CPLEX = "a076750e-1247-5638-91d2-ce28b192dca0" @eval using CPLEX
end

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
                       "the correct directory."
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
            unit_type_data = filter(:unit_type => x -> x == unit_type,
                                    unit_specs
                             )[1, :]

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
                C2N_specs
            )

            unit_specs[(unit_specs.unit_type .== unit_type), :construction_duration] .= size(capex_tl)[1]
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
        "FROM assets WHERE agent_id = ", agent_id,
        " AND completion_pd <= ", pd,
        " AND cancellation_pd > ", pd,
        " AND retirement_pd > ", pd,
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
    beta[2] = ((beta[2] * size(y)[1] 
                + settings["demand"]["historical_demand_growth_rate"] 
                  * (demand_projection_window - size(y)[1])
               ) / demand_projection_window
              )

    # Project future demand
    proj_horiz = fc_pd - size(visible_demand)[1]
    x_proj = hcat(
                 ones(proj_horiz),
                 [i for i=size(visible_demand)[1]+size(demand_history)[1]+1:fc_pd+size(demand_history)[1]]
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


function get_net_demand(
    db,
    pd,
    fc_pd,
    demand_forecast,
    system_portfolios,
    unit_specs,
)
    # Calculate the amount of forecasted net demand in future periods
    installed_cap_forecast =
        DataFrame(period = Int64[], derated_capacity = Float64[])

    vals = (pd, pd)

    total_caps = DataFrame(
                     period = Int64[],
                     total_eff_cap = Float64[]
                 )

    for i = pd:maximum(keys(system_portfolios))
        year_portfolio = innerjoin(
            system_portfolios[i],
            unit_specs,
            on = :unit_type,
        )
        transform!(
            year_portfolio,
            [:num_units, :capacity, :capacity_factor] =>
                ((num_units, cap, CF) -> num_units .* cap .* CF) =>
                    :effective_capacity,
        )

        total_year_cap = sum(year_portfolio[!, :effective_capacity])
        push!(total_caps, [i, total_year_cap])
    end

    demand_forecast = innerjoin(demand_forecast, total_caps, on = :period)
    transform!(
        demand_forecast,
        [:total_demand, :total_eff_cap] =>
            ((demand, eff_cap) -> demand - eff_cap) => :net_demand,
    )

    return demand_forecast
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
    agent_params,
    db,
    current_pd,
    long_econ_results,
    C2N_specs,
)
    PA_uids = create_PA_unique_ids(settings, unit_specs, asset_counts)

    PA_fs_dict = create_PA_pro_formas(PA_uids, fc_pd)

    PA_uids, PA_fs_dict = populate_PA_pro_formas(
        settings,
        PA_uids,
        PA_fs_dict,
        unit_specs,
        fc_pd,
        agent_params,
        db,
        current_pd,
        long_econ_results,
        C2N_specs,
    )

    return PA_uids, PA_fs_dict

end

function create_PA_unique_ids(settings, unit_specs, asset_counts)
    PA_uids = DataFrame(
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
                PA_uids,
                [unit_type project_type lag nothing uid NPV allowed]
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
        )[!, :retirement_pd,]

        for ret_pd in ret_pds
            for lag = 0:settings["agent_opt"]["num_future_periods_considered"]
                row = [unit_type project_type lag ret_pd uid NPV allowed]
                push!(PA_uids, row)
                uid += 1
            end
        end
    end

    return PA_uids

end


function create_PA_pro_formas(PA_uids, fc_pd)
    PA_fs_dict = Dict()
    for uid in PA_uids[!, :uid]
        PA_fs_dict[uid] = DataFrame(
            year = 1:fc_pd,
            capex = zeros(fc_pd),
            remaining_debt_principal = zeros(fc_pd),
            debt_payment = zeros(fc_pd),
            interest_payment = zeros(fc_pd),
            depreciation = zeros(fc_pd),
            gen = zeros(fc_pd),
        )
    end

    return PA_fs_dict

end


function populate_PA_pro_formas(
    settings,
    PA_uids,
    PA_fs_dict,
    unit_specs,
    fc_pd,
    agent_params,
    db,
    current_pd,
    long_econ_results,
    C2N_specs,
)
    for uid in PA_uids[!, :uid]
        # Retrieve current project alternative definition for convenience
        current_PA = filter(:uid => x -> x == uid, PA_uids)[1, :]

        # Retrieve relevant unit type specs for convenience
        unit_type_data = filter(
            :unit_type => 
                x -> x == current_PA[:unit_type],
            unit_specs
            )[1, :]

        # If this is a potential new construction project, compute series
        #   related to construction costs
        if current_PA[:project_type] == "new_xtr"
            # Generate this alternative's expected construction
            #   expenditure profile

            PA_fs_dict[uid][!, :capex] = generate_capex_profile(
                db,
                settings,
                current_pd,
                unit_type_data,
                current_PA[:lag],
                fc_pd,
                C2N_specs,
            )

            # Set up the time-series of outstanding debt based on this
            #   construction expenditure profile
            set_initial_debt_principal_series(
                PA_fs_dict[uid],
                unit_type_data,
                current_PA[:lag],
                agent_params,
            )

            # Generate "prime movers": debt payments and depreciation
            generate_prime_movers(
                unit_type_data,
                PA_fs_dict[uid],
                current_PA[:lag],
                agent_params[1, :cost_of_debt],
            )
        end

        # Forecast unit revenue ($/period) and generation (kWh/period)
        forecast_unit_revenue_and_gen(
            settings,
            unit_type_data,
            PA_fs_dict[uid],
            db,
            current_pd,
            current_PA[:lag],
            long_econ_results,
            mode = current_PA[:project_type],
            orig_ret_pd = current_PA[:ret_pd],
        )

        # Forecast unit operating costs: VOM, fuel cost, and FOM
        forecast_unit_op_costs(
            settings,
            unit_type_data,
            PA_fs_dict[uid],
            current_PA[:lag],
            mode = current_PA[:project_type],
            orig_ret_pd = current_PA[:ret_pd],
        )

        # Propagate the accounting logic flow
        propagate_accounting_line_items(PA_fs_dict[uid], db)

        # If the current project is a retirement project, adjust the FS to show
        #   the difference between the original retirement plan and the 
        #   proposed outcome
        if current_PA[:project_type] == "retirement"
            PA_fs_dict[uid] =
                compute_FS_delta_value(PA_fs_dict[uid], current_PA[:lag], db)
        end

        FCF_NPV, PA_fs_dict[uid] =
            compute_alternative_NPV(PA_fs_dict[uid], agent_params)
        file_name = joinpath(
            "tmp",
            string(
                current_PA[:unit_type], "_",
                current_PA[:project_type], "_lag",
                current_PA[:lag], "_fs.csv"
            )
        )
        CSV.write(file_name, PA_fs_dict[uid])

        # Save the NPV result
        filter(:uid => x -> x == uid, PA_uids, view = true)[1, :NPV] = FCF_NPV

    end

    return PA_uids, PA_fs_dict

end


function create_FS_dict(data, fc_pd, num_lags; mode = "new_xtr")
    fs_dict = Dict()
    num_alts = size(data)[1]
    for i = 1:num_alts
        for j = 0:num_lags
            if mode == "retirement"
                ret_pd = data[i, :retirement_pd]
                project_type = "retirement"
            elseif mode == "new_xtr"
                ret_pd = "X"
                project_type = "new_xtr"
            end

            unit_name = string(
                data[i, :unit_type], "_", project_type, "_", ret_pd, "_lag-", j
            )

            unit_FS = DataFrame(
                year = 1:fc_pd,
                capex = zeros(fc_pd),
                gen = zeros(fc_pd),
                remaining_debt_principal = zeros(fc_pd),
                debt_payment = zeros(fc_pd),
                interest_payment = zeros(fc_pd),
                depreciation = zeros(fc_pd),
            )

            fs_dict[unit_name] = unit_FS
        end
    end
    return fs_dict
end


"""
    generate_capex_profile(unit_type_data, lag)

Creates a uniform expenditure profile for a potential construction project,
and appends appropriate numbers of leading and lagging zeros (no construction 
expenditures outside the construction period).
"""
function generate_capex_profile(
    db,
    settings,
    current_pd,
    unit_type_data,
    lag,
    fc_pd,
    C2N_specs,
)
    # No construction expenditures before the project begins
    head_zeros = zeros(lag)

    if occursin("C2N", unit_type_data[:unit_type])
        capex_tl, activity_schedule = project_C2N_capex(
            db,
            settings,
            unit_type_data,
            lag,
            fc_pd,
            current_pd,
            C2N_specs,
        )
        capex = capex_tl[!, :total_capex]
    else
        # Uniformly distribute projected total project costs over the
        #   construction period
        capex_per_pd = (
            unit_type_data[:overnight_capital_cost] *
            unit_type_data[:capacity] *
            settings["constants"]["MW2kW"] /
            unit_type_data[:construction_duration]
        )
        capex =
            ones(
                convert(
                    Int64,
                    ceil(
                        round(
                            unit_type_data[:construction_duration],
                            digits = 3,
                        ),
                    ),
                ),
            ) .* capex_per_pd
    end

    # No construction expenditures from the end of construction until the end
    #   of the model's forecast period
    tail_zeros = zeros(fc_pd - lag - size(capex)[1])

    # Concatenate the above series into one
    capex_column = vcat(head_zeros, capex, tail_zeros)

    return capex_column
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
            data = C2N_specs["greenfield"]["PWR"]
        elseif occursin("HTGR", unit_type_data[:unit_type])
            rxtr_type = "HTGR"
            data = C2N_specs["greenfield"]["HTGR"]
        else
            rxtr_type = "SFR"
            data = C2N_specs["greenfield"]["SFR"]
        end
    else
        if unit_type_data[:unit_type] == "C2N1"
            conversion_type = "electrical"
            rxtr_type = "PWR"
        elseif unit_type_data[:unit_type] == "C2N2"
            conversion_type = "steam_noTES"
            rxtr_type = "HTGR"
        else
            conversion_type = "steam_TES"
            rxtr_type = "SFR"
        end
        data = C2N_specs[conversion_type][assumption]
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


"""
    set_initial_debt_principal_series(unit_fs, unit_type_data, lag, agent_params)

Set up the record of the accrual of debt during construction.
"""
function set_initial_debt_principal_series(
    unit_fs,
    unit_type_data,
    lag,
    agent_params,
)
    duration = convert(
        Int64,
        ceil(round(unit_type_data[:construction_duration], digits = 3)),
    )

    for i = (lag + 1):(lag + duration)
        unit_fs[i, :remaining_debt_principal] =
            (sum(unit_fs[1:i, :capex]) * agent_params[1, :debt_fraction])
    end
end


"""
    generate_prime_movers(unit_type_data, unit_fs, lag, d)

Generate the "prime mover" time series during the operating life of the plant:
  - debt payments
  - interest due
  - remaining debt principal
  - depreciation
Revenue is calculated in a separate function.
"""
function generate_prime_movers(unit_type_data, unit_fs, lag, cod)
    unit_d_x = convert(
        Int64,
        ceil(round(unit_type_data[:construction_duration], digits = 3)),
    )
    unit_op_life = unit_type_data[:unit_life]
    financing_term = 30
    rev_head_start =
        ceil(Int64, round(unit_type_data[:rev_head_start], digits = 3))

    # Find the end of the capital expenditures period
    capex_end = 0
    for i = 1:size(unit_fs)[1]
        if round(unit_fs[i, :capex], digits = 3) > 0
            capex_end = i
        end
    end

    for i =
        (capex_end + 1 - rev_head_start):(capex_end + 1 - rev_head_start + financing_term)
        # Apply a constant debt payment (sinking fund at cost of debt), based
        #   on the amount of debt outstanding at the end of the xtr project
        unit_fs[i, :debt_payment] = (
            unit_fs[lag + unit_d_x, :remaining_debt_principal] .* cod ./
            (1 - (1 + cod) .^ (-1 * financing_term))
        )

        # Determine the portion of each payment which pays down interest
        #   (instead of principal)
        unit_fs[i, :interest_payment] =
            (unit_fs[i - 1, :remaining_debt_principal] * cod)

        # Update the amount of principal remaining at the end of the period
        unit_fs[i, :remaining_debt_principal] = (
            unit_fs[i - 1, :remaining_debt_principal] -
            (unit_fs[i, :debt_payment] - unit_fs[i, :interest_payment])
        )
    end

    dep_term = 20
    for i =
        (capex_end + 1 - rev_head_start):(capex_end + 1 - rev_head_start + dep_term)
        # Apply straight-line depreciation, based on debt outstanding at the
        #   project's completion: all units have a dep. life of 20 years
        unit_fs[i, :depreciation] =
            (unit_fs[capex_end, :remaining_debt_principal] ./ dep_term)
    end

end


"""
    forecast_unit_revenue_and_gen(settings, unit_type_data, unit_fs, db, current_pd, lag, long_econ_results; mode, ret_pd)

Forecast the revenue which will be earned by the current unit type,
and the total generation (in kWh) of the unit, per time period.
"""
function forecast_unit_revenue_and_gen(
    settings,
    unit_type_data,
    unit_fs,
    db,
    current_pd,
    lag,
    long_econ_results;
    mode = "new_xtr",
    orig_ret_pd,
)
    # Compute the original retirement period
    # Minimum of ret_pd or size of the unit FS (i.e. the forecast period)
    if orig_ret_pd == nothing
        orig_ret_pd = size(unit_fs)[1]
    end

    # Load past years' dispatch results from A-LEAF
    ALEAF_dispatch_results =
        DBInterface.execute(db, "SELECT * FROM ALEAF_dispatch_results") |>
        DataFrame
    ALEAF_results =
        average_historical_ALEAF_results(settings, ALEAF_dispatch_results)

    # Compute the unit's total generation for each period, in kWh
    compute_total_generation(
        settings,
        current_pd,
        unit_type_data,
        unit_fs,
        lag,
        long_econ_results,
        ALEAF_results["wtd_hist_gens"];
        mode = mode,
        orig_ret_pd = orig_ret_pd,
    )

    # Compute total projected revenue, with VRE adjustment if appropriate, and
    #   save to the unit financial statement
    compute_total_revenue(
        settings,
        current_pd,
        unit_type_data,
        unit_fs,
        lag,
        long_econ_results,
        ALEAF_results["wtd_hist_revs"];
        mode = mode,
        orig_ret_pd = orig_ret_pd,
    )

end


function average_historical_ALEAF_results(settings, ALEAF_dispatch_results)
    wtd_hist_revs = nothing
    wtd_hist_gens = nothing

    if size(ALEAF_dispatch_results)[1] != 0
        # Get a list of all unique years in the dataframe
        years = unique(ALEAF_dispatch_results, :period)[!, :period]

        hist_decay = settings["dispatch"]["hist_decay"]

        # Compute scaling factor to normalize to 1
        c = 1 / sum(1 / (1 + hist_decay)^k for k in years)

        # Add a column with the diminishing historical weighting factor
        transform!(
            ALEAF_dispatch_results,
            [:period] => ((pd) -> c ./ (1 + hist_decay) .^ pd) => :hist_wt,
        )

        # Get weighted total revenue and generation
        transform!(
            ALEAF_dispatch_results,
            [:total_rev, :hist_wt] =>
                ((rev, wt) -> rev .* hist_decay) => :wtd_total_rev,
        )
        transform!(
            ALEAF_dispatch_results,
            [:gen_total, :hist_wt] =>
                ((gen, wt) -> gen .* hist_decay) => :wtd_total_gen,
        )

        wtd_hist_revs = unstack(
            select(ALEAF_dispatch_results, [:unit_type, :wtd_total_rev]),
            :unit_type,
            :wtd_total_rev,
            combine = sum,
        )

        wtd_hist_gens = unstack(
            select(ALEAF_dispatch_results, [:unit_type, :wtd_total_gen]),
            :unit_type,
            :wtd_total_gen,
            combine = sum,
        )

    end

    hist_results = Dict(
                       "wtd_hist_revs" => wtd_hist_revs,
                       "wtd_hist_gens" => wtd_hist_gens
                   )


    return hist_results

end


"""
    compute_total_revenue(unit_type_data, unit_fs, lag; mode, orig_ret_pd)

Compute the final projected revenue stream for the current unit type, adjusting
unit availability if it is a VRE type.
"""
function compute_total_revenue(
    settings,
    current_pd,
    unit_type_data,
    unit_fs,
    lag,
    long_econ_results,
    wtd_hist_revs;
    mode,
    orig_ret_pd = 9999,
)
    # Helpful short variables
    unit_d_x = convert(
        Int64,
        ceil(round(unit_type_data[:construction_duration], digits = 3)),
    )
    unit_op_life = unit_type_data[:unit_life]

    # Add a Revenue column to the financial statement dataframe
    unit_fs[!, :Revenue] .= 0.0

    # In "new_xtr" mode, revenues START accruing after the lag plus
    #   construction duration
    # In "retirement" mode, revenues start in the current period and CEASE
    #   accruing after the soonest of {the lag, the mandatory retirement
    #   period, the end of the forecast period}.
    if mode == "new_xtr"
        # Find the end of the capital expenditures period
        capex_end = 0
        for i = 1:size(unit_fs)[1]
            if round(unit_fs[i, :capex], digits = 3) > 0
                capex_end = i
            end
        end
        rev_start = (
            capex_end + 1 -
            ceil(Int64, round(unit_type_data[:rev_head_start], digits = 3))
        )
        rev_end = min(rev_start + unit_op_life, size(unit_fs)[1])
    elseif mode == "retirement"
        rev_start = 1
        # The maximum end of revenue period is capped at the length of the
        #   unit_fs dataframe (as some units with unspecified retirement
        #   periods default to a retirement period of 9999).
        rev_end = min(orig_ret_pd, size(unit_fs)[1] - 1)
    end

    # Compute final projected revenue series
    agg_econ_results = combine(
        groupby(long_econ_results, [:y, :unit_type]),
        [:annualized_rev_per_unit, :annualized_policy_adj_per_unit] .=> sum,
        renamecols = false,
    )

    if wtd_hist_revs == nothing
        hist_rev = 0
        hist_wt = 0
    else
        try
            hist_rev = wtd_hist_revs[1, unit_type_data[:unit_type]]
            hist_wt = settings["dispatch"]["hist_wt"]
        catch
            hist_rev = 0
            hist_wt = 0
        end
    end

    for y = rev_start:rev_end
        row = filter(
            [:y, :unit_type] =>
                ((t, unit_type) -> (t == y) && (unit_type == unit_type_data[:unit_type])),
            agg_econ_results,
        )

        if size(row)[1] != 0
            if hist_wt == 0
                wt = 0
            else
                wt = hist_wt^(y - current_pd)
            end
            unit_fs[y, :Revenue] = (
                (1 - hist_wt) * (
                    row[1, :annualized_rev_per_unit] +
                    row[1, :annualized_policy_adj_per_unit]
                ) + hist_wt * hist_rev
            )
        else
            unit_fs[y, :Revenue] = 0
        end
    end

    avg_cpp_life_rem = 15

    if occursin("C2N", unit_type_data[:unit_type])
        C2N_start = floor(Int64, unit_type_data[:cpp_ret_lead]) + lag + 1
        C2N_end = (
            floor(Int64, unit_type_data[:cpp_ret_lead]) +
            lag +
            avg_cpp_life_rem
        )
        for y = C2N_start:C2N_end
            coal_data = filter(
                [:y, :unit_type] =>
                    ((t, unit_type) -> (t == y) && (unit_type == "Coal")),
                agg_econ_results,
            )
            if (size(coal_data)[1] != 0)
                unit_fs[y, :Revenue] = (
                    unit_fs[y, :Revenue] -
                    unit_type_data["num_cpp_rets"] * (
                        coal_data[1, :annualized_rev_per_unit] +
                        coal_data[1, :annualized_policy_adj_per_unit]
                    )
                )
            end
        end
    end

end

"""
    compute_total_generation(unit_type_data, unit_fs, lag; mode, orig_ret_pd)

Calculate the unit's total generation for the period, in kWh.
"""
function compute_total_generation(
    settings,
    current_pd,
    unit_type_data,
    unit_fs,
    lag,
    long_econ_results,
    wtd_hist_gens;
    mode,
    orig_ret_pd,
)
    # Helpful short variable names
    unit_d_x = convert(
        Int64,
        ceil(round(unit_type_data[:construction_duration], digits = 3)),
    )
    unit_op_life = unit_type_data[:unit_life]

    # In "new_xtr" mode, generation persists from end of construction to end
    #   of life
    # In "retirement" mode, generation occurs from the current period to the
    #   soonest of {the lag, the mandatory retirement period, the end of the
    #   forecast period}.
    if mode == "new_xtr"
        # Find the end of the capital expenditures period
        capex_end = 0
        for i = 1:size(unit_fs)[1]
            if round(unit_fs[i, :capex], digits = 3) > 0
                capex_end = i
            end
        end
        gen_start = (
            capex_end + 1 -
            ceil(Int64, round(unit_type_data[:rev_head_start], digits = 3))
        )

        gen_end = min(gen_start + unit_op_life, size(unit_fs)[1])
    elseif mode == "retirement"
        gen_start = 1
        # The maximum end of generation period is capped at the length of the
        #   unit_fs dataframe (as some units with unspecified retirement
        #   periods default to a retirement period of 9999).
        gen_end = min(orig_ret_pd, size(unit_fs)[1])
    end

    # If the project is a C2N project, initialize a column to record lost coal
    #   generation
    if occursin("C2N", unit_type_data[:unit_type])
        unit_fs[!, :coal_gen] .= 0.0
    end

    transform!(
        long_econ_results,
        [:gen, :Probability, :num_units] =>
            ((gen, prob, num_units) -> gen .* prob .* 365 ./ num_units) =>
                :annualized_gen_per_unit,
    )
    agg_econ_results = combine(
        groupby(long_econ_results, [:y, :unit_type]),
        :annualized_gen_per_unit => sum,
        renamecols = false,
    )

    if wtd_hist_gens == nothing
        hist_gen = 0
        hist_wt = 0
    else
        try
            hist_gen = wtd_hist_gens[1, unit_type_data[:unit_type]]
            hist_wt = settings["dispatch"]["hist_wt"]
        catch
            hist_gen = 0
            hist_wt = settings["dispatch"]["hist_wt"]
        end
    end

    # Distribute generation values time series
    for y = gen_start:gen_end
        row = filter(
            [:y, :unit_type] =>
                (t, unit_type) -> (t == y) && (unit_type == unit_type_data[:unit_type]),
            agg_econ_results,
        )

        if size(row)[1] != 0
            if hist_wt == 0
                wt = 0
            else
                wt = hist_wt^(y - current_pd)
            end
            unit_fs[y, :gen] = (
                (1 - hist_wt) * row[1, :annualized_gen_per_unit] +
                hist_wt * hist_gen
            )
        else
            unit_fs[y, :gen] = 0
        end
    end

    #coal_avg_rem_life = 15

    #if occursin("C2N", unit_type_data[:unit_type])
    #    for y = (floor(Int64, unit_type_data[:cpp_ret_lead])+lag+1):(floor(Int64, unit_type_data[:cpp_ret_lead])+lag+coal_avg_rem_life)
    #        coal_data = filter([:y, :unit_type] => (t, unit_type) -> (t == y) && (unit_type == "coal"), agg_econ_results)
    #        if (size(coal_data)[1] != 0)
    #            unit_fs[y, :coal_gen] = unit_type_data["num_cpp_rets"] * (coal_data[1, :annualized_gen_per_unit])
    #        end
    #    end
    #end


end


"""
    forecast_unit_op_costs(unit_type_data, unit_fs, lag; mode, orig_ret_pd)

Forecast cost line items for the current unit:
 - fuel cost
 - VOM
 - FOM
"""
function forecast_unit_op_costs(
    settings,
    unit_type_data,
    unit_fs,
    lag;
    mode = "new_xtr",
    orig_ret_pd,
)
    # Helpful short variable names
    unit_d_x = convert(
        Int64,
        ceil(round(unit_type_data[:construction_duration], digits = 3)),
    )
    unit_op_life = unit_type_data[:unit_life]
    rev_head_start =
        ceil(Int64, round(unit_type_data[:rev_head_start], digits = 3))
    cpp_ret_lead =
        floor(Int64, round(unit_type_data[:cpp_ret_lead], digits = 3))

    if !occursin("C2N", unit_type_data[:unit_type])
        # Compute total fuel cost
        # Unit conversions:
        #  gen [MWh/year] * FC_per_MWh [$/MWh] = $/year
        transform!(
            unit_fs,
            [:gen] =>
                ((gen) -> gen .* unit_type_data[:FC_per_MWh]) => :Fuel_Cost,
        )

        # Compute total VOM cost incurred during generation
        # Unit conversions:
        #   gen [MWh/year] * VOM [$/MWh] = $/year
        transform!(
            unit_fs,
            [:gen] => ((gen) -> gen .* unit_type_data[:VOM]) => :VOM_Cost,
        )
    else
        transform!(
            unit_fs,
            [:gen, :coal_gen] =>
                (
                    (gen, coal_gen) ->
                        gen .* unit_type_data[:FC_per_MWh] .-
                        coal_gen .* 143.07
                ) => :Fuel_Cost,
        )
        transform!(
            unit_fs,
            [:gen, :coal_gen] =>
                (
                    (gen, coal_gen) ->
                        gen .* unit_type_data[:VOM] .- coal_gen .* 4.4
                ) => :VOM_Cost,
        )
    end

    # Compute total FOM cost for each period
    # Unit conversions:
    #   FOM [$/kW-year] * capacity [MW] * [1000 kW / 1 MW] = $/year
    if mode == "new_xtr"
        # Find the end of the capital expenditures period
        capex_end =
            maximum(filter(:capex => capex -> capex > 0, unit_fs)[!, :year])
        pre_zeros = zeros(capex_end - rev_head_start)
        op_ones = ones(unit_op_life)
    elseif mode == "retirement"
        pre_zeros = zeros(0)
        op_ones = ones(min(size(unit_fs)[1], orig_ret_pd))
    end
    post_zeros = zeros(size(unit_fs)[1] - size(pre_zeros)[1] - size(op_ones)[1])
    unit_fs[!, :FOM_Cost] = vcat(pre_zeros, op_ones, post_zeros)

    unit_fs[!, :FOM_Cost] .= (
        unit_fs[!, :FOM_Cost] .* unit_type_data[:FOM] .*
        unit_type_data[:capacity] .* settings["constants"]["MW2kW"]
    )
    if occursin("C2N", unit_type_data[:unit_type])
        C2N_start = lag + 1 + cpp_ret_lead
        C2N_end = C2N_start + unit_op_life
        for y = C2N_start:C2N_end
            unit_fs[y, :FOM_Cost] = (
                unit_fs[y, :FOM_Cost] - (
                    unit_type_data["num_cpp_rets"] *
                    500 *
                    settings["constants"]["MW2kW"] *
                    39.70
                )
            )
        end
    end

end


"""
    propagate_accounting_line_items(unit_fs, db)

Compute out and save all accounting line items:
 - EBITDA
 - EBIT
 - EBT
 - taxes owed
 - Net Income
 - Free Cash Flow
"""
function propagate_accounting_line_items(unit_fs, db)
    # Compute EBITDA
    transform!(
        unit_fs,
        [:Revenue, :Fuel_Cost, :VOM_Cost, :FOM_Cost] =>
            ((rev, fc, VOM, FOM) -> rev - fc - VOM - FOM) => :EBITDA,
    )

    # Compute EBIT
    transform!(
        unit_fs,
        [:EBITDA, :depreciation] => ((EBITDA, dep) -> EBITDA - dep) => :EBIT,
    )

    # Compute EBT
    transform!(
        unit_fs,
        [:EBIT, :interest_payment] =>
            ((EBIT, interest) -> EBIT - interest) => :EBT,
    )

    # Retrieve the system corporate tax rate from the database
    # Extract the value into a temporary dataframe
    tax_rate =
        DBInterface.execute(
            db,
            string(
                "SELECT value FROM model_params ",
                "WHERE parameter == 'tax_rate'",
            ),
        ) |> DataFrame
    # Pull out the bare value
    tax_rate = tax_rate[1, :value]

    # Compute taxes owed
    transform!(unit_fs, [:EBT] => ((EBT) -> EBT .* tax_rate) => :Tax_Owed)

    # Compute net income
    transform!(
        unit_fs,
        [:EBT, :Tax_Owed] => ((EBT, tax_owed) -> EBT - tax_owed) => :Net_Income,
    )

    # Compute free cash flow (FCF)
    transform!(
        unit_fs,
        [:Net_Income, :interest_payment, :capex] =>
            ((NI, interest, capex) -> NI + interest - capex) => :FCF,
    )

end


function compute_FS_delta_value(unit_fs, lag, db)
    # Make a copy of the baseline fs to use for the early-retirement case
    early_ret_fs = deepcopy(unit_fs)

    # Zero out all values which should be zero due to early retirement
    columns_to_adjust = [:gen, :Revenue, :VOM_Cost, :Fuel_Cost, :FOM_Cost]
    for col in columns_to_adjust
        # Set values to 0, starting at the end of the lag period
        early_ret_fs[(lag + 1):size(unit_fs)[1], col] .= 0
    end

    # Re-propagate the accounting logic for the early-retirement sheet
    propagate_accounting_line_items(early_ret_fs, db)

    # Create a final dataframe for the delta between these two world-states
    final_fs = deepcopy(unit_fs)

    # Columns to adjust: exclude year and discount factor
    cols_to_diff = (col for col in names(unit_fs))
    for col in cols_to_diff
        final_fs[:, col] = early_ret_fs[!, col] - unit_fs[!, col]
    end

    final_fs[!, :year] = unit_fs[!, :year]

    return final_fs

end


"""
    compute_alternative_NPV(unit_fs, agent_params)

For this project alternative (unit type + lag time), compute the project's FCF NPV.
"""
function compute_alternative_NPV(unit_fs, agent_params)
    # Discount rate is WACC
    d = (
        agent_params[1, :debt_fraction] * agent_params[1, :cost_of_debt] +
        (
            (1 - agent_params[1, :debt_fraction]) *
            agent_params[1, :cost_of_equity]
        )
    )

    # Add a column of compounded discount factors to the dataframe
    transform!(
        unit_fs,
        [:year] =>
            ((year) -> (1 + d) .^ (-1 .* (year .- 1))) => :discount_factor,
    )

    # Discount the alternative's FCF NPV
    FCF_NPV = transpose(unit_fs[!, :FCF]) * unit_fs[!, :discount_factor]

    return FCF_NPV, unit_fs

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
    PA_uids,
    PA_fs_dict,
    total_demand,
    asset_counts,
    agent_params,
    unit_specs,
    current_pd,
    system_portfolios,
    db,
    agent_id,
    agent_fs,
    fc_pd,
)
    # Create the model object
    @debug "Setting up model..."

    # Initialize a model with the settings-specified optimizer
    m = create_model_with_optimizer(settings)

    # Parameter names
    num_alternatives = size(PA_uids)[1]
    num_time_periods = size(PA_fs_dict[PA_uids[1, :uid]])[1]

    # Set up variables
    # Number of units of each type to build: must be Integer
    @variable(m, u[1:num_alternatives] >= 0, Int)

    # Compute expected marginal generation and effective nameplate capacity
    #   contribution per alternative type
    marg_gen = zeros(num_alternatives, num_time_periods)
    marg_eff_cap = zeros(num_alternatives, num_time_periods)
    for i = 1:size(PA_uids)[1]
        for j = 1:num_time_periods
            # Marginal generation in kWh
            marg_gen[i, j] = (
                PA_fs_dict[PA_uids[i, :uid]][j, :gen] /
                settings["constants"]["MW2kW"]
            )
            # Convert total anticipated marginal generation to an effective
            #   nameplate capacity and save to the appropriate entry in the
            #   marginal generation array
            unit_type_data = filter(
                :unit_type => x -> x == PA_uids[i, :unit_type],
                unit_specs,
            )

            if PA_uids[i, :project_type] == "new_xtr"
                op_start =
                    PA_uids[i, :lag] + unit_type_data[1, :construction_duration]
                if j > op_start
                    marg_eff_cap[i, j] = (
                        unit_type_data[1, :capacity] *
                        unit_type_data[1, :capacity_factor]
                    )
                end
            elseif PA_uids[i, :project_type] == "retirement"
                if (j > PA_uids[i, :lag]) && (j <= PA_uids[i, :ret_pd])
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
    PA_uids[!, :current] .= 0
    for i = 1:size(PA_uids)[1]
        if PA_uids[i, :project_type] == "retirement"
            PA_uids[i, :current] = 1
        elseif (
            (PA_uids[i, :project_type] == "new_xtr") && (PA_uids[i, :lag] == 0)
        )
            PA_uids[i, :current] = 1
        end
    end

    # Prevent the agent from intentionally causing foreseeable energy shortages
    for i = 1:settings["agent_opt"]["shortage_protection_period"]
        # Get extant capacity sufficiency measures (projected demand and
        #   capacity)
        sufficiency =
            filter(:period => x -> x == current_pd + i - 1, total_demand)
        pd_total_demand = sufficiency[1, :total_demand]
        total_eff_cap = sufficiency[1, :total_eff_cap]

        # Enforce different capacity assurance requirements based on extant
        #   reserve margin
        if (
            total_eff_cap / pd_total_demand >
            settings["agent_opt"]["cap_decrease_threshold"]
        )
            margin = settings["agent_opt"]["cap_decrease_margin"]
        elseif (
            total_eff_cap / pd_total_demand >
            settings["agent_opt"]["cap_maintain_threshold"]
        )
            margin = settings["agent_opt"]["cap_maintain_margin"]
        else
            margin = settings["agent_opt"]["cap_increase_margin"]
        end
        @constraint(
            m,
            transpose(u .* PA_uids[:, :current]) * marg_eff_cap[:, i] >=
            (total_eff_cap - pd_total_demand) * margin
        )

    end

    # Create arrays of expected marginal debt, interest, dividends, and FCF
    #   per unit type
    marg_debt = zeros(num_alternatives, num_time_periods)
    marg_int = zeros(num_alternatives, num_time_periods)
    marg_div = zeros(num_alternatives, num_time_periods)
    marg_FCF = zeros(num_alternatives, num_time_periods)
    for i = 1:size(PA_uids)[1]
        project = PA_fs_dict[PA_uids[i, :uid]]
        for j = 1:num_time_periods
            # Retrieve the marginal value of interest
            # Scale to units of $B
            marg_debt[i, j] = project[j, :remaining_debt_principal] / 1e9
            marg_int[i, j] = project[j, :interest_payment] / 1e9
            marg_div[i, j] = (
                project[j, :remaining_debt_principal] *
                agent_params[1, :cost_of_equity] / 1e9
            )
            marg_FCF[i, j] = project[j, :FCF] / 1e9
        end
    end


    # Prevent the agent from reducing its credit metrics below Moody's Baa
    #   rating thresholds (from the Unregulated Power Companies ratings grid)
    for i = 1:6
        @constraint(
            m,
            (
                agent_fs[i, :FCF] / 1e9 +
                sum(u .* marg_FCF[:, i]) +
                (1 - 4.2) * (
                    agent_fs[i, :interest_payment] / 1e9 +
                    sum(u .* marg_int[:, i])
                )
            ) >= 0
        )
    end

    # Enforce the user-specified maximum number of new construction/retirement
    #   projects by type per period, and the :allowed field in PA_uids
    for i = 1:num_alternatives
        if PA_uids[i, :project_type] == "new_xtr"
            @constraint(
                m,
                u[i] .<=
                convert(Int64, PA_uids[i, :allowed]) .*
                settings["agent_opt"]["max_type_newbuilds_per_pd"]
            )
        elseif PA_uids[i, :project_type] == "retirement"
            unit_type = PA_uids[i, :unit_type]
            ret_pd = PA_uids[i, :ret_pd]
            asset_count = filter(
                [:unit_type, :retirement_pd] =>
                    (x, y) -> x == unit_type && y == ret_pd,
                asset_counts,
            )[
                1,
                :count,
            ]
            max_retirement = (
                convert(Int64, PA_uids[i, :allowed]) .* min(
                    asset_count,
                    settings["agent_opt"]["max_type_rets_per_pd"],
                )
            )
            #            @constraint(m, u[i] .<= max_retirement)
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
    ret_summation_matrix = zeros(size(retireable_asset_counts)[1], size(PA_uids)[1])
    for i = 1:size(retireable_asset_counts)[1]
        for j = 1:size(PA_uids)[1]
            if (
                (PA_uids[j, :project_type] == "retirement") &
                (
                    PA_uids[j, :unit_type] ==
                    retireable_asset_counts[i, :unit_type]
                ) &
                (
                    PA_uids[j, :ret_pd] ==
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

    # Set up the coal retirements of matrix, to ensure that the total "pool"
    #   of coal plants available for retirement is respected across the entire
    #   visible horizon
    coal_retirements = zeros(size(PA_uids)[1], 10)
    planned_coal_units_operating = DataFrame(pd = Int64[], num_units = Int64[])

    for i = 1:size(PA_uids)[1]
        # Mark all of the direct coal retirements in the matrix
        if (PA_uids[i, :project_type] == "retirement") && (PA_uids[i, :unit_type] == "coal")
            coal_retirements[i, PA_uids[i, :lag]+1] = 1
        end

        # Mark all of the C2N-forced coal retirements in the matrix
        if occursin("C2N", PA_uids[i, :unit_type])
            unit_type_data = filter(:unit_type => ((ut) -> ut == PA_uids[i, :unit_type]), unit_specs)
            target_coal_ret_pd = PA_uids[i, :lag] + unit_type_data[1, :cpp_ret_lead] + 1
            if target_coal_ret_pd <= size(coal_retirements)[2]
                coal_retirements[i, target_coal_ret_pd] = 1
            end
        end

    end

    for i=1:size(coal_retirements)[2]
        num_units = 0

        # Get number of planned coal units operating during this period
        coal_ops = DBInterface.execute(db, "SELECT unit_type, COUNT(unit_type) FROM assets WHERE unit_type = 'coal' AND C2N_reserved = 0 AND agent_id = $agent_id AND completion_pd <= $i AND retirement_pd >= $i AND cancellation_pd >= $i GROUP BY unit_type") |> DataFrame
        if size(coal_ops)[1] > 0
            num_units = coal_ops[1, "COUNT(unit_type)"]
        end

        # Record the number of operating coal units in this period
        push!(planned_coal_units_operating, (i, num_units))
    end

    # Constrain rolling total of coal retirements
    for i=1:size(coal_retirements)[2]
        @constraint(m, sum(transpose(u) * coal_retirements[:, i]) <= planned_coal_units_operating[i, :num_units])
        @constraint(m, sum(sum(transpose(u) * coal_retirements[:, j]) for j = 1:i) <= planned_coal_units_operating[1, :num_units])
    end

    # Create the objective function 
    profit_lamda = settings["agent_opt"]["profit_lamda"] / 1e9
    credit_rating_lamda = settings["agent_opt"]["credit_rating_lamda"]
    cr_horizon = settings["agent_opt"]["cr_horizon"]
    int_bound = settings["agent_opt"]["int_bound"]

    @objective(
        m,
        Max,
        (
            profit_lamda * (transpose(u) * PA_uids[!, :NPV]) +
            credit_rating_lamda * (
                sum(agent_fs[1:cr_horizon, :FCF]) / 1e9 +
                sum(transpose(u) * marg_FCF[:, 1:cr_horizon]) +
                sum(agent_fs[1:cr_horizon, :interest_payment]) / 1e9 +
                sum(transpose(u) * marg_int[:, 1:cr_horizon]) -
                (int_bound) * (
                    sum(agent_fs[1:cr_horizon, :interest_payment]) / 1e9 +
                    sum(transpose(u) * marg_int[:, 1:cr_horizon])
                )
            )
        )
    )


    @debug "Optimization model set up."

    return m
end


### Postprocessing
function finalize_results_dataframe(m, PA_uids)
    # Check solve status of model
    status = string(termination_status.(m))

    # If the model solved to optimality, convert the results to type Int64
    if status == "OPTIMAL"
        unit_qty = Int64.(round.(value.(m[:u])))
    else
        # If the model did not solve to optimality, the agent does nothing. Return
        #   a vector of all zeroes instead.
        unit_qty = zeros(Int64, size(PA_uids)[1])
    end

    all_results = hcat(PA_uids, DataFrame(units_to_execute = unit_qty))

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
    for i = 1:size(all_results)[1]
        # Retrieve the individual result row for convenience
        result = all_results[i, :]

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

    # Save the agent's decisions to the database
    save_agent_decisions(db, agent_id, all_results)


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
                    "AND retirement_pd = ? AND agent_id = ?",
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
        new_ret_pd = (
            current_pd +
            result[:lag] +
            convert(
                Int64,
                ceil(round(unit_type_specs[:cpp_ret_lead], digits = 3)),
            )
        )
        C2N_reserved = 0

        # Generate a list of coal units which will still be operational by
        #   the necessary period
        match_vals = (new_ret_pd, new_ret_pd, agent_id, C2N_reserved)
        ret_candidates =
            DBInterface.execute(
                db,
                string(
                    "SELECT asset_id FROM assets WHERE unit_type = 'coal' ",
                    "AND completion_pd <= ? AND retirement_pd > ? ",
                    "AND agent_id = ? AND C2N_reserved = ?",
                ),
                match_vals,
            ) |> DataFrame

        # Set the number of units to execute
        units_to_execute = unit_type_specs["num_cpp_rets"] * result[:units_to_execute]

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

        # Overwrite the original record's retirement period with the current
        #   period
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
    long_econ_results,
)
    # Retrieve horizontally-abbreviated dataframes
    short_econ_results = select(
        long_econ_results,
        [
            :unit_type,
            :y,
            :d,
            :h,
            :gen,
            :annualized_rev_per_unit,
            :annualized_VOM_per_unit,
            :annualized_FC_per_unit,
            :annualized_policy_adj_per_unit,
        ],
    )
    short_unit_specs = select(unit_specs, [:unit_type, :FOM])

    # Retrieve the agent's complete portfolio forecast
    agent_portfolio_forecast =
        get_agent_portfolio_forecast(agent_id, db, current_pd, fc_pd)

    # Inner join the year's portfolio with financial pivot
    fin_results = innerjoin(
        short_econ_results,
        agent_portfolio_forecast,
        on = [:y, :unit_type],
    )

    # Fill in total revenue
    transform!(
        fin_results,
        [
            :annualized_rev_per_unit,
            :annualized_policy_adj_per_unit,
            :num_units,
        ] => ((rev, adj, num_units) -> (rev .+ adj) .* num_units) => :total_rev,
    )

    # Fill in total VOM
    transform!(
        fin_results,
        [:annualized_VOM_per_unit, :num_units] =>
            ((VOM, num_units) -> VOM .* num_units) => :total_VOM,
    )

    # Fill in total fuel costs
    transform!(
        fin_results,
        [:annualized_FC_per_unit, :num_units] =>
            ((FC, num_units) -> FC .* num_units) => :total_FC,
    )

    # Create the annualized results dataframe so far
    results_pivot = combine(
        groupby(fin_results, :y),
        [:total_rev, :total_VOM, :total_FC] .=> sum;
        renamecols = false,
    )

    # If the agent is projected to reach zero installed capacity within the 
    #   forecast period, pad the results_pivot with zeros until the
    #   end of the forecast period
    if size(results_pivot)[1] < fc_pd
        y_start = maximum(results_pivot[!, :y]) + 1
        pad_length = fc_pd - size(results_pivot)[1]

        pad_df = DataFrame(
            y = y_start:(y_start + pad_length - 1),
            total_rev = zeros(Int64, pad_length),
            total_VOM = zeros(Int64, pad_length),
            total_FC = zeros(Int64, pad_length),
        )

        results_pivot = vcat(results_pivot, pad_df)
    end

    fs = DataFrame(
        base_pd = ones(Int64, fc_pd) .* current_pd,
        projected_pd = current_pd:(current_pd + fc_pd - 1),
        revenue = results_pivot[!, :total_rev],
        VOM = results_pivot[!, :total_VOM],
        fuel_costs = results_pivot[!, :total_FC],
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

    FOM_df = select(
        combine(groupby(FOM_df, :y), :total_FOM => sum; renamecols = false),
        [:y, :total_FOM],
    )
    rename!(FOM_df, :y => :projected_pd, :total_FOM => :FOM)

    fs = innerjoin(fs, FOM_df, on = :projected_pd)

    # Fill in total depreciation
    total_pd_depreciation =
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
    total_pd_interest =
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

    # Fill in total capex
    total_pd_capex =
        DBInterface.execute(
            db,
            string(
                "SELECT projected_pd, SUM(capex) ",
                "FROM capex_projections ",
                "WHERE agent_id = $agent_id ",
                "AND base_pd = $current_pd GROUP BY projected_pd",
            ),
        ) |> DataFrame

    # Join the scheduled columns to the financial statement
    fs = leftjoin(fs, total_pd_depreciation, on = :projected_pd)
    fs = leftjoin(fs, total_pd_interest, on = :projected_pd)
    fs = leftjoin(fs, total_pd_capex, on = :projected_pd)

    # Standardize column names
    rename!(
        fs,
        Symbol("SUM(depreciation)") => :depreciation,
        Symbol("SUM(interest_payment)") => :interest_payment,
        Symbol("SUM(capex)") => :capex,
    )

    # Replace missing values with 0s
    fs[!, :depreciation] = coalesce.(fs[!, :depreciation], 0)
    fs[!, :interest_payment] = coalesce.(fs[!, :interest_payment], 0)
    fs[!, :capex] = coalesce.(fs[!, :capex], 0)


    ### Computed FS quantities
    # EBITDA
    transform!(
        fs,
        [:revenue, :VOM, :FOM, :fuel_costs] =>
            ((rev, VOM, FOM, FC) -> (rev - VOM - FOM - FC)) => :EBITDA,
    )

    # EBIT
    transform!(
        fs,
        [:EBITDA, :depreciation] => ((EBITDA, dep) -> (EBITDA - dep)) => :EBIT,
    )

    # EBT
    transform!(
        fs,
        [:EBIT, :interest_payment] =>
            ((EBIT, interest) -> EBIT - interest) => :EBT,
    )

    # Retrieve the system corporate tax rate from the database
    # Extract the value into a temporary dataframe
    tax_rate =
        DBInterface.execute(
            db,
            string(
                "SELECT value FROM model_params ",
                "WHERE parameter == 'tax_rate'",
            ),
        ) |> DataFrame
    # Pull out the bare value
    tax_rate = tax_rate[1, :value]

    # Compute actual tax paid
    transform!(fs, :EBT => ((EBT) -> EBT * tax_rate) => :tax_paid)

    # Net Income
    transform!(
        fs,
        [:EBT, :tax_paid] => ((EBT, tax) -> (EBT - tax)) => :Net_Income,
    )

    # Free Cash Flow
    transform!(
        fs,
        [:Net_Income, :depreciation, :capex] =>
            ((NI, dep, capex) -> NI + dep - capex) => :FCF,
    )

    # Save the dataframe to the database
    save_agent_fs!(fs, agent_id, db)

    return fs


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

    for row in Tuple.(eachrow(fs))
        DBInterface.execute(
            db,
            string(
                "INSERT INTO agent_financial_statements VALUES ",
                "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ),
            row,
        )
    end

end

function save_agent_decisions(db, agent_id, decision_df)
    decision_df[!, :agent_id] .= agent_id
    cols_to_ignore = [:uid, :current]
    select!(decision_df, :agent_id, Not(vcat([:agent_id], cols_to_ignore)))
    for row in Tuple.(eachrow(decision_df))
        DBInterface.execute(
            db,
            string(
                "INSERT INTO agent_decisions ",
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            ),
            row,
        )
    end
end


function display_agent_choice_results(CLI_args, all_results)
    if CLI_args["verbosity"] == 2
        @info "Project alternatives to execute:"
        @info filter(:units_to_execute => u -> u > 0, all_results)
    elseif CLI_args["verbosity"] == 3
        @debug "Alternatives to execute:"
        @debug all_results
    end

end

end
