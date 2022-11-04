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

using SQLite, DataFrames, CSV, JuMP, GLPK, Cbc, Logging, Tables, SCIP, HiGHS

try
    using CPLEX
catch LoadError
    @info string("CPLEX is not available!")
end

include("./dispatch.jl")
using .Dispatch
include("./C2N_projects.jl")
using .C2N

export ProjectAlternative, load_db, get_current_period, set_up_local_paths, get_agent_id, get_agent_params,load_unit_type_data, set_forecast_period, extrapolate_demand, project_demand_flat, project_demand_exponential, allocate_fuel_costs, create_FS_dict, get_unit_specs, get_table, show_table, get_WIP_projects_list, get_demand_forecast, get_net_demand, get_next_asset_id, ensure_projects_not_empty, authorize_anpe, generate_capex_profile, set_initial_debt_principal_series, generate_prime_movers, forecast_unit_revenue_and_gen, forecast_unit_op_costs, propagate_accounting_line_items, compute_alternative_NPV,set_up_model, get_current_assets_list, convert_to_marginal_delta_FS, postprocess_agent_decisions, set_up_project_alternatives, update_agent_financial_statement, save_agent_decisions

#####
# Constants
#####
MW2kW = 1000            # Conversion factor from MW to kW
MMBTU2BTU = 1000000     # Conversion factor from MMBTU to BTu
kW2W = 1000             # Conversion factor from kW to W
hours_per_year = 8760   # Number of hours in a year (without final 0.25 day)


#####
# Setup functions
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
        command = string("SELECT * FROM agent_params WHERE agent_id = ", agent_id)
        df = DBInterface.execute(db, command) |> DataFrame
    catch e
        @error "Could not get agent parameters from file:"
        @error e
        exit()
    end
end


function set_up_local_paths(settings)
    settings["file_paths"]["ABCE_abs_path"] = @__DIR__
    if settings["simulation"]["run_ALEAF"] == true
        try
            settings["file_paths"]["ALEAF_abs_path"] = ENV["ALEAF_DIR"]            
        catch LoadError
            println("The environment variable ALEAF_abs_path does not appear to be set. Please make sure it points to the correct directory.")
        end
    else
        settings["file_paths"]["ALEAF_abs_path"] = "NULL_PATH"
    end

    return settings

end


function load_unit_type_data(unit_data_file)
    unit_data = CSV.read(unit_data_file, DataFrame)
    num_types = size(unit_data)[1]
    return unit_data, num_types
end


function set_forecast_period(unit_specs, num_lags)
    transform!(unit_specs, [:d_x, :unit_life] => ((lead_time, unit_life) -> lead_time + unit_life) => :full_life)
    max_horizon = convert(Int64, ceil(maximum(unit_specs[!, :full_life]) + num_lags))
    return max_horizon
end


function allocate_fuel_costs(unit_data, fuel_costs)
    num_units = size(unit_data)[1]
    unit_data[!, :uc_fuel] = zeros(num_units)
    for i = 1:num_units
        unit_data[i, :uc_fuel] = fuel_costs[fuel_costs[!, :fuel_type] == unit_data[i, :fuel_type], :cost_per_mmbtu]
    end
    return unit_data
end




#####
# Database interaction functions
#####


function get_table(db, table_name)
    command = string("SELECT * FROM ", string(table_name))
    df = DBInterface.execute(db, command) |> DataFrame
    return df
end


function show_table(db, table_name)
    command = string("SELECT * FROM ", string(table_name))
    df = DBInterface.execute(db, command) |> DataFrame
    # @info string("\nTable \'", table_name, "\':")
    # @info df
    return df
end


function get_WIP_projects_list(db, pd, agent_id)
    # Get a list of all WIP (non-complete, non-cancelled) projects for the given agent
    # Time-period convention:
    #   <status>_pd == "the first period where this asset has status <status>"
    SQL_get_proj = SQLite.Stmt(db, string("SELECT asset_id FROM assets WHERE agent_id = ", agent_id, " AND completion_pd > ", pd, " AND cancellation_pd > ", pd))
    project_list = DBInterface.execute(SQL_get_proj) |> DataFrame
    return project_list
end


function get_current_assets_list(db, pd, agent_id)
    # Get a list of all of the agent's currently-operating assets, plus their
    #   unit type and mandatory retirement date
    # Time-period convention:
    #   <status>_pd == "the first period where this asset has status <status>"
    SQL_get_assets = SQLite.Stmt(db, string("SELECT asset_id, unit_type, retirement_pd, C2N_reserved FROM assets WHERE agent_id = ", agent_id, " AND completion_pd <= ", pd, " AND cancellation_pd > ", pd, " AND retirement_pd > ", pd))
    asset_list = DBInterface.execute(SQL_get_assets) |> DataFrame

    # Count the number of assets by type
    asset_counts = combine(groupby(asset_list, [:unit_type, :retirement_pd, :C2N_reserved]), nrow => :count)
    # @info asset_counts

    return asset_list, asset_counts
end


function extrapolate_demand(visible_demand, db, pd, fc_pd, settings)
    mode = settings["demand"]["demand_projection_mode"]
    if mode == "flat"
        demand = project_demand_flat(visible_demand, fc_pd)
    elseif mode == "exp_termrate"
        term_demand_gr = settings["demand"]["terminal_demand_growth_rate"]
        demand = project_demand_exp_termrate(visible_demand, fc_pd, term_demand_gr)
    elseif mode == "exp_fitted"
        demand = project_demand_exp_fitted(visible_demand, db, pd, fc_pd, settings)
    else
        @error string("The specified demand extrapolation mode, ", mode, ", is not implemented.")
        @error "Please use 'flat', 'exp_termrate', or 'exp_fitted' at this time."
        @error "Terminating..."
        exit()
    end

    return demand
end


function project_demand_flat(visible_demand, fc_pd)
    demand = DataFrame(demand = zeros(Float64, convert(Int64, fc_pd)))
    demand[1:size(visible_demand)[1], :total_demand] .= visible_demand[!, :total_demand]
    demand[(size(visible_demand)[1] + 1):fc_pd, :total_demand] .= demand[size(visible_demand)[1], :total_demand] 
    return demand
end


function project_demand_exp_termrate(visible_demand, fc_pd, term_demand_gr)
    total_demand = deepcopy(visible_demand)
    prev_pd = visible_demand[size(visible_demand)[1], :period]
    prev_real_demand = last(visible_demand[!, :real_demand])
    for i = prev_pd+1:fc_pd
        next_real_demand = prev_real_demand * (1 + term_demand_gr) ^ (i - prev_pd - 1)
        push!(total_demand, [i, next_real_demand])
    end
    return total_demand
end


function project_demand_exp_fitted(visible_demand, db, pd, fc_pd, settings)
    # Retrieve historical data, if available
    demand_projection_window = settings["demand"]["demand_projection_window"]
    demand_history_start = max(0, pd - settings["demand"]["demand_visibility_horizon"])
    start_and_end = (demand_history_start, pd)
    demand_history = DBInterface.execute(db, "SELECT demand FROM demand WHERE period >= ? AND period < ? ORDER BY period ASC", start_and_end) |> DataFrame

    # Create suitable arrays for x (including intercept) and y
    num_obs = size(visible_demand)[1] + size(demand_history)[1]
    x = hcat(ones(num_obs), [i for i=0:num_obs-1])
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
    beta[2] = (beta[2] * size(y)[1] + settings["demand"]["historical_demand_growth_rate"] * (demand_projection_window - size(y)[1])) / demand_projection_window

    # Project future demand
    proj_horiz = fc_pd - size(visible_demand)[1]
    x_proj = hcat(ones(proj_horiz), [i for i=size(visible_demand)[1]+size(demand_history)[1]+1:fc_pd+size(demand_history)[1]])
    y_log_proj = x_proj[:, 1] .* beta[1] + x_proj[:, 2] .* beta[2]
    y_proj = exp.(x_proj[:, 1] .* beta[1] + x_proj[:, 2] .* beta[2])

    all_proj = DataFrame(intercept = x_proj[:, 1], x_val = x_proj[:, 2], y_log_proj = y_log_proj, y_val = y_proj)

    # Concatenate input and projected demand data into a single DataFrame
    demand = DataFrame(demand = vcat(visible_demand[!, :demand], y_proj))
    return demand
end


function get_demand_forecast(db, pd, agent_id, fc_pd, settings)
    # Retrieve 
    vals = (pd, pd + settings["demand"]["demand_visibility_horizon"])
    visible_demand = DBInterface.execute(db, "SELECT period, demand FROM demand WHERE period >= ? AND period < ?", vals) |> DataFrame

    rename!(visible_demand, :demand => :real_demand)

    demand_forecast = extrapolate_demand(visible_demand, db, pd, fc_pd, settings)
    prm = DBInterface.execute(db, "SELECT * FROM model_params WHERE parameter = 'PRM'") |> DataFrame

    transform!(demand_forecast, :real_demand => ((dem) -> dem .* (1 + prm[1, :value]) .+ settings["system"]["peak_initial_reserves"]) => :total_demand)

    return demand_forecast
end


function get_net_demand(db, pd, agent_id, fc_pd, demand_forecast, all_year_system_portfolios, unit_specs)
    # Calculate the amount of forecasted net demand in future periods
    installed_cap_forecast = DataFrame(period = Int64[], derated_capacity = Float64[])
    vals = (pd, pd)

    total_caps = DataFrame(period = Int64[], total_eff_cap = Float64[])

    for i = pd:pd+fc_pd-1
        year_portfolio = innerjoin(all_year_system_portfolios[i], unit_specs, on = :unit_type)
        transform!(year_portfolio, [:num_units, :capacity, :CF] => ((num_units, cap, CF) -> num_units .* cap .* CF) => :effective_capacity)
        total_year_cap = sum(year_portfolio[!, :effective_capacity])
        push!(total_caps, [i, total_year_cap])
    end

    demand_forecast = innerjoin(demand_forecast, total_caps, on = :period)
    transform!(demand_forecast, [:total_demand, :total_eff_cap] => ((demand, eff_cap) -> demand - eff_cap) => :net_demand)

    return demand_forecast
end


function get_next_asset_id(db)
    tables_to_check = ["assets", "WIP_projects", "asset_updates", "WIP_updates", "depreciation_projections"]

    id_vals = Vector{Any}()

    for table in tables_to_check
        id_val = DBInterface.execute(db, "SELECT MAX(asset_id) FROM $table") |> DataFrame
        id_val = id_val[1, Symbol("MAX(asset_id)")]
        push!(id_vals, id_val)
    end

    # Convert id_vals into a skipmissing object
    id_vals = skipmissing(id_vals)

    # Return the next available asset ID (one greater than the current largest ID)
    next_id = maximum(id_vals) + 1

    return next_id    
end


function get_unit_specs(db)
    # Retrieve the table of unit specifications from the DB
    df = DBInterface.execute(db, "SELECT * FROM unit_specs") |> DataFrame
    num_types = size(df)[1]
    # Convert the Int-type columns to Int64
    df[!, :unit_life] = convert.(Int64, df[:, :unit_life])
    return df, num_types
end


function ensure_projects_not_empty(db, agent_id, project_list, current_period)
    # FAKE: only exists to ensure the Julia and Python scripts actually have something to do
    try
        if (size(project_list)[1] == 0) && (current_period <= 5)
            # Current period restriction is a testing spoof only to allow
            #    observation of post-completion model behavior

            # Create a new project
            new_asset_id = get_next_asset_id(db)

            # Assign some dummy data to the project
            xtr_vals = (new_asset_id, string(agent_id), 0, 1000, 10, 0)
            asset_vals = (new_asset_id, string(agent_id), "gas", "no", "no", 9999, 0)
            DBInterface.execute(db, "INSERT INTO WIP_projects VALUES (?, ?, ?, ?, ?, ?)", xtr_vals)
            DBInterface.execute(db, "INSERT INTO assets VALUES (?, ?, ?, ?, ?, ?, ?)", asset_vals)
            # @info string("Created project ", new_asset_id)

            # Update the list of WIP construction projects and return it
            project_list = get_WIP_projects_list(db, agent_id)
            return project_list
        else
            # If at least one construction project already exists, there's no
            #    need to do anything.
            return project_list
        end
    catch e
        @error "Could not insert a seed project into the agent's project list"
        @error e
        exit()
    end
end


function authorize_anpe(db, agent_id, current_period, project_list, unit_specs)
    # Loop through each project and authorize $100 of ANPE by setting the anpe value in WIP_projects
    for i = 1:size(project_list[!, :asset_id])[1]
        current_asset = project_list[i, :asset_id]
        asset_type = DBInterface.execute(db, string("SELECT unit_type FROM assets WHERE asset_id = ", current_asset)) |> DataFrame
        unit = filter(row -> row[:unit_type] == asset_type[1, :unit_type], unit_specs)
        # Authorize a uniform expenditure over the life of the project
        anpe_val = unit[1, :uc_x] * unit[1, :capacity] * 1000 / unit[1, :d_x]
#        anpe_val = 100000000   # $1B/period
        vals = (anpe_val, current_period, current_asset)
        DBInterface.execute(db, "UPDATE WIP_updates SET anpe = ? WHERE period = ? AND asset_id = ?", vals)
    end
end


#####
# NPV functions
#####


function set_up_project_alternatives(settings, unit_specs, asset_counts, num_lags, fc_pd, agent_params, db, current_pd, long_econ_results, C2N_specs)
    PA_uids = create_PA_unique_ids(settings, unit_specs, asset_counts, num_lags)

    PA_fs_dict = create_PA_pro_formas(PA_uids, fc_pd)

    PA_uids, PA_fs_dict = populate_PA_pro_formas(settings, PA_uids, PA_fs_dict, unit_specs, fc_pd, agent_params, db, current_pd, long_econ_results, C2N_specs)

    return PA_uids, PA_fs_dict

end

function create_PA_unique_ids(settings, unit_specs, asset_counts, num_lags)
    PA_uids = DataFrame(
                   unit_type = String[],
                   project_type = String[],
                   lag = Int64[],
                   ret_pd = Union{Int64, Nothing}[],
                   uid = Int64[],
                   NPV = Float64[],
                   allowed = Bool[]
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
        allowed=true
        if !(unit_type in settings["scenario"]["allowed_xtr_types"])
            allowed=false
        end

        for lag = 0:num_lags
            push!(PA_uids, [unit_type project_type lag nothing uid NPV allowed])
            uid += 1
        end
    end

    # Retirement
    project_type = "retirement"
    allowed = true   # by default, retirements are always allowed
    for unit_type in unique(asset_counts[!, :unit_type])
        for ret_pd in filter([:unit_type, :C2N_reserved] => ((x, reserved) -> (x == unit_type) && (reserved == 0)), asset_counts)[!, :retirement_pd]
            for lag = 0:num_lags
                push!(PA_uids, [unit_type project_type lag ret_pd uid NPV allowed])
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
                              gen = zeros(fc_pd)
                          )
    end

    return PA_fs_dict

end


function populate_PA_pro_formas(settings, PA_uids, PA_fs_dict, unit_specs, fc_pd, agent_params, db, current_pd, long_econ_results, C2N_specs)
    for uid in PA_uids[!, :uid]
        # Retrieve current project alternative definition for convenience
        current_PA = filter(:uid => x -> x == uid, PA_uids)[1, :]

        # Retrieve relevant unit type specs for convenience
        unit_type_data = filter(:unit_type => x -> x == current_PA[:unit_type], unit_specs)[1, :]

        # If this is a potential new construction project, compute series
        #   related to construction costs
        if current_PA[:project_type] == "new_xtr"
            # Generate this alternative's expected construction expenditure profile
            PA_fs_dict[uid][!, :capex] = generate_capex_profile(db, settings, current_pd, unit_type_data, current_PA[:lag], fc_pd, C2N_specs)

            # Set up the time-series of outstanding debt based on this
            #   construction expenditure profile
            set_initial_debt_principal_series(PA_fs_dict[uid], unit_type_data, current_PA[:lag], agent_params)

            # Generate "prime movers": debt payments and depreciation
            generate_prime_movers(unit_type_data, PA_fs_dict[uid], current_PA[:lag], agent_params[1, :cost_of_debt])
        end

        # Forecast unit revenue ($/period) and generation (kWh/period)
        forecast_unit_revenue_and_gen(
            unit_type_data,
            PA_fs_dict[uid],
            db,
            current_pd,
            current_PA[:lag],
            long_econ_results,
            mode = current_PA[:project_type],
            orig_ret_pd = current_PA[:ret_pd])

        # Forecast unit operating costs: VOM, fuel cost, and FOM
        forecast_unit_op_costs(unit_type_data, PA_fs_dict[uid], current_PA[:lag], mode=current_PA[:project_type], orig_ret_pd=current_PA[:ret_pd])

        # Propagate the accounting logic flow
        propagate_accounting_line_items(PA_fs_dict[uid], db)

        # If the current project is a retirement project, adjust the FS to show
        #   the difference between the original retirement plan and the 
        #   proposed outcome
        if current_PA[:project_type] == "retirement"
            PA_fs_dict[uid] = compute_FS_delta_value(PA_fs_dict[uid], current_PA[:lag], db)
        end

        FCF_NPV, PA_fs_dict[uid] = compute_alternative_NPV(PA_fs_dict[uid], agent_params)

        # save a representative example of each unit type to file (new_xtr only)
        savelag = 2
#        if (current_PA[:project_type] == "new_xtr") && (current_PA[:lag] == savelag)
#            ctype = current_PA[:unit_type]
#            fspath = joinpath(
#                         pwd(),
#                         "tmp",
#                         string(
#                             ctype,
#                             "_",
#                             savelag,
#                             "_fs.csv"
#                         )
#                     )
#            CSV.write(fspath, PA_fs_dict[uid])
#        end

        # Save the NPV result
        filter(:uid => x -> x == uid, PA_uids, view=true)[1, :NPV] = FCF_NPV

    end

    return PA_uids, PA_fs_dict

end


function check_valid_vector_mode(mode)
    if !(mode in ["new_xtr", "retirement"])
        # Invalid mode supplied; alert the user and exit
        @error string("Invalid decision vector type specified: ", mode)
        @error "Please ensure that 'mode' is set to either 'new_xtr' or 'retirement'."
        exit()
    end
end


function create_FS_dict(data, fc_pd, num_lags; mode="new_xtr")
    check_valid_vector_mode(mode)

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
            unit_name = string(data[i, :unit_type], "_",
                               project_type, "_",
                               ret_pd, "_lag-",
                               j)
            unit_FS = DataFrame(year = 1:fc_pd, capex = zeros(fc_pd), gen = zeros(fc_pd), remaining_debt_principal = zeros(fc_pd), debt_payment = zeros(fc_pd), interest_payment = zeros(fc_pd), depreciation = zeros(fc_pd))
            fs_dict[unit_name] = unit_FS
        end
    end
    return fs_dict
end


"""
    generate_xtr_cost_profile(unit_type_data, lag)

Creates a uniform expenditure profile for a potential construction project,
and appends appropriate numbers of leading and lagging zeros (no construction 
expenditures outside the construction period).

Arguments:
  unit_type_data: the appropriate unit specification data row for the current
    alternative's unit type, selected from the set of all unit specifications
    (called 'unit_data' in the main file)
  lag (int): the amount of delay before the start of construction for the 
    current project alternative
  fc_pd: total forecast period for all unit types. Defined by the longest 
    lag + construction duration + economic life of any available project
    alternative.

Returns:
  capex_column: the construction expenditure profile of the project, padded
    with zeros to indicate no construction costs outside the construction
    period.
"""
function generate_capex_profile(db, settings, current_pd, unit_type_data, lag, fc_pd, C2N_specs)
    # No construction expenditures before the project begins
    head_zeros = zeros(lag)

    if occursin("C2N", unit_type_data[:unit_type])
        capex_tl, activity_schedule = project_C2N_capex(db, settings, unit_type_data, lag, fc_pd, current_pd, C2N_specs)
        capex = capex_tl[!, :total_capex]
    else
        # Uniformly distribute projected total project costs over the construction
        #   period
        capex_per_pd = unit_type_data[:uc_x] * unit_type_data[:capacity] * MW2kW / unit_type_data[:d_x]
        capex = ones(convert(Int64, ceil(round(unit_type_data[:d_x], digits=3)))) .* capex_per_pd
    end

    # No construction expenditures from the end of construction until the end
    #   of the model's forecast period
    tail_zeros = zeros(fc_pd - lag - size(capex)[1])

    # Concatenate the above series into one
    capex_column = vcat(head_zeros, capex, tail_zeros)

    return capex_column
end


function project_C2N_capex(db, settings, unit_type_data, lag, fc_pd, current_pd, C2N_specs)
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
    capex_tl, activity_schedule = C2N.create_C2N_capex_timeline(db, conversion_type, rxtr_type, current_pd, lag, fc_pd, data, unit_type_data)

    return capex_tl, activity_schedule

end


"""
    set_initial_debt_principal_series(unit_fs, unit_type_data, lag, agent_params)

Set up the record of the accrual of debt during construction.
"""
function set_initial_debt_principal_series(unit_fs, unit_type_data, lag, agent_params)
    for i = lag+1:lag+convert(Int64, ceil(round(unit_type_data[:d_x], digits=3)))
        unit_fs[i, :remaining_debt_principal] = sum(unit_fs[1:i, :capex]) * agent_params[1, :debt_fraction]
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
    unit_d_x = convert(Int64, ceil(round(unit_type_data[:d_x], digits=3)))
    unit_op_life = unit_type_data[:unit_life]
    financing_term = 30
    rev_head_start = ceil(Int64, round(unit_type_data[:rev_head_start], digits=3))

    # Find the end of the capital expenditures period
    capex_end = 0
    for i = 1:size(unit_fs)[1]
        if round(unit_fs[i, :capex], digits=3) > 0
            capex_end = i
        end
    end

    for i = (capex_end+1-rev_head_start):(capex_end+1-rev_head_start+financing_term)
        # Apply a constant debt payment (sinking fund at cost of debt), based
        #   on the amount of debt outstanding at the end of the xtr project
        unit_fs[i, :debt_payment] = unit_fs[lag + unit_d_x, :remaining_debt_principal] .* cod ./ (1 - (1+cod) .^ (-1*financing_term))

        # Determine the portion of each payment which pays down interest
        #   (instead of principal)
        unit_fs[i, :interest_payment] = unit_fs[i-1, :remaining_debt_principal] * cod

        # Update the amount of principal remaining at the end of the period
        unit_fs[i, :remaining_debt_principal] = unit_fs[i-1, :remaining_debt_principal] - (unit_fs[i, :debt_payment] - unit_fs[i, :interest_payment])
    end

    dep_term = 20
    for i = (capex_end+1-rev_head_start):(capex_end+1-rev_head_start+dep_term)
        # Apply straight-line depreciation, based on debt outstanding at the
        #   project's completion: all units have a dep. life of 20 years
        unit_fs[i, :depreciation] = unit_fs[capex_end, :remaining_debt_principal] ./ dep_term
    end

end


"""
    forecast_unit_revenue_and_gen(unit_type_data, unit_fs, db, pd, lag; mode, ret_pd)

Forecast the revenue which will be earned by the current unit type,
and the total generation (in kWh) of the unit, per time period.

This function assumes the unit will be a pure price-taker: it will generate
  during all hours for which it is a marginal or submarginal unit.
If the unit is of a VRE type (as specified in the A-LEAF inputs), then a flat
  de-rating factor is applied to its availability during hours when it is
  eligible to generate.
"""
function forecast_unit_revenue_and_gen(unit_type_data, unit_fs, db, current_pd, lag, long_econ_results; mode="new_xtr", orig_ret_pd)
    check_valid_vector_mode(mode)

    # If the unit is VRE, assign an appropriate availability derate factor
    availability_derate_factor = compute_VRE_derate_factor(unit_type_data)

    # Compute the original retirement period
    # Minimum of ret_pd or size of the unit FS (i.e. the forecast period)
    if orig_ret_pd == nothing
        orig_ret_pd = size(unit_fs)[1]
    end

    # Compute the unit's total generation for each period, in kWh
    compute_total_generation(current_pd, unit_type_data, unit_fs, availability_derate_factor, lag, long_econ_results; mode=mode, orig_ret_pd=orig_ret_pd)

    # Compute total projected revenue, with VRE adjustment if appropriate, and
    #   save to the unit financial statement
    compute_total_revenue(current_pd, unit_type_data, unit_fs, availability_derate_factor, lag, long_econ_results; mode=mode, orig_ret_pd=orig_ret_pd)

end


"""
    compute_submarginal_hours_revenue(unit_type_data, price_curve)

Compute revenue accrued to the unit during hours when it is a submarginal
bidder into the wholesale market, assuming the unit is a pure price taker.
"""
function compute_submarginal_hours_revenue(unit_type_data, price_curve)
    # Get a list of all hours and their prices during which the unit would 
    #   be sub-marginal, from the price-duration curve
    # Marginal cost unit conversion:
    #   MC [$/MWh] = VOM [$/MWh] + FC_per_MWh [$/MWh]
    submarginal_hours = filter(row -> row.lamda > unit_type_data[:VOM] + unit_type_data[:FC_per_MWh] - unit_type_data[:policy_adj_per_MWh], price_curve)

    # Create a conversion factor from number of periods in the price curve
    #   to hours (price curve may be in hours or five-minute periods)
    convert_to_hours = hours_per_year / size(price_curve)[1]

    # Compute total number of submarginal hours
    num_submarg_hours = size(submarginal_hours)[1] * convert_to_hours

    # Calculate total revenue from submarginal hours, adjusting for net penalty
    #   or subsidy due to the effect of policies
    submarginal_hours_revenue = sum((submarginal_hours[!, :lamda] .+ unit_type_data[:policy_adj_per_MWh]) * unit_type_data[:capacity]) * convert_to_hours

    return num_submarg_hours, submarginal_hours_revenue
end


"""
    compute_marginal_hours_revenue(unit_type_data, price_curve, db)

Compute revenue accrued to the unit during hours when it is the marginal
bidder into the wholesale market, assuming the unit is a pure price taker.
The unit "takes credit" for a percentage of the marginal-hour revenue
corresponding to the percentage of unit-type capacity it comprises.
For example, a 200-MW NGCC unit in a system with 1000 MW of total NGCC capacity
receives 200/1000 = 20% of the revenues during hours when NGCC is marginal.
"""
function compute_marginal_hours_revenue(unit_type_data, price_curve, db, pd)
    # Get a list of all hours and their prices during which the unit would 
    #   be marginal, from the price-duration curve
    # Marginal cost unit convertion:
    #   MC [$/MWh] = VOM [$/MWh] + FC_per_MWh [$/MWh]
    marginal_hours = filter(row -> row.lamda == unit_type_data[:VOM] + unit_type_data[:FC_per_MWh] - unit_type_data[:policy_adj_per_MWh], price_curve)

    # Compute the total capacity of this unit type in the current system
    command = string("SELECT asset_id FROM assets WHERE unit_type == '", unit_type_data[:unit_type], "' AND completion_pd <= ", pd, " AND cancellation_pd > ", pd, " AND retirement_pd > ", pd)
    system_type_list = DBInterface.execute(db, command) |> DataFrame
    system_type_capacity = size(system_type_list)[1]

    # Create a conversion factor from number of periods in the price curve
    #   to hours (price curve may be in hours or five-minute periods)
    convert_to_hours = hours_per_year / size(price_curve)[1]

    # Compute effective number of marginal hours
    if system_type_capacity == 0
        num_marg_hours = 0
    else
        num_marg_hours = size(marginal_hours)[1] / system_type_capacity * convert_to_hours
    end

    if size(marginal_hours)[1] == 0
        marginal_hours_revenue = 0
    else
        marginal_hours_revenue = sum((marginal_hours[!, :lamda]) .+ unit_type_data[:policy_adj_per_MWh]) * unit_type_data[:capacity] / system_type_capacity * convert_to_hours
    end

    return num_marg_hours, marginal_hours_revenue
end


function compute_historical_unit_type_results(unit_specs, price_data, db, current_pd)
    # Create dataframe to store results
    unit_hist_results = DataFrame(
                            unit_type = String[],
                            total_gen = Float64[],
                            total_rev = Float64[]
                        )

    for unit_type in unit_specs[!, :unit_type]
        # Filter unit type data for convenience
        unit_type_data = filter(:unit_type => x -> x == unit_type, unit_specs)[1, :]

        # Get the unit type's number of hours generated and revenue while
        #   submarginal
        unit_submarg_hours, unit_submarg_hours_revenue = compute_submarginal_hours_revenue(unit_type_data, price_data)

        # Get the unit type's number of hours generated and revenue while
        #   marginal
        unit_marg_hours, unit_marg_hours_revenue = compute_marginal_hours_revenue(unit_type_data, price_data, db, current_pd)

        unit_total_rev = unit_submarg_hours_revenue + unit_marg_hours_revenue
        unit_total_gen = (unit_submarg_hours + unit_marg_hours) * unit_type_data[:capacity]

        push!(unit_hist_results, [unit_type unit_total_gen unit_total_rev])
    end

    return unit_hist_results

end


"""
    compute_VRE_derate_factor(unit_type_data)

If the unit is declared as `is_VRE` in the A-LEAF inputs, return an
availability derating factor equal to its capacity credit factor.
"""
function compute_VRE_derate_factor(unit_type_data)
    availability_derate_factor = 1

    # Julia reads booleans as UInt8, so have to convert to something sensible
    if convert(Int64, unit_type_data[:is_VRE][1]) == 1
        availability_derate_factor = unit_type_data[:CF]
    end

    return availability_derate_factor   
end



"""
    compute_total_revenue(unit_type_data, unit_fs, availability_derate_factor, lag; mode, orig_ret_pd)

Compute the final projected revenue stream for the current unit type, adjusting
unit availability if it is a VRE type.
"""
function compute_total_revenue(current_pd, unit_type_data, unit_fs, availability_derate_factor, lag, long_econ_results; mode, orig_ret_pd=9999)
    check_valid_vector_mode(mode)

    # Helpful short variables
    unit_d_x = convert(Int64, ceil(round(unit_type_data[:d_x], digits=3)))
    unit_op_life = unit_type_data[:unit_life]

    # Add a Revenue column to the financial statement dataframe
    unit_fs[!, :Revenue] .= 0.0

    # In "new_xtr" mode, revenues START accruing after the lag plus
    #   construction duration
    # In "retirement" mode, revenues start in the current period and CEASE accruing
    #   after the soonest of {the lag, the mandatory retirement period, the end
    #   of the forecast period}.
    if mode == "new_xtr"
        # Find the end of the capital expenditures period
        capex_end = 0
        for i = 1:size(unit_fs)[1]
            if round(unit_fs[i, :capex], digits=3) > 0
                capex_end = i
            end
        end
        rev_start = capex_end + 1 - ceil(Int64, round(unit_type_data[:rev_head_start], digits=3))
        rev_end = min(rev_start + unit_op_life, size(unit_fs)[1])
    elseif mode == "retirement"
        rev_start = 1
        # The maximum end of revenue period is capped at the length of the
        #   unit_fs dataframe (as some units with unspecified retirement
        #   periods default to a retirement period of 9999).
        rev_end = min(orig_ret_pd, size(unit_fs)[1]-1)
    end

    # If the project is a C2N project, use the AdvancedNuclear results
    # Otherwise, use the unit's own results
    type_filter = unit_type_data[:unit_type]
    if occursin("C2N", unit_type_data[:unit_type])
        type_filter = "AdvancedNuclear"
    end

    # Compute final projected revenue series
    agg_econ_results = combine(groupby(long_econ_results, [:y, :unit_type]), [:annualized_rev_perunit, :annualized_policy_adj_perunit] .=> sum, renamecols=false)
    for y = rev_start:rev_end
        row = filter([:y, :unit_type] => (t, unit_type) -> (t == y) && (unit_type == type_filter), agg_econ_results)

        if size(row)[1] != 0
            unit_fs[y, :Revenue] = row[1, :annualized_rev_perunit] + row[1, :annualized_policy_adj_perunit]
        else
            unit_fs[y, :Revenue] = 0
        end
    end

    avg_cpp_life_rem = 15

    if occursin("C2N", unit_type_data[:unit_type])
        for y = (floor(Int64, unit_type_data[:cpp_ret_lead])+lag+1):(floor(Int64, unit_type_data[:cpp_ret_lead])+lag+avg_cpp_life_rem)
            coal_data = filter([:y, :unit_type] => (t, unit_type) -> (t == y) && (unit_type == "Coal"), agg_econ_results)
            if (size(coal_data)[1] != 0)
                unit_fs[y, :Revenue] = unit_fs[y, :Revenue] - unit_type_data["num_cpp_rets"]*(coal_data[1, :annualized_rev_perunit] + coal_data[1, :annualized_policy_adj_perunit])
            end
        end
    end

end

"""
    compute_total_generation(unit_type_data, unit_fs, availability_derate_factor, lag; mode, orig_ret_pd)

Calculate the unit's total generation for the period, in kWh.
"""
function compute_total_generation(current_pd, unit_type_data, unit_fs, availability_derate_factor, lag, long_econ_results; mode, orig_ret_pd)
    check_valid_vector_mode(mode)

    # Helpful short variable names
    unit_d_x = convert(Int64, ceil(round(unit_type_data[:d_x], digits=3)))
    unit_op_life = unit_type_data[:unit_life]

    # In "new_xtr" mode, generation persists from end of construction to end of life
    # In "retirement" mode, generation occurs from the current period to the
    #   soonest of {the lag, the mandatory retirement period, the end of the
    #   forecast period}.
    if mode == "new_xtr"
        # Find the end of the capital expenditures period
        capex_end = 0
        for i = 1:size(unit_fs)[1]
            if round(unit_fs[i, :capex], digits=3) > 0
                capex_end = i
            end
        end
        gen_start = capex_end + 1 - ceil(Int64, round(unit_type_data[:rev_head_start], digits=3))

        gen_end = min(gen_start + unit_op_life, size(unit_fs)[1])
    elseif mode == "retirement"
        gen_start = 1
        # The maximum end of generation period is capped at the length of the
        #   unit_fs dataframe (as some units with unspecified retirement
        #   periods default to a retirement period of 9999).
        gen_end = min(orig_ret_pd, size(unit_fs)[1])
    end

    # If the project is a C2N project, use the AdvancedNuclear results
    # Otherwise, use the unit's own results
    type_filter = unit_type_data[:unit_type]
    if occursin("C2N", unit_type_data[:unit_type])
        type_filter = "AdvancedNuclear"
        unit_fs[!, :coal_gen] .= 0.0
    end

    transform!(long_econ_results, [:gen, :Probability, :num_units] => ((gen, prob, num_units) -> gen .* prob .* 365 ./ num_units) => :annualized_gen_perunit)
    agg_econ_results = combine(groupby(long_econ_results, [:y, :unit_type]), :annualized_gen_perunit => sum, renamecols=false)

    # Distribute generation values time series
    for y = gen_start:gen_end
        row = filter([:y, :unit_type] => (t, unit_type) -> (t == y) && (unit_type == type_filter), agg_econ_results)

        if size(row)[1] != 0
            unit_fs[y, :gen] = row[1, :annualized_gen_perunit]
        else
            unit_fs[y, :gen] = 0
        end
    end

    #coal_avg_rem_life = 15

    #if occursin("C2N", unit_type_data[:unit_type])
    #    for y = (floor(Int64, unit_type_data[:cpp_ret_lead])+lag+1):(floor(Int64, unit_type_data[:cpp_ret_lead])+lag+coal_avg_rem_life)
    #        coal_data = filter([:y, :unit_type] => (t, unit_type) -> (t == y) && (unit_type == "Coal"), agg_econ_results)
    #        if (size(coal_data)[1] != 0)
    #            unit_fs[y, :coal_gen] = unit_type_data["num_cpp_rets"] * (coal_data[1, :annualized_gen_perunit])
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
function forecast_unit_op_costs(unit_type_data, unit_fs, lag; mode="new_xtr", orig_ret_pd)
    check_valid_vector_mode(mode)

    # Helpful short variable names
    unit_d_x = convert(Int64, ceil(round(unit_type_data[:d_x], digits=3)))
    unit_op_life = unit_type_data[:unit_life]
    rev_head_start = ceil(Int64, round(unit_type_data[:rev_head_start], digits=3))
    cpp_ret_lead = floor(Int64, round(unit_type_data[:cpp_ret_lead], digits=3))

    if !occursin("C2N", unit_type_data[:unit_type])
        # Compute total fuel cost
        # Unit conversions:
        #  gen [MWh/year] * FC_per_MWh [$/MWh] = $/year
        transform!(unit_fs, [:gen] => ((gen) -> gen .* unit_type_data[:FC_per_MWh]) => :Fuel_Cost)

        # Compute total VOM cost incurred during generation
        # Unit conversions:
        #   gen [MWh/year] * VOM [$/MWh] = $/year
        transform!(unit_fs, [:gen] => ((gen) -> gen .* unit_type_data[:VOM]) => :VOM_Cost)
    else
        transform!(unit_fs, [:gen, :coal_gen] => ((gen, coal_gen) -> gen .* unit_type_data[:FC_per_MWh] .- coal_gen .* 143.07) => :Fuel_Cost)
        transform!(unit_fs, [:gen, :coal_gen] => ((gen, coal_gen) -> gen .* unit_type_data[:VOM] .- coal_gen .* 4.4) => :VOM_Cost)
    end

    # Compute total FOM cost for each period
    # Unit conversions:
    #   FOM [$/kW-year] * capacity [MW] * [1000 kW / 1 MW] = $/year
    if mode == "new_xtr"
        # Find the end of the capital expenditures period
        capex_end = maximum(filter(:capex => capex -> capex > 0, unit_fs)[!, :year])
        pre_zeros = zeros(capex_end-rev_head_start)
        op_ones = ones(unit_op_life)
    elseif mode == "retirement"
        pre_zeros = zeros(0)
        op_ones = ones(min(size(unit_fs)[1], orig_ret_pd))
    end
    post_zeros = zeros(size(unit_fs)[1] - size(pre_zeros)[1] - size(op_ones)[1])
    unit_fs[!, :FOM_Cost] = vcat(pre_zeros, op_ones, post_zeros)

    unit_fs[!, :FOM_Cost] .= unit_fs[!, :FOM_Cost] .* unit_type_data[:FOM] .* unit_type_data[:capacity] .* MW2kW
    if occursin("C2N", unit_type_data[:unit_type])
        for y = (lag+1+cpp_ret_lead):(lag+1+cpp_ret_lead+unit_op_life)
            unit_fs[y, :FOM_Cost] = unit_fs[y, :FOM_Cost] - (unit_type_data["num_cpp_rets"] * 632 * MW2kW * 39.70)
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
    transform!(unit_fs, [:Revenue, :Fuel_Cost, :VOM_Cost, :FOM_Cost] => ((rev, fc, VOM, FOM) -> rev - fc - VOM - FOM) => :EBITDA)

    # Compute EBIT
    transform!(unit_fs, [:EBITDA, :depreciation] => ((EBITDA, dep) -> EBITDA - dep) => :EBIT)

    # Compute EBT
    transform!(unit_fs, [:EBIT, :interest_payment] => ((EBIT, interest) -> EBIT - interest) => :EBT)

    # Retrieve the system corporate tax rate from the database
    command = string("SELECT value FROM model_params WHERE parameter == 'tax_rate'")
    # Extract the value into a temporary dataframe
    tax_rate = DBInterface.execute(db, command) |> DataFrame
    # Pull out the bare value
    tax_rate = tax_rate[1, :value]

    # Compute taxes owed
    transform!(unit_fs, [:EBT] => ((EBT) -> EBT .* tax_rate) => :Tax_Owed)

    # Compute net income
    transform!(unit_fs, [:EBT, :Tax_Owed] => ((EBT, tax_owed) -> EBT - tax_owed) => :Net_Income)

    # Compute free cash flow (FCF)
    transform!(unit_fs, [:Net_Income, :interest_payment, :capex] => ((NI, interest, capex) -> NI + interest - capex) => :FCF)

end


function compute_FS_delta_value(unit_fs, lag, db)
    # Make a copy of the baseline fs to use for the early-retirement case
    early_ret_fs = deepcopy(unit_fs)

    # Zero out all values which should be zero due to early retirement
    columns_to_adjust = [:gen, :Revenue, :VOM_Cost, :Fuel_Cost, :FOM_Cost]
    for col in columns_to_adjust
        # Set values to 0, starting at the end of the lag period
        early_ret_fs[lag+1:size(unit_fs)[1], col] .= 0
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
    #d = agent_params[1, :debt_fraction] * agent_params[1, :cost_of_debt] + (1 - agent_params[1, :debt_fraction]) * agent_params[1, :cost_of_equity]
    d = agent_params[1, :cost_of_debt]

    # Add a column of compounded discount factors to the dataframe
    transform!(unit_fs, [:year] => ((year) -> (1+d) .^ (-1 .* (year .- 1))) => :discount_factor)

    # Discount the alternative's FCF NPV
    FCF_NPV = transpose(unit_fs[!, :FCF]) * unit_fs[!, :discount_factor]

    return FCF_NPV, unit_fs

end


function set_baseline_future_years(db, agent_id)
    all_unit_manifest = DBInterface.execute(db, "SELECT *, COUNT(asset_id) FROM assets GROUP BY unit_type, completion_pd, retirement_pd")
    owned_unit_manifest = filter(:agent_id => x -> x == agent_id, all_unit_manifest)

    num_forecast_years = 10

    all_extant_units = Dict()
    owned_extant_units = Dict()

    for i = 1:num_forecast_years
        # Determine the total number of extant units by type during this year
        aeudf = groupby(filter([:completion_pd, :retirement_pd] => (c, r) -> (c <= i) && (r > i), unit_manifest), :unit_type)
        all_extant_units[i] = combine(aeudf, Symbol("COUNT(asset_id)") => sum)

        # Determine the number of owned extant (operational) units by type during this year
        oeudf = groupby(filter([:completion_pd, :retirement_pd] => (c, r) -> (c <= i) && (r > i), unit_manifest), :unit_type)
        owned_extant_units[i] = combine(oeudf, Symbol("COUNT(asset_id)") => sum)
    end

    


end


### JuMP optimization model initialization
"""
    set_up_model(unit_FS_dict, ret_fs_dict, fc_pd, available_demand, NPV_results, ret_NPV_results)

Set up the JuMP optimization model, including variables, constraints, and the
objective function.

Returns:
  m (JuMP model object)
"""
function set_up_model(settings, PA_uids, PA_fs_dict, total_demand, asset_counts, agent_params, unit_specs, current_pd, system_portfolios, db, agent_id, agent_fs, fc_pd)
    # Create the model object
    # @info "Setting up model..."

    solver = lowercase(settings["simulation"]["solver"])
    if solver == "cplex"
        # using CPLEX
        m = Model(CPLEX.Optimizer)
    elseif solver == "glpk"
        m = Model(GLPK.Optimizer)
    elseif solver == "cbc"
        m = Model(Cbc.Optimizer)
    elseif solver == "scip"
        m = Model(SCIP.Optimizer)
    elseif solver == "highs"
        m = Model(HiGHS.Optimizer)
    else
        throw(error("Solver `$solver` not supported. Try `cplex` instead."))
    end
    set_silent(m)

    # Parameter names
    num_alternatives = size(PA_uids)[1]
    num_time_periods = size(PA_fs_dict[PA_uids[1, :uid]])[1]

    # Set up variables
    # Number of units of each type to build: must be Integer
    @variable(m, u[1:num_alternatives] >= 0, Int)

    transform!(total_demand, [:real_demand, :total_eff_cap] => ((dem, eff_cap) -> eff_cap - dem) => :excess_capacity)

    for i = 1:size(total_demand)[1]
        if total_demand[i, :excess_capacity] < 0
            total_demand[i, :excess_capacity] = 0
        end
    end

    # Compute expected marginal generation and effective nameplate capacity
    #   contribution per alternative type
    marg_gen = zeros(num_alternatives, num_time_periods)
    marg_eff_cap = zeros(num_alternatives, num_time_periods)
    for i = 1:size(PA_uids)[1]
        for j = 1:num_time_periods
            marg_gen[i, j] = PA_fs_dict[PA_uids[i, :uid]][j, :gen] / MW2kW    # in kWh
            # Convert total anticipated marginal generation to an effective
            #   nameplate capacity and save to the appropriate entry in the
            #   marginal generation array
            unit_type_data = filter(:unit_type => x -> x == PA_uids[i, :unit_type], unit_specs)
            if PA_uids[i, :project_type] == "new_xtr"
                if j > PA_uids[i, :lag] + unit_type_data[1, :d_x]
                    marg_eff_cap[i, j] = unit_type_data[1, :capacity] * unit_type_data[1, :CF]
                end
            elseif PA_uids[i, :project_type] == "retirement"
                if (j > PA_uids[i, :lag]) && (j <= PA_uids[i, :ret_pd])
                    marg_eff_cap[i, j] = (-1) * unit_type_data[1, :capacity] * unit_type_data[1, :CF]
                end
            end
        end
    end

    # Prevent excessive overbuild in the medium term
#    for i = 3:10
#        @constraint(m, transpose(u) * marg_eff_cap[:, i] + sum(system_portfolios[i][!, :effective_capacity]) <= forecasted_demand[current_pd+i, :total_demand] * 1.1)   # settings["planning_reserve_margin"]
#    end

    # Prevent the agent from intentionally causing foreseeable energy shortages
    shortage_protection_pd = 8
    @info "pd:  peak_demand  eff_capacity"
    for i = 1:shortage_protection_pd
        k = current_pd + i
        pd_total_demand = filter(:period => x -> x == current_pd + i - 1, total_demand)[1, :total_demand]
        pdtd = round(pd_total_demand, digits=1)
        total_eff_cap = filter(:period => x -> x == current_pd + i - 1, total_demand)[1, :total_eff_cap]
        tec = round(total_eff_cap, digits=1)
        @info "$k:     $pdtd      $tec"
        if (total_eff_cap > pd_total_demand)
            margin = -0.3
        else
            margin = 0.0
        end
        @constraint(m, transpose(u) * marg_eff_cap[:, i] >= (total_eff_cap - pd_total_demand) * margin)

    end

    filename = joinpath(
                   pwd(),
                   "tmp",
                   "marg_eff_cap.csv"
               )
    CSV.write(filename, DataFrame(marg_eff_cap, :auto))

    # Create arrays of expected marginal debt, interest, dividends, and FCF per unit type
    marg_debt = zeros(num_alternatives, num_time_periods)
    marg_int = zeros(num_alternatives, num_time_periods)
    marg_div = zeros(num_alternatives, num_time_periods)
    marg_FCF = zeros(num_alternatives, num_time_periods)
    for i = 1:size(PA_uids)[1]
        for j = 1:num_time_periods
            # Retrieve the marginal value of interest
            # Scale to units of $B
            marg_debt[i, j] = PA_fs_dict[PA_uids[i, :uid]][j, :remaining_debt_principal] / 1e9
            marg_int[i, j] = PA_fs_dict[PA_uids[i, :uid]][j, :interest_payment] / 1e9
            marg_div[i, j] = PA_fs_dict[PA_uids[i, :uid]][j, :remaining_debt_principal] * agent_params[1, :cost_of_equity] / 1e9
            marg_FCF[i, j] = PA_fs_dict[PA_uids[i, :uid]][j, :FCF] / 1e9
        end
    end

    filename = joinpath(
                   pwd(),
                   "tmp",
                   "marg_FCF.csv"
               )
    CSV.write(filename, DataFrame(marg_int, :auto))

    agent_fs_path = joinpath(
                        pwd(),
                        "tmp",
                        string(
                            "agent_",
                            agent_id,
                            "_pd_",
                            current_pd,
                            ".csv"
                        )
                    )
    CSV.write(agent_fs_path, agent_fs)

    # Prevent the agent from reducing its credit metrics below Moody's Baa
    #   rating thresholds (from the Unregulated Power Companies ratings grid)
    for i = 1:6
        @constraint(m, agent_fs[i, :FCF] / 1e9 + sum(u .* marg_FCF[:, i]) + (1 - 4.2) * (agent_fs[i, :interest_payment] / 1e9 + sum(u .* marg_int[:, i])) >= 0)
        # These constraints need to be rewritten
        #@constraint(m, (agent_NI / 1e9 + sum((transpose(u) .* marg_FCF[:, i]))) / (0.2) - (agent_principal / 1e9 + sum(transpose(u) .* marg_debt[:, i])) >= 0)
        #@constraint(m, (agent_NI + (transpose(u) * marg_FCF[:, i])) + sum((transpose(u) .* (marg_FCF[:, i] - marg_div[:, i] - 0.15 .* marg_debt[:, i]))) >= 0)
    end

    # Enforce the user-specified maximum number of new construction/retirement
    #   projects by type per period, and the :allowed field in PA_uids
    for i = 1:num_alternatives
        if PA_uids[i, :project_type] == "new_xtr"
            @constraint(m, u[i] .<= convert(Int64, PA_uids[i, :allowed]) .* settings["agent_opt"]["max_type_newbuilds_per_pd"])
        elseif PA_uids[i, :project_type] == "retirement"
            unit_type = PA_uids[i, :unit_type]
            ret_pd = PA_uids[i, :ret_pd]
            asset_count = filter([:unit_type, :retirement_pd] => (x, y) -> x == unit_type && y == ret_pd, asset_counts)[1, :count]
            max_retirement = convert(Int64, PA_uids[i, :allowed]) .* min(asset_count, settings["agent_opt"]["max_type_rets_per_pd"])
#            @constraint(m, u[i] .<= max_retirement)
        end
    end

    # The total number of assets of each type to be retired cannot exceed
    #   the total number of assets of that type owned by the agent

    # Setup
    # Convenient variable for the number of distinct retirement alternatives,
    #   excluding assets reserved for C2N conversion
    retireable_asset_counts = filter([:C2N_reserved] => ((reserved) -> reserved == 0), asset_counts)
    k = size(retireable_asset_counts)[1]
    # Shortened name for the number of lag periods to consider
    #   1 is added, as the user-set value only specifies future periods,
    #   with the "lag = 0" instance being implied
    num_lags = settings["agent_opt"]["num_future_periods_considered"] + 1
    # Create the matrix to collapse lagged options
    # This matrix has one long row for each retiring-asset category,
    #   with 1's in each element where the corresponding element of u[] is
    #   one of these units
    ret_summation_matrix = zeros(k, size(PA_uids)[1])
    for i = 1:k
        for j = 1:size(PA_uids)[1]
            if (PA_uids[j, :project_type] == "retirement") & (PA_uids[j, :unit_type] == retireable_asset_counts[i, :unit_type]) & (PA_uids[j, :ret_pd] == retireable_asset_counts[i, :retirement_pd])
                ret_summation_matrix[i, j] = 1
            end
        end
    end

    filename = joinpath(
                   pwd(),
                   "tmp",
                   "ret_rum_mat_$agent_id.csv"
               )
    CSV.write(filename, DataFrame(ret_summation_matrix, :auto))

    # Specify constraint: the agent cannot plan to retire more units (during
    #   all lag periods) than exist of that unit type
    for i = 1:size(retireable_asset_counts)[1]
        @constraint(m, sum(ret_summation_matrix[i, :] .* u) <= retireable_asset_counts[i, :count])
    end

    # Create the objective function 

    lamda_1 = 1.0 / 1e9
    lamda_2 = 0.0
    lim = 6
    int_bound = 5.0
 
    @objective(
        m,
        Max,
        (
            lamda_1 * (transpose(u) * PA_uids[!, :NPV])
            + lamda_2 * (
                sum(agent_fs[1:lim, :FCF]) / 1e9
                + sum(transpose(u) * marg_FCF[:, 1:lim])
                + sum(agent_fs[1:lim, :interest_payment]) / 1e9
                + sum(transpose(u) * marg_int[:, 1:lim])
                - (int_bound) * (sum(agent_fs[1:lim, :interest_payment]) / 1e9 + sum(transpose(u) * marg_int[:, 1:lim]))
            )
        )
    )


    # @info "Optimization model set up."

    return m
end


### Postprocessing
function postprocess_agent_decisions(all_results, unit_specs, db, current_pd, agent_id)
    for i = 1:size(all_results)[1]
        # Retrieve the individual result row for convenience
        result = all_results[i, :]

        if result[:project_type] == "new_xtr"
        # New construction decisions are binding for this period only
            if (result[:lag] == 0) && (result[:units_to_execute] != 0)
                record_new_construction_projects(result, unit_specs, db, current_pd, agent_id)
                # If the project is a C2N project, retire two coal units at
                #    the appropriate time
                if occursin("C2N", result[:unit_type])
                    record_asset_retirements(result, db, agent_id, unit_specs, current_pd, mode="C2N_newbuild")
                end
            end
        elseif result[:project_type] == "retirement"
        # Retirements are recorded as binding, whether the retirement date is
        #   this period or in the future
            if result[:units_to_execute] != 0
                record_asset_retirements(result, db, agent_id, unit_specs, current_pd)
            end
        else
            @warn "I'm not sure what to do with the following project type:"
            @warn result
            @warn "This decision entry will be skipped."
        end

    end

end


function record_new_construction_projects(result, unit_data, db, current_pd, agent_id)
    # Retrieve unit_specs data for this unit type
    unit_type_specs = filter(:unit_type => x -> x == result[:unit_type], unit_data)

    # Set default initial values
    cum_occ = unit_type_specs[1, :uc_x] * unit_type_specs[1, :capacity] * MW2kW
    rcec = cum_occ
    cum_exp = 0
    cum_d_x = unit_type_specs[1, :d_x]
    rtec = cum_d_x
    start_pd = current_pd
    completion_pd = current_pd + unit_type_specs[1, :d_x]
    cancellation_pd = 9999
    retirement_pd = current_pd + unit_type_specs[1, :d_x] + unit_type_specs[1, :unit_life]
    total_capex = 0
    cap_pmt = 0
    anpe = unit_type_specs[1, :uc_x] * unit_type_specs[1, :capacity] * 1000 / unit_type_specs[1, :d_x]
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
            cum_d_x,
            rtec,
            cum_exp,
            anpe
        )
        DBInterface.execute(
            db,
            "INSERT INTO WIP_updates VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            WIP_projects_vals
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
            C2N_reserved
        )
        DBInterface.execute(
            db,
            "INSERT INTO asset_updates VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            assets_vals
        )
    end

end


function record_asset_retirements(result, db, agent_id, unit_specs, current_pd; mode="standalone")
    if mode == "standalone"
        # Generate a list of assets which match 'result' on agent owner, 
        #   unit type, and mandatory retirement date
        match_vals = (result[:unit_type], result[:ret_pd], agent_id)
        ret_candidates = DBInterface.execute(
            db,
            "SELECT asset_id FROM assets WHERE unit_type = ? AND retirement_pd = ? AND agent_id = ?",
            match_vals
        ) |> DataFrame

        # Set the number of units to execute
        units_to_execute = result[:units_to_execute]

        # Set the new retirement pd
        new_ret_pd = result[:lag]

    elseif mode == "C2N_newbuild"
        # Determine the period in which the coal units must retire
        unit_type_specs = filter(:unit_type => x -> x == result[:unit_type], unit_specs)[1, :]
        new_ret_pd = current_pd + result[:lag] + convert(Int64, ceil(round(unit_type_specs[:cpp_ret_lead], digits=3)))
        C2N_reserved = 0 # Ensure a coal unit cannot be selected twice for C2N decommissioning

        # Generate a list of coal units which will still be operational by
        #   the necessary period
        match_vals = (new_ret_pd, new_ret_pd, agent_id, C2N_reserved)
        ret_candidates = DBInterface.execute(
            db,
            "SELECT asset_id FROM assets WHERE unit_type = 'Coal' AND completion_pd <= ? AND retirement_pd > ? AND agent_id = ? AND C2N_reserved = ?",
            match_vals
        ) |> DataFrame

        # Set the number of units to execute
        units_to_execute = unit_type_specs["num_cpp_rets"]
        # @info "Units to execute: ", units_to_execute

    end

    # Retire as many of these matching assets as is indicated by the agent
    #   optimization result
    for j = 1:units_to_execute
        asset_to_retire = ret_candidates[convert(Int64, j), :asset_id]
        asset_data = DBInterface.execute(
            db,
            "SELECT * FROM assets WHERE asset_id = $asset_to_retire"
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
            replacement_data
        )
    end

end


function update_agent_financial_statement(agent_id, db, unit_specs, current_pd, fc_pd, long_econ_results, all_year_portfolios)
    # Retrieve horizontally-abbreviated dataframes
    short_econ_results = select(long_econ_results, [:unit_type, :y, :d, :h, :gen, :annualized_rev_perunit, :annualized_VOM_perunit, :annualized_FC_perunit, :annualized_policy_adj_perunit])
    short_unit_specs = select(unit_specs, [:unit_type, :FOM])

    # Inner join the year's portfolio with financial pivot
    fin_results = innerjoin(short_econ_results, all_year_portfolios, on = [:y, :unit_type])

    # Fill in total revenue
    transform!(fin_results, [:annualized_rev_perunit, :annualized_policy_adj_perunit, :num_units] => ((rev, adj, num_units) -> (rev .+ adj) .* num_units) => :total_rev)

    # Fill in total VOM
    transform!(fin_results, [:annualized_VOM_perunit, :num_units] => ((VOM, num_units) -> VOM .* num_units) => :total_VOM)

    # Fill in total fuel costs
    transform!(fin_results, [:annualized_FC_perunit, :num_units] => ((FC, num_units) -> FC .* num_units) => :total_FC)

    # Create the annualized results dataframe so far
    results_pivot = combine(groupby(fin_results, :y), [:total_rev, :total_VOM, :total_FC] .=> sum; renamecols=false)

    fs = DataFrame(
             base_pd = ones(Int64, fc_pd) .* current_pd,
             projected_pd = current_pd:current_pd+fc_pd-1,
             revenue = results_pivot[!, :total_rev],
             VOM = results_pivot[!, :total_VOM],
             fuel_costs = results_pivot[!, :total_FC]
         )

    # Fill in total FOM
    # Inner join the FOM table with the agent portfolios
    FOM_df = innerjoin(all_year_portfolios, short_unit_specs, on = :unit_type)
    transform!(FOM_df, [:FOM, :num_units] => ((FOM, num_units) -> FOM .* num_units .* 1000) => :total_FOM)
    FOM_df = select(combine(groupby(FOM_df, :y), :total_FOM => sum; renamecols=false), [:y, :total_FOM])
    rename!(FOM_df, :y => :projected_pd, :total_FOM => :FOM)

    fs = innerjoin(fs, FOM_df, on = :projected_pd)

    # Fill in total depreciation
    total_pd_depreciation = DBInterface.execute(db, "SELECT projected_pd, SUM(depreciation) FROM depreciation_projections WHERE agent_id = $agent_id AND base_pd = $current_pd GROUP BY projected_pd") |> DataFrame

    # Fill in total interest payments
    total_pd_interest = DBInterface.execute(db, "SELECT projected_pd, SUM(interest_payment) FROM financing_schedule WHERE agent_id = $agent_id AND base_pd = $current_pd GROUP BY projected_pd") |> DataFrame

    # Fill in total capex
    total_pd_capex = DBInterface.execute(db, "SELECT projected_pd, SUM(capex) FROM capex_projections WHERE agent_id = $agent_id AND base_pd = $current_pd GROUP BY projected_pd") |> DataFrame

    # Join the scheduled columns to the financial statement, correct column
    #    names, and replace missing values with zeroes
    fs = leftjoin(fs, total_pd_depreciation, on = :projected_pd)
    fs = leftjoin(fs, total_pd_interest, on = :projected_pd)
    fs = leftjoin(fs, total_pd_capex, on = :projected_pd)
    rename!(fs, Symbol("SUM(depreciation)") => :depreciation, Symbol("SUM(interest_payment)") => :interest_payment, Symbol("SUM(capex)") => :capex)
    fs[!, :depreciation] = coalesce.(fs[!, :depreciation], 0)
    fs[!, :interest_payment] = coalesce.(fs[!, :interest_payment], 0)
    fs[!, :capex] = coalesce.(fs[!, :capex], 0)


    ### Computed FS quantities
    # EBITDA
    transform!(fs, [:revenue, :VOM, :FOM, :fuel_costs] => ((rev, VOM, FOM, FC) -> (rev - VOM - FOM - FC)) => :EBITDA)

    # EBIT
    transform!(fs, [:EBITDA, :depreciation] => ((EBITDA, dep) -> (EBITDA - dep)) => :EBIT)

    # EBT
    transform!(fs, [:EBIT, :interest_payment] => ((EBIT, interest) -> EBIT - interest) => :EBT)

    # Tax owed
    tax_rate = 0.21
    transform!(fs, :EBT => ((EBT) -> EBT * tax_rate) => :tax_paid)

    # Net Income
    transform!(fs, [:EBT, :tax_paid] => ((EBT, tax) -> (EBT - tax)) => :Net_Income)

    # Free Cash Flow
    transform!(fs, [:Net_Income, :depreciation, :capex] => ((NI, dep, capex) -> NI + dep - capex) => :FCF)

    # Save the dataframe to the database
    save_agent_fs!(fs, agent_id, db)

    return fs


end


function save_agent_fs!(fs, agent_id, db)
    # Add the agent id to the dataframe
    fs[!, :agent_id] .= agent_id

    # Get the column order from the agent_financial_statements database table
    sql_cmd = "SELECT name FROM PRAGMA_TABLE_INFO('agent_financial_statements')"
    fs_col_order = DBInterface.execute(db, sql_cmd) |> DataFrame
    fs_col_order = collect(fs_col_order[!, "name"]) # convert DF col to vector

    # Reorder the agent fs to match the DB table
    fs = select(fs, fs_col_order)

    try
        tup_fs = Tuple.(eachrow(fs))
        for i = 1:size(fs)[1]
            values = tup_fs[i]
            ins_cmd = "INSERT INTO agent_financial_statements VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
            DBInterface.execute(db, ins_cmd, values)
        end
    catch e
        throw(e)
    end

end

function save_agent_decisions(db, agent_id, decision_df)
    decision_df[!, :agent_id] .= agent_id
    cols_to_ignore = [:uid]
    select!(decision_df, :agent_id, Not(vcat([:agent_id], cols_to_ignore)))
    tup_rows = Tuple.(eachrow(decision_df))
    for i=1:size(decision_df)[1]
        values = tup_rows[i]
        ins_cmd = "INSERT INTO agent_decisions VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
        DBInterface.execute(db, ins_cmd, values)
    end
end



end

