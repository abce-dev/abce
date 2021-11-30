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

using SQLite, DataFrames, CSV

export load_db, get_current_period, get_agent_id, get_agent_params, load_unit_type_data, set_forecast_period, extrapolate_demand, project_demand_flat, project_demand_exponential, allocate_fuel_costs, create_unit_FS_dict, get_unit_specs, get_table, show_table, get_WIP_projects_list, get_demand_forecast, get_net_demand, get_next_asset_id, ensure_projects_not_empty, authorize_anpe, add_xtr_events, create_NPV_results_df, generate_xtr_exp_profile, set_initial_debt_principal_series

#####
# Constants
#####
MW2kW = 1000   # Conversion factor from MW to kW


#####
# Setup functions
#####

function load_db(db_file)
    try
        db = SQLite.DB(db_file)
        return db
    catch e
        println("Couldn't load the database:")
        println(e)
        exit()
    end
end


function get_current_period()
    try
        pd = parse(Int64, ARGS[2])
        return pd
    catch e
        println("Couldn't retrieve a period number from the command line:")
        println(e)
        exit()
    end
end


function get_agent_id()
    try
        agent_id = parse(Int64, ARGS[3])
        return agent_id
    catch e
        println("Couldn't retrieve the agent ID from the command line:")
        println(e)
        exit()
    end
end


function get_agent_params(db, agent_id)
    try
        command = string("SELECT * FROM agent_params WHERE agent_id = ", agent_id)
        df = DBInterface.execute(db, command) |> DataFrame
    catch e
        println("Could not get agent parameters from file:")
        println(e)
        exit()
    end
end


function load_unit_type_data(unit_data_file)
    unit_data = CSV.read(unit_data_file, DataFrame)
    num_types = size(unit_data)[1]
    return unit_data, num_types
end


function set_forecast_period(df, num_lags)
    transform!(df, [:d_x, :unit_life] => ((lead_time, unit_life) -> lead_time + unit_life) => :full_life)
    max_horizon = maximum(df[!, :full_life]) + num_lags
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
    println(string("\nTable \'", table_name, "\':"))
    println(df)
    return df
end


function get_WIP_projects_list(db, pd, agent_id)
    # Get a list of all WIP (non-complete, non-cancelled) projects for the given agent
    SQL_get_proj = SQLite.Stmt(db, string("SELECT asset_id FROM assets WHERE agent_id = ", agent_id, " AND completion_pd > ", pd, " AND cancellation_pd > ", pd))
    project_list = DBInterface.execute(SQL_get_proj) |> DataFrame
    return project_list
end


function extrapolate_demand(visible_demand, db, pd, fc_pd, settings)
    mode = settings["demand_projection_mode"]
    if mode == "flat"
        demand = project_demand_flat(visible_demand, fc_pd)
    elseif mode == "exp_termrate"
        term_demand_gr = settings["terminal_demand_growth_rate"]
        demand = project_demand_exp_termrate(visible_demand, fc_pd, term_demand_gr)
    elseif mode == "exp_fitted"
        demand = project_demand_exp_fitted(visible_demand, db, pd, fc_pd, settings)
    else
        println(string("The specified demand extrapolation mode, ", mode, ", is not implemented."))
        println("Please use 'flat', 'exp_termrate', or 'exp_fitted' at this time.")
        println("Terminating...")
        exit()
    end

    return demand
end


function project_demand_flat(visible_demand, fc_pd)
    demand = DataFrame(demand = zeros(Float64, convert(Int64, fc_pd)))
    demand[1:size(visible_demand)[1], :demand] .= visible_demand[!, :demand]
    demand[(size(visible_demand)[1] + 1):fc_pd, :demand] .= demand[size(visible_demand)[1], :demand] 
    return demand
end


function project_demand_exp_termrate(visible_demand, fc_pd, term_demand_gr)
    future_demand = ones(fc_pd - size(visible_demand)[1])
    for i = 1:size(future_demand)[1]
        future_demand[i] = last(visible_demand)[1] * (1 + term_demand_gr) ^ i
    end
    demand = DataFrame(demand = vcat(visible_demand[!, :demand], future_demand))
    return demand
end


function project_demand_exp_fitted(visible_demand, db, pd, fc_pd, settings)
    # Retrieve historical data, if available
    demand_projection_window = settings["demand_projection_window"]
    demand_history_start = max(0, pd - settings["demand_visibility_horizon"])
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
    beta[2] = (beta[2] * size(y)[1] + settings["historical_demand_growth_rate"] * (demand_projection_window - size(y)[1])) / demand_projection_window

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
    demand_vis_horizon = settings["demand_visibility_horizon"]
    vals = (pd, pd + demand_vis_horizon)
    visible_demand = DBInterface.execute(db, "SELECT demand FROM demand WHERE period >= ? AND period < ?", vals) |> DataFrame
    demand_forecast = extrapolate_demand(visible_demand, db, pd, fc_pd, settings)
    prm = DBInterface.execute(db, "SELECT * FROM model_params WHERE parameter = 'PRM'") |> DataFrame
    demand_forecast = demand_forecast .* (1 + prm[1, :value]) .+ settings["peak_initial_reserves"]
    return demand_forecast
end


function get_net_demand(db, pd, agent_id, fc_pd, demand_forecast)
    # Calculate the amount of forecasted net demand in future periods
    installed_cap_forecast = DataFrame(period = Int64[], derated_capacity = Float64[])
    vals = (pd, pd)
    # Select a list of all current assets, which are not cancelled, retired, or hidden from public view
    current_assets = DBInterface.execute(db, "SELECT * FROM assets WHERE cancellation_pd > ? AND retirement_pd > ? AND revealed = 'true'", vals) |> DataFrame
    if size(current_assets)[1] == 0
        println("There are no currently-active generation assets in the system; unpredictable behavior may occur.")
    end
    current_assets[!, :capacity] = zeros(size(current_assets)[1])
    current_assets[!, :CF] = zeros(size(current_assets)[1])
    unit_specs = DBInterface.execute(db, "SELECT * FROM unit_specs") |> DataFrame
    for i = 1:size(current_assets)[1]
        asset_type = current_assets[i, :unit_type]
        unit = filter(row -> row[:unit_type] == asset_type, unit_specs)
        current_assets[i, :capacity] = unit[1, :capacity]
        current_assets[i, :CF] = unit[1, :CF]
        transform!(current_assets, [:capacity, :CF] => ((cap, cf) -> cap .* cf) => :derated_capacity)
    end
    for i=pd:pd+fc_pd-1
        future_active_assets = filter(row -> (row[:completion_pd] <= i) && (row[:retirement_pd] > i), current_assets)
        total_cap = sum(future_active_assets[!, :derated_capacity])
        df = DataFrame(period = i, derated_capacity = total_cap)
        append!(installed_cap_forecast, df)
    end
    net_demand_forecast = demand_forecast[!, :demand] - installed_cap_forecast[!, :derated_capacity]
    return net_demand_forecast
end


function get_next_asset_id(db)
    # Return the next available asset ID (one greater than the current largest ID)
    SQL_get_ids = SQLite.Stmt(db, string("SELECT asset_id FROM assets"))
    asset_df = DBInterface.execute(SQL_get_ids) |> DataFrame
    #asset_df[!, :asset_id] = asset_df[:, :asset_id]
    next_id = maximum(asset_df[!, :asset_id]) + 1
    return next_id    
end


function get_unit_specs(db)
    # Retrieve the table of unit specifications from the DB
    df = DBInterface.execute(db, "SELECT * FROM unit_specs") |> DataFrame
    num_types = size(df)[1]
    # Convert the Int-type columns to Int64
    df[!, :d_x] = convert.(Int64, df[:, :d_x])
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
            println(string("Created project ", new_asset_id))

            # Update the list of WIP construction projects and return it
            project_list = get_WIP_projects_list(db, agent_id)
            return project_list
        else
            # If at least one construction project already exists, there's no
            #    need to do anything.
            return project_list
        end
    catch e
        println("Could not insert a seed project into the agent's project list")
        println(e)
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
        DBInterface.execute(db, "UPDATE WIP_projects SET anpe = ? WHERE period = ? AND asset_id = ?", vals)
    end
end


#####
# NPV functions
#####

"""
    create_NPV_results_DF(unit_data_df, num_lags)

Create a dataframe to hold the results of NPV calculations for the various
  types, expanded by the number of allowed lags.
"""
function create_NPV_results_df(unit_data, num_lags)
    alternative_names = Vector{String}()
    num_alternatives = size(unit_data)[1] * (num_lags + 1)

    for i = 1:size(unit_data)[1]
        for j = 0:num_lags
            name = string(unit_data[i, :unit_type], "_lag-", j)
            push!(alternative_names, name)
        end
    end

    NPV_results = DataFrame(name=alternative_names, NPV=zeros(num_alternatives))
    return alternative_names, NPV_results

end


function create_unit_FS_dict(unit_data, fc_pd, num_lags)
    fs_dict = Dict()
    num_types = size(unit_data)[1]
    for i = 1:num_types
        for j = 0:num_lags
            short_name = unit_data[i, :unit_type]
            unit_name = string(short_name, "_lag-", j)
            unit_FS = DataFrame(year = 1:fc_pd, xtr_exp = zeros(fc_pd), gen = zeros(fc_pd), remaining_debt_principal = zeros(fc_pd), debt_payment = zeros(fc_pd), interest_due = zeros(fc_pd), depreciation = zeros(fc_pd))
            fs_dict[unit_name] = unit_FS
        end
    end
    return fs_dict
end


"""
    populate_unit_alternative_FS(unit_type, unit_data, lag, fc_pd)

For a given project alternative (unit type + lag duration), calculate out its
marginal contribution to the financial statements.
"""
function populate_unit_alternative_FS(unit_type, unit_data, lag, unit_FS_dict, fc_pd)
    name = string(unit_type, "_lag-", lag)
    fs = unit_FS_dict[name]

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
  xtr_exp_column: the construction expenditure profile of the project, padded
    with zeros to indicate no construction costs outside the construction
    period.
"""
function generate_xtr_exp_profile(unit_type_data, lag, fc_pd)
    # No construction expenditures before the project begins
    head_zeros = zeros(lag)

    # Uniformly distribute projected total project costs over the construction
    #   period
    xtr_exp_per_pd = unit_type_data[1, :uc_x] * unit_type_data[1, :capacity] * MW2kW / unit_type_data[1, :d_x]
    xtr_exp = ones(unit_type_data[1, :d_x]) .* xtr_exp_per_pd

    # No construction expenditures from the end of construction until the end
    #   of the model's forecast period
    tail_zeros = zeros(fc_pd - lag - unit_type_data[1, :d_x])

    # Concatenate the above series into one
    xtr_exp_column = vcat(head_zeros, xtr_exp, tail_zeros)

    return xtr_exp_column
end


"""
    set_initial_debt_principal_series(unit_fs, unit_type_data, lag, agent_params)

Set up the record of the accrual of debt during construction.
"""
function set_initial_debt_principal_series(unit_fs, unit_type_data, lag, agent_params)
    for i = lag+1:lag+unit_type_data[1, :d_x]
        unit_fs[i, :remaining_debt_principal] = sum(unit_fs[1:i, :xtr_exp]) * agent_params[1, :debt_fraction]
    end
end




end
