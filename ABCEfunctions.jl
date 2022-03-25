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

using SQLite, DataFrames, CSV, JuMP, GLPK, Logging, CPLEX

export load_db, get_current_period, get_agent_id, get_agent_params, load_unit_type_data, set_forecast_period, extrapolate_demand, project_demand_flat, project_demand_exponential, allocate_fuel_costs, create_FS_dict, get_unit_specs, get_table, show_table, get_WIP_projects_list, get_demand_forecast, get_net_demand, get_next_asset_id, ensure_projects_not_empty, authorize_anpe, create_NPV_results_df, generate_xtr_exp_profile, set_initial_debt_principal_series, generate_prime_movers, forecast_unit_revenue_and_gen, forecast_unit_op_costs, propagate_accounting_line_items, compute_alternative_NPV, set_up_model, get_current_assets_list, convert_to_marginal_delta_FS, run_scenario_dispatches

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
    @info string("\nTable \'", table_name, "\':")
    @info df
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
    SQL_get_assets = SQLite.Stmt(db, string("SELECT asset_id, unit_type, retirement_pd FROM assets WHERE agent_id = ", agent_id, " AND completion_pd <= ", pd, " AND cancellation_pd > ", pd, " AND retirement_pd > ", pd))
    asset_list = DBInterface.execute(SQL_get_assets) |> DataFrame

    # Count the number of assets by type
    asset_counts = combine(groupby(asset_list, [:unit_type, :retirement_pd]), nrow => :count)

    # Bugfix: ensure all `retirement_pd` values are integers
    asset_counts[!, :retirement_pd] = convert.(Int64, asset_counts[!, :retirement_pd])

    @info asset_counts

    return asset_list, asset_counts
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
        @error string("The specified demand extrapolation mode, ", mode, ", is not implemented.")
        @error "Please use 'flat', 'exp_termrate', or 'exp_fitted' at this time."
        @error "Terminating..."
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
    current_assets = DBInterface.execute(db, "SELECT * FROM assets WHERE cancellation_pd > ? AND retirement_pd > ?", vals) |> DataFrame
    if size(current_assets)[1] == 0
        @warn "There are no currently-active generation assets in the system; unpredictable behavior may occur."
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
    return installed_cap_forecast, net_demand_forecast
end


function get_next_asset_id(db)
    # Return the next available asset ID (one greater than the current largest ID)
    asset_max = DBInterface.execute(db, "SELECT MAX(asset_id) FROM assets") |> DataFrame
    WIP_max = DBInterface.execute(db, "SELECT MAX(asset_id) FROM WIP_projects") |> DataFrame
    pending_asset_max = DBInterface.execute(db, "SELECT MAX(asset_id) FROM asset_updates") |> DataFrame
    pending_WIP_max = DBInterface.execute(db, "SELECT MAX(asset_id) FROM WIP_updates") |> DataFrame
    results = skipmissing([asset_max[1, 1], WIP_max[1, 1], pending_asset_max[1, 1], pending_WIP_max[1, 1]])
    next_id = maximum(results) + 1

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
            @info string("Created project ", new_asset_id)

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
    # ANPE: Authorized Next-Period Expenditures
    # Loop through each project and authorize ANPE by setting the anpe value in WIP_projects
    for i = 1:size(project_list[!, :asset_id])[1]
        current_asset = project_list[i, :asset_id]
        asset_type = DBInterface.execute(db, string("SELECT unit_type FROM assets WHERE asset_id = ", current_asset)) |> DataFrame
        unit = filter(row -> row[:unit_type] == asset_type[1, :unit_type], unit_specs)
        project_data = DBInterface.execute(db, string("SELECT * FROM WIP_projects WHERE asset_id = ", current_asset, " AND period = ", current_period)) |> DataFrame
        # Authorize a uniform expenditure over the life of the project
        # Maximum productive investment rate: the yearly expenditure required
        #   to complete the project to its nominal cost in the nominal
        #   construction duration
        # Minimum investment: $0 (no negative ANPE values allowed)
        max_invest_rate = unit[1, :uc_x] * unit[1, :capacity] * MW2kW / unit[1, :d_x]
        anpe_val = min(max_invest_rate, max(project_data[1, :rcec], 0))

        vals = (anpe_val, current_period, current_asset)
        DBInterface.execute(db, "UPDATE WIP_projects SET anpe = ? WHERE period = ? AND asset_id = ?", vals)
    end
end


#####
# NPV functions
#####

function check_valid_vector_mode(mode)
    if !(mode in ["new_xtr", "retire"])
        # Invalid mode supplied; alert the user and exit
        @error string("Invalid decision vector type specified: ", mode)
        @error "Please ensure that 'mode' is set to either 'new_xtr' or 'retire'."
        exit()
    end
end



"""
    create_NPV_results_DF(unit_data_df, num_lags; mode="new_xtr")

Create a dataframe to hold the results of NPV calculations for the various
  types, expanded by the number of allowed lags.

Arguments:
  unit_data_df (DataFrame): unit data, depending on mode specification:
    new_xtr: unit specification data from DB
    retire: asset counts pivot table (agent's # of existing assets by type
        and retirement pd)

  num_lags (int): number of lags the agent considers (spec'd in settings.yml)

  mode options:
    new_xtr: for new construction, uses unit specification as input
    retire: for retiring existing assets, uses asset counts from DB as input
"""
function create_NPV_results_df(unit_data_df, num_lags; mode="new_xtr")
    check_valid_vector_mode(mode)

    alternative_names = Vector{String}()
    num_entries = size(unit_data_df)[1]
    num_alternatives = num_entries * (num_lags + 1)

    for i = 1:num_entries
        for j = 0:num_lags
            if mode == "new_xtr"
                name = string(unit_data_df[i, :unit_type], "_0_lag-", j)
            elseif mode == "retire"
                name = string(unit_data_df[i, :unit_type], "_", unit_data_df[i, :retirement_pd], "_lag-", j)
            end
            push!(alternative_names, name)
        end
    end

    NPV_results = DataFrame(name=alternative_names, NPV=zeros(num_alternatives))
    return alternative_names, NPV_results

end


function create_FS_dict(data, fc_pd, num_lags; mode="new_xtr")
    check_valid_vector_mode(mode)

    fs_dict = Dict()
    num_alts = size(data)[1]
    for i = 1:num_alts
        for j = 0:num_lags
            if mode == "retire"
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
            unit_FS = DataFrame(year = 1:fc_pd, xtr_exp = zeros(fc_pd), gen = zeros(fc_pd), remaining_debt_principal = zeros(fc_pd), debt_payment = zeros(fc_pd), interest_payment = zeros(fc_pd), depreciation = zeros(fc_pd))
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
    unit_d_x = unit_type_data[1, :d_x]
    unit_op_life = unit_type_data[1, :unit_life]
    for i = (lag + unit_d_x + 1):(lag + unit_d_x + unit_op_life)
        # Apply a constant debt payment (sinking fund at cost of debt), based
        #   on the amount of debt outstanding at the end of the xtr project
        unit_fs[i, :debt_payment] = unit_fs[lag + unit_d_x, :remaining_debt_principal] .* cod ./ (1 - (1+cod) .^ (-1*unit_op_life))

        # Determine the portion of each payment which pays down interest
        #   (instead of principal)
        unit_fs[i, :interest_payment] = unit_fs[i-1, :remaining_debt_principal] * cod

        # Update the amount of principal remaining at the end of the period
        unit_fs[i, :remaining_debt_principal] = unit_fs[i-1, :remaining_debt_principal] - (unit_fs[i, :debt_payment] - unit_fs[i, :interest_payment])

        # Apply straight-line depreciation, based on debt outstanding at the
        #   project's completion
        unit_fs[i, :depreciation] = unit_fs[lag + unit_d_x, :xtr_exp] ./ unit_op_life
    end

end


"""
    forecast_unit_revenue_and_gen(unit_type_data, unit_fs, price_curve, db, pd, lag; mode, ret_pd)

Forecast the revenue which will be earned by the current unit type,
and the total generation (in kWh) of the unit, per time period.

This function assumes the unit will be a pure price-taker: it will generate
  during all hours for which it is a marginal or submarginal unit.
If the unit is of a VRE type (as specified in the A-LEAF inputs), then a flat
  de-rating factor is applied to its availability during hours when it is
  eligible to generate.
"""
function forecast_unit_revenue_and_gen(unit_type_data, unit_fs, price_curve, db, pd, lag; mode="new_xtr", orig_ret_pd=9999)
    check_valid_vector_mode(mode)

    # Compute estimated revenue from submarginal hours
    num_submarg_hours, submarginal_hours_revenue = compute_submarginal_hours_revenue(unit_type_data, price_curve)

    # Compute estimated revenue from marginal hours
    num_marg_hours, marginal_hours_revenue = compute_marginal_hours_revenue(unit_type_data, price_curve, db, pd)

    # If the unit is VRE, assign an appropriate availability derate factor
    availability_derate_factor = compute_VRE_derate_factor(unit_type_data)

    # Compute the original retirement period
    # Minimum of ret_pd or size of the unit FS (i.e. the forecast period)
    ret_pd = min(size(unit_fs)[1], orig_ret_pd)

    # Compute total projected revenue, with VRE adjustment if appropriate, and
    #   save to the unit financial statement
    compute_total_revenue(unit_type_data, unit_fs, submarginal_hours_revenue, marginal_hours_revenue, availability_derate_factor, lag; mode=mode, orig_ret_pd=ret_pd)

    # Compute the unit's total generation for each period, in kWh
    compute_total_generation(unit_type_data, unit_fs, num_submarg_hours, num_marg_hours, availability_derate_factor, lag; mode=mode, orig_ret_pd=ret_pd)

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
    submarginal_hours = filter(row -> row.lamda > unit_type_data[1, :VOM] + unit_type_data[1, :FC_per_MWh], price_curve)

    # Create a conversion factor from number of periods in the price curve
    #   to hours (price curve may be in hours or five-minute periods)
    convert_to_hours = hours_per_year / size(price_curve)[1]

    # Compute total number of submarginal hours
    num_submarg_hours = size(submarginal_hours)[1] * convert_to_hours

    # Calculate total revenue from submarginal hours
    submarginal_hours_revenue = sum(submarginal_hours[!, :lamda] * unit_type_data[1, :capacity]) * convert_to_hours

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
    marginal_hours = filter(row -> row.lamda == unit_type_data[1, :VOM] + unit_type_data[1, :FC_per_MWh], price_curve)

    # Compute the total capacity of this unit type in the current system
    command = string("SELECT asset_id FROM assets WHERE unit_type == '", unit_type_data[1, :unit_type], "' AND completion_pd <= ", pd, " AND cancellation_pd > ", pd, " AND retirement_pd > ", pd)
    system_type_list = DBInterface.execute(db, command) |> DataFrame
    system_type_capacity = size(system_type_list)[1]

    # Create a conversion factor from number of periods in the price curve
    #   to hours (price curve may be in hours or five-minute periods)
    convert_to_hours = hours_per_year / size(price_curve)[1]

    # Compute effective number of marginal hours
    num_marg_hours = size(marginal_hours)[1] * unit_type_data[1, :capacity] / system_type_capacity * convert_to_hours

    if size(marginal_hours)[1] == 0
        marginal_hours_revenue = 0
    else
        marginal_hours_revenue = sum(marginal_hours[!, :lamda]) * unit_type_data[1, :capacity] / system_type_capacity * convert_to_hours
    end

    return num_marg_hours, marginal_hours_revenue
end


"""
    compute_VRE_derate_factor(unit_type_data)

If the unit is declared as `is_VRE` in the A-LEAF inputs, return an
availability derating factor equal to its capacity credit factor.
"""
function compute_VRE_derate_factor(unit_type_data)
    availability_derate_factor = 1

    # Julia reads booleans as UInt8, so have to convert to something sensible
    if convert(Int64, unit_type_data[1, "is_VRE"][1]) == 1
        availability_derate_factor = unit_type_data[1, :CF]
    end

    return availability_derate_factor   
end



"""
    compute_total_revenue(unit_type_data, unit_fs, submarginal_hours_revenue, marginal_hours_revenue, availability_derate_factor, lag; mode, orig_ret_pd)

Compute the final projected revenue stream for the current unit type, adjusting
unit availability if it is a VRE type.
"""
function compute_total_revenue(unit_type_data, unit_fs, submarginal_hours_revenue, marginal_hours_revenue, availability_derate_factor, lag; mode, orig_ret_pd=9999)
    check_valid_vector_mode(mode)

    # Helpful short variables
    unit_d_x = unit_type_data[1, :d_x]
    unit_op_life = unit_type_data[1, :unit_life]

    # Add a Revenue column to the financial statement dataframe
    unit_fs[!, :Revenue] .= 0.0

    # In "new_xtr" mode, revenues START accruing after the lag plus
    #   construction duration
    # In "retire" mode, revenues start in the current period and CEASE accruing
    #   after the soonest of {the lag, the mandatory retirement period, the end
    #   of the forecast period}.
    if mode == "new_xtr"
        rev_start = lag + unit_d_x + 1
        rev_end = lag + unit_d_x + unit_op_life
    elseif mode == "retire"
        rev_start = 1
        # The maximum end of revenue period is capped at the length of the
        #   unit_fs dataframe (as some units with unspecified retirement
        #   periods default to a retirement period of 9999).
        rev_end = min(orig_ret_pd, size(unit_fs)[1])
    end

    # Compute final projected revenue series
    unit_fs[rev_start:rev_end, :Revenue] .= (submarginal_hours_revenue + marginal_hours_revenue) * availability_derate_factor

end



"""
    compute_total_generation(unit_type_data, unit_fs, availability_derate_factor, lag; mode, orig_ret_pd)

Calculate the unit's total generation for the period, in kWh.
"""
function compute_total_generation(unit_type_data, unit_fs, num_submarg_hours, num_marg_hours, availability_derate_factor, lag; mode, orig_ret_pd=9999)
    check_valid_vector_mode(mode)

    # Helpful short variable names
    unit_d_x = unit_type_data[1, :d_x]
    unit_op_life = unit_type_data[1, :unit_life]

    # Compute total generation
    gen = (num_submarg_hours + num_marg_hours) * unit_type_data[1, :capacity] * availability_derate_factor * MW2kW   # in kWh

    # In "new_xtr" mode, generation persists from end of construction to end of life
    # In "retire" mode, generation occurs from the current period to the
    #   soonest of {the lag, the mandatory retirement period, the end of the
    #   forecast period}.
    if mode == "new_xtr"
        gen_start = lag + unit_d_x + 1
        gen_end = lag + unit_d_x + unit_op_life
    elseif mode == "retire"
        gen_start = 1
        # The maximum end of generation period is capped at the length of the
        #   unit_fs dataframe (as some units with unspecified retirement
        #   periods default to a retirement period of 9999).
        gen_end = min(orig_ret_pd, size(unit_fs)[1])
    end

    # Distribute generation values time series
    unit_fs[gen_start:gen_end, :gen] .= gen
end


"""
    forecast_unit_op_costs(unit_type_data, unit_fs, lag; mode, orig_ret_pd)

Forecast cost line items for the current unit:
 - fuel cost
 - VOM
 - FOM
"""
function forecast_unit_op_costs(unit_type_data, unit_fs, lag; mode="new_xtr", orig_ret_pd=9999)
    check_valid_vector_mode(mode)

    # Helpful short variable names
    unit_d_x = unit_type_data[1, :d_x]
    unit_op_life = unit_type_data[1, :unit_life]

    # Compute total fuel cost
    # Unit conversions:
    #  gen [kWh/year] * FC_per_MWh [$/MWh] * [1 MWh / 1000 kWh] = $/year
    transform!(unit_fs, [:gen] => ((gen) -> gen .* unit_type_data[1, :FC_per_MWh] ./ MW2kW) => :Fuel_Cost)

    # Compute total VOM cost incurred during generation
    # Unit conversions:
    #   gen [kWh/year] * VOM [$/MWh] * [1 MWh / 1000 kWh] = $/year
    transform!(unit_fs, [:gen] => ((gen) -> gen .* unit_type_data[1, :VOM] ./ MW2kW) => :VOM_Cost)

    # Compute total FOM cost for each period
    # Unit conversions:
    #   FOM [$/kW-year] * capacity [MW] * [1000 kW / 1 MW] = $/year
    if mode == "new_xtr"
        pre_zeros = zeros(lag + unit_d_x)
        op_ones = ones(unit_op_life)
    elseif mode == "retire"
        pre_zeros = zeros(0)
        op_ones = ones(min(unit_op_life, size(unit_fs)[1], orig_ret_pd))
    end
    post_zeros = zeros(size(unit_fs)[1] - size(pre_zeros)[1] - size(op_ones)[1])
    unit_fs[!, :FOM_Cost] = vcat(pre_zeros, op_ones, post_zeros)

    unit_fs[!, :FOM_Cost] .= unit_fs[!, :FOM_Cost] .* unit_type_data[1, :FOM] .* unit_type_data[1, :capacity] .* MW2kW

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
    transform!(unit_fs, [:Net_Income, :interest_payment, :xtr_exp] => ((NI, interest, xtr_exp) -> NI + interest - xtr_exp) => :FCF)

end


"""
    convert_to_marginal_delta_FS(unit_fs, lag)

For this project alternative, compute the delta between the original projected
FCF (without premature retirement) and the new projected FCF (with premature
retirement).
"""
function convert_to_marginal_delta_FS(unit_fs, lag; mode="new_xtr")
    cols_to_invert = [col for col in names(unit_fs) if col != "year"]
    for col in cols_to_invert
        # Set values to 0 during the lag (no change)
        unit_fs[1:lag, col] .= 0
        # Invert values post-retirement to indicate lost revenues/generation
        #   and avoided costs
        unit_fs[!, col] .= -unit_fs[!, col]
    end

end


"""
    compute_alternative_NPV(unit_fs, agent_params)

For this project alternative (unit type + lag time), compute the project's FCF NPV.
"""
function compute_alternative_NPV(unit_fs, agent_params)
    # Discount rate is WACC
    d = agent_params[1, :debt_fraction] * agent_params[1, :cost_of_debt] + (1 - agent_params[1, :debt_fraction]) * agent_params[1, :cost_of_equity]

    # Add a column of compounded discount factors to the dataframe
    transform!(unit_fs, [:year] => ((year) -> (1+d) .^ (-1 .* (year .- 1))) => :discount_factor)

    # Discount the alternative's FCF NPV
    FCF_NPV = transpose(unit_fs[!, :FCF]) * unit_fs[!, :discount_factor]

    return FCF_NPV

end


### JuMP optimization model initialization
"""
    set_up_model(unit_FS_dict, ret_fs_dict, fc_pd, available_demand, NPV_results, ret_NPV_results)

Set up the JuMP optimization model, including variables, constraints, and the
objective function.

Returns:
  m (JuMP model object)
"""
function set_up_model(settings, agent_params, unit_FS_dict, ret_FS_dict, available_demand, new_xtr_NPV_df, ret_NPV_df, asset_counts, unit_specs, installed_capacity_forecast)
    # Create the model object
    @info "Setting up model..."
    m = Model(GLPK.Optimizer)

    # For debugging, enable the following line to increase verbosity
    set_optimizer_attribute(m, "msg_lev", GLPK.GLP_MSG_ALL)

    # Set the MIP optimality gap parameter to improve solve times
    set_optimizer_attribute(m, "mip_gap", 0.001)

    # Concatenate all results into unified data structures
    all_FS_dict = merge(unit_FS_dict, ret_FS_dict)
    all_NPV_results = vcat(new_xtr_NPV_df, ret_NPV_df)

    # Parameter names
    alternative_names = [item for item in keys(all_FS_dict)]
    num_alternatives = size(all_NPV_results)[1]
    num_time_periods = size(all_FS_dict[alternative_names[1]])[1]

    # Set up variables
    # Number of units of each type to build: must be Integer
    @variable(m, u[1:num_alternatives] >= 0, Int)

    # To prevent unnecessary infeasibility conditions, convert nonpositive
    #   available_demand values to 0
    nonnegative_available_demand = deepcopy(available_demand)
    excess_supply = deepcopy(available_demand) * (-1)
    for i = 1:size(available_demand)[1]
        if nonnegative_available_demand[i] < 0
            nonnegative_available_demand[i] = 0
        end
    end

    # Compute expected marginal generation contribution per alternative type
    marg_gen = zeros(num_alternatives, num_time_periods)
    for i = 1:num_alternatives
        alt_name = string(all_NPV_results[i, :unit_type], "_",
                          all_NPV_results[i, :project_type], "_",
                          all_NPV_results[i, :retirement_pd], "_lag-",
                          all_NPV_results[i, :lag])
        for j = 1:num_time_periods
            # Convert total anticipated marginal generation to an effective
            #   nameplate capacity and save to the appropriate entry in the
            #   marginal generation array
            marg_gen[i, j] = all_FS_dict[alt_name][j, :gen] / (hours_per_year * MW2kW)
        end
    end

    # Restrict total construction to be less than maximum available demand,
    for i = 4:num_time_periods
        if excess_supply[i] > 0
            @constraint(m, transpose(u) * marg_gen[:, i] <= (-1) * excess_supply[i] * agent_params[1, "aggressiveness"] * 0.2)
#        elseif excess_supply[i] > 0
#            @constraint(m, transpose(u) * marg_gen[:, i] <= (-1) * excess_supply[i] * 0.33)
        end
    end

    # Prevent agents from causing foreseeable short-term supply scarcity
    protected_pd = 3
    for i = 1:protected_pd
        if excess_supply[i] < 0
            @constraint(m, transpose(u) * marg_gen[:, i] >= 0)
        elseif (excess_supply[i] >= 0) && (excess_supply[i] < installed_capacity_forecast[i, :derated_capacity] * 0.05)
            @constraint(m, transpose(u) * marg_gen[:, i] >= (-1) * excess_supply[i] * 0.25)
        elseif excess_supply[i] > installed_capacity_forecast[i, :derated_capacity] * 0.05
            @constraint(m, transpose(u) * marg_gen[:, i] >= (-1) * excess_supply[i] * 0.4)
        end
    end

    # Prevent agents from excessively overbuilding in the long term
#    long_term_start = 3
#    for i = long_term_start:num_time_periods
#        @constraint(m, transpose(u) * marg_gen[:, i] <= installed_capacity_forecast[i, :derated_capacity] * 0.1 - excess_supply[i])
#    end

    for i = 1:num_alternatives
        if all_NPV_results[i, :project_type] == "new_xtr"
            @constraint(m, u[i] .<= settings["max_type_newbuilds_per_pd"])
        elseif all_NPV_results[i, :project_type] == "retirement"
            unit_type = all_NPV_results[i, :unit_type]
            ret_pd = all_NPV_results[i, :retirement_pd]
            asset_count = filter([:unit_type, :retirement_pd] => (x, y) -> x == unit_type && y == ret_pd, asset_counts)[1, :count]
            max_retirement = min(asset_count, settings["max_type_rets_per_pd"])
            @constraint(m, u[i] .<= max_retirement)
        end
    end

    # Convenient variable for the number of distinct retirement alternatives
    k = size(asset_counts)[1]
    # Shortened name for the number of lag periods to consider
    #   1 is added, as the user-set value only specifies future periods,
    #   with the "lag = 0" instance being implied
    num_lags = settings["num_future_periods_considered"] + 1
    # Create the matrix to collapse lagged options
    ret_summation_matrix = zeros(k, k*num_lags)
    for i = 1:k
        ret_summation_matrix[i, ((i-1)*num_lags+1):(i*num_lags)] .= 1
    end

    # Specify constraint: the agent cannot plan to retire more units (during
    #   all lag periods) than exist of that unit type
    @constraint(m, ret_summation_matrix * u[size(new_xtr_NPV_df)[1]+1:end, :] .<= asset_counts[!, :count])
        

    # Create the objective function 
    @objective(m, Max, transpose(u) * all_NPV_results[!, :NPV])
    @info "Optimization model set up."

    return m
end


### DISPATCH FUNCTIONS

function initialize_repdays(settings)
    @info "Initializing representative days."
    # Read in all repdays data
    repdays_spec = CSV.read(settings["repdays_file"], DataFrame)
    load_ts_data = CSV.read(settings["load_ts_file"], DataFrame)
    wind_ts_data = CSV.read(settings["wind_ts_file"], DataFrame)
    solar_ts_data = CSV.read(settings["solar_ts_file"], DataFrame)

    # Set up unscaled representative days
    load_repdays = DataFrame()
    wind_repdays = DataFrame()
    solar_repdays = DataFrame()
    for day in repdays_spec[!, :Day]
        load_repdays = vcat(load_repdays, load_ts_data[24*day+1:24*(day+1), :])
        wind_repdays = vcat(wind_repdays, wind_ts_data[24*day+1:24*(day+1), :])
        solar_repdays = vcat(solar_repdays, solar_ts_data[24*day+1:24*(day+1), :])
    end

    return repdays_spec, (load_repdays, wind_repdays, solar_repdays)
end


function initialize_indices(db, repdays_spec)
    distinct_units = DBInterface.execute(db, "SELECT DISTINCT unit_type FROM unit_specs") |> DataFrame
    num_units = size(distinct_units)[1]
    num_days = size(repdays_spec)[1]
    num_hours = 24

    return (num_units, num_days, num_hours)
end


function scale_annual_repdays(y, peak_demand, repdays_data, portfolio, unit_specs)
    # Make copies of repdays data for this year
    year_load_repdays = deepcopy(repdays_data[1])
    year_wind_repdays = deepcopy(repdays_data[2])
    year_solar_repdays = deepcopy(repdays_data[3])

    # Scale the load data
    year_load_repdays[!, :Load] = year_load_repdays[!, :LoadShape] .* peak_demand

    # Scale the wind and solar availability data
    year_wind_repdays[!, :Available_Capacity] .= (
        year_wind_repdays[!, :WindShape]
        .* filter([:unit_type, :y] => (x, y) -> x == "Wind" && y == y, portfolio)[1, :num_units]
        .* filter(:unit_type => x -> x == "Wind", unit_specs)[1, :capacity]
    )

    year_solar_repdays[!, :Available_Capacity] .= (
        year_solar_repdays[!, :SolarShape]
        .* filter([:unit_type, :y] => (x, y) -> x == "Solar" && y == y, portfolio)[1, :num_units]
        .* filter(:unit_type => x -> x == "Solar", unit_specs)[1, :capacity]
    )

    return (year_load_repdays, year_wind_repdays, year_solar_repdays)

end


function set_up_annual_dispatch_model(indices, repdays_data, portfolio, unit_specs)
    @info "Setting up the annual dispatch model."
    num_units, num_days, num_hours = indices
    load_repdays, wind_repdays, solar_repdays = repdays_data
    # Initialize model
    m = Model(with_optimizer(CPLEX.Optimizer))

    # g: quantity generated (MWh) for each unit type
    @variable(m, g[1:num_units, 1:num_days, 1:num_hours] >= 0)

    # c: number of units of each type committed in each hour
    @variable(m, c[1:num_units, 1:num_days, 1:num_hours] >= 0, Int)

    # Total generation per hour must be greater than or equal to demand
    @constraint(
        m,
        mkt_equil[k=1:num_days, j=1:num_hours],
        sum(g[i, k, j] for i = 1:num_units) >= load_repdays[(k-1)*24+j, :Load]
    )

    # Number of committed units per hour must be less than or equal to the
    #   total number of units of that type in the system
    @constraint(
        m,
        commit_balance[i=1:num_units, k=1:num_days, j=1:num_hours],
        c[i, k, j] <= portfolio[i, :num_units]
    )

    # Limit total generation per unit type each hour to the total capacity of
    #   all committed units of this type with committed units subject to 
    #   minimum and maximum power levels
    @constraint(
        m,
        min_pl[i=3:num_units, k=1:num_days, j=1:num_days],
        (
            g[i, k, j]
                <= c[i, k, j] .* unit_specs[i, :capacity]
                .* unit_specs[i, :CF] .* unit_specs[i, :PMAX]
        )
    )

    @constraint(
        m,
        max_pl[i=3:num_units, k=1:num_days, j=1:num_days],
        (
            g[i, k, j]
                <= c[i, k, j] .* unit_specs[i, :capacity]
                .* unit_specs[i, :CF] .* unit_specs[i, :PMIN]
        )
    )

    # Limit solar and wind generation to their actual hourly availability
    @constraint(
        m,
        wind_gen[k=1:num_days, j=1:num_hours],
        g[1, k, j] <= wind_repdays[(k-1)*24+j, :Available_Capacity]
    )

    @constraint(
        m,
        solar_gen[k=1:num_days, j=1:num_hours],
        g[2, k, j] <= solar_repdays[(k-1)*24+j, :Available_Capacity]
    )

    # Ramping constraints
    @constraint(
        m,
        RUL[i=1:num_units, k=1:num_days, j=1:num_hours-1],
        (
            g[i, k, j+1] - g[i, k, j]
                <= (c[i, k, j+1] .* unit_specs[i, :RUL]
                .* unit_specs[i, :capacity] .* unit_specs[i, :CF])
        )
    )

    @constraint(
        m,
        RDL[i=1:num_units, k=1:num_days, j=1:num_hours-1],
        (
            g[i, k, j+1] - g[i, k, j]
                >= (-1) * c[i, k, j+1] * unit_specs[i, :RDL]
                * unit_specs[i, :capacity] * unit_specs[i, :CF]
        )
    )

    # Set up objective: minimize cost of service
    @objective(
        m,
        Min,
        (
            sum(sum(sum(g[i, k, j] for j = 1:num_hours) for k = 1:num_days)
            .* (unit_specs[i, :VOM] + unit_specs[i, :FC_per_MWh])
            for i = 1:num_units)
        )
    )

    return m
end


function solve_annual_dispatch(m, y, indices, unit_specs, all_gc_results, all_price_results)
    # Unpack indices
    num_units, num_days, num_hours = indices

    # Make a copy of m in order to solve the LP version later on, to
    #   get shadow prices
    m_copy = copy(m)
    set_optimizer(m_copy, CPLEX.Optimizer)

    @info "Solving MILP problem..."
    solve_annual_MILP_dispatch(m, y, indices, unit_specs, all_gc_results)
    @info "Solving LP problem..."
    solve_annual_LP_dispatch(m_copy, y, indices, unit_specs, all_price_results)
    @info "Dispatch problems solved; all data saved."

end

function solve_annual_MILP_dispatch(m, y, indices, unit_specs, all_gc_results)
    # Unpack indices
    num_units, num_days, num_hours = indices

    # Solve the optimization problem
    optimize!(m)
    status = termination_status.(m)
    @info "MILP status: $status"

    # Retrieve the generation and commitment time-series
    gen_qty = value.(m[:g])
    c = value.(m[:c])

    # Save generation and commitment results
    for k = 1:num_days
        for j = 1:num_hours
            for i = 1:num_units
                line = (
                    y = y,
                    d = k,
                    h = j,
                    unit_type = unit_specs[i, :unit_type],
                    gen = gen_qty[i, k, j],
                    commit = c[i, k, j]
                )
                push!(all_gc_results, line)
            end
        end
    end

end

function solve_annual_LP_dispatch(m_copy, y, indices, unit_specs, all_price_results)
    # Unpack indices
    num_units, num_days, num_hours = indices

    # Relax the integrality constraints in m_copy (i.e. c[...]) to allow
    #   retrieval of the dual solution containing market prices
    undo = relax_integrality(m_copy)

    # Solve the LP
    optimize!(m_copy)
    status = termination_status.(m_copy)
    @info "LP status: $status"

    # Retrieve the electricity price data
    mkt_px = reshape(shadow_price.(m_copy[:mkt_equil]), (num_days, num_hours))

    # Save price results
    for j = 1:num_hours
        for k = 1:num_days
            line = (
                y = y,
                d = k,
                h = j,
                price = (-1) * mkt_px[k, j]
            )
            push!(all_price_results, line)
        end
    end

end


function handle_annual_dispatch(y, scenario, indices, repdays_data, unit_specs, all_gc_results, all_price_results)
    # Set up year-specific data by selecting from scenario or calculating
    @info "Setting up data for year $y..."
    peak_demand = scenario.PD_series[y]
    year_portfolio = filter(row -> (row[:y] == y), scenario.portfolio)
    year_repdays_data = scale_annual_repdays(y, peak_demand, repdays_data, year_portfolio, unit_specs)
    @info "Data set up."

    m = set_up_annual_dispatch_model(indices, year_repdays_data, year_portfolio, unit_specs)

    # Solve the MILP and LP dispatch problems for this year
    solve_annual_dispatch(m, y, indices, unit_specs, all_gc_results, all_price_results)
end


function compute_financial_results(all_gc_results, all_price_results, unit_specs)
    # Make a copy of all_gc_results for manipulation
    unit_revcost_data = deepcopy(all_gc_results)

    # Inner join with price data to enable straightforward row-wise functions
    unit_revcost_data = innerjoin(unit_revcost_data, all_price_results, on = [:y, :d, :h])

    # Retrieve a reduced set of unit specs and inner join onto the financial
    #   data dataframe
    reduced_unit_specs = select(unit_specs, [:unit_type, :capacity, :FOM, :VOM, :FC_per_MWh, :CF])
    unit_revcost_data = innerjoin(unit_revcost_data, reduced_unit_specs, on = :unit_type)

    # Compute revenue
    transform!(unit_revcost_data, [:gen, :price] => ((gen, price) -> gen .* price) => :revenue)

    # Compute total variable operating cost
    transform!(unit_revcost_data, [:gen, :VOM, :FC_per_MWh] => ((gen, VOM, FC) -> gen .* (VOM + FC)) => :variable_opcost)

    return unit_revcost_data

end


function run_scenario_dispatches(scenario, db, settings, unit_specs)
    # Initialize the representative days
    repdays_spec, repdays_data = initialize_repdays(settings)

    # Initialize standard indices
    indices = initialize_indices(db, repdays_spec)

    # Initialize results dataframes
    all_gc_results = DataFrame(y = Int[], d = Int[], h = Int[], unit_type = String[], gen = Float64[], commit = Int[])
    all_price_results = DataFrame(y = Int[], d = Int[], h = Int[], price = Float64[])

    # Run dispatch for all years in this scenario
    for y = 1:size(scenario.PD_series)[1]
        handle_annual_dispatch(y, scenario, indices, repdays_data, unit_specs, all_gc_results, all_price_results)
    end

    unit_revcost_data = compute_financial_results(all_gc_results, all_price_results, unit_specs)

    return unit_revcost_data

end


end


