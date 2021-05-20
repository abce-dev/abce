module ABCEfunctions

using SQLite, DataFrames, CSV

export load_db, get_current_period, get_agent_id, get_agent_params, load_unit_type_data, load_demand_data, set_forecast_period, extrapolate_demand, allocate_fuel_costs, create_unit_FS_dict, get_unit_specs, get_table, show_table, get_WIP_projects_list, get_demand_forecast, get_net_demand, get_next_asset_id, ensure_projects_not_empty, authorize_anpe, add_xtr_events

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


function load_demand_data(demand_data_file)
    demand_data = CSV.read(demand_data_file, DataFrame)
    return demand_data
end


function set_forecast_period(df, num_lags)
    transform!(df, [:d_x, :unit_life] => ((lead_time, unit_life) -> lead_time + unit_life) => :full_life)
    max_horizon = maximum(df[!, :full_life]) + num_lags
    return max_horizon
end


function extrapolate_demand(available_demand, fc_pd)
    demand = DataFrame(demand = zeros(Float64, convert(Int64, fc_pd)))
    demand[1:size(available_demand)[1], :demand] .= available_demand[!, :demand]
    demand[(size(available_demand)[1] + 1):fc_pd, :demand] .= demand[size(available_demand)[1], :demand] 
    return demand
end


function allocate_fuel_costs(unit_data, fuel_costs)
    num_units = size(unit_data)[1]
    unit_data[!, :uc_fuel] = zeros(num_units)
    for i = 1:num_units
        unit_data[i, :uc_fuel] = fuel_costs[fuel_costs[!, :fuel_type] == unit_data[i, :fuel_type], :cost_per_mmbtu]
    end
    return unit_data
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


function get_demand_forecast(db, pd, demand_vis_horizon, agent_id, fc_pd)
    # Get a list of forecasted future demand amounts
    # Forecast no increase after the end of the future visibility window
    # Hardcoded visibility window of 5
    vals = (pd, pd + demand_vis_horizon)
    demand_forecast = DBInterface.execute(db, "SELECT demand FROM demand WHERE period >= ? AND period < ?", vals) |> DataFrame
    println(demand_forecast)
    demand_forecast = extrapolate_demand(demand_forecast, fc_pd)
    return demand_forecast
end


function get_net_demand(db, pd, agent_id, fc_pd, demand_forecast)
    # Calculate the amount of forecasted net demand in future periods
    installed_cap_forecast = DataFrame(period = Int64[], derated_capacity = Float64[])
    vals = (pd, pd)
    # Select a list of all current assets, which are not cancelled, retired, or hidden from public view
    current_assets = DBInterface.execute(db, "SELECT * FROM assets WHERE cancellation_pd > ? AND retirement_pd > ? AND revealed = 'true'", vals) |> DataFrame
    println("Current assets:")
    println(current_assets)
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
# NPV transformation functions
#####

function add_xtr_events(unit_data, unit_num, unit_FS_dict, agent_params)
    # Generate the events which occur during the construction period:
    #    construction expenditures and the accrual of debt
    unit_num = convert(Int64, unit_num)
    for i = 1:unit_data[unit_num, :d_x]
        # Linearly distribute construction costs over the construction duration
        unit_FS_dict[i, :xtr_exp] = unit_data[unit_num, :uc_x] * unit_data[unit_num, :capacity] * 1000 / unit_data[unit_num, :d_x]
        unit_FS_dict[i, :remaining_debt_principal] = sum(fs[j, :xtr_exp] for j in 1:i) * agent_params[:debt_fraction]
    end
end


#function generate_operation_prime_movers(unit_data, unit_num, unit_FS_dict, agent_params)
    # Generate the prime-mover events which occur during the plant's operation:
    #   capital repayment, total amount of generation, and depreciation
#    for i = (unit_data[unit_num, :d_x] + 1):(unit_data[unit_num, :d_x] + unit_data[unit_num, :unit_life])
        # Apply a constant debt payment (sinking fund at the cost of debt)
#    end
#end





















end
