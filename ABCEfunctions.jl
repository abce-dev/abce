module ABCEfunctions

using SQLite, DataFrames

export load_db, get_current_period, get_agent_id, get_table, show_table, get_active_projects_list, ensure_projects_not_empty, authorize_anpe

#####
# Setup functions
#####

function load_db()
    try
        db = SQLite.DB(ARGS[1])
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


function load_demand_data(demand_data_file):
    demand_data = CSV.read(demand_data_file, DataFrame)
    return demand_data
end


function set_forecast_period(df)
    transform!(df, [:d_x, :unit_life] => ((lead_time, unit_life) -> lead_time + unit_life) => :full_life)
    max_horizon = max(df[!, :full_life])
    return max_horizon
end


function forecast_demand(available_demand, fc_pd)
    demand = zeros(Float64, fc_pd)
    demand[1:size(available_demand)] = available_demand[!, :demand]
    return demand
end


function allocate_fuel_costs(unit_data, fuel_costs):
    num_units = size(unit_data)[1]
    unit_data[!, :uc_fuel] = zeros(num_units)
    for i = 1:num_units
        unit_data[i, :uc_fuel] = c_fuel[df[i, :fuel_type]]
    end
    return unit_data
end


function create_unit_FS_dict(unit_data)
    fs_dict = Dict{String, DataFrame}
    num_types = size(unit_data)[1]
    for i = 1:num_types:
        unit_name = unit_data[i, :name]
        unit_FS = DataFrame(year = 1:fc_pd, xtr_exp = zeros(fc_pd), gen = zeros(fc_pd), remaining_debt_principal = zeros(fc_pd), debt_payment = zeros(fc_pd), interest_due = zeros(fc_pd), depreciation = zeros(fc_pd))
        fs_dict[unit_name] = unit_FS
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


function get_WIP_projects_list(db, agent_id)
    # Get a list of all active (non-complete, non-cancelled) projects for the given agent
    SQL_get_proj = SQLite.Stmt(db, string("SELECT asset_id FROM assets WHERE agent_id = ", agent_id, " AND is_complete = 'no' AND is_cancelled = 'no'"))
    project_list = DBInterface.execute(SQL_get_proj) |> DataFrame
    return project_list
end


function get_next_available_id(db)
    # Return the next available asset ID (one greater than the current largest ID)
    SQL_get_ids = SQLite.Stmt(db, string("SELECT asset_id FROM assets"))
    asset_df = DBInterface.execute(SQL_get_ids) |> DataFrame
    asset_df[!, :asset_id] = tryparse.(Int64, asset_df[:, :asset_id])
    next_id = maximum(asset_df[!, :asset_id]) + 1
    return next_id    
end


function ensure_projects_not_empty(db, agent_id, project_list, current_period)
    # FAKE: only exists to ensure the Julia and Python scripts actually have something to do
    try
        if (size(project_list)[1] == 0) && (current_period <= 5)
            # Current period restriction is a testing spoof only to allow
            #    observation of post-completion model behavior

            # Create a new project
            new_asset_id = get_next_available_id(db)

            # Assign some dummy data to the project
            xtr_vals = (new_asset_id, string(agent_id), 0, 1000, 10, 0)
            asset_vals = (new_asset_id, string(agent_id), "gas", "no", "no", 9999, 0)
            DBInterface.execute(db, "INSERT INTO xtr_projects VALUES (?, ?, ?, ?, ?, ?)", xtr_vals)
            DBInterface.execute(db, "INSERT INTO assets VALUES (?, ?, ?, ?, ?, ?, ?)", asset_vals)
            println(string("Created project ", new_asset_id))

            # Update the list of active construction projects and return it
            project_list = get_active_projects_list(db, agent_id)
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


function authorize_anpe(db, agent_id, current_period, project_list, unit_data)
    # Loop through each project and authorize $100 of ANPE by setting the anpe value in xtr_projects
    for i = 1:size(project_list[!, :asset_id])[1]
        current_asset = project_list[i, :asset_id]
        println("Authorizing expenditures for project ", current_asset)
        anpe_val = 1000000000   # $1B/period
        vals = (anpe_val, current_period, current_asset)
        DBInterface.execute(db, "UPDATE xtr_projects SET anpe = ? WHERE period = ? AND asset_id = ?", vals)
    end
end


#####
# NPV transformation functions
#####

function add_xtr_events(unit_data, unit_num, unit_FS_dict, agent_params)
    # Generate the events which occur during the construction period:
    #    construction expenditures and the accrual of debt
    for i = 1:unit_data[unit_num, :d_x]
        # Linearly distribute construction costs over the construction duration
        unit_FS_dict[i, :xtr_exp] = unit_data[unit_num, :uc_x], * unit_data[unit_num, :capacity] * 1000 / unit_data[unit_num, :d_x]
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
