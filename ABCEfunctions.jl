module ABCEfunctions

using SQLite, DataFrames

export load_db, get_current_period, get_agent_id, get_table, show_table, get_active_projects_list, ensure_projects_not_empty, authorize_anpe

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

function get_active_projects_list(db, agent_id)
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
            asset_vals = (new_asset_id, string(agent_id), "no", "no", 9999, 0)
            DBInterface.execute(db, "INSERT INTO xtr_projects VALUES (?, ?, ?, ?, ?, ?)", xtr_vals)
            DBInterface.execute(db, "INSERT INTO assets VALUES (?, ?, ?, ?, ?, ?)", asset_vals)
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

function authorize_anpe(db, agent_id, current_period, project_list)
    # Loop through each project and authorize $100 of ANPE by setting the anpe value in xtr_projects
    for i = 1:size(project_list[!, :asset_id])[1]
        current_asset = project_list[i, :asset_id]
        println("Authorizing expenditures for project ", current_asset)
        vals = (100, current_period, current_asset)
        DBInterface.execute(db, "UPDATE xtr_projects SET anpe = ? WHERE period = ? AND asset_id = ?", vals)
    end
end

end
