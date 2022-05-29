module C2N

using Logging, CSV, DataFrames

mutable struct C2N_project
    conv_type::String
    rxtr_type::String
    lag::Int64
    uid::Int64
    fs::DataFrame
end


function C2N_project(conv_type, rxtr_type, base_pd, lag, uid, fc_pd)
    # Validate inputs
    check_init_inputs(conv_type, lag, uid, fc_pd)

    capex_tl = project_capex_profile(base_pd, lag, db, fc_pd, C2N_PA, unit_specs; status="new", asset_id=nothing)


    fs = DataFrame(
        year = 1:fc_pd,
        capex = capex_tl[!, :capex],
        remaining_debt_principal = zeros(fc_pd),
        debt_payment = zeros(fc_pd),
        interest_payment = zeros(fc_pd),
        depreciation = zeros(fc_pd)
    )

    return C2N_project(C2N_type, lag, uid, fs)
end


function create_C2N_fs(db, conv_type, rxtr_type, base_pd, lag, fc_pd, C2N_specs)
    check_init_inputs(conv_type, rxtr_type, lag, fc_pd)
    capex_tl = project_capex_profile(base_pd, lag, db, fc_pd, rxtr_type, conv_type, C2N_specs; status="new", asset_id=nothing)

end


function check_init_inputs(conv_type, rxtr_type, lag, fc_pd)
    # Valid type declarations
    valid_conv_types = ["greenfield", "electrical", "steam_noTES", "steam_TES"]
    valid_rxtr_types = ["PWR", "SFR", "HTGR"]

    # Check project categorical descriptors
    if !in(conv_type, valid_conv_types)
        err_msg = string(
            "Conversion type $conv_type is not valid, or is not currently supported.\n",
            "Currently supported conversion types:\n",
            valid_conv_types
        )
        throw(ArgumentError(conv_type, err_msg))
    end

     if !in(rxtr_type, valid_rxtr_types)
        err_msg = string(
            "Reactor type $rxtr_type is not valid, or is not currently supported.\n",
            "Currently supported reactor types:\n",
            valid_rxtr_types
        )
        throw(ArgumentError(rxtr_type, err_msg))
    end

    # Validate lag
    if lag < 0
        throw(DomainError(lag, "Lag must be a nonnegative integer."))
    end

    # Validate fc_pd
    if fc_pd < 1
        throw(DomainError(fc_pd, "fc_pd (forecast period) must be an integer greater than 0."))
    end

end


function project_capex_profile(base_pd, lag, db, fc_pd, rxtr_type, conv_type, C2N_specs; status="new", asset_id=nothing)
    # Determine the project's current status
    if status == "new"
        project_current_status = deepcopy(C2N_specs)
    elseif status == "ongoing"
        # Retrieve current project status from database
        pd_of_interest = base_pd - 1
        project_WIP_update = DBInterface.execute(db, "SELECT * FROM WIP_C2N WHERE asset_id = $asset_id AND period = $pd_of_interest") |> DataFrame
        project_WIP_update = project_WIP_update[1, :]
        project_current_status = Dict(
            :npp_ns_xtr => Dict("cost_rem" => project_WIP_update[:npp_ns_xtr_cost_rem], "time_rem" => project_WIP_update[:npp_ns_xtr_time_rem]),
            :npp_safety_xtr => Dict("cost_rem" => project_WIP_update[:npp_safety_xtr_cost_rem], "time_rem" => project_WIP_update[:npp_safety_xtr_time_rem]),
            :cpp_dnd => Dict("cost_rem" => project_WIP_update[:cpp_dnd_cost_rem], "time_rem" => project_WIP_update[:cpp_dnd_time_rem]),
            :cpp_wr => Dict("cost_rem" => project_WIP_update[:cpp_wr_cost_rem], "time_rem" => project_WIP_update[:cpp_wr_time_rem]),
            :cpp_nrc => Dict("cost_rem" => project_WIP_update[:cpp_nrc_cost_rem], "time_rem" => project_WIP_update[:cpp_nrc_time_rem])
        )
    end

    # Set up the timeline of expenditures
    max_horizon = fc_pd
    capex_tl = DataFrame(
        year = base_pd:base_pd+fc_pd-1,
        cpp_dnd = zeros(fc_pd),
        cpp_nrc = zeros(fc_pd),
        cpp_wr = zeros(fc_pd),
        npp_ns_xtr = zeros(fc_pd),
        npp_safety_xtr = zeros(fc_pd)
    )

    # Generate the table of activities ongoing during project intervals
    activity_schedule = get_C2N_available_activities(conv_type, project_current_status, base_pd, lag)

    # Convert the binary project schedule into a capex schedule
    capex_activity_schedule = allocate_funds_to_activities(activity_schedule, project_current_status, C2N_specs)

    # Convert the interval-based schedule into a yearly schedule
    capex_tl, capex_activity_schedule = convert_to_annual_capex_schedule(capex_activity_schedule)

    println(capex_activity_schedule)
    println(capex_tl)

    # Finalize the CapEx projection
    #transform!(capex_tl, [:npp_ns_xtr, :npp_safety_xtr, :cpp_wr, :cpp_nrc, :cpp_dnd] => ((a, b, c, d, e) -> a + b + c + d + e) => total_capex)

    return capex_tl

end


function get_C2N_available_activities(conv_type, project_current_status, base_pd, lag)
    all_activities = ["cpp_dnd", "cpp_nrc", "cpp_wr", "npp_ns_xtr", "npp_safety_xtr"]

    activity_schedule = DataFrame(
        state_start = Float64[],
        state_end = Float64[],
        cpp_dnd = Int[],
        cpp_nrc = Int[],
        cpp_wr = Int[],
        npp_ns_xtr = Int[],
        npp_safety_xtr = Int[]
    )

    project_ongoing = true

    # The first entry in the dataframe is always the lag period
    test_status = deepcopy(project_current_status)

    while project_ongoing
        available_activities = get_available_activities(test_status, conv_type)
        if size(available_activities)[1] == 0
            project_ongoing = false
            break
        end

        # Determine the expected end date of all activities which are currently
        #   available for work
        end_dates = DataFrame(activity = String[], time_rem = Float64[])
        for activity in available_activities
            time_rem = test_status[activity]["time_rem"]
            push!(end_dates, [activity, time_rem])
        end

        # Find the activity with the shortest time remaining
        sort!(end_dates, :time_rem)
        println(end_dates)
        next_finished = end_dates[1, :activity]
        next_interval_duration = end_dates[1, :time_rem]

        # Any activities eligible for work during this interval receive a 1
        #   in the activity_schedule; otherwise, they get a 0
        activity_bools = []
        for activity in all_activities
            if activity in available_activities
                push!(activity_bools, 1)
            else
                push!(activity_bools, 0)
            end
        end

        # Add the row to the activity_schedule dataframe
        if size(activity_schedule)[1] == 0
            next_start = base_pd + lag
        else
            next_start = last(activity_schedule[!, :state_end])
        end
        next_end = next_start + next_interval_duration
        push!(activity_schedule, hcat(next_start, next_end, transpose(activity_bools)))

        # Subtract this interval's duration from the time remaining for any
        #   available activities
        for activity in available_activities
            test_status[activity]["time_rem"] -= end_dates[1, :time_rem]
        end

    end


    return activity_schedule
end


function get_available_activities(project_current_status, conv_type)
    if conv_type == "greenfield"
        available_activities = get_available_activities_greenfield(project_current_status)
    elseif conv_type == "electrical"
        available_activities = get_available_activities_electrical(project_current_status)
    elseif conv_type == "steam_noTES"
        available_activities = get_available_activities_steam_noTES(project_current_status)
    elseif conv_type == "steam_TES"
        available_activities = get_available_activities_steam_TES(project_current_status)
    else
        println("I don't recognize that conversion project type. Check your inputs and try again.")
        exit()
    end

    return available_activities

end


function get_greenfield_available_activities(project_current_status)
    all_activities = ["npp_ns_xtr", "npp_safety_xtr", "cpp_dnd", "cpp_wr", "cpp_nrc"]
    for activity in all_activities
        act_time_rem = project_current_status[activity]["time_rem"]
        act_cost_rem = project_current_status[activity]["cost_rem"]
        if (round(act_time_rem, digits = 3) == 0) && (round(act_cost_rem, digits = 3) == 0)
            deleteat!(all_activities, findall(x -> x == activity, all_activities))
        end
    end

    # The first remaining activity in all_activities is the sole available
    #   activity
    # If no activities remain, the project is done; return empty list
    if size(all_activities)[1] != 0
        available_activities = [all_activities[1]]
    else
        available_activities = []
    end

    return available_activities
end


function get_available_activities_electrical(project_current_status)
    all_activities = [key for key in keys(project_current_status)]
    completed_activities = []
    for activity in all_activities
        act_time_rem = project_current_status[activity]["time_rem"]
        if (round(act_time_rem, digits = 3) == 0)
            push!(completed_activities, activity)
        end
    end

    incomplete_activities = [all_activities[i] for i = 1:size(all_activities)[1] if !(all_activities[i] in completed_activities)]

    # Add activities which are available to be worked on to the
    #   available_activities list
    # This section will not make much sense without the documentary flowcharts;
    #   see docs for explanation
    # TODO: replace with digraph traversal
    available_activities = []

    # Activities with no prerequisite: cpp_wr, npp_ns_xtr
    if "cpp_wr" in incomplete_activities
        push!(available_activities, "cpp_wr")
    end

    if "npp_ns_xtr" in incomplete_activities
        push!(available_activities, "npp_ns_xtr")
    end

    # CPP NRC license approval prerequisite: cpp_wr
    if ("cpp_nrc" in incomplete_activities) && !("cpp_wr" in incomplete_activities)
        push!(available_activities, "cpp_nrc")
    end

    # NPP safety xtr prerequisite: CPP NRC license
    if ("npp_safety_xtr" in incomplete_activities) && !("cpp_nrc" in incomplete_activities) && ("npp_ns_xtr" in completed_activities)
        push!(available_activities, "npp_safety_xtr")
    end

    # CPP D&D only after all other activities are done
    if !("npp_safety_xtr" in incomplete_activities) && !("cpp_dnd" in completed_activities)
        push!(available_activities, "cpp_dnd")
    end

    return available_activities
end

function get_available_activities_steam_noTES(project_current_status)
    all_activities = [key for key in keys(project_current_status)]
    completed_activities = []
    for activity in all_activities
        act_time_rem = project_current_status[activity]["time_rem"]
        if (round(act_time_rem, digits = 3) == 0)
            push!(completed_activities, activity)
        end
    end

    incomplete_activities = [all_activities[i] for i = 1:size(all_activities)[1] if !(all_activities[i] in completed_activities)]

    # Add activities which are available to be worked on to the
    #   available_activities list
    # This section will not make much sense without the documentary flowcharts;
    #   see docs for explanation
    # TODO: replace with digraph traversal
    available_activities = []

    # Activities with no prerequisites: cpp_wr, cpp_nrc, npp_ns_xtr
    if "cpp_wr" in incomplete_activities
        push!(available_activities, "cpp_wr")
    end

    if "cpp_nrc" in incomplete_activities
        push!(available_activities, "cpp_nrc")
    end

    if "npp_ns_xtr" in incomplete_activities
        push!(available_activities, "npp_ns_xtr")
    end

    # Only start NPP safety xtr after all three of the above are done
    if ("cpp_wr" in completed_activities) && ("cpp_nrc" in completed_activities) && ("npp_ns_xtr" in completed_activities) && ("npp_safety_xtr" in incomplete_activities)
        push!(available_activities, "npp_safety_xtr")
    end

    # Only start cpp_dnd after all other activities are complete
    if !("npp_safety_xtr" in incomplete_activities) && !("cpp_dnd" in completed_activities)
        push!(available_activities, "cpp_dnd")
    end

    println(available_activities)

    return available_activities
end


function get_available_activities_steam_TES(project_current_status)
    all_activities = [key for key in keys(project_current_status)]
    completed_activities = []
    for activity in all_activities
        act_time_rem = project_current_status[activity]["time_rem"]
        if (round(act_time_rem, digits = 3) == 0)
            push!(completed_activities, activity)
        end
    end

    incomplete_activities = [all_activities[i] for i = 1:size(all_activities)[1] if !(all_activities[i] in completed_activities)]

    # Add activities which are available to be worked on to the
    #   available_activities list
    # This section will not make much sense without the documentary flowcharts;
    #   see docs for explanation
    # TODO: replace with digraph traversal
    available_activities = []

    # Activities with no prerequisites: cpp_nrc, npp_ns_xtr
    if "cpp_nrc" in incomplete_activities
        push!(available_activities, "cpp_nrc")
    end

    if "npp_ns_xtr" in incomplete_activities
        push!(available_activities, "npp_ns_xtr")
    end

    # NPP safety xtr prerequisite: NPP NS xtr
    if ("npp_ns_xtr" in completed_activities) && ("npp_safety_xtr" in incomplete_activities)
        push!(available_activities, "npp_safety_xtr")
    end

    # After NPP safety xtr is done, then do cpp_dnd and cpp_wr
    if ("npp_safety_xtr" in completed_activities) && !("cpp_dnd" in completed_activities)
        push!(available_activities, "cpp_dnd")
    end

    if ("npp_safety_xtr" in completed_activities) && !("cpp_wr" in completed_activities)
        push!(available_activities, "cpp_wr")
    end

    return available_activities
end



function allocate_funds_to_activities(activity_schedule, C2N_specs, project_current_status)
    all_activities = [key for key in keys(C2N_specs)]

    capex_activity_schedule = deepcopy(activity_schedule)

    for activity in all_activities
        if C2N_specs[activity]["time_rem"] != 0
            linear_cost = C2N_specs[activity]["cost_rem"] / C2N_specs[activity]["time_rem"]
        else
            linear_cost = 0
        end
        transform!(capex_activity_schedule, [Symbol(activity), :state_start, :state_end] => ((activity, s_start, s_end) -> activity .* linear_cost * 900 * 1000) => Symbol(activity))
    end

    return capex_activity_schedule

end


function convert_to_annual_capex_schedule(capex_activity_schedule)
    start_pd = floor(Int64, minimum(capex_activity_schedule[!, :state_start]))
    end_pd = ceil(Int64, maximum(capex_activity_schedule[!, :state_end]))

    # Get unscaled capex totals
    transform!(capex_activity_schedule, [:cpp_dnd, :cpp_nrc, :cpp_wr, :npp_ns_xtr, :npp_safety_xtr] => ((a, b, c, d, e) -> a .+ b .+ c .+ d .+ e) => :unscaled_capex)

    capex_tl = DataFrame(
        period = start_pd:end_pd,
        total_capex = zeros(end_pd - start_pd + 1)
    )

    for i = start_pd:end_pd
        subset = filter([:state_start, :state_end] => ((s_start, s_end) -> !(s_end <= i) && !(s_start >= i + 1)), capex_activity_schedule)
        for j = 1:size(subset)[1]
            duration_in_i = min(i+1, subset[j, :state_end]) - max(i, subset[j, :state_start])
            capex_tl[i - start_pd + 1, :total_capex] += duration_in_i * subset[j, :unscaled_capex]
        end
    end

    return capex_tl, capex_activity_schedule

end



end
