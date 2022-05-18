module C2N

using Logging, CSV, DataFrames

mutable struct C2N_project
    conversion_type::String
    rxtr_type::String
    initial_state::String
    lag::Int64
    uid::Int64
    fs::DataFrame
end


function C2N_project(conversion_type, rxtr_type, initial_state, lag, uid, fc_pd)
    # Validate inputs
    check_init_inputs(conversion_type, initial_state, lag, uid, fc_pd)


    fs = DataFrame(
        year = 1:fc_pd,
        capex = zeros(fc_pd),
        remaining_debt_principal = zeros(fc_pd),
        debt_payment = zeros(fc_pd),
        interest_payment = zeros(fc_pd),
        depreciation = zeros(fc_pd),
        gen = zeros(fc_pd)
    )
    return C2N_project(C2N_type, lag, uid, fs)
end


function check_init_inputs(conversion_type, rxtr_type, initial_state, lag, uid, fc_pd)
    # Valid type declarations
    #valid_conversion_types = ["none", "electric", "steam_noTES", "steam_TES"]
    valid_conversion_types = ["none", "electric"]
    #valid_rxtr_types = ["pwr", "sfr", "htgr"]
    valid_rxtr_types = ["pwr"]
    valid_initial_states = ["greenfield", "C2N"]

    # Check project categorical descriptors
    if !in(conversion_type, valid_conversion_types):
        err_msg = string(
            "Conversion type $conversion_type is not valid, or is not currently supported.\n",
            "Currently supported conversion types:\n",
            valid_conversion_types
        )
        throw(ArgumentError(conversion_type, err_msg))
    end

     if !in(rxtr, valid_rxtr_types):
        err_msg = string(
            "Reactor type $rxtr_type is not valid, or is not currently supported.\n",
            "Currently supported reactor types:\n",
            valid_rxtr_types
        )
        throw(ArgumentError(rxtr_type, err_msg))
    end

    if !in(initial_state, valid_initial_states):
        err_msg = string(
            "Initial state $initial_state is not valid, or is not currently supported.\n",
            "Currently-supported initial states:\n",
            valid_initial_states
        )
        throw(ArgumentError(initial_state, err_msg))
    end

    # Validate lag
    if lag < 0
        throw(DomainError(lag, "Lag must be a nonnegative integer."))
    end

    # Validate uid
    if uid < 0
        throw(ArgumentError(uid, "uid must be a nonnegative integer."))
    end

    # Validate fc_pd
    if fc_pd < 1
        throw(DomainError(fc_pd, "fc_pd (forecast period) must be an integer greater than 0."))
    end

end


function project_capex_profile(C2N_PA, unit_specs; status="new", db=nothing, asset_id=nothing, current_pd=9999)
    # list of all project activities
    all_activities = ["cpp_dnd", "cpp_nrc", "npp_ns_xtr", "npp_safety_xtr", "npp_comm"]

    # Short parameter names for convenience
    rxtr_type = C2N_PA[:rxtr_type]
    conv_type = C2N_PA[:conversion_type]

    # Retrieve unit specs for convenience
    rxtr_specs = filter(:unit_type => x -> x == rxtr_type, unit_specs)

    # Retrieve specs for this project type for convenience
    C2N_specs = DBInterface.execute(db, "SELECT * FROM C2N_specs WHERE reactor_type = $rxtr_type AND conversion_type = $conv_type") |> DataFrame
    C2N_specs = C2N_specs[1, :]

    # Determine the project's current status
    if status == "new"
        project_current_status = deepcopy(C2N_specs)[1, :]
    elseif status == "ongoing"
        # Retrieve current project status from database
        project_current_status = DBInterface.execute(db, "SELECT * FROM WIP_C2N WHERE asset_id = $asset_id AND period = $current_pd") |> DataFrame
        project_current_status = project_current_status[1, :]

    end

    # Set up the timeline of expenditures
    max_horizon = 15
    capex_tl = DataFrame(
        year = 1:max_horizon,
        cpp_dnd = zeros(max_horizon),
        cpp_nrc = zeros(max_horizon),
        npp_ns_xtr = zeros(max_horizon),
        npp_safety_xtr = zeros(max_horizon),
        npp_comm = zeros(max_horizon)
    )

    activity_schedule = get_C2N_available_activities(conversion_type, initial_state, project_current_status)

    capex_tl = allocate_funds_to_activities(activity_schedule, capex_tl, project_current_status, C2N_specs)

    # Finalize the CapEx projection
    transform!(capex_tl, [:cpp_dnd, :cpp_nrc, :npp_ns_xtr, :npp_safety_xtr, :npp_comm] => ((a, b, c, d, e) -> a + b + c + d + e) => :total_capex)

    return capex_tl

end


function get_C2N_available_activities(conversion_type, initial_state, project_current_status)
    all_activities = ["cpp_dnd", "cpp_nrc", "npp_ns_xtr", "npp_safety_xtr", "npp_comm"]

    activity_schedule = DataFrame(
        state_end = Float64[],
        cpp_dnd = Int[],
        cpp_nrc = Int[],
        npp_ns_xtr = Int[],
        npp_safety_xtr = Int[],
        npp_comm = Int[]
    )

    project_ongoing = true
    while project_ongoing
        # Determine which activities are available immediately
        if initial_state == "greenfield"
            available_activities = get_greenfield_available_activities(project_current_status)
        elseif initial_state == "C2N"
            if conversion_type == "electrical"
                available_activities = get_elec_C2N_available_activities(project_current_status)
            end
        end

        # Determine the expected end date of all activities which are currently
        #   available for work
        end_dates = DataFrame(:activity = String[], :time_rem = Float64)
        for activity in available_activities:
            act_time_rem = Symbol(string(activity, "_time_rem"))
            time_rem = project_current_status[act_time_rem]
            push!(end_dates, (:activity = activity, :time_rem = time_rem))
        end

        # Find the activity with the shortest time remaining
        sort!(end_dates, :time_rem)
        next_finished = end_dates[1, :activity]
        if size(activity_schedule)[1] == 0
            next_time = end_dates[1, :time_rem]
        else
            next_time = last(activity_schedule)[1, :state_end] + end_dates[1, :time_rem]
        end


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
        push!(activity_schedule, transpose(vcat(next_time, activity_bools)))

        # If no activities remain, end the loop
        if size(available_activities)[1] == 0
            project_ongoing = false
        end
    end


    return activity_schedule
end

function get_greenfield_available_activities(project_current_status)
    activity_order = ["npp_ns_xtr", "npp_safety_xtr", "npp_comm", "cpp_dnd", "cpp_nrc"]
    for activity in activity_order
        act_time_rem = Symbol(string(activity, "_time_rem"))
        act_cost_rem = Symbol(string(activity, "_cost_rem"))
        if (round(Int, project_current_status[act_time_rem]) == 0) && (round(Int, project_current_status[act_cost_rem]) == 0)
            deleteat!(activity_order, findall(x -> x == activity, activity_order))
        end
    end

    # The first remaining activities in activity_order is the sole available
    #   activity
    # If no activities remain, the project is done; return empty list
    if size(activity_order)[1] != 0
        available_activities = [activity_order[1]]
    else
        available_activities = []
    end

    return available_activities
end


function get_elec_C2N_available_activities(project_current_status)
    all_activities = ["npp_ns_xtr", "npp_safety_xtr", "npp_comm", "cpp_dnd", "cpp_nrc"]
    complete_activities = []
    for activity in all_activities
        act_time_rem = Symbol(string(activity, "_time_rem"))
        act_cost_rem = Symbol(string(activity, "_cost_rem"))
        if (round(Int, project_current_status[act_time_rem]) == 0) && (round(Int, project_current_status[act_cost_rem]) == 0)
            push!(complete_activities, activity)
        end
    end

    incomplete_activities = [all_activities[i] for i = 1:size(all_activities)[1] if !(all_activities[i] in completed_activities)]

    # Add activities which are available to be worked on to the
    #   available_activities list
    # This section will not make much sense without the documentary flowcharts;
    #   see docs for explanation
    # TODO: replace with digraph traversal
    available_activities = []

    if "cpp_dnd" in incomplete_activities
        push!(available_activities, "cpp_dnd")
    end

    if !("cpp_dnd" in incomplete_activities) && ("npp_ns_xtr" in incomplete_activities)
        push!(available_activities, "npp_ns_xtr")
    end

    if "cpp_nrc" in incomplete_activities
        push!(available_activities, "cpp_nrc")
    end

    if !("cpp_dnd" in incomplete_activities) && !("cpp_nrc" in incomplete_activities) && !("npp_ns_xtr" in incomplete_activities)
        push!(available_activities, "npp_safety_xtr")
    end

    if "npp_safety_xtr" in complete_activities
        push!(available_activities, "npp_comm")
    end

    return available_activities
end


function allocate_funds_to_activities(activity_schedule, capex_tl, project_status, C2N_specs)
    all_activities = ["npp_ns_xtr", "npp_safety_xtr", "npp_comm", "cpp_dnd", "cpp_nrc"]

    capex_tl[!, :duration] = zeros(size(capex_tl)[1])
    capex_tl[1, :duration] = activity_schedule[1, :state_end]
    for i = 2:size(capex_tl)[1]
        capex_tl[i, :duration] = activity_schedule[i, :state_end] - activity_schedule[i-1, :state_end]
    end

    for activity in all_activities
        # Add the activity schedule indicator column into capex_tl
        capex_tl[!, Symbol(string(activity, "_i"))] = activity_schedule[!, activity]

        activity_duration = sum(capex_tl[!, :duration] .* capex_tl[!, Symbol(string(activity, "_i"))])

        # Assume that costs per activity take place linearly over the course of
        #   that activity: prorate total remaining activity costs over total
        #   remaining activity time
        for i = 1:size(capex_tl)[1]
            capex_tl[i, activity] = capex_tl[i, :duration] / activity_duration * project_status[Symbol(string(activity, "_cost_rem"))]
        end
    end   

    return capex_tl

end





end
