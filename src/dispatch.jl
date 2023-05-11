module Dispatch

using Requires, Logging, CSV, DataFrames, JuMP, GLPK, Cbc, XLSX, SQLite, HiGHS


# Initialize this module, with CPLEX as an optional library if available
function __init__()
    @require CPLEX = "a076750e-1247-5638-91d2-ce28b192dca0" @eval using CPLEX
end


function execute_dispatch_economic_projection(
    CLI_args,
    db,
    settings,
    fc_pd,
    total_demand,
    unit_specs,
    system_portfolios,
)
    @debug string(
        "Running the dispatch simulation for ",
        settings["dispatch"]["num_dispatch_years"],
        " years...",
    )

    # Set up all timeseries data
    ts_data = load_ts_data(
        joinpath(
            CLI_args["inputs_path"],
            "ts_data",
        ),
        settings["dispatch"]["num_repdays"],
    )

    system_portfolios = fill_portfolios_missing_units(system_portfolios, unit_specs)

    all_gc_results, all_price_results = handle_annual_dispatch(
        settings,
        CLI_args["current_pd"],
        system_portfolios,
        total_demand,
        ts_data,
        unit_specs,
    )

    long_econ_results = postprocess_results(
        settings,
        all_gc_results,
        all_price_results,
        ts_data[:repdays_data],
        system_portfolios,
        unit_specs,
        CLI_args["current_pd"],
        fc_pd,
    )

    return long_econ_results
end


function get_system_portfolios(db, settings, start_year, unit_specs)
    # Set up portfolio dictionaries
    @debug "Setting up dispatch portfolios..."
    system_portfolios = Dict()

    end_year = (start_year
                + convert(Int64, settings["dispatch"]["num_dispatch_years"])
                - 1
               )

    for y = start_year:end_year
        # Retrieve annual system portfolios
        system_portfolio =
            DBInterface.execute(
                db,
                string(
                    "SELECT unit_type, COUNT(unit_type) ",
                    "FROM assets WHERE completion_pd <= $y ",
                    "AND retirement_pd > $y ",
                    "AND cancellation_pd > $y ",
                    "GROUP BY unit_type",
                ),
            ) |> DataFrame

        # Clean up column names in the portfolio dataframes
        rename!(system_portfolio, Symbol("COUNT(unit_type)") => :num_units)

        system_portfolios[y] = system_portfolio
    end

    return system_portfolios
end


function fill_portfolios_missing_units(system_portfolios, unit_specs)
    # Ensure that at least 1 unit of every type in unit_specs is represented in
    #   every year of the system portfolio, by adding 1 instance of missing
    #   unit types.
    for y=minimum(keys(system_portfolios)):maximum(keys(system_portfolios))
        for unit_type in unit_specs[!, :unit_type]
            if !in(unit_type, system_portfolios[y][!, :unit_type])
                push!(system_portfolios[y], (unit_type, 1))
            end
        end
    end

    return system_portfolios
end


function handle_annual_dispatch(
    settings,
    current_pd,
    system_portfolios,
    total_demand,
    ts_data,
    unit_specs,
)
    all_gc_results = set_up_gc_results_df()
    all_prices = set_up_prices_df()

    # Run the annual dispatch for the user-specified number of dispatch years
    for y = current_pd:current_pd + settings["dispatch"]["num_dispatch_years"] - 1
        @debug "\n\nDISPATCH SIMULATION: YEAR $y"

        # Select the current year's expected portfolio
        year_portfolio = system_portfolios[y]

        # Determine appropriate total demand for this year
        year_demand =
            filter(:period => ((pd) -> pd == y), total_demand)[1, :real_demand]

        # Set up and run the dispatch simulation for this year
        # This function updates all_gc_results and all_prices in-place, and
        #   returns a boolean to determine whether the next year should be run
        results = run_annual_dispatch(
            settings,
            y,
            year_portfolio,
            year_demand,
            ts_data,
            unit_specs,
        )

        # Save new generation, commitment, and price results
        all_gc_results = vcat(all_gc_results, results[:new_gc_results])
        all_prices = vcat(all_prices, results[:new_prices])

        @debug "DISPATCH SIMULATION: YEAR $y COMPLETE."

        if !results[:run_next_year]
            break
        end

        if (
            (results[:total_ENS] > settings["system"]["max_total_ENS"]) &&
            (y == current_pd)
        )
            msg = string(
                "In the upcoming dispatch year, the system is ",
                "projected to experience a total energy shortage ",
                "which exceeds the maximum allowed level of energy ",
                "not served (ENS). The system portfolio is likely ",
                "in an unsalvageable state of supply insufficiency ",
                "in the immediate term. The run will now terminate.",
            )
            throw(ErrorException(msg))
        end

    end

    # If no years ran correctly, throw an error and exit
    if size(all_gc_results)[1] == 0
        msg = string(
            "No dispatch simulations could be run. Re-check ",
            "your inputs.",
        )
        throw(ErrorException(msg))
    end

    return all_gc_results, all_prices

end


function load_ts_data(ts_file_dir, num_repdays)
    # Load the time-series demand and VRE data into dataframes
    load_data =
        CSV.read(joinpath(ts_file_dir, "timeseries_load_hourly.csv"), DataFrame)
    wind_data =
        CSV.read(joinpath(ts_file_dir, "timeseries_wind_hourly.csv"), DataFrame)
    solar_data =
        CSV.read(joinpath(ts_file_dir, "timeseries_pv_hourly.csv"), DataFrame)
    repdays_data =
        CSV.read(joinpath(ts_file_dir, "repDays_$num_repdays.csv"), DataFrame)

    ts_data = Dict(
        :load_data => load_data,
        :wind_data => wind_data,
        :solar_data => solar_data,
        :repdays_data => repdays_data,
    )

    return ts_data
end

function set_up_gc_results_df()
    all_gc_results = DataFrame(
        y = Int[],
        d = Int[],
        h = Int[],
        unit_type = String[],
        gen = Float64[],
        commit = Int[],
    )

    return all_gc_results
end


function set_up_prices_df()
    all_prices = DataFrame(y = Int[], d = Int[], h = Int[], price = Float64[])

    return all_prices
end


function scale_load(ts_data, peak_demand)
    ts_data[:load_data][!, :Load] =
        (ts_data[:load_data][!, :LoadShape] * peak_demand)

    return ts_data
end


function set_up_load_repdays(ts_data)
    load_repdays = DataFrame()
    for day in ts_data[:repdays_data][!, :Day]
        load_repdays[!, Symbol(day)] =
            ts_data[:load_data][(24 * day + 1):(24 * (day + 1)), :Load]
    end

    ts_data[:load_repdays] = load_repdays

    return ts_data
end


function scale_wind_solar_data(ts_data, year_portfolio, unit_specs)
    wind_specs = filter(:unit_type => x -> x == "wind", unit_specs)[1, :]
    solar_specs = filter(:unit_type => x -> x == "solar", unit_specs)[1, :]

    # For non-zero wind and solar capacity, scale the WindShape series by the
    #   installed capacity to get total instantaneous VRE availability
    # If either wind or solar has 0 installed capacity, set its entire time
    #   series to 0.0
    if !isempty(filter(:unit_type => x -> x == "wind", year_portfolio))
        ts_data[:wind_data][!, :wind] = (
            ts_data[:wind_data][!, :WindShape] *
            filter(:unit_type => x -> x == "wind", year_portfolio)[
                1,
                :num_units,
            ] *
            wind_specs[:capacity]
        )
    else
        ts_data[:wind_data][!, :wind] .= 0.0
    end

    if !isempty(filter(:unit_type => x -> x == "solar", year_portfolio))
        ts_data[:solar_data][!, :solar] = (
            ts_data[:solar_data][!, :SolarShape] *
            filter(:unit_type => x -> x == "solar", year_portfolio)[
                1,
                :num_units,
            ] *
            solar_specs[:capacity]
        )
    else
        ts_data[:solar_data][!, :solar] .= 0.0
    end

    return ts_data
end


function set_up_wind_solar_repdays(ts_data)
    wind_repdays = DataFrame()
    solar_repdays = DataFrame()

    for day in ts_data[:repdays_data][!, :Day]
        wind_repdays[!, Symbol(day)] =
            ts_data[:wind_data][(24 * day + 1):(24 * (day + 1)), :wind]
        solar_repdays[!, Symbol(day)] =
            ts_data[:solar_data][(24 * day + 1):(24 * (day + 1)), :solar]
    end

    ts_data[:wind_repdays] = wind_repdays
    ts_data[:solar_repdays] = solar_repdays

    return ts_data
end


function set_up_model(settings, ts_data, year_portfolio, unit_specs)
    # Create joined portfolio-unit_specs dataframe, to ensure consistent
    #   accounting for units which are actually present and consistent
    #   unit ordering
    portfolio_specs = innerjoin(unit_specs, year_portfolio, on = :unit_type)
    portfolio_specs[!, :unit_index] = 1:size(portfolio_specs)[1]

    # Set up wind_index and solar_index for easy filtering later
    wind_index = 0
    solar_index = 0
    if !isempty(filter(:unit_type => x -> x == "wind", portfolio_specs))
        wind_index = filter(:unit_type => x -> x == "wind", portfolio_specs)[
            1,
            :unit_index,
        ]
    end
    if !isempty(filter(:unit_type => x -> x == "solar", portfolio_specs))
        solar_index = filter(:unit_type => x -> x == "solar", portfolio_specs)[
            1,
            :unit_index,
        ]
    end

    # Helpful named constants
    num_units = size(portfolio_specs)[1]
    num_days = size(ts_data[:repdays_data])[1]
    num_hours = 24

    # Break out timeseries data sets for convenience
    load_data = ts_data[:load_data]
    wind_data = ts_data[:wind_data]
    solar_data = ts_data[:solar_data]
    repdays_data = ts_data[:repdays_data]
    load_repdays = ts_data[:load_repdays]
    wind_repdays = ts_data[:wind_repdays]
    solar_repdays = ts_data[:solar_repdays]

    # Initialize JuMP model
    if lowercase(settings["simulation"]["solver"]) == "cplex"
        m = Model(CPLEX.Optimizer)
    elseif lowercase(settings["simulation"]["solver"]) == "glpk"
        m = Model(GLPK.Optimizer)
    elseif lowercase(settings["simulation"]["solver"]) == "cbc"
        m = Model(Cbc.Optimizer)
    elseif lowercase(settings["simulation"]["solver"]) == "highs"
        m = Model(HiGHS.Optimizer)
    else
        throw(error("Solver not supported. Try `cplex` instead."))
    end
    set_silent(m)

    # g: quantity generated (in MWh) for each unit type
    @variable(m, g[1:num_units, 1:num_days, 1:num_hours] >= 0)

    # c: number of units of each type committed in each hour
    @variable(m, c[1:num_units, 1:num_days, 1:num_hours] >= 0, Int)

    # s: penalty variable for energy not served in each hour
    @variable(m, s[1:num_days, 1:num_hours] >= 0)

    # Total generation per hour must be greater than or equal to the demand
    @constraint(
        m,
        mkt_equil[k = 1:num_days, j = 1:num_hours],
        sum(g[i, k, j] for i = 1:num_units) + s[k, j] >= load_repdays[j, k]
    )

    # Number of committed units per hour must be less than or equal to the
    #   total number of units of that type in the system
    for i = 1:num_units
        for k = 1:num_days
            for j = 1:num_hours
                @constraint(m, c[i, k, j] <= portfolio_specs[i, :num_units])
            end
        end
    end

    # Limit total generation per unit type each hour to the total capacity of
    #   all committed units of this type, with committed units subject to
    #   minimum and maximum power levels
    # wind and solar are excluded
    for i = 1:num_units
        if !(portfolio_specs[i, :unit_type] in ["wind", "solar"])
            for k = 1:num_days
                for j = 1:num_hours
                    @constraint(
                        m,
                        g[i, k, j] <= (
                            c[i, k, j] .* portfolio_specs[i, :capacity] .*
                            portfolio_specs[i, :capacity_factor] .*
                            portfolio_specs[i, :max_PL]
                        )
                    )
                    @constraint(
                        m,
                        g[i, k, j] >= (
                            c[i, k, j] .* portfolio_specs[i, :capacity] .*
                            portfolio_specs[i, :capacity_factor] .*
                            portfolio_specs[i, :min_PL]
                        )
                    )
                end
            end
        end
    end

    # Limit solar and wind generation to their actual hourly availability, if
    #   wind and/or solar exist in the system
    for k = 1:num_days
        for j = 1:num_hours
            if wind_index != 0
                @constraint(m, g[wind_index, k, j] <= wind_repdays[j, k])
            end
            if solar_index != 0
                @constraint(m, g[solar_index, k, j] <= solar_repdays[j, k])
            end
        end
    end

    # Ramping constraints
    for i = 1:num_units
        for k = 1:num_days
            for j = 1:(num_hours - 1)
                # Ramp-up constraint
                @constraint(
                    m,
                    (
                        g[i, k, j + 1] - g[i, k, j] <=
                        c[i, k, j + 1] .* portfolio_specs[i, :ramp_up_limit] *
                        portfolio_specs[i, :capacity] *
                        portfolio_specs[i, :capacity_factor]
                    )
                )
                # Ramp-down constraint
                @constraint(
                    m,
                    (
                        g[i, k, j + 1] - g[i, k, j] >=
                        (-1) *
                        c[i, k, j + 1] *
                        portfolio_specs[i, :ramp_down_limit] *
                        portfolio_specs[i, :capacity] *
                        portfolio_specs[i, :capacity_factor]
                    )
                )
            end
        end
    end

    ENS_penalty = settings["constants"]["big_number"]

    @objective(
        m,
        Min,
        sum(
            sum(
                sum(g[i, k, j] + ENS_penalty * s[k, j] for j = 1:num_hours) for
                k = 1:num_days
            ) .* (
                portfolio_specs[i, :VOM] + portfolio_specs[i, :FC_per_MWh] -
                portfolio_specs[i, :policy_adj_per_MWh]
            ) for i = 1:num_units
        )
    )

    return m, portfolio_specs
end


function solve_model(model; model_type = "integral")
    if model_type == "integral"
        # Optimize!
        optimize!(model)
        status = string(termination_status.(model))
        @debug "$model_type model status: $status"
        gen_qty = nothing
        c = nothing
        s = nothing
        if status == "OPTIMAL"
            gen_qty = value.(model[:g])
            c = value.(model[:c])
            s = value.(model[:s])
        end

        returns = model, status, gen_qty, c, s

    elseif model_type == "relaxed_integrality"
        optimize!(model)
        status = string(termination_status.(model))
        @debug "$model_type model status: $status"

        returns = model

    end

    return returns
end


function assemble_gc_results(y, gen_qty, c, portfolio_specs)
    new_gc_results = set_up_gc_results_df()

    # Save generation and commitment results to a dataframe
    # g[i, k, j] and c[i, k, j] are arranged as
    #   [num_units, num_days, num_hours]
    for k = 1:size(c)[2]            # num_days
        for j = 1:size(c)[3]        # num_hours
            for i = 1:size(c)[1]    # num_units
                line = (
                    y = y,
                    d = k,
                    h = j,
                    unit_type = portfolio_specs[i, :unit_type],
                    gen = gen_qty[i, k, j],
                    commit = round(Int, c[i, k, j]),
                )
                push!(new_gc_results, line)
            end
        end
    end

    return new_gc_results
end


function reshape_shadow_prices(shadow_prices, y, settings)
    price_df = DataFrame(y = Int[], d = Int[], h = Int[], price = Float64[])

    # Convert the (repdays, hours) table into a long (y, d, h, price) table
    for k = 1:size(shadow_prices)[1]        # num_days
        for j = 1:size(shadow_prices)[2]    # num_hours
            # Ensure the shadow price is no greater than the system cap
            price = (-1) * shadow_prices[k, j]
            if price > settings["system"]["price_cap"]
                price = settings["system"]["price_cap"]
            end
            line = (y = y, d = k, h = j, price = price)
            push!(price_df, line)
        end
    end

    return price_df
end


function propagate_all_results(all_gc_results, all_prices, current_pd, end_year)
    final_dispatched_year = maximum(all_gc_results[!, :y])

    final_year_gc =
        filter(:y => x -> x == final_dispatched_year, all_gc_results)
    final_year_prices =
        filter(:y => x -> x == final_dispatched_year, all_prices)

    for y = (final_dispatched_year + 1):(current_pd + end_year - 1)
        # Copy the final_year_gc results forward, updating the year
        next_year_gc = deepcopy(final_year_gc)
        next_year_gc[!, :y] .= y
        # Push row-by-row (avoids Julia's strict for-loop scoping)
        for i = 1:size(next_year_gc)[1]
            push!(all_gc_results, next_year_gc[i, :])
        end

        # Copy the final_year_prices results forward, updating the year
        next_year_prices = deepcopy(final_year_prices)
        next_year_prices[!, :y] .= y
        for i = 1:size(next_year_prices)[1]
            push!(all_prices, next_year_prices[i, :])
        end
    end

    return all_gc_results, all_prices

end


function run_annual_dispatch(
    settings,
    y,
    year_portfolio,
    peak_demand,
    ts_data,
    unit_specs,
)
    # Scale the load data to the PD value for this year
    ts_data = scale_load(ts_data, peak_demand)

    # Set up representative days for load
    ts_data = set_up_load_repdays(ts_data)

    # Scale the wind and solar data according to the current year's total
    #   installed capacity
    ts_data = scale_wind_solar_data(ts_data, year_portfolio, unit_specs)

    # Set up representative days for wind and solar
    ts_data = set_up_wind_solar_repdays(ts_data)

    @debug "Setting up optimization model..."
    m, portfolio_specs =
        set_up_model(settings, ts_data, year_portfolio, unit_specs)

    @debug "Optimization model set up."
    @debug string("Solving repday dispatch for year ", y, "...")

    # Create a copy of the model, to use later for the relaxed-integrality
    #   solution
    m_copy = copy(m)
    if lowercase(settings["simulation"]["solver"]) == "cplex"
        set_optimizer(m_copy, CPLEX.Optimizer)
    elseif lowercase(settings["simulation"]["solver"]) == "glpk"
        set_optimizer(m_copy, GLPK.Optimizer)
    elseif lowercase(settings["simulation"]["solver"]) == "cbc"
        set_optimizer(m_copy, Cbc.Optimizer)
    elseif lowercase(settings["simulation"]["solver"]) == "highs"
        set_optimizer(m_copy, HiGHS.Optimizer)
    end
    set_silent(m_copy)

    # Solve the integral optimization problem
    m, status, gen_qty, c, s = solve_model(m, model_type = "integral")

    if status == "OPTIMAL"
        # Save the generation and commitment results from the integral problem
        new_gc_results = assemble_gc_results(y, gen_qty, c, portfolio_specs)

        # Set up a relaxed-integrality version of this model, to allow
        #   retrieval of dual values for the mkt_equil constraint
        @debug "Solving relaxed-integrality problem to get shadow prices..."
        undo = relax_integrality(m_copy)

        # Solve the relaxed-integrality model to compute the shadow prices
        m_copy = solve_model(m_copy, model_type = "relaxed_integrality")
        new_prices = reshape_shadow_prices(
            shadow_price.(m_copy[:mkt_equil]),   # shadow prices
            y,
            settings,
        )

        # Determine the total level of energy not served (ENS) for this
        #   dispatch result
        total_ENS = sum(sum(s[k, j] for k = 1:size(s)[1]) for j = 1:size(s)[2])

        @debug "Year $y dispatch run complete."
        run_next_year = true

    else
        # Dispatcher is unable to solve the current year: assume agents have
        #   not had a chance to address conditions in this year yet
        # Break the loop
        run_next_year = false
        total_ENS = nothing
    end

    results = Dict(
        :new_gc_results => new_gc_results,
        :new_prices => new_prices,
        :run_next_year => run_next_year,
        :total_ENS => total_ENS,
    )

    return results

end


function save_raw_results(all_prices, all_gc_results)
    pfile = joinpath(pwd(), "tmp", "price_results.csv")
    CSV.write(pfile, all_prices)

    gcfile = joinpath(pwd(), "tmp", "./gc_results.csv")
    CSV.write(gcfile, all_gc_results)
end


function combine_and_extend_year_portfolios(system_portfolios, forecast_end_pd)
    # Combine the existing system portfolios into a single dataframe with
    #   indicator column :y
    all_year_portfolios = DataFrame()
    for key in keys(system_portfolios)
        df = system_portfolios[key]
        df[!, :y] .= key
        append!(all_year_portfolios, df)
    end

    # Extend dispatch results by assuming no change after last dispatch year
    last_dispatch_year = maximum([key for key in keys(system_portfolios)])
    for i = last_dispatch_year+1:forecast_end_pd
        df = system_portfolios[last_dispatch_year]
        df[!, :y] .= i
        append!(all_year_portfolios, df)
    end

    return all_year_portfolios

end


function join_results_data_frames(
    all_gc_results,
    all_prices,
    repdays_data,
    all_year_portfolios,
    unit_specs,
)
    # Join in price data to all_gc_results
    long_econ_results = innerjoin(all_gc_results, all_prices, on = [:y, :d, :h])

    # Incorporate repdays probability data
    long_econ_results = innerjoin(
        long_econ_results,
        select(repdays_data, [:index, :Probability]),
        on = [:d => :index],
    )

    # Join in limited unit_specs data to compute some cash flows
    long_econ_results = innerjoin(
        long_econ_results,
        select(
            unit_specs,
            [:unit_type, :capacity, :VOM, :FC_per_MWh, :policy_adj_per_MWh],
        ),
        on = :unit_type,
    )

    # Join in unit number data to long_rev_results
    long_econ_results =
        innerjoin(long_econ_results, all_year_portfolios, on = [:y, :unit_type])

    return long_econ_results

end


function compute_per_unit_cash_flows(long_econ_results)
    # Calculate generation
    transform!(long_econ_results, [:gen, :Probability, :num_units] => ((gen, prob, num_units) -> gen .* prob .* 365 ./ num_units) => :annualized_gen_per_unit)

    # Calculate revenues
    transform!(
        long_econ_results,
        [:gen, :price, :Probability, :num_units] =>
            (
                (gen, price, prob, num_units) ->
                    gen .* price .* prob .* 365 ./ num_units
            ) => :annualized_rev_per_unit,
    )

    # Calculate VOM
    transform!(
        long_econ_results,
        [:gen, :VOM, :Probability, :num_units] =>
            (
                (gen, VOM, prob, num_units) ->
                    gen .* VOM .* prob .* 365 ./ num_units
            ) => :annualized_VOM_per_unit,
    )

    # Calculate fuel cost
    transform!(
        long_econ_results,
        [:gen, :FC_per_MWh, :Probability, :num_units] =>
            (
                (gen, fc, prob, num_units) ->
                    gen .* fc .* prob .* 365 ./ num_units
            ) => :annualized_FC_per_unit,
    )

    # Calculate policy adjustment
    transform!(
        long_econ_results,
        [:gen, :policy_adj_per_MWh, :Probability, :num_units] =>
            (
                (gen, adj, prob, num_units) ->
                    gen .* adj .* prob .* 365 ./ num_units
            ) => :annualized_policy_adj_per_unit,
    )

    for col in eachcol(long_econ_results)
        replace!(col, Inf => 0)
        replace!(col, NaN => 0)
    end

    return long_econ_results

end


function summarize_dispatch_results(settings, unit_specs, long_econ_results)
    dispatch_results = deepcopy(long_econ_results)

    dispatch_results = combine(groupby(dispatch_results, [:y, :unit_type]), [:annualized_gen_per_unit, :annualized_rev_per_unit, :annualized_VOM_per_unit, :annualized_FC_per_unit, :annualized_policy_adj_per_unit] .=> sum)

    rename!(dispatch_results, :annualized_gen_per_unit_sum => :generation, :annualized_rev_per_unit_sum => :revenue, :annualized_VOM_per_unit_sum => :VOM, :annualized_FC_per_unit_sum => :fuel_cost, :annualized_policy_adj_per_unit_sum => :policy_adj)

    # Pivot in FOM data
    FOM_data = select(unit_specs, [:unit_type, :FOM, :capacity])
    dispatch_results = leftjoin(dispatch_results, FOM_data, on = :unit_type)
    transform!(dispatch_results, [:FOM, :capacity] => ((FOM, cap) -> FOM .* cap .* settings["constants"]["MW2kW"]) => :FOM)
    select!(dispatch_results, Not(:capacity))

    # Put the data into a long format for easier filtering
    dispatch_results = stack(dispatch_results, [:generation, :revenue, :VOM, :fuel_cost, :FOM, :policy_adj])
    rename!(dispatch_results, :variable => :dispatch_result, :value => :qty)

    return dispatch_results
end


function postprocess_results(
    settings,
    all_gc_results,
    all_prices,
    repdays_data,
    all_year_system_portfolios,
    unit_specs,
    current_pd,
    fc_pd,
)
    # Propagate the results dataframes out to the end of the projection horizon
    # Assume no change after the last modeled year
    all_gc_results, all_prices =
        propagate_all_results(all_gc_results, all_prices, current_pd, fc_pd)

    # Get a single long dataframe with unit numbers by type for each year
    all_year_portfolios = combine_and_extend_year_portfolios(
        all_year_system_portfolios,
        current_pd + fc_pd,
    )

    # Join in relevant data to create a single master dataframe suitable for
    #   a variety of column-wise operations
    long_econ_results = join_results_data_frames(
        all_gc_results,
        all_prices,
        repdays_data,
        all_year_portfolios,
        unit_specs,
    )

    # Compute cash flows on a per-unit basis: revenue, VOM, fuel cost, and
    #   policy adjustment
    long_econ_results = compute_per_unit_cash_flows(long_econ_results)

    dispatch_results = summarize_dispatch_results(settings, unit_specs, long_econ_results)

    return long_econ_results, dispatch_results

end


end
