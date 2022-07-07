module Dispatch

using Logging, CSV, DataFrames, JuMP, GLPK, XLSX, Logging, CPLEX, SQLite


function execute_dispatch_economic_projection(db, settings, current_pd, fc_pd, total_demand, unit_specs, all_year_system_portfolios)
    # @info string("Running the dispatch simulation for ", settings["num_dispatch_years"], " years...")

    # Set up all timeseries data
    ts_data = load_ts_data(
                  joinpath(settings["ABCE_abs_path"],
                           "inputs",
                           "ALEAF_inputs"
                  ),
                  settings["num_repdays"]
              )

    # Set up dataframes to record all results
    all_prices, all_gc_results = set_up_results_dfs()

    long_econ_results = handle_annual_dispatch(
                            settings,
                            current_pd,
                            fc_pd,
                            all_year_system_portfolios,
                            total_demand,
                            ts_data,
                            unit_specs,
                            all_gc_results,
                            all_prices
                        )

    return long_econ_results
end


function set_up_dispatch_portfolios(db, start_year, fc_pd, agent_id, unit_specs)
    # Set up portfolio dictionaries
    @debug "Setting up dispatch portfolios..."
    all_year_system_portfolios = Dict()
    all_year_agent_portfolios = Dict()

    for y = start_year:start_year+fc_pd
        # Retrieve annual system portfolios
        system_portfolio = DBInterface.execute(db, "SELECT unit_type, COUNT(unit_type) FROM assets WHERE completion_pd <= $y AND retirement_pd > $y AND cancellation_pd > $y GROUP BY unit_type") |> DataFrame

        # Retrieve annual system portfolios for the current agent
        agent_portfolio = DBInterface.execute(db, "SELECT unit_type, COUNT(unit_type) FROM assets WHERE agent_id = $agent_id AND completion_pd <= $y AND retirement_pd > $y AND cancellation_pd > $y GROUP BY unit_type") |> DataFrame

        # Clean up column names in the portfolio dataframes
        rename!(system_portfolio, Symbol("COUNT(unit_type)") => :num_units)
        rename!(agent_portfolio, Symbol("COUNT(unit_type)") => :num_units)

        all_year_system_portfolios[y] = system_portfolio
        all_year_agent_portfolios[y] = agent_portfolio
    end

    return all_year_system_portfolios, all_year_agent_portfolios

end


function handle_annual_dispatch(settings, current_pd, fc_pd, all_year_system_portfolios, total_demand, ts_data, unit_specs, all_gc_results, all_prices)
    # Run the annual dispatch for the user-specified number of dispatch years
    for y = current_pd:current_pd+settings["num_dispatch_years"]
        # @info "\n\nDISPATCH SIMULATION: YEAR $y"

        # Select the current year's expected portfolio
        year_portfolio = all_year_system_portfolios[y]
        # @info year_portfolio

        # Determine appropriate total demand for this year
        year_demand = filter(:period => ((pd) -> pd == y), total_demand)[1, :real_demand]
        # @info year_demand

        # Set up and run the dispatch simulation for this year
        # This function updates all_gc_results and all_prices in-place, and
        #   returns a boolean to determine whether the next year should be run
        run_next_year = run_annual_dispatch(y, year_portfolio, year_demand, ts_data, unit_specs, all_gc_results, all_prices)

        # @info "DISPATCH SIMULATION: YEAR $y COMPLETE."
        if y < current_pd + settings["num_dispatch_years"]
            # @info "RUN NEXT YEAR: $run_next_year"
        end

        if !run_next_year
            break
        end
    end

    # If no years ran correctly, throw an error and exit
    if size(all_gc_results)[1] == 0
        throw(ErrorException("No dispatch simulations could be run. Re-check your inputs."))
    end

    # Propagate the results dataframes out to the end of the projection horizon
    # Assume no change after the last modeled year
    all_gc_results, all_prices = propagate_all_results(fc_pd-1, all_gc_results, all_prices, current_pd)

    # Save the raw results
    save_raw_results(settings, all_prices, all_gc_results)

    # Create the final results set for use by the agent decision algorithm
    # @info "Postprocessing all dispatch results..."
    long_econ_results = postprocess_results(settings, all_year_system_portfolios, all_prices, all_gc_results, ts_data, unit_specs, fc_pd, current_pd)

    # Save a copy of the results to file
    savefile = joinpath(
                   pwd(),
                   "tmp",
                   "long_econ_results_$current_pd.csv"
               )
    # CSV.write(savefile, long_econ_results)

    return long_econ_results    

end


function load_ts_data(ts_file_dir, num_repdays)
    # Load the time-series demand and VRE data into dataframes
    load_data = CSV.read(joinpath(ts_file_dir, "timeseries_load_hourly.csv"), DataFrame)
    wind_data = CSV.read(joinpath(ts_file_dir, "timeseries_wind_hourly.csv"), DataFrame)
    solar_data = CSV.read(joinpath(ts_file_dir, "timeseries_pv_hourly.csv"), DataFrame)
    repdays_data = CSV.read(joinpath(ts_file_dir, "repDays_$num_repdays.csv"), DataFrame)

    ts_data = Dict(:load_data => load_data,
                   :wind_data => wind_data,
                   :solar_data => solar_data,
                   :repdays_data => repdays_data)

    return ts_data
end


function set_up_results_dfs()
    all_prices = DataFrame(y = Int[], d = Int[], h = Int[], price = Float64[])
    all_gc_results = DataFrame(y = Int[], d = Int[], h = Int[], unit_type = String[], gen = Float64[], commit = Int[])

    return all_prices, all_gc_results
end


function scale_load(ts_data, peak_demand)
    ts_data[:load_data][!, :Load] = ts_data[:load_data][!, :LoadShape] * peak_demand

    return ts_data
end


function set_up_load_repdays(ts_data)
    load_repdays = DataFrame()
    for day in ts_data[:repdays_data][!, :Day]
        load_repdays[!, Symbol(day)] = ts_data[:load_data][(24*day + 1):(24*(day + 1)), :Load]
    end

    ts_data[:load_repdays] = load_repdays

    return ts_data
end


function scale_wind_solar_data(ts_data, year_portfolio, unit_specs)
    wind_specs = filter(:unit_type => x -> x == "Wind", unit_specs)[1, :]
    solar_specs = filter(:unit_type => x -> x == "Solar", unit_specs)[1, :]

    # For non-zero wind and solar capacity, scale the WindShape series by the
    #   installed capacity to get total instantaneous VRE availability
    # If either wind or solar has 0 installed capacity, set its entire time
    #   series to 0.0
    if !isempty(filter(:unit_type => x -> x == "Wind", year_portfolio))
        ts_data[:wind_data][!, :Wind] = (ts_data[:wind_data][!, :WindShape]
                                            * filter(:unit_type => x -> x == "Wind", year_portfolio)[1, :num_units]
                                            * wind_specs[:capacity])
    else
        ts_data[:wind_data][!, :Wind] .= 0.0
    end

    if !isempty(filter(:unit_type => x -> x == "Solar", year_portfolio))
        ts_data[:solar_data][!, :Solar] = (ts_data[:solar_data][!, :SolarShape]
                                              * filter(:unit_type => x -> x == "Solar", year_portfolio)[1, :num_units]
                                              * solar_specs[:capacity])
    else
        ts_data[:solar_data][!, :Solar] .= 0.0
    end

    return ts_data
end


function set_up_wind_solar_repdays(ts_data)
    wind_repdays = DataFrame()
    solar_repdays = DataFrame()

    for day in ts_data[:repdays_data][!, :Day]
        wind_repdays[!, Symbol(day)] = ts_data[:wind_data][(24 * day + 1):(24 * (day + 1)), :Wind]
        solar_repdays[!, Symbol(day)] = ts_data[:solar_data][(24 * day + 1):(24 * (day + 1)), :Solar]
    end

    ts_data[:wind_repdays] = wind_repdays
    ts_data[:solar_repdays] = solar_repdays

    return ts_data
end


function set_up_model(ts_data, year_portfolio, unit_specs)
    # Create joined portfolio-unit_specs dataframe, to ensure consistent
    #   accounting for units which are actually present and consistent
    #   unit ordering
    portfolio_specs = innerjoin(unit_specs, year_portfolio, on = :unit_type)
    portfolio_specs[!, :unit_index] = 1:size(portfolio_specs)[1]

    # Set up wind_index and solar_index for easy filtering later
    wind_index = 0
    solar_index = 0
    if !isempty(filter(:unit_type => x -> x == "Wind", portfolio_specs))
        wind_index = filter(:unit_type => x -> x == "Wind", portfolio_specs)[1, :unit_index]
    end
    if !isempty(filter(:unit_type => x -> x == "Solar", portfolio_specs))
        solar_index = filter(:unit_type => x -> x == "Solar", portfolio_specs)[1, :unit_index]
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
    m = Model(CPLEX.Optimizer)

    # Set verbosity to lowest setting
    set_silent(m)

    # g: quantity generated (in MWh) for each unit type
    @variable(m, g[1:num_units, 1:num_days, 1:num_hours] >= 0)

    # c: number of units of each type committed in each hour
    @variable(m, c[1:num_units, 1:num_days, 1:num_hours] >= 0, Int)

    # Total generation per hour must be greater than or equal to the demand
    @constraint(m, mkt_equil[k=1:num_days, j=1:num_hours], sum(g[i, k, j] for i = 1:num_units) >= load_repdays[j, k])

    # Number of committed units per hour must be less than or equal to the total
    #   number of units of that type in the system
    for i = 1:num_units
        for k = 1:num_days
            for j = 1:num_hours
                @constraint(m, c[i, k, j] <= portfolio_specs[i, :num_units])
            end
        end
    end

    # Limit total generation per unit type each hour to the total capacity of all
    #   committed units of this type, with committed units subject to minimum and
    #   maximum power levels
    # Wind and solar are excluded
    for i = 1:num_units
        if !(portfolio_specs[i, :unit_type] in ["Wind", "Solar"])
            for k = 1:num_days
                for j = 1:num_hours
                    @constraint(m, g[i, k, j] <= c[i, k, j] .* portfolio_specs[i, :capacity] .* portfolio_specs[i, :CF] .* portfolio_specs[i, :PMAX])
                    @constraint(m, g[i, k, j] >= c[i, k, j] .* portfolio_specs[i, :capacity] .* portfolio_specs[i, :CF] .* portfolio_specs[i, :PMIN])
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
            for j = 1:num_hours-1
                # Ramp-up constraint
                @constraint(m, g[i, k, j+1] - g[i, k, j] <= c[i, k, j+1] .* portfolio_specs[i, :RUL] * portfolio_specs[i, :capacity] * portfolio_specs[i, :CF])
                # Ramp-down constraint
                @constraint(m, g[i, k, j+1] - g[i, k, j] >= (-1) * c[i, k, j+1] * portfolio_specs[i, :RDL] * portfolio_specs[i, :capacity] * portfolio_specs[i, :CF])
            end
        end
    end

    @objective(m, Min, sum(sum(sum(g[i, k, j] for j = 1:num_hours) for k = 1:num_days) .* (portfolio_specs[i, :VOM] + portfolio_specs[i, :FC_per_MWh] - portfolio_specs[i, :policy_adj_per_MWh]) for i = 1:num_units))

    return m, portfolio_specs
end


function solve_model(model; model_type="integral")
    if model_type == "integral"
        # Optimize!
        optimize!(model)
        status = string(termination_status.(model))
        # @info "$model_type model status: $status"
        gen_qty = nothing
        c = nothing
        if status == "OPTIMAL"
            gen_qty = value.(model[:g])
            c = value.(model[:c])
        end

        returns = model, status, gen_qty, c

    elseif model_type == "relaxed_integrality"
        optimize!(model)
        status = string(termination_status.(model))
        # @info "$model_type model status: $status"
        

        returns = model

    end

    return returns
end


function save_gc_results(ts_data, num_hours, num_units, y, gen_qty, c, all_gc_results, portfolio_specs)
    # Save generation and commitment results
    for k = 1:size(ts_data[:repdays_data])[1]
        for j = 1:num_hours
            for i = 1:num_units
                line = (y = y, d = k, h = j,
                        unit_type = portfolio_specs[i, :unit_type],
                        gen = gen_qty[i, k, j],
                        commit = round(Int, c[i, k, j]))
                push!(all_gc_results, line)
            end
        end
    end

    return all_gc_results
end


function reshape_shadow_price(model, ts_data, y, num_hours, num_days, all_prices)
    market_prices = reshape(shadow_price.(model[:mkt_equil]),
                     (size(ts_data[:repdays_data])[1], num_hours))

    # The result of shadow_price comes out weirdly reshaped, so this operation
    #   orders the data by hour (inner loop) then by day (outer loop)
    prices = zeros(1, num_hours*num_days)
    for j = 1:num_hours
        for k = 1:num_days
            line = (y = y, d = k, h = j,
                    price = (-1) * market_prices[k, j])
            push!(all_prices, line)
        end
    end

    return all_prices
end


function propagate_all_results(end_year, all_gc_results, all_prices, current_pd)
    final_dispatched_year = maximum(all_gc_results[!, :y])

    final_year_gc = filter(:y => x -> x == final_dispatched_year, all_gc_results)
    final_year_prices = filter(:y => x -> x == final_dispatched_year, all_prices)

    for y = final_dispatched_year+1:current_pd+end_year
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


function run_annual_dispatch(y, year_portfolio, peak_demand, ts_data, unit_specs, all_gc_results, all_prices)
    # Constants
    num_hours = 24

    # Scale the load data to the PD value for this year
    ts_data = scale_load(ts_data, peak_demand)

    # Set up representative days for load
    ts_data = set_up_load_repdays(ts_data)

    # Scale the wind and solar data according to the current year's total
    #   installed capacity
    ts_data = scale_wind_solar_data(ts_data, year_portfolio, unit_specs)

    # Set up representative days for wind and solar
    ts_data = set_up_wind_solar_repdays(ts_data)

    # @info "Setting up optimization model..."
    m, portfolio_specs = set_up_model(ts_data, year_portfolio, unit_specs)

    # @info "Optimization model set up."
    # @info string("Solving repday dispatch for year ", y, "...")

    # Create a copy of the model, to use later for the relaxed-integrality
    #   solution
    m_copy = copy(m)
    set_optimizer(m_copy, CPLEX.Optimizer)
    set_silent(m_copy)

    # Solve the integral optimization problem
    m, status, gen_qty, c = solve_model(m, model_type = "integral")

    if status == "OPTIMAL"
        # Save the generation and commitment results from the integral problem
        all_gc_results = save_gc_results(ts_data, num_hours, size(year_portfolio)[1], y, gen_qty, c, all_gc_results, portfolio_specs)

        # Set up a relaxed-integrality version of this model, to allow retrieval
        #   of dual values for the mkt_equil constraint
        # @info "Solving relaxed-integrality problem to get shadow prices..."
        undo = relax_integrality(m_copy)

        # Solve the relaxed-integrality model to compute the shadow prices
        m_copy = solve_model(m_copy, model_type = "relaxed_integrality")
        all_prices = reshape_shadow_price(
                         m_copy,
                         ts_data,
                         y,
                         num_hours,
                         size(ts_data[:repdays_data])[1],
                         all_prices
                     )

        # @info "Year $y dispatch run complete."
        run_next_year = true

    else
        # Dispatcher is unable to solve the current year: assume agents have
        #   not had a chance to address conditions in this year yet
        # Break the loop
        run_next_year = false
    end

    return run_next_year

end


function save_raw_results(settings, all_prices, all_gc_results)
    pfile = joinpath(
                pwd(),
                "tmp",
                "price_results.csv"
            )
    # CSV.write(pfile, all_prices)

    gcfile = joinpath(
                 pwd(),
                 "tmp",
                 "./gc_results.csv"
             )
    # CSV.write(gcfile, all_gc_results)
end


function pivot_gc_results(all_gc_results, all_prices, repdays_data)
    # @info "Postprocessing results..."
    # Pivot generation data by unit type
    # Output format: y, d, h, Wind, Solar, ..., AdvancedNuclear
    g_pivot = select(all_gc_results, Not(:commit))
    g_pivot = unstack(g_pivot, :unit_type, :gen)
    g_pivot = innerjoin(g_pivot, all_prices, on = [:y, :d, :h])
    g_pivot = innerjoin(g_pivot, select(repdays_data, [:index, :Probability]), on = [:d => :index])

    # Pivot commitment data by unit type
    # Output format: y, d, h, Wind, ..., AdvancedNuclear
    c_pivot = select!(all_gc_results, Not(:gen))
    c_pivot = unstack(c_pivot, :unit_type, :commit)

    return g_pivot, c_pivot
end


function set_up_ppx_revenue(g_pivot, ts_data, unit_specs)
    # Add revenue columns to g_pivot
    for unit_type in unit_specs[!, :unit_type]
        # Retrieve this unit's specs for easy reference
        unit_spec = filter(:unit_type => x -> x == unit_type, unit_specs)
        transform!(g_pivot, [Symbol(unit_type), :Probability] => ((gen, prob) -> gen .* prob .* 365) => Symbol(string(unit_type, "total_gen")))
        transform!(g_pivot, [Symbol(string(unit_type, "total_gen")), :price] => ((gen, price) -> gen .* price) => Symbol(string(unit_type, "rev")))
        transform!(g_pivot, [Symbol(string(unit_type, "total_gen"))] => ((gen) -> gen .* (unit_spec[1, :VOM] + unit_spec[1, :FC_per_MWh] .* unit_spec[1, :heat_rate])) => Symbol(string(unit_type, "opcost")))
        transform!(g_pivot, [Symbol(string(unit_type, "rev")), Symbol(string(unit_type, "opcost"))] => ((rev, opcost) -> rev - opcost) => Symbol(string(unit_type, "profit")))
    end

    return g_pivot

end


function reorganize_g_pivot(g_pivot, unit_specs)
    # Reorganize g_pivot in order to retrieve a DataFrame of annual profit values
    #   by unit type
    fields_to_select = [:y, :d, :h]
    for unit_type in unit_specs[!, :unit_type]
        append!(fields_to_select, [Symbol(string(unit_type, :profit))])
    end
    g_unpivot = select(g_pivot, fields_to_select)

    # Rename unit-specific columns from <type>profit to just <type> for clarity
    for unit_type in unit_specs[!, :unit_type]
        rename!(g_unpivot, Symbol(string(unit_type, "profit")) => Symbol(unit_type))
    end

    g_unpivot = stack(g_unpivot, 4:size(g_unpivot)[2])

    # Rename new columns from automatic names to more descriptive ones
    rename!(g_unpivot, :variable => :unit_type, :value => :profit)

    # Pivot g_pivot so that generation results from each year (by unit type) are grouped together
    fields_to_select = [:y, :d, :h]
    for unit_type in unit_specs[!, :unit_type]
        append!(fields_to_select, [Symbol(string(unit_type, "total_gen"))])
    end
    g_gen = select(g_pivot, fields_to_select)

    # Rename unit-specific columns from <type>total_gen to just <type>
    for unit_type in unit_specs[!, :unit_type]
        rename!(g_gen, Symbol(string(unit_type, "total_gen")) => Symbol(unit_type))
    end

    g_gen = stack(g_gen, Not([:y, :d, :h]))

    rename!(g_gen, :variable => :unit_type, :value => :ydh_gen)

    g_gen = combine(groupby(g_gen, [:y, :unit_type]), :ydh_gen => sum => :total_gen)

    return g_unpivot, g_gen

end


function compute_final_profit(g_unpivot, system_portfolios, unit_specs, fc_pd, current_pd)
    # Combine the system portfolios into a single dataframe with indicator
    #   column :y
    all_year_portfolios = create_all_year_portfolios(system_portfolios, fc_pd, current_pd)

    # Calculate operating and net profit per unit per unit type
    g_unpivot = innerjoin(g_unpivot, all_year_portfolios, on = [:unit_type, :y])
    transform!(g_unpivot, [:profit, :num_units] => ((profit, num_units) -> profit ./ num_units) => :ProfitPerUnit)

    g_profitsum = combine(groupby(g_unpivot, [:y, :unit_type]), :ProfitPerUnit => sum => :annual_profit)
    g_profitsum[!, :FOM_per_unit] .= 0.0

    for i = 1:size(g_profitsum)[1]
        unit_type = g_profitsum[i, :unit_type]
        unit_spec = filter(
                        :unit_type => x -> x == unit_type,
                        unit_specs
                    )
        g_profitsum[i, :FOM_per_unit] = unit_spec[1, :FOM] * unit_spec[1, :capacity] * 1000
    end

    # Re-pivot the final result to show annual net profit per unit for each
    #   unit type and year
    transform!(g_profitsum, [:annual_profit, :FOM_per_unit] =>
               ((profit, FOM) -> profit - FOM) =>
               :net_profit_per_unit)

    final_profit_pivot = unstack(
                             select(
                                 g_profitsum,
                                 Not([:annual_profit, :FOM_per_unit])),
                             :unit_type, :net_profit_per_unit
                         )

    return final_profit_pivot

end


function create_all_year_portfolios(system_portfolios, fc_pd, current_pd)
    # Combine the system portfolios into a single dataframe with indicator
    #   column :y

    # Explicit dispatch results
    all_year_portfolios = DataFrame()
    for key in keys(system_portfolios)
        df = system_portfolios[key]
        df[!, :y] .= key
        append!(all_year_portfolios, df)
    end

    # Extend dispatch results by assuming no change after last dispatch year
    last_dispatch_year = maximum([key for key in keys(system_portfolios)])
    for i = last_dispatch_year+1:current_pd+fc_pd
        df = system_portfolios[length(keys(system_portfolios))-1]
        df[!, :y] .= i
        append!(all_year_portfolios, df)
    end

    return all_year_portfolios

end


function compute_final_gen(g_gen, system_portfolios, fc_pd, current_pd)
    all_year_portfolios = create_all_year_portfolios(system_portfolios, fc_pd, current_pd)

    # Calculate operating and net profit per unit per unit type
    g_gen = innerjoin(g_gen, all_year_portfolios, on = [:unit_type, :y])
    transform!(g_gen, [:total_gen, :num_units] => ((total_gen, num_units) -> total_gen ./ num_units) => :GenPerUnit)

    final_gen_pivot = g_gen

    return final_gen_pivot

end


function postprocess_long_results(g_pivot, settings, system_portfolios, unit_specs, fc_pd, current_pd)
    # Organize data
    long_rev_results = deepcopy(g_pivot)

    # Get a single long dataframe with unit numbers by type for each year
    all_year_portfolios = create_all_year_portfolios(system_portfolios, fc_pd, current_pd)

    # Get a list of all long_rev_results columns whose names match entries in
    #   unit_specs
    types_list = unit_specs[!, :unit_type]
    extant_types = [col for col in names(long_rev_results) if col in types_list]
    
    long_rev_results = stack(long_rev_results, extant_types)
    rename!(long_rev_results, :variable => :unit_type, :value => :gen)
    short_unit_specs = select(unit_specs, [:unit_type, :capacity, :VOM, :FC_per_MWh, :policy_adj_per_MWh])
    long_rev_results = innerjoin(long_rev_results, short_unit_specs, on = :unit_type)

    # Create dataframe of subsidy values
    subs = DataFrame(
               unit_type = String[],
               subs_val = Float64[]
                    )
    for unit_type in unique(all_year_portfolios[!, :unit_type])
        subs_val = 1  # multiplier
        if (occursin("C2N", unit_type)) || (unit_type == "AdvancedNuclear")
            subs_val = settings["C2N_subsidy"]
        end
        push!(subs, [unit_type subs_val])
    end

    # Append unit number data to long_rev_results
    long_rev_results = innerjoin(long_rev_results, all_year_portfolios, on = [:y, :unit_type])

    # Append subsidy data to long_rev_results
    long_rev_results = innerjoin(long_rev_results, subs, on = :unit_type)

    # Calculate revenues
    transform!(long_rev_results, [:gen, :price, :subs_val, :Probability, :num_units] => ((gen, price, subs, prob, num_units) -> gen .* price .* subs .* prob .* 365 ./ num_units) => :annualized_rev_perunit)

    # Calculate VOM
    transform!(long_rev_results, [:gen, :VOM, :Probability, :num_units] => ((gen, VOM, prob, num_units) -> gen .* VOM .* prob .* 365 ./ num_units) => :annualized_VOM_perunit)

    # Calculate fuel cost
    transform!(long_rev_results, [:gen, :FC_per_MWh, :Probability, :num_units] => ((gen, fc, prob, num_units) -> gen .* fc .* prob .* 365 ./ num_units) => :annualized_FC_perunit)

    # Calculate policy adjustment
    transform!(long_rev_results, [:gen, :policy_adj_per_MWh, :Probability, :num_units] => ((gen, adj, prob, num_units) -> gen .* adj .* prob .* 365 ./ num_units) => :annualized_policy_adj_perunit)

    return long_rev_results

end 


function postprocess_results(settings, system_portfolios, all_prices, all_gc_results, ts_data, unit_specs, fc_pd, current_pd)
    # Pivot the generation and commitment results
    g_pivot, c_pivot = pivot_gc_results(all_gc_results, all_prices, ts_data[:repdays_data])

    long_econ_results = postprocess_long_results(g_pivot, settings, system_portfolios, unit_specs, fc_pd, current_pd)

    return long_econ_results
end

end
