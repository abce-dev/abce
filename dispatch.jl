module Dispatch

using Logging, CSV, DataFrames, JuMP, GLPK, XLSX, Logging, CPLEX


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
    all_prices_df = DataFrame(y = Int[], d = Int[], h = Int[], price = Float64[])
    all_gc_results_df = DataFrame(y = Int[], d = Int[], h = Int[], unit_type = String[], gen = Float64[], commit = Int[])

    return all_prices_df, all_gc_results_df
end


function scale_load(ts_data, peak_demand)
    ts_data[:load_data][!, :Load] = ts_data[:load_data][!, :LoadShape] * peak_demand

    return load_data
end


function set_up_load_repdays(ts_data)
    load_repdays = DataFrame()
    for day in repdays_data[!, :Day]
        load_repdays[!, Symbol(day)] = ts_data[:load_data][(24*day + 1):(24*(day + 1)), :Load]
    end

    return load_repdays
end


function scale_wind_solar_data(ts_data, year_portfolio, unit_specs)
    wind_specs = filter(:unit_type => x -> x == "Wind", unit_specs)[1, :]
    solar_specs = filter(:unit_type => x -> x == "Solar", unit_specs)[1, :]

    ts_data[:wind_data][!, :Wind] = (ts_data[:wind_data][!, :WindShape]
                                        * filter(:unit_type => x -> x == "Wind", year_portfolio)[1, :num_units]
                                        * wind_specs[:capacity])
    ts_data[:solar_data][!, :Solar] = (ts_data[:solar_data][!, :SolarShape]
                                          * filter(:unit_type => x -> x == "Solar", year_portfolio)[1, :num_units]
                                          * solar_specs[:capacity])

    return ts_data
end


function set_up_wind_solar_repdays(ts_data)
    wind_repdays = DataFrame()
    solar_repdays = DataFrame()

    for day in ts_data[:repdays][!, :Day]
        ts_data[:wind_repdays][!, Symbol(day)] = ts_data[:wind_data][(24 * day + 1):(24 * (day + 1)), :Wind]
        ts_data[:solar_repdays][!, Symbol(day)] = ts_data[:solar_data][(24 * day + 1):(24 * (day + 1)), :Solar]
    end    

    return ts_data
end


function set_up_model(ts_data, year_portfolio, unit_specs)
    # Helpful named constants
    num_units = size(year_portfolio)[1]
    num_days = size(ts_data[:repdays])[1]
    num_hours = 24

    # Break out timeseries data sets for convenience
    load_data = ts_data[:load_data]
    wind_data = ts_data[:wind_data]
    solar_data = ts_data[:solar_data]
    repdays_data = ts_data[:repdays_data]

    # Initialize JuMP model
    m = Model(with_optimizer(CPLEX.Optimizer))

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
                @constraint(m, c[i, k, j] <= year_portfolio[i, :num_units])
            end
        end
    end

    # Limit total generation per unit type each hour to the total capacity of all
    #   committed units of this type, with committed units subject to minimum and
    #   maximum power levels
    for i = 3:num_units
        for k = 1:num_days
            for j = 1:num_hours
                @constraint(m, g[i, k, j] <= c[i, k, j] .* unit_specs[i, :CAP] .* unit_specs[i, :CF] .* unit_specs[i, :PMAX])
                @constraint(m, g[i, k, j] >= c[i, k, j] .* unit_specs[i, :CAP] .* unit_specs[i, :CF] .* unit_specs[i, :PMIN])
            end
        end
    end

    # Limit solar and wind generation to their actual hourly availability
    for k = 1:num_days
        for j = 1:num_hours
            @constraint(m, g[1, k, j] <= wind_repdays[j, k])
            @constraint(m, g[2, k, j] <= solar_repdays[j, k])
        end
    end

    # Ramping constraints
    for i = 1:num_units
        for k = 1:num_days
            for j = 1:num_hours-1
                # Ramp-up constraint
                @constraint(m, g[i, k, j+1] - g[i, k, j] <= c[i, k, j+1] .* unit_specs[i, :RUL] * unit_specs[i, :CAP] * unit_specs[i, :CF])
                # Ramp-down constraint
                @constraint(m, g[i, k, j+1] - g[i, k, j] >= (-1) * c[i, k, j+1] * unit_specs[i, :RDL] * unit_specs[i, :CAP] * unit_specs[i, :CF])
            end
        end
    end

    @objective(m, Min, sum(sum(sum(g[i, k, j] for j = 1:num_hours) for k = 1:num_days) .* (unit_specs[i, :VOM] + unit_specs[i, :FC] .* unit_specs[i, :HR]) for i = 1:num_units))

    return m
end


function solve_model(model; model_type="integral")
    if model_type == "integral"
        # Optimize!
        optimize!(model)
        status = termination_status.(model)
        @info "$model_type model status: $status"
        gen_qty = value.(model[:g])
        c = value.(model[:c])

        returns = model, gen_qty, c

    else if model_type == "relaxed_integrality"
        optimize!(model)
        status = termination_status.(model)
        @info "$model_type model status: $status"
        

        returns = model

    end

    return returns
end


function save_gc_results(ts_data, num_hours, num_units, y, gen_qty, c, all_gc_results)
    # Save generation and commitment results
    for k = 1:size(ts_data[:repdays_data])[1]
        for j = 1:num_hours
            for i = 1:num_units
                line = (y = y, d = k, h = j,
                        unit_type = unit_specs[i, :UNIT_TYPE],
                        gen = gen_qty[i, k, j],
                        commit = c[i, k, j])
                push!(all_gc_results, line)
            end
        end
    end

    return all_gc_results
end


function reshape_shadow_price(model, num_hours, num_days, all_prices)
    market_prices = reshape(shadow_price.(model[:mkt_equil]),
                     (size(ts_data[:repdays_data])[1], num_hours))

    # The result of shadow_price comes out weirdly reshaped, so this operation
    #   orders the data by hour (inner loop) then by day (outer loop)
    prices = zeros(1, num_hours*num_days)
    for j = 1:num_hours
        for k = 1:num_days
            line = (y = y, d = k, h = j,
                    price = (-1) * mkt_prices[k, j])
            push!(all_prices, line)
        end
    end

    return all_prices
end


function run_annual_dispatch(y, year_portfolio, ts_data, unit_specs)
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

    @info "Setting up optimization model..."
    m = set_up_model(ts_data, year_portfolio, unit_specs)

    @info "Optimization model set up."
    @info string("Solving repday dispatch for year ", y, "...")

    # Create a copy of the model, to use later for the relaxed-integrality
    #   solution
    m_copy = copy(m)
    set_optimizer(m_copy, CPLEX.Optimizer)

    # Solve the integral optimization problem
    m, gen_qty, c = solve_model(m, model_type = "integral")

    # Save the generation and commitment results from the integral problem
    all_gc_results = save_gc_results(ts_data, num_hours, num_units, y, gen_qty, c, all_gc_results)

    # Set up a relaxed-integrality version of this model, to allow retrieval
    #   of dual values for the mkt_equil constraint
    @info "Solving relaxed-integrality problem to get shadow prices..."
    undo = relax_integrality(m_copy)

    # Solve the relaxed-integrality model to compute the shadow prices
    m_copy = solve_model(m_copy, model_type = "relaxed_integrality")
    all_prices = save_shadow_price(
                     m_copy,
                     num_hours,
                     size(ts_data[:repdays_data])[1],
                     all_prices
                 )

    @info "Year $y dispatch run complete."

    return all_gc_results, all_prices

end


function save_raw_results(all_prices, all_gc_results)
    pfile = "./price_results.csv"
    CSV.write(pfile, all_prices)

    gcfile = "./gc_results.csv"
    CSV.write(gcfile, all_gc_results)
end


function pivot_gc_results(all_gc_results, all_prices)
    @info "Postprocessing results..."
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
    for unit_type in unit_specs[!, :UNIT_TYPE]
        # Retrieve this unit's specs for easy reference
        unit_spec = filter(:UNIT_TYPE => x -> x == unit_type, unit_specs)
        transform!(g_pivot, [Symbol(unit_type), :Probability] => ((gen, prob) -> gen .* prob .* 365) => Symbol(string(unit_type, "norm")))
        transform!(g_pivot, [Symbol(string(unit_type, "norm")), :price] => ((gen, price) -> gen .* price) => Symbol(string(unit_type, "rev")))
        transform!(g_pivot, [Symbol(string(unit_type, "norm"))] => ((gen) -> gen .* (unit_spec[1, :VOM] + unit_spec[1, :FC] .* unit_spec[1, :HR])) => Symbol(string(unit_type, "opcost")))
        transform!(g_pivot, [Symbol(string(unit_type, "rev")), Symbol(string(unit_type, "opcost"))] => ((rev, opcost) -> rev - opcost) => Symbol(string(unit_type, "profit")))
    end

    return g_pivot

end


function reorganize_g_pivot(g_pivot, unit_specs)
    # Reorganize g_pivot in order to retrieve a DataFrame of annual profit values
    #   by unit type
    fields_to_select = [:y, :d, :h]
    for unit_type in unit_specs[!, :UNIT_TYPE]
        append!(fields_to_select, [Symbol(string(unit_type, :profit))])
    end
    g_unpivot = select(g_pivot, fields_to_select)

    # Rename unit-specific columns from <type>profit to just <type> for clarity
    for unit_type in unit_specs[!, :UNIT_TYPE]
        rename!(g_unpivot, Symbol(string(unit_type, "profit")) => Symbol(unit_type))
    end

    g_unpivot = stack(g_unpivot, 4:size(g_unpivot)[2])

    # Rename new columns from automatic names to more descriptive ones
    rename!(g_unpivot, :variable => :unit_type, :value => :profit)

    return g_pivot

end


function compute_net_profit(g_pivot, year_portfolio, unit_specs)
    # Calculate operating and net profit per unit per unit type
    g_unpivot = innerjoin(g_unpivot, portfolio, on = :unit_type)
    transform!(g_unpivot, [:profit, :num_units] => ((profit, num_units) -> profit ./ num_units) => :ProfitPerUnit)

    g_profitsum = combine(groupby(g_unpivot, [:y, :unit_type]), :ProfitPerUnit => sum => :annual_profit)
    g_profitsum[!, :FOM_per_unit] .= 0.0

    for i = 1:size(g_profitsum)[1]
        unit_type = g_profitsum[i, :unit_type]
        unit_spec = filter(
                        :UNIT_TYPE => x -> x == unit_type,
                        unit_specs
                    )
        num_curr_units = filter(
                             :unit_type => x -> x == unit_type,
                             portfolio
                         )[1, :num_units]
        g_profitsum[i, :FOM_per_unit] = unit_spec[1, :FOM] * unit_spec[1, :CAP] * 1000
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


function postprocess_results(all_prices, all_gc_results, ts_data, unit_specs)
    # Pivot the generation and commitment results
    g_pivot, c_pivot = pivot_gc_results(all_gc_results, all_prices)

    # Set up revenue columns in g_pivot
    g_pivot = set_up_ppx_revenue(g_pivot, ts_data, unit_specs)

    # Reorganize the g_pivot dataframe
    g_pivot = reorganize_g_pivot(g_pivot, unit_specs)

    # Compute operating and net profit, including by unit type
    final_profit_pivot = compute_final_profit(g_pivot, portfolio, unit_specs)

    CSV.write("unit_profit_summary.csv", final_profit_pivot)

end

end
