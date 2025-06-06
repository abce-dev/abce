##########################################################################
# Copyright 2023 Argonne National Laboratory
#
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


module Dispatch

using Requires, Logging, CSV, DataFrames, JuMP, GLPK, Cbc, XLSX, SQLite, HiGHS, CPLEX


# Initialize this module, with CPLEX as an optional library if available
#function __init__()
#    @require CPLEX = "a076750e-1247-5638-91d2-ce28b192dca0" @eval using CPLEX
#end


function execute_dispatch_economic_projection(
    CLI_args,
    db,
    settings,
    fc_pd,
    total_demand,
    unit_specs,
    system_portfolios;
    run_mode="forecast",
    downselection_mode="scenario_reduction"
)
    @debug string(
        "Running the dispatch simulation for ",
        settings["dispatch"]["num_dispatch_years"],
        " years...",
    )

    all_repdays, all_grc_results, all_price_results = handle_annual_dispatch(
        db,
        settings,
        CLI_args,
        system_portfolios,
        total_demand,
        unit_specs;
        run_mode=run_mode,
        downselection_mode=downselection_mode
    )

    long_econ_results, dispatch_results = postprocess_results(
        settings,
        all_grc_results,
        all_price_results,
        all_repdays,
        system_portfolios,
        unit_specs,
        CLI_args["current_pd"],
        fc_pd,
    )

    return long_econ_results, dispatch_results
end


function get_system_portfolios(db, settings, start_year, unit_specs)
    # Set up portfolio dictionaries
    @debug "Setting up dispatch portfolios..."
    system_portfolios = Dict()

    end_year = (
        start_year +
        convert(Int64, settings["dispatch"]["num_dispatch_years"]) - 1
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
    for y = minimum(keys(system_portfolios)):maximum(keys(system_portfolios))
        for unit_type in unit_specs[!, :unit_type]
            if !in(unit_type, system_portfolios[y][!, :unit_type])
                push!(system_portfolios[y], (unit_type, 1))
            end
        end
    end

    return system_portfolios
end


function set_up_ts_data(settings, ts_data_dir, current_pd, fc_pd, downselection_mode)
    if (settings["dispatch"]["downselection"] == "exact") || (downselection_mode == "exact")
        num_days = 365
    else
        num_days = settings["dispatch"]["num_repdays"]
    end

    ts_data = load_ts_data(
        ts_data_dir,
        current_pd,
        fc_pd,
        num_repdays=num_days
    )

    return ts_data
end


function handle_annual_dispatch(
    db,
    settings,
    CLI_args,
    system_portfolios,
    total_demand,
    unit_specs;
    run_mode="forecast",
    downselection_mode="scenario_reduction"
)
    current_pd = CLI_args["current_pd"]

    ts_data_dir = joinpath(
        CLI_args["inputs_path"],
        settings["file_paths"]["timeseries_data_dir"],
    )

    all_grc_results = set_up_grc_results_df()
    all_prices = set_up_prices_df()

    if run_mode == "forecast"
        num_years = settings["dispatch"]["num_dispatch_years"]
    elseif run_mode == "current"
        num_years = 1
    end

    all_repdays = nothing

    # Run the annual dispatch for the user-specified number of dispatch years
    for y = current_pd:(current_pd + num_years - 1)
        @info "\n\nDISPATCH SIMULATION: YEAR $y"

        # Set up the timeseries data for this year
        ts_data = set_up_ts_data(
            settings,
            ts_data_dir,
            CLI_args["current_pd"],
            y,
            downselection_mode
        )

        y_repdays = deepcopy(ts_data[:repdays_data])
        y_repdays[!, :y] .= y

        if all_repdays == nothing
            all_repdays = y_repdays
        else
            all_repdays = vcat(all_repdays, y_repdays)
        end

        # Select the current year's expected portfolio
        year_portfolio = system_portfolios[y]

        # Determine appropriate total demand for this year
        year_demand =
            filter(:period => ((pd) -> pd == y), total_demand)[1, :real_demand]

        # Set up and run the dispatch simulation for this year
        # This function updates all_grc_results and all_prices in-place, and
        #   returns a boolean to determine whether the next year should be run
        results = run_annual_dispatch(
            settings,
            y,
            year_portfolio,
            year_demand,
            ts_data,
            unit_specs;
            run_mode=run_mode,
            downselection_mode=downselection_mode
        )

        # Save new generation, commitment, and price results
        all_grc_results = vcat(all_grc_results, results[:new_grc_results])
        all_prices = vcat(all_prices, results[:new_prices])

        if run_mode == "current"
            save_summary_statistics(
                db,
                CLI_args["current_pd"],
                results[:summary_statistics]
            )
        end

        @debug "DISPATCH SIMULATION: YEAR $y COMPLETE."

    end

    return all_repdays, all_grc_results, all_prices

end


function load_ts_data(ts_file_dir, base_pd, fc_pd; num_repdays=nothing)
    # Load the time-series demand and VRE data into dataframes
    load_data =
        CSV.read(joinpath(ts_file_dir, "timeseries_load_hourly.csv"), DataFrame)
    wind_data =
        CSV.read(joinpath(ts_file_dir, "timeseries_wind_hourly.csv"), DataFrame)
    solar_data =
        CSV.read(joinpath(ts_file_dir, "timeseries_pv_hourly.csv"), DataFrame)
    reg_data = 
        CSV.read(joinpath(ts_file_dir, "timeseries_reg_hourly.csv"), DataFrame)
    spin_data = 
        CSV.read(joinpath(ts_file_dir, "timeseries_spin_hourly.csv"), DataFrame)
    nspin_data = 
        CSV.read(joinpath(ts_file_dir, "timeseries_nspin_hourly.csv"), DataFrame)

    if num_repdays == nothing
        repdays_data = nothing
    elseif num_repdays == 365
        repdays_data = CSV.read(
            joinpath(
                ts_file_dir,
                "repDays_365.csv",
            ),
            DataFrame,
        )
    else
        repdays_data = CSV.read(
            joinpath(
                ts_file_dir,
                string(
                    "bp_",
                    base_pd,
                    "_fp_",
                    fc_pd,
                    "_repDays_$num_repdays.csv",
                ),
            ),
            DataFrame
        )
    end

    ts_data = Dict(
        :load_data => load_data,
        :wind_data => wind_data,
        :solar_data => solar_data,
        :reg_data => reg_data,
        :spin_data => spin_data,
        :nspin_data => nspin_data,
        :repdays_data => repdays_data,
    )

    return ts_data
end

function set_up_grc_results_df()
    all_grc_results = DataFrame(
        y = Int[],
        d = Int[],
        h = Int[],
        unit_type = String[],
        gen = Float64[],
        reg = Float64[],
        spin = Float64[],
        nspin = Float64[],
        commit = Int[],
        su = Int[],
        sd = Int[],
    )

    return all_grc_results
end


function set_up_prices_df()
    prices_df = DataFrame(
        y = Int[], 
        d = Int[], 
        h = Int[], 
        lambda = Float64[], 
        reg_rmp = Float64[], 
        spin_rmp = Float64[], 
        nspin_rmp = Float64[]
    )

    return prices_df
end


function scale_load(ts_data, peak_demand)
    ts_data[:load_data][!, :Load] =
        (ts_data[:load_data][!, :LoadShape] * peak_demand)

    return ts_data
end


function scale_AS_data(ts_data, peak_demand, init_peak_demand)
    ts_data[:reg_data][!, :AS_qty] = ts_data[:reg_data][!, :ReqReg] * peak_demand / init_peak_demand
    ts_data[:spin_data][!, :AS_qty] = ts_data[:spin_data][!, :ReqSR] * peak_demand / init_peak_demand
    ts_data[:nspin_data][!, :AS_qty] = ts_data[:nspin_data][!, :ReqNSR] * peak_demand / init_peak_demand

    return ts_data
end


function set_repdays_params(settings, ts_data)
    if settings["dispatch"]["downselection"] == "exact"
        if size(ts_data[:load_data])[1] <= 24
            num_days = 1
            num_hours = Int(round(size(ts_data[:load_data])[1]))
        else
            num_days = Int(round(size(ts_data[:load_data])[1] / 24))
            num_hours = 24
        end
    else
        num_days = settings["dispatch"]["num_repdays"]
        num_hours = 24
    end

    return num_days, num_hours

end


function set_up_load_repdays(downselection_mode, num_days, num_hours, ts_data)
    load_repdays = DataFrame()

    if downselection_mode == "exact"
        for day = 1:num_days
            load_repdays[!, Symbol(day)] =
                ts_data[:load_data][(num_hours * (day-1) + 1):(num_hours * day), :Load]
        end
    else
        for day in ts_data[:repdays_data][!, :Day]
            load_repdays[!, Symbol(day)] =
                ts_data[:load_data][(num_hours * (day-1) + 1):(num_hours * day), :Load]
        end
    end

    ts_data[:load_repdays] = load_repdays

    return ts_data
end


function set_up_AS_repdays(downselection_mode, num_days, num_hours, ts_data)
    reg_repdays = DataFrame()
    spin_repdays = DataFrame()
    nspin_repdays = DataFrame()

    if downselection_mode == "exact"
        for day = 1:num_days
            reg_repdays[!, Symbol(day)] = ts_data[:reg_data][(num_hours * (day - 1) + 1):(num_hours * day), :AS_qty]
            spin_repdays[!, Symbol(day)] = ts_data[:spin_data][(num_hours * (day - 1) + 1):(num_hours * day), :AS_qty]
            nspin_repdays[!, Symbol(day)] = ts_data[:nspin_data][(num_hours * (day - 1) + 1):(num_hours * day), :AS_qty]
        end
    else
        for day in ts_data[:repdays_data][!, :Day]
            reg_repdays[!, Symbol(day)] = ts_data[:reg_data][(num_hours * day + 1):(num_hours * (day + 1)), :AS_qty]
            spin_repdays[!, Symbol(day)] = ts_data[:spin_data][(num_hours * day + 1):(num_hours * (day + 1)), :AS_qty]
            nspin_repdays[!, Symbol(day)] = ts_data[:nspin_data][(num_hours * day + 1):(num_hours * (day + 1)), :AS_qty]
        end
    end

    ts_data[:reg_repdays] = reg_repdays
    ts_data[:spin_repdays] = spin_repdays
    ts_data[:nspin_repdays] = nspin_repdays

    return ts_data
end


function set_up_wind_solar_repdays(downselection_mode, num_days, num_hours, ts_data)
    wind_repdays = DataFrame()
    solar_repdays = DataFrame()

    if downselection_mode == "exact"
        days = 1:num_days
    else
        days = ts_data[:repdays_data][!, :Day]
    end

    for day in days
        wind_repdays[!, Symbol(day)] =
            ts_data[:wind_data][(num_hours * (day-1) + 1):(num_hours * day), :WindShape]
        solar_repdays[!, Symbol(day)] =
            ts_data[:solar_data][(num_hours * (day-1) + 1):(num_hours * day), :SolarShape]
    end

    ts_data[:wind_repdays] = deepcopy(wind_repdays)
    ts_data[:solar_repdays] = deepcopy(solar_repdays)

    return ts_data
end


function set_up_model(settings, num_days, num_hours, ts_data, year_portfolio, unit_specs; gen_data, commit_data)
    # Create joined portfolio-unit_specs dataframe, to ensure consistent
    #   accounting for units which are actually present and consistent
    #   unit ordering
    portfolio_specs = innerjoin(unit_specs, year_portfolio, on = :unit_type, makeunique=true)
    portfolio_specs[!, :unit_index] = 1:size(portfolio_specs)[1]

    convnuc_index = nothing
    if !isempty(filter(:unit_type => x -> x == "conventional_nuclear", portfolio_specs))
        convnuc_index = filter(:unit_type => x -> x == "conventional_nuclear", portfolio_specs)[1, :unit_index]
    end

    # Helpful named constants
    num_units = size(portfolio_specs)[1]

    # Break out timeseries data sets for convenience
    load_repdays = ts_data[:load_repdays]
    wind_repdays = ts_data[:wind_repdays]
    solar_repdays = ts_data[:solar_repdays]
    reg_repdays = ts_data[:reg_repdays]
    spin_repdays = ts_data[:spin_repdays]
    nspin_repdays = ts_data[:nspin_repdays]

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

    # r: quantity (MWh) reserved for frequency regulation (combined up/down)
    @variable(m, r[1:num_units, 1:num_days, 1:num_hours] >= 0)

    # sr: quantity (MWh) reserved for spinning reserve
    @variable(m, sr[1:num_units, 1:num_days, 1:num_hours] >= 0)

    # nsr: quantity (MWh) reserved for non-spinning reserve
    @variable(m, nsr[1:num_units, 1:num_days, 1:num_hours] >= 0)

    # c: number of units of each type committed in each hour
    @variable(m, c[1:num_units, 1:num_days, 1:num_hours] >= 0, Int)

    # su: number of units of each type started up in each hour
    @variable(m, su[1:num_units, 1:num_days, 1:num_hours] >= 0)

    # sd: number of units of each type shut down in each hour
    @variable(m, sd[1:num_units, 1:num_days, 1:num_hours] >= 0)

    # ens: penalty variable for energy not served in each hour
    @variable(m, ens[1:num_days, 1:num_hours] >= 0)

    # rns: penalty variable for regulation not served in each hour
    @variable(m, rns[1:num_days, 1:num_hours] >= 0)

    # sns: penalty variable for spinning reserve not served in each hour
    @variable(m, sns[1:num_days, 1:num_hours] >= 0)

    # nsns: penalty variable for non-spinning reserve not served in each hour
    @variable(m, nsns[1:num_days, 1:num_hours] >= 0)

    # Total generation per hour plus energy not served must be greater than
    #   or equal to the demand
    @constraint(
        m,
        mkt_equil[k = 1:num_days, j = 1:num_hours],
        sum(g[i, k, j] for i = 1:num_units) + ens[k, j] == load_repdays[j, k]
    )

    # Total regulation per hour plus regulation not served must be greater 
    #   than or equal to the requirement
    @constraint(
        m,
        reg_mkt_equil[k = 1:num_days, j = 1:num_hours],
        sum(r[i, k, j] for i = 1:num_units) + rns[k, j] == reg_repdays[j, k]
    )

    # Total spinning reserve per hour plus spinning reserve not served must be
    #   greater than or equal to the requirement
    @constraint(
        m,
        spin_mkt_equil[k = 1:num_days, j = 1:num_hours],
        sum(sr[i, k, j] for i = 1:num_units) + sns[k, j] == spin_repdays[j, k]
    )

    # Total non-spinning reserve per hour plus non-spinning reserve not served
    #   must be greater than or equal to the requirement
    @constraint(
        m,
        nspin_mkt_equil[k = 1:num_days, j = 1:num_hours],
        sum(nsr[i, k, j] for i = 1:num_units) + nsns[k, j] == nspin_repdays[j, k]
    )

    # Number of committed units per hour must be less than or equal to the
    #   total number of units of that type in the system
    for i = 1:num_units
        if portfolio_specs[i, :is_VRE] == 0
            for k = 1:num_days
                for j = 1:num_hours
                    @constraint(m, c[i, k, j] <= portfolio_specs[i, :esc_num_units])
                end
            end
        end
    end

    # If no gen_data and commit_data are provided, assume that the time series 
    #   starts with the first hour of the year or repday year
    if (gen_data == nothing) || (commit_data == nothing)
        for i = 1:num_units
            if portfolio_specs[i, :is_VRE] == 0
            # on the first hour of the "year", all units committed are counted as being started up this hour
            @constraint(m, su[i, 1, 1] == c[i, 1, 1])
            end
        end
    else
        # If gen_data and commit_data are available, match this slice's
        #   first-hour results to last slice's last-hour results
        for i = 1:num_units
            # Match generation: ramp-up
            @constraint(
                m,
                (
                    g[i, 1, 1] - gen_data[i] <=
                    c[i, 1, 1] .* portfolio_specs[i, :ramp_up_limit] *
                    portfolio_specs[i, :capacity] *
                    portfolio_specs[i, :capacity_factor]
                )
            )

            # Match generation: ramp-down constraint
            @constraint(
                m,
                (
                    g[i, 1, 1] - gen_data[i] >=
                    (-1) *
                    c[i, 1, 1] *
                    portfolio_specs[i, :ramp_down_limit] *
                    portfolio_specs[i, :capacity] *
                    portfolio_specs[i, :capacity_factor]
                )
            )

            # If not VRE, match commitment + start-up - shut_down
            if portfolio_specs[i, :is_VRE] == 0
                @constraint(m, c[i, 1, 1] - commit_data[i] - su[i, 1, 1] + sd[i, 1, 1] == 0)
            end
        end
    end


    # Compute number of units started up and shut down per hour
    for i = 1:num_units
        # Exclude wind and solar
        if portfolio_specs[i, :is_VRE] == 0
            for k = 1:num_days
                # Intraday constraint
                if k > 1
                    @constraint(m, c[i, k, 1] - c[i, k-1, num_hours] - su[i, k, 1] + sd[i, k, 1] == 0)
                end

                for j = 2:num_hours
                    @constraint(m, c[i, k, j] - c[i, k, j-1] - su[i, k, j] + sd[i, k, j] == 0)
                end
            end
        end
    end

    # Limit total energy market participation per unit type each hour to the total capacity of
    #   all committed units of this type, with committed units subject to
    #   minimum and maximum power levels
    # wind and solar are excluded
    for i = 1:num_units
        if portfolio_specs[i, :is_VRE] == 0
            for k = 1:num_days
                for j = 1:num_hours
                    # Constrain total market participation to capacity * capacity_factor
                    @constraint(
                        m,
                        g[i, k, j] + r[i, k, j] + sr[i, k, j] + nsr[i, k, j] <= (
                            c[i, k, j] .* portfolio_specs[i, :capacity] .*
                            portfolio_specs[i, :capacity_factor]
                        )
                    )
                    @constraint(
                        m,
                        g[i, k, j] + r[i, k, j] + sr[i, k, j] + nsr[i, k, j] >= (
                            c[i, k, j] .* portfolio_specs[i, :capacity] .*
                            portfolio_specs[i, :min_PL]
                        )
                    )
                end
            end
        end
    end

    # Limit solar and wind generation and AS participation to their actual 
    #   hourly availability, if wind and/or solar exist in the system
    for i = 1:num_units
        if portfolio_specs[i, :is_VRE] == 1
            # Allocate representative days data based on fundamental unit type
            if occursin("wind", portfolio_specs[i, :unit_type])
                availability = wind_repdays
            elseif occursin("solar", portfolio_specs[i, :unit_type])
                availability = solar_repdays
            else
                continue
            end

            # Set constraint
            for k = 1:num_days
                for j = 1:num_hours
                    @constraint(
                        m,
                        g[i, k, j] + r[i, k, j] + sr[i, k, j] + nsr[i, k, j]
                          <= availability[j, k] * portfolio_specs[i, :esc_num_units] * portfolio_specs[i, :capacity]
                    )
                end
            end
        end
    end

    # Constrain AS participation to specified maximum levels
    for i = 1:num_units
        for k = 1:num_days
            for j = 1:num_hours
                # Frequency regulation
                @constraint(
                    m,
                    r[i, k, j] 
                      <= (c[i, k, j]
                          .* portfolio_specs[i, :capacity] 
                          .* portfolio_specs[i, :capacity_factor] 
                          .* portfolio_specs[i, :max_regulation]
                         )
                )

                # Spinning reserve
                @constraint(
                    m,
                    sr[i, k, j] 
                      <= (c[i, k, j]
                          .* portfolio_specs[i, :capacity] 
                          .* portfolio_specs[i, :capacity_factor] 
                          .* portfolio_specs[i, :max_spinning_reserve]
                         )
                )

                # Non-spinning reserve
                @constraint(
                    m,
                    nsr[i, k, j] 
                      <= (c[i, k, j]
                          .* portfolio_specs[i, :capacity] 
                          .* portfolio_specs[i, :capacity_factor] 
                          .* portfolio_specs[i, :max_nonspinning_reserve]
                         )
                )

            end
        end
    end


    # Force all conventional_nuclear units to be committed at all times
    for k = 1:num_days
        for j = 1:num_hours
            if convnuc_index != nothing
                @constraint(m, c[convnuc_index, k, j] == portfolio_specs[convnuc_index, :num_units])
            end
        end
    end

    # Ramping constraints
    for i = 1:num_units
        if portfolio_specs[i, :is_VRE] == 0
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
    end

    ENS_penalty = settings["constants"]["big_number"]
    ASNS_penalty = ENS_penalty * settings["dispatch"]["ASNS_penalty_ratio"]
    RNS_subpenalty = settings["dispatch"]["rns_subpenalty"]
    SNS_subpenalty = settings["dispatch"]["sns_subpenalty"]
    NSNS_subpenalty = settings["dispatch"]["nsns_subpenalty"]
    gamma_reg = settings["dispatch"]["gamma_reg"]
    gamma_spin = settings["dispatch"]["gamma_spin"]
    gamma_nspin = settings["dispatch"]["gamma_nspin"]

    # Rescale AS subpenalties
    total = RNS_subpenalty + SNS_subpenalty + NSNS_subpenalty
    RNS_subpenalty = RNS_subpenalty / total
    SNS_subpenalty = SNS_subpenalty / total
    NSNS_subpenalty = NSNS_subpenalty / total

    @objective(
        m,
        Min,
        sum(   # sum over d
            sum(   # sum over h
                sum(   # sum over i
                    # Generation price
                    g[i, k, j] .* (
                        portfolio_specs[i, :VOM] 
                        + portfolio_specs[i, :FC_per_MWh] 
                        - portfolio_specs[i, :carbon_tax_per_MWh] 
                        - portfolio_specs[i, :tax_credits_per_MWh]
                    )
                    # Ancillary services prices
                    + (r[i, k, j] .* gamma_reg + sr[i, k, j] .* gamma_spin + nsr[i, k, j] .* gamma_nspin) .* portfolio_specs[i, :VOM]
                    # Commitment/no-load costs
                    + c[i, k, j] .* portfolio_specs[i, :no_load_cost]
                    # Start-up costs
                    + su[i, k, j] .* portfolio_specs[i, :start_up_cost]
                    # Shut-down costs
                    + sd[i, k, j] .* portfolio_specs[i, :shut_down_cost]
                for i = 1:num_units)
                # Penalty for energy not served and ancillary services
                #   not served
                + ens[k, j] .* ENS_penalty + (RNS_subpenalty .* rns[k, j] .+ SNS_subpenalty .* sns[k, j] + NSNS_subpenalty .* nsns[k, j]) .* ASNS_penalty
            for j = 1:num_hours)
        for k = 1:num_days)
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
        r = nothing
        sr = nothing
        nsr = nothing
        c = nothing
        su = nothing
        sd = nothing
        ens = nothing
        rns = nothing
        sns = nothing
        nsns = nothing
        if status == "OPTIMAL"
            gen_qty = value.(model[:g])
            r = value.(model[:r])
            sr = value.(model[:sr])
            nsr = value.(model[:nsr])
            c = value.(model[:c])
            su = value.(model[:su])
            sd = value.(model[:sd])
            ens = value.(model[:ens])
            rns = value.(model[:rns])
            sns = value.(model[:sns])
            nsns = value.(model[:nsns])
        end

        returns = model, status, gen_qty, r, sr, nsr, c, su, sd, ens, rns, sns, nsns

    elseif model_type == "relaxed_integrality"
        optimize!(model)
        status = string(termination_status.(model))
        @debug "$model_type model status: $status"

        returns = model

    end

    return returns
end


function assemble_grc_results(y, gen_qty, r, sr, nsr, c, su, sd, portfolio_specs; starting_index=1)
    new_grc_results = set_up_grc_results_df()

    # Save generation and commitment results to a dataframe
    # g[i, k, j] and c[i, k, j] are arranged as
    #   [num_units, num_days, num_hours]
    for k = 1:size(c)[2]            # num_days
        for j = 1:size(c)[3]        # num_hours
            for i = 1:size(c)[1]    # num_units
                line = (
                    y = y,
                    d = k + starting_index - 1,
                    h = j,
                    unit_type = portfolio_specs[i, :unit_type],
                    gen = gen_qty[i, k, j],
                    reg = r[i, k, j],
                    spin = sr[i, k, j],
                    nspin = nsr[i, k, j],
                    commit = round(Int, c[i, k, j]),
                    su = round(Int, su[i, k, j]),
                    sd = round(Int, sd[i, k, j]),
                )
                push!(new_grc_results, line)
            end
        end
    end

    return new_grc_results
end


function reshape_shadow_prices(gen_shadow_prices, reg_shadow_prices, spin_shadow_prices, nspin_shadow_prices, y, settings; starting_index=1)
    price_df = set_up_prices_df()

    # Convert the (repdays, hours) table into a long (y, d, h, price) table
    for k = 1:size(gen_shadow_prices)[1]        # num_days
        for j = 1:size(gen_shadow_prices)[2]    # num_hours
            gen_price = min(settings["system"]["price_cap"], (-1) * gen_shadow_prices[k, j])
            reg_price = min(settings["system"]["AS_price_cap"], (-1) * reg_shadow_prices[k, j])
            spin_price = min(settings["system"]["AS_price_cap"], (-1) * spin_shadow_prices[k, j])
            nspin_price = min(settings["system"]["AS_price_cap"], (-1) * nspin_shadow_prices[k, j])

            # Add the price data to the table, indexed by y,d,h
            line = (
                y = y,
                d = k + starting_index - 1,
                h = j,
                lambda = gen_price,
                reg_rmp = reg_price,
                spin_rmp = spin_price,
                nspin_rmp = nspin_price
            )
            push!(price_df, line)
        end
    end

    return price_df
end


function propagate_all_results(all_grc_results, all_prices, current_pd, end_year)
    final_dispatched_year = maximum(all_grc_results[!, :y])

    final_year_grc =
        filter(:y => x -> x == final_dispatched_year, all_grc_results)
    final_year_prices =
        filter(:y => x -> x == final_dispatched_year, all_prices)

    for y = (final_dispatched_year + 1):(current_pd + end_year - 1)
        # Copy the final_year_grc results forward, updating the year
        next_year_grc = deepcopy(final_year_grc)
        next_year_grc[!, :y] .= y
        # Push row-by-row (avoids Julia's strict for-loop scoping)
        for i = 1:size(next_year_grc)[1]
            push!(all_grc_results, next_year_grc[i, :])
        end

        # Copy the final_year_prices results forward, updating the year
        next_year_prices = deepcopy(final_year_prices)
        next_year_prices[!, :y] .= y
        for i = 1:size(next_year_prices)[1]
            push!(all_prices, next_year_prices[i, :])
        end
    end

    return all_grc_results, all_prices

end


function run_annual_dispatch(
    settings,
    y,
    year_portfolio,
    peak_demand,
    ts_data,
    unit_specs;
    run_mode="forecast",
    downselection_mode="scenario_reduction"
)
    # Set up the appropriate scenario reduction parameters
    if run_mode == "current"
        num_days = Int(round(size(ts_data[:load_data])[1] / 24))
        num_hours = 24
    else
        num_days, num_hours = set_repdays_params(settings, ts_data)
    end

    # Scale the load data to the PD value for this year
    ts_data = scale_load(ts_data, peak_demand)

    # Set up representative days for load
    ts_data = set_up_load_repdays(downselection_mode, num_days, num_hours, ts_data)

    # Set up representative days for wind and solar
    ts_data = set_up_wind_solar_repdays(downselection_mode, num_days, num_hours, ts_data)

    # Scale the ancillary services to the PD value for this year
    ts_data = scale_AS_data(ts_data, peak_demand, settings["scenario"]["peak_demand"])

    # Set up representative days for AS
    ts_data = set_up_AS_repdays(downselection_mode, num_days, num_hours, ts_data)

    # If running in forecast mode, run the entire "year" at once
    if (run_mode == "forecast") && (settings["dispatch"]["num_repdays"] < 100)
        @debug "Setting up optimization model..."
        m, portfolio_specs =
            set_up_model(settings, num_days, num_hours, ts_data, year_portfolio, unit_specs; gen_data=nothing, commit_data=nothing)

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
        m, status, gen_qty, r, sr, nsr, c, su, sd, ens, rns, sns, nsns = solve_model(m, model_type = "integral")

        # Set up default return values if the dispatch year is infeasible
        new_grc_results = nothing
        new_prices = nothing
        summary_statistics = nothing

        # Save the generation and commitment results from the integral problem
        new_grc_results = assemble_grc_results(y, gen_qty, r, sr, nsr, c, su, sd, portfolio_specs)

        # Set up a relaxed-integrality version of this model, to allow
        #   retrieval of dual values for the mkt_equil constraint
        @debug "Solving relaxed-integrality problem to get shadow prices..."
        undo = relax_integrality(m_copy)

        # Solve the relaxed-integrality model to compute the shadow prices
        m_copy = solve_model(m_copy, model_type = "relaxed_integrality")
        new_prices = reshape_shadow_prices(
            shadow_price.(m_copy[:mkt_equil]),   # generation shadow prices
            shadow_price.(m_copy[:reg_mkt_equil]),  # regulation shadow prices
            shadow_price.(m_copy[:spin_mkt_equil]),  # spinning reserve shadow prices
            shadow_price.(m_copy[:nspin_mkt_equil]), # non-spinning reserve shadow prices
            y,
            settings,
        )

        @debug "Year $y dispatch run complete."

    else
        # If running in 'current' mode (i.e. a full year), separate the total
        #   year into subperiods to reduce memory requirements
        subpd = settings["dispatch"]["annual_dispatch_subperiod"] # in days
        starting_index = 1
        gen_data = nothing
        commit_data = nothing

        all_grc_results = set_up_grc_results_df()
        all_prices = set_up_prices_df()
        all_ens = nothing
        all_rns = nothing
        all_sns = nothing
        all_nsns = nothing

        while starting_index < 365
            # Take slices of the time-series data
            ending_index = min(starting_index + subpd - 1, size(ts_data[:load_repdays])[2])
            @info "  Simulating dispatch for days $starting_index - $ending_index"

            ts_data_slice = Dict(
                :load_repdays => ts_data[:load_repdays][:, starting_index:ending_index],
                :wind_repdays => ts_data[:wind_repdays][:, starting_index:ending_index],
                :solar_repdays => ts_data[:solar_repdays][:, starting_index:ending_index],
                :reg_repdays => ts_data[:reg_repdays][:, starting_index:ending_index],
                :spin_repdays => ts_data[:spin_repdays][:, starting_index:ending_index],
                :nspin_repdays => ts_data[:nspin_repdays][:, starting_index:ending_index],
            )

            num_slice_days = ending_index - starting_index + 1

            @debug "Setting up optimization model..."
            m, portfolio_specs =
                set_up_model(settings, num_slice_days, num_hours, ts_data_slice, year_portfolio, unit_specs; gen_data=gen_data, commit_data=commit_data)

            @debug "Optimization model set up."
            @debug string("Solving annual dispatch for subperiod beginning ", starting_index, "...")

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
            m, status, gen_qty, r, sr, nsr, c, su, sd, ens, rns, sns, nsns = solve_model(m, model_type = "integral")

            # Save last-hour commitment and generation data to set up
            #   next period's continuity constraints
            commit_data = deepcopy(c)[:, num_slice_days, num_hours]
            gen_data = deepcopy(gen_qty)[:, num_slice_days, num_hours]

            # Save the generation, commitment, and ancillary services results
            #    from the integral problem
            all_grc_results = vcat(all_grc_results, assemble_grc_results(y, gen_qty, r, sr, nsr, c, su, sd, portfolio_specs; starting_index=starting_index))

            if starting_index == 1
                # On the first subperiod, overwrite the default `nothing' 
                #   objects with the new Matrix objects
                all_ens = ens
                all_rns = rns
                all_sns = sns
                all_nsns = nsns
            else
                # For subsequent subperiods, concatenate new data to the end
                #   of the existing Matrix objects
                all_ens = vcat(all_ens, ens)
                all_rns = vcat(all_rns, rns)
                all_sns = vcat(all_sns, sns)
                all_nsns = vcat(all_nsns, nsns)
            end

            # Set up a relaxed-integrality version of this model, to allow
            #   retrieval of dual values for the mkt_equil constraint
            @debug "Solving relaxed-integrality problem to get shadow prices..."
            undo = relax_integrality(m_copy)

            # Solve the relaxed-integrality model to compute the shadow prices
            m_copy = solve_model(m_copy, model_type = "relaxed_integrality")
            price_results = reshape_shadow_prices(
                shadow_price.(m_copy[:mkt_equil]),   # generation shadow prices
                shadow_price.(m_copy[:reg_mkt_equil]),  # regulation shadow prices
                shadow_price.(m_copy[:spin_mkt_equil]),  # spinning reserve shadow prices
                shadow_price.(m_copy[:nspin_mkt_equil]), # non-spinning reserve shadow prices
                y,
                settings;
                starting_index=starting_index
            )

            all_prices = vcat(all_prices, price_results)

            @debug "Subperiod starting $starting_index dispatch run complete."

            # Empty the old optimization model
            empty!(m)
            empty!(m_copy)


            # Update the new starting_index for the next slice of data
            starting_index = ending_index + 1
        end

        # Rename data objects to their standard names
        new_grc_results = all_grc_results
        new_prices = all_prices
        ens = all_ens
        rns = all_rns
        sns = all_sns
        nsns = all_nsns

    end

    if run_mode == "current"
        summary_statistics = calculate_summary_statistics(new_grc_results, new_prices, ens, rns, sns, nsns)
    else
        summary_statistics = nothing
    end

    results = Dict(
        :new_grc_results => new_grc_results,
        :new_prices => new_prices,
        :summary_statistics => summary_statistics
    )

    # Set m and m_copy to nothing, to reduce memory utilization
    m = nothing
    m_copy = nothing

    return results

end


function calculate_summary_statistics(new_grc_results, new_prices, ens, rns, sns, nsns)
    econ_res = innerjoin(new_grc_results, new_prices, on = [:y, :d, :h])

    # Compute weighted average generation price
    transform!(econ_res, [:gen, :lambda] => ((gen, lambda) -> gen .* lambda) => :gen_tx)
    wa_gen_price = sum(econ_res[!, :gen_tx]) / sum(econ_res[!, :gen])

    # Compute weighted average regulation price
    transform!(econ_res, [:reg, :reg_rmp] => ((reg, reg_rmp) -> reg .* reg_rmp) => :reg_tx)
    wa_reg_price = sum(econ_res[!, :reg_tx]) / sum(econ_res[!, :reg])

    # Compute weighted average spinning reserve price
    transform!(econ_res, [:spin, :spin_rmp] => ((spin, spin_rmp) -> spin .* spin_rmp) => :spin_tx)
    wa_spin_price = sum(econ_res[!, :spin_tx]) / sum(econ_res[!, :spin])

    # Compute weighted average non-spinning reserve price
    transform!(econ_res, [:nspin, :nspin_rmp] => ((nspin, nspin_rmp) -> nspin .* nspin_rmp) => :nspin_tx)
    wa_nspin_price = sum(econ_res[!, :nspin_tx]) / sum(econ_res[!, :nspin])

    total_ens = sum(sum(ens[k, j] for k = 1:size(ens)[1]) for j = 1:size(ens)[2])
    total_rns = sum(sum(rns[k, j] for k = 1:size(rns)[1]) for j = 1:size(rns)[2])
    total_sns = sum(sum(sns[k, j] for k = 1:size(sns)[1]) for j = 1:size(sns)[2])
    total_nsns = sum(sum(nsns[k, j] for k = 1:size(nsns)[1]) for j = 1:size(nsns)[2])

    summary_statistics = Dict(
        :wa_gen_price => wa_gen_price,
        :wa_reg_price => wa_reg_price,
        :wa_spin_price => wa_spin_price,
        :wa_nspin_price => wa_nspin_price,
        :total_ens => total_ens,
        :total_rns => total_rns,
        :total_sns => total_sns,
        :total_nsns => total_nsns,
    )

    return summary_statistics
end


function save_summary_statistics(db, current_pd, ss)
    vals = (
        current_pd,
        ss[:wa_gen_price],
        ss[:wa_reg_price],
        ss[:wa_spin_price],
        ss[:wa_nspin_price],
        ss[:total_ens],
        ss[:total_rns],
        ss[:total_sns],
        ss[:total_nsns],
    )

    stmt = "INSERT INTO annual_dispatch_summary VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)"

    DBInterface.execute(db, stmt, vals)
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
    for i = (last_dispatch_year + 1):forecast_end_pd
        df = system_portfolios[last_dispatch_year]
        df[!, :y] .= i
        append!(all_year_portfolios, df)
    end

    return all_year_portfolios

end


function join_results_data_frames(
    settings,
    all_grc_results,
    all_prices,
    all_repdays,
    all_year_portfolios,
    unit_specs,
)
    # Join in price data to all_grc_results
    long_econ_results = innerjoin(all_grc_results, all_prices, on = [:y, :d, :h])

    # Incorporate repdays probability data
    if settings["dispatch"]["downselection"] == "exact"
        long_econ_results[!, :Probability] .= 1.0/365
    else
        long_econ_results = innerjoin(
            long_econ_results,
            select(all_repdays, [:index, :Probability, :y]),
            on = [:y, :d => :index],
        )
    end

    # Join in limited unit_specs data to compute some cash flows
    long_econ_results = innerjoin(
        long_econ_results,
        select(
            unit_specs,
            [:unit_type, :capacity, :VOM, :FC_per_MWh, :carbon_tax_per_MWh, :tax_credits_per_MWh],
        ),
        on = :unit_type,
    )

    # Join in unit number data to long_rev_results
    long_econ_results =
        innerjoin(long_econ_results, all_year_portfolios, on = [:y, :unit_type], makeunique=true)

    return long_econ_results

end


function compute_per_unit_cash_flows(long_econ_results)
    # Calculate generation
    transform!(
        long_econ_results,
        [:gen, :Probability, :num_units] =>
            ((gen, prob, num_units) -> gen .* prob .* 365 ./ num_units) =>
                :annualized_gen_per_unit,
    )

    # Calculate regulation service
    transform!(
        long_econ_results,
        [:reg, :Probability, :num_units] =>
            ((reg, prob, num_units) -> reg .* prob .* 365 ./ num_units) =>
                :annualized_reg_per_unit,
    )

    # Calculate spinning reserve service
    transform!(
        long_econ_results,
        [:spin, :Probability, :num_units] =>
            ((spin, prob, num_units) -> spin .* prob .* 365 ./ num_units) =>
                :annualized_spin_per_unit,
    )

    # Calculate non-spinning reserve service
    transform!(
        long_econ_results,
        [:nspin, :Probability, :num_units] =>
            ((nspin, prob, num_units) -> nspin .* prob .* 365 ./ num_units) =>
                :annualized_nspin_per_unit,
    )

    # Calculate generation revenue
    transform!(
        long_econ_results,
        [:gen, :lambda, :Probability, :num_units] =>
            ((gen, lambda, prob, num_units) -> (gen .* lambda .* prob .* 365 ./ num_units))
                => :annualized_gen_rev_per_unit,
    )

    # Calculate frequency regulation revenue
    transform!(
        long_econ_results,
        [:reg, :reg_rmp, :Probability, :num_units]
        => ((reg, reg_rmp, prob, num_units)
        -> (reg .* reg_rmp .* prob .* 365 ./ num_units))
        => :annualized_reg_rev_per_unit,
    )

    # Calculate spinning reserve revenue
    transform!(
        long_econ_results,
        [:spin, :spin_rmp, :Probability, :num_units]
        => ((spin, spin_rmp, prob, num_units)
        -> (spin .* spin_rmp .* prob .* 365 ./ num_units))
        => :annualized_spin_rev_per_unit,
    )

    # Calculate non-spinning reserve revenue
    transform!(
        long_econ_results,
        [:nspin, :nspin_rmp, :Probability, :num_units]
        => ((nspin, nspin_rmp, prob, num_units)
        -> (nspin .* nspin_rmp .* prob .* 365 ./ num_units))
        => :annualized_nspin_rev_per_unit,
    )

    # Calculate total reserves revenue
    transform!(
        long_econ_results,
        [:annualized_reg_rev_per_unit, :annualized_spin_rev_per_unit, :annualized_nspin_rev_per_unit]
        => ((ann_reg_rev, ann_spin_rev, ann_nspin_rev)
        -> (ann_reg_rev .+ ann_spin_rev .+ ann_nspin_rev))
        => :annualized_reserves_rev_per_unit,
    )

    # Calculate total revenue
    transform!(
        long_econ_results,
        [:annualized_gen_rev_per_unit, :annualized_reserves_rev_per_unit]
        => ((gen_rev, res_rev)
        -> (gen_rev .+ res_rev))
        => :annualized_rev_per_unit,
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
        [:gen, :carbon_tax_per_MWh, :Probability, :num_units] =>
            (
                (gen, adj, prob, num_units) ->
                    gen .* adj .* prob .* 365 ./ num_units
            ) => :annualized_carbon_tax_per_unit,
    )

    # Calculate tax credit adjustment
    transform!(
        long_econ_results,
        [:gen, :tax_credits_per_MWh, :Probability, :num_units] =>
            (
                (gen, tc, prob, num_units) ->
                    gen .* tc .* prob .* 365 ./ num_units
            ) => :annualized_tax_credits_per_unit,
    )


    for col in eachcol(long_econ_results)
        replace!(col, Inf => 0)
        replace!(col, NaN => 0)
    end

    return long_econ_results

end


function summarize_dispatch_results(settings, unit_specs, long_econ_results)
    dispatch_results = deepcopy(long_econ_results)

    dispatch_results = combine(
        groupby(dispatch_results, [:y, :unit_type]),
        [
            :annualized_gen_per_unit,
            :annualized_reg_per_unit,
            :annualized_spin_per_unit,
            :annualized_nspin_per_unit,
            :annualized_gen_rev_per_unit,
            :annualized_reg_rev_per_unit,
            :annualized_spin_rev_per_unit,
            :annualized_nspin_rev_per_unit,
            :annualized_reserves_rev_per_unit,
            :annualized_rev_per_unit,
            :annualized_VOM_per_unit,
            :annualized_FC_per_unit,
            :annualized_carbon_tax_per_unit,
            :annualized_tax_credits_per_unit,
        ] .=> sum,
    )

    rename!(
        dispatch_results,
        :annualized_gen_per_unit_sum => :generation,
        :annualized_reg_per_unit_sum => :regulation,
        :annualized_spin_per_unit_sum => :spinning_reserve,
        :annualized_nspin_per_unit_sum => :nonspinning_reserve,
        :annualized_gen_rev_per_unit_sum => :gen_rev,
        :annualized_reg_rev_per_unit_sum => :reg_rev,
        :annualized_spin_rev_per_unit_sum => :spin_rev,
        :annualized_nspin_rev_per_unit_sum => :nspin_rev,
        :annualized_reserves_rev_per_unit_sum => :reserves_rev,
        :annualized_rev_per_unit_sum => :revenue,
        :annualized_VOM_per_unit_sum => :VOM,
        :annualized_FC_per_unit_sum => :fuel_cost,
        :annualized_carbon_tax_per_unit_sum => :carbon_tax,
        :annualized_tax_credits_per_unit_sum => :tax_credits,
    )

    # Pivot in FOM data
    FOM_data = select(unit_specs, [:unit_type, :FOM, :capacity])
    dispatch_results = leftjoin(dispatch_results, FOM_data, on = :unit_type)
    transform!(
        dispatch_results,
        [:FOM, :capacity] =>
            ((FOM, cap) -> FOM .* cap .* settings["constants"]["MW2kW"]) =>
                :FOM,
    )
    select!(dispatch_results, Not(:capacity))

    # Put the data into a long format for easier filtering
    dispatch_results = stack(
        dispatch_results,
        [:generation, :regulation, :spinning_reserve, :nonspinning_reserve, :gen_rev, :reg_rev, :spin_rev, :nspin_rev, :reserves_rev, :revenue, :VOM, :fuel_cost, :FOM, :carbon_tax, :tax_credits],
    )
    rename!(dispatch_results, :variable => :dispatch_result, :value => :qty)

    return dispatch_results
end


function postprocess_results(
    settings,
    all_grc_results,
    all_prices,
    all_repdays,
    all_year_system_portfolios,
    unit_specs,
    current_pd,
    fc_pd,
)
    # Propagate the results dataframes out to the end of the projection horizon
    # Assume no change after the last modeled year
    all_grc_results, all_prices =
        propagate_all_results(all_grc_results, all_prices, current_pd, fc_pd)

    # Get a single long dataframe with unit numbers by type for each year
    all_year_portfolios = combine_and_extend_year_portfolios(
        all_year_system_portfolios,
        current_pd + fc_pd,
    )

    # Join in relevant data to create a single master dataframe suitable for
    #   a variety of column-wise operations
    long_econ_results = join_results_data_frames(
        settings,
        all_grc_results,
        all_prices,
        all_repdays,
        all_year_portfolios,
        unit_specs,
    )

    # Compute cash flows on a per-unit basis: revenue, VOM, fuel cost, and
    #   policy adjustment
    long_econ_results = compute_per_unit_cash_flows(long_econ_results)

    dispatch_results =
        summarize_dispatch_results(settings, unit_specs, long_econ_results)

    return long_econ_results, dispatch_results
end


function finalize_annual_dispatch_results(db, current_pd, long_econ_results, dispatch_results)
    #save_annual_dispatch_summary(db, current_pd, long_econ_results)
    #save_annual_dispatch_hourly_results(db, current_pd, long_econ_results)
    save_annual_dispatch_unit_summary(db, current_pd, dispatch_results)
    #save_annual_dispatch_hourly_unit_results(db, current_pd, long_econ_results)
end


function save_annual_dispatch_unit_summary(db, current_pd, dispatch_results)
    pivot = unstack(dispatch_results, :unit_type, :dispatch_result, :qty)

    pivot[!, :period] .= current_pd

    # Put the columns in the same order as the DB table
    col_order = DBInterface.execute(
        db,
        "SELECT name FROM PRAGMA_TABLE_INFO('annual_dispatch_unit_summary')",
    ) |> DataFrame
    col_order = collect(col_order[!, "name"])
    pivot = select(pivot, col_order)

    # Set up the Sql "(?, ?, ..., ?)" string of correct length
    fill_tuple = string("(", repeat("?, ", size(col_order)[1]-1), "?)")

    # Add each row to the database table
    for row in Tuple.(eachrow(pivot))
        DBInterface.execute(
            db,
            string("INSERT INTO annual_dispatch_unit_summary VALUES $fill_tuple"),
            row,
        )
    end
end

function save_annual_dispatch_hourly_unit_results(db, current_pd, long_econ_results)
    # Get columns of interest
    results = deepcopy(long_econ_results[:, [:y, :d, :h, :unit_type, :gen, :reg, :spin, :nspin]])

    # Rename the columns to match the DB standard
    rename!(
        results,
        :y => :period,
        :d => :day,
        :h => :hour,
        :gen => :generation,
        :reg => :regulation,
        :spin => :spinning_reserve,
        :nspin => :nonspinning_reserve,
    )

    stmt = DBInterface.prepare(db, """INSERT INTO annual_disp_hr_unit_results VALUES (:period, :day, :hour, :unit_type, :generation, :regulation, :spinning_reserve, :nonspinning_reserve)""")

    DBInterface.executemany(
        stmt,
        (period = results[!, :period],
         day = results[!, :day],
         hour = results[!, :hour],
         unit_type = results[!, :unit_type],
         generation = results[!, :generation],
         regulation = results[!, :regulation],
         spinning_reserve = results[!, :spinning_reserve],
         nonspinning_reserve = results[!, :nonspinning_reserve],
        )
    )

end


end
