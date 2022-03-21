using Logging

@info "Loading packages..."

using CSV, DataFrames, JuMP, GLPK, XLSX, Logging, CPLEX, BilevelJuMP

@info "Initializing data..."

# Peak demand values: number of entries determines number of years considered
PD = [80000, 81000, 83000]
repday_ids = [10, 90, 180, 292, 355]

# Load the time-series demand and VRE data into dataframes
ts_data = CSV.read("./inputs/ALEAF_inputs/timeseries_load_hourly.csv", DataFrame)
wind_data = CSV.read("./inputs/ALEAF_inputs/timeseries_wind_hourly.csv", DataFrame)
solar_data = CSV.read("./inputs/ALEAF_inputs/timeseries_pv_hourly.csv", DataFrame)

# Load the unit specs data into a dataframe
unit_specs_data, unit_spec_labels = XLSX.readtable("./inputs/ALEAF_inputs/ALEAF_Master_LC_GEP_dispatch.xlsx", "Gen Technology", header=true)
unit_specs = DataFrame()
for i in 1:size(unit_spec_labels)[1]
    unit_specs[!, unit_spec_labels[i]] = unit_specs_data[i]
end

# Set up system portfolio data
portfolio_data, portfolio_headers = XLSX.readtable("./inputs/ALEAF_inputs/ALEAF_ERCOT_5units.xlsx", "gen", header=true)
portfolio = DataFrame()

# Add a CF column to unit_specs
unit_specs[!, :CF] .= 1.0
# non-VRE uses CAPCRED as CF
for i = 1:size(unit_specs)[1]
    if unit_specs[i, :VRE_Flag] == false
        unit_specs[i, :CF] = unit_specs[i, :CAPCRED]
    end
end

# Set up some convenient parameter names
num_units = size(unit_specs)[1]
num_days = size(repday_ids)[1]
num_hours = 24
@info "Data initialized."

# Iterate over all years
for y = 1:size(PD)[1]
    @info string("Preparing period ", y)

    # Scale the load data to the PD value for this year
    ts_data[!, :Load] = ts_data[!, :LoadShape] * PD[y]

    # Set up representative days for load
    load_repdays = DataFrame()
    for day in repday_ids
        load_repdays[!, Symbol(day)] = ts_data[24*day+1:24*(day+1), :Load]
    end

    # Assign the "current" portfolio
    # TODO: update yearly with new/retired units
    # Retrieving from fixed columns because I had problems filtering on
    #   Symbol-type headers.........
    portfolio[!, :unit_type] = portfolio_data[7]
    portfolio[!, :num_units] = portfolio_data[9]

    # Scale the wind and solar data according to the total installed capacities
    wind_data[!, :Wind] = wind_data[!, :WindShape] * filter(:unit_type => x -> x == "Wind", portfolio)[1, :num_units] * filter(:UNIT_TYPE => x -> x == "Wind", unit_specs)[1, :CAP]
    solar_data[!, :Solar] = solar_data[!, :SolarShape] * filter(:unit_type => x -> x == "Solar", portfolio)[1, :num_units] * filter(:UNIT_TYPE => x -> x == "Solar", unit_specs)[1, :CAP]

    # Set up representative days for wind and solar
    wind_repdays = DataFrame()
    solar_repdays = DataFrame()
    for day in repday_ids
        wind_repdays[!, Symbol(day)] = wind_data[24*day+1:24*(day+1), :Wind]
        solar_repdays[!, Symbol(day)] = solar_data[24*day+1:24*(day+1), :Solar]
    end

    @info "Setting up optimization model..."

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
                @constraint(m, c[i, k, j] <= portfolio[i, :num_units])
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

    @info "Optimization model set up."
    @info string("Solving repday dispatch for year ", y, "...")

    optimize!(m)
    status = termination_status.(m)
    @info "Status: ", status
    gen_qty = value.(m[:g])
    c = value.(m[:c])

    all_g_results = Dict()
    all_c_results = Dict()
    for i = 1:num_days
        g_results = DataFrame(gen_qty[:, i, :])
        c_results = DataFrame(c[:, i, :])
        insertcols!(g_results, 1, :unit_type => unit_specs[!, :UNIT_TYPE])
        insertcols!(c_results, 1, :unit_type => unit_specs[!, :UNIT_TYPE])
        all_g_results[i] = g_results
        all_c_results[i] = c_results
    end

    g_results = all_g_results[1]
    c_results = all_c_results[1]
    for i = 2:num_days
        append!(g_results, all_g_results[i])
        append!(c_results, all_c_results[i])
    end

    gfile = string("./gen_results", y, ".csv")
    cfile = string("./com_results", y, ".csv")
    CSV.write(gfile, g_results)
    CSV.write(cfile, c_results)

    @info "Results written to $gfile and $cfile."
end
