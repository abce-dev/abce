module ABCEfunctions

using DataFrames, CSV
export load_unit_type_data, get_demand_data_file, load_demand_data, set_forecast_period, forecast_demand, create_unit_FS_dfs


function load_unit_type_data(unit_file)
    try
        df = CSV.read(unit_file, DataFrame)
        num_types = size(df)[1]
        return df, num_types
    catch e
        println(string("Could not load data from the unit operational data file", unit_file))
        println(e)
        exit()
    end
end


function get_demand_data_file()
    return "./data/default_demand.csv"
end


function load_demand_data(demand_file)
    try
        df = CSV.read(demand_file, DataFrame)[!, :demand]
        return df
    catch e
        println(string("Could not load data from the demand file ", demand_file))
        println(e)
        exit()
    end
end


function set_forecast_period(df)
    return maximum(df[!, :d_x] + df[!, :unit_life])
end


function forecast_demand(available_demand, fc_pd)
    fill_demand = last(available_demand) * ones(fc_pd - size(available_demand)[1])
    return vcat(available_demand, fill_demand)
end


function create_unit_FS_dfs(unit_data, horizon)
    try
        fs_dict = Dict{String, DataFrame}()
        for i = 1:size(unit_data)[1]
            new_df = DataFrame(year = 1:horizon, xtr_exp = zeros(horizon), gen = zeros(horizon), remaining_debt_principal = zeros(horizon), debt_payment = zeros(horizon), interest_due = zeros(horizon), depreciation = zeros(horizon))
            name = unit_data[i, :name]
            fs_dict[name] = new_df
        end
        return fs_dict
    catch e
        println("Could not create the dataframes for the unit financial statements:")
        println(e)
        exit()
    end
end




end
