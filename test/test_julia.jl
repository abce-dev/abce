include("./test.jl")
include("../src/dispatch.jl")

using .sTest, .Dispatch, YAML, DataFrames, CSV

function test_reshape_shadow_prices(shadow_prices, check_reshaped_shadow_prices, y, settings)
    reshaped_shadow_prices = Dispatch.reshape_shadow_prices(
                                 shadow_prices,
                                 y,
                                 settings
                             )

    sTest.test(reshaped_shadow_prices, check_reshaped_shadow_prices)
end


function test_propagate_all_results(all_gc_results, all_prices, current_pd, end_year, test_gc_results, test_prices)
    propagated_gc_results, propagated_prices = Dispatch.propagate_all_results(all_gc_results, all_prices, current_pd, end_year)

    sTest.test(propagated_gc_results, test_gc_results)
    sTest.test(propagated_prices, test_prices)
end


function test_combine_and_extend_year_portfolios(dict_of_portfolios, comparison_portfolios_df, forecast_end_pd)
    combined_portfolios_df = Dispatch.combine_and_extend_year_portfolios(dict_of_portfolios, forecast_end_pd)

    sTest.test(combined_portfolios_df, comparison_portfolios_df)

end




###########################################
# Loading test data                       #
###########################################

input_shadow_prices = Matrix(
                          CSV.read(
                              "./test_data/raw_shadow_prices.csv",
                              DataFrame,
                              header=false
                          )
                      )

all_prices = CSV.read("./test_data/reshaped_prices.csv", DataFrame)
all_gc_results = CSV.read("./test_data/all_gc_results.csv", DataFrame)
prop_prices = CSV.read("./test_data/prop_prices.csv", DataFrame)
prop_gc_results = CSV.read("./test_data/prop_gc_results.csv", DataFrame)

# Portfolio data
y1_system_portfolio = CSV.read("./test_data/y1_system_portfolio.csv", DataFrame)
y2_system_portfolio = CSV.read("./test_data/y2_system_portfolio.csv", DataFrame)
year_portfolios = Dict(1 => y1_system_portfolio, 2 => y2_system_portfolio)
extended_system_portfolios = CSV.read("./test_data/extended_system_portfolio.csv", DataFrame)

###########################################
# Test runner, with list of all tests     #
###########################################

function run_tests()
    settings = YAML.load_file("../settings.yml")

    # Put list of tests here
    test_reshape_shadow_prices(
        input_shadow_prices,
        all_prices,
        1,                 # current year
        settings
    )

    test_propagate_all_results(
        all_gc_results,
        all_prices,
        1,                 # current year
        3,                 # end year
        prop_gc_results,   # static test data
        prop_prices        # static test data
    )

    test_combine_and_extend_year_portfolios(
        Dict(1 => y1_system_portfolio, 2 => y2_system_portfolio),
        extended_system_portfolios,
        3
    )

end


run_tests()

