include("./test.jl")
include("../src/dispatch.jl")

using .sTest, .Dispatch, YAML, DataFrames

function test_reshape_shadow_prices(shadow_prices, check_reshaped_shadow_prices, y, settings)
    reshaped_shadow_prices = Dispatch.reshape_shadow_prices(
                                 shadow_prices,
                                 y,
                                 settings
                             )

    sTest.test(reshaped_shadow_prices, check_reshaped_shadow_prices)

    return reshaped_shadow_prices
end


function test_propagate_all_results(all_gc_results, all_prices, current_pd, end_year, test_gc_results, test_prices)
    propagated_gc_results, propagated_prices = Dispatch.propagate_all_results(all_gc_results, all_prices, current_pd, end_year)

    sTest.test(propagated_gc_results, test_gc_results)
    sTest.test(propagated_prices, test_prices)
end



settings = YAML.load_file("../settings.yml")
test_shadow_prices = [-0.0 -0.0 -1.0 -5.0
                      -0.0 -2.0 -9.1 -10000.0
                      -1.5 -3.3 -0.0 -0.0]

check_prices = [1 1 1 0.0
                1 1 2 0.0
                1 1 3 1.0
                1 1 4 5.0
                1 2 1 0.0
                1 2 2 2.0
                1 2 3 9.1
                1 2 4 9001
                1 3 1 1.5
                1 3 2 3.3
                1 3 3 0.0
                1 3 4 0.0]
check_prices = DataFrame(check_prices, [:y, :d, :h, :price])

all_gc_results = [1 1 1 "abc" 1.0 1
                  1 1 1 "def" 100 1
                  1 1 2 "abc" 2.0 1
                  1 1 2 "def" 112 2
                  1 2 1 "abc" 75.1 1
                  1 2 1 "def" 0.0 0
                  1 2 2 "abc" 0.1 1
                  1 2 2 "def" 10000 256]
all_gc_results = DataFrame(all_gc_results, [:y, :d, :h, :unit_type, :gen, :commit])

y2_gc_results = copy(all_gc_results)
y2_gc_results[!, :y] .= 2

c2 = copy(check_prices)
c2[!, :y] .= 2

check_prop_gc_results = vcat(all_gc_results, y2_gc_results)
check_prop_prices = vcat(check_prices, c2)


###########################################
# Test runner, with list of all tests     #
###########################################

function run_tests()
    # Put list of tests here
    all_prices = test_reshape_shadow_prices(test_shadow_prices, check_prices, 1, settings)

    test_propagate_all_results(all_gc_results, all_prices, 1, 2, check_prop_gc_results, check_prop_prices)
end


run_tests()

