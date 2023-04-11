include("../src/dispatch.jl")

using Logging, DataFrames, YAML, .Dispatch

function test(thing, value)
    # Gather some information about the function calling test, for reporting
    st = StackTraces.stacktrace()
    outer_func_line, func_name, func_line = get_stacktrace_entries(st)

    # Set up the variables for reporting test results
    pass = true
    test_header_msg = string(
                          "Invoking test ", func_name,
                          " (line ", outer_func_line, 
                          "):\n  Sub-test at line ", func_line, 
                          "\n  Result: "
                      )
    result_msg = "pass"

    if typeof(value) in [DataType, UnionAll]
    # If `value` is a DataType, perform type-checking for type or supertype
        if !(typeof(thing) == value) && !(thing <: value)
            pass = false
            result_msg = string(thing, " is not ", value)
        end
    else
    # If `value` is not a DataType, assume a value comparison is desired
        if thing != value
            pass = false
            result_msg = string(thing, " != ", value)
        end
    end

    msg = string(test_header_msg, result_msg)

    if pass
        global num_passed += 1
        @info msg
    else
        global num_failed += 1
        @warn msg
    end
end


function get_stacktrace_entries(st)
    # Determine the location of the test invocation in the top-level scope
    outer_func_line = st[3].line

    # Determine the test being called and the location of the sub-call
    #   to `test`
    func_name = st[2].func
    func_line = st[2].line

    return outer_func_line, func_name, func_line
end


function report_results()
    # Report summary of results to user
    println("\n\n")
    @info "==============SUMMARY=============="
    @info "Total number of tests passed: $num_passed"
    @info "Total number of tests failed: $num_failed"
end


###########################################
# User-defined tests here                 #
###########################################

function test_reshape_shadow_prices(shadow_prices, check_reshaped_shadow_prices, y, settings)
    reshaped_shadow_prices = Dispatch.reshape_shadow_prices(
                                 shadow_prices,
                                 y,
                                 settings
                             )

    test(reshaped_shadow_prices, check_reshaped_shadow_prices)

    return reshaped_shadow_prices
end


settings = YAML.load_file("../settings.yml")
test_shadow_prices = [-0.0 -0.0 -1.0 -5.0
                      -0.0 -2.0 -9.1 -10000.0
                      -1.5 -3.3 -0.0 -0.0]

check = [1 1 1 0.0
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

check = DataFrame(check, [:y, :d, :h, :price])
check[!, :y] = convert.(Int64, check[:, :y])
check[!, :d] = convert.(Int64, check[:, :d])
check[!, :h] = convert.(Int64, check[:, :h])

###########################################
# Test runner, with list of all tests     #
###########################################

function run_tests()
    # Set up global counters for passed and failed tests
    global num_passed = 0
    global num_failed = 0

    # Put list of tests here
    test_reshape_shadow_prices(test_shadow_prices, check, 1, settings)
    test_reshape_shadow_prices(test_shadow_prices, check, 2, settings)

    report_results()
end


run_tests()

