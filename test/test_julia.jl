include("../src/dispatch.jl")

using Logging, DataFrames, .Dispatch

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

function test_set_up_results_dfs()
    all_prices, all_gc_results = Dispatch.set_up_results_dfs()

    test_all_prices = DataFrame(
                          y = Int[],
                          g = Int[],
                          h = Int[],
                          price = Float64[]
                      )

    test(all_prices, DataFrame)
    test(all_gc_results, DataFrame)
    test(names(all_prices), names(test_all_prices))

end


function test_reshape_shadow_price(shadow_prices, check_reshaped_shadow_prices)
    reshaped_shadow_prices = Dispatch.reshape_shadow_price(
                                 shadow_prices,
                                 y,
                                 num_days,
                                 num_hours
                             )

    test(all_prices, check_reshaped_shadow_prices)

    return all_prices
end



###########################################
# Test runner, with list of all tests     #
###########################################

function run_tests()
    # Set up global counters for passed and failed tests
    global num_passed = 0
    global num_failed = 0

    # Put list of tests here
    test_set_up_results_dfs()

    report_results()
end


run_tests()

