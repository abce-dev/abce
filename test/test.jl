module sTest

using Logging, DataFrames, YAML

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
    extra_info = ""

    if typeof(value) in [DataType, UnionAll]
    # If `value` is a DataType, perform type-checking for type or supertype
        if !(typeof(thing) == value) && !(thing <: value)
            pass = false
            result_msg = "fail"
            extra_info = string("/n", thing, " is not ", value)
        end
    elseif isa(thing, DataFrame) && isa(value, DataFrame)
    # If both objects are DataFrames, use DataFrame-specific equivalence check
        # Check for same columns
        if issetequal(names(thing), names(value))
            # If all columns match, reorder columns in thing to follow value
            thing = select(thing, names(value))

            # Sort both dataframes lexicographically, using defaults
            sort!(thing)
            sort!(value)

            # Check whether the dataframes are now identical
            if !isequal(thing, value)
                pass = false
                result_msg = "fail"
                extra_info = string("\n", thing, "\n is not equal to:\n", value)
            end
        end
    else
    # If `value` is not a DataType, assume a value comparison is desired
        if thing != value
            pass = false
            result_msg = "fail"
            extra_info = string("\n", thing, " != ", value)
        end
    end

    msg = string(test_header_msg, result_msg, extra_info)

    if pass
        @info msg
    else
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


end
