module sTest

using Logging, DataFrames, CSV

function test(thing, value)
    # Set up reports resulting information, including defaults
    result = set_up_result()

    # For any of the following applicable test cases, perform an appropriate test
    # If `value` is a DataType, perform type-checking for type or supertype
    if typeof(value) in [DataType, UnionAll]
        result = test_check_type(thing, value, result)

    # If both objects are DataFrames, use DataFrame-specific equivalence check
    elseif isa(thing, DataFrame) && isa(value, DataFrame)
        result = test_dataframes_identical(thing, value, result)

    # If `value` is not a DataType, assume a value comparison is desired
    else
        result = test_equal(thing, value, result)

    end

    # Report the results of the test to the user
    report_test_result(result)
end


function set_up_result()
    # Gather some information about the function calling test, for reporting
    st = StackTraces.stacktrace()
    outer_func_line, func_name, func_line = get_stacktrace_entries(st)

    # Set up the variables for reporting test results
    result = Dict()
    result["pass"] = true
    result["test_header_msg"] = string(
                                    "Invoking test ", func_name,
                                    " (line ", outer_func_line, 
                                    "):\n  Sub-test at line ", func_line, 
                                    "\n  Result: "
                      )
    result["result_msg"] = "pass"
    result["extra_info"] = ""

    return result
end


function get_stacktrace_entries(st)
    # Determine the location of the test invocation in the top-level scope
    outer_func_line = st[4].line

    # Determine the test being called and the location of the sub-call
    #   to `test`
    func_name = st[3].func
    func_line = st[3].line

    return outer_func_line, func_name, func_line
end


function test_compare_types(thing, value, result)
    if !(typeof(thing) == value) && !(thing <: value)
        result["pass"] = false
        result["result_msg"] = "fail"
        result["extra_info"] = string("/n", thing, " is not ", value)
    end

    return result
end


function test_dataframes_identical(thing, value, result)
    # Check for same columns
    if issetequal(names(thing), names(value))
        # If all columns match, reorder columns in thing to follow value
        thing = select(thing, names(value))

        # Sort both dataframes lexicographically, using defaults
        sort!(thing)
        sort!(value)

        for col_name in names(value)
            if eltype(thing[!, col_name]) <: Union{Integer, Real}
                thing[!, col_name] = round.(thing[!, col_name], digits=5)
            end

            if eltype(value[!, col_name]) <: Union{Integer, Real}
                value[!, col_name] = round.(value[!, col_name], digits=5)
            end
        end

        # Check whether the dataframes are now identical
        if !isequal(thing, value)
            result["pass"] = false
            result["result_msg"] = "fail"
            result["extra_info"] = string("\n", thing, "\n is not equal to:\n", value)
        end
    end

    return result
end


function test_equal(thing, value, result)
    if thing != value
        result["pass"] = false
        result["result_msg"] = "fail"
        result["extra_info"] = string("\n", thing, " != ", value)
    end

    return result
end


function report_test_result(result)
    msg = string(result["test_header_msg"], result["result_msg"], result["extra_info"])

    if result["pass"]
        @info msg
    else
        @warn msg
    end

end

end
