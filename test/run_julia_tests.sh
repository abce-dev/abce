#!/bin/bash
test_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
abce_dir="$( dirname $test_dir )"
cat "$abce_dir/env/Project.toml"
echo $abce_dir
# Use testSysimage.so if available
if [[ -f "$test_dir/testSysimage.so" ]]; then
    julia --project="$abce_dir/env" -J "$test_dir/testSysimage.so" "$test_dir/test_julia.jl"
else
    julia --project="$abce_dir/env" "$test_dir/test_julia.jl"
fi
