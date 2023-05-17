#!/bin/bash

##########################################################################
## Copyright 2023 Argonne National Laboratory
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


test_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
abce_dir="$( dirname $test_dir )"

echo "Running julia tests..."

# Use testSysimage.so if available
if [[ -f "$test_dir/testSysimage.so" ]]; then
    julia --project="$ABCE_ENV" -J "$test_dir/testSysimage.so" "$test_dir/test_julia.jl"
    rc=$?
else
    julia --project="$ABCE_ENV" "$test_dir/test_julia.jl"
    rc=$?
fi

echo "Tests finished."

exit $rc
