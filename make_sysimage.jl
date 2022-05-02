# Package and function precompiler for unit_choice.jl
#
# Loading and compilation can be very slow for some Julia packages
#   (especially JuMP, CSV, and DataFrames), and some package functions require
#   lengthy first-use runtime compilation.
# This script uses the Julia PackageCompiler tool to create a Sysimage file,
#   which contains precompiled binaries for all of the packages used in
#   unit_choice.jl.
# The script first loads default packages into the current environment, and then
#   runs the current version of unit_choice.jl. Due to the command-line
#   flag `--trace-compile=precompile.jl`, it pipes information about runtime
#   function compilation to precompile.jl. Finally, the sysimage is created
#   by specifying the desired packages, and passing in the function precompile
#   data.
#
# After a fresh regeneration of the abceSysimage.so file, running unit_choice.jl
#   takes 3-4 seconds total. You may need to rerun this file occasionally as
#   additional function and package additions slow down execution time.

##########################################################################
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

using PackageCompiler, Pkg

# Load all default Julia packages into the current environment
Pkg.activate

# Run the current version of unit_choice.jl, outputting function-compilation
#   records to `./precompile.jl`
run(`julia --trace-compile=precompile.jl agent_choice.jl --settings_file=./settings.yml --current_pd=1 --agent_id=201`)

# Create `abceSysimage.so` using the specified packages and the newly
#   generated `precompile.jl` file.
create_sysimage([:CPLEX, :CSV, :DataFrames, :GLPK, :JuMP, :SQLite, :XLSX, :YAML]; sysimage_path="abceSysimage.so", precompile_statements_file="./precompile.jl")
