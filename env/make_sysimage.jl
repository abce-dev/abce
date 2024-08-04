# Package and function precompiler for agent_choice.jl
##########################################################################
# Copyright 2023 Argonne National Laboratory
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

using PackageCompiler, Pkg, Requires

abce_abs_path = ENV["ABCE_DIR"]
env_path = joinpath(abce_abs_path, "env")
script_path = joinpath(abce_abs_path, "src", "agent_choice.jl")
settings_file= joinpath(abce_abs_path, "settings.yml")
inputs_path = joinpath(abce_abs_path, "inputs")

# Activate the current environment
Pkg.activate(env_path)

# Run agent_choice.jl as a standalone script, outputting function-compilation
#   records to `./precompile.jl`
run(`julia --project=$env_path --trace-compile=precompile.jl $script_path --current_pd=0 --agent_id=201 --verbosity=2 --settings_file=$settings_file --inputs_path=$inputs_path --abce_abs_path=$abce_abs_path`)

pkg_list = [
    :ArgParse,
    :Cbc,
    :CSV,
    :DataFrames,
    :GLPK,
    :HiGHS,
    :Logging,
    :JuMP,
    :SQLite,
    :Tables,
    :XLSX,
    :YAML
]

# If CPLEX is available, add it to the package list for precompilation
try
    using CPLEX
    push!(pkg_list, :CPLEX)
catch
    println("skipping CPLEX")
end

# Create `abceSysimage.so` using the specified packages and the newly
#   generated `precompile.jl` file.
create_sysimage(
    pkg_list;
    sysimage_path="abceSysimage.so",
    precompile_statements_file="./precompile.jl"
)

