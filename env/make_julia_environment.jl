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


using Pkg, Logging, DelimitedFiles
Pkg.instantiate()

Pkg.add("ArgParse")
Pkg.build("ArgParse")

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--reqs_file"
    help = "absolute path to the julia reqirements file"
    required = false
    default = joinpath(ENV["ABCE_ENV"], "julia_requirements.csv")
end

CL_args = parse_args(s)

# Initialize a dictionary to track any problems arising in the process
problems = Dict()

# Delete preexisting .toml files to avoid contamination
files_to_delete = ["Manifest.toml", "Project.toml"]
for dfile in files_to_delete
    full_file = joinpath(ENV["ABCE_ENV"], dfile)
    if isfile(full_file)
        rm(full_file)
        @info "Removed file $full_file"
   end
end

# Activate local environment
Pkg.activate(ENV["ABCE_ENV"])

# Get the list of packages to set up
julia_pkg_list = readdlm(CL_args["reqs_file"], String)

# Add all Julia packages to the environment, and ensure they are all built
for pkg in julia_pkg_list
    try
        @info "Adding $pkg"
        Pkg.add(pkg)
        Pkg.build(pkg)
    catch e
        msg = string("A problem occurred while trying to install package ",
                     "$pkg. Installation of this package will be skipped, ",
                     "which may cause issues at ABCE runtime.")
        @warn msg
        problems[pkg] = e
    end
end

# Ensure all non-default Python packages are installed via Conda.jl
@info "Ensuring Conda packages are installed..."
Pkg.build("Conda")
using Conda

req_file = joinpath(@__DIR__, "requirements.txt")
open(req_file, "r") do filehandle
    conda_list = readlines(filehandle)
    for cpkg in conda_list
        @info "Adding $cpkg"
        try
            Conda.add(cpkg)
        catch e
            msg = string("A problem occurred while trying to install package ",
                         "$cpkg. Installation of this package will be ",
                         "skipped, which may cause issues at ABCE runtime.")
            @warn msg
            problems[cpkg] = e
        end
    end
end

# Inform the user if there were any problems
if length(problems) == 0
    @info "All packages installed and built successfully."
else
    println("\n\n")
    msg = string("All Julia packages and Python libraries loaded, with the ",
                 "following exceptions:")
    @warn msg

    for key in keys(problems)
        @warn string(key, ": ", problems[key])
    end
end

