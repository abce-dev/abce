using Pkg, Logging, ArgParse

julia_pkg_list = ["ArgParse",
                  "CSV",
                  "Conda",
                  "DataFrames",
                  "FileIO",
                  "GLPK",
                  "HDF5",
                  "Infiltrator",
                  "Ipopt",
                  "JLD2",
                  "JSON",
                  "JuMP",
                  "LinearAlgebra",
                  "MathOptInterface",
                  "Memento",
                  "PackageCompiler",
                  "PowerModels",
                  "PyCall",
                  "SCIP",
                  "SQLite",
                  "Tables",
                  "XLSX",
                  "YAML",
                  "Cbc",
                  "SCIP",
                  "HiGHS"]

problems = Dict()

# Set up command-line parser
s = ArgParseSettings()
@add_arg_table s begin
    "-c"
        help = "make clean, i.e. delete existing *.toml files in current directory"
        action = :store_true
end
CLI_args = parse_args(s)

# If the user specifies the -c flag, try to delete Manifest.toml and
#   Project.toml, if they exist in the current directory
if CLI_args["c"] == true
    files_to_delete = ["Manifest.toml", "Project.toml"]
    for dfile in files_to_delete
        try
            rm(abspath(dfile))
            @info "Removed file $dfile"
        catch e
            if occursin("no such file", e.msg)
                current_dir = @__DIR__
                @info "No extant file $dfile in the current directory, $current_dir."
            else
                throw(e)
            end
        end
    end
end

# Activate local environment
Pkg.activate(".")

# Add Julia libraries needed for ALEAF
# Add and build the optimizer packages
try
    println("Adding CPLEX")
    Pkg.add("CPLEX")
    println("Building CPLEX...")
    Pkg.build("CPLEX")
    println("CPLEX built.")
catch e
    @warn "There was a problem installing or building CPLEX."
    problems["CPLEX"] = e
end

# Add the non-optimizer packages (which don't need to be built)
for i=1:size(julia_pkg_list)[1]
    pkg = julia_pkg_list[i]
    try
        @info string("Adding ", pkg)
        Pkg.add(pkg)
        Pkg.build(pkg)
    catch e
        @warn "There was a problem installing $pkg. Installation of this package will be skipped, which may cause issues at ABCE runtime."
        problems[pkg] = e
    end
end

@info "All libraries added to the environment."

# Ensure all non-default Python packages are installed via Conda.jl
@info "Ensuring Conda packages are installed..."
Pkg.build("Conda")
using Conda

req_file = "requirements.txt"
open(req_file, "r") do filehandle
    conda_list = readlines(filehandle)
    for i=1:size(conda_list)[1]
        cpkg = conda_list[i]
        println(string("Adding ", cpkg))
        try
            Conda.add(cpkg)
        catch e
            @warn string("A problem occurred while trying to install package $cpkg. Installation of this package will be skipped, which may cause issues at ABCE runtime.")
            problems[cpkg] = e
        end
    end
end

@info "========================================================================\n\n"
@info "All packages and Python libraries loaded, with the following exceptions:\n"

for key in keys(problems)
    @warn string(key, ": ", problems[key])
end
