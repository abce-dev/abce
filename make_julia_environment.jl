using Pkg, Logging, ArgParse

julia_pkg_list = [
    "ArgParse",
    "CPLEX",
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
    "HiGHS"
]

# Initialize a dictionary to track any problems arising in the process
problems = Dict()

# Delete preexisting .toml files to avoid contamination
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

# Activate local environment
Pkg.activate(".")

# Add all Julia packages to the environment, and ensure they are all built
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
