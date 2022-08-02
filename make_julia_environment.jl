using Pkg, Logging, ArgParse

julia_pkg_list = ["ArgParse",
                  "CSV",
                  "Cbc",
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
                  "SQLite",
                  "Tables",
                  "XLSX",
                  "YAML"]

conda_list = ["numpy",
              "scipy",
              "matplotlib",
              "pandas"]

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
@info "Adding CPLEX"
Pkg.add("CPLEX")
@info "Building CPLEX..."
Pkg.build("CPLEX")
@info "CPLEX built."

# Add the non-optimizer packages (which don't need to be built)
for i=1:size(julia_pkg_list)[1]
    @info string("Adding ", julia_pkg_list[i])
    Pkg.add(julia_pkg_list[i])
    Pkg.build(julia_pkg_list[i])
end

@info "All libraries added to the environment."

# Ensure all non-default Python packages are installed via Conda.jl
@info "Ensuring Conda packages are installed..."
Pkg.build("Conda")
using Conda
for i=1:size(conda_list)[1]
    @info string("Adding ", conda_list[i])
    Conda.add(conda_list[i])
end

@info "All packages and Python libraries loaded successfully."
