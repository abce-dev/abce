using Pkg

# Activate local environment
Pkg.activate(".")

julia_pkg_list = ["CSV",
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
                  "PowerModels",
                  "PyCall",
                  "SQLite",
                  "XLSX",
                  "YAML"]

conda_list = ["numpy",
              "scipy",
              "matplotlib",
              "pandas"]

# Add Julia libraries needed for ALEAF
# Add and build the optimizer packages
# CPLEX version 0.6.0 is required for compatibility with CPLEX 12.8
println("Adding CPLEX")
Pkg.add(Pkg.PackageSpec(;name="CPLEX", version="0.6.0"))
println("Building CPLEX...")
Pkg.build("CPLEX")
println("CPLEX built.")

println("Adding Gurobi")
Pkg.add("Gurobi")
println("Building Gurobi...")
Pkg.build("Gurobi")
println("Gurobi built.")

# Add the non-optimizer packages (which don't need to be built)
for i=1:size(julia_pkg_list)[1]
    println(string("Adding ", julia_pkg_list[i]))
    Pkg.add(julia_pkg_list[i])
end

println("All libraries added to the environment.")

# Ensure all non-default Python packages are installed via Conda.jl
println("Ensuring Conda packages are installed...")
using Conda
for i=1:size(conda_list)[1]
    println(string("Adding ", conda_list[i]))
    Conda.add(conda_list[i])
end

println("All packages and Python libraries loaded successfully.")
