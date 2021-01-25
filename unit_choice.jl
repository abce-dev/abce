# Simple two-option project choice model

println("Loading packages...")
using JuMP, GLPK, LinearAlgebra, DataFrames, CSV
println("Packages loaded successfully.")

###### DATA ######

df = CSV.read("./h_units.csv", DataFrame)

#unit_names = ["big", "small"]
#op_costs = [8., 5.]
#op_revs = [16., 10.]
#unit_life = [30., 20.]

genco_max_exp = 30
d = 0.05
npv_costs = zeros(Float64, 0)
npv_revs = zeros(Float64, 0)

##################

println("Setting up model...")

num_types = size(df[!, :names])[1]      # Number of available unit types

transform!(df, [:op_costs, :unit_life] => ((op_cost, life) -> op_cost .* (1+d .- (1+d) .^ (-1*life))/d) => :npv_costs)
transform!(df, [:op_revs, :unit_life] => ((op_rev, life) -> op_rev .* (1+d .- (1+d) .^ (-1*life))/d) => :npv_revs)

#for i = 1:num_types
#    append!(npv_costs, op_costs[i] * (1+d - (1+d)^(-1*unit_life[i]))/d)
#    append!(npv_revs, op_revs[i] * (1+d - (1+d)^(-1*unit_life[i]))/d)
#end

m = Model(GLPK.Optimizer)
@variable(m, u[1:num_types] >= 0, Int)

@constraint(m, transpose(u) * df[!, :op_costs] <= genco_max_exp)    # Single-year cost ceiling

@objective(m, Max, transpose(u) * (df[!, :npv_revs] - df[!, :npv_costs]))
println("Model set up.")

#### SOLVE THE MODEL
println("Solving problem...")
optimize!(m)
status = termination_status.(m)
unit_qty = value.(u)
total_rev = transpose(unit_qty) * df[!, :npv_revs]
total_cost = transpose(unit_qty) * df[!, :npv_costs]
net_profit = transpose(unit_qty) * (df[!, :npv_revs] - df[!, :npv_costs])

println(status)
println("The GenCo should build ", unit_qty[1], " A units, and ", unit_qty[2], " B units")
println("Total company revenue = ", total_rev)
println("Total company costs = ", total_cost)
println("Company net profit = ", net_profit)
