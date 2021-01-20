# Simple two-option project choice model

###### DATA ######

op_costs = [10., 5.]
op_revs = [16., 10.]
unit_life = [30., 20.]
genco_max_exp = 20
d = 0.05
npv_costs = zeros(Float64, 0)
npv_revs = zeros(Float64, 0)

##################

for i = 1:2
    append!(npv_costs, op_costs[i] * (1+d - (1+d)^(-1*unit_life[i]))/d)
    append!(npv_revs, op_revs[i] * (1+d - (1+d)^(-1*unit_life[i]))/d)
end

println("Loading packages...")
using JuMP, GLPK, LinearAlgebra
println("Packages loaded successfully.")

println("Setting up model...")
m = Model(GLPK.Optimizer)
@variable(m, u[1:2] >= 0, Int)

@constraint(m, transpose(u) * op_costs <= genco_max_exp)    # Single-year cost ceiling

@objective(m, Max, transpose(u) * (npv_revs - npv_costs))
println("Model set up.")

#### SOLVE THE MODEL
println("Solving problem...")
optimize!(m)
status = termination_status.(m)
unit_qty = value.(u)
total_rev = transpose(unit_qty) * npv_revs
total_cost = transpose(unit_qty) * npv_costs
net_profit = transpose(unit_qty) * (npv_revs - npv_costs)

println(status)
println("The GenCo should build ", unit_qty[1], " A units, and ", unit_qty[2], " B units")
println("Total company revenue = ", total_rev)
println("Total company costs = ", total_cost)
println("Company net profit = ", net_profit)
