# Simple two-option project choice model

println("Loading packages...")
using JuMP, GLPK, LinearAlgebra, DataFrames, CSV, Printf, YAML
# Include localy module of ABCE functions
include("./ABCEfunctions.jl")
using .ABCEfunctions
println("Packages loaded successfully.")

###### Set up inputs
println("Initializing data...")

# Get GenCo number from command line
if size(ARGS)[1] != 0
    gc_id = ARGS[1]
else
    println("No GenCo ID provided. Using default data.")
    gc_id = nothing
end

max_demand = 1500     # Maximum available unserved demand
d = 0.05             # Discount rate
de_ratio = .5
debt_cap = .6
tax_rate = 0.21
int_cap = 9000000     # $/year, max amount of debt interest allowed

#fc_pd = 50    # Number of periods in the forecast horizon
avg_e_price = 0.07   # $/kWh, avg electricity price

# Dictionary of fuel costs
c_fuel = Dict([("nuc", .64), ("ng", 2), ("coal", 5)])  # $/MMBTU
c_fuel_sd = Dict([("nuc", 0.01), ("ng", 0.1), ("coal", 0.12)])

# Sequence of available unserved demand
# The sequence starts with the upcoming period, and is padded after its
#   end to match the length of the fs_dict dataframes.
#available_demand = [100, 500, 1000, 1500, 1500, 1500, 1500, 1500]

# System parameters
# Read in from CSV
df = CSV.read("./h_units.csv", DataFrame)   # Updated unit operational data
if gc_id == nothing
    demand_filename = "./default_demand.csv"
else
    demand_filename = string("./gc", gc_id, "_demand.csv")
end
available_demand = CSV.read(demand_filename, DataFrame)[!, :demand]

# Set number of available unit types based on the vertical size of the df array
num_types = size(df)[1]

# Ensure that forecast horizon is long enough to accommodate the end of life
#   the most long-lived unit
durations = df[!, :d_x] + df[!, :unit_life]
fc_pd = maximum(durations)

# Extend the unserved demand data to match the total forecast period
# Assume the final value for unserved demand remains the same for the entire
#   horizon
fill_demand = last(available_demand) * ones(fc_pd - size(available_demand)[1])
available_demand = vcat(available_demand, fill_demand) 

# Allocate the correct fuel cost to each generator according to its fuel type
# RV (will go inside the MC loop)
df[!, :uc_fuel] = zeros(size(df)[1])
for i = 1:size(df)[1]
    df[i, :uc_fuel] = c_fuel[df[i, :fuel_type]]
end

# Add empty column for project NPVs in df
df[!, :FCF_NPV] = zeros(size(df)[1])

# Create per-unit financial statement tables
fs_dict = Dict{String, DataFrame}()
for i = 1:size(df)[1]
    new_df = DataFrame(year = 1:fc_pd, xtr_exp = zeros(fc_pd), gen = zeros(fc_pd), remaining_debt_principal = zeros(fc_pd), debt_payment = zeros(fc_pd), interest_due = zeros(fc_pd), depreciation = zeros(fc_pd))
    name = df[i, :name]
    fs_dict[name] = new_df
end

# Populate financial statements with top-line data
# This data is deterministic, and is on a per-unit basis (not per-kW or per-kWh)
# Therefore, you can multiply each of these dataframes by its corresponding
#    u[] value to determine the actual impact of the units chosen on the GC's
#    final financial statement
for i = 1:size(df)[1]
    name = df[i, :name]
    fs = fs_dict[name]

    # Generate events during the construction period
    for j = 1:df[i, :d_x]
        # Linearly distribute construction costs over the construction duration
        fs[j, :xtr_exp] = df[i, :uc_x] * df[i, :capacity] * 1000 / df[i, :d_x]
        fs[j, :remaining_debt_principal] = sum(fs[k, :xtr_exp] for k in 1:j) * de_ratio
    end

    # Generate prime-mover events from the end of construction to the end of the unit's life
    for j = (df[i, :d_x] + 1):(df[i, :d_x] + df[i, :unit_life])
        # Apply constant debt payment (sinking fund at cost of debt)
        # This amount is always calculated based on the amount of debt outstanding at
        #    the project's completion
        fs[j, :debt_payment] = fs[df[i, :d_x], :remaining_debt_principal] .* d ./ (1 - (1+d) .^ (-1*df[i, :unit_life]))
        # Determine portion of payment which pays down interest (instead of principal)
        fs[j, :interest_due] = fs[j-1, :remaining_debt_principal] * d
        # Update the amount of principal remaining at the end of the year
        fs[j, :remaining_debt_principal] = fs[j-1, :remaining_debt_principal] - (fs[j, :debt_payment] - fs[j, :interest_due])

        # Set kWh of generation for the year
        fs[j, :gen] = df[i, :capacity] * df[i, :CF] * 8760 * 1000

        # Apply straight-line depreciation
        fs[j, :depreciation] = fs[df[i, :d_x], :xtr_exp] ./ df[i, :unit_life]
    end

    # Apply reactive functions to the rest of the dataframe
    # Compute total revenue for the year
    transform!(fs, [:gen] => ((gen) -> avg_e_price .* gen) => :Revenue)
    # Compute total cost of fuel used
    transform!(fs, [:gen] => ((gen) -> df[i, :uc_fuel] .* df[i, :heat_rate] .* gen ./ 1000000) => :Fuel_Cost)
    # Compute total VOM cost incurred during generation
    transform!(fs, [:gen] => ((gen) -> df[i, :VOM] .* gen) => :VOM_Cost)
    # Compute total FOM cost for the year (non-reactive)
    fs[!, :FOM_Cost] = zeros(size(fs)[1])
    fs[(df[i, :d_x]+1):(df[i, :d_x]+df[i, :unit_life]), :FOM_Cost] = ones(df[i, :unit_life]) .* (df[i, :FOM] * df[i, :capacity] * 1000)
    # Compute EBITDA
    transform!(fs, [:Revenue, :Fuel_Cost, :VOM_Cost, :FOM_Cost] => ((rev, fc, VOM, FOM) -> rev - fc - VOM - FOM) => :EBITDA)
    # Compute EBIT
    transform!(fs, [:EBITDA, :depreciation] => ((EBITDA, dep) -> EBITDA - dep) => :EBIT)
    # Compute EBT
    transform!(fs, [:EBIT, :interest_due] => ((EBIT, interest) -> EBIT - interest) => :EBT)
    # Compute taxes owed
    transform!(fs, [:EBT] => ((EBT) -> EBT .* tax_rate) => :tax_owed)
    # Compute net income
    transform!(fs, [:EBT, :tax_owed] => ((EBT, tax) -> EBT - tax) => :Net_Income)
    # Compute FCF
    transform!(fs, [:Net_Income, :interest_due, :xtr_exp] => ((NI, interest, xtr_exp) -> NI + interest - xtr_exp) => :FCF)

    # Add column of compounded discount factors
    transform!(fs, [:year] => ((year) -> (1+d) .^ (-1 .* year)) => :d_factor)
    # 
    df[i, :FCF_NPV] = transpose(fs[!, :FCF]) * fs[!, :d_factor]
end

println(df)
println("Data initialized.")

###### Set up the model
println("Setting up model...")
m = Model(GLPK.Optimizer)
@variable(m, u[1:num_types] >= 0, Int)
@variable(m, z[1:fc_pd])

# Restrict total construction to be less than maximum available demand (subject to capacity factor)
for i = 1:size(fs_dict[df[1, :name]])[1]
    @constraint(m, sum(u[j] * fs_dict[df[j, :name]][i, :gen] for j=1:3) / (8760*1000) <= available_demand[i])
end
#@constraint(m, transpose(u) * (df[!, :capacity] .* df[!, :CF]) <= max_demand)
# Constraint on max amount of interest payable per year
@constraint(m, transpose(u) * (df[!, :uc_x] .* df[!, :capacity] .* de_ratio .* d ./ (1 .- (1+d) .^ (-1 .* df[!, :unit_life]))) <= int_cap)

@objective(m, Max, transpose(u) * df[!, :FCF_NPV])
println("Model set up.")


###### Solve the model
println("Solving problem...")
optimize!(m)
status = termination_status.(m)
unit_qty = value.(u)


###### Display the results
println(status)
println("Units to build:")
println(hcat(select(df, :name), DataFrame(units = unit_qty)))
println("Total NPV of all built projects = ", transpose(unit_qty) * df[!, :FCF_NPV])
