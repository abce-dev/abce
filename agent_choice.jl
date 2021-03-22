# Agent decision model

println("\n\n-----------------------------------------------------------------")
println("Julia agent choice algorithm: starting")
println("Loading packages...")
using JuMP, GLPK, LinearAlgebra, DataFrames, CSV, Printf, YAML
# Include localy module of ABCE functions
include("./ABCEfunctions.jl")
using .ABCEfunctions
println("Packages loaded successfully.")

###### Set up inputs
println("Initializing data...")

# Load the database
db = load_db()
pd = get_current_period()
agent_id = get_agent_id()

# Set up agent-specific data
# Get a list of all ongoing construction projects for the current agent
agent_projects = get_active_projects_list(db, agent_id)

# Unit type data file name
unit_data_file = "./data/h_units.csv"

# Get agent financial parameters
agent_params = get_agent_params(db, agent_id)


#d = 0.05              # Discount rate
#de_ratio = .5         # Max debt/equity ratio
#tax_rate = 0.21       # Corporate tax rate
#int_cap = 9000000     # $/year, max amount of debt interest allowed
avg_e_price = 0.07    # $/kWh, avg electricity price

# Read in fuel costs from file
c_fuel = CSV.read("./data/fuel_costs.csv", DataFrame)

# System parameters
# Read unit operational data (df) and number of unit types (num_types)
df, num_types = load_unit_type_data(unit_data_file)
# Locate the appropriate demand data file
demand_filename = get_demand_data_file()
# Load the demand data
available_demand = load_demand_data(demand_filename)

# Ensure that forecast horizon is long enough to accommodate the end of life
#   for the most long-lived unit
fc_pd = set_forecast_period(df)

# Extend the unserved demand data to match the total forecast period (constant projection)
available_demand = forecast_demand(available_demand, fc_pd)

# Allocate the correct fuel cost to each generator according to its fuel type
df[!, :uc_fuel] = zeros(size(df)[1])
for i = 1:size(df)[1]
    df[i, :uc_fuel] = c_fuel[df[i, :fuel_type]]
end

# Add empty column for project NPVs in df
df[!, :FCF_NPV] = zeros(size(df)[1])

# Create per-unit financial statement tables
fs_dict = create_unit_FS_dfs(df, fc_pd)

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
println("\n Julia: finishing")
