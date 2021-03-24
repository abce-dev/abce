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

# File names
unit_data_file = "./data/unit_specs.csv"
fuel_cost_file = "./data/fuel_costs.csv"
demand_data_file = "./data/default_demand.csv"

# Load the inputs
db = load_db()
pd = get_current_period()
agent_id = get_agent_id()

# Set up agent-specific data
# Get a list of all ongoing construction projects for the current agent
agent_projects = get_active_projects_list(db, agent_id)

# Get agent financial parameters
agent_params = get_agent_params(db, agent_id)
d = agent_params[1, :discount_rate]

# Set an average e- price
avg_e_price = 0.07    # $/kWh, avg electricity price

# Read in fuel costs from file
c_fuel = CSV.read(fuel_cost_data, DataFrame)

# System parameters
# Read unit operational data (unit_data) and number of unit types (num_types)
unit_data, num_types = load_unit_type_data(unit_data_file)

# Load the demand data
available_demand = load_demand_data(demand_filename)

# Ensure that forecast horizon is long enough to accommodate the end of life
#   for the most long-lived unit
fc_pd = set_forecast_period(unit_data)

# Extend the unserved demand data to match the total forecast period (constant projection)
available_demand = forecast_demand(available_demand, fc_pd)

# Allocate the correct fuel cost to each generator according to its fuel type
unit_data = allocate_fuel_costs(unit_data)

# Add empty column for project NPVs in unit_data
unit_data[!, :FCF_NPV] = zeros(num_types)

# Create per-unit financial statement tables
unit_FS_dict = create_unit_FS_dict(unit_data)

# Populate financial statements with top-line data
# This data is deterministic, and is on a per-unit basis (not per-kW or per-kWh)
# Therefore, you can multiply each of these dataframes by its corresponding
#    u[] value to determine the actual impact of the units chosen on the GC's
#    final financial statement
for i = 1:num_types
    name = unit_data[i, :name]
    fs = unit_FS_dict[name]

    # Generate events during the construction period
    add_xtr_events(unit_data, i, fs, agent_params)

    # Generate prime-mover events from the end of construction to the end of the unit's life
    for j = (unit_data[i, :d_x] + 1):(unit_data[i, :d_x] + unit_data[i, :unit_life])
        # Apply constant debt payment (sinking fund at cost of debt)
        # This amount is always calculated based on the amount of debt outstanding at
        #    the project's completion
        fs[j, :debt_payment] = fs[unit_data[i, :d_x], :remaining_debt_principal] .* d ./ (1 - (1+d) .^ (-1*unit_data[i, :unit_life]))
        # Determine portion of payment which pays down interest (instead of principal)
        fs[j, :interest_due] = fs[j-1, :remaining_debt_principal] * d
        # Update the amount of principal remaining at the end of the year
        fs[j, :remaining_debt_principal] = fs[j-1, :remaining_debt_principal] - (fs[j, :debt_payment] - fs[j, :interest_due])

        # Set kWh of generation for the year
        fs[j, :gen] = unit_data[i, :capacity] * unit_data[i, :CF] * 8760 * 1000

        # Apply straight-line depreciation
        fs[j, :depreciation] = fs[unit_data[i, :d_x], :xtr_exp] ./ unit_data[i, :unit_life]
    end

    # Apply reactive functions to the rest of the dataframe
    # Compute total revenue for the year
    transform!(fs, [:gen] => ((gen) -> avg_e_price .* gen) => :Revenue)
    # Compute total cost of fuel used
    transform!(fs, [:gen] => ((gen) -> unit_data[i, :uc_fuel] .* unit_data[i, :heat_rate] .* gen ./ 1000000) => :Fuel_Cost)
    # Compute total VOM cost incurred during generation
    transform!(fs, [:gen] => ((gen) -> unit_data[i, :VOM] .* gen) => :VOM_Cost)
    # Compute total FOM cost for the year (non-reactive)
    fs[!, :FOM_Cost] = zeros(size(fs)[1])
    fs[(unit_data[i, :d_x]+1):(unit_data[i, :d_x]+unit_data[i, :unit_life]), :FOM_Cost] = ones(unit_data[i, :unit_life]) .* (unit_data[i, :FOM] * unit_data[i, :capacity] * 1000)
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
    unit_data[i, :FCF_NPV] = transpose(fs[!, :FCF]) * fs[!, :d_factor]
end

println(unit_data)
println("Data initialized.")

###### Set up the model
println("Setting up model...")
m = Model(GLPK.Optimizer)
@variable(m, u[1:num_types] >= 0, Int)
@variable(m, z[1:fc_pd])

# Restrict total construction to be less than maximum available demand (subject to capacity factor)
for i = 1:size(unit_FS_dict[unit_data[1, :name]])[1]
    @constraint(m, sum(u[j] * unit_FS_dict[unit_data[j, :name]][i, :gen] for j=1:3) / (8760*1000) <= available_demand[i])
end
# Constraint on max amount of interest payable per year
@constraint(m, transpose(u) * (unit_data[!, :uc_x] .* unit_data[!, :capacity] .* de_ratio .* d ./ (1 .- (1+d) .^ (-1 .* unit_data[!, :unit_life]))) <= agent_params[1, :interest_cap])

@objective(m, Max, transpose(u) * unit_data[!, :FCF_NPV])
println("Model set up.")


###### Solve the model
println("Solving problem...")
optimize!(m)
status = termination_status.(m)
unit_qty = value.(u)


###### Display the results
println(status)
println("Units to build:")
println(hcat(select(unit_data, :name), DataFrame(units = unit_qty)))
println("Total NPV of all built projects = ", transpose(unit_qty) * unit_data[!, :FCF_NPV])


###### Save the new units into the `assets` and `WIP_projects` DB tables
for i = 1:num_units
    for j = 1:unit_qty[i]
        next_id = get_next_available_id(db)
        # Update `WIP_projects` table
        rcec = unit_data[i, :uc_x] * unit_data[i, :capacity]
        rtec = unit_data[i, :d_x]
        WIP_projects_vals = (next_id, agent_id, pd, rcec, rtec, rcec / 10)
        DBInterface.execute(db, "INSERT INTO WIP_projects VALUES (?, ?, ?, ?, ?, ?)", WIP_projects_vals)

        # Update `assets` table
        unit_type = unit_data[i, :name]
        completion_pd = pd + unit_data[i, :d_x]
        cancellation_pd = 9999
        retirement_pd = pd + unit_data[i, :d_x] + unit_data[i, :unit_life]
        cap_pmt = 0
        assets_vals = (next_id, agent_id, unit_type, completion_pd, cancellation_pd, retirement_pd, cap_pmt)
        DBInterface.execute(db, "INSERT INTO assets VALUES (?, ?, ?, ?, ?, ?, ?)", assets_vals)
    end
end

##### Authorize ANPE for all current WIP projects
# Retrieve all WIP projects
WIP_projects = get_WIP_projects_list(db, agent_id)
# Authorize ANPE for the upcoming period (default: $1B/year)
authorize_anpe(db, agent_id, pd, WIP_projects, unit_data)




show_table(db, "assets")
show_table(db, "WIP_projects")








println("\n Julia: finishing")
