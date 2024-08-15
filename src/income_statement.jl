


module ISpf

using DataFrames, CSV, SQLite


function set_forecast_horizon(unit_specs, other_args)
    pf_horiz = maximum(unit_specs[!, :construction_duration]) + maximum(unit_specs[!, :unit_life]) + 1

    return pf_horiz
end


function create_income_statement(agent_id, current_pd, pf_horiz, other_args)
    pf = DataFrame(
        agent_id = Int[],
        base_pd = Int[],
        fc_pd = Int[],
        delta_pd = Int[],
    )

    for y = 0:(pf_horiz-1)
        row = [agent_id, current_pd, current_pd + y, y]
        push!(pf, row)
    end

    return pf
end


function add_revenue(pf, other_args)
    annual_revenue = 1000000000.0
    pf[!, :revenue] .= annual_revenue

    return pf
end


function add_costs(pf, other_args)
    annual_costs = -500000000.0
    pf[!, :VOM] .= annual_costs
    pf[!, :FOM] .= 0.0
    pf[!, :fuel_cost] .= 0.0
    pf[!, :carbon_tax] .= 0.0

    return pf
end


function add_depreciation(pf, other_args)
    annual_depreciation = -80000000.0
    pf[!, :depreciation] .= annual_depreciation

    return pf
end


function add_interest(pf, other_args)
    annual_interest = -120000000.0
    pf[!, :interest_payment] .= annual_interest

    return pf
end


function compute_tax_credits(pf, other_args)
    pf[!, :tax_credits] .= 500000.0

    return pf
end


function account_for_tax(pf, settings, mode="IS")
#    pf = compute_tax_credits(pf, other_args)

    tax_rate = settings["system"]["tax_rate"]
    tax_credits_discount = settings["system"]["tax_credits_discount"]

    if mode == "subproject"
        # Allow negative tax (assume adequate tax burden in the rest of
        #    the firm)
        transform!(pf, :EBT => ((EBT) -> EBT * tax_rate * (-1)) => :nominal_tax)
        transform!(pf, [:nominal_tax, :tax_credits] => ((tax, tc) -> tax + tc) => :realized_tax)
        pf[!, :tax_credit_sale_revenue] .= 0.0
    else
        # If this is the full income statement, enforce a minimum tax owed of $0
        #   and monetize excess tax credits by selling at the "haircut" rate
        transform!(
            pf,
            :EBT
            => ByRow(EBT -> ifelse.(EBT >= 0, EBT * tax_rate * (-1), 0))
            => :nominal_tax
        )

        transform!(
            pf,
            [:nominal_tax, :tax_credits]
            => ByRow((tax, tc) -> ifelse.((-1)*tax >= tc, tax + tc, 0))
            => :realized_tax
        )

        transform!(
            pf,
            [:nominal_tax, :tax_credits]
             => ByRow((tax, tc) -> ifelse.((-1)*tax < tc, (tc + tax) * (1 - tax_credits_discount) * (1 - tax_rate), 0)) 
             => :tax_credit_sale_revenue
        )
    end

    return pf
end


function EBITDA_to_EBT(pf)
    # Compute EBITDA
    transform!(
        pf,
        [:revenue, :VOM, :FOM, :fuel_cost, :carbon_tax] =>
        ((rev, vom, fom, fc, ctax) -> rev + vom + fom + fc + ctax) =>
        :EBITDA
    )

    # Compute EBIT
    transform!(
        pf,
        [:EBITDA, :depreciation] =>
        ((ebitda, dep) -> ebitda + dep) =>
        :EBIT
    )

    # Compute EBT
    transform!(
        pf,
        [:EBIT, :interest_payment] =>
        ((ebit, int) -> ebit + int) =>
        :EBT
    )

    return pf
end


function compute_net_income(pf)
    transform!(
        pf,
        [:EBT, :realized_tax, :tax_credit_sale_revenue]
        => ((EBT, tax, tc_sale_rev) -> EBT + tax + tc_sale_rev)
        => :net_income
    )

    return pf
end


function reorder_pf(pf)
    select!(
        pf,
        :agent_id,
        :base_pd,
        :fc_pd,
        :delta_pd, 
        :revenue, 
        :VOM,
        :FOM,
        :fuel_cost,
        :carbon_tax,
        :EBITDA,
        :depreciation,
        :EBIT,
        :interest_payment,
        :EBT,
        :tax_credits,
        :realized_tax,
        :net_income
    )

    return pf
end




function forecast_income_statement(db, settings, agent_id, current_pd)
    ###### Set up data and settings
    pf_horiz = 8
    tax_rate = settings["system"]["tax_rate"]
    other_args = nothing

    # Retrieve necessary DB tables
    unit_specs = DBInterface.execute(db, "SELECT * FROM unit_specs") |> DataFrame

    # Determine the forecast horizon
    pf_horiz = set_forecast_horizon(unit_specs, other_args)

    # Create the income statement (year columns only)
    pf = create_income_statement(agent_id, current_pd, pf_horiz, other_args)

    ###### Account for all exogenous factors
    # Account for revenues
    pf = add_revenue(pf, other_args)

    # Account for costs
    pf = add_costs(pf, other_args)

    # Account for depreciation
    pf = add_depreciation(pf, other_args)

    # Account for interest expense
    pf = add_interest(pf, other_args)

    ###### Account for all computed factors
    pf = EBITDA_to_EBT(pf)

    # Account for tax owed
    pf = account_for_tax(pf, settings, "IS")

    # Compute net income
    transform!(
        pf,
        [:EBT, :realized_tax] =>
        ((ebt, tax) -> ebt + tax) =>
        :net_income
    )

    ##### Finalize format and save to file
    pf = reorder_pf(pf)

    println(pf)
    CSV.write("./IS.csv", pf)
end


end
