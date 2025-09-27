import pandas as pd
import yaml

default_filename = "../inputs/agent_specifications.yml"


def check_filename(fname):
    status = 0
    try:
        handle = yaml.load(open(fname, "r"), Loader=yaml.FullLoader)
    except:
        status = 1
        handle = None

    return status, handle


def get_file(mode):
    status = 1
    
    # mode settings
    if mode == "portfolio":
        default_filename = "../inputs/agent_specifications.yml"
        thing_to_load = "total portfolio"
    elif mode == "unit_specs":
        default_filename = "../inputs/unit_specs.yml"
        thing_to_load = "unit specs data"
    else:
        print(f"File retrieval mode {mode} is not supported.")
        exit()

    while status != 0:
        fname = input(f"Provide filename from which to load the {thing_to_load}, or press ENTER to use the default of {default_filename}: ")
        if fname == "":
            fname = default_filename
        status, handle = check_filename(fname)
        if status != 0:
            print("File could not be loaded as a valid YAML input file.")

    return handle


def consolidate_portfolio(handle):
    portfolio = {}

    for agent, agent_data in handle.items():
        apf = agent_data["starting_portfolio"]
        for unit_type, num in apf.items():
            if unit_type in portfolio.keys():
                portfolio[unit_type] += num
            else:
                portfolio[unit_type] = num

    return portfolio


def show_options():
    valid_inputs = ["1", "2", "3"]
    selection = ""
    while selection not in valid_inputs:
        selection = input(f"Select an option:\n1) Show total portfolio\n2) Divide evenly\n3) Hand-specify individual agents\n")

    selection = int(selection)

    return selection


def divide_evenly(system_portfolio):
    num_agents = ""
    check = 0

    while not check:
        num_agents = input("Specify number of agents: ")
        try:
            num_agents = int(num_agents)
            check = 1
        except:
            print("Please enter an integer.")

    # Divide the system portfolio evenly for the specified number of agents
    divided_portfolio = {}
    added_units = {}
    total_portfolio = {}

    for unit_type, num_units in system_portfolio.items():
        total_num_units = num_units
        divides_evenly = False
        while not divides_evenly:
            divided_num_units = total_num_units / num_agents
            if divided_num_units.is_integer():
                divides_evenly = True
                divided_num_units = int(divided_num_units)
            else:
                total_num_units += 1

        # Record the number of units in the divided portfolio
        divided_portfolio[unit_type] = divided_num_units

        # Record any units added to the portfolio for divisibility
        if divided_num_units * num_agents != num_units:
            added_units[unit_type] = divided_num_units * num_agents - num_units

        total_portfolio[unit_type] = divided_num_units * num_agents

    return num_agents, divided_portfolio, added_units, total_portfolio


def get_portfolio_capacity(pf, unit_specs):
    pf_caps = {}
    for unit_type, unit_num in pf.items():
        pf_caps[unit_type] = {}
        pf_caps[unit_type]["ICAP"] = unit_num * unit_specs[unit_type]["capacity"]
        pf_caps[unit_type]["UCAP"] = unit_num * unit_specs[unit_type]["capacity"] * unit_specs[unit_type]["capacity_factor"]

    return pf_caps


def compare_portfolio_capacity(pf1, pf2, unit_specs):
    pf1_ICAP = 0
    pf1_UCAP = 0
    pf2_ICAP = 0
    pf2_UCAP = 0

    pf1_caps = get_portfolio_capacity(pf1, unit_specs)
    pf2_caps = get_portfolio_capacity(pf2, unit_specs)

    for unit_type, unit_type_data in pf1_caps.items():
        pf1_ICAP += unit_type_data["ICAP"]
        pf1_UCAP += unit_type_data["UCAP"]

    for unit_type, unit_type_data in pf2_caps.items():
        pf2_ICAP += unit_type_data["ICAP"]
        pf2_UCAP += unit_type_data["UCAP"]

    return pf1_ICAP, pf1_UCAP, pf2_ICAP, pf2_UCAP


def show_divided_portfolio(num_agents, divided_portfolio, added_units, total_portfolio):
    pprint(divided_portfolio, "Divided portfolio among "+str(num_agents)+" agents")
    pprint(added_units, "Units added to ensure even divisibility")
    pprint(total_portfolio, "New total system portfolio")


def show_capacity_comparison(pf1_name, pf1_ICAP, pf1_UCAP, pf2_name, pf2_ICAP, pf2_UCAP):
    print(f"{pf1_name} portfolio installed capacity: {pf1_ICAP} MWe")
    print(f"{pf2_name} portfolio installed capacity: {pf2_ICAP} MWe")
    print(f"Difference: {pf2_ICAP - pf1_ICAP} MWe, {round((pf2_ICAP - pf1_ICAP) / pf1_ICAP * 100, 1)}%")
    print(f"{pf1_name} portfolio unforced capacity: {pf1_UCAP} MWe")
    print(f"{pf2_name} portfolio unforced capacity: {pf2_UCAP} MWe")
    print(f"Difference: {pf2_UCAP - pf1_UCAP} MWe, {round((pf2_UCAP - pf1_UCAP) / pf1_UCAP * 100, 1)}%")


def pprint(pf, name):
    # Pretty print for portfolio dictionaries

    # Determine the longest possible line length
    longest = 0
    for unit_type, unit_num in pf.items():
        if len(unit_type) + len(str(unit_num)) > longest:
            longest = len(unit_type) + len(str(unit_num))

    longest += 3

    # Print the portfolio with consistent line lengths
    print(f"{name}:")
    for unit_type, unit_num in pf.items():
        space = " "*(longest - len(unit_type) - len(str(unit_num)))
        print(f"{unit_type}:{space}{unit_num}")



def run():
    is_end = 0

    agent_specs_data = get_file("portfolio")
    unit_specs = get_file("unit_specs")
    system_portfolio = consolidate_portfolio(agent_specs_data)

    selection = show_options()
    if selection == 1:  # Show total portfolio
        pprint(system_portfolio, "Total system portfolio")
    elif selection == 2:   # Divide evenly
        num_agents, divided_portfolio, added_units, total_portfolio = divide_evenly(system_portfolio)
        show_divided_portfolio(num_agents, divided_portfolio, added_units, total_portfolio)
        init_sys_ICAP, init_sys_UCAP, div_sys_ICAP, div_sys_UCAP = compare_portfolio_capacity(system_portfolio, total_portfolio, unit_specs)
        show_capacity_comparison("initial system", init_sys_ICAP, init_sys_UCAP, "divided system", div_sys_ICAP, div_sys_UCAP)
    else:   # Hand-specify individual agents
        print("This option isn't implemented yet.")


if __name__ == "__main__":
    run()
