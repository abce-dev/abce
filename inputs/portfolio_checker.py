import sys
import pandas as pd
import argparse
from pathlib import Path
import yaml


def cli_args():
    parser = argparse.ArgumentParser(
                 description="Check the total system portfolio composition in a given file."
             )
    parser.add_argument(
        "--portfolio_file",
        type=str,
        help="Relative path to system portfolio file.",
        default=Path(Path.cwd() / "agent_specifications.yml")
    )

    parser.add_argument(
        "--unit_specs_file",
        type=str,
        help="Relative path to unit specifications file.",
        default=Path(Path.cwd() / "unit_specs.yml")
    )

    args = parser.parse_args()

    return args


def load_agent_portfolios(portfolio_file):
    # Load data from an agent specification-type file
    pf_full_data = yaml.load(
                       open(portfolio_file, "r"),
                       Loader=yaml.FullLoader
                   )

    # Initialize a dictionary to collect agent portfolios
    pfs = {}
    for agent_id, agent_data in pf_full_data.items():
        if "starting_portfolio" in agent_data.keys():
            pfs[agent_id] = agent_data["starting_portfolio"]

    return pfs


def load_unit_specifications(unit_specs_file):
    # Load data from the unit specs file
    unit_specs_full_data = yaml.load(
                               open(unit_specs_file, "r"),
                               Loader=yaml.FullLoader
                           )

    # Set up the dictionary of unit capacities
    unit_caps = {}
    for unit_type, unit_type_data in unit_specs_full_data.items():
        unit_caps[unit_type] = unit_type_data["capacity"]

    return unit_caps


def get_total_unit_numbers(pfs):
    # Convert portfolio dictionaries into a dataframe, filling missing
    #   values with 0. Rows are unit types, columns are agent ids.
    unit_nums_table = pd.DataFrame.from_dict(pfs).fillna(0).astype("int")

    # Add a column of totals by unit type.
    unit_nums_table["total"] = (unit_nums_table[
                                    [agent_id for agent_id in pfs.keys()]
                                ].sum(axis=1)
                               )

    return unit_nums_table


def get_installed_capacity(pfs, unit_caps):
    # Convert pfs into a dataframe; rows are unit types, columns are agent ids
    pf_df = (pd.DataFrame.from_dict(pfs)
               .fillna(0)
               .astype("int")
               .reset_index()
               .rename(columns={"index": "unit_type"})
            )

    # Unpivot the dataframe to create a joinable long-style tidy dataset
    pf_df = (pd.melt(
                 pf_df,
                 id_vars="unit_type",
                 value_vars = [agent_id for agent_id in pfs.keys()]
             ).rename(columns={"variable": "agent_id", "value": "num_units"})
            )

    # Convert unit capacity data into a dataframe
    unit_caps = (pd.DataFrame.from_dict(unit_caps, orient="index")
                   .reset_index()
                   .rename(columns={"index": "unit_type", 0: "capacity"})
                )

    # Join the unit capacity data into the portfolio dataframe by unit type
    caps_table = pf_df.merge(unit_caps, on="unit_type", how="left")

    # Create a column of total capacity by agent
    caps_table["total_capacity"] = (caps_table["num_units"] 
                                    * caps_table["capacity"])

    # Pivot caps_table to show installed capacity by unit type and agent
    caps_table = pd.pivot_table(
                     caps_table,
                     values="total_capacity",
                     index="unit_type",
                     columns="agent_id"
                 )

    # Add a column of totals by unit type.
    caps_table["Total capacity (MWe)"] = (caps_table[
                                              [agent_id for agent_id in pfs.keys()]
                                          ].sum(axis=1)
                                         )

    return caps_table


def show_results(unit_nums_table, caps_table):
    print("\nTotal number of units by type:")
    print(unit_nums_table)
    print("\n")

    print("\nTotal installed capacity by unit type:")
    print(caps_table)
    print("\n")

    grand_total = sum(caps_table["Total capacity (MWe)"])
    print(f"Grand total system nameplate capacity = {grand_total} MWe.")



def check_portfolio(args):
    # Read in the agent portfolios
    pfs = load_agent_portfolios(args.portfolio_file)

    # Read in unit capacity data
    unit_caps = load_unit_specifications(args.unit_specs_file)

    # Get number of units by type, per agent and system total
    unit_nums_table = get_total_unit_numbers(pfs)

    # Get installed capacity by type, per agent and system total
    caps_table = get_installed_capacity(pfs, unit_caps)

    show_results(unit_nums_table, caps_table)


if __name__ == "__main__":
    args = cli_args()
    check_portfolio(args)
