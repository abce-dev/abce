##########################################################################
# Copyright 2023 Argonne National Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################

import os
import numpy as np
import pandas as pd
from pathlib import Path
import yaml
import sqlite3
import matplotlib.pyplot as plt
import logging
import argparse

unit_type_colors = {
    "coal": "#555051",
    "NGCC": "#03cec8",
    "NGCT": "#0340c8",
    "wind": "#16b904",
    "legacy wind": "#0c6802",
    "solar": "#ffc800",
    "legacy solar": "#a98503",
    "conventional nuclear": "#ff6800",
    "generic SMR": "#f1a66d",
    "PWR offsite C2N": "#ff9aa7",
    "PWR onsite C2N": "#ff0022",
    "HTGR offsite C2N": "#cb8ae9",
    "HTGR onsite C2N": "#a50eec",
    "SFR offsite C2N": "#ebf697",
    "SFR onsite C2N": "#c6e000",
    "PUN units": "#000000",
}


def write_raw_db_to_excel(settings, db):
    # Get the names of all database tables
    db_tables = pd.read_sql_query(
        "SELECT name FROM sqlite_master WHERE type='table';", db
    )

    # Set up the path to the ultimate outputs directory
    out_file = Path(
        Path.cwd()
        / "outputs"
        / settings["simulation"]["scenario_name"]
        / settings["file_paths"]["output_file"]
    )

    # Write each table to a tab in the excel sheet
    with pd.ExcelWriter(out_file) as writer:
        for i in range(len(db_tables)):
            # Get table name
            table = db_tables.loc[i, "name"]

            # Get table data
            table_data = pd.read_sql_query(f"SELECT * FROM {table}", db)

            # Write table data to excel tab
            try:
                table_data.to_excel(writer, sheet_name=table, engine="openpyxl")
            except:
                logging.info(f"Unable to save table {table} to excel, likely due to excessive length.")


###############################################################################
# Functions for postprocessing and plotting
###############################################################################


def get_agent_list(db):
    # Get a list of all unique agent IDs
    agent_list = pd.read_sql_query("SELECT agent_id FROM agent_params", db)
    agent_list = agent_list.agent_id.values.tolist()

    return agent_list


def get_unit_specs(db):
    unit_specs_full = pd.read_sql_query("SELECT * FROM unit_specs", db)

    unit_specs = unit_specs_full[
        ["unit_type", "construction_duration", "capacity"]
    ]

    return unit_specs


def get_portfolio_profile(db, agent_id, unit_specs):
    # Retrieve a long dataframe showing the number of installed units
    #   by type by year
    portfolios = pd.DataFrame()

    # If the agent_id argument is specified, filter on that specific agent
    if agent_id == None:
        agent_filter = ""
    else:
        agent_filter = f"agent_id = {agent_id} AND "

    # Set the total time horizon to num_years + the construction duration of
    #   the fastest-to-build project
    ads = pd.read_sql_query(
        "SELECT * FROM annual_dispatch_summary",
        db
    )

    horizon = int(
        unit_specs.construction_duration.min(skipna=True)
        + max(ads["period"]) + 1
        + 1
    )

    # Read and process the portfolio for each year in the horizon
    for i in range(horizon):
        year_pf = pd.read_sql_query(
            f"SELECT unit_type, COUNT(unit_type) "
            + f"FROM assets WHERE "
            + agent_filter
            + f"completion_pd <= {i} AND "
            + f"retirement_pd > {i} AND "
            + f"cancellation_pd > {i} "
            + "GROUP BY unit_type",
            db,
        )

        # Pivot the dataframe to long format
        year_pf = pd.pivot_table(
            year_pf,
            values="COUNT(unit_type)",
            index="unit_type",
            aggfunc=np.sum,
        )

        # Merge in the shortened unit_specs data
        year_pf = year_pf.merge(unit_specs, how="outer", on="unit_type")
        year_pf = year_pf.fillna(0)

        # Add a year index column
        year_pf["year"] = i

        # Rename the COUNT column for simplicity
        year_pf = year_pf.rename(columns={"COUNT(unit_type)": "num_units"})

        # Either create or append to the dataframe
        if (len(portfolios) == 0) and (len(year_pf) != 0):
            portfolios = year_pf
        elif len(portfolios) != 0:
            portfolios = pd.concat([portfolios, year_pf])

    # Compute total capacity
    portfolios["total_capacity"] = (
        portfolios["num_units"] * portfolios["capacity"]
    )
    portfolios = pd.pivot_table(
        portfolios, values="total_capacity", index="year", columns=["unit_type"]
    )

    return portfolios


def organize_portfolio_profile(portfolio_profile):
    # Re-order the columns
    col_order = [
        "PUN_unit",
        "conventional_nuclear",
        "coal",
        "PWR_C2N0_single",
        "PWR_C2N1_single",
        "HTGR_C2N0_single",
        "HTGR_C2N2_single",
        "SFR_C2N0_single",
        "SFR_C2N3_single",
        "advanced_nuclear",
        "ngcc",
        "ngct",
        "wind_old",
        "wind",
        "solar_old",
        "solar",
    ]

    cols_to_select = []
    extant_units = portfolio_profile.columns.values.tolist()

    for name in col_order:
        if name in extant_units:
            cols_to_select.append(name)

    portfolio_profile = portfolio_profile[cols_to_select]

    # Rename columns to more readable names
    readable_col_names = {
        "PUN_unit": "PUN units",
        "conventional_nuclear": "conventional nuclear",
        "PWR_C2N0_single": "PWR offsite C2N",
        "PWR_C2N1_single": "PWR onsite C2N",
        "HTGR_C2N0_single": "HTGR offsite C2N",
        "HTGR_C2N2_single": "HTGR onsite C2N",
        "SFR_C2N0_single": "SFR offsite C2N",
        "SFR_C2N3_single": "SFR onsite C2N",
        "advanced_nuclear": "generic SMR",
        "ngcc": "NGCC",
        "ngct": "NGCT",
        "wind_old": "legacy wind",
        "solar_old": "legacy solar",
    }

    portfolio_profile = portfolio_profile.rename(columns=readable_col_names)

    return portfolio_profile


def plot_portfolio_profile(settings, agent_id, portfolio, tag=None, descriptor=None):
    # Set up figure titles and filenames
    if settings is not None:
        tag = settings["simulation"]["scenario_name"]
        descriptor = " ".join(tag.split("_"))

    # Set up figure-specific strings according to the agent_id specified
    if agent_id == None:
        title = f"Total system portfolio evolution\n{descriptor}"
        filename = f"{tag}_total_system_portfolio_evolution.png"
    else:
        title = f"Agent {agent_id} portfolio evolution\n{descriptor}"
        filename = f"{tag}_agent_{agent_id}_portfolio_evolution.png"

    # Remove empty data columns
    for column in list(portfolio.columns):
        if sum(portfolio[column]) == 0:
            portfolio = portfolio.drop(column, axis=1)

    # Set up the figure
    fig, ax = plt.subplots(constrained_layout=True, dpi=250)

    # Rescale MW to GW
    portfolio = portfolio / 1000

    # Add the data
    portfolio.plot.bar(
        stacked=True,
        ax=ax,
        rot=0,
        color=unit_type_colors,
        edgecolor="black",
        width=1,
        linewidth=0.5,
    )

    # Set xticks to every 5 years
    ticks = [i*5 for i in range((len(portfolio) // 5) + 1)]
    ax.set_xticks(ticks)

    # Add titles and axis labels
    fig.suptitle(title)
    ax.set_ylabel("Installed capacity (GWe)")

    # Add the legend
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1], loc="center left", bbox_to_anchor=(1.0, 0.5))

    # Save the figure
    if settings is not None:
        fig.get_figure().savefig(
            Path("outputs", settings["simulation"]["scenario_name"], filename)
        )
    else:
        target_dir = Path(os.getenv("ABCE_DIR")) / "imgs" / tag
        if not target_dir.is_dir():
            os.makedirs(target_dir)
        fig.get_figure().savefig(target_dir / filename)


def plot_portfolios(db, settings, unit_specs, agent_id, tag, descriptor):
    if agent_id == None:
        msg = "total system portfolio"
    else:
        msg = f"agent {agent_id}'s portfolio"

    logging.debug(f"Procesing data for {msg}...")
    portfolio_profile = get_portfolio_profile(
        db, agent_id, unit_specs,
    )

    portfolio_profile = organize_portfolio_profile(portfolio_profile)

    plot_portfolio_profile(settings, agent_id, portfolio_profile, tag, descriptor)
    logging.debug(f"Plot for {msg} saved.")


def postprocess_results(args, abce_model, settings=None, tag=None, descriptor=None):
    logging.info("Postprocessing results...")

    # Save the raw database as an Excel format for easier viewing and manual
    #   postprocessing/debugging
    if settings is not None:
        write_raw_db_to_excel(settings, abce_model.db)

    # Get a list of all agent ids
    agent_list = get_agent_list(abce_model.db)

    # Get an subset of unit_specs columns
    unit_specs = get_unit_specs(abce_model.db)

    # Plot portfolio evolution for all agents, plus the overall system
    if not args.no_plots:
        for agent_id in agent_list + [None]:
            plot_portfolios(abce_model.db, settings, unit_specs, agent_id, tag, descriptor)

    logging.info("Postprocessing complete.")


def cli_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dir",
        type=str,
        required=True
    )

    parser.add_argument(
        "--default",
        "-d",
        action="store_true",
    )

    args = parser.parse_args()
    args.no_plots = False
    return args


class Model(object):
    def __init__(self, db):
        self.db = db


if __name__ == "__main__":
    args = cli_args()

    # Open database and create the fake model object
    db_file = Path(args.dir) / "abce_db.db"
    db = sqlite3.connect(db_file)
    m = Model(db)

    # If not default, ask the user for the scenario name to use in
    #   figure titles
    descriptor = None
    if not args.default:
        descriptor = input("Provide a scenario descriptor to use in the figure titles:\n")

    # Separate out the immediate parent directory name as a tag for filenames
    # This is the same as the scenario name in a live run of ABCE
    tag = db_file.parts[-2]

    postprocess_results(args, m, settings=None, tag=tag, descriptor=descriptor)


