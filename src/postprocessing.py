import numpy as np
import pandas as pd
from pathlib import Path
import yaml
import sqlite3
import matplotlib.pyplot as plt
import logging

unit_type_colors = {
    "coal": "#555051",
    "ngcc": "#03cec8",
    "ngct": "#0340c8",
    "wind": "#16b904",
    "solar": "#ffc800",
    "conventional_nuclear": "#ff6800",
    "advanced_nuclear": "#f1a66d",
    "PWR_C2N0_single": "#ff9aa7",
    "PWR_C2N1_single": "#ff0022",
    "HTGR_C2N0_single": "#cb8ae9",
    "HTGR_C2N2_single": "#a50eec",
    "SFR_C2N0_single": "#ebf697",
    "SFR_C2N3_single": "#c6e000"
}


def write_raw_db_to_excel(settings, db):
    # Get the names of all database tables
    db_tables = pd.read_sql_query(
                    "SELECT name FROM sqlite_master WHERE type='table';",
                     db
                )

    # Set up the path to the ultimate outputs directory
    out_file = Path(
                   Path.cwd() /
                   "outputs" /
                   settings["simulation"]["scenario_name"] /
                   settings["file_paths"]["output_file"]
               )

    # Write each table to a tab in the excel sheet
    with pd.ExcelWriter(out_file) as writer:
        for i in range(len(db_tables)):
            # Get table name
            table = db_tables.loc[i, "name"]

            # Get table data
            table_data = pd.read_sql_query(
                             f"SELECT * FROM {table}",
                             db
                         )

            # Write table data to excel tab
            table_data.to_excel(writer, sheet_name=table, engine="openpyxl")


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


def get_portfolio_profile(settings, db, agent_id, unit_specs):
    # Retrieve a long dataframe showing the number of installed units
    #   by type by year
    portfolios = pd.DataFrame()

    # If the agent_id argument is specified, filter on that specific agent
    if agent_id == None:
        agent_filter = ""
    else:
        agent_filter = f"agent_id = {agent_id} AND "

    # Set the total time horizon to num_steps + the construction duration of
    #   the fastest-to-build project
    horizon = (min(unit_specs.construction_duration) 
               + int(settings["simulation"]["num_steps"])
               + 1
              )

    # Read and process the portfolio for each year in the horizon
    for i in range(horizon):
        year_pf = pd.read_sql_query(
                      f"SELECT unit_type, COUNT(unit_type) " +
                      f"FROM assets WHERE " +
                      agent_filter +
                      f"completion_pd <= {i} AND " +
                      f"retirement_pd > {i} AND " +
                      f"cancellation_pd > {i} " +
                       "GROUP BY unit_type",
                      db
                  )

        # Pivot the dataframe to long format 
        year_pf = pd.pivot_table(
                      year_pf,
                      values="COUNT(unit_type)",
                      index="unit_type",
                      aggfunc=np.sum
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
    portfolios["total_capacity"] = (portfolios["num_units"] 
                                    * portfolios["capacity"]
                                   )
    portfolios = pd.pivot_table(
                     portfolios,
                     values="total_capacity",
                     index="year",
                     columns=["unit_type"]
                 )

    return portfolios


def plot_portfolio_profile(settings, agent_id, portfolio):
    # Set up figure-specific strings according to the agent_id specified
    if agent_id == None:
        title = "Total system portfolio evolution"
        filename = "total_system_portfolio_evolution.png"
    else:
        title = f"Agent {agent_id} portfolio evolution"
        filename = f"agent_{agent_id}_portfolio_evolution.png"

    # Remove empty data columns
    for column in list(portfolio.columns):
        if sum(portfolio[column]) == 0:
            portfolio = portfolio.drop(column, axis=1)

    # Set up the figure
    fig = plt.figure(constrained_layout=True, dpi=250)

    # Add the data
    portfolio.plot.bar(
        stacked=True,
        ax=fig.gca(),
        rot=0,
        color=unit_type_colors,
        edgecolor='black'
    )

    # Add titles and axis labels
    fig.suptitle(title)
    fig.gca().set_ylabel("Installed capacity (MWe)")

    # Move the legend outside of the chart area
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    # Save the figure
    fig.get_figure().savefig(
        Path(
            "outputs",
            settings["simulation"]["scenario_name"],
            filename
        )
    )


def postprocess_portfolios(db, settings, unit_specs, agent_id):
    if agent_id == None:
        msg = "total system portfolio"
    else:
        msg = f"agent {agent_id}'s portfolio"

    logging.debug(f"Procesing data for {msg}...")
    portfolio_profile = get_portfolio_profile(
                            settings,
                            db,
                            agent_id,
                            unit_specs
                        )

    plot_portfolio_profile(settings, agent_id, portfolio_profile)
    logging.debug(f"Plot for {msg} saved.")


def postprocess_results(abce_model, settings):
    logging.info("Postprocessing results...")

    # Save the raw database as an Excel format for easier viewing and manual
    #   postprocessing/debugging
    write_raw_db_to_excel(settings, abce_model.db)

    # Get a list of all agent ids
    agent_list = get_agent_list(abce_model.db)

    # Get an subset of unit_specs columns
    unit_specs = get_unit_specs(abce_model.db)

    # Plot portfolio evolution for all agents, plus the overall system
    for agent_id in agent_list + [None]:
        postprocess_portfolios(abce_model.db, settings, unit_specs, agent_id)

    logging.info("Postprocessing complete.")


if __name__ == "__main__":
    settings = yaml.load(open("../settings.yml", "r"), Loader=yaml.FullLoader)
    db = sqlite3.connect("../outputs/ABCE_ERCOT_allC2N/abce_db.db")
    postprocess_results(db, settings)



