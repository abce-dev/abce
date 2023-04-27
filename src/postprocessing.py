import numpy as np
import pandas as pd
from pathlib import Path
import yaml
import sqlite3
import matplotlib.pyplot as plt

def write_raw_db_to_excel(abce_model, settings):
    # Get the names of all database tables
    db_tables = pd.read_sql_query(
                    "SELECT name FROM sqlite_master WHERE type='table';",
                     abce_model.db
                )

    # Set up the path to the ultimate outputs directory
    out_file = Path(
                   settings["file_paths"]["ABCE_abs_path"] /
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
                             abce_model.db
                         )

            # Write table data to excel tab
            table_data.to_excel(writer, sheet_name=table, engine="openpyxl")


def get_agent_list(db):
    # Get a list of all unique agent IDs
    agent_list = pd.read_sql_query("SELECT agent_id FROM agent_params", db)
    agent_list = agent_list.agent_id.values.tolist()

    return agent_list


def set_horizon(settings, db):
    # Set the total time horizon to num_steps + the construction duration of
    #   the fastest-to-build project
    xtr_durations = pd.read_sql_query(
                        "SELECT construction_duration FROM unit_specs",
                        db
                    )
    offset = min(xtr_durations.construction_duration)

    horizon = int(settings["simulation"]["num_steps"] + offset)

    return horizon


def get_unit_specs(db):
    unit_specs_full = pd.read_sql_query("SELECT * FROM unit_specs", db)

    unit_specs = unit_specs_full[["unit_type", "capacity"]]

    return unit_specs


def get_agent_portfolio_profile(db, agent_id, unit_specs, horizon):
    # Retrieve a long dataframe showing the agent's number of installed units
    #   by type by year
    portfolios = pd.DataFrame()

    horizon = 10

    for i in range(horizon+1):
        year_pf = pd.read_sql_query(
                      f"SELECT unit_type, COUNT(unit_type) " +
                      f"FROM assets " +
                      f"WHERE agent_id = {agent_id} " +
                      f"AND completion_pd <= {i} " +
                      f"AND retirement_pd > {i} " +
                      f"AND cancellation_pd > {i} " +
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


def plot_agent_portfolios(agent_id, agent_portfolios):
    agent_pf_profile = agent_portfolios[agent_id]

    # Remove empty data columns
    for column in list(agent_pf_profile.columns):
        if sum(agent_pf_profile[column]) == 0:
            agent_pf_profile = agent_pf_profile.drop(column, axis=1)

    # Set up the figure
    fig = plt.figure(constrained_layout=True, dpi=250)

    # Add the data
    agent_pf_profile.plot.bar(stacked=True, ax=fig.gca())

    # Add titles and axis labels
    fig.suptitle(f"Agent {agent_id} portfolio evolution: IRA policies only")
    fig.gca().set_ylabel("Installed capacity (MWe)")

    # Move the legend outside of the chart area
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    # Save the figure
    fig.get_figure().savefig(f"agent_{agent_id}_pf_evolution.png")


def postprocess_results(db, settings):
    agent_list = get_agent_list(db)

    horizon = set_horizon(settings, db)

    unit_specs = get_unit_specs(db)

    agent_portfolios = {}

    for agent_id in agent_list:
        print(f"Processing data for agent {agent_id}...")
        agent_portfolios[agent_id] = get_agent_portfolio_profile(
                                         db,
                                         agent_id,
                                         unit_specs,
                                         horizon
                                     )

        plot_agent_portfolios(agent_id, agent_portfolios)
        print(f"Plot for agent {agent_id} saved.")


if __name__ == "__main__":
    settings = yaml.load(open("../settings.yml", "r"), Loader=yaml.FullLoader)
    db = sqlite3.connect("../outputs/ABCE_ERCOT_allC2N/abce_db.db")
    postprocess_results(db, settings)



