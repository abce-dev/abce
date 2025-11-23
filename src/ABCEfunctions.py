##########################################################################
## Copyright 2023 Argonne National Laboratory
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


import sqlite3
import os
import pandas as pd
from . import scenario_reduction as sr


def get_next_asset_id(db, suggested_next_id):
    # Start a list of possible next ids
    next_id_candidates = [suggested_next_id]

    # Check all possible locations for max asset ids
    tables_to_check = [
        "assets",
        "WIP_projects",
        "asset_updates",
        "WIP_updates",
        "depreciation_projections",
    ]

    for table in tables_to_check:
        id_val = pd.read_sql(f"SELECT MAX(asset_id) FROM {table}", db).iloc[
            0, 0
        ]
        next_id_candidates.append(id_val)

    # The largest of all non-None elements becomes the next asset id
    next_id = max([item for item in next_id_candidates if item is not None]) + 1

    return next_id


def execute_scenario_reduction(args, db, current_pd, fc_pd, settings, unit_specs):
    # Get the number of wind and solar units to allow computation of net
    #   demand
    current_portfolio = pd.read_sql_query(
        f"SELECT unit_type FROM assets "
        + f"WHERE completion_pd <= {fc_pd} "
        + f"AND retirement_pd > {fc_pd}",
        db,
    )
    num_wind = len(current_portfolio.loc[current_portfolio.unit_type == "wind"])
    num_wind_old = len(current_portfolio.loc[current_portfolio.unit_type == "wind_old"])
    num_solar = len(
        current_portfolio.loc[current_portfolio.unit_type == "solar"]
    )
    num_solar_old = len(current_portfolio.loc[current_portfolio.unit_type == "solar_old"])

    # Get the capacity of wind and solar units
    wind_cap = unit_specs["wind"]["capacity"]
    wind_old_cap = unit_specs["wind_old"]["capacity"]
    solar_cap = unit_specs["solar"]["capacity"]
    solar_old_cap = unit_specs["solar_old"]["capacity"]

    # Get peak demand for this period
    peak_demand = pd.read_sql_query(
        f"SELECT demand FROM demand " + f"WHERE period = {fc_pd}", db
    ).iloc[0, 0]

    # Set up directory locations
    ts_data_dir = os.path.join(
        args.inputs_path,
        "ts_data"
    )

    temp_data_dir = os.path.join(
        args.inputs_path,
        "ts_data",
        "scenario_reduction_tmp",
    )

    output_dir = os.path.join(
        args.inputs_path,
        "ts_data",
    )

    plot_output_dir = os.path.join(output_dir, "sr_plots")

    sr.run_scenario_reduction(
        time_resolution="hourly",
        num_scenarios_list=[settings["dispatch"]["num_repdays"]],
        generate_input_data_flag=True,
        data_location_timeseries=ts_data_dir,
        data_input_path=temp_data_dir,
        data_output_path=output_dir,
        plot_output_path=plot_output_dir,
        windCapacity=num_wind * wind_cap + num_wind_old * wind_old_cap,
        solarCapacity=num_solar * solar_cap + num_solar_old * solar_old_cap,
        peakDemand=peak_demand,
        current_pd=current_pd,
        fc_pd=fc_pd
    )


def update_DB_table_inplace(db, cur, table, new_data, where):
    """
    A centralized function to update data in any DB table. The final form of
    the command constructed will be:
        UPDATE table SET key1 = new_data[key1], key2 = new_data[key2], ...,
            keyn = new_data[keyn] WHERE col1 = where[1] AND col2 = where[2]
            AND ... AND colk = where[k]

    Arguments:
      - db (sqlite3 database connection)
      - cur (sqlite3 cursor object for db)
      - table (string): name of table in db to update
      - new_data (dict): column headers and new values to set for each
      - where (dict): search/match criteria to determine which rows are updated
    """

    # Initialize the command
    update_cmd = f"UPDATE {table} SET "

    # Add the list of data to update
    val_list = []
    for column, value in new_data.items():
        if isinstance(value, str):
            val_list.append(f"{column} = '{value}'")
        else:
            val_list.append(f"{column} = {value}")
    update_cmd += ", ".join(val_list)

    # Add the list of matching/filtering conditions
    update_cmd += f" WHERE "
    val_list = []
    for column, value in where.items():
        # Protect strings with single quotes
        if isinstance(value, str):
            val_list.append(f"{column} = '{value}'")
        else:
            val_list.append(f"{column} = {value}")
    update_cmd += " AND ".join(val_list)

    # Execute the constructed command
    cur.execute(update_cmd)

