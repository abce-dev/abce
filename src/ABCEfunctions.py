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
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
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


def process_outputs(settings, output_dir, unit_specs):
    """
    A handler function for postprocessing A-LEAF results stored from the
      individual time-steps of an ABCE simulation run.

    Postprocessing tasks:
      - "Capacity expansion" (i.e. portfolio composition) (in # of units
          and MW installed)
      - System-wide results, e.g. total system costs and weighted-average
          prices
      - Generation unit type results: amount of generation and AS, and
          profitability by unit type
      - Price data: sorted by time-stamp and from high to low
    """

    # Postprocessing settings
    scenario_name = settings["simulation"]["scenario_name"]
    file_types = [
        "dispatch_summary_OP",
        "expansion_result",
        "system_summary_OP",
        "system_tech_summary_OP",
    ]

    # Get a list of each type of file from the ABCE A-LEAF outputs dir
    # This will allow postprocessing to work even if some files are missing
    file_lists = {}
    for ftype in file_types:
        file_lists[ftype] = glob.glob(os.path.join(output_dir, f"*{ftype}*"))

    # Postprocess the expansion results
    expansion_results = process_expansion_results(
        file_lists["expansion_result"], output_dir, scenario_name, unit_specs
    )

    # Postprocess the system-level results
    system_summary_results = process_system_summary(
        file_lists["system_summary_OP"], output_dir, scenario_name
    )

    # Postprocess the generation unit type results
    system_tech_results = process_tech_summary(
        file_lists["system_tech_summary_OP"],
        output_dir,
        scenario_name,
        unit_specs,
    )

    # Postprocess the electricity price data
    unsorted_lmp_data, sorted_lmp_data = process_dispatch_data(
        file_lists["dispatch_summary_OP"], output_dir, scenario_name
    )

    # Write results to xlsx
    writer = pd.ExcelWriter(
        Path(
            self.settings["file_paths"]["ABCE_abs_path"]
            / "outputs"
            / self.settings["simulation"]["scenario_name"]
            / "abce_ppx_outputs.xlsx"
        )
    )
    expansion_results.to_excel(writer, sheet_name="exp_results", index=False)
    system_summary_results.to_excel(
        writer, sheet_name="sys_summary_results", index=False
    )
    system_tech_results["gen"].to_excel(
        writer, sheet_name="tech_generation", index=False
    )
    system_tech_results["rr"].to_excel(
        writer, sheet_name="tech_reg", index=False
    )
    system_tech_results["sr"].to_excel(
        writer, sheet_name="tech_spin", index=False
    )
    system_tech_results["nsr"].to_excel(
        writer, sheet_name="tech_nonspin", index=False
    )
    system_tech_results["rev"].to_excel(
        writer, sheet_name="tech_revenue", index=False
    )
    system_tech_results["prof"].to_excel(
        writer, sheet_name="tech_profit", index=False
    )
    unsorted_lmp_data.to_excel(writer, sheet_name="unsorted_LMP", index=False)
    sorted_lmp_data.to_excel(writer, sheet_name="sorted_LMP", index=False)
    writer.save()

    # Plot PDCs
    plot_pdcs(sorted_lmp_data)


def process_expansion_results(
    exp_file_list, output_dir, scenario_name, unit_specs
):
    """
    Process the "capacity expansion" results output by A-LEAF. For ABCE,
      A-LEAF should never be able to actually build or retire new units, so
      this provides both an easy reminder of the portfolio composition but
      also a diagnostic in case of an incorrect setting that allows A-LEAF
      to build or retire units.
    """
    num_files = len(exp_file_list)
    for i in range(num_files):
        file_name = os.path.join(
            output_dir, f"{scenario_name}__expansion_result__step_{i}.csv"
        )
        df = pd.read_csv(file_name)

        # Use the first file as a seed
        if i == 0:
            exp_df = df
            exp_df = exp_df.rename(columns={"u_i": "step_0"})
        # Fill in only the # of units column from subsequent files
        else:
            exp_df[f"step_{i}"] = df["u_i"]

    # Delete unneeded columns and give unit id more helpful names
    exp_df = exp_df.drop(["u_new_i", "u_ret_i"], axis=1)

    # Warning: this currently assumes the units in unit_specs are listed in the
    #   same order as in the A-LEAF outputs
    # unit_specs should inherit its ordering from the master A-LEAF unit specs
    #   file, but this may be fragile.
    exp_df["unit_id"] = unit_specs["unit_type"]
    exp_df = exp_df.rename(columns={"unit_id": "unit_type"})

    return exp_df


def process_system_summary(ss_file_list, output_dir, scenario_name):
    """
    Process and collect all system-summary data outputs from A-LEAF.

    The average LMP and RMP data in the A-LEAF output files are computed by
      straight average, whereas I need a weighted average. This function also
      loads in the dispatch data to quickly compute the correct weighted-
      average prices for each time step.
    """
    num_files = len(ss_file_list)
    for i in range(num_files):
        file_name = os.path.join(
            output_dir, f"{scenario_name}__system_summary_OP__step_{i}.csv"
        )
        df = pd.read_csv(file_name)
        # For this workflow, df should only have one line. Alert the user if
        # otherwise
        if len(df) != 1:
            logging.warn(
                "It looks like you've run multiple ALEAF scenarios ",
                "in the same run. Postprocessing outputs may be ",
                "incorrectly translated.",
            )

        # Give df a more helpful first-column key
        df = df.rename(columns={"Scenario": "step"})
        df.loc[0, "step"] = i

        # Correctly calculate weighted average electricity and AS prices, and
        #   add them to the appropriate columns
        ds_file_name = os.path.join(
            output_dir, f"{scenario_name}__dispatch_summary_OP__step_{i}.csv"
        )
        dsdf = pd.read_csv(ds_file_name)
        # Set up weighted price columns
        dsdf["wtd_LMP"] = dsdf["g_idht"] * dsdf["LMP_dht"]
        dsdf["wtd_RMP_R"] = dsdf["r_r_idht"] * dsdf["RMP_R_dht"]
        dsdf["wtd_RMP_S"] = dsdf["r_s_idht"] * dsdf["RMP_S_dht"]
        dsdf["wtd_RMP_NS"] = dsdf["r_ns_idht"] * dsdf["RMP_NS_dht"]
        # Compute weighted average prices
        elec_wap = dsdf["wtd_LMP"].sum() / dsdf["g_idht"].sum()
        r_wap = dsdf["wtd_RMP_R"].sum() / dsdf["r_r_idht"].sum()
        s_wap = dsdf["wtd_RMP_S"].sum() / dsdf["r_s_idht"].sum()
        ns_wap = dsdf["wtd_RMP_NS"].sum() / dsdf["r_ns_idht"].sum()
        # Insert corrected values into df
        df.loc[0, "LMP"] = elec_wap
        df.loc[0, "RMP_R"] = r_wap
        df.loc[0, "RMP_S"] = s_wap
        df.loc[0, "RMP_NS"] = ns_wap

        # If this is step 0, use the dataframe as a seed
        if i == 0:
            ss_df = df
        # Append each file's data to the main df
        else:
            ss_df = pd.concat([ss_df, df])

    return ss_df


def process_tech_summary(sts_file_list, output_dir, scenario_name, unit_specs):
    """
    Process all sets of A-LEAF unit-type summary statistics. Cleans up unit
      type names and then organizes results into a unified dataframe.
    """
    num_files = len(sts_file_list)
    for i in range(num_files):
        # Read in the system technology summary for time-step i
        file_name = os.path.join(
            output_dir, f"{scenario_name}__system_tech_summary_OP__step_{i}.csv"
        )
        df = pd.read_csv(file_name)

        # Clean up unit type representation
        df = df.rename(columns={"UnitGroup": "unit_type"})

        # Replace default numbers with unit names
        # Same fragility warning as process_expansion_results()
        df["unit_type"] = unit_specs["unit_type"]

        # If this is step 0, use the file to set up DataFrames for all tracked
        #   items
        if i == 0:
            gen_df = (
                df[["unit_type", "Generation"]]
                .copy()
                .rename(columns={"Generation": "step_0"})
            )
            rr_df = (
                df[["unit_type", "Reserve_Reg"]]
                .copy()
                .rename(columns={"Reserve_Reg": "step_0"})
            )
            sr_df = (
                df[["unit_type", "Reserve_Spin"]]
                .copy()
                .rename(columns={"Reserve_Spin": "step_0"})
            )
            nsr_df = (
                df[["unit_type", "Reserve_Non"]]
                .copy()
                .rename(columns={"Reserve_Non": "step_0"})
            )
            rev_df = (
                df[["unit_type", "UnitRevenue"]]
                .copy()
                .rename(columns={"UnitRevenue": "step_0"})
            )
            prof_df = (
                df[["unit_type", "UnitProfit"]]
                .copy()
                .rename(columns={"UnitProfit": "step_0"})
            )
        else:
            gen_df[f"step_{i}"] = df["Generation"].copy()
            rr_df[f"step_{i}"] = df["Reserve_Reg"].copy()
            sr_df[f"step_{i}"] = df["Reserve_Spin"].copy()
            nsr_df[f"step_{i}"] = df["Reserve_Non"].copy()
            rev_df[f"step_{i}"] = df["UnitRevenue"].copy()
            prof_df[f"step_{i}"] = df["UnitProfit"].copy()

    results = {
        "gen": gen_df,
        "rr": rr_df,
        "sr": sr_df,
        "nsr": nsr_df,
        "rev": rev_df,
        "prof": prof_df,
    }
    return results


def process_dispatch_data(ds_file_list, output_dir, scenario_name):
    """
    Process all sets of A-LEAF dispatch output data, in order to create
      dataframes of electricity price (LMP) data. Function filters the raw data
      to remove redundant LMP entries, and produces time-sorted and
      high-to-low sorted dataframes of results for all time steps.

    Returns:
      - unsorted_LMP_data: LMP data for each time step, sorted by timestamp
      - sorted_LMP_data: LMP data for each time step, sorted high to low
          (to create price duration curves)
    """

    num_files = len(ds_file_list)
    for i in range(num_files):
        # Read in the ALEAF dispatch output data for time-step i
        file_name = os.path.join(
            output_dir, f"{scenario_name}__dispatch_summary_OP__step_{i}.csv"
        )
        df = pd.read_csv(file_name)

        # Filter out duplicate price entries: the dispatch data df has one
        #   copy of each time period's price data for each unit type. I'm
        #   filtering for unique time periods by only selecting data for the
        #   "WIND" generator (though any generator type would work the same).
        df = df[df["UnitGroup"] == "WIND"].reset_index()

        # Create dataframes to hold unsorted and sorted LMP data separately
        unsorted_lmp = pd.DataFrame(df["LMP_dht"].copy())
        sorted_lmp = unsorted_lmp.copy().sort_values(
            by="LMP_dht", ascending=False, ignore_index=True
        )

        # Name the new column in each df appropriately
        if i == 0:
            unsorted_lmp_data = unsorted_lmp.rename(
                columns={"LMP_dht": "step_0"}
            )
            sorted_lmp_data = sorted_lmp.rename(columns={"LMP_dht": "step_0"})
        else:
            unsorted_lmp_data[f"step_{i}"] = unsorted_lmp
            sorted_lmp_data[f"step_{i}"] = sorted_lmp

    return unsorted_lmp_data, sorted_lmp_data


def plot_pdcs(sorted_lmp_data):
    """
    A function to plot price duration curves. Currently only plots
      step 0 results.
    """
    # Create log-log plots of the period-0 A-LEAF dispatch results
    x = np.arange(len(sorted_lmp_data)) + 1
    fig, ax = plt.subplots()
    ax.set_title("Log-log plot of period-0 price duration curve")
    ax.set_xlabel("log(period of the year)")
    ax.set_ylabel("log(price)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.plot(x, sorted_lmp_data["step_0"])
    fig.savefig("./abce_pd0_pdc.png")
