##########################################################################
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

def get_next_asset_id(db, first_asset_id):
    max_id = pd.read_sql("SELECT MAX(asset_id) FROM assets", db).iloc[0, 0]
    if max_id == None:
        next_id = first_asset_id
    else:
        next_id = max_id + 1

    return next_id


def process_outputs(settings, output_dir):
    ALEAF_scenario_name = settings["ALEAF_scenario_name"]
    file_types = ["dispatch_summary_OP", "expansion_result", "system_summary_OP", "system_tech_summary_OP"]
    file_lists = {}
    for ftype in file_types:
        file_lists[ftype] = glob.glob(os.path.join(output_dir, f"*{ftype}*"))
    expansion_results, expansion_results_mw = process_expansion_results(file_lists["expansion_result"], output_dir, ALEAF_scenario_name)
    system_summary_results = process_system_summary(file_lists["system_summary_OP"], file_lists["dispatch_summary_OP"], output_dir, ALEAF_scenario_name)
    system_tech_results = process_tech_summary(file_lists["system_tech_summary_OP"], output_dir, ALEAF_scenario_name)
    unsorted_lmp_data, sorted_lmp_data = process_dispatch_data(file_lists["dispatch_summary_OP"], output_dir, ALEAF_scenario_name)

    # Write results to xlsx
    writer = pd.ExcelWriter("abce_ppx_outputs.xlsx")
    expansion_results.to_excel(writer, sheet_name="exp_results", index=False)
    expansion_results_mw.to_excel(writer, sheet_name="exp_results_mw", index=False)
    system_summary_results.to_excel(writer, sheet_name="sys_summary_results", index=False)
    system_tech_results["gen"].to_excel(writer, sheet_name="tech_generation", index=False)
    system_tech_results["rr"].to_excel(writer, sheet_name="tech_reg", index=False)
    system_tech_results["sr"].to_excel(writer, sheet_name="tech_spin", index=False)
    system_tech_results["nsr"].to_excel(writer, sheet_name="tech_nonspin", index=False)
    system_tech_results["rev"].to_excel(writer, sheet_name="tech_revenue", index=False)
    system_tech_results["prof"].to_excel(writer, sheet_name="tech_profit", index=False)
    unsorted_lmp_data.to_excel(writer, sheet_name="unsorted_LMP", index=False)
    sorted_lmp_data.to_excel(writer, sheet_name="sorted_LMP", index=False)
    writer.save()

    # Plot PDCs
    plot_pdcs(sorted_lmp_data)


def process_expansion_results(exp_file_list, output_dir, ALEAF_scenario_name):
    num_files = len(exp_file_list)
    for i in range(num_files):
        file_name = os.path.join(output_dir, f"{ALEAF_scenario_name}__expansion_result__step_{i}.csv")
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
    unit_types = ["Wind", "Solar", "NGCC", "NGCT", "Advanced Nuclear"]
    exp_df["unit_id"] = unit_types
    exp_df = exp_df.rename(columns={"unit_id": "unit_type"})

    # Create separate dataframe in MW instead of # of units
    caps = [100, 100, 200, 50, 300]
    exp_df_mw = exp_df.copy()
    exp_df_mw = exp_df_mw.mul(caps, axis=0)
    exp_df_mw["unit_type"] = unit_types
 
    return exp_df, exp_df_mw


def process_system_summary(ss_file_list, ds_file_list, output_dir, ALEAF_scenario_name):
    num_files = len(ss_file_list)
    for i in range(num_files):
        file_name = os.path.join(output_dir, f"{ALEAF_scenario_name}__system_summary_OP__step_{i}.csv")
        df = pd.read_csv(file_name)
        # For my workflow, df should only have one line. Alert the user if otherwise
        if len(df) != 1:
            print("WARNING: It looks like you've run multiple ALEAF scenarios in the same run. Postprocessing outputs may be incorrectly translated.")

        # Give df a more helpful first-column key
        df = df.rename(columns={"Scenario": "step"})
        df.loc[0, "step"] = i

        # Correctly calculate weighted average electricity and AS prices, and 
        #   add them to the appropriate columns
        ds_file_name = os.path.join(output_dir, f"{ALEAF_scenario_name}__dispatch_summary_OP__step_{i}.csv")
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

        # Use the first file as a seed
        if i == 0:
            ss_df = df
        # Append each file's data to the main df
        else:
            ss_df = ss_df.append(df)

    return ss_df


def process_tech_summary(sts_file_list, output_dir, ALEAF_scenario_name):
    num_files = len(sts_file_list)
    for i in range(num_files):
        file_name = os.path.join(output_dir, f"{ALEAF_scenario_name}__system_tech_summary_OP__step_{i}.csv")
        df = pd.read_csv(file_name)

        # Clean up unit type representation
        df = df.rename(columns={"UnitGroup": "unit_type"})
        unit_types = ["Wind", "Solar", "NGCC", "NGCT", "Advanced Nuclear"]
        df["unit_type"] = unit_types

        # If step 0, use the file to set up DataFrames for all tracked items
        if i == 0:
            gen_df = df[["unit_type", "Generation"]].copy().rename(columns={"Generation": "step_0"})
            rr_df = df[["unit_type", "Reserve_Reg"]].copy().rename(columns={"Reserve_Reg": "step_0"})
            sr_df = df[["unit_type", "Reserve_Spin"]].copy().rename(columns={"Reserve_Spin": "step_0"})
            nsr_df = df[["unit_type", "Reserve_Non"]].copy().rename(columns={"Reserve_Non": "step_0"})
            rev_df = df[["unit_type", "UnitRevenue"]].copy().rename(columns={"UnitRevenue": "step_0"})
            prof_df = df[["unit_type", "UnitProfit"]].copy().rename(columns={"UnitProfit": "step_0"})
        else:
            gen_df[f"step_{i}"] = df["Generation"].copy()
            rr_df[f"step_{i}"] = df["Reserve_Reg"].copy()
            sr_df[f"step_{i}"] = df["Reserve_Spin"].copy()
            nsr_df[f"step_{i}"] = df["Reserve_Non"].copy()
            rev_df[f"step_{i}"] = df["UnitRevenue"].copy()
            prof_df[f"step_{i}"] = df["UnitProfit"].copy()

    results = {"gen": gen_df, "rr": rr_df, "sr": sr_df, "nsr": nsr_df, "rev": rev_df, "prof": prof_df}
    return results


def process_dispatch_data(ds_file_list, output_dir, ALEAF_scenario_name):
    num_files = len(ds_file_list)
    for i in range(num_files):
        file_name = os.path.join(output_dir, f"{ALEAF_scenario_name}__dispatch_summary_OP__step_{i}.csv")
        df = pd.read_csv(file_name)

        # Filter out duplicate price entries
        df = df[df["UnitGroup"] == "WIND"].reset_index()
        unsorted_lmp = pd.DataFrame(df["LMP_dht"].copy())
        sorted_lmp = unsorted_lmp.copy().sort_values(by="LMP_dht", ascending=False, ignore_index=True)

        if i == 0:
            unsorted_lmp_data = unsorted_lmp.rename(columns={"LMP_dht": "step_0"})
            sorted_lmp_data = sorted_lmp.rename(columns={"LMP_dht": "step_0"})
        else:
            unsorted_lmp_data[f"step_{i}"] = unsorted_lmp
            sorted_lmp_data[f"step_{i}"] = sorted_lmp

    return unsorted_lmp_data, sorted_lmp_data


def plot_pdcs(sorted_lmp_data):
    # Create log-log plots of the period-0 A-LEAF dispatch results
    x = np.arange(1, len(sorted_lmp_data)+1, 1)
    fig, ax = plt.subplots()
    ax.set_title("Log-log plot of period-0 price duration curve")
    ax.set_xlabel("log(period of the year")
    ax.set_ylabel("log(price)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.plot(x, sorted_lmp_data["step_0"])
    fig.savefig("./abce_pd0_pdc.png")



















 
