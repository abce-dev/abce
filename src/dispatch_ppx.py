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


import pandas as pd
import logging

###############################################################################
# Constants
###############################################################################

std_names = {
    "day": "day",
    "hour": "hour",
    "time": "time",
    "Unit_Type": "unit_type",
    "c_idh": "c_idh",
    "g_idht": "gen_idht",
    "r_r_idht": "reg_idht",
    "r_s_idht": "spin_idht",
    "r_ns_idht": "nspin_idht",
    "LMP_dht": "LMP_dht",
    "RMP_R_dht": "RMP_reg_dht",
    "RMP_S_dht": "RMP_spin_dht",
    "RMP_NS_dht": "RMP_nspin_dht",
}

###############################################################################
# Housekeeping
###############################################################################


def set_up_logging():
    logging.basicConfig(level=logging.INFO)


###############################################################################
# Preprocessing
###############################################################################


def read_in_dispatch_file(dispatch_file):
    df = pd.read_csv(dispatch_file)

    return df


def select_standard_columns(dsp_df):
    dsp_df = dsp_df.rename(columns=std_names)

    slim_df = dsp_df[list(value for key, value in std_names.items())].copy(
        deep=True
    )

    return slim_df


###############################################################################
# Processing the dispatch data
###############################################################################


def compute_timeseries_revenues(df):
    df["gen_rev"] = df["gen_idht"] * df["LMP_dht"]
    df["reg_rev"] = df["reg_idht"] * df["RMP_reg_dht"]
    df["spin_rev"] = df["spin_idht"] * df["RMP_spin_dht"]
    df["nspin_rev"] = df["nspin_idht"] * df["RMP_nspin_dht"]

    return df


def pivot_dispatch_results(df):
    slim_df = df[
        [
            "unit_type",
            "gen_idht",
            "reg_idht",
            "spin_idht",
            "nspin_idht",
            "gen_rev",
            "reg_rev",
            "spin_rev",
            "nspin_rev",
        ]
    ].copy(deep=True)
    agg_dsp_pivot = slim_df.groupby(["unit_type"]).sum()

    return agg_dsp_pivot


def join_unit_data(dsp_pivot, system_portfolio, unit_specs):
    # Join system_portfolio into the dispatch results
    dsp_ext_pivot = dsp_pivot.join(system_portfolio, how="inner")

    # Convert unit_specs into a dataframe for this join
    unit_specs_df = (
        pd.DataFrame.from_dict(unit_specs, orient="index")
        .reset_index()
        .rename(columns={"index": "unit_type"})
    )

    # Join unit_specs into the dispatch results, for convenience
    dsp_ext_pivot = dsp_ext_pivot.join(
        unit_specs_df.set_index("unit_type"), how="inner"
    )

    return dsp_ext_pivot


def compute_per_unit_results(agg_dsp_pivot):
    # Copy out columns which are already on a per-unit basis
    dsp_pivot_PU = agg_dsp_pivot[
        [
            "num_units",
            "VOM",
            "FOM",
            "FC_per_MWh",
            "carbon_tax_per_MWh",
            "tax_credits_per_MWh",
            "capacity",
        ]
    ].copy(deep=True)

    # Compute per-unit generation/AS quantities and revenues
    services = ["gen", "reg", "spin", "nspin"]
    for service in services:
        dsp_pivot_PU[f"{service}_total"] = (
            agg_dsp_pivot[f"{service}_idht"] / agg_dsp_pivot["num_units"]
        )
        dsp_pivot_PU[f"{service}_rev"] = (
            agg_dsp_pivot[f"{service}_rev"] / agg_dsp_pivot["num_units"]
        )

    dsp_pivot_PU["generation"] = dsp_pivot_PU["gen_total"]

    dsp_pivot_PU["revenue"] = (
        dsp_pivot_PU["gen_rev"]
        + dsp_pivot_PU["reg_rev"]
        + dsp_pivot_PU["spin_rev"]
        + dsp_pivot_PU["nspin_rev"]
    )

    # Get per-unit costs and policy incentive/penalty impacts
    dsp_pivot_PU["var_costs"] = dsp_pivot_PU["VOM"] * dsp_pivot_PU["generation"]
    dsp_pivot_PU["fuel_cost"] = (
        dsp_pivot_PU["FC_per_MWh"] * dsp_pivot_PU["generation"]
    )
    dsp_pivot_PU["FOM"] = dsp_pivot_PU["FOM"] * dsp_pivot_PU["capacity"] * 1000
    dsp_pivot_PU["carbon_tax"] = (
        dsp_pivot_PU["carbon_tax_per_MWh"] * dsp_pivot_PU["generation"]
    )

    dsp_pivot_PU["tax_credits"] = dsp_pivot_PU["tax_credits_per_MWh"] * dsp_pivot_PU["generation"]

    return dsp_pivot_PU


def downselect_dispatch_econ_results(dsp_pivot_PU):
    final_dsp_results = dsp_pivot_PU[
        [
            "generation",
            "reg_total",
            "spin_total",
            "nspin_total",
            "gen_rev",
            "reg_rev",
            "spin_rev",
            "nspin_rev",
            "revenue",
            "VOM",
            "fuel_cost",
            "FOM",
            "carbon_tax",
            "tax_credits",
        ]
    ].copy(deep=True)

    return final_dsp_results


def resolve_nans(final_dsp_results):
    final_dsp_results = final_dsp_results.fillna(0)

    return final_dsp_results


###############################################################################
# Running the script
###############################################################################


def postprocess_dispatch(dispatch_file, system_portfolio, unit_specs):
    logging.info(f"Processing dispatch file at {dispatch_file}.")

    # Preprocess the file
    df = read_in_dispatch_file(dispatch_file)
    df = select_standard_columns(df)

    # Process results
    df = compute_timeseries_revenues(df)
    agg_dsp_pivot = pivot_dispatch_results(df)
    agg_dsp_pivot = join_unit_data(agg_dsp_pivot, system_portfolio, unit_specs)
    dsp_pivot_PU = compute_per_unit_results(agg_dsp_pivot)
    final_dsp_results = downselect_dispatch_econ_results(dsp_pivot_PU)
    final_dsp_results = resolve_nans(final_dsp_results)

    logging.info("Dispatch file processed.")

    return final_dsp_results
