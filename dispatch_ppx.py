import pandas as pd
import logging

###############################################################################
# Constants
###############################################################################

dsp_file = "./ABCE_IL__dispatch_summary_OP__step_0.csv"

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
    "RMP_NS_dht": "RMP_nspin_dht"
}

sys_pf = pd.DataFrame.from_dict(
             {"Wind": [343],
              "Solar": [26],
              "NGCC": [300],
              "NGCT": [500],
              "AdvancedNuclear": [0],
              "ConventionalNuclear": [11],
              "Coal": [26]
             },
             orient="index",
             columns=["num_units"]
         )

unit_specs = pd.DataFrame.from_dict(
                 {"Wind": [0, 10, 0, 25, 100],
                  "Solar": [0, 8, 0, 25, 100],
                  "NGCC": [1.5, 3, 4, 0, 200],
                  "NGCT": [3, 5, 4, 0, 50],
                  "Coal": [4, 12, 10, 0, 500],
                  "ConventionalNuclear": [2, 7, 1, 15, 1000],
                  "AdvancedNuclear": [1, 5, 1, 15, 300]
                 },
                 orient="index",
                 columns=["VOM", "FOM", "FC", "policy_adj", "capacity"]
             )

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

    slim_df = dsp_df[list(value for key, value in std_names.items())].copy(deep=True)

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
    slim_df = (df[[
                  "unit_type",
                  "gen_idht",
                  "reg_idht",
                  "spin_idht",
                  "nspin_idht",
                  "gen_rev",
                  "reg_rev",
                  "spin_rev",
                  "nspin_rev"
               ]].copy(deep=True))
    agg_dsp_pivot = slim_df.groupby(["unit_type"]).sum()

    return agg_dsp_pivot


def join_unit_data(dsp_pivot, sys_pf, unit_specs):
    dsp_ext_pivot = dsp_pivot.join(sys_pf, how="left")
    dsp_ext_pivot = dsp_ext_pivot.join(unit_specs, how="left")

    return dsp_ext_pivot


def compute_perunit_results(agg_dsp_pivot):
    # Copy out columns which are already on a per-unit basis
    dsp_pivot_PU = (agg_dsp_pivot[[
                        "num_units",
                        "VOM",
                        "FOM",
                        "FC",
                        "policy_adj",
                        "capacity"
                    ]].copy(deep=True))

    # Get per-unit generation/AS quantities
    dsp_pivot_PU["gen_total"] = agg_dsp_pivot["gen_idht"] / agg_dsp_pivot["num_units"]
    dsp_pivot_PU["reg_total"] = agg_dsp_pivot["reg_idht"] / agg_dsp_pivot["num_units"]
    dsp_pivot_PU["spin_total"] = agg_dsp_pivot["spin_idht"] / agg_dsp_pivot["num_units"]
    dsp_pivot_PU["nspin_total"] = agg_dsp_pivot["nspin_idht"] / agg_dsp_pivot["num_units"]

    # Get per-unit revenue quantities
    dsp_pivot_PU["gen_rev"] = (
        agg_dsp_pivot["gen_rev"] / agg_dsp_pivot["num_units"]
    )
    dsp_pivot_PU["reg_rev"] = (
        agg_dsp_pivot["reg_rev"] / agg_dsp_pivot["num_units"]
    )
    dsp_pivot_PU["spin_rev"] = (
        agg_dsp_pivot["spin_rev"] / agg_dsp_pivot["num_units"]
    )
    dsp_pivot_PU["nspin_rev"] = (
        agg_dsp_pivot["nspin_rev"] / agg_dsp_pivot["num_units"]
    )
    dsp_pivot_PU["total_rev"] = (
        dsp_pivot_PU["gen_rev"] + dsp_pivot_PU["reg_rev"] 
        + dsp_pivot_PU["spin_rev"] + dsp_pivot_PU["nspin_rev"]
    )

    # Get per-unit costs and policy incentive/penalty impacts
    dsp_pivot_PU["var_costs"] = (
        dsp_pivot_PU["VOM"] + dsp_pivot_PU["capacity"] * dsp_pivot_PU["gen_total"]
    )
    dsp_pivot_PU["fixed_costs"] = (
        dsp_pivot_PU["FOM"] * dsp_pivot_PU["capacity"]
    )
    dsp_pivot_PU["total_policy_adj"] = (
        dsp_pivot_PU["policy_adj"] * dsp_pivot_PU["gen_total"]
    )

    # Get total per-unit costs
    dsp_pivot_PU["total_costs"] = (
        dsp_pivot_PU["var_costs"] + dsp_pivot_PU["fixed_costs"] 
    )
    dsp_pivot_PU["total_rev_w_policy"] = (
        dsp_pivot_PU["total_rev"] + dsp_pivot_PU["total_policy_adj"]
    )

    # Get total operating profit
    dsp_pivot_PU["op_profit"] = (
        dsp_pivot_PU["total_rev_w_policy"] - dsp_pivot_PU["total_costs"]
    )

    return dsp_pivot_PU


###############################################################################
# Running the script
###############################################################################

def postprocess_dispatch(dispatch_file):
    logging.info(f"Processing dispatch file at {dispatch_file}.")

    # Preprocess the file
    df = read_in_dispatch_file(dispatch_file)
    df = select_standard_columns(df)

    # Process results
    df = compute_timeseries_revenues(df)
    agg_dsp_pivot = pivot_dispatch_results(df)
    agg_dsp_pivot = join_unit_data(agg_dsp_pivot, sys_pf, unit_specs)
    dsp_pivot_PU = compute_perunit_results(agg_dsp_pivot)

    print(dsp_pivot_PU)

    logging.info("Dispatch file processed.")


if __name__ == "__main__":
    set_up_logging()

    # Process the dispatch file
    postprocess_dispatch(dsp_file)
