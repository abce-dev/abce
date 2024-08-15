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
import pandas as pd
import numpy as np
import yaml
import openpyxl
from pathlib import Path


def load_data(file_name):
    if Path(file_name).suffix in [".yml", ".yaml"]:
        file_contents = yaml.load(open(file_name, "r"), Loader=yaml.FullLoader)
    elif Path(file_name).suffix == ".csv":
        file_contents = pd.read_csv(file_name)

    return file_contents


def process_system_portfolio(db, current_pd):
    # Retrieve the list of currently-operational assets by type
    system_portfolio = pd.read_sql_query(
        f"SELECT unit_type, COUNT(unit_type) FROM assets "
        + f"WHERE completion_pd <= {current_pd} AND "
        + f"retirement_pd > {current_pd} "
        + f"GROUP BY unit_type",
        db,
    )

    # Rename the aggregation column to match the standard
    system_portfolio = system_portfolio.rename(
        columns={"COUNT(unit_type)": "num_units"}
    )

    return system_portfolio


def update_unit_specs_for_ALEAF(unit_specs_data, system_portfolio):
    # Separate out ATB_search_settings
    for unit_type, unit_type_data in unit_specs_data.items():
        if "ATB_search_settings" in unit_type_data.keys():
            for ATB_key, ATB_value in unit_type_data[
                "ATB_search_settings"
            ].items():
                unit_type_data[f"ATB_{ATB_key}"] = ATB_value

    # Convert unit_specs data from dict to DataFrame
    unit_specs = (
        pd.DataFrame.from_dict(unit_specs_data, orient="index")
        .reset_index()
        .rename(columns={"index": "unit_type"})
    )

    # Sort unit_specs alphabetically to allow definitive TechX id ordering
    unit_specs = unit_specs.sort_values("unit_type")

    # Create the UNIT_x columns
    unit_specs["UNITGROUP"] = unit_specs["unit_type"]
    unit_specs["UNIT_CATEGORY"] = unit_specs["unit_type"]

    # Create constant-value columns
    unit_specs["FOR"] = 0
    unit_specs["FCR"] = 0.05
    unit_specs["Charge_CAP"] = 0
    unit_specs["STOCAP"] = 0
    unit_specs["STOMIN"] = 0
    unit_specs["INVEST_FLAG"] = "FALSE"
    unit_specs["RET_FLAG"] = "FALSE"
    unit_specs["BATEFF"] = 0
    unit_specs["Integrality"] = "TRUE"
    unit_specs["Outages"] = "FALSE"
    unit_specs["bus_i"] = 1
    unit_specs["GenCo ID"] = 1
    unit_specs["MAXINVEST"] = 0
    unit_specs["MININVEST"] = 0
    unit_specs["MINRET"] = 0
    unit_specs["ATB_Setting_ID"] = "ATB_ID_1"

    # Create the computed Emission column
    unit_specs["Emission"] = (
        unit_specs["emissions_per_MMBTU"] * unit_specs["heat_rate"] / 1000
    )

    # Take care of the row-by-row operations
    unit_specs["Tech_ID"] = ""
    unit_specs["Commitment"] = ""
    unit_specs["GEN UID"] = ""
    tech_id = 1
    unit_specs = unit_specs.set_index("unit_type")
    for unit_type in unit_specs.itertuples():
        # Set Commitment to the inverse of is_VRE
        unit_specs.loc[unit_type.Index, "Commitment"] = not unit_type.is_VRE

        # Set the Tech_ID value in alphabetical order
        unit_specs.loc[unit_type.Index, "Tech_ID"] = f"Tech{tech_id}"
        tech_id += 1

        # Set the GEN UID field
        unit_specs.loc[unit_type.Index, "GEN UID"] = f"{unit_type.Index}_1"

    # Pivot in num_units data from the system portfolio for convenience
    unit_specs = unit_specs.reset_index().rename(columns={"index": "unit_type"})
    unit_specs = unit_specs.merge(system_portfolio, how="outer", on="unit_type")

    # Fill any NaNs with zeros
    unit_specs = unit_specs.fillna(0)

    return unit_specs


def create_ALEAF_unit_dataframes(unit_specs):
    gen_technology_cols = {
        "Tech_ID": "Tech_ID",
        "UNITGROUP": "UNITGROUP",
        "UNIT_CATEGORY": "UNIT_CATEGORY",
        "unit_type": "UNIT_TYPE",
        "fuel_type": "FUEL",
        "capacity": "CAP",
        "max_PL": "PMAX",
        "min_PL": "PMIN",
        "FOR": "FOR",
        "overnight_capital_cost": "CAPEX",
        "FCR": "FCR",
        "retirement_cost": "RETC",
        "FOM": "FOM",
        "VOM": "VOM",
        "FC_per_MMBTU": "FC",
        "no_load_cost": "NLC",
        "start_up_cost": "SUC",
        "shut_down_cost": "SDC",
        "heat_rate": "HR",
        "ramp_up_limit": "RUL",
        "ramp_down_limit": "RDL",
        "max_regulation": "MAXR",
        "max_spinning_reserve": "MAXSR",
        "max_nonspinning_reserve": "MAXNSR",
        "capacity_factor": "CAPCRED",
        "emissions_per_MMBTU": "EMSFAC",
        "Charge_CAP": "Charge_CAP",
        "STOCAP": "STOCAP",
        "STOMIN": "STOMIN",
        "INVEST_FLAG": "INVEST_FLAG",
        "RET_FLAG": "RET_FLAG",
        "is_VRE": "VRE_Flag",
        "BATEFF": "BATEFF",
        "Commitment": "Commitment",
        "Integrality": "Integrality",
        "Emission": "Emission",
        "Outages": "Outages",
    }

    gen_technology = (
        unit_specs[list(gen_technology_cols.keys())]
        .copy()
        .rename(columns=gen_technology_cols)
    )

    gen_cols = {
        "GEN UID": "GEN UID",
        "bus_i": "bus_i",
        "Tech_ID": "Tech_ID",
        "UNITGROUP": "UNITGROUP",
        "UNIT_CATEGORY": "UNIT_CATEGORY",
        "unit_type": "UNIT_TYPE",
        "fuel_type": "FUEL",
        "num_units": "EXUNITS",
        "capacity": "CAP",
        "MAXINVEST": "MAXINVEST",
        "MININVEST": "MININVEST",
        "MINRET": "MINRET",
    }

    gen = unit_specs[list(gen_cols.keys())].copy().rename(columns=gen_cols)

    ATB_settings_cols = {
        "ATB_Setting_ID": "ATB_Setting_ID",
        "Tech_ID": "Tech_ID",
        "UNITGROUP": "UNITGROUP",
        "UNIT_CATEGORY": "UNIT_CATEGORY",
        "unit_type": "UNIT_TYPE",
        "fuel_type": "FUEL",
        "ATB_Tech": "Tech",
        "ATB_TechDetail": "TechDetail",
        "ATB_Case": "Case",
        "ATB_CRP": "CRP",
        "ATB_Scenario": "Scenario",
        "ATB_Year": "Year",
        "ATB_ATB_year": "ATB Year",
    }

    ATB_settings = (
        unit_specs[list(ATB_settings_cols.keys())]
        .copy()
        .rename(columns=ATB_settings_cols)
    )

    return gen_technology, gen, ATB_settings


def create_ALEAF_Master_file(ALEAF_data, settings):
    # Dictionary of all tabs and their metadata
    # In final implementation, will be replaced by the ALEAF standard
    #   data schema
    tabs_to_create = {
        "ALEAF Master Setup": {
            "ABCE_tab_name": "ALEAF_Master_setup",
            "data": None,
        },
        "CPLEX Setting": {"ABCE_tab_name": "CPLEX_settings", "data": None},
        "GLPK Setting": {"ABCE_tab_name": "GLPK_settings", "data": None},
        "CBC Setting": {"ABCE_tab_name": "CBC_settings", "data": None},
        "Gurobi Setting": {"ABCE_tab_name": "Gurobi_settings", "data": None},
        "HiGHS Setting": {"ABCE_tab_name": "HiGHS_settings", "data": None},
    }

    # Pull the ALEAF settings data into dictionaries
    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        tabs_to_create[ALEAF_tab_name]["data"] = ALEAF_data["ALEAF_Master"][
            tab_data["ABCE_tab_name"]
        ]

    # Pull the appropriate solver name from the settings file
    tabs_to_create["ALEAF Master Setup"]["solver_name"] = settings[
        "simulation"
    ]["solver"]

    # Finalize the <solver> Setting tab data
    for solver_tab, tab_data in tabs_to_create.items():
        if solver_tab != "ALEAF Master Setup":
            # Set up metadata about solver settings
            invalid_items = [
                "solver_setting_list",
                "num_solver_setting",
                "solver_direct_mode_flag",
            ]
            solver_setting_list = ", ".join(
                [
                    parameter
                    for parameter in tabs_to_create[solver_tab]["data"].keys()
                    if parameter not in invalid_items
                ]
            )

            # Set up solver_direct_mode_flag: TRUE if CPLEX, FALSE otherwise
            mode_flag = "false"
            if "CPLEX" in solver_tab:
                mode_flag = "true"

            # Create the dictionary of extra rows for all solver tabs
            solver_extra_items = {
                "solver_direct_mode_flag": mode_flag,
                "num_solver_setting": len(tabs_to_create[solver_tab]["data"]),
                "solver_setting_list": solver_setting_list,
            }

            for key, value in solver_extra_items.items():
                if key not in tab_data["data"].keys():
                    tab_data["data"].update({key: value})

    # Construct the path to which this file should be written
    output_path = (
        Path(os.environ["ALEAF_DIR"])
        / "setting"
        / settings["ALEAF"]["ALEAF_master_settings_file"]
    )

    # Write this file to the destination
    write_workbook_and_close("ALEAF_Master", tabs_to_create, output_path)


def create_ALEAF_Master_LC_GEP_file(
    ALEAF_data, gen_technology, ATB_settings, settings
):
    tabs_to_create = {
        "LC_GEP Setting": {"ABCE_tab_name": "LC_GEP_settings", "data": None},
        "Planning Design": {"ABCE_tab_name": "planning_design", "data": None},
        "Simulation Setting": {
            "ABCE_tab_name": "simulation_settings",
            "data": None,
        },
        "Simulation Configuration": {
            "ABCE_tab_name": "scenario_settings",
            "data": None,
            "orient": "horizontal",
        },
        "Gen Technology": {
            "ABCE_tab_name": "unit_specs",
            "data": gen_technology,
        },
        "File Path": {
            "ABCE_tab_name": "ALEAF_relative_file_paths",
            "data": None,
        },
        "Scenario Reduction Setting": {
            "ABCE_tab_name": "scenario_reduction_settings",
            "data": None,
        },
        "ATB Setting": {
            "ABCE_tab_name": "ATB_search_settings",
            "data": ATB_settings,
        },
    }

    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        if not isinstance(tab_data["data"], pd.DataFrame):
            tabs_to_create[ALEAF_tab_name]["data"] = ALEAF_data[
                "ALEAF_Master_LC_GEP"
            ][tab_data["ABCE_tab_name"]]

    # Finalize tab data
    # Finalize "Planning Design" tab
    # Add extra items
    pd_data = tabs_to_create["Planning Design"]["data"]
    pd_extra_items = {
        "targetyear_value": (
            pd_data["final_year_value"] - pd_data["current_year_value"]
        ),
        "load_increase_rate_value": 1,
        "num_simulation_per_stage_value": (
            (pd_data["final_year_value"] - pd_data["current_year_value"])
            / pd_data["numstages_value"]
        ),
    }
    tabs_to_create["Planning Design"]["data"].update(pd_extra_items)

    # Rename items
    items_to_rename = {
        "planning_reserve_margin": "planning_reserve_margin_value"
    }
    for key, value in items_to_rename.items():
        if value not in tabs_to_create["Planning Design"]["data"].keys():
            tabs_to_create["Planning Design"]["data"][value] = tabs_to_create[
                "Planning Design"
            ]["data"].pop(key)

    # Finalize "Simulation Setting" tab
    # Add extra items
    ss_extra_items = {
        "test_system_file_name": f"ALEAF_{tabs_to_create['Simulation Setting']['data']['test_system_name']}"
    }
    tabs_to_create["Simulation Setting"]["data"].update(ss_extra_items)

    # Rename items
    items_to_rename = {
        "capex_projection_flag": "capax_projection_flag",
        "network_reduction_flag": "Network_reduction_flag",
    }
    for key, value in items_to_rename.items():
        if value not in tabs_to_create["Simulation Setting"]["data"].keys():
            tabs_to_create["Simulation Setting"]["data"][
                value
            ] = tabs_to_create["Simulation Setting"]["data"].pop(key)

    # Finalize "Simulation Configuration" tab
    # Update existing data items
    # Update peak demand
    tabs_to_create["Simulation Configuration"]["data"][
        "peak_demand"
    ] = settings["scenario"]["peak_demand"]

    # Update policies
    if "policies" in settings["scenario"].keys():
        for policy, policy_data in settings["scenario"]["policies"].items():
            if policy_data["enabled"]:
                qty = policy_data["qty"]
                if (
                    "PTC" in policy
                    and "unit_type" in policy_data["eligibility"].keys()
                ):
                    eligible_types = policy_data["eligibility"]["unit_type"]
                    if any("wind" in unit_type for unit_type in eligible_types):
                        tabs_to_create["Simulation Configuration"]["data"][
                            "wind_PTC"
                        ] = qty
                    if any(
                        "solar" in unit_type for unit_type in eligible_types
                    ):
                        tabs_to_create["Simulation Configuration"]["data"][
                            "solar_PTC"
                        ] = qty
                    if any(
                        "nuclear" in unit_type for unit_type in eligible_types
                    ):
                        tabs_to_create["Simulation Configuration"]["data"][
                            "nuclear_PTC"
                        ] = qty
                if (
                    "ITC" in policy
                    and "unit_type" in policy_data["eligibility"].keys()
                ):
                    eligible_types = policy_data["eligibility"]["unit_type"]
                    if any("wind" in unit_type for unit_type in eligible_types):
                        tabs_to_create["Simulation Configuration"]["data"][
                            "wind_ITC"
                        ] = qty
                    if any(
                        "solar" in unit_type for unit_type in eligible_types
                    ):
                        tabs_to_create["Simulation Configuration"]["data"][
                            "solar_ITC"
                        ] = qty
                    if any(
                        "nuclear" in unit_type for unit_type in eligible_types
                    ):
                        tabs_to_create["Simulation Configuration"]["data"][
                            "nuclear_ITC"
                        ] = qty
                if "CTAX" in policy:
                    tabs_to_create["Simulation Configuration"]["data"][
                        "carbon_tax"
                    ] = qty

    # Add extra items
    sc_extra_items = {
        "planning_reserve_margin_value": tabs_to_create["Planning Design"][
            "data"
        ]["planning_reserve_margin_value"],
        "load_increase_rate_value": tabs_to_create["Planning Design"]["data"][
            "load_increase_rate_value"
        ],
    }
    tabs_to_create["Simulation Configuration"]["data"].update(sc_extra_items)

    # Rename items
    items_to_rename = {
        "scenario_name": "Scenario",
        "peak_demand": "PD",
        "carbon_tax": "CTAX",
        "RPS_percentage": "RPS",
        "wind_PTC": "PTC_W",
        "solar_PTC": "PTC_S",
        "nuclear_PTC": "PTC_N",
        "wind_ITC": "ITC_W",
        "solar_ITC": "ITC_S",
        "nuclear_ITC": "ITC_N",
    }
    for key, value in items_to_rename.items():
        if (
            value
            not in tabs_to_create["Simulation Configuration"]["data"].keys()
        ):
            tabs_to_create["Simulation Configuration"]["data"][
                value
            ] = tabs_to_create["Simulation Configuration"]["data"].pop(key)

    # Finalize "Scenario Reduction Setting" tab
    srs_data = tabs_to_create["Scenario Reduction Setting"]["data"]
    data_sets = [
        srs_data["input_type_load_shape_flag"],
        srs_data["input_type_load_MWh_flag"],
        srs_data["input_type_wind_shape_flag"],
        srs_data["input_type_wind_MWh_flag"],
        srs_data["input_type_solar_shape_flag"],
        srs_data["input_type_solar_MWh_flag"],
        srs_data["input_type_net_load_MWh_flag"],
    ]
    num_data_sets = sum(1 for item in data_sets if item == "TRUE")
    srs_extra_items = {"num_data_set": num_data_sets}
    tabs_to_create["Scenario Reduction Setting"]["data"].update(srs_extra_items)

    # Construct the path to which this file should be written
    output_path = (
        Path(os.environ["ALEAF_DIR"])
        / "setting"
        / settings["ALEAF"]["ALEAF_model_settings_file"]
    )

    # Write this file to the destination
    write_workbook_and_close("ALEAF_Master_LC_GEP", tabs_to_create, output_path)


def create_ALEAF_portfolio_file(ALEAF_data, gen, settings):
    tabs_to_create = {
        "case setting": {"ABCE_tab_name": "grid_settings", "data": None},
        "gen": {"ABCE_tab_name": "system_portfolio", "data": gen},
        "bus": {"ABCE_tab_name": "buses", "data": None, "orient": "horizontal"},
        "branch": {
            "ABCE_tab_name": "branch",
            "data": None,
            "orient": "horizontal",
        },
        "sub_area": {
            "ABCE_tab_name": "sub_area",
            "data": None,
            "orient": "horizontal",
        },
        "sub_area_mapping": {
            "ABCE_tab_name": "sub_area_mapping",
            "data": None,
            "orient": "horizontal",
        },
    }

    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        if not isinstance(tab_data["data"], pd.DataFrame):
            tabs_to_create[ALEAF_tab_name]["data"] = ALEAF_data[
                "ALEAF_portfolio"
            ][tab_data["ABCE_tab_name"]]

    # Construct the path to which this file should be written
    output_path = (
        Path(os.environ["ALEAF_DIR"])
        / "data"
        / settings["ALEAF"]["ALEAF_model_type"]
        / settings["ALEAF"]["ALEAF_region"]
        / settings["ALEAF"]["ALEAF_portfolio_file"]
    )

    # Write this file to the destination
    write_workbook_and_close("ALEAF_ERCOT", tabs_to_create, output_path)


def write_workbook_and_close(base_filename, tabs_to_create, output_file_path):
    # Load all tab data into pandas dataframes
    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        # If the data is already in a dataframe, no need to convert it
        if not isinstance(tab_data["data"], pd.DataFrame):
            orient = "index"
            if "orient" in tab_data.keys():
                if tab_data["orient"] == "horizontal":
                    df = pd.DataFrame.from_dict([tab_data["data"]])
                elif tab_data["orient"] == "multiline_horizontal":
                    df = pd.DataFrame.from_dict(
                        tab_data["data"], orient="index"
                    )
            else:
                df = (
                    pd.DataFrame.from_dict(
                        tab_data["data"], orient="index", columns=["Value"]
                    )
                    .reset_index()
                    .rename(columns={"index": "Setting"})
                )

            tabs_to_create[ALEAF_tab_name]["data"] = df

    # Create an ExcelWriter object to contain all tabs and save file
    writer_object = pd.ExcelWriter(output_file_path, engine="openpyxl")

    # Write all tabs to file
    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        tab_data["data"].to_excel(
            writer_object, sheet_name=ALEAF_tab_name, index=False
        )

    # Write all changes and close the file writer
    writer_object.close()


def set_unit_type_policy_adjustment(unit_type, unit_type_data, settings):
    # Initialize all units with zero policy adjustment
    carbon_tax_per_MWh = 0
    tax_credits_per_MWh = 0
    tax_credits_per_MW = 0

    if "policies" in settings["scenario"].keys():
        # If some policies are specified, determine their effect
        for policy, policy_specs in settings["scenario"]["policies"].items():
            if policy_specs["enabled"]:
                is_eligible = False
                if "eligibility" not in policy_specs.keys():
                    # If there are no eligibility criteria specified, the
                    #   policy applies to all units
                    is_eligible = True
                elif (
                    "unit_type" in policy_specs["eligibility"].keys()
                    and unit_type in policy_specs["eligibility"]["unit_type"]
                ):
                    # Check to see if the unit's type is included in a list of
                    #   eligible unit types
                    is_eligible = True
                else:
                    # Check for any other eligibility criteria
                    for criterion, value in policy_specs["eligibility"].items():
                        if (
                            criterion != "unit_type"
                            and unit_type_data[criterion] == value
                        ):
                            is_eligible = True

                if is_eligible:
                    if "CTAX" in policy:
                        carbon_tax_per_MWh -= (
                            unit_type_data["heat_rate"]
                            * unit_type_data["emissions_per_MMBTU"]
                            * policy_specs["qty"]
                        )
                    elif "PTC" in policy:
                        tax_credits_per_MWh += policy_specs["qty"]
                    elif "ITC" in policy:
                        tax_credits_per_MW += policy_specs["qty"]

    policy_results = {
        "carbon_tax_per_MWh": carbon_tax_per_MWh,
        "tax_credits_per_MWh": tax_credits_per_MWh,
        "tax_credits_per_MW": tax_credits_per_MW,
    }

    return carbon_tax_per_MWh, tax_credits_per_MWh, tax_credits_per_MW,



def compute_unit_specs_cols(unit_specs, settings):
    for unit_type, unit_type_data in unit_specs.items():
        # Ensure all units have fuel cost elements
        if not unit_type_data["uses_fuel"]:
            unit_type_data["FC_per_MMBTU"] = 0
            unit_type_data["fuel_type"] = "none"

        # Add fuel cost per MWh column
        unit_type_data["FC_per_MWh"] = (
            unit_type_data["FC_per_MMBTU"] * unit_type_data["heat_rate"]
        )

        # Add policy adjustment per MWh column
        unit_type_data["carbon_tax_per_MWh"], unit_type_data["tax_credits_per_MWh"], unit_type_data["tax_credits_per_MW"] = set_unit_type_policy_adjustment(
            unit_type, unit_type_data, settings
        )

        # Add C2N-related factors, with a default value of 0
        C2N_vals = ["cpp_ret_lead", "num_cpp_rets", "rev_head_start"]
        for val in C2N_vals:
            if val not in unit_type_data.keys():
                unit_type_data[val] = 0

    return unit_specs


def initialize_unit_specs(settings, args):
    unit_specs = load_data(
        Path(args.inputs_path)
        / settings["file_paths"]["unit_specs_data_file"]
    )

    unit_specs = compute_unit_specs_cols(unit_specs, settings)

    return unit_specs


def update_ALEAF_data(ALEAF_data, settings):
    # Update ALEAF_data with settings data
    ALEAF_data["ALEAF_Master_LC_GEP"]["scenario_settings"][
        "scenario_name"
    ] = settings["simulation"]["scenario_name"]

    return ALEAF_data


def create_ALEAF_files(settings, ALEAF_data, unit_specs_data, db, current_pd):
    # Process the system portfolio
    system_portfolio = process_system_portfolio(db, current_pd)

    # Process the dictionary unit_specs into the A-LEAF-style format
    unit_specs = update_unit_specs_for_ALEAF(unit_specs_data, system_portfolio)

    # Process the unit_specs data into A-LEAF-ready dataframes
    gen_technology, gen, ATB_settings = create_ALEAF_unit_dataframes(unit_specs)

    # Update ALEAF data
    ALEAF_data = update_ALEAF_data(ALEAF_data, settings)

    # Create the ALEAF_Master.xlsx file
    create_ALEAF_Master_file(ALEAF_data, settings)

    # Create the ALEAF_Master_LC_GEP.xlsx file
    create_ALEAF_Master_LC_GEP_file(
        ALEAF_data, gen_technology, ATB_settings, settings
    )

    # Create the ALEAF_portfolio.xlsx file
    create_ALEAF_portfolio_file(ALEAF_data, gen, settings)
