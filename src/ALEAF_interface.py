# A module allowing interactions between ABCE and A-LEAF

import os
import pandas as pd
import openpyxl
import sqlite3
import yaml
from pathlib import Path

# import local modules
from . import seed_creator as sc


def prepare_xlsx_data(ref_data_file, destination_file):
    book = openpyxl.load_workbook(ref_data_file)
    writer = pd.ExcelWriter(destination_file, engine="openpyxl")
    writer.book = book
    writer.sheets = dict((sheet.title, sheet) for sheet in book.worksheets)

    return book, writer


def update_ALEAF_system_portfolio(ALEAF_portfolio_ref, ALEAF_portfolio_remote, db, current_pd):
    # Get a list of all possible unit types, by pulling from the unit_specs table
    all_unit_types = pd.read_sql_query(f"SELECT unit_type FROM unit_specs GROUP BY unit_type", db)
    all_unit_types = all_unit_types["unit_type"].to_list()

    # Get the updated list of currently-operating units by unit type
    unit_type_count = pd.read_sql_query(f"SELECT unit_type, COUNT(unit_type) FROM assets WHERE completion_pd <= {current_pd} AND retirement_pd > {current_pd} AND cancellation_pd > {current_pd} GROUP BY unit_type", db)
    unit_type_count = unit_type_count.set_index("unit_type")

    db.commit()

    # Retrieve and organize the A-LEAF system portfolio
    book, writer = prepare_xlsx_data(ALEAF_portfolio_ref, ALEAF_portfolio_remote)
    df = organize_ALEAF_portfolio(writer)

    # Update unit type numbers
    for unit_type in all_unit_types:
        if unit_type in unit_type_count.index:
            new_count = unit_type_count.loc[unit_type, "COUNT(unit_type)"]
        else:
            new_count = 0
        df.loc[(df["bus_i"] == 1) & (df["UNIT_TYPE"] == unit_type), "EXUNITS"] = new_count

    df.to_excel(writer, sheet_name="gen", header=True, index=False)
    writer.save()


def organize_ALEAF_portfolio(writer):
    df = pd.DataFrame(writer.sheets["gen"].values)
    df.columns = df.iloc[0]
    df = df.drop(df.index[0]).reset_index(drop=True)
    df["EXUNITS"] = df["EXUNITS"].astype("int64")
    df["UNIT_TYPE"] = df["UNIT_TYPE"].astype("string")

    return df


def update_ALEAF_model_settings(ALEAF_model_settings_ref, ALEAF_model_settings_remote, db, settings, period=0):
    # Read in ALEAF simulation settings from ALEAF_Master_LC_GEP.xlsx-style file
    book, writer = prepare_xlsx_data(ALEAF_model_settings_ref, ALEAF_model_settings_remote)
    sim_settings = pd.DataFrame(writer.sheets["Simulation Configuration"].values)
    sim_settings.columns = sim_settings.iloc[0]
    sim_settings = sim_settings.drop(sim_settings.index[0]).reset_index().drop("index", axis=1)

    # Read in demand data from DB
    demand = pd.read_sql_query(f"SELECT demand FROM demand WHERE period = {period}", db)

    # If this is the first step of the run, update top-level settings
    if period == 0:
        # Update ALEAF scenario name
        sim_settings.loc[0, "Scenario"] = settings["simulation"]["ALEAF_scenario_name"]


    # Update periodic data
    # Update peak demand (PD, MW)
    sim_settings.loc[sim_settings["Scenario"] == settings["simulation"]["ALEAF_scenario_name"], "PD"] = demand.iloc[0, 0]

    # Write the updated sheet to file
    sim_settings.to_excel(writer, sheet_name = "Simulation Configuration", header=True, index=False)

    # Update key columns from unit_specs table
    gen_tech = pd.DataFrame(writer.sheets["Gen Technology"].values)
    gen_tech.columns = gen_tech.iloc[0]
    gen_tech = gen_tech.drop(gen_tech.index[0]).reset_index().drop("index", axis=1)

    ALEAF_header_converter = {"UNIT_TYPE": "unit_type",
                              "FC": "FC_per_MWh",
                              "VOM": "VOM"
                             }

    cols_to_update = [key for key in ALEAF_header_converter.keys() if key != "UNIT_TYPE"]

    # Get the final unit_specs data from the database
    unit_specs = pd.read_sql_query("SELECT * FROM unit_specs", db)

    # Update all data in cols_to_update
    for i in range(len(gen_tech)):
        unit_type = gen_tech.loc[i, "UNIT_TYPE"]
        for col in cols_to_update:
            gen_tech.loc[i, col] = unit_specs.loc[unit_specs["unit_type"] == unit_type, ALEAF_header_converter[col]].values[0]

    gen_tech.to_excel(writer, sheet_name = "Gen Technology", header=True, index=False)

    writer.save()


def update_ALEAF_policy_settings(ALEAF_model_settings_ref, ALEAF_model_settings_remote, policies_dict, unit_specs):
    valid_CTAX_names = ["CTAX", "ctax", "carbon_tax", "carbontax"]
    valid_PTC_names = ["PTC", "ptc", "production_tax_credit", "productiontaxcredit"]
    valid_wind_names = ["Wind", "wind"]
    valid_solar_names = ["Solar", "PV", "Solar PV", "solar", "pv", "solar pv", "Solar_PV", "solar_pv"]
    valid_nuclear_names = ["Nuclear", "nuclear", "AdvancedNuclear", "advancednuclear", "Advanced_Nuclear", "advanced_nuclear"]

    policy_cols = ["CTAX", "RPS", "PTC_W", "PTC_S", "ITC_S", "ITC_W", "CAPPMT"]

    # Read in ALEAF simulation settings from ALEAF_Master_LC_GEP.xlsx-style file
    book, writer = prepare_xlsx_data(ALEAF_model_settings_ref, ALEAF_model_settings_remote)
    sim_settings = pd.DataFrame(writer.sheets["Simulation Configuration"].values)
    sim_settings.columns = sim_settings.iloc[0]
    sim_settings = sim_settings.drop(sim_settings.index[0]).reset_index().drop("index", axis=1)

    # Zero out all policy columns
    for col in policy_cols:
        sim_settings[col] = 0
    for key, val in policies_dict.items():
        if (key in valid_CTAX_names) and (val["enabled"]):
            sim_settings["CTAX"] = val["qty"]
        elif (key in valid_PTC_names) and (val["enabled"]):
            for unit_type in val["eligible"]:
                if unit_type in valid_wind_names:
                    sim_settings["PTC_W"] = val["qty"]
                elif unit_type in valid_solar_names:
                    sim_settings["PTC_S"] = val["qty"]
                elif unit_type in valid_nuclear_names:
                    sim_settings["PTC_Nuc"] = val["qty"]

    # Write the updated sheet to file
    sim_settings.to_excel(writer, sheet_name = "Simulation Configuration", header=True, index=False)
    writer.save()


# New functions start here


def read_ALEAF_schema():
    ALEAF_schema = yaml.load(
                       open("./inputs/ALEAF_settings_schema.yml", "r"),
                       Loader=yaml.FullLoader
                   )

    return ALEAF_schema



def read_ABCE_data():
    # To replace with ABCE data setup integration
    unit_specs = yaml.load(
                     open("./inputs/unit_specs.yml", "r"),
                     Loader=yaml.FullLoader
                 )
    settings = yaml.load(open("./settings.yml", "r"), Loader=yaml.FullLoader)

    # Package all ABCE data into a dictionary
    ABCE_data = {
        "settings": settings,
        "unit_specs": unit_specs
    }

    return ABCE_data


def initialize_ALEAF_data(ALEAF_schema, ABCE_data):
    # Set up the top-level container dictionary for ALEAF data
    ALEAF_data_dict = {}

    # Iterate through all files
    for file_name, file_data in ALEAF_schema.items():
        # Set up the file data dictionary
        file_data_dict = {}

        # Iterate through each file's tabs
        # Filter out metadata (non-tab-name) keys
        tab_types = [
            key for key in file_data.keys() if key != "default_filename"
        ]

        for tab_type in tab_types:
            # Get all tab data
            tab_data = file_data[tab_type]

            # Set up the final tab data dictionary
            tab_data_dict = {}

            # Save the tab metadata to the dictionary
            tab_data_dict["ALEAF_name"] = tab_data["tab_name"]
            tab_data_dict["field_order"] = tab_data["field_order"]
            tab_data_dict["field_orientation"] = tab_data["field_orientation"]

            # Filter out metadata (non-field-name) keys
            field_names = [
                key for key in tab_data.keys() 
                if key not in [
                    "tab_name",
                    "field_order",
                    "field_orientation"
                ]
            ]

            for field in field_names:
                field_data = tab_data[field]
                print(field_data["description"])

            # Save the tab data to the file data dictionary
            file_data_dict[tab_type] = tab_data_dict

        # Save the file data dictionary to the ALEAF data dictionary
        ALEAF_data_dict[file_name] = file_data_dict

    return ALEAF_data_dict


def finalize_dynamic_ALEAF_fields(ABCE_data, ALEAF_settings_data):

    return ALEAF_settings_data


def write_ALEAF_settings_file(ALEAF_settings_data, ALEAF_target_dir, file_name):
    if file_name in ["ALEAF_Master", "ALEAF_Master_LC_GEP"]:
        file_path = Path(ALEAF_target_dir) / "setting" / f"{file_name}.xlsx"
        data = ALEAF_settings_data[file_name]

        logging.info(f"ALEAF input file {file_name} written to {file_path}.")
    elif file_name == "ALEAF_portfolio":
        file_path = Path(ALEAF_target_dir) / "data" / ALEAF_region / f"ALEAF_{ALEAF_region}.xlsx"
        data = ALEAF_settings_data[file_name]

        logging.info("ALEAF input file {file_name} written to {file_path}.")
    else:
        logging.error(f"Unknown ALEAF input file type specified: {file_name}")
        logging.error(f"Please provide either ALEAF_Master, ALEAF_Master_LC_GEP, or ALEAF_portfolio")
        raise ValueError
    






if __name__ == "__main__":
    ALEAF_schema = read_ALEAF_schema()
    ABCE_data = read_ABCE_data()
    ALEAF_data_dict = initialize_ALEAF_data(ALEAF_schema, ABCE_data)
    print(ALEAF_data_dict)
    #ALEAF_settings_data = finalize_ALEAF_data(ABCE_data, ALEAF_schemas)
    #write_ALEAF_input_files(ALEAF_settings_data)


