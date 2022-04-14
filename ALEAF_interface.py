# A module allowing interactions between ABCE and A-LEAF

import os
import pandas as pd
import openpyxl
import sqlite3

# import local modules
import seed_creator as sc


def prepare_xlsx_data(ref_data_file, destination_file):
    book = openpyxl.load_workbook(ref_data_file)
    writer = pd.ExcelWriter(destination_file, engine="openpyxl")
    writer.book = book
    writer.sheets = dict((sheet.title, sheet) for sheet in book.worksheets)

    return book, writer


def update_ALEAF_system_portfolio(ALEAF_portfolio_ref, ALEAF_portfolio_remote, db, current_pd):
    # Get the updated list of currently-operating units by unit type
    unit_type_count = pd.read_sql_query(f"SELECT unit_type, COUNT(unit_type) FROM assets WHERE completion_pd <= {current_pd} AND retirement_pd > {current_pd} AND cancellation_pd > {current_pd} GROUP BY unit_type", db)
    unit_type_count = unit_type_count.set_index("unit_type")

    # Retrieve and organize the A-LEAF system portfolio
    book, writer = prepare_xlsx_data(ALEAF_portfolio_ref, ALEAF_portfolio_remote)
    df = organize_ALEAF_portfolio(writer)

    # Update unit type numbers
    for unit_type in list(unit_type_count.index):
        df.loc[(df["bus_i"] == 1) & (df["Unit Type"] == unit_type), "EXUNITS"] = unit_type_count.loc[unit_type, "COUNT(unit_type)"]
    df.to_excel(writer, sheet_name="gen", header=True, index=False)
    writer.save()


def organize_ALEAF_portfolio(writer):
    df = pd.DataFrame(writer.sheets["gen"].values)
    df.columns = df.iloc[0]
    df = df.drop(df.index[0]).reset_index(drop=True)
    df["EXUNITS"] = df["EXUNITS"].astype("int64")
    df["Unit Type"] = df["Unit Type"].astype("string")

    return df


def update_ALEAF_model_settings(ALEAF_model_settings_ref, ALEAF_model_settings_remote, db, settings, period=0):
    # Read in ALEAF simulation settings from ALEAF_Master_LC_GEP.xlsx-style file
    book, writer = prepare_xlsx_data(ALEAF_model_settings_ref, ALEAF_model_settings_remote)
    sim_config = pd.DataFrame(writer.sheets["Simulation Configuration"].values)
    sim_config.columns = sim_config.iloc[0]
    sim_config = sim_config.drop(sim_config.index[0]).reset_index().drop("index", axis=1)

    # Read in demand data from DB
    demand = pd.read_sql_query(f"SELECT demand FROM demand WHERE period = {period}", db)

    # If this is the first step of the run, update top-level settings
    if period == 0:
        # Update ALEAF scenario name
        sim_config.loc[0, "Scenario"] = settings["ALEAF_scenario_name"]


    # Update periodic data
    # Update peak demand (PD, MW)
    sim_config.loc[sim_config["Scenario"] == settings["ALEAF_scenario_name"], "PD"] = demand.iloc[0, 0]

    # Write the updated sheet to file
    sim_config.to_excel(writer, sheet_name = "Simulation Configuration", header=True, index=False)
    writer.save()


def update_ALEAF_policy_settings(ALEAF_model_settings_ref, ALEAF_model_settings_remote, policies_dict, unit_specs):
    valid_CTAX_names = ["CTAX", "ctax", "carbon_tax", "carbontax"]
    valid_PTC_names = ["PTC", "ptc", "production_tax_credit", "productiontaxcredit"]

    # Read in ALEAF simulation settings from ALEAF_Master_LC_GEP.xlsx-style file
    book, writer = prepare_xlsx_data(ALEAF_model_settings_ref, ALEAF_model_settings_remote)
    sim_config = pd.DataFrame(writer.sheets["Simulation Configuration"].values)
    sim_config.columns = sim_config.iloc[0]
    sim_config = sim_config.drop(sim_config.index[0]).reset_index().drop("index", axis=1)

    # Zero out all policy columns
    policy_cols = ["CTAX", "RPS", "PTC_W", "PTC_S", "ITC_S", "ITC_W", "CAPPMT"]
    for col in policy_cols:
        sim_config[col] = 0
    for key, val in policies_dict.items():
        if (key in valid_CTAX_names) and (val["enabled"]):
            sim_config["CTAX"] = val["qty"]
        elif (key in valid_PTC_names) and (val["enabled"]):
            sim_config["PTC_W"] = val["qty"]
            sim_config["PTC_S"] = val["qty"]

    print(sim_config)

    # Write the updated sheet to file
    sim_config.to_excel(writer, sheet_name = "Simulation Configuration", header=True, index=False)
    writer.save()



