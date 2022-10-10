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
    sim_config = pd.DataFrame(writer.sheets["Simulation Configuration"].values)
    sim_config.columns = sim_config.iloc[0]
    sim_config = sim_config.drop(sim_config.index[0]).reset_index().drop("index", axis=1)

    # Zero out all policy columns
    for col in policy_cols:
        sim_config[col] = 0
    for key, val in policies_dict.items():
        if (key in valid_CTAX_names) and (val["enabled"]):
            sim_config["CTAX"] = val["qty"]
        elif (key in valid_PTC_names) and (val["enabled"]):
            for unit_type in val["eligible"]:
                if unit_type in valid_wind_names:
                    sim_config["PTC_W"] = val["qty"]
                elif unit_type in valid_solar_names:
                    sim_config["PTC_S"] = val["qty"]
                elif unit_type in valid_nuclear_names:
                    sim_config["PTC_Nuc"] = val["qty"]

    # Write the updated sheet to file
    sim_config.to_excel(writer, sheet_name = "Simulation Configuration", header=True, index=False)
    writer.save()



