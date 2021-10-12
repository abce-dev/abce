# A module allowing interactions between ABCE and A-LEAF

import os
import pandas as pd
import openpyxl
import sqlite3

# import local modules
import seed_creator as sc


def get_new_units(db, current_pd):
    # Retrieve the unit types of all units completed during the current period
    new_assets = pd.read_sql_query("SELECT unit_type FROM assets " +
                                   f"WHERE completion_pd = {current_pd}", db)
    # Count the number of units of each type using groupby()
    new_assets["num_units"] = 1
    new_assets = new_assets.groupby("unit_type").sum()
    return new_assets


def prepare_xlsx_data(ref_data_file, destination_file):
    book = openpyxl.load_workbook(ref_data_file)
    writer = pd.ExcelWriter(destination_file, engine="openpyxl")
    writer.book = book
    writer.sheets = dict((sheet.title, sheet) for sheet in book.worksheets)

    return book, writer


def update_ALEAF_system_portfolio(ALEAF_portfolio_ref, ALEAF_portfolio_remote, db, current_pd):
    new_assets = get_new_units(db, current_pd)
    book, writer = prepare_xlsx_data(ALEAF_portfolio_ref, ALEAF_portfolio_remote)

    df = organize_ALEAF_portfolio(writer)
    for unit_type in list(new_assets.index):
        df.loc[(df["bus_i"] == 1) & (df["Unit Type"] == unit_type), "EXUNITS"] += new_assets.loc[unit_type, "num_units"]
    df.to_excel(writer, sheet_name="gen", header=True, index=False)
    writer.save()


def organize_ALEAF_portfolio(writer):
    df = pd.DataFrame(writer.sheets["gen"].values)
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])

    return df


def update_ALEAF_demand(ALEAF_model_settings_ref, ALEAF_model_settings_remote, db, period=0):
    demand = pd.read_sql_query(f"SELECT demand FROM demand WHERE period = {period}", db)
    book, writer = prepare_xlsx_data(ALEAF_model_settings_ref, ALEAF_model_settings_remote)
    sim_config = pd.DataFrame(writer.sheets["Simulation Configuration"].values)
    sim_config.columns = sim_config.iloc[0]
    sim_config = sim_config.drop(sim_config.index[0])

    sim_config.loc[sim_config["Scenario"] == "ABCE_base", "PD"] = demand.iloc[0, 0]
    sim_config.to_excel(writer, sheet_name = "Simulation Configuration", header=True, index=False)
    writer.save()


