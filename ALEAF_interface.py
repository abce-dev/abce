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


def load_excel_workbook(filename):
    book = openpyxl.load_workbook(filename)
    writer = pd.ExcelWriter(filename, engine="openpyxl")
    writer.book = book
    writer.sheets = dict((sheet.title, sheet) for sheet in book.worksheets)

    return book, writer


def set_ALEAF_pwd(ALEAF_master_settings_file, ALEAF_absolute_path):
    """
    Inplace operation to update the ALEAF pwd setting to the desired location

    TODO: Fix magic number [4] in update.to_excel() line (should choose row
      number based on a search for "pwd" in the first column).
    """
    book, writer = load_excel_workbook(ALEAF_master_settings_file)
    update = pd.DataFrame({"Setting": "pwd_location", "Value": f"{ALEAF_absolute_path}"}, index=[0])
    update.to_excel(writer, sheet_name="ALEAF Master Setup", startrow=4, startcol=0, header=False, index=False)
    writer.save()


def update_ALEAF_system_portfolio(ALEAF_sys_portfolio_path, db, current_pd):
    new_assets = get_new_units(db, current_pd)
    book, writer = load_excel_workbook(ALEAF_sys_portfolio_path)

    df = get_organized_ALEAF_portfolio(writer)
    for unit_type in list(new_assets.index):
        df.loc[(df["bus_i"] == 1) & (df["Unit Type"] == unit_type), "EXUNITS"] += new_assets.loc[unit_type, "num_units"]
    df.to_excel(writer, sheet_name="gen", header=True, index=False)
    writer.save()


def get_organized_ALEAF_portfolio(writer):
    df = pd.DataFrame(writer.sheets["gen"].values)
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])

    return df


def update_ALEAF_demand(ALEAF_model_settings_path, db, current_pd):
    demand = pd.read_sql_query(f"SELECT demand FROM demand WHERE period = {current_pd}", db)
    book, writer = load_excel_workbook(ALEAF_model_settings_path)
    sim_config = pd.DataFrame(writer.sheets["Simulation Configuration"].values)
    sim_config.columns = sim_config.iloc[0]
    sim_config = sim_config.drop(sim_config.index[0])

    sim_config.loc[sim_config["PLOTORDER"] == 1, "PD"] = demand.iloc[0, 0]
    sim_config.to_excel(writer, sheet_name = "Simulation Configuration", header=True, index=False)
    writer.save()



