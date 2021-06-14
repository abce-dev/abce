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
    print(new_assets)
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
    ALEAF_portfolio = pd.DataFrame(writer.sheets["gen"].values)

    df = pd.DataFrame(writer.sheets["gen"].values)
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])
    for unit_type in list(new_assets.index):
        df.loc[(df["bus_i"] == 1) & (df["Unit Type"] == unit_type), "EXUNITS"] += new_assets.loc[unit_type, "num_units"]
    df.to_excel(writer, sheet_name="gen", header=True, index=False)
    writer.save()







