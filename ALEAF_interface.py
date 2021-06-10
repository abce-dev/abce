# A module allowing interactions between ABCE and A-LEAF

import os
import pandas as pd
import openpyxl
import sqlite3

# import local modules
import seed_creator as sc


def get_new_units(db, current_pd):
    # Retrieve the unit types of all units completed during the current period
    new_assets = pd.read_sql_query("SELECT unit_type FROM assets" +
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
    update.to_excel(writer, "ALEAF Master Setup", startrow=4, startcol=0, header=False, index=False)
    writer.save()


def load_ALEAF_system_portfolio(ALEAF_sys_portfolio_path):
    ALEAF_system_portfolio = pd.read_excel(ALEAF_sys_portfolio_path,
                                           engine="openpyxl",
                                           sheet_name="gen")
    return ALEAF_system_portfolio


def load_non_portfolio_sheets(ALEAF_sys_portfolio_path):
    ALEAF_xlsx_file = pd.ExcelFile(ALEAF_sys_portfolio_path, engine="openpyxl")
    sheet_names = ALEAF_xlsx_file.sheet_names
    sheet_names.remove("gen")   # The "gen" sheet holds the system portfolio data
    sheet_data = dict()
    for sheet in sheet_names:
        sheet_data[sheet] = pd.read_excel(ALEAF_sys_portfolio_path,
                                          engine="openpyxl",
                                          sheet_name=sheet)
    return sheet_data


def split_ALEAF_system_portfolio_by_bus(ALEAF_system_portfolio, used_bus=1):
    other_buses_portfolio = ALEAF_system_portfolio[ALEAF_system_portfolio["bus_i"] != used_bus]
    ALEAF_system_portfolio = ALEAF_system_portfolio[ALEAF_system_portfolio["bus_i"] == used_bus]
    return ALEAF_system_portfolio, other_buses_portfolio


def recombine_ALEAF_system_portfolio_by_bus(ALEAF_system_portfolio, other_buses_portfolio):
    ALEAF_system_portfolio = pd.concat([ALEAF_system_portfolio, other_buses_portfolio],
                                       axis=0,
                                       ignore_index=True)
    return ALEAF_system_portfolio


def update_ALEAF_system_portfolio(new_assets, ALEAF_system_portfolio):
    # Set the index to "Unit Type" for easier referencing
    ALEAF_system_portfolio = ALEAF_system_portfolio.set_index("Unit Type")
    # Add the appropriate number of units to the system portfolio
    for unit in list(new_assets.index):
        ALEAF_system_portfolio.loc[unit, "EXUNITS"] += new_assets.loc[unit, "num_units"]
    # Reset "Unit Type" to be a column again
    ALEAF_system_portfolio = ALEAF_system_portfolio.reset_index()


def overwrite_ALEAF_system_file(ALEAF_sys_portfolio_path, system_portfolio, non_portfolio_sheets):
    """
    Recombine the system portfolio dataframe and the dictionary of
    non-portfolio dataframes into a single excel sheet under the original
    sheet name

    Arguments:
      system_portfolio (DataFrame): the dataframe containing the ALEAF system
        portfolio
      non_portfolio_sheets (dict): a dictionary containing all of the
        non-portfolio Excel file sheets as individual DataFrames
    """

    with pd.ExcelWriter(ALEAF_sys_portfolio_path) as writer:
        system_portfolio.to_excel(writer, sheet_name="gen", index=False)
        for sheet in list(non_portfolio_sheets.keys()):
            non_portfolio_sheets[sheet].to_excel(writer, sheet_name=sheet, index=False) 












