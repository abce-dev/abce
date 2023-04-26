import pandas as pd
from pathlib import Path

def write_raw_db_to_excel(abce_model, settings):
    # Get the names of all database tables
    db_tables = pd.read_sql_query(
                    "SELECT name FROM sqlite_master WHERE type='table';",
                     abce_model.db
                )

    # Set up the path to the ultimate outputs directory
    out_file = Path(
                   settings["file_paths"]["ABCE_abs_path"] /
                   "outputs" /
                   settings["simulation"]["scenario_name"] /
                   settings["file_paths"]["output_file"]
               )

    # Write each table to a tab in the excel sheet
    with pd.ExcelWriter(out_file) as writer:
        for i in range(len(db_tables)):
            # Get table name
            table = db_tables.loc[i, "name"]

            # Get table data
            table_data = pd.read_sql_query(
                             f"SELECT * FROM {table}",
                             abce_model.db
                         )

            # Write table data to excel tab
            table_data.to_excel(writer, sheet_name=table, engine="openpyxl")
