# A script to create (or delete and recreate) the base ABCE database file
# Specification:
                                        # INIT    STEP
abce_tables = {"WIP_projects": 
                 [("asset_id", "integer"), # Julia, Python
                  ("agent_id", "text"),
                  ("period", "real"),
                  ("rcec", "real"),
                  ("rtec", "real"),
                  ("anpe", "real")
                 ],

               "assets":
                 [("asset_id", "integer"),
                  ("agent_id", "text"),
                  ("unit_type", "text"),
                  ("revealed", "text"),
                  #("start_pd", "real"),  #TODO: implement start period record
                  ("completion_pd", "real"),
                  ("cancellation_pd", "real"),
                  ("retirement_pd", "real"),
                  ("total_capex", "real"),
                  ("cap_pmt", "real")
                 ],

               "agent_params":
                 [("agent_id", "text"),
                  ("discount_rate", "real"),
                  ("tax_rate", "real"),
                  ("term_growth_rate", "real"),
                  ("debt_fraction", "real"),
                  ("debt_cost", "real"),
                  ("equity_cost", "real"),
                  ("interest_cap", "real")
                 ],

               "unit_specs":
                 [("unit_type", "text"),
                  ("fuel_type", "text"),
                  ("capacity", "real"),
                  ("uc_x", "real"),
                  ("d_x", "real"),
                  ("heat_rate", "real"),
                  ("VOM", "real"),
                  ("FOM", "real"),
                  ("unit_life", "real"),
                  ("CF", "real"),
                  ("fuel_cost", "real")
                 ],

               "demand":
                 [("period", "real"),
                  ("demand", "real")
                 ],

               "price_curve":
                 [("lamda", "real")]
              }


import sqlite3
import os
import sys
import pandas as pd

def clear_db_file(abce_db, replace):
    if os.path.exists(abce_db):
        if replace:
            os.remove(abce_db)
            print(f"Existing file at {abce_db} deleted.")
        else:
            print(f"Okay, please remove the file at {abce_db}, or specify --replace on the command line to avoid this message.")
            print("Terminating...")
            exit()


def create_db_file(abce_db):
    print(f"Creating a new database file at {abce_db}.")
    db = sqlite3.connect(abce_db)
    cur = db.cursor()
    return db, cur


def make_table(cur, table_name):
    sql_cols = []
    for column in abce_tables[table_name]:
        sql_cols.append(f"{column[0]} {column[1]}")
    cmd = f"CREATE TABLE {table_name} (" + ", ".join(sql_cols) + ")"
    cur.execute(cmd)


def create_all_tables(cur):
    for table in abce_tables:
        make_table(cur, table)


def create_database(db_file_name, replace=False):
    # Check whether the specified file already exists and delete it if allowed
    clear_db_file(db_file_name, replace)
    # Create a database seed file and associated cursor object
    db, cur = create_db_file(db_file_name)
    # Create all tables in the database
    create_all_tables(cur)
    # Commit changes and close the connection to the database
    db.commit()
    print(f"Database created in file '{db_file_name}'.")
    return db, cur


# Set name and path for ABCE database file
if __name__ == "__main__":
    create_database(sys.argv[1])
    # TODO: add argparse here
