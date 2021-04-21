# A script to create (or delete and recreate) the base ABCE database file
# Specification:
                                        # INIT    STEP
abce_tables = {"WIP_projects": 
                 [("asset_id", "text"), # Julia, Python
                  ("agent_id", "text"),
                  ("period", "real"),
                  ("rcec", "real"),
                  ("rtec", "real"),
                  ("anpe", "real")
                 ],

               "assets":
                 [("asset_id", "text"),
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

def clear_db_file(abce_db):
    if os.path.exists(abce_db):
        user_resp = get_user_consent_to_delete(abce_db)
        if user_resp in ["Y", "y", "Yes", "yes", "YES"]:
            os.remove(abce_db)
            print(f"Existing file at {abce_db} deleted.")
        else:
            print(f"Okay, please remove the file at {abce_db} before running this script again.")
            print("Terminating...")
            exit()


def get_user_consent_to_delete(abce_db):
    print(f"A file already exists at {abce_db}.")
    user_resp = ""
    valid_responses = ["Y", "y", "N", "n", "yes", "no", "Yes", "No"]
    while user_resp not in valid_responses:
        user_resp = input("Is it OK to delete it? [y/n] ")
    return user_resp


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


def create_database(db_file_name):
    # Check whether the specified file already exists; if so, ask the user's
    #    permission to delete it
    clear_db_file(db_file_name)
    # Create a database seed file and associated cursor object
    db, cur = create_db_file(db_file_name)
    # Create all tables in the database
    create_all_tables(cur)
    # Commit changes and close the connection to the database
    db.commit()
    db.close()
    print(f"Database created in file '{db_file_name}'.")



# Set name and path for ABCE database file
if __name__ == "__main__":
    create_database(sys.argv[1])
