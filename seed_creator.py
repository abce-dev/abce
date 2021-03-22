# A script to create (or delete and recreate) the base ABCE database file
# Specification:
# table 1: xtr_projects
#       HEADER         TYPE   SET BY
#    - asset_id      : text : Julia (init), Python (step)
#    - agent_id      : text : Julia (init), Python (step)
#    - period        : real : Julia (init), Python (step)
#    - rcec          : real : Julia (init), Python (step)
#    - rtec          : real : Julia (init), Python (step)
#    - anpe          : real : Julia (init), Julia (step)
#
# table 2: assets
#       HEADER            TYPE   SET BY
#    - asset_id         : text : Python
#    - agent_id         : text : Python
#    - completion_pd    : text : Python
#    - cancellation_pd  : text : Julia
#    - retirement_pd    : real : Python
#    - cap_pmt          : real : Python


import sqlite3
import os

def set_database_filename():
    abce_db = "./abce_db.db"
    return abce_db

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
        user_resp = input("Is it OK to delete it? [y/n]")
    return user_resp

def create_db_file(abce_db):
    print(f"Creating a new database file at {abce_db}.")
    db = sqlite3.connect(abce_db)
    cur = db.cursor()
    return db, cur

def create_assets_table(cur):
    cur.execute("""CREATE TABLE assets
                 (asset_id text, agent_id text, completion_pd real,
                 cancellation_pd real, retirement_pd real, cap_pmt real)""")

def create_WIP_projects_table(cur):
    cur.execute("""CREATE TABLE WIP_projects
                 (asset_id text, agent_id text, period real, rcec real,
                 rtec real, anpe real)""")

# Set name and path for ABCE database file
if __name__ == "__main__":
    # Get the desired filename for the database
    abce_db = set_database_filename()

    # Check whether the file already exists, and ask user's permission to delete
    #    if an existing file is found
    clear_db_file(abce_db)

    # Create a database seed file and associated cursor object
    db, cur = create_db_file(abce_db)

    # Create the `assets` and `WIP_projects` tables, with associated headers
    create_assets_table(cur)
    create_WIP_projects_table(cur)

    # Commit changes and close the connection to the database
    db.commit()
    db.close()

    print(f"Database created in file '{abce_db}'.")
