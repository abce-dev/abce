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
#       HEADER         TYPE   SET BY
#    - asset_id      : text : Python
#    - agent_id      : text : Python
#    - is_complete  : text : Python
#    - is_cancelled  : text : Python
#    - ret_pd        : real : Python
#    - cap_pmt       : real : Python


import sqlite3
import os

# Set name and path for ABCE database file
abce_db = "./abce_db.db"

# If file already exists, ask if it's OK to delete it
if os.path.exists(abce_db):
    print(f"A file already exists at {abce_db}.")
    user_resp = ""
    valid_responses = ["Y", "y", "N", "n", "yes", "no", "Yes", "No"]
    yes = ["Y", "y", "yes", "Yes"]
    while user_resp not in valid_responses:
        user_resp = input("Is it OK to delete it? [y/n]")
    if user_resp in yes:
        os.remove(abce_db)
        print(f"Existing file at {abce_db} deleted.")
    else:
        print(f"Please move the file at {abce_db} before running this script again.")
        print("Terminating...")
        exit()
else:
    print(f"Creating new database file at {abce_db}...")


# Create a header-only database seed for projects
con = sqlite3.connect(abce_db)
cur = con.cursor()
cur.execute('''CREATE TABLE xtr_projects
               (asset_id text, agent_id text, period real, rcec real, rtec real, anpe real)''')
cur.execute('''CREATE TABLE assets
               (asset_id text, agent_id text, is_complete text, is_cancelled text, ret_pd real, cap_pmt real)''')

con.commit()
con.close()

print(f"File {abce_db} created.")
