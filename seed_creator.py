# A script to create (or delete and recreate) the base ABCE database file
# Specification:
# table 1: WIP_projects
#       HEADER         TYPE   INIT BY  STEPPED BY
#    - asset_id      : text : Julia,   Python
#    - agent_id      : text : Julia,   Python
#    - period        : real : Julia,   Python
#    - rcec          : real : Julia,   Python
#    - rtec          : real : Julia,   Python
#    - anpe          : real : Julia,   Julia 
#
# table 2: assets
#       HEADER            TYPE   INIT BY        STEPPED BY
#    - asset_id         : text : Python/Julia,  n/a
#    - agent_id         : text : Python/Julia,  n/a
#    - completion_pd    : text : Python/Julia,  Python
#    - cancellation_pd  : text : Julia,         Julia
#    - retirement_pd    : real : Python/Julia,  Python
#    - cap_pmt          : real : Python/Julia,  Python
#    - unit_type        : text : Python/Julia,  n/a
#
# table 3: agent_params
#       HEADER             TYPE   SET BY
#    - agent_id          : text : model
#    - discount_rate     : real : model
#    - tax_rate          : real : model
#    - term_growth_rate  : real : model
#    - debt_fraction     : real : model
#    - interest_cap      : real : model
#
# table 4: unit_specs
#
#    - unit_type    : text : seed_creator
#    - fuel_type    : text : seed_creator
#    - capacity     : real : seed_creator
#    - uc_x         : real : seed_creator
#    - d_x          : real : seed_creator
#    - heat_rate    : real : seed_creator
#    - VOM          : real : seed_creator
#    - FOM          : real : seed_creator
#    - unit_life    : real : seed_creator
#    - CF           : real : seed_creator
#    - fuel_cost    : real : seed_creator
# table 5: demand
#
#    - period       : real : model
#    - demand       : real : model




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
        user_resp = input("Is it OK to delete it? [y/n] ")
    return user_resp


def create_db_file(abce_db):
    print(f"Creating a new database file at {abce_db}.")
    db = sqlite3.connect(abce_db)
    cur = db.cursor()
    return db, cur


def create_assets_table(cur):
    cur.execute("""CREATE TABLE assets
                 (asset_id text, agent_id text, unit_type text,
                  completion_pd real, cancellation_pd real, retirement_pd real,
                  cap_pmt real)""")


def create_WIP_projects_table(cur):
    cur.execute("""CREATE TABLE WIP_projects
                 (asset_id text, agent_id text, period real, rcec real,
                  rtec real, anpe real)""")


def create_agent_params_table(cur):
    cur.execute("""CREATE TABLE agent_params
                   (agent_id text, discount_rate real, tax_rate real,
                    term_growth_rate real, debt_fraction real,
                    interest_cap real)""")

def create_unit_specs_table(cur):
    cur.execute("""CREATE TABLE	unit_specs
                   (unit_type text, fuel_type text, capacity real, uc_x real,
                    d_x real, heat_rate real, VOM real, FOM real,
                    unit_life real, CF real, fuel_cost real)""")

def create_demand_table(cur):
    cur.execute("CREATE TABLE demand (period real, demand real)")




# Set name and path for ABCE database file
if __name__ == "__main__":
    # Get the desired filename for the database
    abce_db = set_database_filename()

    # Check whether the file already exists, and ask user's permission to delete
    #    if an existing file is found
    clear_db_file(abce_db)

    # Create a database seed file and associated cursor object
    db, cur = create_db_file(abce_db)

    # Create the `assets`, `WIP_projects`, and `agent_params` tables, 
    #    including headers
    create_assets_table(cur)
    create_WIP_projects_table(cur)
    create_agent_params_table(cur)
    create_unit_specs_table(cur)
    create_demand_table(cur)

    # Commit changes and close the connection to the database
    db.commit()
    db.close()

    print(f"Database created in file '{abce_db}'.")
