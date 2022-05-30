# A script to create (or delete and recreate) the base ABCE database file

##########################################################################
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################

import sqlite3
import os
import sys
import pandas as pd

# Database Specification:
abce_tables = {"WIP_projects": 
                 [("asset_id", "integer", "PRIMARY KEY"),
                  ("agent_id", "integer"),
                  ("period", "real"),
                  ("cum_occ", "real"),
                  ("rcec", "real"),
                  ("cum_d_x", "real"),
                  ("rtec", "real"),
                  ("cum_exp", "real"),
                  ("anpe", "real")
                 ],

               "assets":
                 [("asset_id", "integer", "PRIMARY KEY"),
                  ("agent_id", "text"),
                  ("unit_type", "text"),
                  ("start_pd", "integer"),
                  ("completion_pd", "integer"),
                  ("cancellation_pd", "integer"),
                  ("retirement_pd", "integer"),
                  ("total_capex", "real"),
                  ("cap_pmt", "real"),
                  ("C2N_reserved", "integer")
                 ],

               # Table to temporarily hold updates on construction progress
               #   during decision rounds
               # Agents append their own updates to this table at the end of
               #   their decision turn. After all decision turns in each round,
               #   the Model object executes these updates into the
               #   'WIP_projects' table, and empties this table.
               "WIP_updates": 
                 [("asset_id", "integer", "PRIMARY KEY"),
                  ("agent_id", "integer"),
                  ("period", "real"),
                  ("cum_occ", "real"),
                  ("rcec", "real"),
                  ("cum_d_x", "real"),
                  ("rtec", "real"),
                  ("cum_exp", "real"),
                  ("anpe", "real")
                 ],

               # Table to temporarily hold agent decisions about asset status
               #   updates (e.g. retirements) during decision rounds
               # Agents append their own updates to this table at the end of
               #   their decision turn. After all decision turns in this round,
               #   the Model object executes these updates into the 'assets'
               #   table, and empties this table.
               "asset_updates":
                 [("asset_id", "integer", "PRIMARY KEY"),
                  ("agent_id", "integer"),
                  ("unit_type", "text"),
                  ("start_pd", "integer"),
                  ("completion_pd", "integer"),
                  ("cancellation_pd", "integer"),
                  ("retirement_pd", "integer"),
                  ("total_capex", "real"),
                  ("cap_pmt", "real"),
                  ("C2N_reserved", "integer")
                 ],

               "agent_params":
                 [("agent_id", "text", "PRIMARY KEY"),
                  ("tax_rate", "real"),
                  ("term_growth_rate", "real"),
                  ("debt_fraction", "real"),
                  ("cost_of_debt", "real"),
                  ("cost_of_equity", "real"),
                  ("starting_fcf", "real"),
                  ("starting_debt", "real"),
                  ("starting_PPE", "real")
                 ],

               "financing_schedule":
                 [("agent_id", "text", "PRIMARY KEY"),
                  ("period", "integer"),
                  ("type", "text"),
                  ("asset", "integer"),
                  ("original_principal", "real"),
                  ("term", "integer"),
                  ("outstanding_principal", "real")
                 ],

               "agent_debt":
                 [("agent_id", "text", "PRIMARY KEY"),
                  ("period", "integer"),
                  ("outstanding_principal", "real")
                 ],

               "unit_specs":
                 [("unit_type", "text", "PRIMARY KEY"),
                  ("fuel_type", "text"),
                  ("capacity", "real"),        # MW
                  ("uc_x", "real"),            # $/kW
                  ("d_x", "real"),             # years
                  ("heat_rate", "real"),       # BTU/Wh
                  ("VOM", "real"),             # $/MWh
                  ("FOM", "real"),             # $/kW-yr
                  ("unit_life", "real"),       # years
                  ("CF", "real"),              # frac
                  ("PMAX", "real"),            # frac
                  ("PMIN", "real"),            # frac
                  ("RUL", "real"),             # frac PL change/hr
                  ("RDL", "real"),             # frac PL change/hr
                  ("ATB_FC", "real"),          # $/MWh or $/MMBTU
                  ("ATB_FC_units", "text"),    # unit for ATB FC
                  ("FC_per_MWh", "real"),      # $/MWh
                  ("is_VRE", "text"),          # boolean
                  ("emissions_rate", "real"),  # tCO2 / MWh
                  ("policy_adj_per_MWh", "real"), # net $/MWh subsidy or penalty due to carbon tax, PTC, etc.,
                  ("cpp_ret_lead", "real"),    # lead time between xtr start and cpp retirement
                  ("rev_head_start", "real")   # time before xtr finish when revenues begin
                 ],


               "demand":
                 [("period", "real"),
                  ("demand", "real")
                 ],

               "price_curve":
                 [("lamda", "real")
                 ],

               "model_params":
                 [("parameter", "text"),
                  ("value", "real")
                 ],

               "WIP_C2N":
                 [("asset_id", "integer"),
                   ("C2N_type", "text"),
                   ("pd", "integer"),
                   ("license_issued", "boolean"),
                   ("cpp_dnd_cost_rem", "real"),
                   ("cpp_dnd_time_rem", "real"),
                   ("cpp_wr_cost_rem", "real"),
                   ("cpp_wr_time_rem", "real"),
                   ("cpp_nrc_cost_rem", "real"),
                   ("cpp_nrc_time_rem", "real"),
                   ("npp_ns_xtr_cost_rem", "real"),
                   ("npp_ns_xtr_time_rem", "real"),
                   ("npp_safety_xtr_cost_rem", "real"),
                   ("npp_safety_xtr_time_rem", "real"),
                   ("total_cost_rem", "real"),
                   ("total_time_rem", "real")
                 ],

               "WIP_C2N_updates":
                 [("asset_id", "integer"),
                   ("C2N_type", "text"),
                   ("pd", "integer"),
                   ("license_issued", "boolean"),
                   ("cpp_dnd_cost_rem", "real"),
                   ("cpp_dnd_time_rem", "real"),
                   ("cpp_wr_cost_rem", "real"),
                   ("cpp_wr_time_rem", "real"),
                   ("cpp_nrc_cost_rem", "real"),
                   ("cpp_nrc_time_rem", "real"),
                   ("npp_ns_xtr_cost_rem", "real"),
                   ("npp_ns_xtr_time_rem", "real"),
                   ("npp_safety_xtr_cost_rem", "real"),
                   ("npp_safety_xtr_time_rem", "real"),
                   ("total_cost_rem", "real"),
                   ("total_time_rem", "real")
                 ],


               "financial_instrument_manifest":
                 [("agent_id", "integer", "PRIMARY_KEY"),
                  ("instrument_id", "integer"),
                  ("instrument_type", "text"),
                  ("asset_id", "integer"),
                  ("pd_issued", "integer"),
                  ("initial_principal", "real"),
                  ("maturity_pd", "integer"),
                  ("rate", "real")
                 ],

               "agent_financing_schedule":
                 [("instrument_id", "integer", "PRIMARY KEY"),
                  ("agent_id", "integer"),
                  ("base_pd", "integer"),
                  ("projected_pd", "integer"),
                  ("total_payment", "real"),
                  ("interest_payment", "real"),
                  ("principal_payment", "real")
                 ],

               "capex_projections":
                 [("agent_id", "integer", "PRIMARY KEY"),
                  ("asset_id", "integer"),
                  ("base_pd", "integer"),
                  ("projected_pd", "integer"),
                  ("capex", "real")
                 ],

               "depreciation_projections":
                 [("agent_id", "integer", "PRIMARY KEY"),
                  ("asset_id", "integer"),
                  ("completion_pd", "integer"),
                  ("base_pd", "integer"),
                  ("projected_pd", "integer"),
                  ("depreciation", "real"),
                  ("beginning_book_value", "real")
                 ],

               "agent_financial_statements":
                 [("agent_id", "integer", "PRIMARY KEY"),
                  ("base_pd", "integer"),
                  ("projected_pd", "integer"),
                  ("revenue", "real"),
                  ("VOM", "real"),
                  ("FOM", "real"),
                  ("fuel_costs", "real"),
                  ("EBITDA", "real"),
                  ("depreciation", "real"),
                  ("EBIT", "real"),
                  ("interest_payment", "real"),
                  ("EBT", "real"),
                  ("tax_paid", "real"),
                  ("Net_Income", "real"),
                  ("capex", "real"),
                  ("FCF", "real")
                 ]

              }


def ask_user_permission_to_delete(abce_db):
    user_resp = ""
    acceptable_responses = ["Y", "y", "Yes", "yes", "N", "n", "No", "no"]
    agree_responses = ["Y", "y", "Yes", "yes"]
    reply = False
    while user_resp not in acceptable_responses:
        user_resp = input("There is already a database file at {abce_db}. Can I delete it? [y/n] ")
    if user_resp in agree_responses:
        reply = True
    return reply
    

def clear_db_file(abce_db, force):
    if os.path.exists(abce_db):
        if force:
            os.remove(abce_db)
            print(f"Existing file at {abce_db} deleted.")
        else:
            user_response = ask_user_permission_to_delete(abce_db)
            if user_response:
                os.remove(abce_db)
                print(f"Existing file at {abce_db} deleted.")
                print("(Hint: you can specify --force or -f on the command line to automatically delete an existing DB file.)")
            else:
                print("DB file at {abce_db} not deleted. Please move it or specify a different file name.")
                print("Terminating...")
                exit()


def create_db_file(abce_db):
    print(f"Creating a new database file at {abce_db}.")
    db = sqlite3.connect(abce_db, timeout=10)
    cur = db.cursor()
    return db, cur


def make_table(cur, table_name):
    sql_cols = []
    for column in abce_tables[table_name]:
        sql_cols.append(f"{column[0]} {column[1]}")
    cmd = f"CREATE TABLE {table_name} (" + ", ".join(sql_cols) + ")"
    cur.execute(cmd)


def create_database(db_file_name, replace=False):
    # Check whether the specified file already exists and delete it if allowed
    clear_db_file(db_file_name, replace)
    # Create a database seed file and associated cursor object
    db, cur = create_db_file(db_file_name)
    # Create all tables in the database
    for table in abce_tables:
        make_table(cur, table)
    # Commit changes and close the connection to the database
    db.commit()
    print(f"Database created in file '{db_file_name}'.")
    return db, cur


# Set name and path for ABCE database file
if __name__ == "__main__":
    create_database(sys.argv[1])
    # TODO: add argparse here
