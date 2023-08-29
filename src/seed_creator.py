##########################################################################
# Copyright 2023 Argonne National Laboratory
#
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

# A script to create (or delete and recreate) the base ABCE database file

import sqlite3
import os
import sys
import logging
import pandas as pd

# Database Specification:
abce_tables = {
    "WIP_projects": [
        ("asset_id", "integer", "PRIMARY KEY"),
        ("agent_id", "integer"),
        ("period", "real"),
        ("cum_occ", "real"),
        ("rcec", "real"),
        ("cum_construction_duration", "real"),
        ("rtec", "real"),
        ("cum_construction_exp", "real"),
        ("anpe", "real"),
    ],
    "assets": [
        ("asset_id", "integer", "PRIMARY KEY"),
        ("agent_id", "text"),
        ("unit_type", "text"),
        ("start_pd", "integer"),
        ("completion_pd", "integer"),
        ("cancellation_pd", "integer"),
        ("retirement_pd", "integer"),
        ("total_capex", "real"),
        ("cap_pmt", "real"),
        ("C2N_reserved", "integer"),
    ],
    # Table to temporarily hold updates on construction progress
    #   during decision rounds
    # Agents append their own updates to this table at the end of
    #   their decision turn. After all decision turns in each round,
    #   the Model object executes these updates into the
    #   'WIP_projects' table, and empties this table.
    "WIP_updates": [
        ("asset_id", "integer", "PRIMARY KEY"),
        ("agent_id", "integer"),
        ("period", "real"),
        ("cum_occ", "real"),
        ("rcec", "real"),
        ("cum_construction_duration", "real"),
        ("rtec", "real"),
        ("cum_construction_exp", "real"),
        ("anpe", "real"),
    ],
    # Table to temporarily hold agent decisions about asset status
    #   updates (e.g. retirements) during decision rounds
    # Agents append their own updates to this table at the end of
    #   their decision turn. After all decision turns in this round,
    #   the Model object executes these updates into the 'assets'
    #   table, and empties this table.
    "asset_updates": [
        ("asset_id", "integer", "PRIMARY KEY"),
        ("agent_id", "integer"),
        ("unit_type", "text"),
        ("start_pd", "integer"),
        ("completion_pd", "integer"),
        ("cancellation_pd", "integer"),
        ("retirement_pd", "integer"),
        ("total_capex", "real"),
        ("cap_pmt", "real"),
        ("C2N_reserved", "integer"),
    ],
    "agent_params": [
        ("agent_id", "text", "PRIMARY KEY"),
        ("debt_fraction", "real"),
        ("cost_of_debt", "real"),
        ("cost_of_equity", "real"),
        ("starting_debt", "real"),
        ("starting_PPE", "real"),
    ],
    "unit_specs": [
        ("unit_type", "text", "PRIMARY KEY"),
        ("fuel_type", "text"),
        ("capacity", "real"),  # MW
        ("overnight_capital_cost", "real"),  # $/kW
        ("construction_duration", "real"),  # years
        ("heat_rate", "real"),  # MMBTU/MWh
        ("VOM", "real"),  # $/MWh
        ("FOM", "real"),  # $/kW-yr
        ("unit_life", "real"),  # years
        ("capacity_factor", "real"),  # frac
        ("max_PL", "real"),  # frac
        ("min_PL", "real"),  # frac
        ("ramp_up_limit", "real"),  # frac PL change/hr
        ("ramp_down_limit", "real"),  # frac PL change/hr
        ("max_regulation", "real"), # max. regulation participation, % PL
        ("max_spinning_reserve", "real"),
        ("max_nonspinning_reserve", "real"),
        ("no_load_cost", "real"),
        ("FC_per_MMBTU", "real"),  # $/MMBTU
        ("FC_per_MWh", "real"),  # $/MWh
        ("is_VRE", "text"),  # boolean
        ("emissions_per_MMBTU", "real"),  # tCO2 / MMBTU
        # net $/MWh subsidy or penalty due to carbon tax, PTC, etc.,
        ("policy_adj_per_MWh", "real"),
        # lead time between xtr start and cpp retirement
        ("tax_credits_per_MWh", "real"),
        ("tax_credits_per_MW", "real"),
        ("cpp_ret_lead", "real"),
        # number of coal units which must be retired
        ("num_cpp_rets", "integer"),
        # time before xtr finish when revenues begin
        ("rev_head_start", "real"),
    ],
    "demand": [("period", "real"), ("demand", "real")],
    "model_params": [("parameter", "text"), ("value", "real")],
    "WIP_C2N": [
        ("asset_id", "integer"),
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
        ("total_time_rem", "real"),
    ],
    "WIP_C2N_updates": [
        ("asset_id", "integer"),
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
        ("total_time_rem", "real"),
    ],
    "financial_instrument_manifest": [
        ("agent_id", "integer", "PRIMARY_KEY"),
        ("instrument_id", "integer"),
        ("instrument_type", "text"),
        ("asset_id", "integer"),
        ("pd_issued", "integer"),
        ("initial_principal", "real"),
        ("maturity_pd", "integer"),
        ("rate", "real"),
    ],
    "financing_schedule": [
        ("instrument_id", "integer", "PRIMARY KEY"),
        ("agent_id", "integer"),
        ("base_pd", "integer"),
        ("projected_pd", "integer"),
        ("total_payment", "real"),
        ("interest_payment", "real"),
        ("principal_payment", "real"),
    ],
    "capex_projections": [
        ("agent_id", "integer", "PRIMARY KEY"),
        ("asset_id", "integer"),
        ("base_pd", "integer"),
        ("projected_pd", "integer"),
        ("capex", "real"),
    ],
    "depreciation_projections": [
        ("agent_id", "integer", "PRIMARY KEY"),
        ("asset_id", "integer"),
        ("completion_pd", "integer"),
        ("base_pd", "integer"),
        ("projected_pd", "integer"),
        ("depreciation", "real"),
        ("beginning_book_value", "real"),
    ],
    "agent_financial_statements": [
        ("agent_id", "integer", "PRIMARY KEY"),
        ("base_pd", "integer"),
        ("projected_pd", "integer"),
        ("capex", "real"),
        ("remaining_debt_principal", "real"),
        ("depreciation", "real"),
        ("revenue", "real"),
        ("VOM", "real"),
        ("FOM", "real"),
        ("fuel_cost", "real"),
        ("policy_adj", "real"),
        ("EBITDA", "real"),
        ("EBIT", "real"),
        ("interest_payment", "real"),
        ("EBT", "real"),
        ("tax_owed", "real"),
        ("tax_credits", "real"),
        ("net_income", "real"),
        ("FCF", "real"),
        ("dividends", "real"),
        ("retained_earnings", "real"),
        ("ICR", "real"),
        ("FCF_debt_ratio", "real"),
        ("RE_debt_ratio", "real"),
        ("moodys_score", "real"),
    ],
    "agent_decisions": [
        ("agent_id", "integer", "PRIMARY KEY"),
        ("base_pd", "integer"),
        ("unit_type", "text"),
        ("project_type", "text"),
        ("lag", "integer"),
        ("ret_pd", "integer"),
        ("NPV", "real"),
        ("allowed", "text"),
        ("units_to_execute", "integer"),
    ],
    "ALEAF_dispatch_results": [
        ("period", "integer"),
        ("unit_type", "text"),
        ("generation", "real"),
        ("reg_total", "real"),
        ("spin_total", "real"),
        ("nspin_total", "real"),
        ("gen_rev", "real"),
        ("reg_rev", "real"),
        ("spin_rev", "real"),
        ("nspin_rev", "real"),
        ("revenue", "real"),
        ("VOM", "real"),
        ("fuel_cost", "real"),
        ("FOM", "real"),
        ("policy_adj", "real"),
        ("tax_credits", "real"),
    ],
}


def ask_user_permission_to_delete(abce_db):
    user_resp = ""
    acceptable_responses = ["Y", "y", "Yes", "yes", "N", "n", "No", "no"]
    agree_responses = ["Y", "y", "Yes", "yes"]
    reply = False
    while user_resp not in acceptable_responses:
        user_resp = input(
            f"There is already a database file at {abce_db}. "
            + "Can I delete it? [y/n] "
        )
    if user_resp in agree_responses:
        reply = True
    return reply


def clear_db_file(abce_db, force):
    if os.path.exists(abce_db):
        if force:
            os.remove(abce_db)
            logging.info(f"Existing file at {abce_db} deleted.")
        else:
            user_response = ask_user_permission_to_delete(abce_db)
            if user_response:
                os.remove(abce_db)
                logging.log(45, f"Existing file at {abce_db} deleted.")
                logging.log(
                    45,
                    (
                        "(Hint: you can specify --force or -f on the command "
                        + "line to automatically delete an existing DB file.)"
                    ),
                )
            else:
                logging.log(
                    45,
                    (
                        f"DB file at {abce_db} not deleted. Please move it or "
                        + "specify a different file name."
                    ),
                )
                logging.log(45, "Terminating...")
                exit()


def create_db_file(abce_db):
    logging.info(f"Creating a new database file at {abce_db}.")
    db = sqlite3.connect(str(abce_db), timeout=10)
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
    logging.info(f"Database created in file '{db_file_name}'.")
    return db, cur


# Set name and path for ABCE database file
if __name__ == "__main__":
    create_database(sys.argv[1])
