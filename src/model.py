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

import math
import os
import shutil
import subprocess
import yaml
import numpy as np
import pandas as pd
import sqlite3
import logging
from pathlib import Path
from mesa import Agent, Model
from mesa.time import RandomActivation
# import nrelpy.atb as ATB

# import local modules
from .agent import GenCo
from . import ABCEfunctions as ABCE
from . import seed_creator as sc
from . import dispatch_ppx as dsp
from . import input_data_management as idm

import warnings
warnings.filterwarnings("ignore")


class GridModel(Model):
    ''' A model with some number of GenCos. '''

    def __init__(self, settings, args):
        # Copy the command-line arguments as member data
        self.args = args
        self.settings = settings

        # If verbosity is 2 or 3, show the ABCE splash header
        if self.args.verbosity >= 2:
            self.show_abce_header()

        # Check ./outputs/ dir and clear out old files
        self.prepare_outputs_directory()

        # Initialize database for storing and managing all simulation data
        self.db_file = (Path.cwd() / 
                        "outputs" / 
                        settings["simulation"]["ALEAF_scenario_name"] / 
                        settings["file_paths"]["db_file"]
                       )
        self.db, self.cur = sc.create_database(self.db_file, self.args.force)

        # Initialize all data
        self.load_all_data()

        # Create the local tmp/ directory inside the current working directory,
        #   if it doesn't already exist
        tmp_dir_location = (Path.cwd() / "tmp")
        Path(tmp_dir_location).mkdir(exist_ok=True)

        # If running A-LEAF, set up any necessary file paths
        if self.settings["simulation"]["run_ALEAF"]:
            self.set_ALEAF_file_paths()

        # Initialize the model one time step before the true start date
        self.current_pd = self.settings["constants"]["time_before_start"]

        # Define a dictionary to hold all agents for easier indexing
        self.agents = {}

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Create agents
        for agent_id, agent_params in self.agent_specs.items():
            agent = GenCo(
                agent_id,
                self,
                agent_params
            )

            # Initialize assets owned by this agent, and add them to
            #   the database
            self.initialize_agent_assets(agent)

            # Add this agent to the dictionary of agents
            self.agents[agent_id] = agent

            # Add this agent to the schedule
            self.schedule.add(agent)

        # Make sure all pending database transactions are resolved
        self.db.commit()


    def load_all_data(self):
        # Retrieve the input data
        self.unit_specs = idm.initialize_unit_specs(self.settings)

        # Save unit_specs to the database
        self.add_unit_specs_to_db()

        # Load all-period demand data into the database
        self.load_demand_data_to_db()

        # Add model parameters to the database
        self.load_model_parameters_to_db()

        # Read agent specification data
        self.load_agent_specifications()


    def load_agent_specifications(self):
        # Read in the GenCo parameters data from file
        agent_specs_file_name = Path(
            self.settings["file_paths"]["ABCE_abs_path"] /
            self.settings["file_paths"]["agent_specifications_file"]
        )
        self.agent_specs = yaml.load(
            open(
                agent_specs_file_name,
                'r'),
            Loader=yaml.FullLoader)


    def prepare_outputs_directory(self):
        # Set up the output data path
        self.ABCE_output_data_path = (
            Path(os.getcwd()) /
            "outputs" /
            self.settings["simulation"]["ALEAF_scenario_name"]
        )

        if not Path(self.ABCE_output_data_path).is_dir():
            # If the desired output directory doesn't already exist, create it
            Path(self.ABCE_output_data_path).mkdir(exist_ok=True, parents=True)
        else:
            # Otherwise, delete any existing files in the directory
            for existing_file in Path(self.ABCE_output_data_path).iterdir():
                logging.debug(f"Deleting file {existing_file} from the output directory")
                (Path(self.ABCE_output_data_path) / existing_file).unlink()


    def show_abce_header(self):
        logo_file = Path(
                        self.settings["file_paths"]["ABCE_abs_path"] /
                        "misc" /
                        self.settings["file_paths"]["logo"]
                    )
        with open(logo_file, "r") as logo:
            for line in logo.read().splitlines():
                logging.log(self.settings["constants"]["vis_lvl"], line)


    def set_ALEAF_file_paths(self):
        """ Set up all absolute paths to ALEAF and its input files, and
              save them as member data.
        """
        # Set path to ALEAF outputs
        self.ALEAF_output_data_path = Path(
            Path(os.environ["ALEAF_DIR"]) /
            "output" /
            self.settings["ALEAF"]["ALEAF_model_type"] /
            self.settings["ALEAF"]["ALEAF_region"] /
            f"scenario_1_{self.settings['simulation']['ALEAF_scenario_name']}"
        )


    def add_unit_specs_to_db(self):
        """
        This function selects only the unit_specs columns needed for the
          database, then saves the table to the database.
        """
        # Reshape unit_specs into a dataframe
        unit_specs = pd.DataFrame.from_dict(self.unit_specs, orient="index").reset_index().rename(columns={"index": "unit_type"})

        # Retrieve the column-header schema for the 'unit_specs' table
        self.cur.execute("SELECT * FROM unit_specs")
        unit_specs_db_cols = [element[0] for element in self.cur.description]

        # Select the necessary DB columns in order
        unit_specs = unit_specs[unit_specs_db_cols]

        # Save the finalized unit specs data to the DB, and set the member data
        unit_specs.to_sql(
            "unit_specs",
            self.db,
            if_exists="replace",
            index=False)



    def initialize_agent_retirement_schedule(self, agent):
        # Validate the number of retirements versus number of owned units
        for unit_type, num_units in agent.starting_portfolio.items():
            total_ret_units = 0
            if unit_type not in agent.scheduled_retirements.keys():
                # If this unit type is not found in the agent's retirement
                #   schedule, initialize an empty dictionary
                agent.scheduled_retirements[unit_type] = {}

            # Determine the total number of units of this type scheduled
            #   scheduled to be retired by this agent
            for retirement_pd, ret_num_units in agent.scheduled_retirements[unit_type].items():
                total_ret_units += ret_num_units

            # Reconcile total number of scheduled retirements with total number
            #   of owned assets of this type
            if total_ret_units > num_units:
                raise ValueError(f"Portfolio specification mismatch for agent #{agent_id}: total owned {unit_type} units = {num_units}, but total specified {unit_type} retirements = {total_ret_units}.")

            # If there are any units with unspecified retirement dates, set
            #   their retirement date to a large number
            elif total_ret_units < num_units:
                agent.scheduled_retirements[unit_type][self.settings["constants"]["distant_time"]] = num_units - total_ret_units



    def initialize_agent_assets(self, agent):
        # Set up the agent's retirement schedule, and validate that the 
        #   schedule of retirements is compatible with the agent's portfolio
        self.initialize_agent_retirement_schedule(agent)

        # Retrieve the column-header schema for the 'assets' table
        self.cur.execute("SELECT * FROM assets")
        assets_col_names = [element[0] for element in self.cur.description]

        # Create a master dataframe to hold all asset records for this
        #   agent-unit type combination (to reduce the frequency of saving
        #   to the database)
        master_assets_df = pd.DataFrame(columns=assets_col_names)

        # Assign units to this agent as specified in the portfolio file,
        #   and record each in the master_asset_df dataframe
        for unit_type, unit_ret_data in agent.scheduled_retirements.items():
            # Create assets in blocks, according to the number of units per
            #   retirement period
            for retirement_pd, num_units in unit_ret_data.items():
                # Compute unit capex according to the unit type spec
                unit_capex = self.compute_total_capex_preexisting(unit_type)
                cap_pmt = 0

                # Set up all data values except for the asset id
                asset_dict = {"agent_id": agent.unique_id,
                              "unit_type": unit_type,
                              "start_pd": self.settings["constants"]["time_before_start"],
                              "completion_pd": 0,
                              "cancellation_pd": self.settings["constants"]["distant_time"],
                              "retirement_pd": int(retirement_pd),
                              "total_capex": unit_capex,
                              "cap_pmt": cap_pmt,
                              "C2N_reserved": 0
                             }

                # For each asset in this block, create a dataframe record and
                #   store it to the master_assets_df
                for i in range(num_units):
                    # Find the largest extant asset id, and set the current
                    #   asset id 1 higher
                    asset_dict["asset_id"] = max(
                        ABCE.get_next_asset_id(
                            self.db,
                            self.settings["constants"]["first_asset_id"]
                        ),
                        max(
                            master_assets_df["asset_id"],
                            default=self.settings["constants"]["first_asset_id"]
                        ) + 1
                    )

                    # Convert the dictionary to a dataframe format and save
                    new_record = pd.DataFrame(asset_dict, index=[0])
                    master_assets_df = master_assets_df.append(new_record)

        # Once all assets from all unit types for this agent have had records
        #   initialized, save the dataframe of all assets into the 'assets'
        #   DB table
        master_assets_df.to_sql(
            "assets",
            self.db,
            if_exists="append",
            index=False
        )

        self.db.commit()


    def compute_total_capex_preexisting(self, unit_type):
        unit_cost_per_kW = self.unit_specs[unit_type]["overnight_capital_cost"]
        unit_capacity = self.unit_specs[unit_type]["capacity"]

        total_capex = unit_cost_per_kW * unit_capacity * self.settings["constants"]["MW2kW"]

        return total_capex


    def create_unit_type_retirement_df(
            self, unit_type, agent_id, starting_asset_id):
        """
        Create the step-function mapping to determine which units are assigned
          which mandatory retirement periods.
        """
        # Filter the df of retirement period data for the current unit type
        unit_type_rets = self.ret_data[(self.ret_data["unit_type"] == unit_type) & (
            self.ret_data["agent_id"] == agent_id)].copy()

        # Sort from soonest to furthest-away retirement period
        unit_type_rets = unit_type_rets.sort_values(
            by="retirement_pd", axis=0, ascending=True, ignore_index=True)

        # Generate the thresholds for each retirement period, starting with
        #   the next available asset id. These thresholds indicate the END
        #   (i.e. NON-INCLUSIVE) of each assignment interval. A row with an
        #   `rp_threshold` of 2005 and a `retirement_pd` of 12 means that any
        #   assets with IDs strictly less than 2005 will receive a
        #   `retirement_pd` of 12 (unless they qualify for a lower threshold
        #   category).
        #unit_type_rets["rp_threshold"] = np.cumsum(unit_type_rets["num_copies"].to_numpy(), axis=0) + starting_asset_id

        return unit_type_rets

    def load_demand_data_to_db(self):
        # Load all-period demand data into the database
        demand_data_file = (
            Path(self.settings["file_paths"]["ABCE_abs_path"]) /
            self.settings["file_paths"]["demand_data_file"]
        )
        demand_df = pd.read_csv(demand_data_file) * self.settings["scenario"]["peak_demand"]

        # Create an expanded range of periods to backfill with demand_df data
        new_index = list(range(self.settings["demand"]["total_forecast_horizon"]))
        demand_df = demand_df.reindex(new_index, method="ffill")

        # Save data to DB
        demand_df.to_sql(
            "demand",
            self.db,
            if_exists="replace",
            index_label="period")
        self.db.commit()


    def load_model_parameters_to_db(self):
        # Load specific parameters specified in settings.yml to the database
        prm = self.settings["system"]["planning_reserve_margin"]
        self.cur.execute(f"INSERT INTO model_params VALUES ('PRM', {prm})")

        tax_rate = self.settings["system"]["tax_rate"]
        self.cur.execute(
            f"INSERT INTO model_params VALUES ('tax_rate', {tax_rate})")
        self.db.commit()



    def step(self, demo=False):
        """
        Advance the model by one step.
        """
        self.current_pd += 1

        if self.current_pd == 0:
            self.has_ABCE_sysimage, self.has_dispatch_sysimage = self.check_for_sysimage_files()

        self.display_step_header()

        # Advance the status of all WIP projects to the current period
        self.update_WIP_projects()

        # Update financial statements and financial projections for all agents
        self.update_agent_financials()

        # Compute the scenario reduction results for this year
        ABCE.execute_scenario_reduction(
            self.db,
            self.current_pd,
            self.settings,
            self.unit_specs
        )

        # Close the database to avoid access problems in the Julia scope
        self.db.commit()
        self.db.close()

        # Iterate through all agent turns
        self.schedule.step()

        logging.log(
            self.settings["constants"]["vis_lvl"],
            "\nAll agent turns are complete.\n"
        )

        # Reopen the database connection, now that Julia execution is complete
        self.db = sqlite3.connect(self.db_file, timeout=10)
        self.cur = self.db.cursor()

        # Show update tables in the terminal
        self.show_round_updates()

        # Transfer all decisions and updates from the 'asset_updates' and
        #   'WIP_updates' tables into their respective public-information
        #   equivalents
        self.execute_all_status_updates()

        if demo:
            logging.log(self.settings["constants"]["vis_lvl"], "\n")
            user_response = input("Press Enter to continue: ")

        if self.settings["simulation"]["run_ALEAF"]:
            # Re-load the baseline A-LEAF data
            ALEAF_data = idm.load_data(Path(self.settings["ALEAF"]["ALEAF_data_file"]))

            # Generate all three A-LEAF input files and save them to the 
            #   appropriate subdirectories in the A-LEAF top-level directory
            idm.create_ALEAF_files(self.settings, ALEAF_data, self.unit_specs, self.db, self.current_pd)

            # Run A-LEAF
            logging.log(self.settings["constants"]["vis_lvl"], "Running A-LEAF...")
            run_script_path = Path(os.environ["ALEAF_DIR"]) / "execute_ALEAF.jl"
            ALEAF_env_path = Path(os.environ["ALEAF_DIR"]) / "."
            ALEAF_sysimage_path = Path(os.environ["ALEAF_DIR"]) / "aleafSysimage.so"
            aleaf_cmd = f"julia --project={ALEAF_env_path} -J {ALEAF_sysimage_path} {run_script_path} {self.settings['ALEAF']['ALEAF_abs_path']}"

            if self.args.verbosity < 2:
                sp = subprocess.check_call(aleaf_cmd,
                                           shell=True,
                                           stdout=open(os.devnull, "wb"))
            else:
                sp = subprocess.check_call(aleaf_cmd, shell=True)

            self.save_ALEAF_outputs()
            self.process_ALEAF_dispatch_results()


    def display_step_header(self):
        if self.current_pd != 0:
            logging.log(self.settings["constants"]["vis_lvl"], "\n\n\n")
        logging.log(self.settings["constants"]["vis_lvl"], "=" * 60)
        logging.log(self.settings["constants"]["vis_lvl"], f"   Simulation step: {self.current_pd}")
        logging.log(self.settings["constants"]["vis_lvl"], "=" * 60)


    def show_round_updates(self):
        logging.debug("Updates to assets:")
        logging.debug(pd.read_sql("SELECT * FROM asset_updates", self.db))
        logging.debug("\nConstruction project updates (last 10 entries):")
        logging.debug(pd.read_sql("SELECT * FROM WIP_updates", self.db).tail(n=10))


    def check_for_sysimage_files(self):
        ABCE_sysimage_path = (Path(
            self.settings["file_paths"]["ABCE_abs_path"]) /
            "env" /
            self.settings["file_paths"]["ABCE_sysimage_file"]
        )

        dispatch_sysimage_path = (Path(
            self.settings["file_paths"]["ABCE_abs_path"]) /
            "env" /
            self.settings["file_paths"]["dispatch_sysimage_file"]
        )

        has_ABCE_sysimage = True
        has_dispatch_sysimage = True

        if not Path(ABCE_sysimage_path).exists():
            msg = (
                f"No sysimage file found at {ABCE_sysimage_path}. " +
                "Execution will proceed, but Julia may run extremely slowly. " +
                "If you already have a sysimage file, please move it to " +
                "the filename {sysimage_path}. If you do not have a " +
                "sysimage file, please run 'julia make_sysimage.jl --mode=abce' in this " +
                "directory.")
            logging.warn(msg)
            has_ABCE_sysimage = False

        if not Path(dispatch_sysimage_path).exists():
            msg = (
                f"No sysimage file found at {dispatch_sysimage_path}. " +
                "Execution will proceed, but the dispatch sub-module may run extremely slowly. " +
                "If you already have a dispatch sysimage file, please move it to " +
                "the filename {dispatch_sysimage_path}. If you do not have a " +
                "dispatch sysimage file, please run 'julia make_sysimage.jl --mode=dispatch' in this " +
                "directory.")
            logging.warn(msg)
            has_dispatch_sysimage = False

        return has_ABCE_sysimage, has_dispatch_sysimage

    def update_agent_financials(self):
        # Update the following database tables for all agents:
        #   - capex projections
        #   - financing instrument manifest
        #   - financing schedule
        #   - depreciation projections
        #   - agent financial statements
        self.update_capex_projections()
        self.update_financial_instrument_manifest()
        self.update_financing_schedule()
        # self.update_PPE_projections()
        self.update_depreciation_projections()
        # self.update_agent_financial_statements()

    def update_capex_projections(self):
        # Based on the current status of any WIP projects, update projections
        #   of capital expenditures for upcoming periods

        # Retrieve a list of all ongoing WIP projects
        WIP_projects = pd.read_sql_query(
            f"SELECT WIP_projects.*, assets.* " +
            "FROM WIP_projects " +
            "INNER JOIN assets " +
            "ON WIP_projects.asset_id = assets.asset_id " +
            "WHERE " +
            f"assets.completion_pd > {self.current_pd} " +
            f"AND assets.retirement_pd > {self.current_pd} " +
            f"AND assets.cancellation_pd > {self.current_pd} " +
            f"AND WIP_projects.period = {self.current_pd}",
            self.db)

        # Create dataframe to hold all new capex_projections entries
        capex_cols = [
            "agent_id",
            "asset_id",
            "base_pd",
            "projected_pd",
            "capex"]
        capex_updates = pd.DataFrame(columns=capex_cols)

        # Iterate through all WIP projects and project capex through the
        #   project's end
        for row in WIP_projects.itertuples():
            asset_id = getattr(row, "asset_id")
            agent_id = getattr(row, "agent_id")

            projected_capex = self.project_capex(
                getattr(row, "unit_type"),
                getattr(row, "cum_construction_exp"),
                getattr(row, "cum_construction_duration"),
                getattr(row, "rcec"),
                getattr(row, "rtec")
            )

            for i in range(len(projected_capex)):
                new_row = [agent_id,
                           asset_id,
                           self.current_pd,
                           self.current_pd + i - 1,
                           projected_capex[i]]
                capex_updates.loc[len(capex_updates.index)] = new_row

        # Write the entire capex_cols dataframe to the capex_projections
        #   DB table
        capex_updates.to_sql(
            "capex_projections",
            self.db,
            if_exists="append",
            index=False)

        self.db.commit()

    def project_capex(self, unit_type, cum_construction_exp, cum_construction_duration, rcec, rtec):
        # For a given WIP project, project the sequence of annual capital
        #   expenditures until the project's expected completion
        # TODO: update with more specific methods for different project types
        # For now, assume the project proceeds linearly
        projected_capex = []
        for i in range(math.ceil(round(rtec + 1, 3))):
            projected_capex.append(rcec / rtec)

        return projected_capex


    def initialize_preexisting_instruments(self, fin_insts_updates):
        # Initialize the instrument ids
        inst_id = self.settings["financing"]["starting_instrument_id"]

        for agent_id, agent in self.agents.items():
            if not agent.inactive:
                # Compute starting level of extant equity
                if agent.debt_fraction != 0:
                    starting_equity = float(agent.starting_debt) / agent.debt_fraction * (1 - agent.debt_fraction)
                else:
                    starting_equity = agent.starting_debt

                # Instantiate a debt record
                debt_row = [agent.unique_id,    # agent_id
                            inst_id,   # instrument_id
                            "debt",             # instrument_type
                            agent.unique_id,    # asset_id (agent_id for starting instruments)
                            self.settings["constants"]["time_before_start"],   # pd_issued
                            float(agent.starting_debt),      # initial_principal
                            self.settings["financing"]["default_debt_term"],   # maturity_pd
                            agent.cost_of_debt  # rate
                           ]

                fin_insts_updates.loc[len(fin_insts_updates.index)] = debt_row
                inst_id += 1

                # Instantiate an equity record
                equity_row = [agent.unique_id,
                              inst_id,
                              "equity",
                              agent.unique_id,
                              self.settings["constants"]["time_before_start"],
                              starting_equity,
                              self.settings["financing"]["default_equity_horizon"],
                              agent.cost_of_equity
                             ]

                fin_insts_updates.loc[len(fin_insts_updates.index)] = equity_row
                inst_id += 1

        return fin_insts_updates


    def update_financial_instrument_manifest(self):
        # Based on projected capital expenditures, project the total set of
        #   active financial instruments which will exist at any point in the
        #   future
        # Assumptions:
        #   - the agent tries to maintain a fixed debt/equity ratio by
        #       amortizing all instruments with a sinking fund at their cost
        #   - the "maturity" life of all instruments is 30 years

        # Pull out all financial instruments which already exist (as of
        #   last period). All other entries will be overwritten
        fin_insts_updates = pd.read_sql_query(
            f"SELECT * FROM financial_instrument_manifest WHERE pd_issued < {self.current_pd - 1}",
            self.db)

        self.db.commit()

        # Delete all other contents of this table
        self.db.cursor().execute("DELETE FROM financial_instrument_manifest")

        # On the first period, add instruments representing preexisting
        #   debt and equity for the agents
        if self.current_pd == 0:
            fin_insts_updates = self.initialize_preexisting_instruments(fin_insts_updates)

        # Get a list of all capex projections
        new_capex_instances = pd.read_sql_query(
            f"SELECT * FROM capex_projections WHERE base_pd >= {self.current_pd} AND projected_pd >= {self.current_pd-1}",
            self.db)
        inst_id = max(fin_insts_updates["instrument_id"]) + 1
        for i in range(len(new_capex_instances.index)):
            agent_id = new_capex_instances.loc[i, "agent_id"]
            asset_id = new_capex_instances.loc[i, "asset_id"]
            total_qty = float(new_capex_instances.loc[i, "capex"])
            pd_issued = new_capex_instances.loc[i, "projected_pd"]
            agent_debt_frac = self.agents[agent_id].debt_fraction
            agent_debt_cost = self.agents[agent_id].cost_of_debt
            agent_equity_cost = self.agents[agent_id].cost_of_equity

            # Set up debt issuance for this project for this year
            debt_row = [agent_id,                         # agent_id
                        inst_id,                          # instrument_id
                        "debt",                           # instrument_type
                        asset_id,                         # asset_id
                        pd_issued,                        # pd_issued
                        total_qty * agent_debt_frac,      # initial_principal
                        pd_issued + self.settings["financing"]["default_debt_term"],             # maturity_pd
                        agent_debt_cost                   # rate
                        ]
            fin_insts_updates.loc[len(fin_insts_updates.index)] = debt_row


            # Set up equity issuance for this project for this year
            equity_row = [agent_id,
                          inst_id + 1,
                          "equity",
                          asset_id,
                          pd_issued,
                          total_qty * (1 - agent_debt_frac),
                          pd_issued + self.settings["financing"]["default_equity_horizon"],
                          agent_equity_cost
                          ]
            fin_insts_updates.loc[len(fin_insts_updates.index)] = equity_row

            inst_id = max(fin_insts_updates["instrument_id"]) + 1

        # Overwrite the financial_instrument_manifest table with the new data
        fin_insts_updates.to_sql(
            "financial_instrument_manifest",
            self.db,
            if_exists="append",
            index=False)
        self.db.commit()

    def update_financing_schedule(self):
        # Retrieve the current list of all forecasted financial instruments
        #   which are not past their maturity date
        all_fin_insts = pd.read_sql_query(
            f"SELECT * FROM financial_instrument_manifest WHERE maturity_pd > {self.current_pd}",
            self.db)

        fin_sched_cols = [
            "instrument_id",
            "agent_id",
            "base_pd",
            "projected_pd",
            "total_payment",
            "interest_payment",
            "principal_payment"]
        fin_sched_updates = pd.DataFrame(columns=fin_sched_cols)

        for i in range(len(all_fin_insts.index)):
            inst_id = all_fin_insts.loc[i, "instrument_id"]
            agent_id = all_fin_insts.loc[i, "agent_id"]
            pd_issued = all_fin_insts.loc[i, "pd_issued"]
            initial_principal = all_fin_insts.loc[i, "initial_principal"]
            rate = all_fin_insts.loc[i, "rate"]
            maturity_pd = all_fin_insts.loc[i, "maturity_pd"]
            total_payment = rate * initial_principal / \
                (1 - (1 + rate) ** (-1 * (maturity_pd - pd_issued)))
            remaining_principal = initial_principal
            for projected_pd in range(pd_issued, maturity_pd):
                interest_payment = rate * remaining_principal
                principal_payment = total_payment - interest_payment

                if projected_pd >= self.current_pd:
                    # Organize data and add to fin_sched_updates
                    new_row = [inst_id,
                               agent_id,
                               self.current_pd,
                               projected_pd,
                               total_payment,
                               interest_payment,
                               principal_payment
                               ]
                    fin_sched_updates.loc[len(
                        fin_sched_updates.index)] = new_row

                # Update the amount of remaining principal
                remaining_principal = remaining_principal - principal_payment

        fin_sched_updates.to_sql(
            "financing_schedule",
            self.db,
            if_exists="append",
            index=False)


    def update_depreciation_projections(self):
        # Currently uses straight-line depreciation only
        # Future: allow user selection of SLD or DDBD
        if self.current_pd == 0:
            # Add depreciation schedule for initial PPE for each agent
            dep_cols = [
                "agent_id",
                "asset_id",
                "completion_pd",
                "base_pd",
                "projected_pd",
                "depreciation",
                "beginning_book_value"]
            dep_projections = pd.DataFrame(columns=dep_cols)
            for agent_id, agent in self.agents.items():
                summary_asset_id = agent.unique_id
                init_PPE = agent.starting_PPE
                dep_horiz = self.settings["financing"]["depreciation_horizon"]
                pd_dep = init_PPE / dep_horiz
                for i in range(dep_horiz):
                    beginning_book_value = init_PPE * \
                        (dep_horiz - i) / dep_horiz
                    new_row = [
                        agent_id,
                        summary_asset_id,
                        0,
                        self.current_pd,
                        i,
                        pd_dep,
                        beginning_book_value]
                    dep_projections.loc[len(dep_projections.index)] = new_row
        else:
            # Start by copying forward all entries related to previously-
            #   completed assets (i.e. where total capex is not in question,
            #   so depreciation values are fixed)
            # From last period, copy forward all entries where:
            #   - the asset was complete as of last period at the latest
            #   - the projected period is this period or later
            # Then simply change all `base_pd` values to the current period
            dep_projections = pd.read_sql_query(
                f"SELECT * FROM depreciation_projections WHERE base_pd = {self.current_pd-1} AND completion_pd < {self.current_pd} AND projected_pd >= {self.current_pd}",
                self.db)
            dep_projections["base_pd"] = self.current_pd

            # Then, recompute expected depreciation schedules for all relevant
            #   WIP projects, including:
            #   - WIPs which are ongoing
            #   - new WIPs since last period
            # WIPs completed last period are NOT included in this section
            WIP_projects = pd.read_sql_query(
                f"SELECT * FROM WIP_projects WHERE period = {self.current_pd-1} AND rtec > 0",
                self.db)

            for row in WIP_projects.itertuples():
                asset_id = getattr(row, "asset_id")
                agent_id = getattr(row, "agent_id")
                starting_pd = getattr(row, "period") + \
                    math.ceil(round(getattr(row, "rtec"), 3))
                asset_PPE = getattr(row, "cum_construction_exp") + getattr(row, "rcec")
                dep_horiz = self.settings["financing"]["depreciation_horizon"]
                pd_dep = asset_PPE / dep_horiz
                for i in range(dep_horiz):
                    book_value = asset_PPE * (dep_horiz - i) / dep_horiz
                    new_row = [
                        agent_id,
                        asset_id,
                        starting_pd,
                        self.current_pd,
                        starting_pd + i,
                        pd_dep,
                        book_value]
                    dep_projections.loc[len(dep_projections.index)] = new_row

        dep_projections.to_sql(
            "depreciation_projections",
            self.db,
            if_exists="append",
            index=False)
        self.db.commit()


    def save_ALEAF_outputs(self):
        # Copy all ALEAF output files to the output directory, with
        #   scenario and step-specific names
        files_to_save = [
            "dispatch_summary_OP",
            "expansion_result",
            "system_summary_OP",
            "system_tech_summary_OP"]
        for outfile in files_to_save:
            old_filename = f"{self.settings['simulation']['ALEAF_scenario_name']}__{outfile}.csv"
            old_filepath = Path(self.ALEAF_output_data_path) / old_filename
            new_filename = f"{self.settings['simulation']['ALEAF_scenario_name']}__{outfile}__step_{self.current_pd}.csv"
            new_filepath = Path(self.ABCE_output_data_path) / new_filename
            shutil.copy2(old_filepath, new_filepath)


    def process_ALEAF_dispatch_results(self):
        # Find and load the ALEAF dispatch file corresponding to the current
        #   simulation period
        fname_pattern = f"dispatch_summary_OP__step_{self.current_pd}"
        ALEAF_dsp_file = None
        for fname in os.listdir(Path(self.ABCE_output_data_path)):
            if "dispatch" in fname and "OP" in fname and f"{self.current_pd}" in fname:
                ALEAF_dsp_file = Path(self.ABCE_output_data_path) / fname

        # Get the number of units which are currently operational (needed for
        #   scaling dispatch results to a per-unit basis)
        sql_query = (f"SELECT unit_type, COUNT(unit_type) FROM assets " +
                         f"WHERE completion_pd <= {self.current_pd} " +
                         f"AND retirement_pd > {self.current_pd} " +
                         f"AND cancellation_pd > {self.current_pd} " +
                         f"GROUP BY unit_type")
        num_units = pd.read_sql_query(sql_query, self.db, index_col="unit_type").rename(columns={"COUNT(unit_type)": "num_units"})

        # Postprocess ALEAF dispatch results
        ALEAF_dsp_results = dsp.postprocess_dispatch(ALEAF_dsp_file, num_units, self.unit_specs)
        ALEAF_dsp_results["period"] = self.current_pd
        ALEAF_dsp_results = ALEAF_dsp_results.reset_index().rename(columns={"index": "unit_type"})

        # Get list of column names for ordering
        cursor = self.db.cursor().execute("SELECT * FROM ALEAF_dispatch_results")
        col_names = [description[0] for description in cursor.description]

        # Reorder ALEAF_dsp_results to match database
        ALEAF_dsp_results = ALEAF_dsp_results[col_names]

        logging.debug(ALEAF_dsp_results)

        ALEAF_dsp_results.to_sql("ALEAF_dispatch_results", self.db, if_exists="append", index=False)


    def update_WIP_projects(self):
        """
        Update the status and projections-to-completion for all current WIP
        projects.

        This function iterates over all WIP projects which are currently
        ongoing (not cancelled or completed).
          1. A stochastic escalation generator extends the projected cost
               and/or schedule to completion according to the project type.
          2. The current ANPE amount is applied to the outstanding RCEC.
          3. If this is enough to complete the project, the project is marked
               as completed ('completion_pd' in the 'assets' table is set to
               the current period).
          4. The new project status, including updated RCEC and RTEC, is saved
               to the 'WIP_projects' table.
        """

        # Get a list of all currently-active construction projects
        if self.current_pd == 0:
            WIP_projects = pd.read_sql(
                "SELECT asset_id FROM assets WHERE " +
                f"completion_pd > {self.current_pd} AND " +
                f"cancellation_pd >= {self.current_pd}",
                self.db)
        else:
            WIP_projects = pd.read_sql(
                "SELECT asset_id FROM assets WHERE " +
                f"completion_pd >= {self.current_pd} AND " +
                f"cancellation_pd > {self.current_pd}",
                self.db)

        # Update each project one at a time
        for asset_id in WIP_projects.asset_id:
            # Select this project's most recent data record
            project_data = pd.read_sql_query(
                f"SELECT * FROM WIP_projects WHERE asset_id = {asset_id} AND period = {self.current_pd-1}",
                self.db)

            # Record the effects of authorized construction expenditures, and
            #   advance the time-remaining estimate by one year
            project_data = self.advance_project_to_current_period(project_data)

            # If this period's authorized expenditures (ANPE) clear the RCEC,
            #   then the project is complete
            if project_data.loc[0, "rcec"] <= self.settings["constants"]["large_epsilon"]:
                # Record the project's completion period as the current period
                self.record_completed_construction_project(project_data)

            # Record updates to the WIP project's status
            self.record_WIP_project_updates(project_data)

            # Record updates to the project's expected completion date in the
            #   database 'assets' table
            self.update_expected_completion_period(project_data)


    def record_completed_construction_project(self, project_data):
        asset_id = project_data.loc[0, "asset_id"]

        # Get asset record from assets
        asset_data = pd.read_sql_query(
            f"SELECT * FROM assets WHERE asset_id = {asset_id}", self.db)

       # Compute periodic sinking fund payments
        unit_type = asset_data.loc[0, "unit_type"]
        unit_life = int(math.ceil(self.unit_specs[unit_type]["unit_life"]))
        capex_payment = 0

        to_update = {"completion_pd": self.current_pd,
                     "retirement_pd": self.current_pd + unit_life,
                     "total_capex": project_data.loc[0, "cum_construction_exp"],
                     "cap_pmt": capex_payment}
        filters = {"asset_id": asset_id}

        ABCE.update_DB_table_inplace(
            self.db,
            self.cur,
            "assets",
            to_update,
            filters
        )

        # Commit changes to database
        self.db.commit()

    def record_WIP_project_updates(self, project_data):
        # Set new_anpe = 0 to avoid inter-period contamination
        #new_anpe = 0
        # Update the 'WIP_projects' table with new RCEC/RTEC/ANPE values
        self.cur.execute("INSERT INTO WIP_projects VALUES " +
                         f"({project_data.asset_id.values[0]}, " +
                         f"{project_data.agent_id.values[0]}, " +
                         f"{self.current_pd}, " +
                         f"{project_data.cum_occ.values[0]}, " +
                         f"{project_data.rcec.values[0]}, " +
                         f"{project_data.cum_construction_duration.values[0]}, " +
                         f"{project_data.rtec.values[0]}, " +
                         f"{project_data.cum_construction_exp.values[0]}, " +
                         f"{project_data.anpe.values[0]})")

        # Save changes to the database
        self.db.commit()

    def compute_total_capex_newbuild(self, asset_id):
        """
        For assets which are newly completed during any period except period 0:
        Retrieve all previously-recorded capital expenditures (via 'anpe', i.e.
        "Authorized Next-Period Expenditures") for the indicated asset ID from
        the database, and sum them to return the total capital expenditure (up
        to the current period).

        Args:
          asset_id (int): asset for which to compute total capex

        Returns:
          total_capex (float): total capital expenditures up to the present period
        """

        total_capex = pd.read_sql("SELECT SUM(anpe) FROM WIP_projects WHERE " +
                                  f"asset_id = {asset_id}", self.db).iloc[0, 0]
        return total_capex

    def compute_sinking_fund_payment(self, agent_id, total_capex, term):
        """
        Compute a constant sinking-fund payment based on a total capital-
        expenditures amount and the life of the investment, using the specified
        agent's financial parameters.

        Args:
          agent_id (int): the unique ID of the owning agent, to retrieve
            financial data
          total_capex (float): total capital expenditures on the project
          term (int or float): term over which to amortize capex

        Returns:
          cap_pmt (float): equal capital repayments to make over the course
            of the indicated amortization term
        """

        agent_params = pd.read_sql("SELECT * FROM agent_params WHERE " +
                                   f"agent_id = {agent_id}", self.db)

        wacc = (agent_params.debt_fraction * agent_params.cost_of_debt +
                (1 - agent_params.debt_fraction) * agent_params.cost_of_equity)
        cap_pmt = total_capex * wacc / (1 - (1 + wacc)**(-term))

        # Convert cap_pmt from a single-entry Series to a scalar
        cap_pmt = cap_pmt[0]

        return cap_pmt

    def advance_project_to_current_period(self, project_data):
        # Escalated RTEC is reduced by one
        project_data.loc[0, "rtec"] -= 1

        # Escalated RCEC is reduced by ANPE
        project_data.loc[0, "rcec"] -= project_data.loc[0, "anpe"]

        # Cumulative capital expenditures are increased by ANPE
        project_data.loc[0, "cum_construction_exp"] += project_data.loc[0, "anpe"]

        # Project ANPE is reset to 0 after being expended
        #project_data.loc[0, "anpe"] = 0

        return project_data

    def update_expected_completion_period(self, project_data):
        asset_id = project_data.loc[0, "asset_id"]
        new_completion_pd = project_data.loc[0, "rtec"] + self.current_pd

        ABCE.update_DB_table_inplace(
            self.db,
            self.cur,
            "assets",
            {"completion_pd": new_completion_pd},
            {"asset_id": asset_id}
        )

        self.db.commit()

    def execute_all_status_updates(self):
        # Record newly-started WIP projects from the agents' decisions
        WIP_updates = pd.read_sql_query("SELECT * FROM WIP_updates", self.db)
        WIP_updates.to_sql(
            "WIP_projects",
            self.db,
            if_exists="append",
            index=False)
        self.db.commit()

        # Record newly-started C2N WIP projects from the agents' decisions
        C2N_updates = pd.read_sql_query(
            "SELECT * FROM WIP_C2N_updates", self.db)
        C2N_updates.to_sql("WIP_C2N", self.db, if_exists="append", index=False)

        # Record status updates to existing assets (i.e. retirements)
        # Convert asset_updates to a dict of dicts for convenience
        asset_updates = pd.read_sql_query(
            "SELECT * FROM asset_updates", self.db)
        for row_num in asset_updates.index:
            new_record = asset_updates.loc[[
                row_num]].copy().reset_index(drop=True)

            orig_record = pd.read_sql_query(
                f"SELECT * FROM assets WHERE asset_id = {new_record.loc[0, 'asset_id']}", self.db)

            if len(orig_record) == 0:
                # The asset does not already exist and an entry must be added
                # Ensure that completion and retirement periods are integers
                new_record.at[0, "completion_pd"] = int(
                    math.ceil(new_record.at[0, "completion_pd"]))
                new_record.at[0, "retirement_pd"] = int(
                    math.ceil(new_record.at[0, "retirement_pd"]))
                pd.DataFrame(new_record).to_sql(
                    "assets", self.db, if_exists="append", index=False)
            else:
                # The prior record must be overwritten
                # Move the filtering data (asset and agent ids) into a separate
                #   dictionary
                new_record = new_record.to_dict(orient="records")[0]
                filters = {}
                col_filters = ["asset_id", "agent_id"]
                for col in col_filters:
                    filters[col] = new_record.pop(col)

                # Update the record in the DB table
                ABCE.update_DB_table_inplace(
                    self.db,
                    self.cur,
                    "assets",
                    new_record,
                    filters
                )
        self.db.commit()

        # Delete all contents of the WIP_updates and asset_updates tables
        self.db.cursor().execute("DELETE FROM WIP_updates")
        self.db.cursor().execute("DELETE FROM asset_updates")
        self.db.commit()
