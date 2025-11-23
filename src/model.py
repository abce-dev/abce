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

# import local modules
from .agent import GenCo
from . import ABCEfunctions as ABCE
from . import seed_creator as sc
from . import dispatch_ppx as dsp
from . import input_data_management as idm



class GridModel(Model):
    """A model with some number of GenCos."""

    def __init__(self, settings, args):
        # Copy the command-line arguments as member data
        self.args = args
        self.settings = settings

        # If verbosity is 2 or 3, show the ABCE splash header
        if self.args.verbosity >= 2:
            self.show_abce_header()

        # Check the outputs dir and clear out old files
        self.prepare_outputs_directory()

        # Initialize the database in the outputs directory
        self.initialize_database()

        # Initialize all data
        self.load_all_data()

        # Ensure a tmp directory exists inside the current working directory
        if self.args.verbosity > 2:
            self.ensure_tmp_dir_exists()

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Create agents
        for agent_id, agent_params in self.agent_specs.items():
            agent = GenCo(agent_id, self, agent_params)

            # Initialize assets owned by this agent, and add them to
            #   the database
            self.initialize_agent_assets(agent)

            # Add this agent to the dictionary of agents
            self.agents[agent_id] = agent

            # Add this agent to the schedule
            self.schedule.add(agent)

        # Make sure all pending database transactions are resolved
        self.db.commit()

        # If the simulation contains a balance-of-system agent, set up its
        #   expansion schedule
        for agent, agent_data in self.agents.items():
            if agent_data.balance_of_system:
                self.set_up_BOS_expansion(agent, agent_data)


    def load_all_data(self):
        # Retrieve the input data
        self.unit_specs = idm.initialize_unit_specs(self.settings, self.args)

        # Save unit_specs to the database
        self.add_unit_specs_to_db()

        # Load all-period demand data into the database
        self.load_demand_data_to_db()

        # Add model parameters to the database
        self.load_model_parameters_to_db()

        # Read agent specification data
        self.load_agent_specifications()

        # Initialize the model one time step before the true start date
        self.current_pd = self.settings["constants"]["time_before_start"]

        # Define a dictionary to hold all agents for easier indexing
        self.agents = {}


    def load_agent_specifications(self):
        # Read in the GenCo parameters data from file
        agent_specs_file_name = Path(
            Path(self.args.inputs_path)
            / self.settings["file_paths"]["agent_specifications_file"]
        )
        self.agent_specs = yaml.load(
            open(agent_specs_file_name, "r"), Loader=yaml.FullLoader
        )

    def prepare_outputs_directory(self):
        # Set up the primary output data path
        self.primary_output_data_path = (
            Path(os.getcwd())
            / "outputs"
            / self.settings["simulation"]["scenario_name"]
        )

        # Set up the output data path for logging and diagnostics
        if self.settings["file_paths"]["output_logging_dir"] == "tmp":
            self.logging_output_data_path = (
                Path(os.getcwd())
                / "tmp"
                / self.settings["simulation"]["scenario_name"]
            )
        else:
            self.logging_output_data_path = (
                Path(self.settings["file_paths"]["output_logging_dir"])
                / self.settings["simulation"]["scenario_name"]
            )

        output_paths = [Path(self.primary_output_data_path), Path(self.logging_output_data_path)]

        for outpath in output_paths:
            if not outpath.is_dir():
                # If the desired output directory doesn't already exist, create it
                outpath.mkdir(exist_ok=True, parents=True)
            else:
                # Otherwise, delete any existing files and subdirectories 
                # in the directory
                for root, dirs, files in os.walk(outpath, topdown=False):
                    for name in files:
                        (Path(root) / name).unlink()
                    for name in dirs:
                        (Path(root) / name).rmdir()

        # Save input data files to the output directory
        # settings.yml (dump yml contents from memory to file)
        with open(Path(self.primary_output_data_path) / "settings.yml", "w") as s_file:
            yaml.dump(self.settings, s_file)

        # inputs/agent_specification.yml (copy file directly)
        as_fname = self.settings["file_paths"]["agent_specifications_file"]
        shutil.copyfile(
            Path(self.args.inputs_path) / as_fname,
            Path(self.primary_output_data_path) / as_fname,
        )

        # inputs/demand_data.yml (copy file directly)
        dd_fname = self.settings["file_paths"]["demand_data_file"]
        shutil.copyfile(
            Path(self.args.inputs_path) / dd_fname,
            Path(self.primary_output_data_path) / dd_fname,
        )

        # inputs/unit_specs.yml (copy file directly)
        us_fname = self.settings["file_paths"]["unit_specs_data_file"]
        shutil.copyfile(
            Path(self.args.inputs_path) / us_fname,
            Path(self.primary_output_data_path) / us_fname,
        )


    def initialize_database(self):
        # Initialize database for storing and managing all simulation data
        self.db_file = (
            self.primary_output_data_path
            / self.settings["file_paths"]["db_file"]
        )
        self.db, self.cur = sc.create_database(self.db_file, self.args.force)


    def ensure_tmp_dir_exists(self):
        # Create the local tmp/ directory inside the current working directory,
        #   if it doesn't already exist
        tmp_dir_location = Path.cwd() / "tmp"
        Path(tmp_dir_location).mkdir(exist_ok=True)


    def show_abce_header(self):
        logo_file = Path(
            self.settings["file_paths"]["ABCE_abs_path"]
            / "misc"
            / self.settings["file_paths"]["logo"]
        )
        with open(logo_file, "r") as logo:
            for line in logo.read().splitlines():
                logging.log(self.settings["constants"]["vis_lvl"], line)

    def add_unit_specs_to_db(self):
        """
        This function selects only the unit_specs columns needed for the
          database, then saves the table to the database.
        """
        # Reshape unit_specs into a dataframe
        unit_specs = (
            pd.DataFrame.from_dict(self.unit_specs, orient="index")
            .reset_index()
            .rename(columns={"index": "unit_type"})
        )

        # Retrieve the column-header schema for the 'unit_specs' table
        self.cur.execute("SELECT * FROM unit_specs")
        unit_specs_db_cols = [element[0] for element in self.cur.description]

        # Select the necessary DB columns in order
        unit_specs = unit_specs[unit_specs_db_cols]

        # Save the finalized unit specs data to the DB, and set the member data
        unit_specs.to_sql(
            "unit_specs", self.db, if_exists="replace", index=False
        )


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
            retirements = agent.scheduled_retirements[unit_type]
            for ret_pd, num_ret_units in retirements.items():
                total_ret_units += num_ret_units

            # Reconcile total number of scheduled retirements with total number
            #   of owned assets of this type
            if total_ret_units > num_units:
                raise ValueError(
                    f"Portfolio specification mismatch for agent "
                    + f"#{agent.unique_id}: total owned {unit_type} units = "
                    + f"{num_units}, but total specified {unit_type} "
                    + f"retirements = {total_ret_units}."
                )

            # If there are any units with unspecified retirement dates, set
            #   their retirement date to a large number
            elif total_ret_units < num_units:
                dist_t = self.settings["constants"]["distant_time"]
                agent.scheduled_retirements[unit_type][dist_t] = (
                    num_units - total_ret_units
                )


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
                asset_dict = {
                    "agent_id": agent.unique_id,
                    "unit_type": unit_type,
                    "start_pd": self.settings["constants"]["time_before_start"],
                    "completion_pd": 0,
                    "cancellation_pd": self.settings["constants"][
                        "distant_time"
                    ],
                    "retirement_pd": int(retirement_pd),
                    "total_capex": unit_capex,
                    "cap_pmt": cap_pmt,
                    "C2N_reserved": 0,
                }

                # For each asset in this block, create a dataframe record and
                #   store it to the master_assets_df
                for i in range(num_units):
                    # Find the largest extant asset id, and set the current
                    #   asset id 1 higher
                    asset_dict["asset_id"] = max(
                        ABCE.get_next_asset_id(
                            self.db,
                            self.settings["constants"]["first_asset_id"],
                        ),
                        max(
                            master_assets_df["asset_id"],
                            default=self.settings["constants"][
                                "first_asset_id"
                            ],
                        )
                        + 1,
                    )

                    # Convert the dictionary to a dataframe format and save
                    new_record = pd.DataFrame(asset_dict, index=[0])
                    master_assets_df = pd.concat([master_assets_df, new_record])

        # Once all assets from all unit types for this agent have had records
        #   initialized, save the dataframe of all assets into the 'assets'
        #   DB table
        master_assets_df.to_sql(
            "assets", self.db, if_exists="append", index=False
        )

        self.db.commit()


    def set_up_BOS_expansion(self, agent, agent_data):
        if agent_data.expansion_strategy == "proportional_expansion":
            pf = pd.read_sql_query(f"SELECT * FROM assets WHERE agent_id = {agent}", self.db)
            pf["agent_id"] = agent

            sys_pf = pd.read_sql_query(f"SELECT * FROM assets", self.db)
            sys_pf["agent_id"] = "system"

            pf = pd.concat([pf, sys_pf])
 
            unit_specs = pd.read_sql_query("SELECT * FROM unit_specs", self.db)

            pf = pf.merge(
                unit_specs[["unit_type", "capacity", "capacity_factor"]],
                on = "unit_type",
            )

            pf["der_cap"] = pf["capacity"] * pf ["capacity_factor"]

            pf = pd.pivot_table(
                pf,
                values = "der_cap",
                index="unit_type",
                columns="agent_id",
                aggfunc="sum",
            ).fillna(0)

            # Re-add the per-unit derated capacity column
            pf = pf.merge(
                unit_specs[["unit_type", "capacity", "capacity_factor"]],
                on = "unit_type",
            )
            pf["der_cap"] = pf["capacity"] * pf ["capacity_factor"]
            pf = pf.set_index("unit_type")
            
            pf = pf[[agent, "system", "der_cap"]]

            pf["agent_pct"] = pf[agent] / sum(pf[agent])
            pf["carryover"] = 0.0

            # Compute the balance-of-system agent's starting market share
            MS_init = sum(pf[agent]) / sum(pf["system"])

            # Retrieve all demand data
            demand = pd.read_sql_query("SELECT * FROM demand", self.db)

            # Get the column schema for the assets table
            cols = sys_pf.columns.tolist()

            new_units = pf[[agent]].copy(deep=True)

            # Add units each year of the demand forecast horizon
            for y in range(min(demand["period"])+1, max(demand["period"])+1):
                # Calculate absolute change in demand for this year
                demand_delta = demand.loc[y, "demand"] - demand.loc[y-1, "demand"]

                # Calculate change in capacity this agent should take care of
                cap_delta = demand_delta * MS_init

                # Total new capacity this agent should add
                total_new = cap_delta * pf["agent_pct"] / pf["der_cap"] + pf["carryover"]

                # Record full new units
                new_units[y] = np.floor(total_new)

                # Update unit carryover
                pf["carryover"] = total_new % 1

            new_units = new_units.drop(columns=[agent])

            added_assets_df = pd.DataFrame(columns=cols)

            # Set up all data values except for the asset id and unit type
            asset_dict = {
                "agent_id": agent,
                "cancellation_pd": self.settings["constants"]["distant_time"],
                "retirement_pd": self.settings["constants"]["distant_time"],
                "total_capex": 0,
                "cap_pmt": 0,
                "C2N_reserved": 0,
            }

            # For each year in the dataframe:
            for y in range(1, max(new_units.columns.tolist())+1):
                # For each unit type in the dataframe:
                for unit_type in new_units.index.values.tolist():
                    if new_units.loc[unit_type, y] != 0:
                        # Create a dataframe record and store it in the
                        #   added_assets_df dataframe
                        asset_dict["unit_type"] = unit_type
                        asset_dict["start_pd"] = y
                        asset_dict["completion_pd"] = y

                        # Find the largest extant asset id, and set the current
                        #   asset id 1 higher
                        asset_dict["asset_id"] = max(
                            ABCE.get_next_asset_id(
                                self.db,
                                self.settings["constants"]["first_asset_id"],
                            ),
                            max(
                                added_assets_df["asset_id"],
                                default=self.settings["constants"][
                                    "first_asset_id"
                                ],
                            )
                            + 1,
                        )

                        # Convert the dictionary to a dataframe format and save
                        new_record = pd.DataFrame(asset_dict, index=[0])
                        added_assets_df = pd.concat([added_assets_df, new_record])

            # Once all additional assets planned for this agent have had records
            #   initialized, save the dataframe of all assets into the 'assets'
            #   DB table
            added_assets_df.to_sql(
                "assets", self.db, if_exists="append", index=False
            )

            self.db.commit()


    def compute_total_capex_preexisting(self, unit_type):
        unit_cost_per_kW = self.unit_specs[unit_type]["overnight_capital_cost"]
        unit_capacity = self.unit_specs[unit_type]["capacity"]

        total_capex = (
            unit_cost_per_kW
            * unit_capacity
            * self.settings["constants"]["MW2kW"]
        )

        return total_capex


    def load_demand_data_to_db(self):
        # Load all-period demand data into the database
        demand_data_file = (
            Path(self.args.inputs_path)
            / self.settings["file_paths"]["demand_data_file"]
        )
        demand_df = (
            pd.read_csv(demand_data_file)
            * self.settings["scenario"]["peak_demand"]
        )

        num_steps = self.settings["simulation"]["num_steps"]
        dem_horiz = self.settings["demand"]["demand_visibility_horizon"]

        if len(demand_df) < num_steps + dem_horiz:
            msg = (
                f"The peak demand time series data supplied in " +
                f"{demand_data_file} is too short for the specified number " +
                f"of simulation years to be run.\n" +
                f"  Number of years given: {len(demand_df)}.\n" +
                f"  Number of years needed: {num_steps + dem_horiz} = " +
                f"{num_steps} simulated year(s) + {dem_horiz} year(s) " +
                f"lookahead.\n"
                f"Supply additional years of peak demand data, or run ABCE " +
                f"for fewer time-steps."
            )
            logging.error(msg)
            exit()

       # Save data to DB
        demand_df.to_sql(
            "demand", self.db, if_exists="replace", index_label="period"
        )
        self.db.commit()

    def load_model_parameters_to_db(self):
        # Load specific parameters specified in settings.yml to the database
        prm = self.settings["system"]["planning_reserve_margin"]
        self.cur.execute(f"INSERT INTO model_params VALUES ('PRM', {prm})")


        tax_rate = self.settings["system"]["tax_rate"]
        self.cur.execute(
            f"INSERT INTO model_params VALUES ('tax_rate', {tax_rate})"
        )

        tax_credits_discount = self.settings["system"]["tax_credits_discount"]
        self.cur.execute(f"INSERT INTO model_params VALUES ('tax_credits_discount', {tax_credits_discount})")

        self.db.commit()

    def step(self, demo=False):
        """
        Advance the model by one step.
        """
        self.current_pd += 1

        # Perform some step-0 setup
        if self.current_pd == 0:
            self.has_ABCE_sysimage = self.check_for_sysimage_file()

        self.display_step_header()

        # Advance the status of all WIP projects to the current period
        self.update_WIP_projects()

        # If any agent will have an empty portfolio this year, delete them
        #   from the schedule
        units_this_year = pd.read_sql_query(
            f"SELECT agent_id, COUNT(asset_id) FROM assets WHERE completion_pd <= {self.current_pd} AND retirement_pd > {self.current_pd} GROUP BY agent_id",
            self.db
        )
        for agent_id, agent in self.agents.items():
            if (str(agent_id) not in units_this_year.agent_id.values) and (self.agents[agent_id] in self.schedule.agents):
                logging.info(f"Removing agent {agent_id} from the simulation due to a size-zero portfolio.")
                self.schedule.remove(self.agents[agent_id])

        # Compute the scenario reduction results for the dispatch forecast
        #   horizon, starting with this year
        fc_horiz = self.settings["dispatch"]["num_dispatch_years"]
        for fc_y in range(self.current_pd, self.current_pd + fc_horiz):
            ABCE.execute_scenario_reduction(
                self.args,
                self.db,
                self.current_pd,
                fc_y,
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
            "\nAll agent turns are complete.\n",
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

        # Set up the command to run dispatch.jl in annual exact mode
        ABCE_ENV = Path(os.environ["ABCE_ENV"])
        annual_disp_script_path = Path(self.settings["file_paths"]["ABCE_abs_path"]) / "src" / "annual_dispatch.jl"
            
        dispatch_cmd = (
            f"julia --project={ABCE_ENV} {annual_disp_script_path} " 
            + f"--ABCE_dir={self.settings['file_paths']['ABCE_abs_path']} "
            + f"--current_pd={self.current_pd} "
            + f"--settings_file={self.args.settings_file} "
        )

        # Run the dispatch simulation
        if self.args.verbosity < 2:
            sp = subprocess.check_call(
                dispatch_cmd, shell=True, stdout=open(os.devnull, "wb")
            )
        else:
            sp = subprocess.check_call(dispatch_cmd, shell=True)


    def display_step_header(self):
        if self.current_pd != 0:
            logging.log(self.settings["constants"]["vis_lvl"], "\n\n\n")
        logging.log(self.settings["constants"]["vis_lvl"], "=" * 60)
        logging.log(
            self.settings["constants"]["vis_lvl"],
            f"   Simulation step: {self.current_pd}",
        )
        logging.log(self.settings["constants"]["vis_lvl"], "=" * 60)


    def show_round_updates(self):
        logging.debug("Updates to assets:")
        logging.debug(pd.read_sql("SELECT * FROM asset_updates", self.db))
        logging.debug("\nConstruction project updates (last 10 entries):")
        logging.debug(
            pd.read_sql("SELECT * FROM WIP_updates", self.db).tail(n=10)
        )


    def check_for_sysimage_file(self):
        ABCE_sysimage_path = (
            Path(os.getenv("ABCE_ENV"))
            / self.settings["file_paths"]["ABCE_sysimage_file"]
        )

        has_ABCE_sysimage = True

        if not Path(ABCE_sysimage_path).exists():
            msg = (
                f"No sysimage file found at {ABCE_sysimage_path}. "
                + "Execution will proceed, but Julia may run "
                + "slowly. If you already have a sysimage file, please move "
                + f"it to the filename {ABCE_sysimage_path}. If you do not have a "
                + "sysimage file, please run 'julia make_sysimage.jl "
                + "--mode=abce' in this directory."
            )
            logging.warn(msg)
            has_ABCE_sysimage = False

        return has_ABCE_sysimage


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
        # Get a list of any balance-of-system agents to exclude
        BOS_agents = []
        not_in = ""
        for agent in self.agents:
            if self.agents[agent].balance_of_system:
                BOS_agents.append(str(agent))
        if len(BOS_agents) != 0:
            not_in = ",".join(BOS_agents)
            not_in = f"AND agent_id NOT IN ({not_in})"

        # Get a list of all currently-active construction projects
        if self.current_pd == 0:
            WIP_projects = pd.read_sql(
                "SELECT asset_id FROM assets WHERE "
                + f"completion_pd > {self.current_pd} AND "
                + f"cancellation_pd >= {self.current_pd} "
                + not_in,
                self.db,
            )
        else:
            WIP_projects = pd.read_sql(
                "SELECT asset_id FROM assets WHERE "
                + f"completion_pd >= {self.current_pd} AND "
                + f"cancellation_pd > {self.current_pd} "
                + not_in,
                self.db,
            )

        # Update each project one at a time
        for asset_id in WIP_projects.asset_id:
            # Select this project's most recent data record
            project_data = pd.read_sql_query(
                (
                    f"SELECT * FROM WIP_projects "
                    + f"WHERE asset_id = {asset_id} "
                    + f"AND period = {self.current_pd-1}"
                ),
                self.db,
            )

            # Record the effects of authorized construction expenditures, and
            #   advance the time-remaining estimate by one year
            project_data = self.advance_project_to_current_period(project_data)

            # If this period's authorized expenditures (ANPE) clear the RCEC,
            #   then the project is complete
            if (
                project_data.loc[0, "rcec"]
                <= self.settings["constants"]["large_epsilon"]
            ):
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
            f"SELECT * FROM assets WHERE asset_id = {asset_id}", self.db
        )

        # Compute periodic sinking fund payments
        unit_type = asset_data.loc[0, "unit_type"]
        unit_life = int(math.ceil(self.unit_specs[unit_type]["unit_life"]))
        capex_payment = 0

        updates = {
            "completion_pd": self.current_pd,
            "retirement_pd": self.current_pd + unit_life,
            "total_capex": project_data.loc[0, "cum_construction_exp"],
            "cap_pmt": capex_payment,
        }
        filters = {"asset_id": asset_id}

        ABCE.update_DB_table_inplace(
            self.db, self.cur, "assets", updates, filters
        )

        # Commit changes to database
        self.db.commit()

    def record_WIP_project_updates(self, project_data):
        # Update the 'WIP_projects' table with new RCEC/RTEC/ANPE values
        self.cur.execute(
            "INSERT INTO WIP_projects VALUES "
            + f"({project_data.asset_id.values[0]}, "
            + f"{project_data.agent_id.values[0]}, "
            + f"{self.current_pd}, "
            + f"{project_data.cum_occ.values[0]}, "
            + f"{project_data.rcec.values[0]}, "
            + f"{project_data.cum_construction_duration.values[0]}, "
            + f"{project_data.rtec.values[0]}, "
            + f"{project_data.cum_construction_exp.values[0]}, "
            + f"{project_data.anpe.values[0]})"
        )

        # Save changes to the database
        self.db.commit()


    def advance_project_to_current_period(self, project_data):
        # Escalated RTEC is reduced by one
        project_data.loc[0, "rtec"] -= 1

        # Escalated RCEC is reduced by ANPE
        project_data.loc[0, "rcec"] -= project_data.loc[0, "anpe"]

        # Cumulative capital expenditures are increased by ANPE
        project_data.loc[0, "cum_construction_exp"] += project_data.loc[
            0, "anpe"
        ]

        # Project ANPE is reset to 0 after being expended
        # project_data.loc[0, "anpe"] = 0

        return project_data


    def update_expected_completion_period(self, project_data):
        asset_id = project_data.loc[0, "asset_id"]
        new_completion_pd = project_data.loc[0, "rtec"] + self.current_pd

        ABCE.update_DB_table_inplace(
            self.db,
            self.cur,
            "assets",
            {"completion_pd": new_completion_pd},
            {"asset_id": asset_id},
        )

        self.db.commit()


    def convert_WIP_updates_to_FIs(self, WIP_updates):
        """
        This function handles the recording of ALL financial information
          ensuing from selected project alternatives by ALL agents.

        This function is run after all agent decision rounds have been
          completed. This data is saved once and is never altered or deleted.
        """

        # Retrieve the list of all extant financial instruments
        extant_FIs = pd.read_sql_query(
            "SELECT * FROM financial_instrument_manifest",
            self.db,
        )

        # Retrieve the schema for capex forecasts
        extant_capex = pd.read_sql_query(
            "SELECT * FROM capex_projections",
            self.db,
        )
        capex_cols = extant_capex.columns.tolist()

        # Retrieve the schema for the financing schedules
        extant_fin_sched = pd.read_sql_query(
            "SELECT * FROM financing_schedule",
            self.db,
        )
        fin_sched_cols = extant_fin_sched.columns.tolist()

        # Retrieve the schema for depreciation forecasts
        extant_dep = pd.read_sql_query(
            "SELECT * FROM depreciation_projections",
            self.db,
        )
        dep_cols = extant_dep.columns.tolist()

        # Create empty dataframes with the same structure as each table,
        #   to store new financial information
        fin_insts_updates = pd.DataFrame(columns=extant_FIs.columns.tolist())
        capex_updates = pd.DataFrame(columns=capex_cols)
        fin_sched_updates = pd.DataFrame(columns=fin_sched_cols)
        dep_updates = pd.DataFrame(columns=dep_cols)

        # Determine the next available instrument id
        next_id = max(extant_FIs["instrument_id"]) + 1

        for index, row in WIP_updates.iterrows():
            # Compute values that don't change per year
            agent_id = row["agent_id"]
            asset_id = row["asset_id"]
            inst_type = "debt"
            debt_term = self.settings["financing"]["default_debt_term"]
            annual_capex = row["cum_occ"] / row["cum_construction_duration"]
            initial_principal = annual_capex * self.agents[agent_id].debt_fraction
            rate = self.agents[agent_id].cost_of_debt

            total_payment = (
                rate * initial_principal 
                / (1 - (1 + rate) ** (-debt_term))
            )

            # Project out annual depreciation
            dep_horiz = self.settings["financing"]["depreciation_horizon"]
            dep_per_pd = - row["cum_occ"] / dep_horiz
            book_value = row["cum_occ"]

            for i in range(dep_horiz):
                y = i + self.current_pd + row["cum_construction_duration"]

                dep_data = [
                    agent_id,
                    asset_id,
                    y,
                    dep_per_pd,
                    book_value
                ]

                dep_updates.loc[len(dep_updates.index)] = dep_data

                book_value += dep_per_pd

            # Iterate over all construction years for this project
            for i in range(int(row["cum_construction_duration"])):
                pd_issued = i + self.current_pd
                maturity_pd = pd_issued + debt_term - 1

                # Add the instrument data for this year
                inst_data = [
                    agent_id,
                    next_id,
                    inst_type,
                    asset_id,
                    pd_issued,
                    initial_principal,
                    maturity_pd,
                    rate
                ]

                fin_insts_updates.loc[len(fin_insts_updates.index)] = inst_data

                # Add the capex data for this year
                capex_data = [
                    agent_id,
                    asset_id,
                    pd_issued,
                    annual_capex
                ]

                capex_updates.loc[len(capex_updates.index)] = capex_data

                # Project out the scheduled payment series
                remaining_principal = initial_principal
                for j in range(debt_term):
                    interest_payment = rate * remaining_principal
                    principal_payment = total_payment - interest_payment

                    # Add the financial schedule data for this instrument
                    inst_sched_data = [
                        next_id,
                        agent_id,
                        pd_issued + j,
                        total_payment,
                        interest_payment,
                        principal_payment,
                    ]

                    fin_sched_updates.loc[len(fin_sched_updates.index)] = inst_sched_data

                    remaining_principal -= principal_payment

                next_id += 1

        # Once all new financial instruments have been set up, add the complete
        #   dataframe to the DB table
        fin_insts_updates.to_sql(
            "financial_instrument_manifest",
            self.db,
            if_exists = "append",
            index = False,
        )

        fin_sched_updates.to_sql(
            "financing_schedule",
            self.db,
            if_exists = "append",
            index=False,
        )

        capex_updates.to_sql(
            "capex_projections",
            self.db,
            if_exists = "append",
            index = False,
        )

        dep_updates.to_sql(
            "depreciation_projections",
            self.db,
            if_exists = "append",
            index = False,
        )

        self.db.commit()


    def execute_all_status_updates(self):
        # Record newly-started WIP projects from the agents' decisions
        WIP_updates = pd.read_sql_query("SELECT * FROM WIP_updates", self.db)
        WIP_updates.to_sql(
            "WIP_projects", self.db, if_exists="append", index=False
        )
        self.db.commit()

        # Set up all financial instruments corresponding to the new projects
        self.convert_WIP_updates_to_FIs(WIP_updates)

        # Record newly-started C2N WIP projects from the agents' decisions
        C2N_updates = pd.read_sql_query(
            "SELECT * FROM WIP_C2N_updates", self.db
        )
        C2N_updates.to_sql("WIP_C2N", self.db, if_exists="append", index=False)

        # Record status updates to existing assets (i.e. retirements)
        # Convert asset_updates to a dict of dicts for convenience
        asset_updates = pd.read_sql_query(
            "SELECT * FROM asset_updates", self.db
        )
        for row_num in asset_updates.index:
            new_record = (
                asset_updates.loc[[row_num]].copy().reset_index(drop=True)
            )

            orig_record = pd.read_sql_query(
                (
                    f"SELECT * FROM assets "
                    + f"WHERE asset_id = {new_record.loc[0, 'asset_id']}"
                ),
                self.db,
            )

            if len(orig_record) == 0:
                # The asset does not already exist and an entry must be added
                # Ensure that completion and retirement periods are integers
                new_record.at[0, "completion_pd"] = int(
                    math.ceil(new_record.at[0, "completion_pd"])
                )
                new_record.at[0, "retirement_pd"] = int(
                    math.ceil(new_record.at[0, "retirement_pd"])
                )
                pd.DataFrame(new_record).to_sql(
                    "assets", self.db, if_exists="append", index=False
                )
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
                    self.db, self.cur, "assets", new_record, filters
                )
        self.db.commit()

        # Delete all contents of the WIP_updates and asset_updates tables
        self.db.cursor().execute("DELETE FROM WIP_updates")
        self.db.cursor().execute("DELETE FROM asset_updates")
        self.db.commit()

