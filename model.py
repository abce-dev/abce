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

import subprocess
import yaml
import numpy as np
import pandas as pd
from mesa import Agent, Model
from mesa.time import RandomActivation

# import local modules
from agent import GenCo
import ABCEfunctions as ABCE
import seed_creator as sc
import price_curve as pc

class GridModel(Model):
    ''' A model with some number of GenCos. '''
    def __init__(self, settings, args):
        # Get input file locations from the settings dictionary
        unit_specs_file = settings["unit_specs_file"]
        fuel_data_file = settings["fuel_data_file"]
        demand_data_file = settings["demand_data_file"]
        price_curve_data_file = settings["price_curve_data_file"]
        time_series_data_file = settings["time_series_data_file"]
        db_file = settings["db_file"]
        # Get parameters from the settings dictionary
        self.num_agents = settings["num_agents"]
        self.first_agent_id = settings["first_agent_id"]
        self.first_asset_id = settings["first_asset_id"]
        self.total_forecast_horizon = settings["total_forecast_horizon"]

        # Copy the command-line arguments as member data
        self.args = args

        # Determine setting for use of a precomputed price curve
        self.use_precomputed_price_curve = True
        if "use_precomputed_price_curve" in settings:
            self.use_precomputed_price_curve = settings["use_precomputed_price_curve"]

        # Initialize the model one time step before the true start date
        self.current_step = -1

        # Initialize database for managing asset and WIP construction project data
        self.db_file = db_file
        self.db, self.cur = sc.create_database(self.db_file, self.args.force)

        # Load unit type specifications and fuel costs
        self.unit_specs = pd.read_csv(unit_specs_file)
        self.fuel_costs = pd.read_csv(fuel_data_file)
        self.add_units_to_db()

        # Load all-period demand data into the database
        self.load_demand_data_to_db(settings)

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Create agents
        for i in range(self.first_agent_id, self.first_agent_id + self.num_agents):
            gc = GenCo(i, self, settings, self.args)
            self.schedule.add(gc)

        # Check whether a market price subsidy is in effect, and its value
        self.set_market_subsidy(settings)
        # Create an appropriate price duration curve
        self.create_price_duration_curve(settings)

        # Save price duration data to the database
        self.price_duration_data.to_sql("price_curve",
                                        con = self.db,
                                        if_exists = "replace")


    def add_units_to_db(self):
        for i in range(len(self.unit_specs)):
            # Get relevant unit specs from file
            unit_type = self.unit_specs.loc[i, "unit_type"]
            fuel_type = self.unit_specs.loc[i, "fuel_type"]
            capacity = self.unit_specs.loc[i, "capacity"]
            uc_x = self.unit_specs.loc[i, "uc_x"]
            d_x = self.unit_specs.loc[i, "d_x"]
            heat_rate = self.unit_specs.loc[i, "heat_rate"]
            VOM = self.unit_specs.loc[i, "VOM"]
            FOM = self.unit_specs.loc[i, "FOM"]
            unit_life = self.unit_specs.loc[i, "unit_life"]
            CF = self.unit_specs.loc[i, "CF"]
            # Incorporate unit fuel cost data from the fuel costs file
            fuel_type_mask = self.fuel_costs.fuel_type == fuel_type
            fuel_cost = (self.fuel_costs[fuel_type_mask]
                         .reset_index()
                         .loc[0, "cost_per_mmbtu"])

            # Insert the values into the unit_specs DB table
            self.cur.execute(f"""INSERT INTO unit_specs VALUES
                            ('{unit_type}', '{fuel_type}', {capacity}, {uc_x},
                             {d_x}, {heat_rate}, {VOM}, {FOM}, {unit_life},
                             {CF}, {fuel_cost})""")


    def load_demand_data_to_db(self, settings):
        # Load all-period demand data into the database
        demand_df = pd.read_csv(settings["demand_data_file"]) * settings["peak_demand"]
        # Create an expanded range of periods to backfill with demand_df data
        new_index = list(range(self.total_forecast_horizon))
        demand_df = demand_df.reindex(new_index, method="ffill")
        # Save data to DB
        demand_df.to_sql("demand", self.db, if_exists="replace", index_label="period")


    def set_market_subsidy(self, settings):
        """
        If a subsidy is enabled, set the model's `subsidy_amount` to the
        quantity specified by the user in 'settings.yml'.

        Args:
          settings (dict): model settings from 'settings.yml'
        """
        self.subsidy_enabled = False
        self.subsidy_amount = 0

        if "enable_subsidy" in settings:
            self.subsidy_enabled = settings["enable_subsidy"]
            if self.subsidy_enabled and "subsidy_amount" in settings:
                self.subsidy_amount = settings["subsidy_amount"]

    def create_price_duration_curve(self, settings):
        # Set up the price curve according to specifications in settings
        if self.use_precomputed_price_curve:
            self.price_duration_data = pc.load_time_series_data(
                                             settings["price_curve_data_file"],
                                             file_type="price",
                                             subsidy=self.subsidy_amount,
                                             output_type = "dataframe")
        else:
            # Create the systemwide merit order curve
            self.merit_curve = pc.create_merit_curve(self.db, self.current_step)
            pc.plot_curve(self.merit_curve, plot_name="merit_curve.png")
            # Load five-minute demand data from file
            self.demand_data = pc.load_time_series_data(
                                     settings["time_series_data_file"],
                                     file_type="load",
                                     peak_demand=settings["peak_demand"])
            pc.plot_curve(self.demand_data, plot_name="demand_curve.png")
            # Create the final price duration curve
            self.price_duration_data = pc.compute_price_duration_curve(
                                              self.demand_data,
                                              self.merit_curve,
                                              settings["price_cap"])
            self.price_duration_data = pd.DataFrame({"lamda": self.price_duration_data})
            # Save a plot of the price duration curve
            pc.plot_curve(self.price_duration_data, plot_name="price_duration.png")


    def step(self):
        """
        Advance the model by one step.
        """
        self.current_step += 1
        if not self.args.quiet:
            print("\n\n\n")
        print("\n=========================================================================")
        print(f"Model step: {self.current_step}")
        print("==========================================================================")
        self.schedule.step()
        if not self.args.quiet:
            print("\nAll agent turns are complete.\n")
        # Reveal new information to all market participants
        projects_to_reveal = self.get_projects_to_reveal()
        self.reveal_decisions(projects_to_reveal)
        if not self.args.quiet:
            print("Table of all assets:")
            print(pd.read_sql("SELECT * FROM assets", self.db))
            print("Table of construction project updates:")
            print(pd.read_sql("SELECT * FROM WIP_projects", self.db).tail(n=8))

        print("Running A-LEAF...")
        aleaf_cmd = "julia /home/kbiegel/kb_aleaf/run.jl"
        if self.args.quiet:
            sp = subprocess.check_call([aleaf_cmd],
                                       shell=True,
                                       stdout=open(os.devnull, "wb"))
        else:
            sp = subprocess.check_call([aleaf_cmd], shell=True)


    def get_projects_to_reveal(self):
        """
        Set which (if any) unrevealed project construction decisions should
           be revealed (`revealed = True`).
        
        This function is a placeholder, which sets all unrevealed decisions
           to be revealed at the end of each decision round, after all
           agents have taken their turns. This behavior may be updated in
           the future.

        Returns:
           projects_to_reveal (pd DataFrame): one-column dataframe of ints,
             listing asset IDs to reveal
        """
        projects_to_reveal = pd.read_sql("SELECT asset_id FROM assets WHERE " +
                                         "revealed = 'false'", self.db)
        return projects_to_reveal

    def reveal_decisions(self, projects_to_reveal):
        """
        Set the indicated set of assets to have the `revealed = True` status
        in the database 'assets' table.

        Args:
           projects_to_reveal (pd DataFrame): one-column dataframe of ints,
             listing asset IDs to reveal
        """
        for asset_id in projects_to_reveal["asset_id"]:
            self.cur.execute(f"UPDATE assets SET revealed = 'true' WHERE " +
                             f"asset_id = '{asset_id}'")
        self.db.commit()

