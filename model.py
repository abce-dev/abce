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

import os
import shutil
import subprocess
import yaml
import numpy as np
import pandas as pd
import logging
from mesa import Agent, Model
from mesa.time import RandomActivation

# import local modules
from agent import GenCo
import ABCEfunctions as ABCE
import seed_creator as sc
import price_curve as pc
import ALEAF_interface as ALI

class GridModel(Model):
    ''' A model with some number of GenCos. '''
    def __init__(self, settings, args):
        self.settings = settings
        # Get input file locations from the settings dictionary
        unit_specs_file = settings["unit_specs_file"]
        fuel_data_file = settings["fuel_data_file"]
        demand_data_file = settings["demand_data_file"]
        price_curve_data_file = settings["seed_dispatch_data_file"]
        time_series_data_file = settings["time_series_data_file"]
        db_file = settings["db_file"]
        self.port_file_1a = settings["port_file_1a"]
        # Get agent parameters from the settings dictionary
        self.num_agents = settings["num_agents"]
        self.first_agent_id = settings["first_agent_id"]
        self.first_asset_id = settings["first_asset_id"]
        self.total_forecast_horizon = settings["total_forecast_horizon"]
        # Get ALEAF parameters from the settings dictionary
        self.ALEAF_abs_path = settings["ALEAF_abs_path"]
        self.ALEAF_master_settings_file = settings["ALEAF_master_settings_file"]
        self.ALEAF_model_type = settings["ALEAF_model_type"]
        self.ALEAF_region = settings["ALEAF_region"]
        self.ALEAF_model_settings_file_name = settings["ALEAF_model_settings_file"]
        self.ALEAF_portfolio_file = settings["ALEAF_portfolio_file"]
        self.ALEAF_scenario_name = settings["ALEAF_scenario_name"]
        # Get model/system parameters from the settings dictionary
        self.planning_reserve_margin = settings["planning_reserve_margin"]

        # Copy the command-line arguments as member data
        self.args = args

        # Initialize the model one time step before the true start date
        self.current_step = -1

        # Initialize database for managing asset and WIP construction project data
        self.db_file = db_file
        self.db, self.cur = sc.create_database(self.db_file, self.args.force)

        # Load all-period demand data into the database
        self.load_demand_data_to_db(settings)

        # Add model parameters to the database
        self.load_model_parameters_to_db(settings)

        # Set up all ALEAF file paths
        self.set_ALEAF_file_paths()

        # If ALEAF is enabled, re-initialize all input data based on local
        #    reference copies
        self.reinitialize_ALEAF_input_data()

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Create agents
        for i in range(self.first_agent_id, self.first_agent_id + self.num_agents):
            gc = GenCo(i, self, settings, self.args)
            self.schedule.add(gc)


        # Determine setting for use of a precomputed price curve
        self.use_precomputed_price_curve = True
        if "use_precomputed_price_curve" in settings:
            self.use_precomputed_price_curve = settings["use_precomputed_price_curve"]

        # Set the source for price data
        self.price_curve_data_file = settings["seed_dispatch_data_file"]

        # Load unit type specifications and fuel costs
        #self.unit_specs = pd.read_csv(unit_specs_file)
        #self.fuel_costs = pd.read_csv(fuel_data_file)
        self.add_units_to_db(from_ALEAF=True)

        # Reset the A-LEAF system portfolio by overwriting "ALEAF_ERCOT.xlsx"
        #    with the copy of "ALEAF_ERCOT_original.xlsx" which is stored
        #    in abce/inputs/ALEAF_inputs
        self.ALEAF_portfolio_new_path = os.path.join(self.ALEAF_abs_path,
                                                     "data",
                                                     self.ALEAF_model_type,
                                                     self.ALEAF_region,
                                                     f"ALEAF_{self.ALEAF_region}.xlsx")
        shutil.copyfile(f"./inputs/ALEAF_inputs/{self.port_file_1a}", self.ALEAF_portfolio_new_path)

        # Reset the A-LEAF model settings file by overwriting "ALEAF_Master_{model_type}.xlsx"
        #    with the copy stored in abce/inputs/aleaf_inputs
        self.ALEAF_model_settings_original_path = f"./inputs/ALEAF_inputs/ALEAF_Master_{self.ALEAF_model_type}_original.xlsx"
        self.ALEAF_model_settings_new_path = os.path.join(self.ALEAF_abs_path,
                                                          "setting",
                                                          f"ALEAF_Master_{self.ALEAF_model_type}.xlsx")
        shutil.copyfile(self.ALEAF_model_settings_original_path, self.ALEAF_model_settings_new_path)

        # Check whether a market price subsidy is in effect, and its value
        self.set_market_subsidy(settings)
                                               
        # Create an appropriate price duration curve
        self.create_price_duration_curve(settings)

        # Save price duration data to the database
        self.price_duration_data.to_sql("price_curve",
                                        con = self.db,
                                        if_exists = "replace")

        # Check ./outputs/ dir and clear out old files
        self.ALEAF_output_path = os.path.join(self.ALEAF_abs_path,
                                              "output",
                                              self.ALEAF_model_type,
                                              self.ALEAF_region,
                                              f"scenario_1_{self.ALEAF_scenario_name}")
        self.ABCE_output_path = os.path.join(os.getcwd(), "outputs", self.ALEAF_scenario_name)
        if not os.path.isdir(self.ABCE_output_path):
            # If the desired output directory doesn't already exist, create it
            os.makedirs(self.ABCE_output_path, exist_ok=True)
        else:
            # Otherwise, delete any existing files in the directory
            for existing_file in os.listdir(self.ABCE_output_path):
                os.remove(os.path.join(self.ABCE_output_path, existing_file))

            


    def set_ALEAF_file_paths(self):
        """ Set up all absolute paths to ALEAF and its input files, and
              save them as member data.
        """
        # Set file paths of local reference copies of ALEAF input data
        ALEAF_inputs_path = "./inputs/ALEAF_inputs"
        self.ALEAF_master_settings_ref = os.path.join(ALEAF_inputs_path,
                                                      "ALEAF_Master_original.xlsx")
        self.ALEAF_model_settings_ref = os.path.join(ALEAF_inputs_path,
                                                     f"ALEAF_Master_{self.ALEAF_model_type}_original.xlsx")
        self.ALEAF_portfolio_ref = os.path.join(ALEAF_inputs_path, 
                                                self.port_file_1a)

        # Set the paths to where settings are stored in the ALEAF directory
        ALEAF_settings_path = os.path.join(self.ALEAF_abs_path, "setting")
        self.ALEAF_master_settings_remote = os.path.join(ALEAF_settings_path,
                                                         "ALEAF_Master.xlsx")
        self.ALEAF_model_settings_remote = os.path.join(ALEAF_settings_path,
                                                        f"ALEAF_Master_{self.ALEAF_model_type}.xlsx")
        self.ALEAF_portfolio_remote = os.path.join(self.ALEAF_abs_path,
                                                   "data",
                                                   f"ALEAF_{self.ALEAF_region}.xlsx")
        self.ATB_remote = os.path.join(self.ALEAF_abs_path,
                                       "data",
                                       self.ALEAF_model_type,
                                       self.ALEAF_region,
                                       "ATBe.csv")

        # Set path to ALEAF outputs
        self.ALEAF_output_data_path = os.path.join(self.ALEAF_abs_path,
                                                  "output",
                                                  self.ALEAF_model_type,
                                                  self.ALEAF_region,
                                                  f"scenario_1_{self.ALEAF_scenario_name}",
                                                  f"{self.ALEAF_scenario_name}__dispatch_summary_OP.csv")


    def reinitialize_ALEAF_input_data(self):
        """ Setting ALEAF inputs requires overwriting fixed xlsx input files.
              This function re-initializes the following ALEAF input files from
              reference copies stored in ./inputs/ALEAF_inputs/.
        """
        # Update the ALEAF_Master_LC_GEP.xlsx model settings file:
        #  - Update the peak demand value in the 'Simulation Configuration' tab
        ALI.update_ALEAF_demand(self.ALEAF_model_settings_ref,
                                        self.ALEAF_model_settings_remote,
                                        self.db,
                                        period=0)

        # Update the ALEAF_ERCOT.xlsx system portfolio data:
        ALI.update_ALEAF_system_portfolio(self.ALEAF_portfolio_ref,
                                          self.ALEAF_portfolio_remote,
                                          self.db,
                                          self.current_step)

    def add_units_to_db(self, from_ALEAF=False):
        if from_ALEAF:
            # Set up header converter from A-LEAF to ABCE unit_spec format
            ALEAF_header_converter = {"UNIT_CATEGORY": "unit_type",
                                      "FUEL": "fuel_type",
                                      "CAP": "capacity",
                                      "INVC": "uc_x",
                                      "HR": "heat_rate",
                                      "CAPCRED": "CF"}

            # Set up the header converter from ATBe to ABCE unit_spec format
            ATB_header_converter = {"CAPEX": "uc_x",
                                    "Variable O&M": "VOM",
                                    "Fixed O&M": "FOM",
                                    "Fuel": "FC_per_MMBTU"}

            # Set up the unit name converter from ALEAF to ABCE
            ALEAF_unit_name_converter = {"Solar PV": "PV",
                                         "Gas CC": "NGCC",
                                         "Gas CT": "NGCT"}

            # Load the unit specs sheet from the settings file
            us_df = pd.read_excel(self.ALEAF_model_settings_ref, engine="openpyxl", sheet_name="Gen Technology")
            # Rename columns and make unit_type the row index
            us_df = us_df.rename(mapper=ALEAF_header_converter, axis=1)
            us_df = us_df.set_index("unit_type")
            # Select only the needed columns by getting the list of headers
            #   required for the unit_specs DB table
            self.cur.execute("SELECT * FROM unit_specs")
            columns_to_select = [item[0] for item in self.cur.description]
            columns_to_select.remove("unit_type")
            # Add columns of zeros for columns which will be computed later
            for column in columns_to_select:
                if column not in us_df.columns:
                    us_df[column] = 0
            # Create the final DataFrame for the unit specs data
            unit_specs_data = us_df[columns_to_select].copy()

            # Load the ATB search settings sheet from ALEAF
            ATB_settings = pd.read_excel(self.ALEAF_model_settings_ref, engine="openpyxl", sheet_name="ATB Setting")
            ATB_settings = ATB_settings.set_index("UNIT_CATEGORY")
            # Remove duplicates (related to multiple scenarios specified in the
            #   Simulation Configuration tab, which are unneeded here)
            ATB_settings = ATB_settings.loc[ATB_settings["ATB_Setting_ID"] == "ATB_ID_1", :]

            # Load the ATB database sheet
            ATB_data = pd.read_csv(self.ATB_remote)

            # Known names of columns which are filled with "ATB" in the
            #   ALEAF unit_spec sheet, names in the ATB format, plus
            #   capacity factor (CF) which is not given in the ALEAF sheet
            to_fill = ["CAPEX", "Variable O&M", "Fixed O&M", "Fuel"]

            for unit_type in list(us_df.index):
                # Retrieve the ATB search/matching settings for this unit type
                unit_settings = dict(ATB_settings.loc[unit_type, :])

                for datum in to_fill:
                    if datum in ATB_header_converter.keys():
                        datum_name = ATB_header_converter[datum]
                    else:
                        datum_name = datum
                    # TODO: set up an A-LEAF -> ATBe search term converter
                    # Attempt to match all search terms
                    mask = ((ATB_data["technology"] == unit_settings["Tech"]) &
                            (ATB_data["techdetail"] == unit_settings["TechDetail"]) &
                            (ATB_data["core_metric_parameter"] == datum) &
                            (ATB_data["core_metric_case"] == unit_settings["Case"]) &
                            (ATB_data["crpyears"] == unit_settings["CRP"]) &
                            (ATB_data["scenario"] == unit_settings["Scenario"]) &
                            (ATB_data["core_metric_variable"] == unit_settings["Year"]))

                    if sum(mask) != 1:
                        # If the mask matches nothing in ATBe, assume that
                        #   the appropriate value is 0 (e.g. battery VOM cost)
                        logging.debug(f"No match (or multiple matches) found for unit type {unit_type}; setting unit_specs value for {datum} to 0.")
                        unit_specs_data.loc[unit_type, datum_name] = 0
                    else:
                        unit_specs_data.loc[unit_type, datum_name] = ATB_data.loc[mask, "value"].values[0]

                # Retrieve the CRP data to fill the 'unit_life' data field
                unit_specs_data.loc[unit_type, "unit_life"] = unit_settings["CRP"]

            # Turn 'unit_type' back into a column from the index of unit_specs_data
            unit_specs_data = unit_specs_data.reset_index()
            # Compute fuel cost per kWh; conversion factor of 1e6 is for BTU -> MMBTU
            unit_specs_data["FC"] = unit_specs_data["FC_per_MMBTU"] * unit_specs_data["heat_rate"] / 1000000
            # FAKE: give all units a construction duration of 5 years
            unit_specs_data["d_x"] = 5
            # Cast the VOM column as Float64 (fixing specific bug)
            unit_specs_data["VOM"] = unit_specs_data["VOM"].astype("float64")

            # Change unit_specs_data unit_type entries to match ABCE standard
            for key in ALEAF_unit_name_converter.keys():
                mask = unit_specs_data["unit_type"] == key
                unit_specs_data.loc[mask, "unit_type"] = ALEAF_unit_name_converter[key]
            unit_specs_data.to_sql("unit_specs", self.db, if_exists = "replace", index = False)
            self.unit_specs = unit_specs_data

        else:
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


    def load_model_parameters_to_db(self, settings):
        # Load specific parameters specified in settings.yml to the database
        prm = settings["planning_reserve_margin"]
        self.cur.execute(f"INSERT INTO model_params VALUES ('PRM', {prm})")


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

    def create_price_duration_curve(self, settings, dispatch_data=None):
        # Set up the price curve according to specifications in settings
        if self.use_precomputed_price_curve:
            if self.current_step <= 0:
                price_curve_data_file = self.price_curve_data_file
            else:
                price_curve_data_file = dispatch_data
            self.price_duration_data = pc.load_time_series_data(
                                             price_curve_data_file,
                                             file_type="price",
                                             subsidy=self.subsidy_amount,
                                             output_type = "dataframe")
        else:
            # Create the systemwide merit order curve
            self.merit_curve = pc.create_merit_curve(self.db, self.current_step)
            pc.plot_curve(self.merit_curve, plot_name="merit_curve.png")
            # Load demand data from file
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

        # Update price data from ALEAF
        if self.current_step != 0:
            new_dispatch_data_filename = f"{self.ALEAF_scenario_name}__dispatch_summary_OP__step_{self.current_step - 1}.csv"
            new_dispatch_data = os.path.join(self.ABCE_output_path, new_dispatch_data_filename)
            print(f"Creating price duration curve using file {new_dispatch_data}")
            self.create_price_duration_curve(self.settings, new_dispatch_data)

            # Save price duration data to the database
            self.price_duration_data.to_sql("price_curve",
                                            con = self.db,
                                            if_exists = "replace")

        # Iterate through all agent turns
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

        if not self.args.no_aleaf:
            # Update the A-LEAF system portfolio based on any new units completed
            #    this round
            # If the current period is 0, then do not update the portfolio
            #    (already done in self.init())
            if self.current_step != 0:
                new_units = ALI.get_new_units(self.db, self.current_step)
                ALEAF_sys_portfolio_path = os.path.join(self.ALEAF_abs_path,
                                                        "data",
                                                        self.ALEAF_model_type,
                                                        self.ALEAF_region,
                                                        self.ALEAF_portfolio_file)
                ALI.update_ALEAF_system_portfolio(ALEAF_sys_portfolio_path, ALEAF_sys_portfolio_path, self.db, self.current_step)

            # Update ALEAF peak demand
            ALI.update_ALEAF_demand(self.ALEAF_model_settings_new_path,
                                    self.ALEAF_model_settings_new_path,
                                    self.db,
                                    self.current_step)

            # Run A-LEAF
            print("Running A-LEAF...")
            run_script_path = os.path.join(self.ALEAF_abs_path, "run.jl")
            ALEAF_sysimage_path = os.path.join(self.ALEAF_abs_path, "aleafSysimage.so")
            aleaf_cmd = f"julia -J{ALEAF_sysimage_path} {run_script_path} {self.ALEAF_abs_path}"
            if self.args.quiet:
                sp = subprocess.check_call([aleaf_cmd],
                                           shell=True,
                                           stdout=open(os.devnull, "wb"))
            else:
                sp = subprocess.check_call([aleaf_cmd], shell=True)

            # Copy all ALEAF output files to the output directory, with
            #   scenario and step-specific names
            files_to_save = ["dispatch_summary_OP", "expansion_result", "system_summary_OP", "system_tech_summary_OP"]
            for outfile in files_to_save:
                old_filename = f"{self.ALEAF_scenario_name}__{outfile}.csv"
                old_filepath = os.path.join(self.ALEAF_output_path, old_filename)
                new_filename = f"{self.ALEAF_scenario_name}__{outfile}__step_{self.current_step}.csv"
                new_filepath = os.path.join(self.ABCE_output_path, new_filename)
                shutil.copy2(old_filepath, new_filepath)

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

