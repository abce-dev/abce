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
    def __init__(self, settings_file_name, settings, args):
        self.settings_file_name = settings_file_name
        self.settings = settings
        # Get agent parameters from the settings dictionary
        self.num_agents = settings["num_agents"]
        self.first_agent_id = settings["first_agent_id"]
        self.first_asset_id = settings["first_asset_id"]
        self.total_forecast_horizon = settings["total_forecast_horizon"]
        # Get ALEAF parameters from the settings dictionary
        self.ALEAF_abs_path = settings["ALEAF_abs_path"]
        self.ALEAF_master_settings_file_name = settings["ALEAF_master_settings_file"]
        self.ALEAF_model_type = settings["ALEAF_model_type"]
        self.ALEAF_region = settings["ALEAF_region"]
        self.ALEAF_model_settings_file_name = settings["ALEAF_model_settings_file"]
        self.ALEAF_portfolio_file_name = settings["ALEAF_portfolio_file"]
        self.ALEAF_scenario_name = settings["ALEAF_scenario_name"]
        # Get model/system parameters from the settings dictionary
        self.planning_reserve_margin = settings["planning_reserve_margin"]
        # Unit conversions
        self.MW2kW = 1000          # Converts MW to kW
        self.MMBTU2BTU = 1000000   # Converts MMBTU to BTU

        # Copy the command-line arguments as member data
        self.args = args

        # Initialize the model one time step before the true start date
        self.current_step = -1

        # Initialize database for managing asset and WIP construction project data
        self.db_file = os.path.join(settings["ABCE_abs_path"],
                                    settings["db_file"])
        self.db, self.cur = sc.create_database(self.db_file, self.args.force)

        # Load all-period demand data into the database
        self.load_demand_data_to_db(self.settings)

        # Add model parameters to the database
        self.load_model_parameters_to_db(self.settings)

        # Set up all ALEAF file paths
        self.set_ALEAF_file_paths(settings)

        # If ALEAF is enabled, re-initialize all input data based on local
        #    reference copies
        self.reinitialize_ALEAF_input_data()

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Add unit type specifications to database, including all parameters
        #   to be loaded from the ATB database and the ABCE supplemental
        #   data file
        self.add_unit_specs_to_db()

        # Create agents
        for i in range(self.first_agent_id, self.first_agent_id + self.num_agents):
            gc = GenCo(i, self, settings, self.args)
            self.schedule.add(gc)

        # Determine setting for use of a precomputed price curve
        self.use_precomputed_price_curve = True
        if "use_precomputed_price_curve" in settings:
            self.use_precomputed_price_curve = settings["use_precomputed_price_curve"]

        # Check whether a market price subsidy is in effect, and its value
        self.set_market_subsidy(self.settings)
                                               
        # Create an appropriate price duration curve
        self.create_price_duration_curve(self.settings)

        # Save price duration data to the database
        self.price_duration_data.to_sql("price_curve",
                                        con = self.db,
                                        if_exists = "replace")

        # Check ./outputs/ dir and clear out old files
        # TODO: replace getcwd() with a command-line argument to specify a 
        #   non-cwd ABCE absolute path
        self.ABCE_output_data_path = os.path.join(settings["ABCE_abs_path"], "outputs", self.ALEAF_scenario_name)
        if not os.path.isdir(self.ABCE_output_data_path):
            # If the desired output directory doesn't already exist, create it
            os.makedirs(self.ABCE_output_data_path, exist_ok=True)
        else:
            # Otherwise, delete any existing files in the directory
            for existing_file in os.listdir(self.ABCE_output_data_path):
                os.remove(os.path.join(self.ABCE_output_data_path, existing_file))

            


    def set_ALEAF_file_paths(self, settings):
        """ Set up all absolute paths to ALEAF and its input files, and
              save them as member data.
        """
        # Set file paths of local reference copies of ALEAF input data
        ALEAF_inputs_path = os.path.join(settings["ABCE_abs_path"], "inputs/ALEAF_inputs")
        self.ALEAF_master_settings_ref = os.path.join(ALEAF_inputs_path,
                                                      "ALEAF_Master_original.xlsx")
        self.ALEAF_model_settings_ref = os.path.join(ALEAF_inputs_path,
                                                     f"ALEAF_Master_{self.ALEAF_model_type}_original.xlsx")
        self.ALEAF_portfolio_ref = os.path.join(ALEAF_inputs_path, 
                                                self.settings["port_file_1a"])

        # Set the paths to where settings are stored in the ALEAF directory
        ALEAF_settings_path = os.path.join(self.ALEAF_abs_path, "setting")
        self.ALEAF_master_settings_remote = os.path.join(ALEAF_settings_path,
                                                         self.ALEAF_master_settings_file_name)
        self.ALEAF_model_settings_remote = os.path.join(ALEAF_settings_path,
                                                        f"ALEAF_Master_{self.ALEAF_model_type}.xlsx")
        self.ALEAF_portfolio_remote = os.path.join(self.ALEAF_abs_path,
                                                   "data",
                                                   self.ALEAF_model_type,
                                                   self.ALEAF_region,
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
                                                  f"scenario_1_{self.ALEAF_scenario_name}")


    def reinitialize_ALEAF_input_data(self):
        """ Setting ALEAF inputs requires overwriting fixed xlsx input files.
              In order to avoid cross-run contamination, these inputs should
              be reset to the user-provided baseline (the settings files in
              ./inputs/ALEAF_inputs/) before the simulation starts.
            This function re-initializes the following ALEAF input files from
              reference copies stored in ./inputs/ALEAF_inputs/.
        """
        # Update the ALEAF_Master_LC_GEP.xlsx model settings file:
        #  - Update the peak demand value in the 'Simulation Configuration' tab
        ALI.update_ALEAF_model_settings(self.ALEAF_model_settings_ref,
                                        self.ALEAF_model_settings_remote,
                                        self.db,
                                        self.settings,
                                        period=0)

        # Update the ALEAF_ERCOT.xlsx system portfolio data:
        ALI.update_ALEAF_system_portfolio(self.ALEAF_portfolio_ref,
                                          self.ALEAF_portfolio_remote,
                                          self.db,
                                          self.current_step)


    def add_unit_specs_to_db(self):
        """
        Load in the A-LEAF unit specification data, including looking up data
          from the NREL Annual Technology Baseline (2020) file as specified by
          the user in ALEAF_Master_LC_GEP.xlsx.
        Data is saved to the database table `unit_specs`.
        """
        # Set up header converter from A-LEAF to ABCE unit_spec format
        ALEAF_header_converter = {"UNIT_TYPE": "unit_type",
                                  "FUEL": "fuel_type",
                                  "CAP": "capacity",
                                  "CAPEX": "uc_x",
                                  "HR": "heat_rate",
                                  "CAPCRED": "CF"}

        # Set up the header converter from ATBe to ABCE unit_spec format
        ATB_header_converter = {"CAPEX": "uc_x",
                                "Variable O&M": "VOM",
                                "Fixed O&M": "FOM",
                                "Fuel": "FC_per_MWh"}

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
        us_df["FC_per_MWh"] = "ATB"
        # Create the final DataFrame for the unit specs data
        unit_specs_data = us_df[columns_to_select].copy()

        # Load the ATB search settings sheet from ALEAF
        ATB_settings = pd.read_excel(self.ALEAF_model_settings_ref, engine="openpyxl", sheet_name="ATB Setting")
        ATB_settings = ATB_settings.set_index("UNIT_TYPE")
        # Remove duplicates (related to multiple scenarios specified in the
        #   Simulation Configuration tab, which are unneeded here)
        ATB_settings = ATB_settings.loc[ATB_settings["ATB_Setting_ID"] == "ATB_ID_1", :]

        # Load the ATB database sheet
        ATB_data = pd.read_csv(self.ATB_remote)

        for unit_type in list(unit_specs_data.index):
            # Retrieve the ATB search/matching settings for this unit type
            unit_settings = dict(ATB_settings.loc[unit_type, :])

            for datum_name in ATB_header_converter.keys():
                if unit_specs_data.loc[unit_type, ATB_header_converter[datum_name]] == "ATB":
                    mask = ((ATB_data["technology"] == unit_settings["Tech"]) &
                            (ATB_data["techdetail"] == unit_settings["TechDetail"]) &
                            (ATB_data["core_metric_parameter"] == datum_name) &
                            (ATB_data["core_metric_case"] == unit_settings["Case"]) &
                            (ATB_data["crpyears"] == unit_settings["CRP"]) &
                            (ATB_data["scenario"] == unit_settings["Scenario"]) &
                            (ATB_data["core_metric_variable"] == unit_settings["Year"]))

                    if sum(mask) != 1:
                        # If the mask matches nothing in ATBe, assume that
                        #   the appropriate value is 0 (e.g. battery VOM cost)
                        logging.warn(f"No match (or multiple matches) found for unit type {unit_type}; setting unit_specs value for {datum_name} to 0.")
                        unit_specs_data.loc[unit_type, ATB_header_converter[datum_name]] = 0
                    else:
                        unit_specs_data.loc[unit_type, ATB_header_converter[datum_name]] = ATB_data.loc[mask, "value"].values[0]

            # Retrieve the units' is_VRE status
            unit_specs_data.loc[unit_type, "is_VRE"] = us_df[us_df.index == unit_type]["VRE_Flag"].values[0]

        # Convert the uc_x column to numeric
        unit_specs_data["uc_x"] = pd.to_numeric(unit_specs_data["uc_x"])
        unit_specs_data["FC_per_MWh"] = pd.to_numeric(unit_specs_data["FC_per_MWh"])

        # Turn 'unit_type' back into a column from the index of unit_specs_data
        unit_specs_data = unit_specs_data.reset_index()

        # Retrieve non-ALEAF parameters from the ABCE supplemental unit
        #   specification file
        unit_specs_ABCE = pd.read_csv(os.path.join(self.settings["ABCE_abs_path"],
                                                   self.settings["unit_specs_abce_supp_file"]))

        # Set unit baseline construction duration and life from supplemental data
        for i in range(len(unit_specs_data)):
            unit_type = unit_specs_data.loc[i, "unit_type"]
            # Set construction duration for this unit
            unit_specs_data.loc[i, "d_x"] = unit_specs_ABCE[unit_specs_ABCE["unit_type"] == unit_type]["d_x"].values[0]
            # Set unit useful life for this unit
            unit_specs_data.loc[i, "unit_life"] = unit_specs_ABCE[unit_specs_ABCE["unit_type"] == unit_type]["unit_life"].values[0]

        # Cast the VOM and FOM columns as Float64 (fixing specific bug)
        unit_specs_data["VOM"] = unit_specs_data["VOM"].astype("float64")
        unit_specs_data["FOM"] = unit_specs_data["FOM"].astype("float64")

        unit_specs_data.to_sql("unit_specs", self.db, if_exists = "replace", index = False)
        self.unit_specs = unit_specs_data


    def load_demand_data_to_db(self, settings):
        # Load all-period demand data into the database
        demand_data_file = os.path.join(settings["ABCE_abs_path"],
                                        settings["demand_data_file"])
        demand_df = pd.read_csv(demand_data_file) * settings["peak_demand"]
        # Create an expanded range of periods to backfill with demand_df data
        new_index = list(range(self.total_forecast_horizon))
        demand_df = demand_df.reindex(new_index, method="ffill")
        # Save data to DB
        demand_df.to_sql("demand", self.db, if_exists="replace", index_label="period")


    def load_model_parameters_to_db(self, settings):
        # Load specific parameters specified in settings.yml to the database
        prm = settings["planning_reserve_margin"]
        self.cur.execute(f"INSERT INTO model_params VALUES ('PRM', {prm})")

        tax_rate = settings["tax_rate"]
        self.cur.execute(f"INSERT INTO model_params VALUES ('tax_rate', {tax_rate})")


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
                price_curve_data_file = os.path.join(settings["ABCE_abs_path"],
                                                     settings["seed_dispatch_data_file"])
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
            time_series_data_file = os.path.join(settings["ABCE_abs_path"],
                                                 settings["time_series_data_file"])
            self.demand_data = pc.load_time_series_data(
                                     time_series_data_file,
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
            new_dispatch_data = os.path.join(self.ABCE_output_data_path, new_dispatch_data_filename)
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

        # Update the A-LEAF system portfolio based on any new units completed
        #    this round
        # If the current period is 0, then do not update the portfolio
        #    (already done in self.init())
        if self.current_step != 0:
            ALI.update_ALEAF_system_portfolio(self.ALEAF_portfolio_remote, self.ALEAF_portfolio_remote, self.db, self.current_step)

        # Update ALEAF peak demand
        ALI.update_ALEAF_model_settings(self.ALEAF_model_settings_remote,
                                        self.ALEAF_model_settings_remote,
                                        self.db,
                                        self.settings,
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

        self.save_ALEAF_outputs()

    def save_ALEAF_outputs(self):
        # Copy all ALEAF output files to the output directory, with
        #   scenario and step-specific names
        files_to_save = ["dispatch_summary_OP", "expansion_result", "system_summary_OP", "system_tech_summary_OP"]
        for outfile in files_to_save:
            old_filename = f"{self.ALEAF_scenario_name}__{outfile}.csv"
            old_filepath = os.path.join(self.ALEAF_output_data_path, old_filename)
            new_filename = f"{self.ALEAF_scenario_name}__{outfile}__step_{self.current_step}.csv"
            new_filepath = os.path.join(self.ABCE_output_data_path, new_filename)
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

