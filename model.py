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


    def initialize_unit_specs_df(self):
        """
        Initialize the unit_specs dataframe, using data from the A-LEAF input
        file. Convert column headers to the ABCE standard names. Set up blank
        columns for values which will be added or computed later.

        Returns:
          unit_specs_data (DataFrame): dataframe containing the initial unit
            specs values loaded from the A-LEAF input. Still contains "ATB"
            entries for values to be searched, and several 0-columns.
        """
        # Set up header converter from A-LEAF to ABCE unit_spec format
        ALEAF_header_converter = {"UNIT_TYPE": "unit_type",
                                  "FUEL": "fuel_type",
                                  "CAP": "capacity",
                                  "CAPEX": "uc_x",
                                  "HR": "heat_rate",
                                  "FC": "FC_per_MWh",
                                  "CAPCRED": "CF",
                                  "VRE_Flag": "is_VRE"}

        # Load the unit specs sheet from the settings file
        us_df = pd.read_excel(self.ALEAF_model_settings_ref, engine="openpyxl", sheet_name="Gen Technology")

        # Rename columns from the A-LEAF standard to ABCE standard, and make
        #   "unit_type" the row index
        us_df = us_df.rename(mapper=ALEAF_header_converter, axis=1)
        us_df = us_df.set_index("unit_type")

        # Now that column names match the ABCE standard, select all columns
        #   from the df which have a matching column in the unit_specs DB table
        # This ensures that only columns which will be used are selected.
        # Load the list of column headers from the DB into the cursor object
        self.cur.execute("SELECT * FROM unit_specs")
        # Convert this result into a list
        columns_to_select = [item[0] for item in self.cur.description]
        # "unit_type" is not needed as it is the index
        columns_to_select.remove("unit_type")

        # Some columns are not pulled from the A-LEAF input data, and will
        #   be computed or added later on. Initialize these columns as zeroes.
        for column in columns_to_select:
            if column not in us_df.columns:
                us_df[column] = 0

        # Create the final DataFrame for the unit specs data
        unit_specs_data = us_df[columns_to_select].copy()

        return unit_specs_data


    def fill_unit_data_from_ATB(self, unit_specs_data):
        """
            In the A-LEAF unit specification, some data values may be
            initialized as "ATB", indicating that appropriate values should be
            retrieved from the ATB data sheet.

            This function matches the relevant unit data type with ATB search
            terms from the A-LEAF "ATB Setting" tab, filters the ATB data, and
            overwrites "ATB" entries with the appropriate data.

            Arguments:
              unit_specs_data (DataFrame)

            Returns:
              unit_specs_data (DataFrame): with all "ATB" values filled
        """
        # Set up the header converter from ATBe to ABCE unit_spec format
        ATB_header_read_converter = {"CAPEX": "uc_x",
                                     "Variable O&M": "VOM",
                                     "Fixed O&M": "FOM",
                                     "Fuel": "FC_per_MWh"}
        ATB_header_write_converter = {"CAPEX": "uc_x",
                                      "Variable O&M": "VOM",
                                      "Fixed O&M": "FOM",
                                      "Fuel": "ATB_FC"}

        # Set up the converter between A-LEAF input sheet search terms and
        #   ATB column headers
        ATB_search_terms_map = {"technology": "Tech",
                                "techdetail": "TechDetail",
                                "core_metric_case": "Case",
                                "crpyears": "CRP",
                                "scenario": "Scenario",
                                "core_metric_variable": "Year"}

        # Load the ATB search settings sheet from ALEAF
        ATB_settings = pd.read_excel(self.ALEAF_model_settings_ref, engine="openpyxl", sheet_name="ATB Setting")
        ATB_settings = ATB_settings.set_index("UNIT_TYPE")

        # Remove extraneous scenarios from the loaded A-LEAF ATB settings.
        #   Some versions of the A-LEAF input sheet may include multiple
        #   scenarios; all but the first will be ignored. ABCE does not allow
        #   multiple scenarios, so this simply prevents issues arising from
        #   copy-and-pasting from example A-LEAF inputs.
        ATB_settings = ATB_settings.loc[ATB_settings["ATB_Setting_ID"] == "ATB_ID_1", :]

        # Load the ATB database sheet
        ATB_data = pd.read_csv(self.ATB_remote)

        # Fill values for each unit type
        for unit_type in list(unit_specs_data.index):
            # Retrieve the ATB search/matching settings for this unit type
            unit_settings = dict(ATB_settings.loc[unit_type, :])

            # Fill all necessary columns for this unit type
            for datum_name in ATB_header_read_converter.keys():
                if unit_specs_data.loc[unit_type, ATB_header_read_converter[datum_name]] == "ATB":
                    # Construct the mask using the ATB search terms map
                    #   defined earlier
                    mask = (ATB_data["core_metric_parameter"] == datum_name)
                    for ATB_key, ALEAF_key in ATB_search_terms_map.items():
                        mask = mask & (ATB_data[ATB_key] == unit_settings[ALEAF_key])

                    if sum(mask) != 1:
                        # If the mask matches nothing in ATBe, assume that
                        #   the appropriate value is 0 (e.g. battery VOM cost),
                        #   but warn the user.
                        logging.warn(f"No match (or multiple matches) found for unit type {unit_type}; setting unit_specs value for {datum_name} to 0.")
                        unit_specs_data.loc[unit_type, ATB_header_write_converter[datum_name]] = 0
                    else:
                        unit_specs_data.loc[unit_type, ATB_header_write_converter[datum_name]] = ATB_data.loc[mask, "value"].values[0]
                        if datum_name == "Fuel":
                            # If the current datum is Fuel, also record its
                            #   associated units
                            unit_specs_data.loc[unit_type, "ATB_FC_units"] = ATB_data.loc[mask, "units"].values[0]

        # Set newly-filled ATB data columns to numeric data types
        #   Columns initialized with "ATB" will be a non-numeric data type,
        #   which causes problems later if not converted.
        for ATB_key, ABCE_key in ATB_header_write_converter.items():
            unit_specs_data[ABCE_key] = pd.to_numeric(unit_specs_data[ABCE_key])

        return unit_specs_data


    def finalize_unit_specs_data(self, unit_specs_data):
        """
        Fill in supplemental unit specification data from
          ABCE files, ensure all fuel cost data is in $/MWh,
          and finalize the layout of the dataframe.

        Arguments:
          unit_specs_data (DataFrame)

        Returns:
          unit_specs_data (DataFrame)
        """
        # Turn 'unit_type' back into a column from the index of unit_specs_data
        unit_specs_data = unit_specs_data.reset_index()

        # Retrieve non-ALEAF parameters from the ABCE supplemental unit
        #   specification file
        unit_specs_ABCE = pd.read_csv(os.path.join(self.settings["ABCE_abs_path"],
                                                   self.settings["unit_specs_abce_supp_file"]))

        # Some generators' (currently NG and Coal) fuel cost is given in the
        #   ATB data in units of $/MMBTU. Convert these values to a $/MWh basis
        #   for consistency.
        for i in range(len(unit_specs_data)):
            if unit_specs_data.loc[i, "ATB_FC_units"] == "$/MWh":
                unit_specs_data.loc[i, "FC_per_MWh"] = unit_specs_data.loc[i, "ATB_FC"]
            elif unit_specs_data.loc[i, "ATB_FC_units"] == "$/MMBTU":
                unit_specs_data.loc[i, "FC_per_MWh"] = unit_specs_data.loc[i, "ATB_FC"] * unit_specs_data.loc[i, "heat_rate"]
            else:
                if not unit_specs_data.loc[i, "is_VRE"]:
                    # If units are VRE, FC units are irrelevant
                    # Otherwise, throw an error
                    raise ValueError(
                        f"The unit {unit_specs_data.loc[i, 'unit_type']}
                        has its fuel cost specified in units of
                        {unit_specs_data.loc[i, 'ATB_FC_units']}. I'm not sure
                        how to convert these units. Double-check your inputs or
                        edit model.py to handle this case."
                    )

        # Set unit baseline construction duration and life from supplemental data
        for i in range(len(unit_specs_data)):
            unit_type = unit_specs_data.loc[i, "unit_type"]
            # Set construction duration for this unit
            unit_specs_data.loc[i, "d_x"] = unit_specs_ABCE[unit_specs_ABCE["unit_type"] == unit_type]["d_x"].values[0]
            # Set unit useful life for this unit
            unit_specs_data.loc[i, "unit_life"] = unit_specs_ABCE[unit_specs_ABCE["unit_type"] == unit_type]["unit_life"].values[0]

        return unit_specs_data




    def add_unit_specs_to_db(self):
        """
        This function loads all unit specification data and saves it to the
          database, as well as to the member data `self.unit_specs`.
        The following steps are performed:
          - Initial data is loaded from the A-LEAF input file.
          - The initial data is updated to follow ABCE format conventions.
          - Data values initialized as "ATB" are searched for and filled from
              the NREL Annual Technology Baseline (2020) (ATB) file.
          - The final unit specification data is loaded from ABCE supplemental
              files.
          - Final formatting steps.
          - Data is saved to the database table `unit_specs` and to
              self.unit_specs.
        """

        # Read data from A-LEAF input file; convert column headers to ABCE
        #   standard; set up placeholder columns for values to be computed
        #   later
        unit_specs_data = self.initialize_unit_specs_df()

        # Match any data values initialized as "ATB" to appropriate entries in
        #   the ATB data sheet, and overwrite them.
        unit_specs_data = self.fill_unit_data_from_ATB(unit_specs_data)

        # Finalize the unit_specs_data dataframe, including ABCE supplemental
        #   data and layout change
        unit_specs_data = self.finalize_unit_specs_data(unit_specs_data)

        # Save the finalized unit specs data to the DB, and set the member data
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

