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

        # Load the unit-type ownership specification for all agents
        self.portfolio_specification = pd.read_csv(os.path.join(self.settings["ABCE_abs_path"], self.settings["portfolios_file"]))
        # Check the portfolio specification to ensure the ownership totals 
        #   equal the total numbers of available units
        self.check_portfolio_specification()

        # Load the mandatory unit retirement data
        ret_data_file = os.path.join(self.settings["ABCE_abs_path"],
                                     settings["retirement_period_specs_file"])
        self.ret_data = pd.read_csv(ret_data_file, comment="#")

        # Create agents
        for agent_id in range(self.first_agent_id, self.first_agent_id + self.num_agents):
            gc = GenCo(agent_id, self, settings, self.args)
            self.schedule.add(gc)
            self.initialize_agent_assets(agent_id)

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
            for datum_name, ALEAF_read_col in ATB_header_read_converter.items():
                if unit_specs_data.loc[unit_type, ALEAF_read_col] == "ATB":
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
        unit_problems = dict()
        unit_specs_data["FC_per_MWh"] = unit_specs_data.apply(
            lambda x:
                x["ATB_FC"] if x["ATB_FC_units"] == "$/MWh" else
                    (x["ATB_FC"] * x["heat_rate"] if x["ATB_FC_units"] == "$/MMBTU" else
                        (0 if x["is_VRE"] == True else
                            (unit_problems.update({x["unit_type"]: x["ATB_FC_units"]})))),
            axis = 1
        )

        if len(unit_problems) != 0:
            raise ValueError(f"I don't recognize the fuel cost units provided in the following cases: {unit_problems}. Check your inputs, or update model.py to handle these additional cases.")

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


    def check_portfolio_specification(self):
        """
        Ensure that the user-supplied portfolio ownership specification is
        correctly set up.
          - Check number of agents versus self.settings["num_agents"]
          - Ensure that the number of extant units of each type equals the
            total number of units of that type owned by all agents

        If either of this criteria is not completely satisfied, give the user
        an error message and end the program.
        """
        # Tracker for portfolio state
        agents_OK = True
        portfolio_OK = True
        num_agents_msg = "Number of agents matches."
        units_owned_msg = "Unit type ownership matches."

        # Read in the A-LEAF total system portfolio
        book, writer = ALI.prepare_xlsx_data(self.ALEAF_portfolio_ref, self.ALEAF_portfolio_ref)
        pdf = ALI.organize_ALEAF_portfolio(writer)

        # Check step 1: ensure number of agents matches
        if self.settings["num_agents"] != self.portfolio_specification["agent_id"].nunique():
            agents_OK = False
            num_agents_msg = f"\nNumber of agents specified differs between settings and portfolio specification.\nSettings file num_agents: {self.settings['num_agents']}  |  Portfolio spec num_agents: {self.portfolio_specification['agent_id'].nunique()}."


        # Check step 2: ensure total assets owned by agent by type equals
        #   total assets by type in the system
        owned_units = self.portfolio_specification.groupby("unit_type")["num_units"].sum()
        incorrect_units = dict()
        for unit_type in list(pdf["Unit Type"]):
            if pdf[pdf["Unit Type"] == unit_type]["EXUNITS"].values[0] != owned_units[unit_type]:
                portfolio_OK = False
                incorrect_units[unit_type] = (pdf[pdf["Unit Type"] == unit_type]["EXUNITS"].values[0], owned_units[unit_type])

        # If unit type numbers don't match, construct a helpful error message for the user
        if not portfolio_OK:
            msg_list = []
            for unit_type in incorrect_units.keys():
                msg = f"{unit_type}:: System Portfolio Total: {incorrect_units[unit_type][0]}  |  Ownership Specification Total: {incorrect_units[unit_type][1]}"
                msg_list.append(msg)
            preamble = "\nThe number of units by generator type does not match between the overall system portfolio specification and the individual agents' portfolios:\n"
            postamble = "\n\nEnsure that the total number of units for each type matches between the system portfolio file and the agent ownership file. \nTerminating..."
            units_owned_msg = preamble + "\n".join(msg_list) + postamble

        # If something went wrong in the checks, display an error message and quit
        if not (agents_OK and portfolio_OK):
            print("\nError: the agent portfolio specification does not match other simulation settings.")
            if not agents_OK:
                print(num_agents_msg)
            if not portfolio_OK:
                print(units_owned_msg)
            exit()


    def initialize_agent_assets(self, agent_id):
        # Get the agent-specific portfolio by unit type: filter the manifest
        #   of agent unit ownership for the current agent's unique ID
        pdf = self.portfolio_specification[self.portfolio_specification["agent_id"] == agent_id]

        # Set the initial asset ID
        asset_id = ABCE.get_next_asset_id(self.db, self.settings["first_asset_id"])

        # Set up the master asset dataframe to hold all asset data for saving
        #   into the database
        master_asset_df = pd.DataFrame()

        # Assign units to this agent as specified in the portfolio file,
        #   and record each in the master_asset_df dataframe
        for unit_type in list(pdf["unit_type"]):
            # Retrieve the list of retirement period data for this unit type
            #   and agent, and create the cumulative-sum threshold mapping
            unit_rets = self.create_unit_type_retirement_df(unit_type, asset_id)

            for j in range(pdf.loc[pdf["unit_type"] == unit_type, "num_units"].values[0]):
                # Assign the appropriate retirement period
                retirement_pd = self.assign_retirement_pd(unit_rets, asset_id)

                # Compute unit capex according to the unit type spec
                unit_capex = self.compute_total_capex_preexisting(unit_type)

                # Default: assume all pre-existing assets have $0 outstanding
                #   financing balance. To add obligations for capital payments
                #   on pre-existing assets, uncomment the code lines below
                # Compute the asset's annual capital payment, if the asset is
                #   not paid off. The asset's unit_life value is used as the
                #   repayment term for financing (i.e. useful life = financing
                #   maturity period).
                #unit_life = self.unit_specs.loc[unit_type, "unit_life"]
                #unit_cap_pmt = self.compute_sinking_fund_payment(agent_id, unit_capex, unit_life)
                cap_pmt = 0

                # Create the dictionary of asset characteristics
                asset_dict = {"asset_id": asset_id,
                              "agent_id": agent_id,
                              "unit_type": unit_type,
                              "start_pd": -1,
                              "completion_pd": 0,
                              "cancellation_pd": 9999,
                              "retirement_pd": retirement_pd,
                              "total_capex": unit_capex,
                              "cap_pmt": cap_pmt}
                new_asset = pd.DataFrame(asset_dict, index=[0])

                # use the first asset's DataFrame as the template
                if len(master_asset_df) == 0:
                    master_asset_df = new_asset
                else:
                    master_asset_df = master_asset_df.append(new_asset)

                # Increment the asset id for the next asset
                asset_id = max(ABCE.get_next_asset_id(self.db, asset_id), max(master_asset_df["asset_id"])+1)

        # Save the dataframe of all assets into the 'assets' DB table
        master_asset_df.to_sql("assets", self.db, if_exists="append", index=False)


    def compute_total_capex_preexisting(self, unit_type):
        unit_cost_per_kW = self.unit_specs.loc[unit_type, "uc_x"]
        unit_capacity = self.unit_specs.loc[unit_type, "capacity"]

        total_capex = unit_cost_per_kW * unit_capacity * self.MW2kW

        return total_capex


    def create_unit_type_retirement_df(self, unit_type, starting_asset_id):
        """
        Create the step-function mapping to determine which units are assigned
          which mandatory retirement periods.
        """
        # Filter the df of retirement period data for the current unit type
        unit_type_rets = self.ret_data[self.ret_data["unit_type"] == unit_type].copy()

        # Sort from soonest to furthest-away retirement period
        unit_type_rets = unit_type_rets.sort_values(by="retirement_pd", axis=0, ascending=True, ignore_index=True)

        # Generate the thresholds for each retirement period, starting with
        #   the next available asset id. These thresholds indicate the END
        #   (i.e. NON-INCLUSIVE) of each assignment interval. A row with an
        #   `rp_threshold` of 2005 and a `retirement_pd` of 12 means that any
        #   assets with IDs strictly less than 2005 will receive a 
        #   `retirement_pd` of 12 (unless they qualify for a lower threshold
        #   category).
        unit_type_rets["rp_threshold"] = np.cumsum(unit_type_rets["num_copies"].to_numpy(), axis=0) + starting_asset_id

        return unit_type_rets


    def assign_retirement_pd(self, unit_type_rets, current_asset_id):
        """Compare the current asset ID against the classification thresholds
           in unit_type_rets, to determine the asset's mandatory retirement
           period.
        """
        # If no retirement period is specified, set the asset to never retire
        #    (i.e. retirement period is a Large Number)
        retirement_pd = 9999

        for i in range(len(unit_type_rets)):
            if current_asset_id < unit_type_rets.loc[i, "rp_threshold"]:
                retirement_pd = unit_type_rets.loc[i, "retirement_pd"]
                break   # Once the correct bin is found, exit the loop

        return retirement_pd


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

        # Transfer all decisions and updates from the 'asset_updates' and
        #   'WIP_updates' tables into their respective public-information
        #   equivalents
        self.execute_all_status_updates()

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
        WIP_projects = pd.read_sql("SELECT * FROM assets WHERE " +
                                   f"completion_pd >= {self.current_step} AND " +
                                   f"cancellation_pd > {self.current_step}",
                                   self.db)

        # Update each project one at a time
        for asset_id in WIP_projects.asset_id:
            # Select this project's most recent data record
            project_data = self.retrieve_project_data(asset_id)

            # If construction cost escalation is enabled in the settings file,
            #   update project_data to reflect cost and schedule escalation
            #   which occurred during the previous period
            if self.settings["xtr_escalation_enabled"]:
                project_data = self.escalate_cost_and_schedule(project_data)

            # Update total capital expenditures to date and current RCEC/RTEC
            project_data = self.advance_project_to_current_period(project_data)

            # If this period's authorized expenditures (ANPE) clear the RCEC,
            #   then the project is complete
            if project_data.loc[0, "rcec"] <= self.settings["large_epsilon"]:
                # Record the project's completion period as the current period
                self.record_completed_xtr_project(project_data)

            # Record updates to the WIP project's status
            self.record_WIP_project_updates(project_data)

            # Record updates to the project's expected completion date in the
            #   database 'assets' table
            self.update_expected_completion_period(project_data)


    def retrieve_project_data(self, asset_id):
        project_data = pd.read_sql("SELECT * FROM WIP_projects WHERE " +
                                   f"asset_id = {asset_id} AND " +
                                   f"period = {self.current_step - 1}",
                                   self.db)

        temp_asset_data = pd.read_sql("SELECT * FROM assets WHERE " +
                                      f"asset_id = {asset_id}", self.db)

        project_data = project_data.merge(temp_asset_data,
                                          how = "inner",
                                          on = ["asset_id", "agent_id"])

        return project_data


    def record_completed_xtr_project(self, project_data):
        asset_id = project_data.loc[0, "asset_id"]

        # Record the project's completion period as the current period
        self.cur.execute(f"UPDATE assets " +
                         f"SET completion_pd = {self.current_step} " +
                         f"WHERE asset_id = {asset_id}")

        # Compute and record the total CapEx for the project
        self.cur.execute(f"UPDATE assets " +
                         f"SET total_capex = {project_data.cum_exp.values[0]} " +
                         f"WHERE asset_id = {asset_id}")

        # Compute periodic sinking fund payments
        unit_type = project_data.loc[0, "unit_type"]
        unit_life = self.unit_specs.loc[unit_type, "unit_life"]
        capex_payment = self.compute_sinking_fund_payment(project_data.loc[0, "agent_id"], project_data.loc[0, "cum_exp"], unit_life)
        self.cur.execute(f"UPDATE assets SET cap_pmt = {capex_payment} " +
                         f"WHERE asset_id = {asset_id}")

        # Commit changes to database
        self.db.commit()


    def record_WIP_project_updates(self, project_data):
        # Set new_anpe = 0 to avoid inter-period contamination
        new_anpe = 0
        # Update the 'WIP_projects' table with new RCEC/RTEC/ANPE values
        self.cur.execute("INSERT INTO WIP_projects VALUES " +
                         f"({project_data.asset_id.values[0]}, " +
                         f"{project_data.agent_id.values[0]}, " +
                         f"{self.current_step}, " +
                         f"{project_data.cum_occ.values[0]}, " +
                         f"{project_data.rcec.values[0]}, " +
                         f"{project_data.cum_d_x.values[0]}, " +
                         f"{project_data.rtec.values[0]}, " +
                         f"{project_data.cum_exp.values[0]}, " +
                         f"{new_anpe})")

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

        wacc = (agent_params.debt_fraction * agent_params.cost_of_debt
                + (1 - agent_params.debt_fraction) * agent_params.cost_of_equity)
        cap_pmt = total_capex * wacc / (1 - (1 + wacc)**(-term))

        # Convert cap_pmt from a single-entry Series to a scalar
        cap_pmt = cap_pmt[0]

        return cap_pmt


    def escalate_cost_and_schedule(self, project_data):
        """
        Generate escalation for a unit's remaining cost expected to completion
        (RCEC) and remaining time expected to completion (RTEC).
        """
        # Retrieve unit type data for easy access
        unit_type = project_data.loc[0, "unit_type"]
        unit_type_data = self.unit_specs.loc[unit_type, :]

        # Compute time delay
        escalated_rtec = project_data.loc[0, "rtec"]

        # Compute cost overrun
        cost_esc_factor = random.uniform(-1.0, 1.0) * unit_type_data["occ_variance"]
        new_cum_occ = project_data.loc[0, "cum_occ"] * (1 + cost_esc_factor)

        # Compute forced schedule delay due to cost increase beyond maximum
        #   productive investment rate
        # Temporary assumption: the max production rate is equivalent to
        #   building the nominal-cost unit in the nominal expected construction
        #   time
        max_invest_rate = unit_type_data["uc_x"] * unit_type_data["capacity"] * self.MW2kW / unit_type_data["d_x"]
        # RTEC is always an integer, due to discrete simulation time periods
        forced_rtec = (new_cum_occ - project_data.loc[0, "cum_exp"]) / max_invest_rate
        # Overall RTEC is the greater of the two generated RTEC values
        schedule_slip = (escalated_rtec - project_data.loc[0, "rtec"]) + max(forced_rtec - project_data.loc[0, "rtec"], 0)
        project_data.loc[0, "cum_d_x"] += schedule_slip
        project_data.loc[0, "rtec"] = math.ceil(project_data.loc[0, "cum_d_x"] - (self.current_step - project_data.loc[0, "start_pd"] - 1))

        # Cover as much cost escalation as possible within the period's ANPE
        #   without exceeding the maximum productive investment rate
        project_data.loc[0, "anpe"] = min(max_invest_rate, (project_data.loc[0, "anpe"] + (new_cum_occ - project_data.loc[0, "cum_occ"])))

        # Compute final new RCEC and cumulative overnight capital cost
        project_data.loc[0, "rcec"] += (new_cum_occ - project_data.loc[0, "cum_occ"])
        project_data.loc[0, "cum_occ"] = new_cum_occ

        return project_data

    def advance_project_to_current_period(self, project_data):
        # Escalated RTEC is reduced by one
        project_data.loc[0, "rtec"] -= 1

        # Escalated RCEC is reduced by ANPE
        project_data.loc[0, "rcec"] -= project_data.loc[0, "anpe"]

        # Cumulative capital expenditures are increased by ANPE
        project_data.loc[0, "cum_exp"] += project_data.loc[0, "anpe"]

        # Project ANPE is reset to 0 after being expended
        project_data.loc[0, "anpe"] = 0

        return project_data


    def update_expected_completion_period(self, project_data):
        asset_id = project_data.loc[0, "asset_id"]
        new_completion_pd = project_data.loc[0, "rtec"] + self.current_step
        
        self.cur.execute(f"UPDATE assets SET completion_pd = {new_completion_pd} " +
                         f"WHERE asset_id = {asset_id}")
        self.db.commit()


    def execute_all_status_updates(self):
        WIP_updates = pd.read_sql_query("SELECT * FROM WIP_updates", self.db)
        asset_updates = pd.read_sql_query("SELECT * FROM asset_updates", self.db)

        # Record newly-started WIP projects from the agents' decisions
        WIP_updates.to_sql("WIP_projects", self.db, if_exists="append", index=False)
        self.db.commit()

        # Record status updates to existing assets (i.e. retirements)
        for i in range(len(asset_updates)):
            new_data = asset_updates.iloc[[i]].copy().reset_index().drop("index", axis=1)
            orig_data = pd.read_sql_query(f"SELECT * FROM assets WHERE asset_id = {new_data.loc[0, 'asset_id']}", self.db)
            if len(orig_data) == 0:
                # The asset does not already exist and an entry must be added
                new_data.to_sql("assets", self.db, if_exists="append", index=False)
            else:
                # Update any columns for which the new_data values do not match
                #   the orig_data values
                for header in new_data.columns:
                    if new_data.loc[0, header] != orig_data.loc[0, header] and header == "retirement_pd":
                        self.cur.execute(f"UPDATE assets SET {header} = {new_data.loc[0, header]} WHERE asset_id = {new_data.loc[0, 'asset_id']}")
        self.db.commit()






