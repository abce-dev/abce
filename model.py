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
from agent import GenCo
import ABCEfunctions as ABCE
import seed_creator as sc
import price_curve as pc
import ALEAF_interface as ALI

import warnings
warnings.filterwarnings("ignore")


class GridModel(Model):
    ''' A model with some number of GenCos. '''

    def __init__(self, settings_file_name, settings, args):
        self.settings_file_name = settings_file_name
        self.settings = settings
        self.solver = settings['solver'].lower()

        # Get agent parameters from the settings dictionary
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
        self.policies = settings["policies"]
        # Unit conversions
        self.MW2kW = 1000          # Converts MW to kW
        self.MMBTU2BTU = 1000000   # Converts MMBTU to BTU

        if 'natural_gas_price' in settings:
            self.natgas_price = settings['natural_gas_price']
        if 'conv_nuclear_FOM' in settings:
            self.conv_nuclear_FOM = settings['conv_nuclear_FOM']

        if 'ATB_year' in settings:
            self.ATB_year = settings['ATB_year']
        else:
            self.ATB_year = 2020
        print(f"Using ATB Year {self.ATB_year}")

        # Copy the command-line arguments as member data
        self.args = args

        # Initialize the model one time step before the true start date
        self.current_pd = -1

        # Initialize database for managing asset and WIP construction project
        # data

        self.db_file = (Path.cwd() / settings["db_file"])
        self.db, self.cur = sc.create_database(self.db_file, self.args.force)

        # Load all-period demand data into the database
        self.load_demand_data_to_db(self.settings)

        # Add model parameters to the database
        self.load_model_parameters_to_db(self.settings)

        # Create the local tmp/ directory inside the current working directory,
        #   if it doesn't already exist
        tmp_dir_location = (Path.cwd() / "tmp")

        Path(tmp_dir_location).mkdir(exist_ok=True)

        # Set up all ALEAF file paths
        self.set_ALEAF_file_paths(settings)

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Add unit type specifications to database, including all parameters
        #   to be loaded from the ATB database and the ABCE supplemental
        #   data file
        self.add_unit_specs_to_db()

        if self.settings["run_ALEAF"]:
            # Initialize the ALEAF model settings and generation technologies
            self.reinitialize_ALEAF_input_data()
            # Initialize the correct policy adjustments by unit type
            ALI.update_ALEAF_policy_settings(
                self.ALEAF_model_settings_remote,
                self.ALEAF_model_settings_remote,
                self.settings["policies"],
                self.unit_specs)

        # Read in the GenCo parameters data from file
        gc_params_file_name = Path(self.settings["ABCE_abs_path"] /
                                   self.settings["gc_params_file"])
        self.gc_params = yaml.load(
            open(
                gc_params_file_name,
                'r'),
            Loader=yaml.FullLoader)

        # Load the unit-type ownership specification for all agents
        self.portfolio_specification = pd.read_csv(
            Path(self.settings["ABCE_abs_path"]) / self.settings["portfolios_file"])

        # Check the portfolio specification to ensure the ownership totals
        #   equal the total numbers of available units
        self.check_num_agents()
        # self.check_total_assets()

        # Load the mandatory unit retirement data
        ret_data_file = Path(self.settings["ABCE_abs_path"] /
                             settings["retirement_period_specs_file"])
        self.ret_data = pd.read_csv(ret_data_file, comment="#")

        # Create agents
        num_agents = len(self.gc_params.keys())
        for agent_id in list(self.gc_params.keys()):
            gc = GenCo(
                agent_id,
                self,
                settings,
                self.gc_params[agent_id],
                self.args)
            self.schedule.add(gc)
            self.initialize_agent_assets(agent_id)

        # Determine setting for use of a precomputed price curve
        self.use_precomputed_price_curve = settings["use_precomputed_price_curve"]

        self.db.commit()

        # Check ./outputs/ dir and clear out old files
        self.ABCE_output_data_path = Path(
            os.getcwd()) / "outputs" / self.ALEAF_scenario_name
        if not Path(self.ABCE_output_data_path).is_dir():
            # If the desired output directory doesn't already exist, create it
            Path(self.ABCE_output_data_path).mkdir(exist_ok=True, parents=True)
        else:
            # Otherwise, delete any existing files in the directory
            for existing_file in Path(self.ABCE_output_data_path).iterdir():
                print("trying to remove directory")
                (Path(self.ABCE_output_data_path) / existing_file).unlink()

    def set_ALEAF_file_paths(self, settings):
        """ Set up all absolute paths to ALEAF and its input files, and
              save them as member data.
        """

        self.ALEAF_remote_path = Path(self.ALEAF_abs_path)
        self.ALEAF_remote_data_path = (Path(self.ALEAF_abs_path) /
                                       "data" /
                                       self.ALEAF_model_type /
                                       self.ALEAF_region)
        # Set file paths of local reference copies of ALEAF input data
        ALEAF_inputs_path = Path(
            settings["ABCE_abs_path"]) / "inputs" / "ALEAF_inputs"
        self.ALEAF_master_settings_ref = (
            Path(ALEAF_inputs_path) /
            settings["ALEAF_master_settings_file"])
        self.ALEAF_model_settings_ref = (Path(ALEAF_inputs_path) /
                                         settings["ALEAF_model_settings_file"])
        self.ALEAF_portfolio_ref = (Path(ALEAF_inputs_path) /
                                    settings["ALEAF_portfolio_file"])

        # Set the paths to where settings are stored in the ALEAF directory
        ALEAF_settings_path = self.ALEAF_remote_path / "setting"
        self.ALEAF_master_settings_remote = (
            Path(ALEAF_settings_path) /
            self.ALEAF_master_settings_file_name)
        self.ALEAF_model_settings_remote = (
            Path(ALEAF_settings_path) /
            f"ALEAF_Master_{self.ALEAF_model_type}.xlsx")
        self.ALEAF_portfolio_remote = (self.ALEAF_remote_data_path /
                                                   f"ALEAF_{self.ALEAF_region}.xlsx")
        # self.ATB_remote = (self.ALEAF_remote_data_path /
        #                                     "ATBe.csv")  
        self.ATB_remote = (ALEAF_inputs_path /
                                            "ATBe.csv")                                                 
        # Set path to ALEAF outputs
        self.ALEAF_output_data_path = (self.ALEAF_remote_path/
                                                  "output"/
                                                  self.ALEAF_model_type/
                                                  self.ALEAF_region/
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
                                          self.current_pd)

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
                                  "FC": "original_FC",
                                  "CAPCRED": "CF",
                                  "VRE_Flag": "is_VRE",
                                  "EMSFAC": "emissions_rate"}

        # Load the unit specs sheet from the settings file
        us_df = pd.read_excel(
            self.ALEAF_model_settings_ref,
            engine="openpyxl",
            sheet_name="Gen Technology")

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
        self.db.commit()

        # Some columns are not pulled from the A-LEAF input data, and will
        #   be computed or added later on. Initialize these columns as zeroes.
        for column in columns_to_select:
            if column not in us_df.columns:
                us_df[column] = 0

        # The original EMSFAC column from A-LEAF is in strange units (100*tCO2/MWh);
        #   convert to tCO2/MWh
        us_df["emissions_rate"] = us_df["emissions_rate"] / 100

        try:
            ng_fuel = us_df['fuel_type'] == 'Gas'
            us_df.loc[ng_fuel, 'original_FC'] = self.natgas_price
            print(f'using specified value: {self.natgas_price}')
        except AttributeError:
            print('Using standard value.')
        try:
            nuke_fuel = us_df['UNITGROUP'] == 'ConventionalNuclear'
            us_df.loc[nuke_fuel, 'FOM'] = self.conv_nuclear_FOM
            print(f'using specified value: {self.conv_nuclear_FOM}')
        except AttributeError:
            print('Using standard value.')

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
                                     "Fuel": "original_FC"}
        ATB_header_write_converter = {"CAPEX": "uc_x",
                                      "Variable O&M": "VOM",
                                      "Fixed O&M": "FOM",
                                      "Fuel": "original_FC"}

        # Set up the converter between A-LEAF input sheet search terms and
        #   ATB column headers
        ATB_search_terms_map = {"technology": "Tech",
                                "techdetail": "TechDetail",
                                "core_metric_case": "Case",
                                "crpyears": "CRP",
                                "scenario": "Scenario",
                                "core_metric_variable": "Year"}

        # Load the ATB search settings sheet from ALEAF
        ATB_settings = pd.read_excel(
            self.ALEAF_model_settings_ref,
            engine="openpyxl",
            sheet_name="ATB Setting")
        ATB_settings = ATB_settings.set_index("UNIT_TYPE")

        # Remove extraneous scenarios from the loaded A-LEAF ATB settings.
        #   Some versions of the A-LEAF input sheet may include multiple
        #   scenarios; all but the first will be ignored. ABCE does not allow
        #   multiple scenarios, so this simply prevents issues arising from
        #   copy-and-pasting from example A-LEAF inputs.
        ATB_settings = ATB_settings.loc[ATB_settings["ATB_Setting_ID"]
                                        == "ATB_ID_1", :]

        # Load the ATB database sheet
        ATB_data = pd.read_csv(self.ATB_remote)
        # print(unit_specs_data.iloc[:, 1:17])

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
                        mask = mask & (
                            ATB_data[ATB_key] == unit_settings[ALEAF_key])

                    if sum(mask) != 1:
                        # If the mask matches nothing in ATBe, assume that
                        #   the appropriate value is 0 (e.g. battery VOM cost),
                        #   but warn the user.
                        logging.warn(
                            f"No match (or multiple matches) found for unit type {unit_type}; setting unit_specs value for {datum_name} to 0.")
                        unit_specs_data.loc[unit_type,
                                            ATB_header_write_converter[datum_name]] = 0
                    else:
                        unit_specs_data.loc[unit_type,
                                            ATB_header_write_converter[datum_name]] = ATB_data.loc[mask,
                                                                                                   "value"].values[0]
                        if datum_name == "Fuel":
                            # If the current datum is Fuel, also record its
                            #   associated units
                            unit_specs_data.loc[unit_type,
                                                "original_FC_units"] = ATB_data.loc[mask,
                                                                                    "units"].values[0]
                elif (ALEAF_read_col == "original_FC") and (unit_specs_data.loc[unit_type, ALEAF_read_col] != "ATB"):
                    unit_specs_data.loc[unit_type,
                                        "original_FC_units"] = "$/MWh"

        # Set newly-filled ATB data columns to numeric data types
        #   Columns initialized with "ATB" will be a non-numeric data type,
        #   which causes problems later if not converted.
        for ATB_key, ABCE_key in ATB_header_write_converter.items():
            unit_specs_data[ABCE_key] = pd.to_numeric(unit_specs_data[ABCE_key])

        return unit_specs_data

    def set_up_policies(self, unit_specs_data):
        """
        Allocate any policy impacts (carbon tax or PTC) to the various unit types

        Sign convention:
          POSITIVE values: subsidy
          NEGATIVE values: penalty
        The value stored in unit_specs_data[unit_type, "policy_adj_per_MWh"]
          is later added to the MARKET price (NOT the unit's bid price)
        """
        unit_specs_data = unit_specs_data.reset_index()

        valid_CTAX_names = ["CTAX", "ctax", "carbon_tax", "carbontax"]
        valid_PTC_names = [
            "PTC",
            "ptc",
            "production_tax_credit",
            "productiontaxcredit"]

        if self.policies is not None:
            for key, val in self.policies.items():
                if val["enabled"]:
                    if key in valid_CTAX_names:
                        unit_specs_data["policy_adj_per_MWh"] = unit_specs_data.apply(
                            lambda x: x["policy_adj_per_MWh"] - x["emissions_rate"] * val["qty"], axis=1)
                    elif key in valid_PTC_names:
                        unit_specs_data["policy_adj_per_MWh"] = unit_specs_data.apply(
                            lambda x: x["policy_adj_per_MWh"] +
                            val["qty"] if x["unit_type"] in val["eligible"] else x["policy_adj_per_MWh"],
                            axis=1)
                    else:
                        err_msg = f"Sorry: the system policy {key} is not implemented, or might be misspelled."
                        raise ValueError(err_msg)

        unit_specs_data = unit_specs_data.set_index("unit_type")

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
        unit_specs_ABCE = pd.read_csv(
            Path(
                self.settings["ABCE_abs_path"]) /
            self.settings["unit_specs_abce_supp_file"])

        # Some generators' (currently NG and Coal) fuel cost is given in the
        #   ATB data in units of $/MMBTU. Convert these values to a $/MWh basis
        #   for consistency.
        # Logic flow:
        #   1. if the unit already has a numeric entry in FC_per_MWh, keep it. else:
        #   2. if the unit's ATB_FC is in $/MWh, copy that value to FC_per_MWh. else:
        #   3. if the unit's ATB_FC is in $/MMBTU, multiply its ATB_FC by its heat rate and copy that value to FC_per_MWh. else:
        #   4. if the unit is flagged as is_VRE == True, use the value of 0. else:
        # 5. record the current unit as having an unresolvable unit problem, and
        # throw an error to the user.
        unit_problems = dict()
        unit_specs_data["FC_per_MWh"] = unit_specs_data.apply(
            lambda x:
                x["original_FC"] if x["original_FC_units"] == "$/MWh" else
            (x["original_FC"] * x["heat_rate"] if x["original_FC_units"] == "$/MMBTU" else
                    (0 if x["is_VRE"] else
                     (unit_problems.update({x["unit_type"]: x["original_FC_units"]})))),
            axis=1
        )

        if len(unit_problems) != 0:
            raise ValueError(
                f"I don't recognize the fuel cost units provided in the following cases: {unit_problems}. Check your inputs, or update model.py to handle these additional cases.")

        # Set unit baseline construction duration and life from supplemental
        # data
        for i in range(len(unit_specs_data)):
            unit_type = unit_specs_data.loc[i, "unit_type"]
            current_unit_ABCE_data = unit_specs_ABCE[unit_specs_ABCE["unit_type"] == unit_type]
            # Set construction duration for this unit
            unit_specs_data.loc[i,
                                "d_x"] = current_unit_ABCE_data["d_x"].values[0]
            # Set unit useful life for this unit
            unit_specs_data.loc[i,
                                "unit_life"] = current_unit_ABCE_data["unit_life"].values[0]
            # Set unit lead time until a coal unit must be retired, if any
            unit_specs_data.loc[i,
                                "cpp_ret_lead"] = current_unit_ABCE_data["cpp_ret_lead"].values[0]
            # Set number of mandatory coal unit retirements
            unit_specs_data.loc[i,
                                "num_cpp_rets"] = current_unit_ABCE_data["num_cpp_rets"].values[0]
            # Set unit revenue head start vs end of xtr period
            unit_specs_data.loc[i,
                                "rev_head_start"] = current_unit_ABCE_data["rev_head_start"].values[0]

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

        unit_specs_data = self.set_up_policies(unit_specs_data)

        # Finalize the unit_specs_data dataframe, including ABCE supplemental
        #   data and layout change
        unit_specs_data = self.finalize_unit_specs_data(unit_specs_data)

        # Save the finalized unit specs data to the DB, and set the member data
        unit_specs_data.to_sql(
            "unit_specs",
            self.db,
            if_exists="replace",
            index=False)
        self.unit_specs = unit_specs_data

    def check_num_agents(self):
        """
        Ensure that all sources of agent data include the same number of
          agents. If not, raise a ValueError.
        """
        # Ensure that any agents specified in the portfolio spec file have
        #   corresponding parameters entries in the gc_params file:
        #   - gc_params.yml => gc_params
        #   - portfolios.csv => self.portfolio_specification
        if not all(agent_id in self.gc_params.keys()
                   for agent_id in self.portfolio_specification.agent_id.unique()):
            num_agents_msg = f"Agents specified in the portfolio specification file {self.settings['portfolios_file']} do not all have corresponding entries in the gc_params file {self.settings['gc_params_file']}. Check your inputs and try again."
            raise ValueError(num_agents_msg)

    def check_total_assets(self):
        """
        Ensure that the total of all assets by type owned by the agents
          matches the master total in the A-LEAF system portfolio file.
          If not, raise a ValueError.
        """
        # Temporarily read in the starting A-LEAF total system portfolio from
        #   the database
        book, writer = ALI.prepare_xlsx_data(
            self.ALEAF_portfolio_ref,
            self.ALEAF_portfolio_ref
        )
        full_pdf = ALI.organize_ALEAF_portfolio(writer)
        system_portfolio = (full_pdf[["Unit Type", "EXUNITS"]]
                            .rename(columns={"Unit Type": "unit_type"})
                            .set_index("unit_type"))

        agent_owned_units = (self.portfolio_specification
                             .groupby("unit_type")["num_units"].sum())
        # Inner join the agent ownership data into the system portfolio data
        system_portfolio = (system_portfolio.join(
            agent_owned_units,
            on="unit_type",
            how="inner")
            .reset_index())

        # Set up a dictionary to log unit number mismatches
        incorrect_units = dict()
        # For any unit_types where the num_units (from agent ownership) do not
        #   match the EXUNITS (from system portfolio), save that type's info to
        #   the incorrect_units dict
        system_portfolio = system_portfolio.apply(
            lambda x:
                incorrect_units.update(
                    {x["unit_type"]: (x["EXUNITS"], x["num_units"])}
                ) if x["num_units"] != x["EXUNITS"] else None,
            axis=1
        )

        # If unit type numbers don't match, construct a helpful error message
        # for the user
        if len(incorrect_units) != 0:
            msg_list = []
            for unit_type in incorrect_units.keys():
                msg = f"{unit_type}:: System Portfolio Total: {incorrect_units[unit_type][0]}  |  Ownership Specification Total: {incorrect_units[unit_type][1]}"
                msg_list.append(msg)
            preamble = "\nThe number of units by generator type does not match between the overall system portfolio specification and the individual agents' portfolios:\n"
            postamble = "\n\nEnsure that the total number of units for each type matches between the system portfolio file and the agent ownership file. \nTerminating..."
            units_owned_msg = preamble + "\n".join(msg_list) + postamble
            raise ValueError(units_owned_msg)

    def initialize_agent_assets(self, agent_id):
        # Get the agent-specific portfolio by unit type: filter the manifest
        #   of agent unit ownership for the current agent's unique ID
        pdf = self.portfolio_specification[self.portfolio_specification["agent_id"] == agent_id]
        # Ensure that all units of a given type are collated into a single
        #   dataframe row
        pdf = pdf.drop("agent_id", axis=1)
        pdf = pdf.groupby(by=["unit_type"]).sum().reset_index()

        # Set the initial asset ID
        asset_id = ABCE.get_next_asset_id(
            self.db, self.settings["first_asset_id"])

        # Retrieve the column-header schema for the 'assets' table
        self.cur.execute("SELECT * FROM assets")
        assets_col_names = [element[0] for element in self.cur.description]
        self.db.commit()

        # Create a master dataframe to hold all asset records for this
        #   agent-unit type combination (to reduce the frequency of saving
        #   to the database)
        master_assets_df = pd.DataFrame(columns=assets_col_names)

        # Assign units to this agent as specified in the portfolio file,
        #   and record each in the master_asset_df dataframe
        for unit_record in pdf.itertuples():
            unit_type = unit_record.unit_type
            num_units = unit_record.num_units

            # Retrieve the list of retirement period data for this unit type
            #   and agent
            unit_rets = self.create_unit_type_retirement_df(
                unit_type, agent_id, asset_id)

            # Out of the total number of assets of this type belonging to this
            #   agent, any "left over" units with no specified retirement date
            #   should all retire at period 9999
            assets_remaining = num_units
            num_not_specified = num_units - sum(unit_rets["num_copies"])
            unit_rets = unit_rets.append({
                "agent_id": agent_id,
                "unit_type": unit_type,
                "retirement_pd": 9999,
                "num_copies": num_not_specified
            }, ignore_index=True)

            # Create assets in blocks, according to the number of units per
            #   retirement period
            for rets_row in unit_rets.itertuples():
                retirement_pd = rets_row.retirement_pd
                num_copies = rets_row.num_copies

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

                # Save all values which don't change for each asset in this
                #   block, i.e. everything but the asset_id
                asset_dict = {"agent_id": agent_id,
                              "unit_type": unit_type,
                              "start_pd": -1,
                              "completion_pd": 0,
                              "cancellation_pd": 9999,
                              "retirement_pd": retirement_pd,
                              "total_capex": unit_capex,
                              "cap_pmt": cap_pmt,
                              "C2N_reserved": 0
                              }

                # For each asset in this block, create a dataframe record and
                #   store it to the master_assets_df
                for i in range(num_copies):
                    # Find the largest extant asset id, and set the current
                    #   asset id 1 higher
                    asset_dict["asset_id"] = max(
                        ABCE.get_next_asset_id(
                            self.db, self.settings["first_asset_id"]), max(
                            master_assets_df["asset_id"], default=self.settings["first_asset_id"]) + 1)

                    # Convert the dictionary to a dataframe format and save
                    new_record = pd.DataFrame(asset_dict, index=[0])
                    master_assets_df = master_assets_df.append(new_record)

            # For any leftover units in assets_remaining with no specified
            #   retirement date, initialize them with the default retirement date
            #   of 9999
            asset_dict["retirement_pd"] = 9999

        # Once all assets from all unit types for this agent have had records
        #   initialized, save the dataframe of all assets into the 'assets'
        #   DB table
        master_assets_df.to_sql(
            "assets",
            self.db,
            if_exists="append",
            index=False)

        self.db.commit()

        financing_row = (agent_id,
                         -1,
                         "debt",
                         0,
                         self.gc_params[agent_id]["starting_debt"],
                         20,
                         self.gc_params[agent_id]["starting_debt"])
        self.cur.execute(
            "INSERT INTO financing_schedule VALUES (?, ?, ?, ?, ?, ?, ?)",
            financing_row)

        debt_row = (agent_id, -1, self.gc_params[agent_id]["starting_debt"])
        self.cur.execute("INSERT INTO agent_debt VALUES (?, ?, ?)", debt_row)
        self.db.commit()

    def compute_total_capex_preexisting(self, unit_type):
        unit_cost_per_kW = self.unit_specs[self.unit_specs.unit_type ==
                                           unit_type]["uc_x"].values[0]
        unit_capacity = self.unit_specs[self.unit_specs.unit_type ==
                                        unit_type]["capacity"].values[0]

        total_capex = unit_cost_per_kW * unit_capacity * self.MW2kW

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

    def load_demand_data_to_db(self, settings):
        # Load all-period demand data into the database
        demand_data_file = (Path(settings["ABCE_abs_path"]) /
                            settings["demand_data_file"])
        demand_df = pd.read_csv(demand_data_file) * settings["peak_demand"]
        # Create an expanded range of periods to backfill with demand_df data
        new_index = list(range(self.total_forecast_horizon))
        demand_df = demand_df.reindex(new_index, method="ffill")
        # Save data to DB
        demand_df.to_sql(
            "demand",
            self.db,
            if_exists="replace",
            index_label="period")
        self.db.commit()

    def load_model_parameters_to_db(self, settings):
        # Load specific parameters specified in settings.yml to the database
        prm = settings["planning_reserve_margin"]
        self.cur.execute(f"INSERT INTO model_params VALUES ('PRM', {prm})")

        tax_rate = settings["tax_rate"]
        self.cur.execute(
            f"INSERT INTO model_params VALUES ('tax_rate', {tax_rate})")
        self.db.commit()

    def create_price_duration_curve(self, settings, dispatch_data=None):
        # Set up the price curve according to specifications in settings
        if self.use_precomputed_price_curve:
            if (self.current_pd <= 0) or (self.settings["run_ALEAF"] == False):
                price_curve_data_file = (Path(settings["ABCE_abs_path"]) /
                                         settings["seed_dispatch_data_file"])
            else:
                price_curve_data_file = dispatch_data
            self.price_duration_data = pc.load_time_series_data(
                price_curve_data_file,
                self.current_pd,
                file_type="price",
                output_type="dataframe")
        else:
            # Create the systemwide merit order curve
            self.merit_curve = pc.create_merit_curve(self.db, self.current_pd)
            pc.plot_curve(self.merit_curve, plot_name="merit_curve.png")
            # Load demand data from file
            time_series_data_file = (Path(settings["ABCE_abs_path"]) /
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
            self.price_duration_data = pd.DataFrame(
                {"lamda": self.price_duration_data})
            # Save a plot of the price duration curve
            pc.plot_curve(
                self.price_duration_data,
                plot_name="price_duration.png")

    def step(self, demo=False):
        """
        Advance the model by one step.
        """
        self.current_pd += 1

        if self.current_pd == 0:
            self.has_ABCE_sysimage, self.has_dispatch_sysimage = self.check_for_sysimage_files()

        if not self.args.quiet:
            pass
            # print("\n\n\n")
        # print("\n=========================================================================")
        # print(f"Model step: {self.current_pd}")
        # print("==========================================================================")

        # Update price data from ALEAF
        if (self.current_pd == 0) or (self.settings["run_ALEAF"] == False):
            self.create_price_duration_curve(self.settings)
        else:
            new_dispatch_data_filename = f"{self.ALEAF_scenario_name}__dispatch_summary_OP__step_{self.current_pd - 1}.csv"
            new_dispatch_data = Path(
                self.ABCE_output_data_path) / new_dispatch_data_filename
            # print(f"Creating price duration curve using file {new_dispatch_data}")
            self.create_price_duration_curve(self.settings, new_dispatch_data)

        # Save price duration data to the database
        self.price_duration_data.to_sql("price_curve",
                                        con=self.db,
                                        index=False,
                                        if_exists="append")
        self.db.commit()

        # Advance the status of all WIP projects to the current period
        self.update_WIP_projects()

        # Update financial statements and financial projections for all agents
        self.update_agent_financials()

        # Compute the scenario reduction results for this year
        ABCE.execute_scenario_reduction(
            self.db,
            self.current_pd,
            self.settings,
            self.unit_specs,
            self.settings["num_repdays"])

        self.db.commit()
        self.db.close()

        # Iterate through all agent turns
        self.schedule.step()
        if not self.args.quiet:
            print("\nAll agent turns are complete.\n")

        self.db = sqlite3.connect(str(Path.cwd() / self.settings["db_file"]))
        self.cur = self.db.cursor()

        # Transfer all decisions and updates from the 'asset_updates' and
        #   'WIP_updates' tables into their respective public-information
        #   equivalents
        self.execute_all_status_updates()

        if demo:
            # print("\n")
            user_response = input("Press Enter to continue: ")

        if not self.args.quiet:
            print("Table of all assets:")
            # print(pd.read_sql("SELECT * FROM assets", self.db))
            # print("Table of construction project updates:")
            # print(pd.read_sql("SELECT * FROM WIP_projects", self.db).tail(n=8))

        if self.settings["run_ALEAF"]:
            # Update the A-LEAF system portfolio based on any new units completed
            #   or units retired this period
            ALI.update_ALEAF_system_portfolio(
                self.ALEAF_portfolio_remote,
                self.ALEAF_portfolio_remote,
                self.db,
                self.current_pd)

            # Update ALEAF peak demand
            ALI.update_ALEAF_model_settings(self.ALEAF_model_settings_remote,
                                            self.ALEAF_model_settings_remote,
                                            self.db,
                                            self.settings,
                                            self.current_pd)

            # Run A-LEAF
            # print("Running A-LEAF...")
            run_script_path = self.ALEAF_remote_path / "execute_ALEAF.jl"
            ALEAF_env_path = self.ALEAF_remote_path / "."
            ALEAF_sysimage_path = self.ALEAF_remote_path / "aleafSysimage.so"
            aleaf_cmd = f"julia --project={ALEAF_env_path} -J {ALEAF_sysimage_path} {run_script_path} {self.ALEAF_abs_path}"

            if self.args.quiet:
                sp = subprocess.check_call(aleaf_cmd,
                                           shell=True,
                                           stdout=open(os.devnull, "wb"))
            else:
                sp = subprocess.check_call(aleaf_cmd, shell=True)

            self.save_ALEAF_outputs()

    def check_for_sysimage_files(self):
        ABCE_sysimage_path = (Path(
            self.settings["ABCE_abs_path"]) /
            self.settings["ABCE_sysimage_file"]
        )

        dispatch_sysimage_path = (Path(
            self.settings["ABCE_abs_path"]) /
            self.settings["dispatch_sysimage_file"]
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
                getattr(row, "cum_exp"),
                getattr(row, "cum_d_x"),
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

    def project_capex(self, unit_type, cum_exp, cum_d_x, rcec, rtec):
        # For a given WIP project, project the sequence of annual capital
        #   expenditures until the project's expected completion
        # TODO: update with more specific methods for different project types
        # For now, assume the project proceeds linearly
        projected_capex = []
        for i in range(math.ceil(round(rtec + 1, 3))):
            projected_capex.append(rcec / rtec)

        return projected_capex

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

        # On the first period, add instruments representing preexisting
        #   debt and equity for the agents
        if self.current_pd < 1:
            inst_id = 1000
            for agent_id, agent_params in self.gc_params.items():
                starting_debt = float(agent_params["starting_debt"])
                debt_frac = agent_params["debt_fraction"]
                starting_equity = float(
                    agent_params["starting_debt"]) / debt_frac * (1 - debt_frac)
                agent_debt_cost = agent_params["cost_of_debt"]
                agent_equity_cost = agent_params["cost_of_equity"]
                debt_row = [agent_id,           # agent_id
                            inst_id,            # instrument_id
                            "debt",             # instrument_type
                            agent_id,
                            # asset_id (agent_id for starting instruments)
                            -1,                 # pd_issued
                            starting_debt,      # initial_principal
                            30,                 # maturity_pd
                            agent_debt_cost     # rate
                            ]
                equity_row = [agent_id,
                              inst_id + 1,
                              "equity",
                              agent_id,
                              -1,
                              starting_equity,
                              30,
                              agent_equity_cost
                              ]
                fin_insts_updates.loc[len(fin_insts_updates.index)] = debt_row
                fin_insts_updates.loc[len(fin_insts_updates.index)] = equity_row

                inst_id += 2

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
            agent_debt_frac = self.gc_params[agent_id]["debt_fraction"]
            agent_debt_cost = self.gc_params[agent_id]["cost_of_debt"]
            agent_equity_cost = self.gc_params[agent_id]["cost_of_equity"]
            amort_pd = 30
            debt_row = [agent_id,                         # agent_id
                        inst_id,                          # instrument_id
                        "debt",                           # instrument_type
                        asset_id,                         # asset_id
                        pd_issued,                        # pd_issued
                        total_qty * agent_debt_frac,      # initial_principal
                        pd_issued + amort_pd,             # maturity_pd
                        agent_debt_cost                   # rate
                        ]
            equity_row = [agent_id,
                          inst_id + 1,
                          "equity",
                          asset_id,
                          pd_issued,
                          total_qty * (1 - agent_debt_frac),
                          pd_issued + amort_pd,
                          agent_equity_cost
                          ]
            fin_insts_updates.loc[len(fin_insts_updates.index)] = debt_row
            fin_insts_updates.loc[len(fin_insts_updates.index)] = equity_row
            inst_id = max(fin_insts_updates["instrument_id"]) + 1

        # Overwrite the financial_instrument_manifest table with the new data
        fin_insts_updates.to_sql(
            "financial_instrument_manifest",
            self.db,
            if_exists="replace",
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
            "agent_financing_schedule",
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
            for agent_id, agent_params in self.gc_params.items():
                summary_asset_id = agent_id
                init_PPE = agent_params["starting_PPE"]
                dep_horiz = 30
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
                asset_PPE = getattr(row, "cum_exp") + getattr(row, "rcec")
                dep_horiz = 20
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
            old_filename = f"{self.ALEAF_scenario_name}__{outfile}.csv"
            old_filepath = Path(self.ALEAF_output_data_path) / old_filename
            new_filename = f"{self.ALEAF_scenario_name}__{outfile}__step_{self.current_pd}.csv"
            new_filepath = Path(self.ABCE_output_data_path) / new_filename
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
            if project_data.loc[0, "rcec"] <= self.settings["large_epsilon"]:
                # Record the project's completion period as the current period
                self.record_completed_xtr_project(project_data)

            # Record updates to the WIP project's status
            self.record_WIP_project_updates(project_data)

            # Record updates to the project's expected completion date in the
            #   database 'assets' table
            self.update_expected_completion_period(project_data)

    def update_agent_debt(self):
        total_new_debt = {201: 0.0, 202: 0.0}

        for agent_id, new_debt in total_new_debt.items():
            if self.current_pd == 0:
                agent_WIP_projects = pd.read_sql_query(
                    f"SELECT asset_id FROM assets WHERE agent_id = {agent_id} AND completion_pd > {self.current_pd} AND retirement_pd > {self.current_pd} AND cancellation_pd > {self.current_pd}",
                    self.db)
            else:
                agent_WIP_projects = pd.read_sql_query(
                    f"SELECT asset_id FROM assets WHERE agent_id = {agent_id} AND completion_pd >= {self.current_pd} AND retirement_pd > {self.current_pd} AND cancellation_pd > {self.current_pd}",
                    self.db)

            for asset_id in agent_WIP_projects.asset_id:
                new_anpe = pd.read_sql_query(
                    f"SELECT anpe FROM WIP_projects WHERE asset_id = {asset_id} AND period = {self.current_pd}",
                    self.db)
                total_new_debt[agent_id] += new_anpe.iloc[0, 0]

            # Update each agent's total debt
            existing_principal = pd.read_sql_query(
                f"SELECT outstanding_principal FROM agent_debt WHERE agent_id = {agent_id} AND period = {self.current_pd-1}",
                self.db)
            starting_debt = self.gc_params[agent_id]["starting_debt"]
            if self.current_pd < 20:
                paid_down = starting_debt * (20. - self.current_pd) / 20.
            else:
                paid_down = starting_debt
            existing_principal = existing_principal.iloc[0, 0]
            total_current_principal = existing_principal - paid_down + \
                total_new_debt[agent_id] * self.gc_params[agent_id]["debt_fraction"]
            debt_row = (agent_id, self.current_pd, total_current_principal)
            self.cur.execute(
                "INSERT INTO agent_debt VALUES (?, ?, ?)", debt_row)
            self.db.commit()

    def record_completed_xtr_project(self, project_data):
        asset_id = project_data.loc[0, "asset_id"]

        # Get asset record from assets
        asset_data = pd.read_sql_query(
            f"SELECT * FROM assets WHERE asset_id = {asset_id}", self.db)

       # Compute periodic sinking fund payments
        unit_type = asset_data.loc[0, "unit_type"]
        unit_life = int(math.ceil(
            self.unit_specs.loc[self.unit_specs.unit_type == unit_type, "unit_life"].values[0]))
        #capex_payment = self.compute_sinking_fund_payment(asset_data.loc[0, "agent_id"], asset_data.loc[0, "cum_exp"], unit_life)
        capex_payment = 0  # to be replaced by capex and financial instrument tracking

        to_update = {"completion_pd": self.current_pd,
                     "retirement_pd": self.current_pd + unit_life,
                     "total_capex": project_data.loc[0, "cum_exp"],
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
                         f"{project_data.cum_d_x.values[0]}, " +
                         f"{project_data.rtec.values[0]}, " +
                         f"{project_data.cum_exp.values[0]}, " +
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
        project_data.loc[0, "cum_exp"] += project_data.loc[0, "anpe"]

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
