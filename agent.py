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

import numpy as np
from mesa import Agent, Model
from mesa.time import RandomActivation
import yaml
import math
import pandas as pd
import subprocess
import csv
import os

# import local modules
import ABCEfunctions as ABCE
import ALEAF_interface as ALI

class GenCo(Agent):
    """ 
    A utility company with a certain number of generation assets.
    """

    def __init__(self, genco_id, model, settings, cli_args):
        """
        Initialize a GenCo class object.

        Detailed Description:
          Initializes a GenCo object based on the provided parameters. Uses
          baseline mesa capabilities to set up an Agent-type object.

          Sets up the unique identifier number and the link to the invoking
          GridModel (`model`). During execution, invokes the Julia
          agent_choice.jl script.

        Args:
           genco_id (int): A unique ID number for the agent, assigned
             sequentially by the Model which owns the agents.
           model (mesa GridModel): Container object which creates, updates,
             and passes data to all agents.
           settings (dict): Runtime/model parameters, loaded and passed
             in by run.py.
           quiet (bool): Set from CLI; sets verbosity level.
        """

        super().__init__(genco_id, model)
        self.model = model
        self.quiet = cli_args.quiet
        self.settings = settings
        self.assign_parameters(settings)
        self.add_initial_assets_to_db(settings)


    def assign_parameters(self, settings):
        # Retrieve parameters from file
        gc_params_file_name = os.path.join(settings["ABCE_abs_path"],
                                           settings["gc_params_file"])
        params = yaml.load(open(gc_params_file_name, 'r'), Loader=yaml.FullLoader)
        # Assign all parameters from params as member data
        for key, val in params.items():
            setattr(self, key, val)

        # Save parameters to DB table `agent_params`
        self.db = self.model.db
        self.cur = self.db.cursor()
        self.cur.execute(f"""INSERT INTO agent_params VALUES ({self.unique_id},
                        {self.discount_rate}, {self.tax_rate},
                        {self.terminal_growth_rate}, {self.debt_fraction},
                        {self.cost_of_debt}, {self.cost_of_equity},
                        {self.interest_cap})""")

        # Miscellaneous parameters
        self.MW2kW = 1000   # Convert MW to kW


    def add_initial_assets_to_db(self, settings):
        # Get supplemental ABCE unit data
        unit_spec_ABCE = pd.read_csv(os.path.join(settings["ABCE_abs_path"], settings["unit_specs_abce_supp_file"]))
        # This currently assumes only one agent.
        # Read in the ALEAF system portfolio
        book, writer = ALI.prepare_xlsx_data(self.model.ALEAF_portfolio_ref, self.model.ALEAF_portfolio_ref)
        pdf = ALI.organize_ALEAF_portfolio(writer)
        # Set the initial asset ID
        asset_id = settings["first_asset_id"]

        # Retrieve the dataframe of mandatory unit retirement periods for
        #   this agent
        ret_df = self.get_agent_retirement_data(settings)

        # Assign all units to this agent, and record each individually in the database
        for unit_type in list(pdf["Unit Type"]):
            # Retrieve and the list of retirement period data for this
            #   unit type and agent, and create the cumulative-sum threshold
            #   mapping
            unit_rets = self.create_unit_type_retirement_df(ret_df, unit_type, asset_id)

            for j in range(pdf.loc[pdf["Unit Type"] == unit_type, "EXUNITS"].values[0]):
                # Assign the appropriate retirement period
                retirement_pd = self.assign_retirement_pd(unit_rets, asset_id)

                # Compute unit capex according to its unit type specification
                unit_capex = self.compute_total_capex_preexisting(unit_type)
                # Default: assume all existing assets have 0 outstanding
                #   financing balance. To add obligations for capital payments
                #   on pre-existing assets, uncomment the code lines below.
                # Compute the asset's annual capital payment, as if the asset
                #   were not paid off. Use the asset's unit_life value as the
                #   repayment term for financing.
                #unit_life = self.model.unit_specs.loc[unit_type, "unit_life"]
                #unit_cap_pmt = self.model.compute_sinking_fund_payment(self.unique_id, unit_capex, unit_life)

                # Create the dictionary of asset characteristics
                asset_dict = {"asset_id": asset_id,
                              "agent_id": self.unique_id,
                              "unit_type": unit_type,
                              "revealed": "true",
                              "completion_pd": 0,
                              "cancellation_pd": 9999,
                              "retirement_pd": retirement_pd,
                              "total_capex": unit_capex,
                              "cap_pmt": 0}
                new_asset = pd.DataFrame(asset_dict, index=[0])

                # Use the first asset's DataFrame as the master
                if asset_id == settings["first_asset_id"]:
                    master_asset_df = new_asset
                else:
                    master_asset_df = master_asset_df.append(new_asset)
                # Increment to get the asset_id
                asset_id += 1
        master_asset_df.to_sql("assets", self.db, if_exists="replace", index=False)


    def step(self):
        """
        Controller function to activate all agent behaviors at each time step.
        """
        # Set the current model step
        self.current_step = self.model.current_step
        print(f"Agent #{self.unique_id} is taking its turn...")

        # Get lists of all assets, all WIP projects, and all operating assets
        self.all_assets = self.get_current_asset_list()
        self.op_assets = self.get_operating_asset_list()

        # Run the agent behavior choice algorithm
        agent_choice_path = os.path.join(self.settings["ABCE_abs_path"],
                                         "agent_choice.jl")
        sysimage_path = os.path.join(self.settings["ABCE_abs_path"],
                                     "abceSysimage.so")
        julia_cmd = (f"julia -J{sysimage_path} {agent_choice_path} " +
                     f"--settings_file={self.model.settings_file_name.name} " +
                     f"--current_pd={self.current_step} " +
                     f"--agent_id={self.unique_id}")
        if self.quiet:
            sp = subprocess.check_call([julia_cmd],
                                       shell = True,
                                       stdout=open(os.devnull, "wb"))
        else:
            sp = subprocess.check_call([julia_cmd], shell = True)

        print(f"Agent #{self.unique_id}'s turn is complete.\n")


    def get_current_asset_list(self):
        """
        Get a list of all assets meeting the following conditions:
          - belongs to the currently-active agent
          - `cancellation_pd` is in the future (not currently cancelled)
          - `retirement_pd` is in the future (not currently retired)

        Returns:
           all_asset_list (list of ints): all asset IDs meeting the above criteria
        """

        all_asset_list = pd.read_sql(f"SELECT asset_id FROM assets WHERE " +
                                     f"agent_id = {self.unique_id} AND " +
                                     f"cancellation_pd > {self.current_step} " +
                                     f"AND retirement_pd > {self.current_step}",
                                      self.model.db)
        all_asset_list = list(all_asset_list["asset_id"])
        return all_asset_list


    def get_WIP_project_list(self):
        """
        Get a list of all assets meeting the following conditions:
          - belongs to the currently-active agent
          - completion_pd is in the future (not yet completed)
          - cancellation_pd is in the future (not currently cancelled)

        Returns:
           WIP_project_list (list of ints): all asset IDs meeting the
             above criteria
        """

        cur = self.model.db.cursor()
        WIP_project_list = pd.read_sql(f"SELECT asset_id FROM assets WHERE " +
                                       f"agent_id = {self.unique_id} AND " +
                                       f"completion_pd >= {self.current_step} " +
                                       f"AND cancellation_pd > {self.current_step}",
                                       self.model.db)
        WIP_project_list = list(WIP_project_list["asset_id"])
        return WIP_project_list


    def get_operating_asset_list(self):
        """
        Get a list of all assets meeting the following conditions:
          - belongs to the currently-active agent
          - completion_pd is in the past (already completed)
          - cancellation_pd is in the future (never cancelled)
          - retirement_pd is in the future (not currently retired)

        Returns:
           op_asset_list (list of ints): all asset IDs meeting the above criteria
        """

        cur = self.model.db.cursor()
        op_asset_list = pd.read_sql(f"SELECT asset_id FROM assets WHERE " +
                                    f"agent_id = {self.unique_id} AND " +
                                    f"completion_pd <= {self.current_step} " +
                                    f"AND cancellation_pd > {self.current_step} "
                                    f"AND retirement_pd > {self.current_step}",
                                    self.model.db)
        op_asset_list = list(op_asset_list["asset_id"])
        return op_asset_list


    def compute_total_capex_preexisting(self, unit_type):
        """
        Compute total capex for units which already exist at the start of the
          run. These units are currently assumed to cost exactly their
          type's estimated construction cost, per the unit_specs data.

        Args:
          unit_type (str): unit type corresponding with a type given in the
            A-LEAF specification

        Returns:
          total_capex (float): total capital expenditures for the asset
        """
        unit_data = self.model.unit_specs.loc[unit_type, :]
        unit_cost_per_kW = unit_data["uc_x"]    # $/kW
        unit_capacity = unit_data["capacity"]   # MW

        total_capex = unit_cost_per_kW * unit_capacity * self.MW2kW

        return total_capex


    def get_agent_retirement_data(self, settings):
        """
        Retrieve the dataframe of mandatory unit retirement periods for
          this agent

        Args:
          settings: settings dictionary

        Returns:
          ret_df: dataframe of all unit-type retirement data for this agent
        """
        ret_data_file = os.path.join(settings["ABCE_abs_path"],
                                     settings["retirement_period_specs_file"])
        ret_df = pd.read_csv(ret_data_file, comment="#")
        ret_df = ret_df[ret_df["agent_id"] == self.unique_id].copy()

        return ret_df


    def create_unit_type_retirement_df(self, ret_df, unit_type, next_asset_id):
        """
        Create the step-function mapping to determine which units receive
        which mandatory retirement periods.
        """
        # Filter the df of retirement period data for the current unit type
        unit_rets = ret_df[ret_df["unit_type"] == unit_type].copy()

        # Sort from lowest to highest retirement period
        unit_rets = unit_rets.sort_values(by="retirement_pd", axis=0, ascending=True, ignore_index=True)

        # Generate the thresholds for each retirement period, starting with
        #   the next available asset id. These thresholds indicate the
        #   END (NON-INCLUSIVE) of each assignment interval. That is, a row
        #   with an rp_threshold of 2005 and a retirement_pd of 12 means that
        #   any assets with IDs strictly less than 2005 will receive a
        #   retirement_pd of 12 (unless they qualify for a lower threshold
        #   category).
        unit_rets["rp_threshold"] = np.cumsum(unit_rets["num_copies"].to_numpy(), axis=0) + next_asset_id

        return unit_rets


    def assign_retirement_pd(self, unit_rets, current_asset_id):
        """
        Compare the current asset ID against the classification thresholds in
          unit_rets, to determine the asset's mandatory retirement period.

        Args:
          unit_rets (DataFrame): retirement threshold mapping for current agent and
            unit type
          current_asset_id (int): unique identifier for the current asset

        Returns:
          retirement_pd (int): mandatory retirement period for this asset
        """

        # If no retirement period is specified, set the asset to never retire
        #    (i.e. retirement period is a Large Number)
        retirement_pd = 9999

        for i in range(len(unit_rets)):
            if current_asset_id < unit_rets.loc[i, "rp_threshold"]:
                retirement_pd = unit_rets.loc[i, "retirement_pd"]
                break   # Once the correct bin is found, exit the loop

        return retirement_pd

