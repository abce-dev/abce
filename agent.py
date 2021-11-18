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
                        {self.debt_cost}, {self.equity_cost},
                        {self.interest_cap})""")


    def add_initial_assets_to_db(self, settings):
        # Converter for non-ABCE-matching ALEAF unit type names
        ALEAF_unit_name_converter = {"Solar PV": "PV",
                                     "Steam": "Coal",
                                     "CC": "NGCC",
                                     "CT": "NGCT"}
        # This currently assumes only one agent.
        # Read in the ALEAF system portfolio
        book, writer = ALI.prepare_xlsx_data(self.model.ALEAF_portfolio_ref, self.model.ALEAF_portfolio_ref)
        pdf = ALI.organize_ALEAF_portfolio(writer)
        # Set the initial asset ID
        asset_id = settings["first_asset_id"]
        # Assign all units to this agent, and record each individually in the database
        for unit_type in list(pdf["Unit Type"]):
            for j in range(pdf.loc[pdf["Unit Type"] == unit_type, "EXUNITS"].values[0]):
                abce_unit_type = unit_type
                if unit_type in ALEAF_unit_name_converter.keys():
                    abce_unit_type = ALEAF_unit_name_converter[unit_type]
                asset_dict = {"asset_id": asset_id,
                              "agent_id": self.unique_id,
                              "unit_type": abce_unit_type,
                              "revealed": "true",
                              "completion_pd": 0,
                              "cancellation_pd": 9999,
                              "retirement_pd": 9999,
                              "total_capex": 12345,
                              "cap_pmt": 6789}
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

        # Update the status of each current WIP project
        if self.current_step > 0:
            self.WIP_projects = self.get_WIP_project_list()
            self.update_WIP_projects()

        # Run the agent behavior choice algorithm
        julia_cmd = ("julia -JabceSysimage.so agent_choice.jl ./settings.yml " +
                     f"{self.current_step} {self.unique_id}")
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


    def update_WIP_projects(self):
        db = self.model.db
        cur = db.cursor()
        if self.WIP_projects:
            # If there are active construction projects under way:
            for asset_id in self.WIP_projects:
                # Select the current project's most recent data record
                WIP_project = pd.read_sql(f"SELECT * FROM WIP_projects WHERE " +
                                          f"asset_id = {asset_id} AND " +
                                          f"period = {self.current_step - 1}",
                                          self.model.db)

                # TODO: implement stochastic RCEC and RTEC escalation
                # Currently: no escalation

                if WIP_project.loc[0, "rcec"] - WIP_project.loc[0, "anpe"] <= 0:
                    # This period's authorized expenditures close out the
                    #    remainder of the project. Update the `assets` table to
                    #    reflect completion status.
                    cur.execute(f"UPDATE assets SET completion_pd = {self.current_step} " +
                                f"WHERE asset_id = {asset_id}")

                    # Compute amount of CapEx from the project (in as-spent, inflation-unadjusted $)
                    total_capex = self.compute_total_capex(asset_id)
                    cur.execute(f"UPDATE assets SET total_capex = {total_capex} " +
                                f"WHERE asset_id = {asset_id}")

                    # Compute periodic sinking fund payments
                    unit_type = (pd.read_sql(f"SELECT unit_type FROM assets " +
                                             f"WHERE asset_id = {asset_id}",
                                             self.model.db).loc[0, "unit_type"])
                    unit_type_mask = self.model.unit_specs["unit_type"] == unit_type
                    unit_life = (self.model.unit_specs
                                 .loc[unit_type_mask, "unit_life"]
                                 .values[0])
                    capex_payment = self.compute_sinking_fund_payment(total_capex, unit_life)
                    cur.execute(f"UPDATE assets SET cap_pmt = {capex_payment} " +
                                f"WHERE asset_id = {asset_id}")
                    db.commit()

                # Set values to update the WIP_project dataframe with new completion data
                period = self.current_step
                rcec = max(WIP_project.loc[0, "rcec"] - WIP_project.loc[0, "anpe"], 0)
                rtec = WIP_project.loc[0, "rtec"] - 1
                anpe = 0    # Reset to 0 to avoid inter-period contamination

                # Update the `WIP_projects` database table
                cur.execute(f"INSERT INTO WIP_projects VALUES " +
                            f"({asset_id}, {self.unique_id}, {period}, " +
                            f"{rcec}, {rtec}, {anpe})")

        else:
            # If there are no construction projects in progress:
            print(f"Agent {self.unique_id} has no ongoing construction projects.")

        # Commit the changes to the database
        db.commit()


    def compute_total_capex(self, asset_id):
        """
        Retrieve all previously-recorded capital expenditures for the indicated
        asset from the database, and sum them to return the total capital
        expenditure (up to the current period).

        Args:
           asset_id (int): asset for which to compute capex

        Returns:
           total_capex (float): total capital expenditures up to the present period
        """

        capex_list = pd.read_sql(f"SELECT anpe FROM WIP_projects WHERE " +
                                 f"asset_id = {asset_id}", self.model.db)
        total_capex = sum(capex_list["anpe"])
        return total_capex


    def compute_sinking_fund_payment(self, total_capex, term):
        """
        Compute a constant sinking-fund payment based on a total capital
        expenditures amount and the life of the unit, using the agent's
        financial parameters.

        Args:
           total_capex (float): total capital expenditures on the project
           term (int or float): term over which to amortize capital
             expenditures

        Returns:
           cap_pmt (float): equal capital repayments to make over the course of
             the indicated amortization term
        """

        wacc = (self.debt_fraction * self.debt_cost
                + (1 - self.debt_fraction) * self.equity_cost)
        cap_pmt = total_capex * wacc / (1 - (1 + wacc)**(-term))
        return cap_pmt






