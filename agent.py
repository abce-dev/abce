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
import logging
from mesa import Agent, Model
from mesa.time import RandomActivation
import yaml
import math
import pandas as pd
import subprocess
import csv
import os
from pathlib import Path

# import local modules
import ABCEfunctions as ABCE
import ALEAF_interface as ALI


class GenCo(Agent):
    """
    A utility company with a certain number of generation assets.
    """

    def __init__(self, genco_id, model, settings, gc_params, cli_args):
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
        """

        super().__init__(genco_id, model)
        self.model = model
        self.settings = settings
        self.args = cli_args
        self.assign_parameters(gc_params)

    def assign_parameters(self, gc_params):
        # Assign all parameters from agent_params as member data
        for key, val in gc_params.items():
            setattr(self, key, val)

        # Save parameters to DB table `agent_params`
        self.db = self.model.db
        self.cur = self.db.cursor()
        self.cur.execute(f"""INSERT INTO agent_params VALUES ({self.unique_id},
                        {self.tax_rate}, {self.terminal_growth_rate},
                        {self.debt_fraction}, {self.cost_of_debt},
                        {self.cost_of_equity}, {self.starting_fcf},
                        {self.starting_debt}, {self.starting_PPE})""")
        self.model.db.commit()

        # Miscellaneous parameters
        self.MW2kW = 1000   # Convert MW to kW

    def step(self):
        """
        Controller function to activate all agent behaviors at each time step.
        """
        # Set the current model step
        self.current_pd = self.model.current_pd
        logging.log(
            self.model.settings["vis_lvl"],
            f"Agent #{self.unique_id} is taking its turn..."
        )

        # Run the agent behavior choice algorithm
        agent_choice_path = (Path(self.settings["ABCE_abs_path"]) /
                             "agent_choice.jl")
        sysimage_cmd = ""

        if self.model.has_ABCE_sysimage:
            sysimage_path = (Path(self.settings["ABCE_abs_path"]) /
                                         self.settings["ABCE_sysimage_file"])
            sysimage_cmd = f"-J {sysimage_path}"

        julia_cmd = (
            f"julia --project={self.settings['ABCE_abs_path']} " + 
            f"{sysimage_cmd} {agent_choice_path} " +
            f"--current_pd={self.current_pd} " +
            f"--agent_id={self.unique_id} " +
            f"--verbosity={self.args.verbosity}"
        )

        sp = subprocess.check_call(julia_cmd, shell=True)

        logging.log(
            self.model.settings["vis_lvl"],
            f"Agent #{self.unique_id}'s turn is complete.\n"
        )

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
                                     f"cancellation_pd > {self.current_pd} " +
                                     f"AND retirement_pd > {self.current_pd}",
                                     self.model.db)
        all_asset_list = list(all_asset_list["asset_id"])
        self.model.db.commit()
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

        WIP_project_list = pd.read_sql(
            f"SELECT asset_id FROM assets WHERE " +
            f"agent_id = {self.unique_id} AND " +
            f"completion_pd >= {self.current_pd} " +
            f"AND cancellation_pd > {self.current_pd}",
            self.model.db)
        WIP_project_list = list(WIP_project_list["asset_id"])
        self.model.db.commit()
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

        op_asset_list = pd.read_sql(f"SELECT asset_id FROM assets WHERE " +
                                    f"agent_id = {self.unique_id} AND " +
                                    f"completion_pd <= {self.current_pd} " +
                                    f"AND cancellation_pd > {self.current_pd} "
                                    f"AND retirement_pd > {self.current_pd}",
                                    self.model.db)
        op_asset_list = list(op_asset_list["asset_id"])
        return op_asset_list
