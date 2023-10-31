##########################################################################
## Copyright 2023 Argonne National Laboratory
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


class GenCo(Agent):
    """
    A utility company with a certain number of generation assets.
    """

    def __init__(self, genco_id, model, gc_params):
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
        self.assign_parameters(gc_params)

    def assign_parameters(self, gc_params):
        # Retrieve all parameter names, the union of the database table column
        #   headers and the items provided in the agent_specifications file
        cur = self.model.db.cursor()
        cur.execute("SELECT * FROM agent_params")
        agent_params_db_cols = set([element[0] for element in cur.description])
        agent_params_fields = set(gc_params.keys())

        # Finalize the list of all possible parameter names
        all_agent_params = list(agent_params_db_cols | agent_params_fields)

        # Remove "agent_id" from this list, as it's already set in __init__()
        all_agent_params.remove("agent_id")

        # Set up default values for fields which are optional to specify
        #   for any agent
        universal_optional_fields = {
            "starting_portfolio": {},
            "scheduled_retirements": {},
            "inactive": False,
        }

        # Start by ensuring that the agent has attributes set for all optional
        #   fields, using the defaults in optional_fields if necessary
        for key, default_value in universal_optional_fields.items():
            if key in gc_params.keys():
                setattr(self, key, gc_params[key])
            else:
                setattr(self, key, default_value)

        # Iterate through all remaining parameters in all_agent_params.
        # An inactive agent is not required to have financial parameters
        #   specified, so these can be filled in with zeros.
        # Other than these two cases, unavailable parameters should raise
        #   an error.
        for param in all_agent_params:
            if param in gc_params.keys():
                setattr(self, param, gc_params[param])
            elif self.inactive:
                setattr(self, param, 0)
            else:
                raise ValueError(
                    f"Agent #{self.unique_id} is missing required parameter "
                    + f"{param} from its specification."
                )

        # Save parameters to DB table `agent_params`
        self.db = self.model.db
        self.cur = self.db.cursor()
        self.cur.execute(
            f"""INSERT INTO agent_params VALUES ({self.unique_id},
                        {self.debt_fraction}, {self.cost_of_debt},
                        {self.cost_of_equity}, {self.starting_debt},
                        {self.starting_PPE}, {self.starting_RE})"""
        )
        self.model.db.commit()

    def step(self):
        """
        Controller function to activate all agent behaviors at each time step.
        """
        # If this agent is designated as inactive, this function should return
        #   immediately
        if self.inactive:
            return

        logging.log(
            self.model.settings["constants"]["vis_lvl"],
            f"Agent #{self.unique_id} is taking its turn...",
        )

        # Run the agent behavior choice algorithm
        agent_choice_path = (
            Path(self.model.settings["file_paths"]["ABCE_abs_path"])
            / "src"
            / "agent_choice.jl"
        )

        sysimage_cmd = ""

        if self.model.has_ABCE_sysimage:
            sysimage_path = (
                Path(self.model.settings["file_paths"]["ABCE_abs_path"])
                / "env"
                / self.model.settings["file_paths"]["ABCE_sysimage_file"]
            )
            sysimage_cmd = f"-J {sysimage_path}"

        julia_env = Path(os.getenv("ABCE_ENV"))

        julia_cmd = (
            f"julia --project={julia_env} "
            + f"{sysimage_cmd} {agent_choice_path} "
            + f"--current_pd={self.model.current_pd} "
            + f"--agent_id={self.unique_id} "
            + f"--verbosity={self.model.args.verbosity} "
            + f"--settings_file={self.model.args.settings_file} "
            + f"--inputs_path={self.model.args.inputs_path} "
            + f"--abce_abs_path={self.model.settings['file_paths']['ABCE_abs_path']} "
        )

        sp = subprocess.check_call(julia_cmd, shell=True)

        logging.log(
            self.model.settings["constants"]["vis_lvl"],
            f"Agent #{self.unique_id}'s turn is complete.\n",
        )

    def get_current_asset_list(self):
        """
        Get a list of all assets meeting the following conditions:
          - belongs to the currently-active agent
          - `cancellation_pd` is in the future (not currently cancelled)
          - `retirement_pd` is in the future (not currently retired)

        Returns:
           all_asset_list (list of ints): all asset IDs meeting the above
               criteria
        """

        all_asset_list = pd.read_sql(
            f"SELECT asset_id FROM assets WHERE "
            + f"agent_id = {self.unique_id} AND "
            + f"cancellation_pd > {self.model.current_pd} "
            + f"AND retirement_pd > {self.model.current_pd}",
            self.model.db,
        )
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
            f"SELECT asset_id FROM assets WHERE "
            + f"agent_id = {self.unique_id} AND "
            + f"completion_pd >= {self.model.current_pd} "
            + f"AND cancellation_pd > {self.model.current_pd}",
            self.model.db,
        )
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
           op_asset_list (list of ints): all asset IDs meeting the
               above criteria
        """

        op_asset_list = pd.read_sql(
            f"SELECT asset_id FROM assets WHERE "
            + f"agent_id = {self.unique_id} AND "
            + f"completion_pd <= {self.model.current_pd} "
            + f"AND cancellation_pd > {self.model.current_pd} "
            f"AND retirement_pd > {self.model.current_pd}",
            self.model.db,
        )
        op_asset_list = list(op_asset_list["asset_id"])
        return op_asset_list
