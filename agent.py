from mesa import Agent, Model
from mesa.time import RandomActivation
import yaml
import math
import pandas as pd
import subprocess
import csv

# import local modules
import generator as gen
import financial_statement as fs
from ABCEfunctions import *

class GenCo(Agent):
    """ A utility company with a certain number of generation assets.
    """
    def __init__(self, genco_id, model):
        """ Initialize a GenCo class object.

            Detailed Description
            --------------------
            Initializes a GenCo object based on the provided parameters. Uses
            baseline mesa capabilities to set up an Agent-type object.

            Sets up the unique identifier number and the link to the invoking
            GridModel (`model`). During execution, invokes the Julia
            agent_choice.jl script.

            Parameters
            ----------
            genco_id : int
                A unique ID number for the agent, assigned sequentially
                by the Model which owns the agents.
            model : GridModel (mesa)
                A mesa GridModel object which creates and passes data to
                all agents.

        """
        super().__init__(genco_id, model)
        self.assign_parameters('./data/gc_params.yml')
        self.model = model
#        self.fs = fs.AgentFS(model = self.model, agent = self)


    def assign_parameters(self, params_file):
        # Retrieve parameters from file
        params = yaml.load(open(params_file, 'r'), Loader=yaml.FullLoader)
        debt_fraction = params['debt_fraction']
        term_growth_rate = params['terminal_growth_rate']
        discount_rate = params['discount_rate']
        tax_rate = params['tax_rate']
        interest_cap = params['interest_cap']

        # Save parameters to DB table `agent_params`
        db = self.model.db
        cur = db.cursor()
        cur.execute(f"""INSERT INTO agent_params VALUES ({self.unique_id}, 
                        {discount_rate}, {tax_rate}, {term_growth_rate}, 
                        {debt_fraction}, {interest_cap})""")


    def step(self):
        """Controller function to activate all agent behaviors at each time step.
        """
        # Set the current model step
        self.current_step = self.get_current_step()

        # Get lists of all assets, all WIP projects, and all operating assets
        self.all_assets = self.get_current_asset_list()
        self.WIP_projects = self.get_WIP_project_list()
        self.op_assets = self.get_operating_asset_list()

        # Update the status of each current WIP project
        self.update_WIP_assets()

        # Retrieve the demand forecast for the upcoming visible periods
        # TODO: integrate with DB
#        self.get_demand_forecast()

        # Write the excess unserved demand to a csv
#        demand_series = pd.DataFrame({'demand': self.available_demand}) * (-1)
#        demand_file = f"./data/gc{self.unique_id}_demand.csv"

        # Run the agent behavior choice algorithm
#        subprocess.run(["/bin/bash", "-c", "julia -JabceSysimage.so agent_choice.jl"], start_new_session=True)
        # Newer invocation
        sp = subprocess.run([f"julia agent_choice.jl"], shell = True)


        # TODO: integrate FS operations with DB
        # Financial statements currently disabled
#        self.fs.step()


    def get_current_step(self):
        """Obtain the current step number from the model.

        """
        current_step = self.model.current_step
        return current_step


    def get_current_asset_list(self):
        # Get a list of all assets belonging to the current agent where:
        #  - cancellation_pd is in the future (not currently cancelled)
        #  - retirement_pd is in the future (not currently retired)
        db = self.model.db
        cur = db.cursor()
        cur.execute(f"SELECT asset_id FROM assets WHERE agent_id = {self.unique_id} AND cancellation_pd > {self.current_step} AND retirement_pd > {self.current_step}")
        all_asset_list = list(cur.fetchall()[0])
        return all_asset_list


    def get_WIP_project_list(self):
        # Get a list of all assets belonging to the current agent where:
        #  - cancellation_pd is in the future (not currently cancelled)
        #  - completion_pd is in the future (not currently completed)
        db = self.model.db
        cur = db.cursor()
        cur.execute(f"SELECT asset_id FROM assets WHERE agent_id = {self.unique_id} AND completion_pd > {self.current_step} AND cancellation_pd > {self.current_step}")
        WIP_project_list = list(cur.fetchall()[0])
        return WIP_project_list


    def get_operating_asset_list(self):
        # Get a list of all assets belonging to the current agent where:
        #  - cancellation_pd is in the future (not currently cancelled)
        #  - completion_pd is in the past (already completed)
        #  - retirement_pd is in the future (not currently retired)
        db = self.model.db
        cur = db.cursor()
        cur.execute(f"SELECT asset_id FROM assets WHERE agent_id = {self.unique_id} AND completion_pd <= {self.current_step} AND cancellation_pd > {self.current_step} AND retirement_pd > {self.current_step}")
        op_asset_list = list(cur.fetchall()[0])
        return op_asset_list


    def update_WIP_assets(self):
        db = self.model.db
        cur = db.cursor()
        if self.WIP_projects:
            # If there are active construction projects under way:
            print(f"Updating ongoing construction projects for agent {self.unique_id}.")
            for asset_id in self.WIP_projects:
                # Select the current project's most recent data record
                cur.execute(f"SELECT * FROM WIP_projects WHERE asset_id = {asset_id} AND period = {self.current_step - 1}")
                WIP_project = pd.DataFrame([list(cur.fetchall()[0])], columns = ["asset_id", "agent_id", "period", "rcec", "rtec", "anpe"])

            # TODO: implement stochastic RCEC and RTEC escalation
            # Currently: no escalation

            if WIP_project.loc[0, "rcec"] - WIP_project.loc[0, "anpe"] <= 0:
                # This period's authorized expenditures close out the
                #    remainder of the project. Update the `assets` table to
                #    reflect completion status.
                cur.execute(f"UPDATE assets SET completion_pd = {self.current_step} WHERE asset_id = {asset_id}")

            # Update the WIP_project dataframe with new completion data
            WIP_project.loc[0, "period"] = self.current_step
            WIP_project.loc[0, "rcec"] = max(WIP_project.loc[0, "rcec"] - WIP_project.loc[0, "anpe"], 0)
            WIP_project.loc[0, "rtec"] = -= 1
            WIP_project.loc[0, "anpe"] = 0    # Reset to 0 to avoid inter-period contamination

            # Update the `WIP_projects` database table
            vals = (WIP_project.loc[0, "asset_id"], WIP_project.loc[0, "agent_id"], 
                    WIP_project.loc[0, "period"], WIP_project.loc[0, "rcec"],
                    WIP_project.loc[0, "rtec"], WIP_project.loc[0, "anpe"])
            cur.execute(f"INSERT INTO WIP_projects VALUES (?, ?, ?, ?, ?, ?)", vals)

        else:
            # If there are no construction projects in progress:
            print(f"Agent {self.unique_id} hase no ongoing construction projects.")

        # Commit the changes to the database
        db.commit()
















    def get_demand_forecast(self):
        # DEPRECATED : will pull from a demand table
        """ Obtain the current demand visibility window from the model.
        """
        if self.model.demand_NTF != []:
            self.demand_forecast = self.model.demand_NTF
            self.current_demand = self.demand_forecast[0]

