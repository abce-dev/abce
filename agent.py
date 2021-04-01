from mesa import Agent, Model
from mesa.time import RandomActivation
import yaml
import math
import pandas as pd
import subprocess
import csv

# import local modules
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
        self.debt_fraction = params['debt_fraction']
        self.debt_cost = params['debt_cost']
        self.equity_cost = params['equity_cost']
        self.term_growth_rate = params['terminal_growth_rate']
        self.discount_rate = params['discount_rate']
        self.tax_rate = params['tax_rate']
        self.interest_cap = params['interest_cap']

        # Save parameters to DB table `agent_params`
        self.db = self.model.db
        self.cur = self.db.cursor()
        self.cur.execute(f"""INSERT INTO agent_params VALUES ({self.unique_id}, 
                        {self.discount_rate}, {self.tax_rate}, {self.term_growth_rate},
                        {self.debt_fraction}, {self.debt_cost}, {self.equity_cost},
                        {self.interest_cap})""")


    def step(self):
        """Controller function to activate all agent behaviors at each time step.
        """
        # Set the current model step
        self.current_step = self.model.current_step

        # Get lists of all assets, all WIP projects, and all operating assets
        self.all_assets = self.get_current_asset_list()
        self.op_assets = self.get_operating_asset_list()

        # Update the status of each current WIP project
        if self.current_step > 0:
            self.WIP_projects = self.get_WIP_project_list()
            print(get_table(self.model.db, self.model.cur, "assets"))
            self.update_WIP_projects()

        # Run the agent behavior choice algorithm
#        subprocess.run(["/bin/bash", "-c", "julia -JabceSysimage.so agent_choice.jl"], start_new_session=True)
        sp = subprocess.check_call([f"julia -JabceSysimage.so agent_choice.jl ./abce_db.db {self.current_step} {self.unique_id}"], shell = True)
        print(get_table(self.model.db, self.model.cur, "WIP_projects"))


    def get_current_asset_list(self):
        # Get a list of all assets belonging to the current agent where:
        #  - cancellation_pd is in the future (not currently cancelled)
        #  - retirement_pd is in the future (not currently retired)
        cur = self.model.db.cursor()
        cur.execute(f"SELECT asset_id FROM assets WHERE agent_id = {self.unique_id} AND cancellation_pd > {self.current_step} AND retirement_pd > {self.current_step}")
        all_asset_list = [int(item[0]) for item in list(cur.fetchall())]
        return all_asset_list


    def get_WIP_project_list(self):
        # Get a list of all assets belonging to the current agent where:
        #  - cancellation_pd is in the future (not currently cancelled)
        #  - completion_pd is in the future (not currently completed)
        cur = self.model.db.cursor()
        cur.execute(f"SELECT asset_id FROM assets WHERE agent_id = {self.unique_id} AND completion_pd > {self.current_step - 1} AND cancellation_pd > {self.current_step}")
        WIP_project_list = list(cur.fetchall())
        WIP_project_list = [int(item[0]) for item in WIP_project_list]
        return WIP_project_list


    def get_operating_asset_list(self):
        # Get a list of all assets belonging to the current agent where:
        #  - cancellation_pd is in the future (not currently cancelled)
        #  - completion_pd is in the past (already completed)
        #  - retirement_pd is in the future (not currently retired)
        cur = self.model.db.cursor()
        cur.execute(f"SELECT asset_id FROM assets WHERE agent_id = {self.unique_id} AND completion_pd <= {self.current_step} AND cancellation_pd > {self.current_step} AND retirement_pd > {self.current_step}")
        op_asset_list = list(cur.fetchall())
        op_asset_list = [int(item[0]) for item in op_asset_list]
        return op_asset_list


    def update_WIP_projects(self):
        db = self.model.db
        cur = db.cursor()
        if self.WIP_projects:
            # If there are active construction projects under way:
            print(f"Updating ongoing construction projects for agent {self.unique_id}.")
            for asset_id in self.WIP_projects:
                # Select the current project's most recent data record
                cur.execute(f"SELECT * FROM WIP_projects WHERE asset_id = {asset_id} AND period = {self.current_step - 1}")
                WIP_project = pd.DataFrame([list(cur.fetchall()[0])], columns = ["asset_id", "agent_id", "period", "rcec", "rtec", "anpe"])
                db.commit()

                # TODO: implement stochastic RCEC and RTEC escalation
                # Currently: no escalation

                if WIP_project.loc[0, "rcec"] - WIP_project.loc[0, "anpe"] <= 0:
                    # This period's authorized expenditures close out the
                    #    remainder of the project. Update the `assets` table to
                    #    reflect completion status.
                    cur.execute(f"UPDATE assets SET completion_pd = {self.current_step} WHERE asset_id = {asset_id}")

                    # Compute amount of CapEx from the project (in as-spent, inflation-unadjusted $)
                    total_capex = self.compute_total_capex(asset_id)
                    cur.execute(f"UPDATE assets SET total_capex = {total_capex} WHERE asset_id = {asset_id}")

                    # Compute periodic sinking fund payments
                    capex_repayment = self.compute_sinking_fund_payment(total_capex)
                    cur.execute(f"UPDATE assets SET cap_pmt = {capex_repayment} WHERE asset_id = {asset_id}")
                    db.commit()

                # Set values to update the WIP_project dataframe with new completion data
                period = self.current_step
                rcec = max(WIP_project.loc[0, "rcec"] - WIP_project.loc[0, "anpe"], 0)
                rtec = WIP_project.loc[0, "rtec"] - 1
                anpe = 0    # Reset to 0 to avoid inter-period contamination

                # Update the `WIP_projects` database table
                cur.execute(f"INSERT INTO WIP_projects VALUES ({asset_id}, {self.unique_id}, {period}, {rcec}, {rtec}, {anpe})")

        else:
            # If there are no construction projects in progress:
            print(f"Agent {self.unique_id} has no ongoing construction projects.")

        # Commit the changes to the database
        db.commit()


    def compute_total_capex(self, asset_id):
        self.cur.execute(f"SELECT anpe FROM WIP_projects WHERE asset_id = {asset_id}")
        capex_list = [item[0] for item in list(self.cur.fetchall())]
        total_capex = sum(capex_list)
        return total_capex


    def compute_sinking_fund_payment(self, total_capex):
        t = 40    # Temporary assumption: new assets paid for with 40-year debt life
        wacc = self.debt_fraction * self.debt_cost + (1 - self.debt_fraction) * self.equity_cost
        cap_pmt = total_capex * wacc / ((1 + wacc)**t - 1)
        return cap_pmt






