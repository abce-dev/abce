from mesa import Agent, Model
from mesa.time import RandomActivation
import yaml
import math
import pandas as pd
import subprocess
import csv

# import local modules
import ABCEfunctions as ABCE

class GenCo(Agent):
    """ A utility company with a certain number of generation assets.
    """
    def __init__(self, genco_id, model, settings):
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
            settings : dictionary
                A dictionary of runtime/model parameters, loaded and passed
                in by run.py.


        """
        super().__init__(genco_id, model)
        self.model = model
        self.gc_params_file = settings["gc_params_file"]
        self.portfolios_file = settings["portfolios_file"]
        self.assign_parameters(self.gc_params_file)
        self.add_initial_assets_to_db(settings)


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


    def add_initial_assets_to_db(self, settings):
        initial_assets = pd.read_csv(self.portfolios_file, skipinitialspace=True)
        for i in range(len(initial_assets)):
            if initial_assets.loc[i, "agent_id"] == self.unique_id:
                agent_id = initial_assets.loc[i, "agent_id"]
                revealed = "true"
                unit_type = initial_assets.loc[i, "unit_type"]
                completion_pd = 0
                cancellation_pd = 9999
                retirement_pd = initial_assets.loc[i, "useful_life"]
                for j in range(initial_assets.loc[i, "num_copies"]):
                    asset_id = ABCE.get_next_asset_id(self.db, settings["first_asset_id"])
                    total_capex = self.model.unit_specs.loc[self.model.unit_specs["unit_type"] == unit_type, "capacity"].values[0] * self.model.unit_specs.loc[self.model.unit_specs["unit_type"] == unit_type, "uc_x"].values[0] * 1000
                    capital_payment = self.compute_sinking_fund_payment(total_capex, self.model.unit_specs.loc[i, "unit_life"])
                    self.cur.execute(f"""INSERT INTO assets VALUES
                                       ({asset_id}, {agent_id}, '{unit_type}', '{revealed}',
                                        {completion_pd}, {cancellation_pd},
                                        {retirement_pd}, {total_capex}, {capital_payment})""")


    def step(self):
        """Controller function to activate all agent behaviors at each time step.
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
#        subprocess.run(["/bin/bash", "-c", "julia -JabceSysimage.so agent_choice.jl"], start_new_session=True)
        sp = subprocess.check_call([f"julia -JabceSysimage.so agent_choice.jl ./settings.yml {self.current_step} {self.unique_id}"], shell = True)

        print(f"Agent #{self.unique_id}'s turn is complete.\n")


    def get_current_asset_list(self):
        # Get a list of all assets belonging to the current agent where:
        #  - cancellation_pd is in the future (not currently cancelled)
        #  - retirement_pd is in the future (not currently retired)
        all_asset_list = pd.read_sql(f"SELECT asset_id FROM assets WHERE agent_id = {self.unique_id} AND cancellation_pd > {self.current_step} AND retirement_pd > {self.current_step} AND revealed = 'true'", self.model.db)
        all_asset_list = list(all_asset_list["asset_id"])
        return all_asset_list


    def get_WIP_project_list(self):
        # Get a list of all assets belonging to the current agent where:
        #  - cancellation_pd is in the future (not currently cancelled)
        #  - completion_pd is in the future (not currently completed)
        cur = self.model.db.cursor()
        WIP_project_list = pd.read_sql(f"SELECT asset_id FROM assets WHERE agent_id = {self.unique_id} AND completion_pd >= {self.current_step} AND cancellation_pd > {self.current_step} AND revealed = 'true'", self.model.db)
        WIP_project_list = list(WIP_project_list["asset_id"])
        return WIP_project_list


    def get_operating_asset_list(self):
        # Get a list of all assets belonging to the current agent where:
        #  - cancellation_pd is in the future (not currently cancelled)
        #  - completion_pd is in the past (already completed)
        #  - retirement_pd is in the future (not currently retired)
        cur = self.model.db.cursor()
        op_asset_list = pd.read_sql(f"SELECT asset_id FROM assets WHERE agent_id = {self.unique_id} AND completion_pd <= {self.current_step} AND cancellation_pd > {self.current_step} AND retirement_pd > {self.current_step}", self.model.db)
        op_asset_list = list(op_asset_list["asset_id"])
        return op_asset_list


    def update_WIP_projects(self):
        db = self.model.db
        cur = db.cursor()
        if self.WIP_projects:
            # If there are active construction projects under way:
            print(f"Updating ongoing construction projects for agent {self.unique_id}.")
            for asset_id in self.WIP_projects:
                # Select the current project's most recent data record
                WIP_project = pd.read_sql(f"SELECT * FROM WIP_projects WHERE asset_id = {asset_id} AND period = {self.current_step - 1}", self.model.db)

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
                    unit_type = pd.read_sql(f"SELECT unit_type FROM assets WHERE asset_id = {asset_id}", self.model.db).loc[0, "unit_type"]
                    unit_life = self.model.unit_specs.loc[self.model.unit_specs["unit_type"] == unit_type, "unit_life"].values[0]
                    capex_payment = self.compute_sinking_fund_payment(total_capex, unit_life)
                    cur.execute(f"UPDATE assets SET cap_pmt = {capex_payment} WHERE asset_id = {asset_id}")
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
        capex_list = pd.read_sql(f"SELECT anpe FROM WIP_projects WHERE asset_id = {asset_id}", self.model.db)
        total_capex = sum(capex_list["anpe"])
        return total_capex


    def compute_sinking_fund_payment(self, total_capex, unit_life):
        wacc = self.debt_fraction * self.debt_cost + (1 - self.debt_fraction) * self.equity_cost
        cap_pmt = total_capex * wacc / (1 - (1 + wacc)**(-unit_life))
        return cap_pmt






