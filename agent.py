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
    def __init__(self, genco_id, model, existing_portfolio):
        """ Initialize a GenCo class object.

            Detailed Description
            --------------------
            Initializes a GenCo object based on the provided parameters. Uses
            baseline mesa capabilities to set up an Agent-type object.

            Sets up the unique identifier number and the link to the invoking
            GridModel (`model`). Also sets up the agent's starting portfolio.

            Parameters
            ----------
            genco_id : int
                A unique ID number for the agent, assigned sequentially
                by the Model which owns the agents.
            model : GridModel (mesa)
                A mesa GridModel object which creates and passes data to
                all agents.
            existing_portfolio : dict (of dicts)
                A dictionary specifying all assets comprising the agent's
                existing portfolio at time t=0.

            Returns
            -------
            None
        """
        super().__init__(genco_id, model)
        self.assign_parameters('./data/gc_params.yml')
        self.model = model
        self.fs = fs.AgentFS(model = self.model, agent = self)
        self.tax_rate = self.model.tax_rate


    def assign_parameters(self, params_file):
        self.params = yaml.load(open(params_file, 'r'), Loader=yaml.FullLoader)
        self.debt_fraction = self.params['debt_fraction']
        self.term_growth_rate = self.params['terminal_growth_rate']
        self.discount_rate = self.params['discount_rate']
        self.tax_rate = self.params['tax_rate']


    def step(self):
        """Controller function to activate all agent behaviors at each time step.
        """
        self.set_current_step()
        for id_num in self.portfolio:
            self.portfolio[id_num].step()
        num_units = len(self.portfolio.keys())
        print(f'I\'m GenCo # {self.unique_id}, and I have {num_units} units.')
        self.get_demand_forecast()
        self.evaluate_current_capacity()
#        self.build_new_units(self.assess_supply_adequacy())

        # Write the excess unserved demand to a csv
#        demand_series = pd.DataFrame({'demand': self.available_demand}) * (-1)
#        demand_file = f"./data/gc{self.unique_id}_demand.csv"
        subprocess.run(["/bin/bash", "-c", "julia -JabceSysimage.so unit_choice.jl"], start_new_session=True)
        self.fs.step()

    def evaluate_current_capacity(self):
        # DEPRECATED : will be a SELECT capacity FROM assets
        """Sum up all current generation capacity owned by this agent.
        """
        self.total_capacity = 0
        for unit in self.portfolio.keys():
            self.total_capacity += self.portfolio[unit].capacity
        print(f'Total capacity currently installed = {self.total_capacity}.')

    def assess_supply_adequacy(self):
        # DEPRECATED : will be a DB operation
        # Also has some dummy behavior which needs to be cleared
        """Determine whether future demand exceeds future supply.
        """
        supply_surplus = list(self.fs.fsdata['capacity'].iloc[self.current_step+1:self.current_step + len(self.demand_forecast)+1] - self.demand_forecast)
        self.available_demand = supply_surplus
        if not all(s > 0 for s in supply_surplus):
            num_new_units = int(math.ceil((-min(supply_surplus) / float(self.model.unit_data.loc['unit_1', 'capacity']))))
        else:
            num_new_units = 0
        return num_new_units

    def build_new_units(self, new_units):
        # DEPRECATED : will be handled by Julia
        """Start a new construction project, and add to agent's portfolio.
        """
        for i in range(new_units):
            unit_id = self.model.id_register.register_unit(agent_id = self.unique_id)
            new_unit = gen.Generator(world_model=self.model, agent=self, id_num=unit_id, gtype='unit_1')
            self.portfolio[unit_id] = new_unit

    def set_current_step(self):
        """Obtain the current step number from the model.

        """
        self.current_step = self.model.current_step
        print(f'\n\n Current step = {self.current_step}')

    def get_demand_forecast(self):
        # DEPRECATED : will pull from a demand table
        """ Obtain the current demand visibility window from the model.
        """
        if self.model.demand_NTF != []:
            self.demand_forecast = self.model.demand_NTF
            self.current_demand = self.demand_forecast[0]

