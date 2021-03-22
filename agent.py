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
        self.fs = fs.AgentFS(model = self.model, agent = self)


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

        # Write the excess unserved demand to a csv
#        demand_series = pd.DataFrame({'demand': self.available_demand}) * (-1)
#        demand_file = f"./data/gc{self.unique_id}_demand.csv"
        subprocess.run(["/bin/bash", "-c", "julia -JabceSysimage.so unit_choice.jl"], start_new_session=True)
        self.fs.step()

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

