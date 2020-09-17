from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
import yaml

portfolio = {100: {'gtype': 'unit_1'}, 101: {'gtype': 'unit_1'}}

class GridModel(Model):
    ''' A model with some number of GenCos. '''
    def __init__(self, n, unit_data_file, demand_data_file):
        self.num_agents = n
        self.current_step = -1
        self.future_vis = 3  # Number of time steps agents can see into the future
        self.schedule = RandomActivation(self)
        self.load_unit_data(unit_data_file)
        self.set_true_demand_profile(demand_data_file)
        self.set_demand_visibility_window()
        # Create agents
        for i in range(self.num_agents):
            gc = GenCo(i, self, portfolio)
            self.schedule.add(gc)

    def load_unit_data(self, filename):
        unit_file = open(filename)
        self.unit_data = yaml.load(unit_file, Loader=yaml.FullLoader)
        self.unit_types = list(self.unit_data.keys())
        self.unit_data_types = list(self.unit_data[self.unit_types[0]].keys())

    def set_true_demand_profile(self, filename):
        demand_file = open(filename)
        demand = yaml.load(demand_file, Loader=yaml.FullLoader)
        self.true_demand_profile = demand['demand']

    def set_demand_visibility_window(self):
        self.demand_NTF = self.true_demand_profile[self.current_step:self.current_step + self.future_vis]

    def step(self):
        ''' Advance the model by one step. '''
        self.current_step += 1
        self.set_demand_visibility_window()
        self.schedule.step()

