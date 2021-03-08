from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
import yaml
import pandas as pd

# import local modules
import id_register
#import demand_history

class GridModel(Model):
    ''' A model with some number of GenCos. '''
    def __init__(self, n, unit_data_file, demand_data_file):
        # Parameters
        self.num_agents = n
        self.current_step = -1
        self.future_vis = 4  # Number of time steps agents can see into the future
        self.tax_rate = 0.265 # Corporate tax rate

        # Load default unit data and demand profile
        self.load_unit_data(unit_data_file)
        self.get_true_market_data(demand_data_file)
        self.set_demand_visibility_window()

        # Set up a public register of unit IDs
        self.id_register = id_register.ID_register()

        # Set up a public demand-history data series
        #self.demand_history = demand_history.DemandHistory(self.prev_demand[0])

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Load the agents' starting portfolios
        c_portfolio = pd.read_csv('./data/portfolios.csv', index_col='agent', skipinitialspace=True)

        # Create agents
        for i in range(self.num_agents):
            gc = GenCo(i, self, c_portfolio.loc[i+1])
            self.schedule.add(gc)

    def load_unit_data(self, filename):
        unit_file = open(filename)
        self.unit_data = yaml.load(unit_file, Loader=yaml.FullLoader)
        self.unit_data = pd.DataFrame.from_dict(self.unit_data, orient='index')
        self.unit_types = self.unit_data.index

    def get_true_market_data(self, filename):
        market_file = open(filename)
        market_data = yaml.load(market_file, Loader=yaml.FullLoader)
        self.prev_demand = market_data['past_demand']
        self.true_demand_profile = market_data['demand']
        self.eprices = market_data['prices']

    def set_demand_visibility_window(self):
        self.demand_NTF = self.true_demand_profile[self.current_step:self.current_step + self.future_vis]

    def step(self):
        ''' Advance the model by one step. '''
        # Move next future demand data point to the historical register
        #if self.current_step != -1:
        #    self.demand_history.add_data(self.demand_NTF[0])
        self.current_step += 1
        self.set_demand_visibility_window()
        self.schedule.step()

