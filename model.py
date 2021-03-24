from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
import yaml
import pandas as pd
import subprocess

# import local modules
import id_register
from ABCEfunctions import *

class GridModel(Model):
    ''' A model with some number of GenCos. '''
    def __init__(self, n, db_file, unit_data_file, demand_data_file):
        # Parameters
        self.num_agents = n
        self.current_step = -1
        self.future_vis = 4  # Number of time steps agents can see into the future

        # Initialize database for managing asset and WIP construction project data
        self.db_file = db_file
        sp = subprocess.run(["python3 seed_creator.py"], shell=True)
        self.db, self.cur = load_database(self.db_file)

        # Load default unit data and demand profile
        self.unit_types, self.unit_data = self.load_unit_data(unit_data_file)
        self.get_true_market_data(demand_data_file)
        self.set_demand_visibility_window()

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Create agents
        for i in range(self.num_agents):
            gc = GenCo(i, self)
            self.schedule.add(gc)

        # Load the agents' starting portfolios
        initial_assets = pd.read_csv('./data/portfolios.csv', skipinitialspace=True)

        # Add all initial assets to the database
        add_initial_assets_to_db(initial_assets, self.cur, self.db)


    def load_unit_data(self, filename):
        unit_file = open(filename)
        unit_data = yaml.load(unit_file, Loader=yaml.FullLoader)
        unit_data = pd.DataFrame.from_dict(unit_data, orient='index')
        unit_types = self.unit_data.index
        return unit_types, unit_data


    def get_true_market_data(self, filename):
        market_file = open(filename)
        market_data = yaml.load(market_file, Loader=yaml.FullLoader)
        self.prev_demand = market_data['past_demand']
        self.true_demand_profile = market_data['demand']
        self.eprices = market_data['prices']


    def set_demand_visibility_window(self):
        self.demand_NTF = self.true_demand_profile[self.current_step:self.current_step + self.future_vis]


    def add_initial_assets_to_db(initial_assets, cur, db):
        for i in range(len(initial_assets)):
            for j in range(initial_assets.loc[i, "num_copies"]):
                asset_id = get_next_asset_id(self.db, self.cur)
                agent_id = initial_assets.loc[i, "agent_id"]
                unit_type = initial_assets.loc[i, "unit_type"]
                completion_pd = 0
                cancellation_pd = 9999
                retirement_pd = initial_assets.loc[i, "useful_life"]
                capital_payment = self.unit_data.loc[i, "overnight_cost"] / initial_assets.loc[i, "useful_life"]
                vals = (asset_id, agent_id, completion_pd, cancellation_pd, retirement_pd, capital_payment)
                cur.execute("""INSERT INTO assets VALUES
                               ({asset_id}, {agent_id}, {unit_type},
                                {completion_pd}, {cancellation_pd},
                                {retirement_pd}, {capital_payment})""")


    def step(self):
        ''' Advance the model by one step. '''
        self.current_step += 1
        self.set_demand_visibility_window()
        self.schedule.step()

