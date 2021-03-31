from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml

# User inputs
unit_data_file = '/home/biegelk/abce/data/unit_specs.csv'
market_data_file = '/home/biegelk/abce/data/market.yml'
db_file = "/home/biegelk/abce/abce_db.db"

# Set default values from which to start assigning agent and asset ID numbers
first_agent_id = 201
first_asset_id = 2001

# Run the model
if __name__ == '__main__':
    abce_model = GridModel(1, db_file, unit_data_file, market_data_file, first_agent_id, first_asset_id)
    for i in range(25):
        abce_model.step()
