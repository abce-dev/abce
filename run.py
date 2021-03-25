from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml

# User inputs
unit_data_file = '/home/biegelk/abce/data/unit_specs.csv'
market_data_file = '/home/biegelk/abce/data/market.yml'
db_file = "/home/biegelk/abce/abce_db.db"

# Run the model
if __name__ == '__main__':
    abce_model = GridModel(1, db_file, unit_data_file, market_data_file)
    for i in range(3):
        abce_model.step()
