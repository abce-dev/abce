from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml

# User inputs
unit_data_file = '/home/biegelk/abce/data/simplified_units.yml'
market_data_file = '/home/biegelk/abce/data/market.yml'

# Run the model
if __name__ == '__main__':
    abce_model = GridModel(1, unit_data_file, market_data_file)
    for i in range(4):
        abce_model.step()
