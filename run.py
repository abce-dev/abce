from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml

# User inputs
unit_data_file = '/home/kbiegel/abce/simplified_units.yml'
market_data_file = '/home/kbiegel/abce/market.yml'

# Run the model
if __name__ == '__main__':
    abce_model = GridModel(1, unit_data_file, market_data_file)
    for i in range(4):
        abce_model.step()
