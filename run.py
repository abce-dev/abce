from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml

# User inputs
unit_data_file = '/home/kbiegel/gamba/simplified_units.yml'
demand_data_file = '/home/kbiegel/gamba/simplified_demand.yml'

# Run the model
if __name__ == '__main__':
    abce_model = GridModel(1, unit_data_file, demand_data_file)
    for i in range(5):
        abce_model.step()
