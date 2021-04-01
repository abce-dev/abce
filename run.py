from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml

# User inputs
unit_data_file = "./data/unit_specs.csv"
fuel_data_file = "./data/fuel_costs.csv"
demand_data_file = "./data/demand_data.csv"
db_file = "./abce_db.db"

# Set default values from which to start assigning agent and asset ID numbers
first_agent_id = 201
first_asset_id = 2001

# Run the model
if __name__ == '__main__':
    abce_model = GridModel(1, db_file, unit_data_file, fuel_data_file, demand_data_file, first_agent_id, first_asset_id)
    for i in range(25):
        abce_model.step()
