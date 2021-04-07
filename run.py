from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml

# User inputs
unit_data_file = "./inputs/unit_specs.csv"
fuel_data_file = "./inputs/fuel_costs.csv"
demand_data_file = "./inputs/demand_data.csv"
#price_curve_data_file = "./inputs/output_DISPATCH_aleafbase.csv"
price_curve_data_file = "./data/RTM_ORDC_REL_DPLY_PRC_ADDR_RSRV_2019.xlsx"
db_file = "./abce_db.db"

# Set default values from which to start assigning agent and asset ID numbers
first_agent_id = 201
first_asset_id = 2001

# Run the model
if __name__ == '__main__':
    abce_model = GridModel(1, db_file, unit_data_file, fuel_data_file, demand_data_file, price_curve_data_file, first_agent_id, first_asset_id)
    for i in range(25):
        abce_model.step()
