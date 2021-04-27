from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml
import pandas as pd

settings_file = "./settings.yml"
with open(settings_file) as setfile:
    settings = yaml.load(setfile, Loader=yaml.FullLoader)

num_steps = settings["num_steps"]
output_file = settings["results_file"]

# Run the model
if __name__ == '__main__':
    abce_model = GridModel(settings_file)
    for i in range(num_steps):
        abce_model.step()

db_tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table';", abce_model.db)
with pd.ExcelWriter(output_file) as writer:
    for i in range(len(db_tables)):
        table = db_tables.loc[i, "name"]
        df = pd.read_sql_query(f"SELECT * FROM {table}", abce_model.db)
        df.to_excel(writer, sheet_name=f"{table}", engine="openpyxl")
