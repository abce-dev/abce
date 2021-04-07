from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml

settings_file = "./settings.yml"
with open(settings_file) as setfile:
    settings = yaml.load(setfile, Loader=yaml.FullLoader)

num_steps = settings["num_steps"]

# Run the model
if __name__ == '__main__':
    abce_model = GridModel(settings_file)
    for i in range(num_steps):
        abce_model.step()
