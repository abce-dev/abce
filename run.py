from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml
import pandas as pd
import argparse


def read_settings(settings_file):
    """
    Read in settings from the settings file.
    """
    settings = yaml.load(settings_file, Loader=yaml.FullLoader)
    return settings


def cli_args():
    """
    Set up the command-line argument parser. Then, read and parse whatever
    arguments are provided via sys.argv into the argument types defined below.

    New command-line options can be specified here.

    Returns:
       args (argparse object): populated namespace, with argument strings as
         attributes. Retrieve values with args.<argument_name>.
    """
    parser = argparse.ArgumentParser(description='Run an ABCE simulation.')
    parser.add_argument("--force", "-f",
                          action="store_true",
                          help="Agree to overwrite any existing DB files.")
    parser.add_argument("--settings_file",
                          type=argparse.FileType("r"),
                          help="Simulation settings file name.",
                          default="./settings.yml")
    parser.add_argument("--silent", "-s",
                          action="store_true",
                          help="Suppress all output except the turn and period counters.")
    args = parser.parse_args()
    return args


def run_model():
    
    args = cli_args()

    settings = read_settings(args.settings_file)

    abce_model = GridModel(settings, args)
    for i in range(settings["num_steps"]):
        abce_model.step()

    db_tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table';", abce_model.db)
    with pd.ExcelWriter(settings["output_file"]) as writer:
        for i in range(len(db_tables)):
            table = db_tables.loc[i, "name"]
            final_db = pd.read_sql_query(f"SELECT * FROM {table}", abce_model.db)
            final_db.to_excel(writer, sheet_name=f"{table}", engine="openpyxl")


# Run the model
if __name__ == '__main__':
    run_model()
