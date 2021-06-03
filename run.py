##########################################################################
# Copyright 2021 Argonne National Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################

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
    parser.add_argument("--quiet", "-q",
                          action="store_true",
                          help="Suppress all output except the turn and period counters.")
    parser.add_argument("--demo", "-d",
                          action="store_true",
                          help="Pause the simulation after each step until user presses a key.")
    args = parser.parse_args()
    return args


def wait_for_user(is_demo):
    """
    If the user specifies the --demo or -d command-line flag, then pause and
      prompt the user to hit 'enter' after every time-step of the simulation.
    """
    if is_demo:
        print("\n")
        user_response = input("Press Enter to continue: ")


def run_model():
    """
    Run the model:
      - process command-line arguments
      - read in settings
      - run the model
      - pull the completed DB into a pandas DataFrame and save it to xlsx
    """    
    args = cli_args()

    settings = read_settings(args.settings_file)

    abce_model = GridModel(settings, args)
    for i in range(settings["num_steps"]):
        abce_model.step()
        wait_for_user(args.demo)

    db_tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE " +
                                  "type='table';", abce_model.db)
    with pd.ExcelWriter(settings["output_file"]) as writer:
        for i in range(len(db_tables)):
            table = db_tables.loc[i, "name"]
            final_db = pd.read_sql_query(f"SELECT * FROM {table}", abce_model.db)
            final_db.to_excel(writer, sheet_name=f"{table}", engine="openpyxl")


# Run the model
if __name__ == '__main__':
    run_model()
