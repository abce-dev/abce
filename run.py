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

import os
import subprocess
from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
from model import GridModel
import yaml
import pandas as pd
import argparse
import ABCEfunctions
from pathlib import Path


def read_settings(settings_file):
    """
    Read in settings from the settings file.
    """
    with open(settings_file, "r") as setfile:
        settings = yaml.load(setfile, Loader=yaml.FullLoader)
    return settings


def set_up_local_paths(args, settings):
    # Set the path for ABCE files to the directory where run.py is saved
    # settings["ABCE_abs_path"] = os.path.realpath(os.path.dirname(__file__))
    settings["ABCE_abs_path"] = Path(__file__).parent
    if settings["run_ALEAF"]:
    # Try to locate an environment variable to specify where A-LEAF is located
        try:
            settings["ALEAF_abs_path"] = Path(os.environ["ALEAF_DIR"])
        except KeyError:
            print("The environment variable ALEAF_abs_path does not appear to be set. Please make sure it points to the correct directory.")
            raise
    else:
        settings["ALEAF_abs_path"] = Path("NULL_PATH")

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
                        type=str,
                        help="Simulation settings file name.",
                        default=Path(Path.cwd()) / "settings.yml")
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress all output except the turn and period counters.")
    parser.add_argument(
        "--demo",
        "-d",
        action="store_true",
        help="Pause the simulation after each step until user presses a key.")
    args = parser.parse_args()
    return args


def check_julia_environment(ABCE_abs_path):
    """
    Check whether Manifest.toml and Project.toml files exist (necessary for
      Julia to correctly load all dependencies).
    If either one is not found, run `make_julia_environment.jl` to
      automatically generate valid .toml files.
    """
    if not ((Path(ABCE_abs_path) / "Manifest.toml").exists()
            and (Path(ABCE_abs_path) / "Project.toml").exists()):
        julia_cmd = (
            f"julia {Path(ABCE_abs_path) / 'make_julia_environment.jl'}")
        try:
            sp = subprocess.check_call([julia_cmd], shell=True)
            # print("Julia environment successfully created.\n\n")
        except subprocess.CalledProcessError:
            # print("Cannot proceed without a valid Julia environment. Terminating...")
            quit()


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

    settings = set_up_local_paths(args, settings)

    check_julia_environment(settings["ABCE_abs_path"])

    abce_model = GridModel(args.settings_file, settings, args)
    for i in range(settings["num_steps"]):
        abce_model.step(demo=args.demo)

    db_tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE " +
                                  "type='table';", abce_model.db)
    with pd.ExcelWriter(settings["output_file"]) as writer:
        for i in range(len(db_tables)):
            table = db_tables.loc[i, "name"]
            final_db = pd.read_sql_query(
                f"SELECT * FROM {table}", abce_model.db)
            final_db.to_excel(writer, sheet_name=f"{table}", engine="openpyxl")

    if abce_model.settings["run_ALEAF"]:
        # Postprocess A-LEAF results
        ABCEfunctions.process_outputs(
            settings,
            abce_model.ABCE_output_data_path,
            abce_model.unit_specs)


# Run the model
if __name__ == '__main__':
    run_model()
