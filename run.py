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
import sys
import subprocess
import logging
from mesa import Agent, Model
from mesa.time import RandomActivation
from src.agent import GenCo
from src.model import GridModel
import yaml
import pandas as pd
import argparse
from pathlib import Path


def read_settings(settings_file):
    """
    Read in settings from the settings file.
    """
    with open(settings_file, "r") as setfile:
        settings = yaml.load(setfile, Loader=yaml.FullLoader)

    return settings


def set_up_local_paths(settings):
    # Set the path for ABCE files to the directory where run.py is saved
    settings["file_paths"]["ABCE_abs_path"] = Path(__file__).parent

    if settings["simulation"]["annual_dispatch_engine"] == "ALEAF":
    # Try to locate an environment variable to specify where A-LEAF is located
        try:
            settings["ALEAF"]["ALEAF_abs_path"] = Path(os.environ["ALEAF_DIR"])
        except KeyError:
            msg = ("The environment variable ALEAF_abs_path does not appear " +
                   "to be set. Please make sure it points to the correct " +
                   "directory.")
            logging.error(msg)
            raise
    else:
        settings["ALEAF"]["ALEAF_abs_path"] = Path("NULL_PATH")

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
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Agree to overwrite any existing DB files."
    )
    parser.add_argument(
        "--settings_file",
        type=str,
        help="Absolute path to simulation settings file.",
        default=Path(Path.cwd() / "settings.yml")
    )
    parser.add_argument(
        "--verbosity",
        choices=[0, 1, 2, 3],
        type=int,
        help="Verbosity of output during runtime. 0 = totally silent; 1 = minimal output; 2 = default output; 3 = full/debug output",
        default=2
    )
    parser.add_argument(
        "--demo",
        "-d",
        action="store_true",
        help="Pause the simulation after each step until user presses a key."
    )
    args = parser.parse_args()
    return args


def initialize_logging(args, vis_lvl):
    # Python logging levels:
    #   CRITICAL = 50
    #   ERROR =    40
    #   WARNING =  30
    #   INFO =     20
    #   DEBUG =    10
    #   NOTSET =    0

    # Default verbosity (CL setting 2) = INFO
    lvl = 20
    if args.verbosity == 0:
        # Do not show logging messages (no messages are generated with
        #   a level of 60 or above)
        lvl = 60
    elif args.verbosity == 1:
        # Only show logging messages of severity ERROR (level = 40) and above
        lvl = 40
    elif args.verbosity == 3:
        # Show all logging messages (level 0 and greater)
        lvl = 0

    fmt = ABCEFormatter(vis_lvl)
    hdlr = logging.StreamHandler(sys.stdout)

    hdlr.setFormatter(fmt)
    logging.root.addHandler(hdlr)

    logging.root.setLevel(lvl)


def check_julia_environment(ABCE_abs_path):
    """
    Check whether Manifest.toml and Project.toml files exist (necessary for
      Julia to correctly load all dependencies).
    If either one is not found, run `make_julia_environment.jl` to
      automatically generate valid .toml files.
    """
    if not ((Path(ABCE_abs_path) / "env" / "Manifest.toml").exists()
            and (Path(ABCE_abs_path) / "env" / "Project.toml").exists()):
        julia_cmd = (
            f"julia {Path(ABCE_abs_path) / 'env' / 'make_julia_environment.jl'}")
        try:
            sp = subprocess.check_call([julia_cmd], shell=True)
            logging.info("Julia environment successfully created.\n\n")
        except subprocess.CalledProcessError:
            logging.error("Cannot proceed without a valid Julia environment. Terminating...")
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

    initialize_logging(args, settings["constants"]["vis_lvl"])

    settings = set_up_local_paths(settings)

    check_julia_environment(settings["file_paths"]["ABCE_abs_path"])

    # Run the simulation
    abce_model = GridModel(settings, args)

    for i in range(settings["simulation"]["num_steps"]):
        abce_model.step(demo=args.demo)

    # Write the raw database to xlsx
    db_tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE " +
                                  "type='table';", abce_model.db)
    with pd.ExcelWriter(
             Path(
                 settings["file_paths"]["ABCE_abs_path"] /
                 "outputs" /
                 settings["simulation"]["scenario_name"] /
                 settings["file_paths"]["output_file"]
             )
        ) as writer:
        for i in range(len(db_tables)):
            table = db_tables.loc[i, "name"]
            final_db = pd.read_sql_query(
                f"SELECT * FROM {table}", abce_model.db)
            final_db.to_excel(writer, sheet_name=f"{table}", engine="openpyxl")


class ABCEFormatter(logging.Formatter):
    """ A custom log formatter for non-standard logging levels.

        This logger will handle graphical elements which should
          be printed to the console on verbosity levels 1, 2, and 3,
          but which shouldn't have any logger-style prefixes.

        This replicates the behavior of print statements, with the
          verbosity-aware handling of logging.
    """

    vis_fmt = "%(msg)s"

    def __init__(self, vis_lvl):
        super().__init__(fmt="%(levelname)s: %(msg)s", datefmt=None, style="%")
        self.vis_lvl = vis_lvl

    def format(self, record):
        # Save the original user-settingsured formatter settings for
        #   later retrieval
        format_orig = self._style._fmt

        # For records with a level of vis_lvl (visual element, specified in
        #   the settings.yml file), use custom format
        if record.levelno == self.vis_lvl:
            self._style._fmt = ABCEFormatter.vis_fmt

        # Call the original Formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format settingsured by the user
        self._style._fmt = format_orig

        return result


# Run the model
if __name__ == '__main__':
    run_model()
