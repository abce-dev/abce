import os
import argparse
import ABCEfunctions as af
import yaml


def cli_args():
    """
    Set up the command-line argument parser. Then, read and parse whatever
      arguments are provided via sys.argv into the argument types defined below.

    Returns:
      args (argparse object): populated namespace, with argument strings as
        attributes. Retrieve values with args.<argument_name>.
    """
    parser = argparse.ArgumentParser(
        description="Postprocess A-LEAF outputs for ABCE.")
    parser.add_argument(
        "--toplevel_dir",
        "-t",
        help="Top-level ABCE directory, containing settings file and output directory.",
        default=".")
    parser.add_argument(
        "--settings_file",
        help="Path to the simulation settings file name INSIDE the top-level ABCE directory.",
        default="settings.yml")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # Retrieve the location of the settings file from argparse
    args = cli_args()

    # Read in settings
    settings_file_name = os.path.join(args.toplevel_dir, args.settings_file)
    try:
        settings = yaml.load(
            open(
                settings_file_name,
                "r"),
            Loader=yaml.FullLoader)
    except BaseException:
        print(
            f"Could not load settings file {settings_file_name}. Check the filename and try again.")

    # Get the location of the output directory
    output_dir = os.path.join(
        args.toplevel_dir,
        "outputs",
        settings["ALEAF_scenario_name"])

    af.process_outputs(settings, output_dir)
