import sys

sys.path.append("../")

from run import *

# Files and parameters
settings_file_path = "./testing_settings.yml"
all_keys = ["simulation", "scenario", "constants", "file_paths", "system", "demand", "dispatch", "agent_opt", "ALEAF"]


# Test the read_settings() function
settings = read_settings(settings_file_path)

def test_read_settings_has_all_keys():
    assert all(key in all_keys for key in settings.keys())

def test_read_settings_no_loose_keys():
    assert all(key in settings.keys() for key in all_keys)

def test_read_settings_scenario_name():
    assert settings["simulation"]["ALEAF_scenario_name"] == "test_scenario"

def test_read_settings_peak_demand():
    assert settings["scenario"]["peak_demand"] == 29000


# Test the set_up_local_paths() function
settings = set_up_local_paths(settings)

def test_set_up_local_paths():
    assert settings["ABCE_abs_path"] == Path(__file__).parent
