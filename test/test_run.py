import pytest
import sys

sys.path.append("../")

from run import *

# Files and parameters
settings_file_path = "./settings_test.yml"
all_keys = ["simulation", "scenario", "constants", "file_paths", "system", "demand", "dispatch", "agent_opt", "ALEAF"]

# Create a settings fixture
@pytest.fixture
def settings():
    return read_settings(settings_file_path)


# Test the read_settings() function
#settings = read_settings(settings_file_path)

def test_read_settings_has_all_keys(settings):
    assert all(key in all_keys for key in settings.keys())

def test_read_settings_no_loose_keys(settings):
    assert all(key in settings.keys() for key in all_keys)

def test_read_settings_scenario_name(settings):
    assert settings["simulation"]["ALEAF_scenario_name"] == "test_scenario"

def test_read_settings_peak_demand(settings):
    assert settings["scenario"]["peak_demand"] == 29000


# Test the set_up_local_paths() function
@pytest.fixture
def fsettings(settings):
    return set_up_local_paths(settings)

def test_set_up_local_paths_ABCE_abs_path(fsettings):
    # First convert the ABCE_abs_path saved in settings["file_paths"] into an
    #   absolute Posix path
    # Then compare it to the directory one level up from this script
    assert Path(fsettings["file_paths"]["ABCE_abs_path"]).resolve() == Path(__file__).parent.parent

def test_set_up_local_paths_ALEAF_abs_path_run_true(settings):
    # Save the original run_ALEAF value to avoid contamination
    orig_run_ALEAF_value = settings["simulation"]["run_ALEAF"]

    # Set run_ALEAF value to True for testing
    settings["simulation"]["run_ALEAF"] = True
    settings = set_up_local_paths(settings)

    assert settings["ALEAF"]["ALEAF_abs_path"] == Path(os.environ["ALEAF_DIR"])

    # Reset run_ALEAF to its original value
    settings["simulation"]["run_ALEAF"] = orig_run_ALEAF_value

def test_set_up_local_paths_ALEAF_abs_path_run_false(settings):
    # Save the original run_ALEAF value to avoid contamination
    orig_run_ALEAF_value = settings["simulation"]["run_ALEAF"]

    # Set run_ALEAF value to False for testing
    settings["simulation"]["run_ALEAF"] = False
    settings = set_up_local_paths(settings)

    assert settings["ALEAF"]["ALEAF_abs_path"] == Path("NULL_PATH")

    # Reset run_ALEAF to its original value
    settings["simulation"]["run_ALEAF"] = orig_run_ALEAF_value
