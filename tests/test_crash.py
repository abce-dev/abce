# Simple tests for ABCE
#
# These tests require the following file tree structure:
#
#. abce/
#├── run.py
#├── settings.yml
#├── Manifest.toml
#├── Project.toml
#└── tests/
#    ├── test_crash.py <--- THIS FILE
#    ├── Manifest.toml (optional: will be copied from abce/ if not found)
#    ├── Project.toml  (optional: will be copied from abce/ if not found)


import os
import shutil
import subprocess as sp
import pytest
import logging

@pytest.fixture
def current_dir():
    return os.getcwd()

def test_crash(current_dir):
    # Set up some file paths.
    # Check current dir and one level above for run.py; if found, use that
    #   dir as the ABCE toplevel dir
    abce_dir = None
    if os.path.exists(os.path.join(current_dir, "run.py")):
        abce_dir = current_dir
    elif os.path.exists(os.path.join(current_dir, "..", "run.py")):
        abce_dir = os.path.join(current_dir, "..")

    # Try to set up the relevant run-script and settings-file paths; if not
    #   found, alert the user
    try:
        run_script = os.path.join(abce_dir, "run.py")
        settings_file = os.path.join(current_dir, "settings.yml")
    except(TypeError):
        logging.error("Could not find the ABCE home directory.")
        logging.error("Be sure the testing script is in (or one directory")
        logging.error("below) the abce toplevel directory.")

    # Make sure environment files with dependencies are copied into the
    #   tests/ dir
    env_files = ["Manifest.toml", "Project.toml"]
    for env_file in env_files:
        if not os.path.exists(os.path.join(current_dir, env_file)):
            orig_toml = os.path.join(abce_dir, env_file)
            new_toml = os.path.join(current_dir, env_file)
            shutil.copy2(orig_toml, new_toml)

    # Run ABCE, and capture the process's return code
    return_code = sp.check_call(f"python3 {run_script} -f --settings_file='{settings_file}'",
                         shell=True)

    # Test passes if the return code is 0
    assert return_code == 0
