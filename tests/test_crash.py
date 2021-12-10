# Simple tests for ABCE
#
# These tests require the following file tree structure:
#
#. abce
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

@pytest.fixture
def current_dir():
    return os.getcwd()

def test_crash(current_dir):
    # Set up some file paths.
    abce_dir = os.path.join(current_dir, "..")
    run_script = os.path.join(abce_dir, "run.py")
    settings_file = os.path.join(abce_dir, "settings.yml")

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
