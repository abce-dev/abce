import pytest
import sys
import subprocess as sp
import sqlite3
import pandas as pd
from pathlib import Path


# Files and parameters
test_settings_file_path = Path(__file__).parent / "settings_test.yml"
run_file_path = Path(__file__).parent.parent / "run.py"

test_db_file = Path(__file__).parent / "test_db.db"
check_db_file = Path(__file__).parent / "check_db.db"

ABCE_run_cmd = [
    "python3",
    str(run_file_path),
    "-f",
    f"--settings_file={test_settings_file_path}"
]

###############################################################################
# Tests
###############################################################################

def test_crash():
    # Run ABCE and check whether the process crashes
    proc = sp.check_call(ABCE_run_cmd)
    assert proc == 0

## Set up the test run's output database as a pytest fixture
#@pytest.fixture
#def test_db():
#    return sqlite3.connect(test_db_file)

# Set up the standard check database as a pytest fixture
#@pytest.fixture
#def check_db():
#    return sqlite3.connect("./check_db.db")

#t_db = sqlite3.connect(test_db_file)
#c_db = sqlite3.connect("./check_db.db")

#tables_list = pd.read_sql_query("SELECT name FROM sqlite_schema WHERE type='table' AND name NOT LIKE 'sqlite_%'", c_db)["name"].tolist()

#@pytest.fixture(params = tables_list)
#def table_name(request):
#    return request.param

#def test_table_identical(check_db, test_db, table_name):
#    sql_stmt = f"SELECT * FROM {table_name}"
#    c_table = pd.read_sql_query(sql_stmt, check_db)
#    t_table = pd.read_sql_query(sql_stmt, test_db)
#    assert c_table.equals(t_table)

