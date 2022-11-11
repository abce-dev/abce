# conftest.py
from pathlib import Path


def pytest_addoption(parser):
    parser.addoption("--settings_file", action="store", default=str(Path(__file__).parent / "settings_test.yml"))
    parser.addoption("--test_db_file", action="store", default=str(Path(__file__).parent / "test_db.db"))
    parser.addoption("--check_db_file", action="store", default=str(Path(__file__).parent / "check_db.db"))
