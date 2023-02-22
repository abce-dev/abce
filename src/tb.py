import os
import pandas as pd
import sqlite3

db_file = "./abce_db.db"

def start():
    db = sqlite3.connect(db_file)
    print(f"connected database at {db_file}.")

    return db
