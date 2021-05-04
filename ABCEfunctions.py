import sqlite3
import sys
import pandas as pd

def connect_to_database(db_file):
    try:
        db = sqlite3.connect(db_file)
        cur = db.cursor()
    except:
        print(sys.exc_info()[0])
        raise
    print("hello")
    print(pd.read_sql_table("assets", db))
    return db, cur


def get_table(db, table_name):
    try:
        df = pd.read_sql_table(table_name, db)
    except:
        raise
    return df


def get_next_asset_id(db, first_asset_id):
    next_id = first_asset_id

    asset_list = pd.read_sql("SELECT asset_id FROM assets", db)
    if len(asset_list) != 0:
        print(type(asset_list.loc[0, "asset_id"]))
        next_id = max(asset_list["asset_id"]) + 1

    return next_id



