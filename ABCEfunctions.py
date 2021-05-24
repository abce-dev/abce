import sqlite3
import sys
import pandas as pd


def get_next_asset_id(db, first_asset_id):
    max_id = pd.read_sql("SELECT MAX(asset_id) FROM assets", db).iloc[0, 0]
    if max_id == None:
        next_id = first_asset_id
    else:
        next_id = max_id + 1

    return next_id



