import sqlite3
import sys
import pandas as pd


def get_next_asset_id(db, first_asset_id):
    max_id = first_asset_id

    asset_list = pd.read_sql("SELECT asset_id FROM assets", db)
    if len(asset_list) != 0:
        max_id = pd.read_sql("SELECT MAX(asset_id) FROM assets", db).iloc[0, 0] + 1

    return max_id



