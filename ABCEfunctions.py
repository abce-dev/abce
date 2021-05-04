import sqlite3
import sys
import pandas as pd


def get_next_asset_id(db, first_asset_id):
    next_id = first_asset_id

    asset_list = pd.read_sql("SELECT asset_id FROM assets", db)
    if len(asset_list) != 0:
        next_id = max(asset_list["asset_id"]) + 1

    return next_id



