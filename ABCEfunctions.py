import sqlite3
import sys
import pandas as pd

def load_database(db_file):
    try:
        db = sqlite3.connect(db_file)
        cur = db.cursor()
        return db, cur
    except:
        print("Could not load database:")
        print(sys.exc_info()[0])
        raise
        quit()

def show_table(db, cur, table_name):
    try:
        names = list(map(lambda x: x[0], db.execute(f"SELECT * FROM {table_name}").description))
        row_list = list()
        for row in cur.execute(f"SELECT * FROM {table_name}"):
            row_list.append(row)
        df = pd.DataFrame(row_list, columns = names)
        print(df)
        return df
    except:
        print(f"Could not load table {table_name}:")
        print(sys.exc_info()[0])
        raise
        quit()

def get_next_asset_id(db, cur):
    cur.execute("SELECT asset_id FROM assets")
    asset_list = list(cur.fetchall())
    if len(asset_list) != 0:
        asset_list = [int(id_num[0]) for id_num in asset_list]
        next_id = max(asset_list) + 1
    else:
        # There are no existing assets
        print("There are no existing assets")
        next_id = 2001
    return next_id



