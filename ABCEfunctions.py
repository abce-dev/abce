import sqlite3
import sys
import pandas as pd

def load_database():
    try:
        db = sqlite3.connect(sys.argv[1])
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


