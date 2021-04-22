import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def get_file_name(filename = None):
    if filename is not None:
        price_file_name = filename
    else:
        try:
            price_file_name = sys.argv[1]
        except:
            print("No price data file specified.")
            print(sys.exc_info()[0])
            raise
    return price_file_name


def load_time_series_data(data_file_name, file_type, subsidy=0, peak_demand=0):
    file_name, file_ext = os.path.splitext(data_file_name)

    if "csv" in file_ext:
        ts_df = pd.read_csv(data_file_name)
    elif "xls" in file_ext and file_type == "price":
        ts_df = pd.read_excel(data_file_name, engine="openpyxl", sheet_name="Jan")
    elif "xls" in file_ext and file_type == "load":
        # timeseriesParams.xlsx only has one sheet
        ts_df = pd.read_excel(data_file_name, engine="openpyxl")
    else:
        # An unsupported file format has been provided; alert the user and end
        print("The file specified for the price curve data is not .csv or .xls/.xlsx.")
        print("Please provide a file in one of those formats.")
        print("Terminating...")
        sys.exit()

    # Invoke an appropriate organization function, depending on file type
    if file_type == "price":
        ts_df = organize_price_data(file_name, ts_df, subsidy)
    elif file_type == "load":
        if peak_demand == 0:
            print(f"Using default peak demand value of {peak_demand}; a value must be specified if this is incorrect.")
        ts_df = organize_load_data(ts_df, peak_demand)
    return ts_df


def organize_price_data(file_name, price_df, subsidy):
    if "output" in file_name or "DISPATCH" in file_name:
        # ALEAF output file
        lamda = price_df.filter(["LMP"], axis=1).iloc[::7].reset_index().drop(labels=["index"], axis=1)
        lamda = lamda.rename(columns={"LMP": "lamda"})
    else:
        # Assume it's an ERCOT file
        lamda = price_df.filter(["Total electricity price"], axis=1).rename(columns={"Total electricity price": "lamda"})
        lamda = lamda.reset_index().drop(labels=["index"], axis=1)
    lamda = lamda.sort_values(by = ["lamda"], ascending = False).reset_index().drop(labels=["index"], axis=1)
    lamda["lamda"] = lamda["lamda"].apply(lambda x: min(9001, x + subsidy))
    return lamda


def organize_load_data(load_df, peak_demand):
    load_duration = load_df.filter(["LoadShape"], axis=1).rename(columns={"LoadShape": "load"})
    load_duration = load_duration.sort_values(by = ["load"], ascending = False).reset_index().drop(labels=["index"], axis=1)
    load_duration = load_duration * peak_demand
    return load_duration


def create_dispatch_curve(active_assets_ids, db):
    cur = db.cursor()
    for asset_id in active_assets_ids:
        cur.execute(f"SELECT asset_id, unit_type, capacity, VOM, fuel_cost FROM assets WHERE assed_id = {asset_id}")


def plot_price_duration_curve(lamda, origin, year, plot_name="price_curve.png"):
    x_vals = np.arange(0, len(lamda), 1)

    fig, ax = plt.subplots()
    ax.plot(x_vals, lamda)
    ax.set_xscale("log")
    ax.set_yscale("log")
    fig.savefig(plot_name)


def compute_unit_revenue(price_duration_data, unit_VOM, unit_capacity, unit_CF, hours_per_year):
    active_period_prices = price_duration_data[price_duration_data["lamda"] > unit_VOM]["lamda"]
    total_revenue = sum(active_period_prices) * unit_capacity * unit_CF * hours_per_year / len(price_duration_data)
    return active_period_prices, total_revenue

























