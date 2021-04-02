import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def get_file_name():
    try:
        price_file_name = sys.argv[1]
    except:
        print("No price data file specified.")
        print(sys.exc_info()[0])
        raise
    return price_file_name


def load_original_data(price_file_name):
    file_name, file_ext = os.path.splitext(price_file_name)

    if "csv" in file_ext:
        price_df = pd.read_csv(price_file_name)
    elif "xls" in file_ext:
        price_df = pd.read_excel(price_file_name, engine="openpyxl", sheet_name="Jan")
    return price_df


def organize_price_data(file_name, price_df):
    if "output" in file_name or "DISPATCH" in file_name:
        # ALEAF output file
        lamda = price_df.filter(["LMP"], axis=1).iloc[::7].reset_index().drop(labels=["index"], axis=1)
        lamda = lamda.rename(columns={"LMP": "lamda"})
        origin = "ALEAF"
        year = 2019
    else:
        # Assume it's an ERCOT file
        lamda = price_df.filter(["System Lamda"], axis=1).rename(columns={"System Lamda": "lamda"})
        origin = "ERCOT_historical"
        year = file_name.split("_")[-1]
    lamda = lamda.sort_values(by = ["lamda"], ascending = False)
    return lamda, origin, year


def plot_price_duration_curve(lamda, origin, year):
    x_vals = np.arange(0, len(lamda), 1)

    fig, ax = plt.subplots()
    ax.plot(x_vals, lamda)
    ax.set_xscale("log")
    ax.set_yscale("log")

    plot_name = f"price_curve_{origin}_{year}.png"
    fig.savefig(plot_name)


def write_price_data_to_file(lamda):
    lamda.to_csv("./price_duration_data.csv", index=False)


def compute_unit_revenue(price_duration_data, unit_VOM, unit_capacity, unit_CF):
    active_period_prices = price_duration_data[price_duration_data["lamda"] > unit_VOM]["lamda"]
    total_revenue = sum(active_period_prices) * unit_capacity * unit_CF * 8760 / len(price_duration_data)
    return active_period_prices, total_revenue

























