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
    lamda = lamda.to_numpy().transpose()[0]
    print(lamda)
    return lamda


def organize_load_data(load_df, peak_demand):
    load_duration = load_df.filter(["LoadShape"], axis=1).rename(columns={"LoadShape": "load"})
    load_duration = load_duration.sort_values(by = ["load"], ascending = False).reset_index().drop(labels=["index"], axis=1)
    load_duration = load_duration * peak_demand
    return load_duration


def create_merit_curve(db, current_pd):
    system_portfolio = pd.read_sql_query(f"SELECT asset_id, unit_type FROM assets WHERE retirement_pd > {current_pd} AND completion_pd <= 0", db)
    unit_specs = pd.read_sql_query("SELECT * FROM unit_specs", db)
    system_portfolio["capacity"] = 0
    system_portfolio["VOM"] = 0
    system_portfolio["FC_per_MMBTU"] = 0
    system_portfolio["heat_rate"] = 0
    for i in range(len(system_portfolio)):
        unit_type = system_portfolio.loc[i, "unit_type"]
        system_portfolio.loc[i, "capacity"] = unit_specs.loc[unit_specs["unit_type"] == unit_type, "capacity"].values[0]
        system_portfolio.loc[i, "VOM"] = unit_specs.loc[unit_specs["unit_type"] == unit_type, "VOM"].values[0]
        system_portfolio.loc[i, "FC_per_MMBTU"] = unit_specs.loc[unit_specs["unit_type"] == unit_type, "fuel_cost"].values[0]
        system_portfolio.loc[i, "heat_rate"] = unit_specs.loc[unit_specs["unit_type"] == unit_type, "heat_rate"].values[0]
    system_portfolio["MC"] = system_portfolio.apply(lambda df: df["heat_rate"] * df["FC_per_MMBTU"]/1000 + df["VOM"], axis=1)
    system_portfolio = system_portfolio.sort_values(by = ["MC"], ascending = True).reset_index().drop(labels=["index"], axis=1)

    x = np.arange(0, sum(system_portfolio["capacity"]), 1)
    y = np.zeros(len(x))

    starting_index = 0
    for i in range(len(system_portfolio)):
        for j in range(int(system_portfolio.loc[i, "capacity"])):
            y[starting_index + j] = system_portfolio.loc[i, "MC"]
        starting_index += int(system_portfolio.loc[i, "capacity"])

    return y


def compute_price_duration_curve(demand, merit_curve):
    prices = np.zeros(len(demand))
    print("Computing PDC")
    for i in range(len(demand)):
        if int(np.around(demand.loc[i, "load"], 0)) >= len(merit_curve):
            # Demand exceeds all available capacity; set price to administrative maximum
            prices[i] = 9001
        else:
            prices[i] = merit_curve[int(np.around(demand.loc[i, "load"], 0))]
    prices = -np.sort(-prices)
    return prices


def plot_curve(data, plot_name="price_curve.png"):
    x_vals = np.arange(0, len(data), 1)

    fig, ax = plt.subplots()
    ax.plot(x_vals, data)
    #ax.set_xscale("log")
    ax.set_yscale("log")
    plt.ylim(1, 9500)
    fig.savefig(plot_name)


def compute_unit_revenue(price_duration_data, unit_VOM, unit_capacity, unit_CF, hours_per_year):
    active_period_prices = price_duration_data[price_duration_data["lamda"] > unit_VOM]["lamda"]
    total_revenue = sum(active_period_prices) * unit_capacity * unit_CF * hours_per_year / len(price_duration_data)
    return active_period_prices, total_revenue

























