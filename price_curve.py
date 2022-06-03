##########################################################################
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


def load_time_series_data(data_file_name, file_type, peak_demand=0, output_type="np.array"):
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
        ts_df = organize_price_data(file_name, ts_df, output_type)
    elif file_type == "load":
        if peak_demand == 0:
            print(f"Using default peak demand value of {peak_demand}; " +
                   "a value must be specified if this is incorrect.")
        ts_df = organize_load_data(ts_df, peak_demand, output_type)
    return ts_df


def organize_price_data(file_name, price_df, output_type):
    if "output_DISPATCH" in file_name:
        # Old ALEAF output file format
        orig_col_name = "LMP"
        row_freq = 7
    elif "dispatch_summary" in file_name:
        # New ALEAF output file format
        orig_col_name = "LMP_dht"
        row_freq = 7
    else:
        # Assume it's an ERCOT file
        orig_col_name = "Total electricity price"
        row_freq = 1
    # Filter and organize the price data
    lamda = pd.DataFrame({"lamda": price_df[orig_col_name].iloc[::row_freq]})
    lamda = (lamda.sort_values(by = ["lamda"])
             .reset_index().drop(labels=["index"], axis=1))
    lamda["lamda"] = lamda["lamda"].apply(lambda x: min(9001, x))
    if output_type == "np.array":
        lamda = lamda.to_numpy().transpose()[0]
    return lamda


def organize_load_data(load_df, peak_demand, output_type):
    load_duration = (load_df.filter(["LoadShape"], axis=1)
                     .rename(columns={"LoadShape": "load"}))
    load_duration = (load_duration.sort_values(by=["load"], ascending=False)
                     .reset_index().drop(labels=["index"], axis=1))
    load_duration = load_duration * peak_demand
    if output_type == "np.array":
        load_duration = load_duration.to_numpy().transpose()[0]
    return load_duration


def create_merit_curve(db, current_pd):
    system_portfolio = pd.read_sql_query(
                         "SELECT assets.asset_id, assets.unit_type, " +
                         "capacity, VOM, FC_per_MWh, heat_rate FROM assets " +
                         "INNER JOIN unit_specs ON assets.unit_type " +
                         "= unit_specs.unit_type WHERE retirement_pd > 0 " +
                         "AND completion_pd <= 0", db)
    system_portfolio["MC"] = (system_portfolio.apply(
                                lambda df: df["FC_per_MWh"] + df["VOM"],
                                axis=1))
    system_portfolio = (system_portfolio.sort_values(by=["MC"], ascending=True)
                        .reset_index().drop(labels=["index"], axis=1))

    y = np.zeros(int(sum(system_portfolio["capacity"])))

    starting_index = 0
    for i in range(len(system_portfolio)):
        for j in range(int(system_portfolio.loc[i, "capacity"])):
            y[starting_index + j] = system_portfolio.loc[i, "MC"]
        starting_index += int(system_portfolio.loc[i, "capacity"])

    return y


def compute_price_duration_curve(demand, merit_curve, price_cap):
    prices = np.zeros(len(demand))
    print("Computing PDC")
    max_supply = len(merit_curve)
    for i, dem in enumerate(np.around(demand).astype(int)):
        if dem >= max_supply:
            # Demand exceeds all available capacity; set price to admin. max
            prices[i] = price_cap
        else:
            prices[i] = merit_curve[dem]
    prices = -np.sort(-prices)
    return prices


def plot_curve(data, plot_name="price_curve.png"):
    x_vals = np.arange(0, len(data), 1)

    fig, ax = plt.subplots()
    ax.plot(x_vals, data)
    ax.set_yscale("log")
    fig.savefig(plot_name)


def compute_unit_revenue(price_duration_data, unit_VOM, unit_capacity, unit_CF, hours_per_year):
    price_mask = price_duration_data["lambda"] > unit_VOM
    active_period_prices = price_duration_data[price_mask]["lamda"]
    total_revenue = (sum(active_period_prices) * unit_capacity
                     * unit_CF * hours_per_year
                     / len(price_duration_data))
    return active_period_prices, total_revenue

























