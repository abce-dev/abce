import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

try:
    price_file_name = sys.argv[1]
except:
    print("No price data file specified.")
    print(sys.exc_info()[0])
    raise

file_name, file_ext = os.path.splitext(price_file_name)

if "csv" in file_ext:
    price_df = pd.read_csv(price_file_name)
elif "xls" in file_ext:
    # ERCOT Excel file
    price_df = pd.read_excel(price_file_name, engine="openpyxl", sheet_name="Jan")
    print(price_df)

if "output" in file_name or "DISPATCH" in file_name:
    # ALEAF output file
    lamda = price_df.filter(["LMP"], axis=1).iloc[::7].reset_index().drop(labels=["index"], axis=1)
    lamda = lamda.rename(columns={"LMP": "lamda"})
    origin = "ALEAF"
    year = 2019
else:
    # Assume it's an ERCOT file
    lamda = price_df.filter(["System Lamda"], axis=1).rename(columns={"System Lamda": "lamda"})
    # ERCOT data can be negative, so ensure there are no zero or negative values
    lamda_adder = -min(lamda["lamda"]) + 0.001
    lamda += lamda_adder
    origin = "ERCOT_historical"
    year = file_name.split("_")[-1]

print(lamda)

lamda = lamda.sort_values(by = ["lamda"], ascending = False)

x_vals = np.arange(0, len(lamda), 1)

fig, ax = plt.subplots()
ax.plot(x_vals, lamda)
ax.set_xscale("log")
ax.set_yscale("log")

plot_name = f"price_curve_{origin}_{year}.png"
fig.savefig(plot_name)

lamda.to_csv("./price_duration_data.csv", index=False)
