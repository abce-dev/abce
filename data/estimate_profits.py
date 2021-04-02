import pandas as pd
import sys
import subprocess as sp

unit_VOM = 35    # $/MWh
unit_FOM = 7     # $/MWh
unit_cap = 100   # MW
unit_CF  = .92

price_duration_data = pd.read_csv("./price_duration_data.csv")
num_periods = len(price_duration_data)

active_period_prices = price_duration_data[price_duration_data["lamda"] > unit_VOM]["lamda"]
total_revenue = sum(active_period_prices) * unit_cap * unit_CF
print("Total revenue = ${:,.2f}".format(total_revenue))

total_VOM = len(active_period_prices) * unit_VOM * unit_cap
total_FOM = unit_FOM * 8760 * unit_cap

print("Total VOM = ${:,.2f}".format(total_VOM))
print("Total VOM = ${:,.2f}".format(total_FOM))

operating_income = total_revenue - total_FOM - total_VOM

print("Operating income = ${:,.2f}".format(operating_income))





