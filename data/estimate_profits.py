import pandas as pd
import sys
import subprocess as sp

from price_curve import *

unit_VOM = 35         # $/MWh
unit_FOM = 10100     # $/MW-yr
unit_capacity = 100   # MW
unit_CF  = .92

# Load the price data from the file specified in the command line
price_data_file = get_file_name()
price_duration_data = load_original_data(price_data_file)

# Organize the price data
price_duration_data = organize_price_data(price_data_file, price_duration_data)
num_periods = len(price_duration_data)

print(price_duration_data)

# Plot the price duration curve
# plot_price_duration_curve(lamda, origin, year)

write_price_data_to_file(price_duration_data)

# Compute total revenue
active_period_prices, total_revenue = compute_unit_revenue(price_duration_data, unit_VOM, unit_capacity, unit_CF)
print("Total revenue = ${:,.2f}".format(total_revenue))

total_VOM = len(active_period_prices) * unit_VOM * unit_capacity * 8760 / len(price_duration_data)
total_FOM = unit_FOM * unit_capacity

print("Total VOM = ${:,.2f}".format(total_VOM))
print("Total FOM = ${:,.2f}".format(total_FOM))

operating_income = total_revenue - total_FOM - total_VOM

print("Operating income = ${:,.2f}".format(operating_income))





