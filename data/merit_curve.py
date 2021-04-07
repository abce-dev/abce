import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

units = pd.read_csv("./unitParams.csv")

# Compute the total capacity for each type of unit
units["total_capacity"] = units.loc[:, "EXUNITS"] * units.loc[:, "CAP"]
units = units.sort_values(by="VOM", ignore_index=True)

# Create an x-axis variable for plotting
x = np.arange(0, sum(units.loc[:, "total_capacity"]))

# Allocate appropriate y-points
y = list()
for j in range(len(units)):
    new_y = np.ones(units.loc[j, "total_capacity"]) * units.loc[j, "VOM"]
    y = [*y, *new_y]

# Plot the merit order curve data
fig, ax = plt.subplots()
ax.plot(x, y)
plot_name = "merit_curve.png"
fig.savefig(plot_name)
