import sys
import yaml
import copy
import subprocess
import pandas as pd
from pathlib import Path

sread = "settings.yml"
asread = "2a1b.yml"

sout = "settings_test.yml"
asout = "2a1b_test.yml"
autorun_params_file = "autorun_params.csv"
inputs_dir = "inputs"



# Load the settings file as a dictionary
sread = Path(".") / sread
sout = Path(".") / sout
settings = yaml.load(
    open(sread, "r"),
    Loader=yaml.FullLoader
)

print(settings)

# Load the agent specs file as a dictionary
asread = Path(inputs_dir) / asread
asout = Path(inputs_dir) / asout
agent_specs = yaml.load(
    open(asread, "r"),
    Loader=yaml.FullLoader
)
as_definitive = copy.deepcopy(agent_specs)

# Load the params file as a pandas dataframe
autorun_params_file = Path(".") / autorun_params_file
runs = pd.read_csv(autorun_params_file)

for r in range(len(runs)):
    # Iterate over the length of the current parameter's list
    sc = f"cost_of_debt_{float(runs.at[r, 'cost_of_debt'])}__" + \
         f"starting_debt_{float(runs.at[r, 'starting_debt'])}__" + \
         f"starting_RE_{float(runs.at[r, 'starting_RE'])}__" + \
         f"k_{float(runs.at[r, 'k'])}"

    print(sc)

    # Update settings file
    settings["simulation"]["scenario_name"] = sc
    settings["file_paths"]["agent_specifications_file"] = str(asout.name)
    with open(sout, "w") as outfile:
        yaml.dump(settings, outfile, default_flow_style=False)

    # Update agent params file
    # cost of debt: identical
    agent_specs[202]["cost_of_debt"] = float(runs.loc[r, "cost_of_debt"].item())
    # starting debt: multiplicative
    agent_specs[202]["starting_debt"] = float(runs.at[r, "starting_debt"]) * as_definitive[202]["starting_debt"]
    # starting RE: multiplicative
    agent_specs[202]["starting_RE"] = float(runs.at[r, "starting_RE"]) * as_definitive[202]["starting_RE"]
    # k: identical
    agent_specs[202]["k"] = float(runs.at[r, "k"])

    with open(asout, "w") as outfile:
        yaml.dump(agent_specs, outfile, default_flow_style=False)

    # Run ABCE with these parameters
    subprocess.run(["python", "run.py", f"--settings_file={sout}"], stdout=sys.stdout, stderr=subprocess.STDOUT)
 
