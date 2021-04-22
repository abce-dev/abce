from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
import yaml
import pandas as pd
import numpy as np
import subprocess

# import local modules
import ABCEfunctions as ABCE
import seed_creator as sc
import price_curve as pc

class GridModel(Model):
    ''' A model with some number of GenCos. '''
    def __init__(self, settings_file):
        with open(settings_file) as setfile:
            self.settings = yaml.load(setfile, Loader=yaml.FullLoader)
        # Get input file locations from the settings dictionary
        unit_specs_file = self.settings["unit_specs_file"]
        fuel_data_file = self.settings["fuel_data_file"]
        demand_data_file = self.settings["demand_data_file"]
        price_curve_data_file = self.settings["price_curve_data_file"]
        db_file = self.settings["db_file"]
        # Get parameters from the settings dictionary
        self.num_agents = self.settings["num_agents"]
        self.first_agent_id = self.settings["first_agent_id"]
        self.first_asset_id = self.settings["first_asset_id"]
        self.total_forecast_horizon = self.settings["total_forecast_horizon"]

        self.current_step = -1

        # Initialize database for managing asset and WIP construction project data
        self.db_file = db_file
        self.db, self.cur = sc.create_database(self.db_file)

        # Load unit type specifications and fuel costs
        self.unit_types, self.unit_specs = self.load_unit_specs(unit_specs_file)
        self.fuel_costs = pd.read_csv(fuel_data_file)
        self.add_units_to_db(self.db, self.cur, self.unit_specs, self.fuel_costs)

        # Load demand data into the database
        demand_df = pd.read_csv(demand_data_file)
        demand_fill = pd.DataFrame(np.ones(self.total_forecast_horizon - len(demand_df)), columns = ["demand"]) * demand_df.iloc[-1]["demand"]
        demand_df = demand_df.append(demand_fill, ignore_index=True)
        for period in list(demand_df.index):
            demand = demand_df.loc[period, "demand"]
            self.cur.execute(f"INSERT INTO demand VALUES ({period}, {demand})")
            self.db.commit()

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Create agents
        for i in range(self.first_agent_id, self.first_agent_id + self.num_agents):
            gc = GenCo(i, self, settings_file)
            self.schedule.add(gc)

        # Check whether a market price subsidy is in effect, and its value
        self.set_market_subsidy()
        # Load and organize the price duration data
        price_duration_data = pc.load_time_series_data(price_curve_data_file, file_type="price", subsidy=self.subsidy_amount)

        # Save price duration data to the database
        for i in range(len(price_duration_data)):
            price = price_duration_data.loc[i, "lamda"]
            self.cur.execute(f"INSERT INTO price_curve VALUES ({price})")
        self.db.commit()
        print(ABCE.get_table(self.db, self.cur, "assets"))


    def load_unit_specs(self, filename):
        unit_specs = pd.read_csv(open(filename))
        unit_types = unit_specs.index
        return unit_types, unit_specs


    def add_units_to_db(self, db, cur, unit_specs, fuel_costs):
        for i in range(len(unit_specs)):
            # Get unit data from file
            unit_type = unit_specs.loc[i, "unit_type"]
            fuel_type = unit_specs.loc[i, "fuel_type"]
            capacity = unit_specs.loc[i, "capacity"]
            uc_x = unit_specs.loc[i, "uc_x"]
            d_x = unit_specs.loc[i, "d_x"]
            heat_rate = unit_specs.loc[i, "heat_rate"]
            VOM = unit_specs.loc[i, "VOM"]
            FOM = unit_specs.loc[i, "FOM"]
            unit_life = unit_specs.loc[i, "unit_life"]
            CF = unit_specs.loc[i, "CF"]
            # Incorporate unit fuel cost data from the fuel costs file
            fuel_cost = fuel_costs[fuel_costs.fuel_type == fuel_type].reset_index().loc[0, "cost_per_mmbtu"]

            # Insert the values into the unit_specs DB table
            cur.execute(f"""INSERT INTO unit_specs VALUES
                            ('{unit_type}', '{fuel_type}', {capacity}, {uc_x},
                             {d_x}, {heat_rate}, {VOM}, {FOM}, {unit_life},
                             {CF}, {fuel_cost})""")
            db.commit()


    def set_market_subsidy(self):
        try:
            self.subsidy_enabled = self.settings["enable_subsidy"]
            if self.subsidy_enabled == True:
                self.subsidy_amount = self.settings["subsidy_amount"]
            else:
                self.subsidy_amount = 0
        except:
            # If "enable_subsidy" is not specified in settings.yml, set to False
            print("No subsidy enabled")
            self.subsidy_enabled = False
            self.subsidy_amount = 0


    def step(self):
        ''' Advance the model by one step. '''
        self.current_step += 1
        print("\n\n\n\n==========================================================================")
        print(f"Model step: {self.current_step}")
        print("==========================================================================")
        self.schedule.step()
        print("\nAll agent turns are complete.\n")
        # Reveal new information to all market participants
        self.reveal_decisions()
        print("Table of all assets:")
        print(ABCE.get_table(self.db, self.cur, "assets"))
        print("Table of construction project updates:")
        print(ABCE.get_table(self.db, self.cur, "WIP_projects").tail(n=8))




    def reveal_decisions(self):
        self.cur.execute("SELECT asset_id FROM assets WHERE revealed = 'false'")
        projects_to_reveal = [i[0] for i in list(self.cur.fetchall())]
        for asset_id in projects_to_reveal:
            self.cur.execute(f"""UPDATE assets SET revealed = 'true' WHERE asset_id = '{asset_id}'""")
        self.db.commit()

