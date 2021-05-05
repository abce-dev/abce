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
    def __init__(self, settings, args):
        # Get input file locations from the settings dictionary
        unit_specs_file = settings["unit_specs_file"]
        fuel_data_file = settings["fuel_data_file"]
        demand_data_file = settings["demand_data_file"]
        price_curve_data_file = settings["price_curve_data_file"]
        time_series_data_file = settings["time_series_data_file"]
        db_file = settings["db_file"]
        # Get parameters from the settings dictionary
        self.num_agents = settings["num_agents"]
        self.first_agent_id = settings["first_agent_id"]
        self.first_asset_id = settings["first_asset_id"]
        self.total_forecast_horizon = settings["total_forecast_horizon"]
        try:
            self.use_precomputed_price_curve = settings["use_precomputed_price_curve"]
        except:
            # Default to `True` but warn the user
            print("No setting specified for price curve generation; will use precomputed curve.")
            self.use_precomputed_price_curve = True


        self.current_step = -1

        # Initialize database for managing asset and WIP construction project data
        self.db_file = db_file
        self.db, self.cur = sc.create_database(self.db_file, args.replace)


        # Load unit type specifications and fuel costs
        self.unit_types, self.unit_specs = self.load_unit_specs(unit_specs_file)
        self.fuel_costs = pd.read_csv(fuel_data_file)
        self.add_units_to_db(self.db, self.cur, self.unit_specs, self.fuel_costs)

        # Load all-period demand data into the database
        self.load_demand_data_to_db(settings)

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Create agents
        for i in range(self.first_agent_id, self.first_agent_id + self.num_agents):
            gc = GenCo(i, self, settings)
            self.schedule.add(gc)

        # Check whether a market price subsidy is in effect, and its value
        self.set_market_subsidy(settings)
        # Set up the price curve according to specifications in settings
        if self.use_precomputed_price_curve:
            price_duration_data = pc.load_time_series_data(price_curve_data_file, file_type="price", subsidy=self.subsidy_amount)
        else:
            # Create the systemwide merit order curve
            merit_curve = pc.create_merit_curve(self.db, self.current_step)
            pc.plot_curve(merit_curve, plot_name="merit_curve.png")
            # Load five-minute demand data from file
            demand_data = pc.load_time_series_data(time_series_data_file, file_type="load", peak_demand=settings["peak_demand"])
            pc.plot_curve(demand_data, plot_name="demand_curve.png")
            # Create the final price duration curve
            price_duration_data = pc.compute_price_duration_curve(demand_data, merit_curve, settings["price_cap"])
            # Save a plot of the price duration curve
            pc.plot_curve(price_duration_data, plot_name="price_duration.png")


        # Save price duration data to the database
        for i in range(len(price_duration_data)):
            self.cur.execute(f"INSERT INTO price_curve VALUES ({price_duration_data[i]})")
        self.db.commit()


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


    def load_demand_data_to_db(self, settings):
        # Load all-period demand data into the database
        demand_df = pd.read_csv(settings["demand_data_file"]) * settings["peak_demand"]
        # Create an expanded range of periods to backfill with demand_df data
        new_index = list(range(self.total_forecast_horizon))
        demand_df = demand_df.reindex(new_index, method="ffill")
        # Save data to DB
        demand_df.to_sql("demand", self.db, if_exists="replace", index_label="period")


    def set_market_subsidy(self, settings):
        self.subsidy_enabled = False
        self.subsidy_amount = 0

        if "enable_subsidy" in settings:
            self.subsidy_enabled = settings["enable_subsidy"]
            if self.subsidy_enabled and "subsidy_amount" in settings:
                self.subsidy_amount = settings["subsidy_amount"]


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
        print(pd.read_sql("SELECT * FROM assets", self.db))
        print("Table of construction project updates:")
        print(pd.read_sql("SELECT * FROM WIP_projects", self.db).tail(n=8))




    def reveal_decisions(self):
        self.cur.execute("SELECT asset_id FROM assets WHERE revealed = 'false'")
        projects_to_reveal = [i[0] for i in list(self.cur.fetchall())]
        for asset_id in projects_to_reveal:
            self.cur.execute(f"""UPDATE assets SET revealed = 'true' WHERE asset_id = '{asset_id}'""")
        self.db.commit()

