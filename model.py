from mesa import Agent, Model
from mesa.time import RandomActivation
from agent import GenCo
import yaml
import pandas as pd
import numpy as np
import subprocess

# import local modules
from ABCEfunctions import *
from seed_creator import *

class GridModel(Model):
    ''' A model with some number of GenCos. '''
    def __init__(self, n, db_file, unit_data_file, fuel_data_file, demand_data_file, first_agent_id, first_asset_id):
        # Parameters
        self.num_agents = n
        self.current_step = -1
        self.future_vis = 4  # Number of time steps agents can see into the future
        self.first_agent_id = first_agent_id  # Starting number for valid agent IDs
        self.first_asset_id = first_asset_id

        # Initialize database for managing asset and WIP construction project data
        self.db_file = db_file
        clear_db_file(self.db_file)
        self.db, self.cur = create_db_file(self.db_file)
        # Create the five DB tables (see `seed_creator.py` for table specifications)
        create_assets_table(self.cur)
        create_WIP_projects_table(self.cur)
        create_agent_params_table(self.cur)
        create_unit_specs_table(self.cur)
        create_demand_table(self.cur)
        self.db.commit()
        print(f"Database created in file '{self.db_file}'.")

        # Load unit type specifications and fuel costs
        self.unit_types, self.unit_data = self.load_unit_data(unit_data_file)
        self.fuel_costs = pd.read_csv(fuel_data_file)
        self.add_units_to_db(self.db, self.cur, self.unit_data, self.fuel_costs)

        # Load demand data into the database
        demand_df = pd.read_csv(demand_data_file)
        demand_fill = pd.DataFrame(np.ones(100 - len(demand_df)), columns = ["demand"]) * demand_df.iloc[-1]["demand"]
        demand_df = demand_df.append(demand_fill, ignore_index=True)
        for period in list(demand_df.index):
            demand = demand_df.loc[period, "demand"]
            self.cur.execute(f"INSERT INTO demand VALUES ({period}, {demand})")
            self.db.commit()

        # Define the agent schedule, using randomly-ordered agent activation
        self.schedule = RandomActivation(self)

        # Create agents
        for i in range(self.first_agent_id, self.first_agent_id + self.num_agents):
            gc = GenCo(i, self)
            self.schedule.add(gc)

        # Load the agents' starting portfolios
        initial_assets = pd.read_csv('./data/portfolios.csv', skipinitialspace=True)

        # Add all initial assets to the database
        self.add_initial_assets_to_db(initial_assets, self.cur, self.db)
        print(get_table(self.db, self.cur, "assets"))


    def load_unit_data(self, filename):
        unit_file = open(filename)
        unit_data = pd.read_csv(unit_file)
        unit_types = unit_data.index
        return unit_types, unit_data


    def add_initial_assets_to_db(self, initial_assets, cur, db):
        for i in range(len(initial_assets)):
            for j in range(initial_assets.loc[i, "num_copies"]):
                asset_id = get_next_asset_id(self.db, self.cur, self.first_asset_id)
                agent_id = initial_assets.loc[i, "agent_id"]
                unit_type = initial_assets.loc[i, "unit_type"]
                completion_pd = 0
                cancellation_pd = 9999
                retirement_pd = initial_assets.loc[i, "useful_life"]
                total_capex = 0    # Dummy value
                capital_payment = self.unit_data.loc[i, "capacity"] * self.unit_data.loc[i, "uc_x"] * 1000 / initial_assets.loc[i, "useful_life"]
                cur.execute(f"""INSERT INTO assets VALUES
                                ({asset_id}, {agent_id}, '{unit_type}', 
                                 {completion_pd}, {cancellation_pd},
                                 {retirement_pd}, {total_capex}, {capital_payment})""")
                db.commit()


    def add_units_to_db(self, db, cur, unit_data, fuel_costs):
        for i in range(len(unit_data)):
            # Get unit data from file
            unit_type = unit_data.loc[i, "unit_type"]
            fuel_type = unit_data.loc[i, "fuel_type"]
            capacity = unit_data.loc[i, "capacity"]
            uc_x = unit_data.loc[i, "uc_x"]
            d_x = unit_data.loc[i, "d_x"]
            heat_rate = unit_data.loc[i, "heat_rate"]
            VOM = unit_data.loc[i, "VOM"]
            FOM = unit_data.loc[i, "FOM"]
            unit_life = unit_data.loc[i, "unit_life"]
            CF = unit_data.loc[i, "CF"]
            # Incorporate unit fuel cost data from the fuel costs file
            fuel_cost = fuel_costs[fuel_costs.fuel_type == fuel_type].cost_per_mmbtu[0]

            # Insert the values into the unit_specs DB table
            cur.execute(f"""INSERT INTO unit_specs VALUES
                            ('{unit_type}', '{fuel_type}', {capacity}, {uc_x},
                             {d_x}, {heat_rate}, {VOM}, {FOM}, {unit_life},
                             {CF}, {fuel_cost})""")
            db.commit()


    def step(self):
        ''' Advance the model by one step. '''
        self.current_step += 1
        print("\n\n====================================================")
        print(f"Model step: {self.current_step}")
        #self.set_demand_visibility_window()
        self.schedule.step()

