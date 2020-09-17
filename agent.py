from mesa import Agent, Model
from mesa.time import RandomActivation
import yaml

# import local modules
import generator as gen

class GenCo(Agent):
    ''' A utility company with a certain number of generation assets. '''
    def __init__(self, genco_id, model, existing_portfolio):
        super().__init__(genco_id, model)
        self.set_up_portfolio(existing_portfolio)
        self.model = model

    def set_up_portfolio(self, existing_portfolio):
        self.portfolio = dict()
        for unit_num in existing_portfolio.keys():
            gtype = existing_portfolio[unit_num]['gtype']
            unit_id = self.model.id_register.get_next_available_id()
            new_unit = gen.Generator(id_num=unit_id, gtype=gtype)
            new_unit.completion = [1]
            self.portfolio[unit_id] = new_unit
            self.model.id_register.add_unit(self.unique_id, unit_id)
            

    def step(self):
        self.set_current_step()
        for id_num in self.portfolio:
            self.portfolio[id_num].step()
        num_units = len(self.portfolio.keys())
        print(f'I\'m GenCo # {self.unique_id}, and I have {num_units} units.')
        self.forecast_demand()
        self.evaluate_current_capacity()
        self.evaluate_new_project()

    def evaluate_current_capacity(self):
        self.total_capacity = 0
        for unit in self.portfolio.keys():
            self.total_capacity += self.portfolio[unit].capacity
        print(f'Total capacity currently installed = {self.total_capacity}.')

    def evaluate_new_project(self):
        for pd in self.demand_forecast:
            if pd > self.total_capacity:
                # Additional capacity needed!
                print("Demand will be higher than capacity. Building a new unit...")
                unit_id = self.model.id_register.get_next_available_id()
                new_unit = gen.Generator(id_num=unit_id, gtype='unit_1')
                self.portfolio[unit_id] = new_unit
                return

    def set_current_step(self):
        self.current_step = self.model.current_step
        print(f'\n\n Current step = {self.current_step}')

    def get_demand_forecast(self):
        if self.model.demand_NTF != []:
            self.demand_forecast = self.model.demand_NTF
            self.current_demand = self.demand_forecast[0]

    def forecast_demand(self):
        # Directly pull demand visibility window from the model
        # Model has complete visibility into a fixed sequence of demand levels
        #    from an input file; only provides next 3 periods to agents
        if self.model.demand_NTF != []:
            self.demand_forecast = self.model.demand_NTF
            self.current_demand = self.demand_forecast[0]

        # Currently can only detect a single upcoming step-change in demand
        for i in range(len(self.demand_forecast)):
            if self.demand_forecast[i] == self.current_demand:
                continue
            if self.demand_forecast[i] > self.current_demand:
                print(f'Demand will increase in {i} periods, from {self.current_demand} to {self.demand_forecast[i]}.')
                self.next_demand_increase_period = self.current_step + i
                return
            elif self.demand_forecast[i] < self.current_demand:
                print(f'Demand will decrease in {i} periods, from {self.current_demand} to {self.demand_forecast[i]}.')
                self.next_demand_increase_period = self.current_step + i
                return

