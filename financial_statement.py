# An object to contain past records and future projections of operational and
# financial data
import pandas as pd
import numpy as np


class FinancialStatement(object):
    """An object to contain past records and future projections of operational
       and financial data for Generators and Agents.
    """
    def __init__(self, model, generator=None, agent=None):
        # Attach the parent model to the FS object
        self.model = model
#        self.current_step = self.model.current_step
        self.column_names = ['capacity', 'capex', 'revenue', 'op_cost', 'EBITDA']
        init_zeroes = np.zeros((40, len(self.column_names)))
        self.fsdata = pd.DataFrame(data = init_zeroes, columns = self.column_names)




class AgentFS(FinancialStatement):
    """The Agent's version of the financial statement. Aggregates the 
       individual financial statements for each of the Agent's owned units,
       and enables Agent decision-making.
    """
    def __init__(self, model, agent):
        super().__init__(self, model)
        self.agent = agent


    def step(self):
        """Current behavior: overwrites Agent's FS each step with new
             aggregated data.
           Planned future behavior: keeps records of past FSs to track
             performance over time.
        """
        self.fsdata = pd.DataFrame(data = np.zeros((40, len(self.column_names))), columns = self.column_names)
#        for unit in self.agent.portfolio.keys():
#            self.agent.portfolio[unit].fs.step()
        self.aggregate_unit_FSs()
        print(self.fsdata.head(10))


    def aggregate_unit_FSs(self):
        for unit in self.agent.portfolio.keys():
            self.fsdata += self.agent.portfolio[unit].fs.fsdata




class GeneratorFS(FinancialStatement):
    """The Generator's version of the financial statement. Tracks capacity
       and the usual financial statement items.
    """
    def __init__(self, model, generator):
        super().__init__(self, model)
        self.generator = generator
        self.fsdata['capacity'].iloc[self.generator.proj_completion_date:self.generator.proj_completion_date + self.generator.asset_life_remaining] = self.generator.capacity


    def step(self):
        self.update_capacity_projections()
        self.update_capex()


    def update_capacity_projections(self):
        # Update projections out to max time horizon
        self.current_step = self.generator.model.current_step
        # Project capacity
        self.fsdata['capacity'].iloc[int(self.generator.proj_completion_date):] = self.generator.capacity


    def update_capex(self):
        if self.generator.status == 'wip':
            if self.generator.xtr_expenditures is not []:
                self.fsdata['capex'].iloc[self.current_step] = self.generator.xtr_expenditures[-1]


