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
        fs_columns = ['capacity', 'capex', 'revenue', 'op_cost', 'EBITDA', 'EBIT', 'EBT', 'tax', 'net_income', 'fcf']
        init_zeroes = np.zeros((40, len(fs_columns)))
        self.fsdata = pd.DataFrame(data = init_zeroes, coluns = self.column_names)
        dep_columns = ['capex', 'PPE', 'inc_depreciation']
        dep_zeroes = np.zeros((40, len(dep_columns)))
        self.depreciation_schedule = pd.DataFrame(data = dep_zeroes, columns = dep_columns)
        debt_columns = ['principal', 'interest']
        debt_zeroes = np.zeros((40,len(self.dep_zeroes)))
        self.debt_schedule = pd.DataFrame(data = debt_zeroes, columns = debt_columns)


    def update_revenue(self):
        """Using the past period's weighted-average electricity price, compute
           the revenue earned by this generator.
        """
        current_price = self.model.eprices[self.current_step]
        self.fsdata['revenue'].iloc[self.current_step] = current_price * self.fsdata['capacity'].iloc[self.current_step]


    def update_opcost(self):
        self.fsdata['op_cost'].iloc[self.current_step] = self.generator.variable_cost + self.generator.fixed_cost

    def update_EBITDA(self):
        self.fsdata['EBITDA'] = self.fsdata['revenue'] - self.fsdata['op_cost']

    def update_EBIT(self):
        self.fsdata['EBIT'] = self.fsdata['EBITDA'] - self.depreciation_schedule['inc_depreciation']

    def update_EBT(self):
        self.fsdata['EBT'] = self.fsdata['EBIT'] - self.debt_schedule['interest']

    def update_tax(self):
        self.fsdata['tax'] = self.fsdata['EBT'] * self.generator.agent.tax_rate

    def update_net_income(self):
        self.fsdata['net_income'] = self.fsdata['EBT'] - self.fsdata['tax']

    def update_fcf(self):
        self.fsdata['fcf'] = self.fsdata['net_income'] + self.depreciation_schedule['inc_depreciation']



class AgentFS(FinancialStatement):
    """The Agent's version of the financial statement. Aggregates the 
       individual financial statements for each of the Agent's owned units,
       and enables Agent decision-making.
    """
    def __init__(self, model, agent):
        super().__init__(self, model)
        self.agent = agent
        self.aggregate_unit_FSs()


    def step(self):
        """Current behavior: overwrites Agent's FS each step with new
             aggregated data.
           Planned future behavior: keeps records of past FSs to track
             performance over time.
        """
        self.fsdata = pd.DataFrame(data = np.zeros((40, len(self.column_names))), columns = self.column_names)
        self.aggregate_unit_FSs()
        if self.agent.current_step < 4:
            print(self.fsdata.head())
        else:
            print(self.fsdata.iloc[self.agent.current_step-3:self.agent.current_step+2])


    def aggregate_unit_FSs(self):
        for unit in self.agent.portfolio.keys():
            self.fsdata += self.agent.portfolio[unit].fs.fsdata


    def update_capacity_projections(self):
        # Project capacity
        for unit in self.agent.portfolio.keys():
            current_unit = self.agent.portfolio[unit]
            self.fsdata['capacity'].iloc[int(current_unit.proj_completion_date):] = current_unit.capacity


    def update_capex(self):
        for unit in self.agent.portfolio.keys():
            current_unit = self.agent.portfolio.keys()
            if current_unit.status == 'wip':
                if current_unit.xtr_expenditures is not []:
                    self.fsdata['capex'].iloc[self.current_step] = current_unit.xtr_expenditures[-1]


    def update_tax(self):
        self.fsdata['tax'] = self.fsdata['EBT'] * self.agent.tax_rate








class GeneratorFS(FinancialStatement):
    """The Generator's version of the financial statement. Tracks capacity
       and the usual financial statement items.
    """
    def __init__(self, model, generator):
        super().__init__(self, model)
        self.model = model
        self.generator = generator
        self.fsdata['capacity'].iloc[self.generator.proj_completion_date:self.generator.proj_completion_date + self.generator.asset_life_remaining] = self.generator.capacity


    def step(self):
        self.current_step = self.generator.model.current_step
        self.update_capacity_projections()
        self.update_capex()
        self.update_revenue()
        self.update_opcost()
        self.update_EBITDA()
        self.update_EBIT()
        self.update_EBT()
        self.update_tax()
        self.update_net_income()


    def update_capacity_projections(self):
        # Project capacity
        self.fsdata['capacity'].iloc[int(self.generator.proj_completion_date):] = self.generator.capacity


    def update_capex(self):
        if self.generator.status == 'wip':
            if self.generator.xtr_expenditures is not []:
                self.fsdata['capex'].iloc[self.current_step] = self.generator.xtr_expenditures[-1]


    def update_tax(self):
        self.fsdata['tax'] = self.fsdata['EBT'] * self.generator.agent.tax_rate










