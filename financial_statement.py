# An object to contain past records and future projections of operational and
# financial data
import pandas as pd
import numpy as np


class FinancialStatement(object):
    """An object to contain past records and future projections of operational
       and financial data for Generators and Agents.
    """
    def __init__(self, model, agent, generator=None):
        # Attach the parent model to the FS object
        self.model = model
        self.agent = agent
        self.fs_columns = ['capacity', 'revenue', 'op_cost', 'EBITDA', 'EBIT', 'EBT', 'tax', 'net_income', 'fcf']
        init_zeroes = np.zeros((40, len(self.fs_columns)))
        self.fsdata = pd.DataFrame(data = init_zeroes, columns = self.fs_columns)
        self.dep_columns = ['capex', 'PPE', 'inc_depreciation']
        dep_zeroes = np.zeros((40, len(self.dep_columns)))
        self.depreciation_schedule = pd.DataFrame(data = dep_zeroes, columns = self.dep_columns)
        self.debt_columns = ['principal', 'interest']
        debt_zeroes = np.zeros((40,len(self.debt_columns)))
        self.debt_schedule = pd.DataFrame(data = debt_zeroes, columns = self.debt_columns)


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
        self.fsdata['tax'] = self.fsdata['EBT'] * self.agent.tax_rate

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
        self.model = model
        self.fs_history = list()
        self.debt_history = list()
        self.dep_history = list()
        self.aggregate_unit_FSs()


    def step(self):
        """Empties last period's working financial statement sheet data;
           prompts each unit to update its financial statement; and aggregates
           all data into a central master sheet. The current step's sheet is
           saved into the AgentFS.history list.
        """
        # Empty last period's data from the FS DataFrames
        self.fsdata = pd.DataFrame(data = np.zeros((40, len(self.fs_columns))), columns = self.fs_columns)
        self.debt_schedule = pd.DataFrame(data = np.zeros((40, len(self.debt_columns))), columns = self.debt_columns)
        self.depreciation_schedule = pd.DataFrame(data = np.zeros((40, len(self.dep_columns))), columns = self.dep_columns)
        # Add up each generator's financial statements to create master records
        self.aggregate_unit_FSs()
        # Display some recent data rows
#        if self.agent.current_step < 4:
#            print(self.fsdata.head())
#        else:
#            print(self.fsdata.iloc[self.agent.current_step-3:self.agent.current_step+2])
        # Save the current period's financial sheet data to the self.history list.
        self.fs_history.append(self.fsdata)
        self.dep_history.append(self.depreciation_schedule)
        self.debt_history.append(self.debt_schedule)
        print(self.debt_schedule.head())


    def aggregate_unit_FSs(self):
        for unit in self.agent.portfolio.keys():
            self.fsdata += self.agent.portfolio[unit].fs.fsdata
            self.depreciation_schedule += self.agent.portfolio[unit].fs.depreciation_schedule
            self.debt_schedule += self.agent.portfolio[unit].fs.debt_schedule





class GeneratorFS(FinancialStatement):
    """The Generator's version of the financial statement. Tracks capacity
       and the usual financial statement items.
    """
    def __init__(self, model, agent, generator):
        super().__init__(self, model)
        self.model = model
        self.generator = generator
        self.fsdata['capacity'].iloc[self.generator.proj_completion_date:self.generator.proj_completion_date + self.generator.asset_life_remaining] = self.generator.capacity


    def step(self):
        self.current_step = self.generator.model.current_step
        self.update_capacity_projections()
        self.update_capex()
        self.update_debt()
        self.update_revenue()
        self.update_opcost()
        self.update_EBITDA()
        self.update_EBIT()
        self.update_EBT()
        self.update_tax()
        self.update_net_income()
        self.update_fcf()


    def update_capacity_projections(self):
        # Project capacity
        self.fsdata['capacity'].iloc[int(self.generator.proj_completion_date):] = self.generator.capacity


    def update_capex(self):
        if self.generator.status == 'wip':
            if self.generator.xtr_expenditures is not []:
                self.depreciation_schedule['capex'].iloc[self.current_step] = self.generator.xtr_expenditures[-1]


    def update_debt(self):
        if self.generator.status == 'wip':
            self.debt_schedule['principal'].iloc[self.current_step] = self.debt_schedule['principal'].iloc[self.current_step - 1] + self.depreciation_schedule['capex'].iloc[self.current_step] * self.generator.agent.debt_fraction


    def update_depreciation(self):
        pass


    def update_tax(self):
        self.fsdata['tax'] = self.fsdata['EBT'] * self.generator.agent.tax_rate










