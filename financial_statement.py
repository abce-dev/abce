# An object to contain past records and future projections of operational and
# financial data
import pandas as pd
import numpy as np


class FinancialStatement(object):
    """An object to contain past records and future projections of operational
       and financial data for Generators and Agents.
    """
    def __init__(self, generator):
        # Attach the owning Generator unit to the FS object
        self.generator = generator
        self.current_step = self.generator.model.current_step
        column_names = ['capacity', 'revenue', 'op_cost', 'EBITDA']
        init_zeroes = np.zeros((40,4))
        self.fsdata = pd.DataFrame(data = init_zeroes, columns = column_names)
        self.fsdata['capacity'].iloc[self.generator.proj_completion_date:] = self.generator.capacity


    def step(self):
        self.update_projections()


    def update_projections(self):
        # Update projections out to max time horizon
        self.current_step = self.generator.model.current_step
        # Project capacity
        self.fsdata['capacity'].iloc[int(self.generator.proj_completion_date):] = self.generator.capacity
        if self.generator.id == 103 or self.generator.id == 104:
            print(self.fsdata.iloc[self.current_step:self.current_step + 5])
