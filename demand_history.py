# A class to create objects which record public information about historical
# electricity demand

import pandas as pd

class DemandHistory(object):
    def __init__(self, init_data):
        if not isinstance(init_data, list):
            init_data = [init_data]
        period = list(range(len(init_data)))
        self.register = pd.DataFrame(data={'period': period, 'demand': init_data})
        self.register = self.register.set_index('period')


    def add_data(self, new_demand_point):
        #new_period = self.register.index[-1] + 1
        #print(new_period)
        self.register = self.register.append({'demand': new_demand_point},
                                             ignore_index=True)
        print(self.register)
