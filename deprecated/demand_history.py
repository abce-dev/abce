# A class to create objects which record public information about historical
# electricity demand

import pandas as pd

class DemandHistory(object):
    """A publicly-available history of current and past demand quantities.
    """
    def __init__(self, init_data):
        """Create an instance of the DemandHistory class.

           Parameters
           ----------
           init_data : list of floats or ints
               Prior history of demand which agents can know.

           Returns
           -------
           new DemandHistory object
        """
        if not isinstance(init_data, list):
            init_data = [init_data]
        period = list(range(len(init_data)))
        self.register = pd.DataFrame(data={'period': period, 'demand': init_data})
        self.register = self.register.set_index('period')


    def add_data(self, new_demand_point):
        """Add a new data point to the time series.

           Parameters
           ----------
           new_demand_point: int or float
               A new data point to append to the end of the current data
               series; most often used to transfer the upcoming "future"
               demand to the most recent past period.

           Returns
           -------
           None
        """
        self.register = self.register.append({'demand': new_demand_point},
                                             ignore_index=True)
