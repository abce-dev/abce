# A class for generation units
import yaml
import mesa

class Generator(object):
    def __init__(self, id_num, gtype, completion=0):
        # Universal identification number, public knowledge
        self.id = id_num
        # Set up unit type, and get parameters from the units.yml file
        self.type = gtype
        units = yaml.load(open('./simplified_units.yml', 'r'), Loader=yaml.FullLoader)
        unit = units[self.type]
        self.capacity = unit['capacity']
        self.overnight_cost = unit['overnight_cost']
        self.xtr_lead_time = unit['xtr_lead_time']
        self.variable_cost = unit['variable_cost']
        self.fixed_cost = unit['fixed_cost']
        self.asset_life_remaining = unit['useful_life']

        # Set up parameters for WIP projects, if needed
        self.completion = list([completion])
        if self.completion[-1] == 1:
            self.status = 'in_service'
        else:
            self.status = 'wip'
        self.xtr_expenditures = list()
        self.total_overnight_cost = 10.0

    def step(self):
        if self.status == 'in_service':
            self.decrement_remaining_life()
        if self.status == 'wip':
            self.update_xtr_progress()
            self.update_xtr_expenses()

    def update_xtr_progress(self):
        """Update the total % completion to date of the project.

            Detailed Description
            --------------------
            Increment the unit's construction completion progress for the most
            recent period. If this additional progress would put completion
            at or over 100%, set the project to completed ('in_service' status)

            The unit adds incremental progress linearly, at a rate of
            (1 / self.xtr_lead_time) per time period.

            Parameters
            ----------
            none

            Returns
            -------
            None
        """
        new_completion = self.completion[-1] + (1/self.xtr_lead_time)
        self.completion.append(new_completion)
        # Cleanup, and status update when project finishes
        if self.completion[-1] >= 1:
            self.completion[-1] = 1
            self.status = 'in_service'
            print(f'Project {self.id} completed and entering service.')
        else:
            print(f'Project {self.id} at {self.completion[-1]} completion')

    def update_xtr_expenses(self):
        """Incur construction expenses related to most recent project progress.

            Current behavior:
             - Incremental cost is a linear function of total cost and incremental completion
                   for previous period.
             - Assumes no schedule slip or cost
        """
        new_expenditures = (self.total_overnight_cost) * (self.completion[-1] - self.completion[-2])
        self.xtr_expenditures.append(new_expenditures)

    def decrement_remaining_life(self):
        if self.status == 'in_service':
            self.asset_life_remaining -= 1



