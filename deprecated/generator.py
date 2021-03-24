# A class for generation units
import yaml
import mesa
import financial_statement as fs

class Generator(object):
    # DEPRECATED
    def __init__(self, world_model, agent, id_num, gtype, completion=0):
        # Universal identification number, public knowledge
        self.id = id_num
        self.model = world_model # System model to which agents and generators belong
        self.agent = agent
        # Set up unit type, and get initial parameters from the units.yml file
        self.type = gtype
        unit = world_model.unit_data.loc[self.type]
        self.capacity = unit['capacity']
        self.overnight_cost = unit['overnight_cost']
        self.xtr_lead_time = unit['xtr_lead_time']
        self.completion_rate_type = unit['completion_rate_type']
        self.variable_cost = unit['variable_cost']
        self.fixed_cost = unit['fixed_cost']
        self.asset_life_remaining = unit['useful_life']

        # Set up parameters for WIP projects, if needed
        self.completion = list([completion])
        if self.completion[-1] == 1:
            self.status = 'in_service'
            self.proj_completion_date = 0
        else:
            self.status = 'wip'
            self.proj_completion_date = self.model.current_step + self.xtr_lead_time
            if self.completion_rate_type == 'linear':
                pass
        self.xtr_expenditures = list() # Blank list of periodic construction capital expenses

        # Initialize own FinancialStatement object
        self.fs = fs.GeneratorFS(model=self.model, agent=self.agent, generator=self)


    def step(self):
        if self.status == 'in_service':
            self.decrement_remaining_life()
        if self.status == 'wip':
            self.update_xtr_progress()
            self.update_xtr_expenses()
        self.fs.step()
        self.update_xtr_status()


    def update_xtr_progress(self):
        """Update the total % completion to date of the project.

            Detailed Description
            --------------------
            Increment the unit's construction completion progress for the most
            recent period. If this additional increment would put the
            completion quantity over 1, set the final completion to 1.

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
        # Simple prediction of construction time remaining:
        # Divide total quantity of progress left to complete (1 - current_completion) by
        #   difference between current completion and last period's completion
        #   (i.e. simple linear extrapolation of most recent completion rate)
        self.proj_completion_date = (1 - self.completion[-1]) / (self.completion[-1] - self.completion[-2]) + self.model.current_step
        # Clean up self.completion if progress overshoots 100%
        if self.completion[-1] >= 1:
            self.completion[-1] = 1


    def update_xtr_status(self):
        """Check whether the project's construction completion status should
           be updated from 'wip' to 'in_service'.
        """
        if self.completion[-1] == 1 and self.status == 'wip':
            self.status = 'in_service'
            print(f'Project {self.id} completed and entering service.')
        elif self.status == 'wip':
            print(f'Project {self.id} at {self.completion[-1]} completion.')


    def update_xtr_expenses(self):
        """Incur construction expenses related to most recent project progress.

            Current behavior:
             - Incremental cost is a linear function of total cost and incremental completion
                   for previous period.
             - Assumes no schedule slip or cost
        """
        new_expenditures = (float(self.overnight_cost)) * (self.completion[-1] - self.completion[-2])
        self.xtr_expenditures.append(new_expenditures)

    def decrement_remaining_life(self):
        if self.status == 'in_service':
            self.asset_life_remaining -= 1


