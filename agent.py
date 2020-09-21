from mesa import Agent, Model
from mesa.time import RandomActivation
import yaml
import math

# import local modules
import generator as gen

class GenCo(Agent):
    """ A utility company with a certain number of generation assets.
    """
    def __init__(self, genco_id, model, existing_portfolio):
        """ Initialize a GenCo class object.

            Detailed Description
            --------------------
            Initializes a GenCo object based on the provided parameters. Uses
            baseline mesa capabilities to set up an Agent-type object.

            Sets up the unique identifier number and the link to the invoking
            GridModel (`model`). Also sets up the agent's starting portfolio.

            Parameters
            ----------
            genco_id : int
                A unique ID number for the agent, assigned sequentially
                by the Model which owns the agents.
            model : GridModel (mesa)
                A mesa GridModel object which creates and passes data to
                all agents.
            existing_portfolio : dict (of dicts)
                A dictionary specifying all assets comprising the agent's
                existing portfolio at time t=0.

            Returns
            -------
            None
        """

        super().__init__(genco_id, model)
        self.set_up_portfolio(existing_portfolio)
        self.model = model

    def set_up_portfolio(self, existing_portfolio):
        """Set up Generator units for each asset in GenCo's assigned portfolio

           Detailed Description
           --------------------
           Sets up Generator objects for each unit specified in the GenCo's
           initial portfolio. The structure of the starting portfolio is still
           under development.

           For each unit detected in the `existing_portfolio` dictionary,
           this function determines its generator type, and sets its unit ID
           according to the next available ID number as provided by the
           model's id_register object. A new unit of the appropriate type is
           initialized, and added to the agent's `self.portfolio` dictionary.

           Parameters
           ----------
           existing_portfolio : dict of dicts
               A dict of dicts (as extracted from a yaml file) describing the
               number, types, and attributes of assets in the GenCo's
               starting portfolio at time t=0.

           Returns
           -------
           None
        """
        self.portfolio = dict()
        for unit_num in existing_portfolio.keys():
            gtype = existing_portfolio[unit_num]['gtype']
            unit_id = self.model.id_register.get_next_available_id()
            new_unit = gen.Generator(id_num=unit_id, gtype=gtype, completion=1)
            new_unit.completion = [1]
            self.portfolio[unit_id] = new_unit
            self.model.id_register.add_unit(self.unique_id, unit_id)

    def step(self):
        """Controller function to activate all agent behaviors at each time step.

           Detailed Description
           --------------------
           This function is automatically invoked by the GridModel, once per
           agent per time step. This function does the following:
             - Updates the current step number
             - Updates the status of all WIP and in-service generation units
                 owned by the agent
             - Prints an identifying statement to track agent activity
             - Updates the agent's demand forecast
             - Totals up the agent's current generation capacity
             - Evaluates whether additional assets are needed, and cues up new
               construction projects if so

           Parameters
           ----------
           none

           Returns
           -------
           None
        """
        self.set_current_step()
        for id_num in self.portfolio:
            self.portfolio[id_num].step()
        num_units = len(self.portfolio.keys())
        print(f'I\'m GenCo # {self.unique_id}, and I have {num_units} units.')
        self.get_demand_forecast()
        #self.forecast_demand()
        self.evaluate_current_capacity()
        self.assess_supply_adequacy()

    def evaluate_current_capacity(self):
        """Sum up all current generation capacity owned by this agent.

        Detailed Description
        --------------------
        This function sums up the total capacity of all generation units
        owned by the utility, including both in-service units and units
        under construction ('wip' status). This value can then be used to
        determine whether additional capacity investments are needed.

        Parameters
        ----------
        none

        Returns
        -------
        None

        """
        self.total_capacity = 0
        for unit in self.portfolio.keys():
            self.total_capacity += self.portfolio[unit].capacity
        print(f'Total capacity currently installed = {self.total_capacity}.')

    def assess_supply_adequacy(self):
        """Determine whether future demand exceeds future supply.

        Detailed Description
        --------------------
        Compares current total capacity to each demand forecast quantity in
        the agent's demand visibility window. If a period is observed where
        demand will be higher than current supply, the agent will start new
        construction projects in order to compensate.

        Once a single period with forecasted unmet demand is detected, this
        function will not examine further future periods.

        Parameters
        ----------
        none

        Returns
        -------
        None

        """
        for pd in self.demand_forecast:
            if pd > self.total_capacity:
                # Additional capacity needed!
                print("Demand will be higher than current capacity. Building additional capacity...")
                # Calculate number of new units required
                # Current: divide predicted unserved demand (pd - current capacity)
                #   by the per-unit capacity of the unit_1 type.
                # To be replaced when unit choice behaviors are modeled in
                #   more detail.
                num_new_units = int(math.ceil((pd - self.total_capacity) / float(self.model.unit_data['unit_1']['capacity'])))
                self.build_new_units(num_new_units)
                return

    def build_new_units(self, num_new_units):
        """Start a new construction project, and add to agent's portfolio.

        Detailed Description
        --------------------
        

        Parameters
        ----------
        num_new_units : int
            Number of new units to build.

        Returns
        -------
        None

        """
        for i in range(num_new_units):
            unit_id = self.model.id_register.get_next_available_id()
            new_unit = gen.Generator(id_num=unit_id, gtype='unit_1')
            self.portfolio[unit_id] = new_unit
            self.model.id_register.add_unit(self.unique_id, unit_id)

    def set_current_step(self):
        """Obtain the current step number from the model.

        """
        self.current_step = self.model.current_step
        print(f'\n\n Current step = {self.current_step}')

    def get_demand_forecast(self):
        """ Obtain the current demand visibility window from the model.

            Detailed Description
            --------------------
            The model provides a fixed look-ahead window of the deterministic
            demand data provided in an input file. The agent retrieves this
            information and uses it for planning.

        """
        if self.model.demand_NTF != []:
            self.demand_forecast = self.model.demand_NTF
            self.current_demand = self.demand_forecast[0]

    def forecast_demand(self):
        """Determine lead time until next unmanageable demand increase.

           Detailed Description
           --------------------
           Loops over the current demand visibility window, and determines
           whether (and when) a capacity shortfall will occur.

           The time step at which this will happen is saved in the variable
           `self.next_demand_increase_period`.

           This function is not currently used, but its capability will be
           needed to allow agents to decide which units can be built in time
           to meet the additional demand.

           Parameters
           ----------
           none

           Returns
           -------
           None

        """
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

