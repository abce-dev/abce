# A class to hold the register of all generation unit IDs currently in use

class ID_register(object):
    """A class to provide a public register of all generation unit ID numbers
       which are currently claimed, and to assign ID numbers to new units.
    """
    def __init__(self):
        """Creates an ID_register object, and sets the starting ID number
           to 100.
        """
        self.register = list()
        self.next_id = 100

    def register_unit(self, agent_id, unit_id = None):
        """Register a new unit, and update the next unclaimed ID.

           Detailed Description
           --------------------
           Add a unit to the register, using either a user-specified unit ID
           number or an automatically-generated ID from the update_next_id()
           method in the ID_register class.

           Parameters
           ----------
           agent_id: (int)
               The unique ID of the agent which owns the generator.

           unit_id: (int, optional)
               A unique ID to assign to the unit. If not specified, or if the
               given ID is already taken, automatically generate a valid
               unclaimed ID based on the existing unit register.

           Returns
           -------
           assigned_id: (int)
               A valid unclaimed ID, recorded in the ID_register object
        """
        if unit_id is None or unit_id in self.register:
            assigned_id = self.next_id
        elif unit_id is not None:
            assigned_id = unit_id
        else:
            print('Something went wrong when trying to assign ID {unit_id} to agent {agent_id}.')
            exit()
        self.register.append(assigned_id)
        self.update_next_id()
        return assigned_id

    def update_next_id(self):
        """Update the next available project ID, ensuring no collisions with
           previously-assigned unit IDs..
        """
        self.next_id = max(self.register) + 1


