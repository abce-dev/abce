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

    def add_unit(self, agent_id, unit_id = None):
        """Add a unit to the register, and change the next available ID.
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
        """Update the next available project ID.

           Detailed Description
           --------------------
           Update the next available project ID, ensuring no collisions with
           existing units.

           This function is called by both the add_unit() and the
           get_next_available_id() functions, so the next available unit ID
           is always guaranteed to be a valid unclaimed ID after those
           operations.

           Parameters
           ----------
           none

           Returns
           -------
           None
        """
        self.next_id = max(self.register) + 1

#        if len(self.register) == 0:
#            self.next_id = 100
#        else:
#            self.next_id = max(self.register) + 1

    def get_next_available_id(self):
        """Return the next unclaimed unit ID, and then update the next claimable
           unit ID number.

           Parameters
           ----------
           None

           Returns
           -------
           assigned_id : int
              The next (sequentially) valid unclaimed unit ID number, to be
              assigned to a new generation unit.
        """
        assigned_id = self.next_id
        self.update_next_id()
        return assigned_id
