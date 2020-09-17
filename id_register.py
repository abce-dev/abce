# A class to hold the register of all generation unit IDs currently in use

class ID_register(object):
    def __init__(self):
        self.register = list()
        self.next_id = 100

    def add_unit(self, agent_id, unit_id):
        self.register.append(unit_id)
        self.update_next_id()

    def update_next_id(self):
        ''' Update the next available project ID, ensuring no collisions with
              existing units. If the next ID is less than any entries in
              the register, reset it to be one greater than the largest
              ID number currently in the register.
        '''
        while self.next_id in self.register:
            self.next_id += 1
        if self.next_id > min(self.register):
            self.next_id = max(self.register) + 1

    def get_next_available_id(self):
        assigned_id = self.next_id
        self.update_next_id()
        return assigned_id
