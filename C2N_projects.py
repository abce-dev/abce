import pandas as pd
import sqlite3

import ABCEfunctions


class C2NProject(object):
    """ A generic C2N project object.
    """
    def __init__(self, agent_id, settings, agent_params, status, C2N_project_defs, db = None, WIP_id = None, current_pd = None, C2N_type = None):
        self.agent_id = agent_id
        self.settings = settings
        self.agent_params = agent_params
        self.status = status
        self.C2N_project_defs = C2N_project_defs
        self.C2N_type = C2N_type

        if status == "new_xtr":
            self.project_status = self.initialize_new_project()
        elif status == "ongoing_xtr":
            project_status = pd.read_sql_query(f"SELECT * FROM C2N_WIP WHERE asset_id = {WIP_id} AND pd = {current_pd}", db)
            self.project_status = project_status.to_dict(orient='records')[0]


    def initialize_new_project(self):
        if self.C2N_type not in self.C2N_project_defs.keys():
            msg = f"Unrecognized type of C2N project: {C2N_type}. Double-check your spelling in the function call, and ensure a corresponding project is specified in the C2N input file."
            raise KeyError(msg)
        else:
            C2N_project = self.C2N_project_defs[self.C2N_type]
            cpp_dnd = {"cost_remaining": C2N_project["cpp_dnd"]["cost"],
                       "time_remaining": C2N_project["cpp_dnd"]["time"]}
            cpp_nrc = {"cost_remaining": C2N_project["cpp_nrc"]["cost"],
                       "time_remaining": C2N_project["cpp_nrc"]["time"]}
            npp_ns_xtr = {"cost_remaining": C2N_project["npp_ns_xtr"]["cost"],
                          "time_remaining": C2N_project["npp_ns_xtr"]["time"]}
            npp_safety_xtr = {"cost_remaining": C2N_project["npp_safety_xtr"]["cost"],
                              "time_remaining": C2N_project["npp_safety_xtr"]["time"]}
            npp_comm = {"cost_remaining": C2N_project["npp_comm"]["cost"],
                        "time_remaining": C2N_project["npp_comm"]["time"]}


        return {cpp_dnd, cpp_nrc, npp_ns_xtr, npp_safety_xtr}


    def allocate_anpe(self, anpe):
        """ Distribute this period's authorized expenditures across all available activities.
        """
        # Track total available funds left to disburse
        remaining_anpe = anpe

        # If no activities have outstanding costs, the project is finished
        if all([activity["cost_remaining"] == 0 for activity in self.project_status.keys()]):
            # Project is done!
            pass
        elif all([activity["cost_remaining"] == 0 for activity in self.project_status.keys() if activity != "npp_comm"]):
            # Allocate all ANPE to commissioning
            
        if self.C2N_type == "greenfield_nuclear":

        elif self.C2N_type == "C2N_elec_only":

        elif self.C2N_type == "C2N_steam_with_TES":

        elif self.C2N_type == "C2N_steam_no_TES":





    def create_capex_timeline(self):
        if self.status == "new_xtr":

        elif self.status == "ongoing_xtr":




if __name__ == "__main__":
    proj = C2NProject(201, None, None)
