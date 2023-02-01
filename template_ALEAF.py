import pandas as pd
import yaml
import openpyxl

ALEAF_data_file = "./inputs/ALEAF_settings.yml"
unit_specs_data_file = "./inputs/unit_specs.yml"
settings_file = "./settings.yml"


def load_data(file_name):
    file_contents = yaml.load(
                        open(file_name, "r"),
                        Loader=yaml.FullLoader
                    )
    return file_contents


def create_ALEAF_Master_file(ALEAF_data, settings):
    # Dictionary of all tabs and their metadata
    # In final implementation, will be replaced by the ALEAF standard
    #   data schema
    tabs_to_create = {
        "ALEAF Master Setup": {
            "ABCE_tab_name": "ALEAF_Master_setup",
            "data": None
        },

        "CPLEX Setting": {
            "ABCE_tab_name": "CPLEX_settings",
            "data": None
        },

        "GLPK Setting": {
            "ABCE_tab_name": "GLPK_settings",
            "data": None
        },

        "CBC Setting": {
            "ABCE_tab_name": "CBC_settings",
            "data": None
        },

        "Gurobi Setting": {
            "ABCE_tab_name": "Gurobi_settings",
            "data": None
        },

        "HiGHS Setting": {
            "ABCE_tab_name": "HiGHS_settings",
            "data": None
        }
    }

    # Pull the ALEAF settings data into dictionaries
    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        tabs_to_create[ALEAF_tab_name]["data"] = ALEAF_data["ALEAF_Master"][tab_data["ABCE_tab_name"]]

    # Finalize the <solver> Setting tab data
    for solver_tab, tab_data in tabs_to_create.items():
        if solver_tab != "ALEAF Master Setup":
            # Set up metadata about solver settings
            solver_setting_list = ", ".join([parameter for parameter in tabs_to_create[solver_tab]["data"].keys()])

            # Set up solver_direct_mode_flag: TRUE if CPLEX, FALSE otherwise
            mode_flag = "false"
            if "CPLEX" in solver_tab:
                mode_flag = "true"

            # Create the dictionary of extra rows for all solver tabs
            solver_extra_items = {
                "solver_direct_mode_flag": mode_flag,
                "num_solver_setting": len(tabs_to_create[solver_tab]["data"]),
                "solver_setting_list": solver_setting_list
            }

            tab_data["data"].update(solver_extra_items)

    write_workbook_and_close("ALEAF_Master", tabs_to_create)


def create_ALEAF_Master_LC_GEP_file(ALEAF_data, settings):
    tabs_to_create = {
        "LC_GEP Setting": {
            "ABCE_tab_name": "LC_GEP_settings",
            "data": None
        },

        "Planning Design": {
            "ABCE_tab_name": "planning_design",
            "data": None
        },

        "Simulation Setting": {
            "ABCE_tab_name": "simulation_settings",
            "data": None
        },

        "Simulation Configuration": {
            "ABCE_tab_name": "scenario_settings",
            "data": None,
            "orient": "horizontal"
        },

        "File Path": {
            "ABCE_tab_name": "ALEAF_relative_file_paths",
            "data": None
        },

        "Scenario Reduction Setting": {
            "ABCE_tab_name": "scenario_reduction_settings",
            "data": None
        },
    }


    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        tabs_to_create[ALEAF_tab_name]["data"] = ALEAF_data["ALEAF_Master_LC_GEP"][tab_data["ABCE_tab_name"]]

    # Finalize tab data
    # Finalize "Planning Design" tab
    # Add extra items
    pd_data = tabs_to_create["Planning Design"]["data"]
    pd_extra_items = {
        "targetyear_value": pd_data["final_year_value"] - pd_data["current_year_value"],
        "load_increase_rate_value": 1, #TODO: determine A-LEAF's expected calculation
        "num_simulation_per_stage_value": (pd_data["final_year_value"] - pd_data["current_year_value"]) / pd_data["numstages_value"]
    }
    tabs_to_create["Planning Design"]["data"].update(pd_extra_items)

    # Rename items
    items_to_rename = {
        "planning_reserve_margin": "planning_reserve_margin_value"
    }
    for key, value in items_to_rename.items():
        tabs_to_create["Planning Design"]["data"][value] = tabs_to_create["Planning Design"]["data"].pop(key)


    # Finalize "Simulation Setting" tab
    # Add extra items
    ss_extra_items = {
        "scenario_name": f"ALEAF_{tabs_to_create['Simulation Setting']['data']['test_system_name']}"
    }
    tabs_to_create["Simulation Setting"]["data"].update(ss_extra_items)

    # Rename items
    items_to_rename = {
        "capex_projection_flag": "capax_projection_flag",
    }
    for key, value in items_to_rename.items():
        tabs_to_create["Simulation Setting"]["data"][value] = tabs_to_create["Simulation Setting"]["data"].pop(key)


    # Finalize "Simulation Configuration" tab
    # Add extra items
    sc_extra_items = {
        "planning_reserve_margin_value": tabs_to_create["Planning Design"]["data"]["planning_reserve_margin_value"],
        "load_increase_rate_value": tabs_to_create["Planning Design"]["data"]["load_increase_rate_value"]
    }
    tabs_to_create["Simulation Configuration"]["data"].update(sc_extra_items)

    # Rename items
    items_to_rename = {
        "scenario_name": "Scenario",
        "peak_demand": "PD",
        "carbon_tax": "CTAX",
        "RPS_percentage": "RPS",
        "wind_PTC": "PTC_W",
        "solar_PTC": "PTC_S",
        "nuclear_PTC": "PTC_N",
        "wind_ITC": "ITC_W",
        "solar_ITC": "ITC_S",
        "nuclear_ITC": "ITC_N"
    }
    for key, value in items_to_rename.items():
        tabs_to_create["Simulation Configuration"]["data"][value] = tabs_to_create["Simulation Configuration"]["data"].pop(key)


    # Finalize "Scenario Reduction Setting" tab
    srs_data = tabs_to_create["Scenario Reduction Setting"]["data"]
    data_sets = [srs_data["input_type_load_shape_flag"],
                 srs_data["input_type_load_MWh_flag"],
                 srs_data["input_type_wind_shape_flag"],
                 srs_data["input_type_wind_MWh_flag"],
                 srs_data["input_type_solar_shape_flag"],
                 srs_data["input_type_solar_MWh_flag"],
                 srs_data["input_type_net_load_MWh_flag"]
                ]
    num_data_sets = sum(1 for item in data_sets if item == "TRUE")
    srs_extra_items = {
        "num_data_set": num_data_sets
    }
    tabs_to_create["Scenario Reduction Setting"]["data"].update(srs_extra_items)

    write_workbook_and_close("ALEAF_Master_LC_GEP", tabs_to_create)


def write_workbook_and_close(base_filename, tabs_to_create):
    # Load all tab data into pandas dataframes
    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        orient = "index"
        if "orient" in tab_data.keys():
            if tab_data["orient"] == "horizontal":
                df = (pd.DataFrame.from_dict(
                         [tab_data["data"]]
                     ))
        else:
            df = (pd.DataFrame.from_dict(
                      tab_data["data"],
                      orient="index",
                      columns=["Value"]
                  ).reset_index().rename(columns={"index": "Setting"}))

        tabs_to_create[ALEAF_tab_name]["data"] = df

    # Create an ExcelWriter object to contain all tabs and save file
    writer_object = pd.ExcelWriter(f"test{base_filename}.xlsx", engine="openpyxl")

    # Write all tabs to file
    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        tab_data["data"].to_excel(
            writer_object,
            sheet_name = ALEAF_tab_name,
            index=False
        )

    # Write all changes and close the file writer
    writer_object.close()


def create_ALEAF_files():
    ALEAF_data = load_data(ALEAF_data_file)
    unit_specs_data = load_data(unit_specs_data_file)
    settings = load_data(settings_file)

    # Create the ALEAF_Master.xlsx file
    create_ALEAF_Master_file(ALEAF_data, settings)

    # Create the ALEAF_Master_LC_GEP.xlsx file
    create_ALEAF_Master_LC_GEP_file(ALEAF_data, settings)

if __name__ == "__main__":
    create_ALEAF_files()
