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
    # Create an ExcelWriter object to contain all tabs and save file
    writer = pd.ExcelWriter("testALEAF_Master.xlsx", engine="openpyxl")

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

    # Load all tab data into pandas dataframes
    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        df = (pd.DataFrame.from_dict(
                  ALEAF_data["ALEAF_Master"][tab_data["ABCE_tab_name"]],
                  orient="index",
                  columns=["Value"]
              ).reset_index().rename(columns={"index": "Setting"}))

        tabs_to_create[ALEAF_tab_name]["data"] = df

    # Finalize the <solver> Setting tab data
    for solver_tab, tab_data in tabs_to_create.items():
        if solver_tab != "ALEAF Master Setup":
            # Set up metadata about solver settings
            solver_setting_list = ", ".join([parameter for parameter in tabs_to_create[solver_tab]["data"].iloc[:, 0]])

            # Set up solver_direct_mode_flag: TRUE if CPLEX, FALSE otherwise
            mode_flag = "FALSE"
            if "CPLEX" in solver_tab:
                mode_flag = "TRUE"

            # Create the dictionary of extra rows for all solver tabs
            solver_extra_rows = {
                "solver_direct_mode_flag": mode_flag,
                "num_solver_setting": len(tabs_to_create[solver_tab]["data"]),
                "solver_setting_list": solver_setting_list
            }

            # Add extra rows to all solver tabs
            for setting, value in solver_extra_rows.items():
                tabs_to_create[solver_tab]["data"].loc[len(tabs_to_create[solver_tab]["data"])] = [setting, value]

    # Write all tabs to file
    for ALEAF_tab_name, tab_data in tabs_to_create.items():
        tab_data["data"].to_excel(
            writer,
            sheet_name = ALEAF_tab_name,
            index=False
        )

    # Write all changes and close the file writer
    writer.close()


def create_ALEAF_files():
    ALEAF_data = load_data(ALEAF_data_file)
    unit_specs_data = load_data(unit_specs_data_file)
    settings = load_data(settings_file)

    # Create the ALEAF_Master.xlsx file
    create_ALEAF_Master_file(ALEAF_data, settings)

if __name__ == "__main__":
    create_ALEAF_files()
