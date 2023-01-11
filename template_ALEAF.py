import pandas as pd
import yaml

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
    # Create the "ALEAF Master Setup" tab
    ALEAF_Master_setup = (pd.DataFrame.from_dict(
                              ALEAF_data["ALEAF_Master"]["ALEAF_Master_setup"],
                              orient="index",
                              columns=["Value"]
                          ).reset_index()
                           .rename(columns={"index": "Setting"})
                         )
    print(ALEAF_Master_setup)



def create_ALEAF_files():
    ALEAF_data = load_data(ALEAF_data_file)
    unit_specs_data = load_data(unit_specs_data_file)
    settings = load_data(settings_file)

    # Create the ALEAF_Master.xlsx file
    create_ALEAF_Master_file(ALEAF_data, settings)

if __name__ == "__main__":
    create_ALEAF_files()
