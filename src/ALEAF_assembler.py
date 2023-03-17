import yaml
import logging

ALEAF_schema_file = "/home/biegelk/abce/inputs/ALEAF_settings_schema.yml"
ALEAF_settings_file = "/home/biegelk/abce/inputs/ALEAF_settings.yml"
settings_file = "/home/biegelk/abce/settings.yml"

def load_data(file_name):
    data = yaml.load(
               open(file_name, "r"),
               Loader=yaml.FullLoader
           )

    return data


def levelize_data(data):
    data_dict = {}
    for key, value in data.items():
        # Skip any parts of the hierarchy encoding tabular data or otherwise
        #   irreducible hierarchies; due to repetition of leaf node keys, 
        #   these need to be processed specially
        if key in ["unit_specs", "ATB_search_settings", "policies"]:
            continue
        if type(value) is not dict:
            data_dict[key] = value
        else:
            data_dict.update(levelize_data(value))

    return data_dict



def initialize_file(file_data, ALEAF_settings, ABCE_settings):
    # Pop the file's exact name from the dictionary, allowing systematic
    #   processing of the other keys, which all indicate excel tabs
    file_name = file_data.pop("default_filename")
    logging.info(f"Initializing file: {file_name}")

    # Set up each of the tabs in turn
    tabs = {}
    for tab_name, tab_data in file_data.items():
        tabs[tab_name] = initialize_tab(tab_data, ALEAF_settings, ABCE_settings)

    # Set up final file data structure
    file_data = {
        "file_name": file_name,
        "tabs": tabs
    }

    return file_data


def initialize_tab(tab_data, ALEAF_settings, ABCE_settings):
    # Pop the tab's ALEAF-style name from the dictionary
    tab_ALEAF_name = tab_data.pop("tab_name")

    logging.info(f"Initializing {tab_ALEAF_name}")

    # Pop the order of the fields in this tab, if available
    tab_data, field_order = get_field_order(tab_data)

    # Pop the orientation of the fields in this tab (vertical is default)
    tab_data, tab_orientation = get_tab_orientation(tab_data)

    # Initialize a dictionary to track any problems setting up the data
    problems = {}

    # Iterate over the remaining fields in tab_data
    for field_name, field_data in tab_data.items():
        # Set the field's value, if possible
        field_data["field_value"] = set_field_value(
                                        field_name,
                                        field_data,
                                        ALEAF_settings,
                                        ABCE_settings,
                                        problems
                                    )

        # Complete initial validation of field value, except for fields
        #   with a current value of None (these are either already broken,
        #   or awaiting actual values with dynamic fill routines, to be applied
        #   later)
        field_problems = validate_field_value(field_name, field_data)
        if len(field_problems) > 0:
            problems[field_name] = field_problems

    # Return the data from this tab
    interior_tab_data = {
        "tab_ALEAF_name": tab_ALEAF_name,
        "field_order": field_order,
        "tab_orientation": tab_orientation,
        "tab_data": tab_data,
        "problems": problems
    }

    return interior_tab_data


def get_field_order(tab_data):
    # Pop the tab's field order from the dictionary, if found
    try:
        field_order = tab_data.pop("field_order")
    except KeyError:
        field_order = None

    return tab_data, field_order


def get_tab_orientation(tab_data):
    # Pop the tab's orientation from the dictionary, if found. If not
    #   specified, assume vertical
    try:
        tab_orientation = tab_data.pop("field_orientation")
    except KeyError:
        tab_orientation = "vertical"

    return tab_data, tab_orientation


def set_field_value(field_name, field_data, ALEAF_settings, ABCE_settings, problems):
    # Set the field's value
    # Settings file takes the highest precedence
    if field_name in ABCE_settings.keys():
        field_value = ABCE_settings[field_name]
    elif field_name in ALEAF_settings.keys():
        field_value = ALEAF_settings[field_name]
    elif "default_value" in field_data.keys():
        field_value = field_data["default_value"]
    elif "input_type" in field_data.keys() and field_data["input_type"] == "dynamic":
        logging.info(field_name)
        field_value = None
    else:
        # Log a problem: this field cannot be assigned a value
        problems[field_name] = "no_value_assignable"
        field_value = None

    return field_value


def validate_field_value(field_name, field_data):
    problems = []
    field_value = field_data["field_value"]  # for ease of reference

    # Create a flag variable to ensure inapplicable checks are not performed
    skip = False

    # If the value is None, validation is either impossible or will be handled
    #   in a later iteration
    if field_value is not None:
        # Check the value's type
        if "types" in field_data.keys():
            # Make a list of type objects, converted from the strs specified
            valid_types = [getattr(__builtins__, type_name) for type_name in field_data["types"]]
            if type(field_value) not in valid_types:
                # If int is an allowed type and the value is convertible to
                #   int, overwrite field_value with an int of itself
                if int in valid_types:
                    try:
                        field_value = int(field_value)
                    except ValueError:
                        pass
                elif float in valid_types and type(field_value) not in valid_types:
                    # If float is an allowed type and the value is convertible
                    #   to float, overwrite field_value with a float of itself
                    try:
                        field_value = float(field_value)
                    except ValueError:
                        pass
                elif str in valid_types and type(field_value) not in valid_types:
                    # Try str, if allowed
                    try:
                        field_value = str(field_value)
                    except ValueError:
                        pass
                elif bool in valid_types and type(field_value) not in valid_types:
                    # Try bool, if allowed
                    try:
                        field_value = bool(field_value)
                    except:
                        pass
                else:
                    # If the type can't be automatically fixed, mark it as a 
                    #   problem and continue to the next field_value
                    problems.append("type")
                    skip = True

        # Check for a list of allowed values
        if "allowed_values" in field_data.keys() and not skip:
            if field_value not in field_data["allowed_values"]:
                problems.append("not_allowed")
                skip = True

        # Check for specified min and max values
        if "min_value" in field_data.keys() and not skip:
            if field_value < field_data["min_value"]:
                problems.append("below_min")
        if "max_value" in field_data.keys() and not skip:
            if field_value > field_data["max_value"]:
                problems.append("above_max")

    return problems


def initialize_ALEAF_schema():
    ALEAF_schema = load_data(ALEAF_schema_file)
    ALEAF_settings_data = load_data(ALEAF_settings_file)
    ABCE_settings_data = load_data(settings_file)

    # Retrieve the non-tabular ALEAF and ABCE settings data as partially
    #    flattened dictionaries
    ALEAF_settings = levelize_data(ALEAF_settings_data)
    ABCE_settings = levelize_data(ABCE_settings_data)

    # Iterate through the files in the schema
    ALEAF_files = {}
    for file_type, file_data in ALEAF_schema.items():
        ALEAF_files[file_type] = initialize_file(file_data, ALEAF_settings, ABCE_settings)
        ALEAF_files[file_type] = instantiate_dynamic_field_values(file_type, file_data)

    return ALEAF_files


def instantiate_dynamic_field_values(file_type, file_data):
    if file_type == "ALEAF_Master":
        for tab_name, tab_data in file_data.items():
            if "solver_direct_mode_flag" in tab_data.keys():
                if "CPLEX" in tab_name:
                    file_data[tab_name]["solver_direct_mode_flag"]["field_value"] = "True"
                else:
                    file_data[tab_name]["solver_direct_mode_flag"]["field_value"] = "False"

            if "solver_setting_list" in tab_data.keys():
                # Get all settings aside from those which are ignored
                ignore = ["solver_direct_mode_flag", "solver_setting_list", "num_solver_setting"]
                all_solver_settings = [key for key in tab_data.keys() if key not in ignore]

                file_data[tab_name]["solver_setting_list"]["field_value"] = ", ".join(all_solver_settings)
                file_data[tab_name]["num_solver_setting"]["field_value"] = len(all_solver_settings)

    elif file_type == "ALEAF_Master_LC_GEP":
        print(file_data["simulation_settings"])

    elif file_type == "ALEAF_portfolio":
        pass

    return file_data






if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    initialize_ALEAF_schema()