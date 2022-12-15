import pandas as pd
import openpyxl
import yaml
import logging

unit_specs_schema_file = "./unit_specs_schema.yml"
unit_sample_data_file = "./inputs/unit_specs.yml"

def read_specification_schema(unit_specs_spec_file):
    # Read in the schema describing all allowable data values in the unit_spec
    #   and important validation information about each
    unit_specs_schema = yaml.load(
                            open(unit_specs_spec_file, "r"),
                            Loader=yaml.FullLoader
                        )

    # Convert all allowable data types into their Python data type instead
    #   of strings
    for value_type in unit_specs_schema.keys():
        unit_specs_schema[value_type]["types"] = (
            list(map(lambda x: getattr(__builtins__, x),
                     unit_specs_schema[value_type]["types"]
            ))
        )

    return unit_specs_schema


def read_unit_data(unit_data_file):
    # Read in the unit specification data from the input file
    unit_specs = yaml.load(
                     open(unit_data_file, "r"),
                     Loader=yaml.FullLoader
                 )

    return unit_specs


def validate_input_specs(unit_specs_schema, unit_specs):
    """
    Validate the user's unit specification file. Currently includes the
      following checks:

       * Checks whether all mandatory unit specification data types are
           provided
       * Checks whether any unknown/nonstandard data value names are given in
           the unit's specification block
       * Checks whether input data values are of an allowed type
       * Checks whether input data values adhere to all value range
           restrictions (for values of allowed types only)
       * Checks whether fuel information is provided for non-VRE units
       * Checks whether a final fuel cost value in $/MWh is computable based
           on values provided

    Function displays issues to the user as they are found, on a unit-by-unit
      basis. If a nonzero number of issues are found, ends program execution.
    """
    specs_ok = True

    for unit_type, unit_type_specs in unit_specs.items():
        # Check whether any universally mandatory specification values
        #   are missing from the spec for this unit type
        universal_missing_values = [value_type for value_type in unit_specs_schema.keys()
                                    if "required" in unit_specs_schema[value_type].keys()
                                    and unit_specs_schema[value_type]["required"]
                                    and value_type not in unit_type_specs.keys()
                                   ]
        if len(universal_missing_values) > 0:
            logging.error(
                f"Unit type {unit_type} is missing the following unit " +
                "specification value(s), which are required for all generation " +
                "unit types:"
            )
            logging.error(universal_missing_values)
            logging.error(
                "Please add these values to this unit's specification. \n"
            )
            specs_ok = False


        # Check for provided data values whose names don't match anything in 
        #   the unit specs schema
        unknown_data_values = [value_type for value_type in unit_type_specs.keys()
                               if value_type not in unit_specs_schema.keys()
                              ]
        if len(unknown_data_values) > 0:
             logging.error(
                 f"Unit type {unit_type} has the following data values in " +
                 f"its unit specification which do not correspond to " +
                 f"allowable data values:"
             )
             logging.error(unknown_data_values)
             logging.error("Check spelling and standard data value names. \n")
             specs_ok = False


        # Check whether provided data values are of an allowed type
        wrong_type_values = [value_type for value_type in unit_type_specs.keys()
                             if value_type not in unknown_data_values
                             and type(unit_type_specs[value_type]) not in unit_specs_schema[value_type]["types"]
                            ]
        if len(wrong_type_values) > 0:
            logging.error(
                f"Unit type {unit_type} has data of an incorrect type " +
                "provided for the following values:"
            )
            for value_type in wrong_type_values:
                logging.error(f"{value_type}: ")
                logging.error(f"Types allowed: {unit_specs_schema[value_type]['types']}")
                logging.error(f"Type provided: {type(unit_type_specs[value_type])} \n")


        # Check whether provided values meet the allowed value range
        for value_type in unit_type_specs.keys():
            if value_type not in wrong_type_values and value_type not in unknown_data_values:
                if "allowed_values" in unit_specs_schema[value_type].keys():
                    if unit_type_specs[value_type] not in unit_specs_schema[value_type]["allowed_values"]:
                        logging.error(
                            f"Unit type {unit_type} has an invalid data value " +
                            f"provided for {value_type}."
                        )
                        logging.error(
                            f"Allowed values for this {value_type}: " +
                            f"{unit_specs_schema[value_type]['allowed_values']}"
                        )
                        logging.error(
                            f"User provided value: {unit_type_specs[value_type]}\n"
                        )
                        specs_ok = False
                if "min_value" in unit_specs_schema[value_type].keys():
                    if unit_type_specs[value_type] < unit_specs_schema[value_type]["min_value"]:
                        logging.error(
                            f"Unit type {unit_type} has an invalid data value " +
                            f"provided for {value_type}."
                        )
                        logging.error(
                            f"Minimum allowed value for {value_type}: " +
                            f"{unit_specs_schema[value_type]['min_value']}"
                        )
                        logging.error(
                            f"User provided value: {unit_type_specs[value_type]}\n"
                        )
                        specs_ok = False
                if "max_value" in unit_specs_schema[value_type].keys():
                    if unit_type_specs[value_type] > unit_specs_schema[value_type]["max_value"]:
                        logging.error(
                            f"Unit type {unit_type} has an invalid data value " +
                            f"provided for {value_type}."
                        )
                        logging.error(
                            f"Maximum allowed value for {value_type}: " +
                            f"{unit_specs_schema[value_type]['max_value']}"
                        )
                        logging.error(
                            f"User provided value: {unit_type_specs[value_type]}\n"
                        )
                        specs_ok = False


        # Check for missing fuel information for fuel-using generators
        if unit_type_specs["uses_fuel"]:
            fuel_missing_values = [value_type for value_type in unit_specs_schema.keys()
                                   if "fuel_related" in unit_specs_schema[value_type].keys()
                                   and unit_specs_schema[value_type]["fuel_related"]
                                   and value_type not in unit_type_specs.keys()]
            if len(fuel_missing_values) > 0:
                logging.error(
                    f"Unit type {unit_type} is marked as a fuel-using " +
                    "generator, but it is missing the following fuel-related " +
                    "specification(s):"
                )
                logging.error(fuel_missing_values)
                logging.error(
                    "Please add these data values to this unit's " +
                    "specification. \n"
                )
                specs_ok = False


        # Check for missing heat rate when fuel is given in $/MMBTU
        if ("fuel_cost_units" in unit_type_specs.keys()
            and unit_type_specs["fuel_cost_units"] == "$/MMBTU"
            and "heat_rate" not in unit_type_specs.keys()):
            logging.error(
                f"Unit type {unit_type} has fuel cost specified in $/MMBTU, " +
                "but has no specified heat rate (MWh/MMBTU)."
            )
            specs_ok = False


    if not specs_ok:
        logging.error(
            "Please rectify these unit specification issues before " +
            "re-running ABCE."
        )
        raise SystemExit


def fill_unit_spec_defaults(unit_specs_schema, unit_specs):
    # For any data value not given in the specification: if that value type has
    #   an available default value in the schema, fill in the default value
    for unit_type, unit_type_specs in unit_specs.items():
        for value_type in unit_specs_schema.keys():
            if value_type not in unit_type_specs.keys() and "default_value" in unit_specs_schema[value_type].keys():
                unit_type_specs[value_type] = unit_specs_schema[value_type]["default_value"]

    return unit_specs


def finalize_unit_specs_data(unit_specs_schema, unit_specs):
    # If any unit specification values are given as "ATB", fill those values
    #   from the ATB data file


    # Compute all final fuel costs in units of $/MWh
    unit_specs = compute_fuel_costs_dpMWh(unit_specs)

    # Set up any data values which are propagated from other values
    unit_specs = set_propagated_values(unit_specs_schema, unit_specs)

    # Set up any A-LEAF-only data values with only one valid value
    unit_specs = set_immutable_ALEAF_values(unit_specs)

    return unit_specs


def compute_fuel_costs_dpMWh(unit_specs):
    # Ensure that all units have a valid fuel cost value in $/MWh
    for unit_type, unit_type_specs in unit_specs.items():
        if unit_type_specs["fuel_cost_units"] == "$/MWh":
            unit_type_specs["FC_per_MWh"] = unit_type_specs["fuel_cost"]
        else:
            unit_type_specs["FC_per_MWh"] = unit_type_specs["fuel_cost"] * unit_type_specs["heat_rate"]

    return unit_specs


def set_propagated_values(unit_specs_schema, unit_specs):
    """
    Certain A-LEAF-specific data values should be propagated from other input
      values, if they are not set directly by the user, including:
        - unit_group (UNITGROUP)
        - unit_category (UNIT_CATEGORY)
        - dispatchable (Commitment)
        - emissions_per_MWh (Emission)
    """
    for unit_type, unit_type_specs in unit_specs.items():
        if "unit_group" not in unit_type_specs.keys():
            unit_type_specs["unit_group"] = unit_type
        if "unit_category" not in unit_type_specs.keys():
            unit_type_specs["unit_category"] = unit_type
        if "dispatchable" not in unit_type_specs.keys():
            if unit_type_specs["is_VRE"]:
                unit_type_specs["dispatchable"] = "FALSE"
            else:
                unit_type_specs["dispatchable"] = "TRUE"
        if "emissions_per_MWh" not in unit_type_specs.keys():
            unit_type_specs["emissions_per_MWh"] = unit_type_specs["emissions_per_MMBTU"] * unit_type_specs["heat_rate"]

    return unit_specs


def set_immutable_ALEAF_values(unit_specs):
    """
    Certain A-LEAF-specific data values must always take on the same values:
        - INVEST_FLAG: FALSE
        - RET_FLAG: FALSE
        - Integrality: TRUE
        - Outages: FALSE

    These data values are not settable by the user, because they have only one
      valid value in the ABCE context.
    """

    for unit_type, unit_type_specs in unit_specs.items():
        unit_type_specs["INVEST_FLAG"] = "FALSE"
        unit_type_specs["RET_FLAG"] = "FALSE"
        unit_type_specs["Integrality"] = "TRUE"
        unit_type_specs["Outages"] = "FALSE"

    return unit_specs






def initialize_unit_specifications(unit_specs_schema_file, unit_data_file):
    """
    Load the unit specifications from the input file, validate all input
      values, and finalize all specification values.
    """

    # Read in specification schema for all unit specs
    unit_specs_schema = read_specification_schema(unit_specs_schema_file)

    # Read in yaml data from input file
    unit_specs = read_unit_data(unit_data_file)

    # Check to ensure that required unit specification values are all provided
    validate_input_specs(unit_specs_schema, unit_specs)

    # Fill in unspecified values with appropriate defaults
    unit_specs = fill_unit_spec_defaults(unit_specs_schema, unit_specs)

    # Finalize unit_specs data
    unit_specs = finalize_unit_specs_data(unit_specs_schema, unit_specs)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    initialize_unit_specifications(unit_specs_schema_file, unit_sample_data_file)
