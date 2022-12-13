import pandas as pd
import openpyxl
import yaml

unit_sample_data_file = "./inputs/unit_specs.yml"
optional_spec_defaults = {
    "fuel_type": "none",
    "fuel_cost": 0,
    "fuel_cost_units": "$/MWh",
    "heat_rate": 0,
    "Charge_CAP": 0,
    "STOCAP": 0,
    "STOMIN": 0,
    "retirement_cost": 0
    "no_load_cost": 0,
    "start_up_cost": 0,
    "shut_down_cost": 0
    "occ_variance": 0,
    "FOR": 0
}

universally_required_specs = [
    "capacity",
    "max_PL",
    "min_PL",
    "overnight_capital_cost",
    "FOM",
    "VOM",
    "ramp_up_limit",
    "ramp_down_limit",
    "capacity_credit",
    "emissions_rate",
    "VRE_flag",
    "uses_fuel",
    "unit_life",
    "construction_duration"
]

fuel_data_values = ["fuel_type", "fuel_cost", "fuel_cost_units"]
valid_fuel_units = ["$/MWh", "$/MMBTU"]


def read_unit_data(unit_data_file):
    # Read in the unit specification data from the input file
    unit_specs = yaml.load(open(unit_data_file, "r"), Loader=yaml.FullLoader)

    return unit_specs


def check_input_specs(unit_specs):
    """
    Validate the user's unit specification file. Currently includes the
      following checks:

       * Checks whether all mandatory unit specification data types are
           provided
       * Checks whether fuel information is provided for non-VRE units
       * Checks whether a valid fuel cost unit is provided
       * Checks whether a final fuel cost value in $/MWh is computable based
           on values provided

    Function displays issues to the user as they are found, on a unit-by-unit
      basis. If a nonzero number of issues are found, ends program execution.
    """
    specs_ok = True
    for unit_type, unit_type_specs in unit_specs:
        # Check whether any universally mandatory specification values
        #   are missing from the spec for this unit type
        universal_missing_values = [value_name for value_name in universally_required_specs
                                    if value_name not in unit_type_specs.keys()]
        if len(universal_missing_values) > 0:
            logging.error(
                f"Unit type {unit_type} is missing the following unit " +
                "specification value(s), which are required for all generation " +
                "unit types:"
            )
            logging.error(universal_missing_values)
            specs_ok = False

        # Check for missing fuel information for fuel-using generators
        if unit_type_specs["uses_fuel"]:
            fuel_missing_values = [value_name for value_name in fuel_data_values
                                   if value_name not in unit_type_specs.keys()]
            if len(fuel_missing_values) > 0:
                logging.error(
                    f"Unit type {unit_type} is marked as a fuel-using " +
                    "generator, but it is missing the following fuel-related " +
                    "specification(s):"
                )
                logging.error(fuel_missing_values)
                specs_ok = False

        # Check for unknown fuel cost unit type
        if unit_type_specs["fuel_cost_units"] not in valid_fuel_units:
            logging.error(
                f"Unit type {unit_type} has an invalid fuel cost unit type: " +
                f"{unit_type_specs['fuel_cost_units']}. Please use one of " +
                f"the following: {valid_fuel_units}."
            )
            specs_ok = False

        # Check for missing heat rate when fuel is given in $/MMBTU
        if (unit_type_specs["fuel_cost_units"] == "$/MMBTU" and
            "heat_rate" not in unit_type_specs.keys()):
            logging.error(
                f"Unit type {unit_type} has fuel cost specified in $/MMBTU, " +
                "but has no specified heat rate (MWh/MMBTU)."
            )
            specs_ok = False

    if not specs_ok:
        logging.info(
            "Please rectify the issues listed before re-running ABCE."
        )
        raise ValueError


def fill_unit_spec_defaults(unit_specs):
    for unit_type, unit_type_specs in unit_specs:
        for data_value in optional_spec_defaults.keys():
            if data_value not in unit_type_specs.keys():
                unit_type_specs[data_value] = optional_spec_defaults[data_value]


def finalize_unit_spec_data(unit_specs):
    # Compute all final fuel costs in units of $/MWh
    unit_specs = compute_fuel_costs_dpMWh(unit_specs)


def compute_fuel_costs_dpMWh(unit_specs):
    # Ensure that all units have a valid fuel cost value in $/MWh
    for unit_type, unit_type_specs in unit_specs:
        if unit_type_specs["fuel_cost_units"] == "$/MWh":
            unit_type_specs["FC_per_MWh"] = unit_type_specs["fuel_cost"]
        else:
            unit_type_specs["FC_per_MWh"] = unit_type_specs["fuel_cost"] * unit_type_specs["heat_rate"]


def create_simulation_specification_file(unit_data_file):
    # Creates A-LEAF "LC_GEP" style input file

    # Read in yaml data from input file
    unit_specs = read_unit_data(unit_data_file)

    # Check to ensure that required unit specification values are all provided
    check_input_specs(unit_specs)

    # Fill in unspecified values with appropriate defaults
    unit_specs = fill_unit_spec_defaults(unit_specs)


if __name__ == "__main__":
    create_simulation_specification_file(unit_sample_data_file)
