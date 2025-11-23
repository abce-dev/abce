##########################################################################
# Copyright 2023 Argonne National Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################

import pandas as pd
import yaml
from pathlib import Path


def load_data(file_name):
    if Path(file_name).suffix in [".yml", ".yaml"]:
        file_contents = yaml.load(open(file_name, "r"), Loader=yaml.FullLoader)
    elif Path(file_name).suffix == ".csv":
        file_contents = pd.read_csv(file_name)

    return file_contents


def set_unit_type_policy_adjustment(unit_type, unit_type_data, settings):
    # Initialize all units with zero policy adjustment
    carbon_tax_per_MWh = 0
    tax_credits_per_MWh = 0
    tax_credits_per_MW = 0

    if "policies" in settings["scenario"].keys():
        # If some policies are specified, determine their effect
        for policy, policy_specs in settings["scenario"]["policies"].items():
            if policy_specs["enabled"]:
                is_eligible = False
                if "eligibility" not in policy_specs.keys():
                    # If there are no eligibility criteria specified, the
                    #   policy applies to all units
                    is_eligible = True
                elif (
                    "unit_type" in policy_specs["eligibility"].keys()
                    and unit_type in policy_specs["eligibility"]["unit_type"]
                ):
                    # Check to see if the unit's type is included in a list of
                    #   eligible unit types
                    is_eligible = True
                else:
                    # Check for any other eligibility criteria
                    for criterion, value in policy_specs["eligibility"].items():
                        if (
                            criterion != "unit_type"
                            and unit_type_data[criterion] == value
                        ):
                            is_eligible = True

                if is_eligible:
                    if "CTAX" in policy:
                        carbon_tax_per_MWh -= (
                            unit_type_data["heat_rate"]
                            * unit_type_data["emissions_per_MMBTU"]
                            * policy_specs["qty"]
                        )
                    elif "PTC" in policy:
                        tax_credits_per_MWh += policy_specs["qty"]
                    elif "ITC" in policy:
                        tax_credits_per_MW += policy_specs["qty"]

    policy_results = {
        "carbon_tax_per_MWh": carbon_tax_per_MWh,
        "tax_credits_per_MWh": tax_credits_per_MWh,
        "tax_credits_per_MW": tax_credits_per_MW,
    }

    return carbon_tax_per_MWh, tax_credits_per_MWh, tax_credits_per_MW,


def compute_unit_specs_cols(unit_specs, settings):
    for unit_type, unit_type_data in unit_specs.items():
        # Ensure all units have fuel cost elements
        if not unit_type_data["uses_fuel"]:
            unit_type_data["FC_per_MMBTU"] = 0
            unit_type_data["fuel_type"] = "none"

        # Add fuel cost per MWh column
        unit_type_data["FC_per_MWh"] = (
            unit_type_data["FC_per_MMBTU"] * unit_type_data["heat_rate"]
        )

        # Add policy adjustment per MWh column
        unit_type_data["carbon_tax_per_MWh"], unit_type_data["tax_credits_per_MWh"], unit_type_data["tax_credits_per_MW"] = set_unit_type_policy_adjustment(
            unit_type, unit_type_data, settings
        )

        # Add C2N-related factors, with a default value of 0
        C2N_vals = ["cpp_ret_lead", "num_cpp_rets", "rev_head_start"]
        for val in C2N_vals:
            if val not in unit_type_data.keys():
                unit_type_data[val] = 0

    return unit_specs


def initialize_unit_specs(settings, args):
    unit_specs = load_data(
        Path(args.inputs_path)
        / settings["file_paths"]["unit_specs_data_file"]
    )

    unit_specs = compute_unit_specs_cols(unit_specs, settings)

    return unit_specs

