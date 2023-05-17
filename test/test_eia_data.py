##########################################################################
## Copyright 2023 Argonne National Laboratory
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


from eia_data import *
import pytest
from datetime import date
import pandas as pd


def test_get_date():
    """
    Tests the get_date function.
    """

    today = date.today().strftime("%B %d, %Y")

    today = today.split(' ')

    m, d, y = get_date()

    assert m == today[0]
    assert d == today[1]
    assert y == int(today[2])

    return


def test_get_eia_generators_none_passed():
    """
    Tests the data download function where
    neither the month nor year are specified.
    """

    df = get_eia_generators()

    assert isinstance(df, pd.DataFrame)

    return


def test_get_eia_generators_month_passed():
    """
    Tests the data download function where
    only the month is passed.
    """

    with pytest.raises(ValueError):
        df = get_eia_generators(month='January')

    return



def test_get_eia_generators_year_passed():
    """
    Tests the data download function where
    only the year is passed.
    """

    with pytest.raises(ValueError):
        df = get_eia_generators(year=2022)

    return


def test_get_eia_generators_bad_year_month():
    """
    Tests the data download function where
    an invalid month and year are passed.
    """
    month = 'thermidor'
    year = 2

    with pytest.raises(ValueError):
        df = get_eia_generators(month=month, year=year)

    return


df = get_eia_generators()


def test_get_region_techs_state():
    """
    Tests the regional filter function using a
    state abbreviation.
    """

    state = 'IL'
    region_techs = get_region_techs(df, region=state)

    assert len(region_techs) > 0

    return


def test_get_region_techs_bad_state():
    """
    Tests the regional filter function using a
    non-existent state abbreviation.
    """

    state = 'QU'
    with pytest.raises(ValueError):
        region_techs = get_region_techs(df, region=state)

    return


def test_get_region_techs_county():
    """
    Tests the regional filter function using a 
    county name.
    """

    county = 'Champaign'
    region_techs = get_region_techs(df, region=county)

    assert len(region_techs) > 0

    return

def test_get_region_techs_bad_county():
    """
    Tests the regional filter function using a 
    non-existent county name.
    """

    county = 'France'
    with pytest.raises(ValueError):
        region_techs = get_region_techs(df, region=county)

    return
