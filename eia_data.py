"""
BSD 3-Clause License

Copyright (c) 2021, Advanced Reactors and Fuel Cycles
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=================================================================================

This code is borrowed from `pygenesys` 
https://github.com/arfc/pygenesys/blob/main/pygenesys/data/eia_data.py

author: Sam Dotson
"""


import pandas as pd
import numpy as np
import os
from datetime import date

months = ['january',
          'february',
          'march',
          'april',
          'may',
          'june',
          'july',
          'august',
          'september',
          'october',
          'november',
          'december']


def get_date():
    """
    Returns the current day, month, and year.
    """
    today = date.today().strftime("%B %d, %Y")

    today = today.split(' ')

    month = today[0]
    day = today[1]
    year = today[2]

    return month, day, int(year)


def get_eia_generators(month=None, year=None):
    """
    This function returns a pandas dataframe containing information on
    all electric generators in the United States from a recent EIA form 
    860M. Note: EIA forms are often delayed therefore this function
    subtracts four months from the current month to guarantee the file
    exists to be downloaded.

    Parameters
    ----------
    month : string
        The month of interest
    year : string
        The year of interest

    Returns
    -------
    df : pandas dataframe
        Holds the data from EIA form 860M
    """
    columns = [
        'Entity ID',
        'Entity Name',
        'Plant Name',
        'Sector',
        'Plant State',
        'Nameplate Capacity (MW)',
        'Technology',
        'Operating Year',
        'Status',
        'Balancing Authority Code',
        'County'
    ]

    # initialize with invalid options
    m = 'thermidor'
    y = 2

    if (month is None) and (year is None):
        m, d, y = get_date()
        month_idx = months.index(m.lower())
        month_idx -= 4
        if month_idx < 0:
            y -= 1
        m = months[month_idx]

    elif (month is not None) and (year is not None):
        print(f"Retrieving EIA Form 860m for {month.capitalize()} {year}")
        m = month
        y = year

    elif ((month is None) or (year is None)):
        
        print(f"Month {month} / Year {year}")
        raise ValueError(("Please specify a month and a year."))


    url = (f"https://www.eia.gov/electricity/data/eia860m/archive/xls/" +
           f"{m}_generator{y}.xlsx")

    try:
        print(f'Downloading from {url}\n')
        df = pd.read_excel(url,
                           sheet_name='Operating',
                           skipfooter=2,
                           skiprows=2,
                           usecols=columns,
                           index_col='Entity ID')
        print('Download successful.')
    except BaseException:
        print('Download failed. Trying different sheet format.')
        try:
            df = pd.read_excel(url,
                               sheet_name='Operating',
                               skipfooter=2,
                               skiprows=1,
                               usecols=columns,
                               index_col='Entity ID')
            print('Download successful.')
        except ValueError:
            fail_str = (f'Download failed. File not found' +
                       f' for Month: {month} and Year: {year}')
            raise ValueError(fail_str)

    return df


def get_region_techs(df, region):
    """
    Gets the existing capacity for a particular technology and region.
    Currently, only electric capacity is considered.

    If the user has an internet connection, get_region_techs will collect
    the most up-to-date electric generator data. 

    Users can also specify a month and year to obtain a specific historical
    dataset from EIA.

    Parameters
    ----------
    df : pandas dataframe
        The dataframe for EIA form 860M
    region : string
        The region of interest. Region may be state or county. The state must
        be given as a two letter abbreviation and a county must be provided as 
        a full name.
    """
        # filter by region
    if len(region) == 2:
        valid_state = (region.upper() in df['Plant State'].values)
        if valid_state:
            region_mask = df['Plant State'] == region.upper()
        else:
            raise ValueError(f'Detected state abbreviation. Abbreviation {region} not found.')
    else:
        valid_county = (region.capitalize() in df['County'].values)
        if valid_county:
            region_mask = df['County'] == region.capitalize()
        else:
            raise ValueError(f'Detected county name. County name {region} not found.')

    df = df[region_mask]
    
    return df