"""
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

    # initialize 
    m = 'thermidor'
    y = 0

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
        print("Please specify a month and a year.")

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
        except FileNotFoundError:
            print('Download failed. File not found.')

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
        be given as an abbreviation a county must be provided as a full name.
    """
        # filter by region
    if len(region) == 2:
        try:
            region_mask = df['Plant State'] == region.upper()
        except BaseException:
            print(f'Detected state abbreviation. Abbreviation {region} not found.')
    else:
        try:
            region_mask = df['County'] == region.upper()
        except BaseException:
            print(f'Detected county name. County name {region} not found.')
    df = df[region_mask]
    
    return df