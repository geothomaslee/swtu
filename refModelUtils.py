#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 10:21:58 2024

@author: thomaslee
"""

import pandas as pd
import numpy as np




ref_file = '/Users/thomaslee/Documents/GitHub/swtu/resources/R05_USANT15.txt'

columns = ['Station_1','Station_2','Station_1_Latitude',
           'Station_1_Longitude','Station_2_Latitude',
           'Station_2_Longitude','Distance','Number_of_Correlations',
           'Period','Wave_Type','Phase_Velocity','Phase_Velocity_Uncertainty']


for period in ['05','06','08','10','12','15','20','25','30','35','40']:
    ref_file = f'/Users/thomaslee/Documents/GitHub/swtu/resources/R{period}_USANT15.txt'
    print(ref_file)

    bounds = [46.0796,47.8242,-122.8625,-120.25]


    df = pd.read_csv(ref_file,names=columns,delimiter='\s+')

    print(f'TOTAL PATHS: {len(df)}')

    df = df.drop(df[(df.Station_1_Latitude > bounds[1])].index)
    df = df.drop(df[(df.Station_1_Latitude < bounds[0])].index)
    df = df.drop(df[(df.Station_2_Latitude > bounds[1])].index)
    df = df.drop(df[(df.Station_2_Latitude < bounds[0])].index)

    df = df.drop(df[(df.Station_1_Longitude > bounds[3])].index)
    df = df.drop(df[(df.Station_1_Longitude < bounds[2])].index)
    df = df.drop(df[(df.Station_2_Longitude > bounds[3])].index)
    df = df.drop(df[(df.Station_2_Longitude < bounds[2])].index)

    phvels = df['Phase_Velocity'].to_list()
    if period == '08':
        print(df['Station_1'])
        print(df['Station_2'])
        print(df['Phase_Velocity'])
    print(phvels)
    avgphvel = np.mean(phvels)
    print(f'Average Phase Velocity for {period}s: {avgphvel}')

