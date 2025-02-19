#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:47:16 2025

@author: Thomas Lee
Rice University
"""

#periods = [5.0,6.0,7.0,8.0,9.0,10.0]
periods = [5.0]

count_dict = {}
ch_count_dict = {}
for period in periods:
    file = f'/Users/thomaslee/fmst_v1.1/Rainier_{period}s_ZZ/otimes.dat'
    chfile = f'/Users/thomaslee/fmst_v1.1/Checkerboard_Rainier_{period}s_ZZ/otimes.dat'

    count = 0
    with open(file, 'r') as otimes:
        lines = otimes.readlines()
        for line in lines:
            if str(line[0]) == '1':
                count += 1

    count_dict[str(period)] = count

    count = 0
    with open(chfile, 'r') as otimes:
        lines = otimes.readlines()
        for line in lines:
            if str(line[0]) == '1':
                count += 1

    ch_count_dict[str(period)] = count

print('Original:')
print(count_dict)
print('Checkerboard:')
print(ch_count_dict)
