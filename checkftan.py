#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:58:33 2025

@author: Thomas Lee
Rice University
"""
import os

import pandas as pd
import matplotlib.pyplot as plt
import obspy

pd.set_option('display.max_columns', None)

ftanDirectory = '/Users/thomaslee/FTAN'
dataDirectory = '/Volumes/NewHDant/RainierAmbient2'
snr = 5
min_wavelengths = 2
stat1 = 'UW.HSR'
stat2 = 'CC.ARAT'

df = pd.read_csv(f'{ftanDirectory}/PhaseVelocities_SNR{snr}_WL{min_wavelengths}.csv')

mask = (df['Station1'] == stat1) & (df['Station2'] == stat2)
if df[mask].empty is True:
    mask = (df['Station2'] == stat1) & (df['Station1'] == stat2)
    if df[mask].empty is True:
        raise ValueError('Could not find measurements for this station pair')

result = df[mask]

columns = result.columns.tolist()
columns = columns[7:]

periods = [8.0]
for period in periods:
    vel = result[f'{period}s'].tolist()[0]

    """
    if os.path.isfile(f'{ftanDirectory}/Folded_Clean/{stat1}_{stat2}_Folded.sac'):
        tr = obspy.read(f'{ftanDirectory}/Folded_Clean/{stat1}_{stat2}_Folded.sac')[0]
    else:
        tr = obspy.read(f'{ftanDirectory}/Folded_Clean/{stat2}_{stat1}_Folded.sac')[0]
    """

    if os.path.isfile(f'{dataDirectory}/Stacks/ZZ/{stat1}_{stat2}.sac'):
        tr = obspy.read(f'{dataDirectory}/Stacks/ZZ/{stat1}_{stat2}.sac')[0]
    else:
        tr = obspy.read(f'{dataDirectory}/Stacks/ZZ/{stat2}_{stat1}.sac')[0]

    pulse_width = 1


    from swtUtils.ftan import foldTrace

    folded_tr = foldTrace(tr)
    print(tr.stats)
    print(folded_tr.stats)

    print(folded_tr.stats.sampling_rate)


    tr.filter('bandpass',freqmin=1/(period+pulse_width),freqmax=1/(period-pulse_width))
    stats = tr.stats.sac
    expected_travel_time = stats['dist'] / vel
    times = [x - 300 for x in tr.times()]
    print(tr.stats.delta)

    heights = [min(tr.data),max(tr.data)]
    travel_times_pos = [expected_travel_time,expected_travel_time]
    travel_times_neg = [-expected_travel_time,-expected_travel_time]

    fig, ax = plt.subplots()
    ax.plot(times,tr.data)
    ax.plot(travel_times_pos,heights,'-')
    ax.plot(travel_times_neg,heights,'-')
    ax.set_xlabel('Lag Time(s)')
    ax.set_ylabel('Correlation')
    plt.show()





