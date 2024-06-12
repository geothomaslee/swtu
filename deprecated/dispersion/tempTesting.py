#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 08:05:31 2024

@author: Thomas Lee
University of New Mexico
Department of Earth and Planetary Science
"""

from glob import glob
from math import floor
import sys
from os.path import basename
import obspy
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

sys.path.append('/users/thomaslee/Documents/pyaftan-master')
import pyaftan #pylint: disable=import-error

testFile = '/Volumes/NewHDant/RainierAmbient/Stacks/ZZ/UW.GSM_XD.MB05.sac'

tr = obspy.read(testFile)[0]

def updateSacHeaders(tr):
    """
    Updates the SAC headers of each trace in the given stream based on the trace stats.
    Parameters:
    - stream: ObsPy Stream object containing one or more traces.
    Returns:
    - The modified stream with updated SAC headers for each trace.
    """

    sachdr = {
        'b': (tr.stats.starttime - tr.stats.endtime)/2,
        'e': (tr.stats.endtime - tr.stats.starttime)/2,
        'stla': tr.stats.sac['stla'],
        'stlo': tr.stats.sac['stlo'],
        'stel': tr.stats.sac['stel'],
        'evla': tr.stats.sac['evla'],
        'evlo': tr.stats.sac['evlo'],
        'kstnm': tr.stats.station,
        'knetwk': tr.stats.network,
        'kcmpnm': tr.stats.channel,
        'delta': tr.stats.delta,
        'dist': tr.stats.sac['dist'],
        'npts': tr.stats.npts,
        'nzyear': tr.stats.starttime.year,
        'nzjday': tr.stats.starttime.julday,
        'nzhour': tr.stats.starttime.hour,
        'nzmin': tr.stats.starttime.minute,
        'nzsec': tr.stats.starttime.second,
        'nzmsec': int(tr.stats.starttime.microsecond / 1000),
    }

    tr.stats.sac = sachdr

    return tr

"""
def plotftan(fparam, dist):
    v1      = dist/(fparam.tamp_1+np.arange(fparam.ncol_1)*dt)
    ampo_1  = fparam.ampo_1[:fparam.ncol_1,:fparam.nrow_1]
    obper1_1= fparam.arr1_1[1,:fparam.nfout1_1]
    gvel1_1 = fparam.arr1_1[2,:fparam.nfout1_1]
    phvel1_1= fparam.arr1_1[3,:fparam.nfout1_1]
    plt.figure()
    ax      = plt.subplot()
    p       = plt.pcolormesh(obper1_1, v1, ampo_1, cmap=cmap,shading='gouraud')
    ax.plot(obper1_1, gvel1_1, '--k', lw=3) #
    if fparam.preflag:
        ax.plot(obper1_1, phvel1_1, '--w', lw=3) #
    if (fparam.nfout2_1!=0):
        obper2_1    = fparam.arr2_1[1,:fparam.nfout2_1]
        gvel2_1     = fparam.arr2_1[2,:fparam.nfout2_1]
        phvel2_1    = fparam.arr2_1[3,:fparam.nfout2_1]
        ax.plot(obper2_1, gvel2_1, '-k', lw=3) #
        if fparam.preflag:
            ax.plot(obper2_1, phvel2_1, '-w', lw=3) #
    cb      = plt.colorbar(p, ax=ax)
    Tmin1   = obper1_1[0]
    Tmax1   = obper1_1[fparam.nfout1_1-1]
    vmin1   = v1[fparam.ncol_1-1]
    vmax1   = v1[0]
    plt.axis([Tmin1, Tmax1, vmin1, vmax1])
    plt.xlabel('Period(s)')
    plt.ylabel('Velocity(km/s)')
    plt.title('Basic FTAN Diagram '+sacname,fontsize=15)
"""




tr = updateSacHeaders(tr)

atr1=pyaftan.aftantrace(tr.data.astype(np.float32),tr.stats)
atr1.makesym()
# aftan analysis using pyaftan
atr1.aftanf77(tmin=1, tmax=25, vmin=1.5,vmax=5,snr=5,tresh=20,phvelname='ak135.disp')
atr1.plotftan(plotflag=1)

arr2_2 = atr1.ftanparam.arr2_2

cp = arr2_2[0,:] #-  central periods, s (real*8)
obper = arr2_2[1,:] #-  apparent periods, s (real*8)
gvel = arr2_2[2,:] #-  group velocities, km/s (real*8)
phvel = arr2_2[3,:] #-  phase velocities, km/s (real*8)
ampdb = arr2_2[4,:] #-  amplitudes, Db (real*8)
snr_db = arr2_2[5,:] #-  signal/noise ratio, Db (real*8)
maxhw = arr2_2[6,:] #-  maximum half width, s (real*8)
amp = arr2_2[7,:] #-  amplitudes
#snr = arr2_2[8,:] #-  signal/noise ratio (optional)

print(obper)
print(phvel)

plt.plot(obper, gvel)
plt.title('gvel')
plt.show()

plt.plot(obper, phvel)
plt.title('phvel')
plt.show()

plt.plot(obper, ampdb)
plt.title('ampdb')
plt.show()

plt.plot(obper, snr_db)
plt.title('snr db')
plt.show()

plt.plot(obper,maxhw)
plt.title('maxhw')
plt.show()

plt.plot(obper,amp)
plt.title('amp')
plt.show()

"""
plt.plot(obper,snr)
plt.show()
"""