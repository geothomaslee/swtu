#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:02:10 2024

@author: Thomas Lee
University of New Mexico
Department of Earth and Planetary Science
"""

import sys
import os
from math import floor
from glob import glob

import numpy as np
import obspy
import matplotlib.pyplot as plt
sys.path.append('/users/thomaslee/Documents/pyaftan-master')
import pyaftan #pylint: disable=import-error

class infDependentError(Exception):
    """Custom class for catching inf values from  aftantrace._aftanipg"""

class omdomError(Exception):
    """Custom class for catching an erorr in _ftanipg that I don't get"""

class sliceError(Exception):
    """Custom class for catching slice indices error"""

def _getTimes(tr):
    """Returns times from -lag to positive lag where lag = 0.5*signal length"""
    npts = tr.stats.npts
    delta = tr.stats.delta
    if (npts % 2) == 0:
        lag = (npts / 2) * delta
    if (npts % 2) == 1:
        lag = (floor(npts / 2)) * delta

    times = np.arange(-lag,lag+delta,delta)
    return times

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

def getStationNames(filepath):
    """Pulls the name of the two stations from the name of the file"""
    basename = os.path.basename(filepath)
    name_no_sac = basename[:-4]
    stat1, stat2 = name_no_sac.split('_')

    return stat1, stat2

def _FTANHelper(filepath,tmin,tmax,nfin):
    """Runs ftan and returns an ftanTrace.fparam object"""
    try:
        tr = obspy.read(filepath)[0]
        tr = updateSacHeaders(tr)
        atr1=pyaftan.aftantrace(tr.data,tr.stats)
        atr1.makesym()
        # aftan analysis using pyaftan
        atr1.aftan(tmin=tmin, tmax=tmax, phvelname='ak135.disp',nfin=nfin)
        fparam = atr1.ftanparam

        return fparam, tr.stats

    except ValueError as e:
        if str(e) == "`y` must contain only finite values.":
            raise infDependentError from e
    except IndexError as e:
        if str(e) == "index 0 is out of bounds for axis 0 with size 0":
            raise omdomError from e
    except TypeError as e:
        if str(e) == "slice indices must be integers or None or have an __index__ method":
            raise sliceError from e

    return None, None

def FTAN(filepath,tmin,tmax,nfin):
    """Wrapper for _FTAN that does a better job handling errors"""
    try:
        fparam, stats = _FTANHelper(filepath,tmin=tmin,tmax=tmax,nfin=nfin)
        if fparam == None:
            return None, None
        return fparam, stats
    except infDependentError:
        return 1
    except omdomError:
        return 2
    except ZeroDivisionError:
        return 3
    except sliceError:
        return 4
    except TypeError as e:
        # This happens when FTAN doesn't work and returns None, we want to
        # just skip through this iteration
        if str(e) == 'cannot unpack non-iterable NoneType object':
            return 5
    except Exception as e:
        print(f'OTHER EXCEPTION RESULTING IN NONE?\n {e}')
        return None, None

def pullRelevantFTANInfo(fparam):
    """Pulls the info from fparam that's relevant for inversions"""
    obper1  = fparam.arr2_2[0,:fparam.nfout2_2]
    gvel1   = fparam.arr2_2[2,:fparam.nfout2_2]
    phvel1  = fparam.arr2_2[3,:fparam.nfout2_2]
    snr = fparam.arr2_2[5,:fparam.nfout2_2]

    return obper1, gvel1, phvel1, snr

def checkAutoCorrelation(filepath):
    """Returns True if it's an autocorrelation"""
    stat1, stat2 = getStationNames(filepath)

    if stat1 == stat2:
        return True
    return False

def plotFTAN(filepath,tmin,tmax):
    """Plots FTAN analysis"""
    tr = obspy.read(filepath)[0]
    times = _getTimes(tr)
    data = tr.data
    stat1, stat2 = getStationNames(filepath)
    fparam = FTAN(filepath,tmin,tmax)

    fparamresult = pullRelevantFTANInfo(fparam)

    obper1 = fparamresult[0]
    gvel1 = fparamresult[1]
    phvel1 = fparamresult[2]

    fig, ax = plt.subplots(2,1)
    ax1 = ax[0]
    ax1.plot(times,data)
    ax1.set_title(f"Cross-Correlation Function for {stat1} to {stat2}")
    ax1.set_xlabel('Lag (s)')

    ax2 = ax[1]
    ax2.plot(obper1,gvel1,label='Group Velocity')
    ax2.plot(obper1,phvel1,label='Phase Velocity')
    ax2.set_ylabel('Velocity (km/s)')
    ax2.set_xlabel('Period (s)')
    ax2.set_title('FTAN Analysis with Phase-Matched Filtering')
    fig.tight_layout()
    plt.legend()
    plt.show()

def getFileList(componentDirectory):
    """Returns a list of all files in the component directory"""
    fileList = glob(os.path.expanduser(componentDirectory) + '/*')

    return fileList
