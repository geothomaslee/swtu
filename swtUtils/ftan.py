#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:29:49 2024

@author: Thomas Lee
University of New Mexico
Department of Earth and Planetary Science

A package for preparing ambient noise cross-correlation functions for FTAN analysis
as implemented by Bensen et. al.

Citations:
    Bensen, G. D., et al. (2007). Processing seismic ambient noise data to
    obtain reliable broad-band surface wave dispersion measurements.
        Geophysical Journal International, 169(3), 1239-1260.
"""
#pylint: disable=invalid-name

import os
from glob import glob
from math import floor
from typing import Tuple

import obspy
from obspy import Trace
from tqdm import tqdm
import pandas as pd
from pandas import DataFrame

def getStackDirectory(dataDirectory: str) -> str:
    """Returns the stack directory given the data directoy"""
    return os.path.expanduser(dataDirectory + '/Stacks')

def getComponentDirectory(dataDirectory: str, component: str) -> str:
    """Returns the component directory given the data directory"""
    if not isinstance(component, str):
        raise TypeError('Component must be 2 character string. Ex. "ZZ"')
    if len(component) != 2:
        raise ValueError('Component must be 2 character string. Ex. "ZZ"')

    stackDirectory = getStackDirectory(dataDirectory)

    return stackDirectory + f'/{component}'

def getSACFileList(componentDirectory: str) -> list:
    """Returns a list of SAC files in the component directory"""
    fileList = glob(componentDirectory + '/*.sac')
    return fileList

def getFoldDirectory(dataDirectory: str, component: str) -> str:
    """Returns the directory containing folded traces"""
    componentDirectory = getComponentDirectory(dataDirectory,component)
    foldDirectory = componentDirectory + '/Folded'

    if not os.path.isdir(foldDirectory):
        os.mkdir(foldDirectory)

    return foldDirectory

def getSACDict(tr: Trace) -> dict:
    """SAC Dict properly formatted for applying headers after folding"""
    sacDict = tr.stats.sac

    del sacDict['npts']
    del sacDict['e']
    del sacDict['b']
    del sacDict['depmin']
    del sacDict['depmax']

    return sacDict

def foldTrace(tr: Trace,delta=0.2) -> Trace:
    """Returns a folded version of the trace"""
    if (tr.stats.npts %2) != 1:
        raise ValueError('Trace must have an odd number of points')

    middle_index = floor(tr.stats.npts / 2)

    left = tr.data[0:middle_index]
    reversed_left = left[::-1]
    right = tr.data[middle_index+1:]

    data = (reversed_left + right) / 2
    new_tr = Trace()
    new_tr.data = data

    new_tr.stats.sac = getSACDict(tr)
    new_tr.stats.delta = delta

    return new_tr

def foldAllTraces(dataDirectory: str,component: str):
    """Folds all traces and puts them in separate folded directory"""
    foldDirectory = getFoldDirectory(dataDirectory,component)
    componentDirectory = getComponentDirectory(dataDirectory,component)

    fileList = getSACFileList(componentDirectory)

    print('WARNING: Must manually define delta in foldTrace. This was the only ' +
          'workaround I could find for a really annoying floating point error')

    for file in tqdm(fileList):
        basename= os.path.basename(file)[0:-4]
        tr = obspy.read(file)[0]
        folded_tr = foldTrace(tr)
        savePath = foldDirectory + f'/{basename}_Folded.sac'
        folded_tr.write(savePath)

def dispOutputToDF(file: str) -> DataFrame:
    """Turns disp file into df"""
    columns=['nf','cper','obper','gvel','phvel','ampl','snr']
    df = pd.read_csv(file,sep='\s+',names=columns)
    return df

def findDispFiles(dispDirectory: str) -> list:
    """Finds all final FTAN output files"""
    fileList = glob(dispDirectory + '/*2_DISP.1*')
    return fileList

def getRelevantInfo(df: DataFrame) -> Tuple[list,list,list]:
    """Returns observed period, phase velocity, and snr info"""
    obper = df['obper'].to_list()
    phvel = df['phvel'].to_list()
    snr = df['snr'].to_list()
    return obper,phvel,snr

def getPeriodRange(df: DataFrame) -> Tuple[float, float]:
    """Returns the range of periods measured by FTAN"""
    vals = df['obper'].to_list()
    return min(vals),max(vals)

def returnSuccessRate(ftanDirectory: str) -> dict:
    """Reports success rate of FTAN after SNR filtering"""
    successDict = {'raw' : 0,
                   'snr' : 0,
                   '1_DISP' : 0,
                   '2_DISP' : 0}

    dirpath = f'{ftanDirectory}/Folded'

    files = glob(dirpath + '/*.sac*')

    for file in files:
        if 'snr' in file:
            successDict['snr'] += 1
        elif '1_DISP' in file:
            successDict['1_DISP'] += 1
        elif '2_DISP' in file:
            successDict['2_DISP'] += 1
        else:
            successDict['raw'] += 1

    return successDict