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

import os
from glob import glob
from math import floor

import obspy
from obspy import Trace
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def getStackDirectory(dataDirectory):
    """Returns the stack directory given the data directoy"""
    return os.path.expanduser(dataDirectory + '/Stacks')

def getComponentDirectory(dataDirectory, component):
    """Returns the component directory given the data directory"""
    if not isinstance(component, str):
        raise TypeError('Component must be 2 character string. Ex. "ZZ"')
    if len(component) != 2:
        raise ValueError('Component must be 2 character string. Ex. "ZZ"')

    stackDirectory = getStackDirectory(dataDirectory)

    return stackDirectory + f'/{component}'

def getSACFileList(componentDirectory):
    """Returns a list of SAC files in the component directory"""
    fileList = glob(componentDirectory + '/*.sac')
    return fileList

def getFoldDirectory(dataDirectory, component):
    componentDirectory = getComponentDirectory(dataDirectory,component)
    foldDirectory = componentDirectory + '/Folded'

    if not os.path.isdir(foldDirectory):
        os.mkdir(foldDirectory)

    return foldDirectory

def getSACDict(tr):
    """SAC Dict properly formatted for applying headers after folding"""
    sacDict = tr.stats.sac

    del sacDict['npts']
    del sacDict['e']
    del sacDict['b']
    del sacDict['depmin']
    del sacDict['depmax']

    return sacDict


def foldTrace(tr):
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

    return new_tr

def foldAllTraces(dataDirectory,component):
    """Folds all traces and puts them in separate folded directory"""
    foldDirectory = getFoldDirectory(dataDirectory,component)
    componentDirectory = getComponentDirectory(dataDirectory,component)

    fileList = getSACFileList(componentDirectory)

    for file in tqdm(fileList):
        basename= os.path.basename(file)[0:-4]
        tr = obspy.read(file)[0]
        folded_tr = foldTrace(tr)
        savePath = foldDirectory + f'/{basename}_Folded.sac'
        folded_tr.write(savePath)


def dispOutputToDF(file):
    columns=['nf','cper','obper','gvel','phvel','ampl','snr']
    df = pd.read_csv(file,delim_whitespace=True,names=columns)
    return df

def findDispFiles(dispDirectory):
    fileList = glob(dispDirectory + '/*2_DISP.1*')
    return fileList

def getRelevantInfo(df):
    obper = df['obper'].to_list()
    phvel = df['phvel'].to_list()
    snr = df['snr'].to_list()
    return obper,phvel,snr

def getPeriodRange(df):
    vals = df['obper'].to_list()
    return min(vals),max(vals)
