#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 07:55:48 2024

@author: Thomas Lee
University of New Mexico
Department of Earth and Planetary Science

A module for creating properly formatted inputs for Fast-Marching Surface
Tomography (FMST) by Nick Rawlinson using outputs from pyaftan by Lili Feng.
Assumes the stacked cross-correlation functions are in the format and file
structure outputted by ambient2.
"""

import os
import pickle
from glob import glob
from tqdm import tqdm
import shutil

import numpy as np
import obspy
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt

import swtUtils.ftan as ft


def saveObj(obj, filename):
    """Quick function for pickling a file"""
    with open(filename, 'wb') as f:
        pickle.dump(obj, f)

def loadObj(filename):
    """Quick function for pickling a file"""
    with open(filename, 'rb') as f:
        return pickle.load(f)

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

def getTomoDirectory(dataDirectory,component):
    if not isinstance(component, str):
        raise TypeError('Component must be 2 character string. Ex. "ZZ"')
    if len(component) != 2:
        raise ValueError('Component must be 2 character string. Ex. "ZZ"')

    return os.path.expanduser(dataDirectory + f'/Tomography/{component}')

def getStationNames(filepath):
    """Pulls the name of the two stations from the name of the file"""
    basename = os.path.basename(filepath)
    name_no_sac = basename[:-4]
    stat1, stat2 = name_no_sac.split('_')

    return stat1, stat2

def getLocalStations(dataDirectory,component):
    """Returns a list of stations that we have cross-correlations for"""
    componentDirectory = getComponentDirectory(dataDirectory,component)

    if not os.path.exists(componentDirectory +'/UniqueStations.pkl'):
        crossCorrelations = glob(componentDirectory + '/*.sac')
        station_list = []
        for file in crossCorrelations:
            stat1, stat2 = getStationNames(file)

            station_list = [s for s in set(station_list + [stat1, stat2]) if s]
            saveObj(station_list,componentDirectory +'/UniqueStations.pkl')
    else:
        station_list = loadObj(componentDirectory +'/UniqueStations.pkl')

    return station_list

def getValidStations(network,bounds,channel,stationList):
    client = Client("IRIS")
    inventory = client.get_stations(network = network,
                                    station = '*',
                                    channel = channel,
                                    minlatitude = bounds[0],
                                    maxlatitude = bounds[1],
                                    minlongitude = bounds[2],
                                    maxlongitude = bounds[3])

    inventoryCopy = inventory.copy()

    for network in inventory:
        for station in network:
            formattedName = f'{network.code}.{station.code}'
            if formattedName not in stationList:
                inventoryCopy = inventoryCopy.remove(network=network.code,station=station.code)

    stationDict = {}
    for network in inventoryCopy:
        for station in network:
            stationDict[f'{network.code}.{station.code}'] = (station.latitude,station.longitude)

    return stationDict

def getInventoryLength(inventory):
    length = 0
    for network in inventory:
        length += len(network)

    return(length)

def makeTomoDirectory(dataDirectory,periods,component):
    """Periods must be a list of all periods, or if only 2 elements, is assumed
    to want every integer period between element 1 and element 2.

    Ex. [1, 5, 7, 8] would use those periods, [1, 8] would use every integer from
    1 to 8 as the periods.
    """
    if not isinstance(periods, list):
        raise TypeError('Periods must be a list, see doc for details')
    if len(periods) < 2:
        raise ValueError('Periods must be a list of length >=2, see docs')
    if len(periods) == 2:
        periods = range(periods[0],periods[1]+1,1)

    tomoDirectory = os.path.expanduser(dataDirectory) + '/Tomography'
    if not os.path.isdir(tomoDirectory):
        os.mkdir(tomoDirectory)

    componentDirectory = tomoDirectory + f'/{component}'
    if not os.path.isdir(componentDirectory):
        os.mkdir(componentDirectory)

    for period in periods:
        periodFolder = componentDirectory + f'/{period}s'
        if not os.path.isdir(periodFolder):
            os.mkdir(periodFolder)

    return tomoDirectory

def makeIssueDict():
    """Makes the empty issue dict"""
    issue_dict = {'autocorrelation' : 0,
                  'filepath not exist' : 0,
                  'fparam_out_none' : 0,
                  'interpOut none' : 0,
                  'bad phvel' : 0,
                  'low snr' : 0,
                  'too close' : 0,
                  'good' : 0}
    return issue_dict

def makefpDict():
    """Makes error dict for fparam"""
    fpDict = {'infDepen' : 0,
              'omdom' : 0,
              'zero' : 0,
              'slice' : 0,
              'fparam returning none' : 0}

    return fpDict

def makeInterpErrorDict():
    """Makes an error dictionary for interpolating the period"""
    interpErrorDict = {'lowest observed period greater than desired' : 0,
                       'highest observed period lower than desired' : 0,
                       'could not find left or right' : 0}
    return interpErrorDict

def makeFMSTInputs(stationDict,dataDirectory,FTANDirectory,period,component,minSNR,minWavelengths,detailedError=True):
    # Get directories
    tomoDirectory = getTomoDirectory(dataDirectory,component)
    periodDirectory = tomoDirectory + f'/{period}s'
    componentDirectory = getComponentDirectory(dataDirectory,component)

    # Get filepaths for output files
    receiverFile = periodDirectory +'/receivers.dat'
    sourcesFile = periodDirectory +'/sources.dat'
    timesFile = periodDirectory + '/otimes.dat'

    # Creating issue dicts
    issue_dict = makeIssueDict()
    fpDict = makefpDict()
    interpErrorDict = makeInterpErrorDict()

    # Create receiver file
    open(receiverFile,'w',encoding='utf-8').close()
    with open(receiverFile, 'a',encoding='utf-8') as file:
        file.write(f'{int(len(list(stationDict.keys())))}\n')
        for station, coords in stationDict.items():
            file.write(f'{coords[0]} {coords[1]}\n')
        file.close()

    # Create source file
    open(sourcesFile,'w',encoding='utf-8').close()
    with open(sourcesFile, 'a',encoding='utf-8') as file:
        file.write(f'{int(len(list(stationDict.keys())))}\n')
        for station, coords in stationDict.items():
            file.write(f'{coords[0]} {coords[1]}\n')
        file.close()

    stationList = list(stationDict.keys())

    open(timesFile,'w',encoding='utf-8').close()

    phvels = []
    with open(timesFile, 'a',encoding='utf-8') as outfile:
        #outfile.write(f'{len(stationList)**2}\n')
        for i, stat1 in tqdm(enumerate(stationList),total=len(stationList)):
            for j, stat2 in enumerate(stationList):
                if stat1 == stat2:
                    issue_dict['autocorrelation'] += 1
                    outfile.write('0 0.0000 1.0\n')
                    continue

                filepath = checkIfFTANExists(stat1,stat2,FTANDirectory)
                if filepath is None:
                    issue_dict['filepath not exist'] += 1
                    outfile.write('0 0.0000 1.0\n')
                    continue

                dist = getDist(stat1,stat2,componentDirectory)
                if dist is None:
                    issue_dict['filepath not exist'] += 1
                    outfile.write('0 0.0000 1.0\n')
                    continue

                df = ft.dispOutputToDF(filepath)
                obper,phvel,snr = ft.getRelevantInfo(df)

                interpOut = interpPeriod(period,obper,phvel,snr)
                interpOut, interpErrorDict = _interpPeriodErrorHandler(interpOut,interpErrorDict)
                if interpOut is None:
                    issue_dict['interpOut none'] += 1
                    outfile.write('0 0.0000 1.0\n')
                    continue

                phvel = interpOut[0]
                snr = interpOut[1]

                if phvel < 1.5 or phvel > 5:
                    outfile.write('0 0.0000 1.0\n')
                    issue_dict['bad phvel'] += 1
                    continue

                if snr < minSNR:
                    issue_dict['low snr'] += 1
                    outfile.write('0 0.0000 1.0\n')
                    continue

                minDist = minWavelengths * getWavelength(period,phvel)
                if dist < minDist:
                    issue_dict['too close'] += 1
                    outfile.write('0 0.0000 1.0\n')
                    continue

                phvels.append(phvel)

                travelTime = round(float(getTravelTime(phvel,dist)), 4)
                issue_dict['good'] += 1
                outfile.write(f'1 {travelTime} 1.0\n')

    if detailedError is True:
        saveObj(issue_dict, f'{periodDirectory}/issueDict.pkl')
        saveObj(fpDict, f'{periodDirectory}/fpDict.pkl')
        saveObj(interpErrorDict,f'{periodDirectory}/interpErrorDict.pkl')

    return phvels

def _fparamErrorCodeHandler(fparam_out,fpDict,outfile):
    """Counts the types of errors from fparam"""
    if isinstance(fparam_out,tuple):
        if fparam_out[0] is None:
            fpDict['fparam returning none'] += 1
            return None, fpDict
        return fparam_out, fpDict

    if fparam_out == 1:
        fpDict['infDepen'] += 1
        outfile.write('0 0.0000 1.0\n')
    if fparam_out == 2:
        fpDict['omdom'] += 1
        outfile.write('0 0.0000 1.0\n')
    if fparam_out == 3:
        fpDict['zero'] += 1
        outfile.write('0 0.0000 1.0\n')
    if fparam_out == 4:
        fpDict['slice'] += 1
        outfile.write('0 0.0000 1.0\n')
    if fparam_out == 5:
        fpDict['cannot unpack none'] += 1
        outfile.write('0 0.0000 1.0\n')

    return None, fpDict

def interpPeriod(period,obper,phvel,snr):
    """Interpolates the phase velocity and snr to the period of interest"""
    if obper[0] > period:
        return 1

    if obper[-1] < period:
        return 2

    left = None
    right = None
    for i, per in enumerate(obper):
        if per < period:
            left = i
        if per > period:
            right = i
            break

    if left is None or right is None:
        return 3

    interpPhvel = np.interp(period,[obper[left],obper[right]],[phvel[left],phvel[right]])
    interpSNR = np.interp(period,[obper[left],obper[right]],[snr[left],snr[right]])

    return interpPhvel,interpSNR

def _interpPeriodErrorHandler(interpPeriodOut,interpPeriodDict):
    """Wrapper for handling errors from interpPeriod"""
    if isinstance(interpPeriodOut, tuple):
        return interpPeriodOut, interpPeriodDict

    if interpPeriodOut == 1:
        interpPeriodDict['lowest observed period greater than desired'] += 1
    if interpPeriodOut == 2:
        interpPeriodDict['highest observed period lower than desired'] += 1
    if interpPeriodOut == 3:
        interpPeriodDict['could not find left or right'] += 1

    return None, interpPeriodDict

def checkIfFTANExists(stat1,stat2,FTANDirectory):
    """Checks if an FTAN output from AFTAN (Bensen) exists"""
    if os.path.exists(FTANDirectory + f'/{stat1}_{stat2}_Folded.sac_2_DISP.1'):
        return FTANDirectory + f'/{stat1}_{stat2}_Folded.sac_2_DISP.1'

    if os.path.exists(FTANDirectory + f'/{stat2}_{stat1}_Folded.sac_2_DISP.1'):
        return FTANDirectory + f'/{stat2}_{stat1}_Folded.sac_2_DISP.1'

    return None

def getDist(stat1,stat2,componentDirectory):
    if os.path.exists(componentDirectory + f'/{stat1}_{stat2}.sac'):
        return obspy.read(componentDirectory + f'/{stat1}_{stat2}.sac')[0].stats.sac['dist']

    if os.path.exists(componentDirectory + f'/{stat2}_{stat1}.sac'):
        return obspy.read(componentDirectory + f'/{stat2}_{stat1}.sac')[0].stats.sac['dist']

    return None

def getTravelTime(vel, dist):
    """Returns travel time given velocity and distance"""
    return dist / vel

def getWavelength(per,vel):
    """Returns the wavelength given period and velocity"""
    return per * vel

def plotIssueDict(issue_dict, label):
    """Plots one of the issue dicts as a histogram"""
    fig, ax = plt.subplots()

    issues = list(issue_dict.keys())
    count = list(issue_dict.values())

    ax.bar(issues, count)

    plt.xticks(rotation=65,ha='right')
    plt.title(label)
    plt.show()

def setupFTANDirectory(FMSTDirectory,period,projectCode,component,_overwrite=False):
    dirName = f'{projectCode}_{period}s_{component}'
    fmstPath = FMSTDirectory + f'/{dirName}'
    fmstMasterPath = FMSTDirectory + f'/{projectCode}_Master'

    if not os.path.isdir(FMSTDirectory):
        raise ValueError('Could not find ftan directory')

    if not isinstance(projectCode,str):
        raise TypeError('Project code must be str. See doc')

    if not isinstance(component, str):
        raise TypeError('Component must be str. See doc')
    if len(component) != 2:
        raise ValueError('Component must be two char str. Ex. "ZZ","EE"')

    if os.path.isdir(fmstPath) is True and _overwrite is True:
        shutil.rmtree(fmstPath)
    elif os.path.isdir(fmstPath) is True and _overwrite is False:
        print(f'FMST path for {projectCode} {component} {period}s already exists')
        return fmstPath

    if not os.path.isdir(fmstMasterPath):
        raise ValueError('Could not find master path for project')

    shutil.copytree(src=fmstMasterPath,
                    dst=fmstPath)

    return fmstPath
