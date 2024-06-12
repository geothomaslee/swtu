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

import numpy as np
import pyaftan_tools as pt
import obspy
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt


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

def getLocalStations(dataDirectory,component):
    """Returns a list of stations that we have cross-correlations for"""
    componentDirectory = getComponentDirectory(dataDirectory,component)

    if not os.path.exists(componentDirectory +'/UniqueStations.pkl'):
        crossCorrelations = glob(componentDirectory + '/*.sac')
        station_list = []
        for file in crossCorrelations:
            stat1, stat2 = pt.getStationNames(file)

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

def makeFMSTInputs(stationDict,dataDirectory,period,component,minSNR,minWavelengths):
    tomoDirectory = getTomoDirectory(dataDirectory,component)
    componentDirectory = getComponentDirectory(dataDirectory,component)
    periodDirectory = tomoDirectory + f'/{period}s'

    receiverFile = periodDirectory +'/receivers.dat'
    sourcesFile = periodDirectory +'/sources.dat'
    timesFile = periodDirectory + '/otimes.dat'

    open(receiverFile,'w').close()
    with open(receiverFile, 'a') as file:
        for station, coords in stationDict.items():
            file.write(f'{coords[0]} {coords[1]}\n')
        file.close()

    open(sourcesFile,'w').close()
    with open(sourcesFile, 'a') as file:
        for station, coords in stationDict.items():
            file.write(f'{coords[0]} {coords[1]}\n')
        file.close()

    stationList = list(stationDict.keys())

    issue_dict = {'autocorrelation' : 0,
                  'filepath not exist' : 0,
                  'fparam_out_none' : 0,
                  'interpOut none' : 0,
                  'bad phvel' : 0,
                  'low snr' : 0,
                  'too close' : 0,
                  'good' : 0}

    open(timesFile,'w').close()

    fpDict = {'infDepen' : 0,
              'omdom' : 0,
              'zero' : 0,
              'slice' : 0,
              'fparam returning none' : 0}

    interpErrorDict = {'lowest observed period greater than desired' : 0,
                       'highest observed period lower than desired' : 0,
                       'could not find left or right' : 0,
                       'unreasonable phvel from interp' : 0,
                       'unreasonable snr' : 0,
                       'bad left phvel' : 0,
                       'bad right phvel' : 0}
    good_phvel_list = []
    good_snr_list = []

    with open(timesFile, 'a') as outfile:
        for i, stat1 in tqdm(enumerate(stationList),total=len(stationList)):
            for j, stat2 in enumerate(stationList):
                if stat1 == stat2:
                    issue_dict['autocorrelation'] += 1

                    continue

                filepath = checkIfCorrelationExists(stat1,stat2,componentDirectory)
                if filepath is None:
                    issue_dict['filepath not exist'] += 1
                    outfile.write('0 0.0000 1.0\n')
                    continue

                fparam_out = pt.FTAN(filepath,tmin=1,tmax=40,nfin=41)
                fparam_out, fpDict = _fparamErrorCodeHandler(fparam_out,fpDict,outfile)
                if fparam_out is None:
                    issue_dict['fparam_out_none'] += 1
                    continue

                fparam = fparam_out[0]
                dist = fparam_out[1]
                obper1, gvel1, phvel1, snr = pt.pullRelevantFTANInfo(fparam)

                interpOut = interpPeriod(period,obper1,phvel1,snr)
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

                good_phvel_list.append(phvel)
                good_snr_list.append(snr)

                travelTime = round(float(getTravelTime(phvel,dist)), 4)
                issue_dict['good'] += 1
                outfile.write(f'1.0 {travelTime} 1.0\n')

    return issue_dict, fpDict, interpErrorDict, good_phvel_list, good_snr_list

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

def interpPeriod(period,obper1,phvel1,snr):
    """Interpolates the phase velocity and snr to the period of interest"""
    if obper1[0] > period:
        return 1

    if obper1[-1] < period:
        return 2

    left = None
    right = None
    for i, obper in enumerate(obper1):
        if obper < period:
            left = i
        if obper > period:
            right = i
            break

    if left is None or right is None:
        return 3

    if phvel1[left] < 1.5:
        left -= 1

    if phvel1[left] < 1.5:
        return 6

    if phvel1[right] < 1.5 or phvel1[right] > 5:
        return 7

    if snr[left] < 2 or snr[right] < 2:
        return 5

    interpPhvel = np.interp(period,[obper1[left],obper1[right]],[phvel1[left],phvel1[right]])
    interpSNR = np.interp(period,[obper1[left],obper1[right]],[snr[left],snr[right]])

    if interpPhvel < 1.5 or interpPhvel > 5:
        return 4

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
    if interpPeriodOut == 4:
        interpPeriodDict['unreasonable phvel from interp'] += 1
    if interpPeriodOut == 5:
        interpPeriodDict['unreasonable snr'] += 1
    if interpPeriodOut == 6:
        interpPeriodDict['bad left phvel'] += 1
    if interpPeriodOut == 7:
        interpPeriodDict['bad right phvel'] += 1

    return None, interpPeriodDict

def checkIfCorrelationExists(stat1,stat2,componentDirectory):
    if os.path.exists(componentDirectory + f'/{stat1}_{stat2}.sac'):
        return componentDirectory + f'/{stat1}_{stat2}.sac'

    if os.path.exists(componentDirectory + f'/{stat2}_{stat1}.sac'):
        return componentDirectory + f'/{stat2}_{stat1}.sac'

    return None

def getTravelTime(vel, dist):
    """Returns travel time given velocity and distance"""
    return dist / vel

def getWavelength(per,vel):
    """Returns the wavelength given period and velocity"""
    return per * vel

def plotIssueDict(issue_dict, label):
    fig, ax = plt.subplots()

    issues = list(issue_dict.keys())
    count = list(issue_dict.values())

    ax.bar(issues, count)

    plt.xticks(rotation=65,ha='right')
    plt.title(label)
    plt.show()


network = 'UW,CC,XU,XD,TA,YH'
stations = 'ALL'
channel='BH*,HH*,EH*'
bound_box = [46.1,47.3,-122.5,-120.9]
dataDirectory = '/Volumes/NewHDant/RainierAmbient'

stationList = getLocalStations(dataDirectory,'ZZ')
stationDict = getValidStations(network,bound_box,channel,stationList)

shortenedDict = {k: v for k, v in list(stationDict.items())[:20]}

periods = [5]
makeTomoDirectory(dataDirectory,periods=periods,component='ZZ')

avg_list = []
for period in periods:
    issue_dict, fpDict, interpErrorDict, good_phvel_list, good_snr_list = good_count = makeFMSTInputs(stationDict,dataDirectory,
                                             period=period,
                                             component='ZZ',
                                             minSNR=5,
                                             minWavelengths=2)
    avg_list.append(np.mean(good_phvel_list))

plt.plot(periods,avg_list)
plt.title('Avg Phase Velocity Over Study Area')
plt.xlabel('Period (s)')
plt.ylabel('Phase Velocity (km/s)')

"""
plotIssueDict(issue_dict,'Overall Issues')
plotIssueDict(fpDict,'FTAN Issues')
plotIssueDict(interpErrorDict,'Interpolation Issues')

plt.hist(good_phvel_list)
plt.title('Measured phase velocities at 5s')
plt.show()
plt.hist(good_snr_list)
plt.show()

print(f'mean of phase vels: {np.mean(good_phvel_list)}')
print(f'std of phase vels: {np.std(good_phvel_list)}')

#46.875 -121.625    0.250 -0.322976E+01 2.96 = 2.009 km/s
"""
