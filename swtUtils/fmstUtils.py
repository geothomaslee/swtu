#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 07:55:48 2024

@author: Thomas Lee
University of New Mexico
Department of Earth and Planetary Science

A module for creating properly formatted inputs for Fast-Marching Surface
Tomography (FMST) by Nick Rawlinson using outputs from AFTAN by Bensen et. al.
Assumes the stacked cross-correlation functions are in the format and file
structure outputted by ambient2.

Citations:
    [1] Bensen, G. D., et al. (2007). Processing seismic ambient noise data to
    obtain reliable broad-band surface wave dispersion measurements.
        Geophysical Journal International, 169(3), 1239-1260.
    [2] Rawlinson, N. and Sambridge M., 2005. "The fast marching method: An
        effective tool for tomographic imaging and tracking multiple phases in
        complex layered media", Explor. Geophys., 36, 341-350.
"""

import os
import pickle
import shutil
from glob import glob
from typing import List,Tuple,Union
from math import cos, atan2, radians, sin, sqrt

from tqdm import tqdm
import numpy as np
import obspy
from obspy.core.inventory import Inventory
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
from deprecated import deprecated
import pandas as pd
from pandas import DataFrame

import swtUtils.ftan as ft #pylint: disable=import-error
#pylint: disable=invalid-name

def saveObj(obj, filename: str):
    """Quick function for pickling a file"""
    with open(filename, 'wb') as f:
        pickle.dump(obj, f)

def loadObj(filename: str):
    """Quick function for pickling a file"""
    with open(filename, 'rb') as f:
        return pickle.load(f)

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

def getTomoDirectory(dataDirectory: str,component: str) -> str:
    """Returns the tomography directory for the component of interest"""
    if not isinstance(component, str):
        raise TypeError('Component must be 2 character string. Ex. "ZZ"')
    if len(component) != 2:
        raise ValueError('Component must be 2 character string. Ex. "ZZ"')

    tomoDirectory = f'{dataDirectory}/Tomography/{component}'
    return tomoDirectory

def getStationNames(filepath: str) -> Tuple[str, str]:
    """Pulls the name of the two stations from the name of the file"""
    basename = os.path.basename(filepath)
    name_no_sac = basename[:-4]
    stat1, stat2 = name_no_sac.split('_')

    return stat1, stat2

def getLocalStations(
        dataDirectory : str,
        component : str,
        forceOverwrite: bool = True) -> List[str]:
    """Returns a list of stations that we have cross-correlations for"""
    componentDirectory = getComponentDirectory(dataDirectory,component)

    if forceOverwrite is True:
        os.remove(componentDirectory + '/UniqueStations.pkl')

    try:
        if not os.path.exists(componentDirectory +'/UniqueStations.pkl'):
            crossCorrelations = glob(componentDirectory + '/*.sac')
            station_list = []
            print('Searching for local stations...')
            for file in tqdm(crossCorrelations):
                stat1, stat2 = getStationNames(file)

                station_list = [s for s in set(station_list + [stat1, stat2]) if s]
                saveObj(station_list,componentDirectory +'/UniqueStations.pkl')
        else:
            print('Searching for local stations...')
            print('Using existing version of UniqueStations.pkl - Is this intentional?')
            station_list = loadObj(componentDirectory +'/UniqueStations.pkl')
    except EOFError:
        print('EOFError often means that UniqueStations.pkl exists but was written improperly, delete it and try again')
    except Exception as e:
        raise e

    return station_list

def getValidStations(
        network : str,
        bounds : Tuple[float,float,float,float],
        channel : str,
        stationList : List[str],
        client: str='IRIS') -> dict:
    """
    Cross-references our list of stations that we have data for in the
    cross-correlation folder with a search of stations in the IRIS database,
    and builds a nice dictionary out of it that is used by makeFMSTInputs.

    Parameters
    ----------
    network : str
        Network(s).
    bounds : list of ints
        [minlat,maxlat,minlon,maxlon]. Study area
    channel : str
        Channels of interest.
    stationList : list of str
        List of all station names that have been found within the stack directory.

    Returns
    -------
    stationDict : dict
        {network.station : (lat,lon).

    """
    _client = Client(client)
    inventory = _client.get_stations(network = network,
                                    station = '*',
                                    channel = channel,
                                    minlatitude = bounds[0],
                                    maxlatitude = bounds[1],
                                    minlongitude = bounds[2],
                                    maxlongitude = bounds[3])

    inventoryCopy = inventory.copy()

    for net in inventory:
        for station in net:
            formattedName = f'{net.code}.{station.code}'
            if formattedName not in stationList:
                inventoryCopy = inventoryCopy.remove(network=net.code,station=station.code)

    stationDict = {}
    for net in inventoryCopy:
        for station in net:
            stationDict[f'{net.code}.{station.code}'] = (station.latitude,station.longitude)

    return stationDict

def getInventoryLength(inventory : Inventory) -> int:
    """Returns the number of stations in your inventory"""
    length = 0
    for network in inventory:
        length += len(network)

    return length

def makeTomoDirectory(
        dataDirectory : str,
        periods : str,
        component : str) -> str:
    """Periods must be a list of all periods, or if only 2 elements, is assumed
    to want every integer period between element 1 and element 2.

    Ex. [1, 5, 7, 8] would use those periods, [1, 8] would use every integer from
    1 to 8 as the periods.
    """
    if not isinstance(periods, list):
        raise TypeError('Periods must be a list, see doc for details')
    if len(periods) < 2:
        raise ValueError('Periods must be a list of length >=2, see docs')
    if len(periods) == 1:
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

def get_all_ftan_outputs(FTANDirectory: str) -> list:
    all_ftan_outputs = glob(f'{FTANDirectory}/*_1_DISP.1')
    return set(all_ftan_outputs)

def getReferenceVelocity(
        stationDict : dict,
        dataDirectory : str,
        FTANDirectory : str,
        period : Union[int,float],
        component : str,
        minSNR : Union[int,float],
        minWavelengths : Union[int,float]) -> Tuple[float,list]:
    """
    For every available FTAN file, finds the velocity for a given period. Will
    return the average of all these measurements to be used as the starting
    velocity for phase velocity inversions, as well as a list of every phase
    velocity measurement for plotting histograms.

    Parameters
    ----------
    stationDict : dict
        Dictionary containing all station information in the format
        {'network.station' : (lat, lon)}
    dataDirectory : str
        Path to data directory.
    FTANDirectory : str
        Path to FTAN directory.
    period : Union[int,float]
        Period of interest.
    component : str
        Cross-component. Ex: 'ZZ', 'NE'
    minSNR : Union[int,float]
        Minimum signal-to-noise ratio.
    minWavelengths : Union[int,float]
        Minimum number of wavelengths between stations to consider. If too close
        phase velocity measurements are inaccurate.

    Returns
    -------
    avgPhvel : float
        Average phase velocity across the entire study area for a given period.
    phvels : list[float]
        List of every measured phase velocity for that period

    """
    componentDirectory = getComponentDirectory(dataDirectory,component)
    stationList = list(stationDict.keys())
    interpErrorDict = makeInterpErrorDict()

    ftan_output_list = get_all_ftan_outputs(FTANDirectory)

    phvels = []
    for stat1 in tqdm(stationList):
        for stat2 in stationList:
            if stat1 == stat2:
                continue

            filepath = checkIfFTANExists(stat1,stat2,FTANDirectory,ftan_output_list)
            if filepath is None:
                continue

            dist = getDist(stat1,stat2,componentDirectory)
            if dist is None:
                continue

            df = ft.dispOutputToDF(filepath)
            obper,phvel,snr = ft.getRelevantInfo(df)

            interpOut = interpPeriod(period,obper,phvel,snr)
            interpOut, interpErrorDict = _interpPeriodErrorHandler(interpOut,interpErrorDict)
            if interpOut is None:
                continue

            phvel = interpOut[0]
            snr = interpOut[1]

            if phvel < 1.5 or phvel > 5:
                continue

            if snr < minSNR:
                continue

            minDist = minWavelengths * getWavelength(period,phvel)
            if dist < minDist:
                continue

            phvels.append(phvel)

    avgPhvel = round(float(np.mean(phvels)),4)
    return avgPhvel, phvel

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

def makeFMSTInputs(stationDict,
                   dataDirectory,
                   FTANDirectory,
                   period,
                   component,
                   minSNR,
                   minWavelengths,
                   writeDists=True,
                   detailedError=True,
                   included_networks=None):
    """
    =====OUTPUT FILE INFO======
    Will makes the 3 input files needed by FMST - sources.dat,receivers.dat,
    and otimes.dat.

    sources.dat - list of all sources
    recivers.dat - list of all receivers
    otimes.dat - list of travel times between combinations

    otimes.dat is expected to be the travel time between every source-receiver
    combination, even if no travel time exists. A 0 is written into the first column
    in that case because FMST still expects a line for every single S-R pair.

    receivers.dat and sources.dat are identical for ambient noise tomography.

    This script will iterate through every possible station-station pair, and will
    give a 0 to otimes in the following cases:

    1. If the station-station pair is an autocorrelation
    2. If no corresponding output from FTAN can be found
    3. If the original SAC file cannot be found to pull a distance from.
    4. If the interpolator returns None.
    5. If the measured phase velocity is below 1.5 km/s or greater than 5 km/s.
    6. If the signal-to-noise ratio is less than minSNR, a function parameter.
    7. If the interstation distance is less than minWavelengths
        - The wavelength is calculated using the phase velocity and period

    ======CALCULATING PHASE VELOCITIES=======
    This is assumed to have already been done by the original AFTAN fortran
    script from Bensen et. al,. I have written this script to assume the file
    structure that is inherited from a script called runFTAN.csh, originally by
    Justin Wilgus as a PhD student in Brandon Schmandt's lab. Chances are if
    you're reading this, you're inheriting this code and therefore will also
    have that script.

    AFTAN automatically adjust the periods it tests based on the SNR, so the outputs
    must be interpolated if you want an exact integer period. This also has to
    interpolate the SNR. See the interpolator itself to see how it works, but it
    does just a linear interpolation based on the nearest points that bound
    the integer period of interest.

    You may encounter a situation where AFTAN actually backtracks at short periods especially,
    and therefore a linear interpolation will break. I haven't had this issue but
    this has been reported to me by others.

    ========ERROR DICTIONARIES============
    With cross-correlation functions its expected that most of your data won't
    pass quality control, but within reason. It's good to know why you're getting
    a 0 value in otimes instead of just moving on. This script creates 3 different
    error dictionaries based on the returned error codes of certain functions. The error
    codes themselves are just integers from the script that spits them out, but
    these functions are themselves wrapped in an error handler that can properly add
    to the dictionary. These are saved inside of the tomography directory as .pkl
    files which need to be read by the pickle module. I have this done already
    inside of swtu.main and it prints them when running main.

    Parameters
    ----------
    stationDict : dict
        Dict in the format from getValidStations.
    dataDirectory : str
        Full filepath to the directory where your stacked cross-correlations
        are kept, with a file structure inherited from
        github/thomaslee/ambient2-tlee-fork.
    FTANDirectory : str
        Full filepath to the directory where FTAN is located. This should be
        the directory that CONTAINS runFTAN.csh
    period : int
        Period of interest.
    component : str
        Component.
    minSNR : int or float
        Minimum signal to noise ratio for phase velocity measurements.
    minWavelengths : int or float
        Minimum number of wavelengths, below which the interstation distance
        is considered too short.
    detailedError : bool, optional
        If true, will pickle the error dictionaries. The default is True.

    Returns
    -------
    phvels : list
        List of every measured phase velocity, generally used for plotting
        a histogram to see if your measurements are generally correct.
    """

    if included_networks is not None:
        print('WARNING: included networks is defined, meaning certain networks may be excluded. Set to None to use all')
        if isinstance(included_networks,list) != True:
            raise TypeError('included_networks must be a list of strings, where each string is a station')
        else:
            for _net_ in included_networks:
                if isinstance(_net_, str) != True:
                    raise TypeError('included_networks must be a list of strings, where each string is a station')

    # Get directories
    tomoDirectory = getTomoDirectory(dataDirectory,component)
    periodDirectory = tomoDirectory + f'/{period}s'
    componentDirectory = getComponentDirectory(dataDirectory,component)

    # Get filepaths for output files
    receiverFile = periodDirectory +'/receivers.dat'
    sourcesFile = periodDirectory +'/sources.dat'
    timesFile = periodDirectory + '/otimes.dat'

    ftanMainDirectory = os.path.dirname(FTANDirectory)
    distFile = ftanMainDirectory + '/distances.dat'

    # Creating issue dicts
    issue_dict = makeIssueDict()
    fpDict = makefpDict()
    interpErrorDict = makeInterpErrorDict()

    # Create receiver file
    open(receiverFile,'w',encoding='utf-8').close()
    with open(receiverFile, 'a',encoding='utf-8') as file:
        file.write(f'{int(len(list(stationDict.keys())))}\n')
        for coords in stationDict.values():
            file.write(f'{coords[0]} {coords[1]}\n')
        file.close()

    # Create source file
    open(sourcesFile,'w',encoding='utf-8').close()
    with open(sourcesFile, 'a',encoding='utf-8') as file:
        file.write(f'{int(len(list(stationDict.keys())))}\n')
        for coords in stationDict.values():
            file.write(f'{coords[0]} {coords[1]}\n')
        file.close()

    stationList = list(stationDict.keys())

    open(timesFile,'w',encoding='utf-8').close()

    ftan_output_list = get_all_ftan_outputs(FTANDirectory)

    phvels = []
    dists = []
    stat1_list = []
    stat2_list = []

    with open(timesFile, 'a',encoding='utf-8') as outfile:
        #outfile.write(f'{len(stationList)**2}\n')
        for stat1 in tqdm(stationList):
            for stat2 in stationList:
                if stat1 == stat2:
                    issue_dict['autocorrelation'] += 1
                    outfile.write('0 0.0000 1.0\n')
                    continue

                filepath = checkIfFTANExists(stat1,stat2,FTANDirectory,ftan_output_list)

                if included_networks is not None:
                    if stat1.split('.')[0] not in included_networks:
                        outfile.write('0 0.0000 1.0\n')
                        continue

                    if stat2.split('.')[0] not in included_networks:
                        outfile.write('0 0.0000 1.0\n')
                        continue

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
                dists.append(dist)
                stat1_list.append(stat1)
                stat2_list.append(stat2)

                travelTime = round(float(getTravelTime(phvel,dist)), 4)
                issue_dict['good'] += 1
                outfile.write(f'1 {travelTime} 1.0\n')


    if writeDists:
        with open(distFile, "w") as distances_file:
            for stat1 in stationList:
                for stat2 in stationList:
                    lat1 = stationDict[stat1][0]
                    lon1 = stationDict[stat1][1]
                    lat2 = stationDict[stat2][0]
                    lon2 = stationDict[stat2][1]
                    dist = haversine_distance(lat1,lon1,lat2,lon2)
                    dist_str = f'{dist:.3f} {lat1} {lon1} {lat2} {lon2}\n'
                    distances_file.write(dist_str)


    if detailedError is True:
        saveObj(issue_dict, f'{periodDirectory}/issueDict.pkl')
        saveObj(fpDict, f'{periodDirectory}/fpDict.pkl')
        saveObj(interpErrorDict,f'{periodDirectory}/interpErrorDict.pkl')

    return phvels, dists, stat1_list, stat2_list

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

def checkIfFTANExists(stat1,stat2,FTANDirectory,ftan_output_list):
    """Checks if an FTAN output from AFTAN (Bensen) exists"""
    if f'{FTANDirectory}/{stat1}_{stat2}_Folded.sac_1_DISP.1' in ftan_output_list:
        return FTANDirectory + f'/{stat1}_{stat2}_Folded.sac_1_DISP.1'

    if f'{FTANDirectory}/{stat2}_{stat1}_Folded.sac_1_DISP.1' in ftan_output_list:
        return FTANDirectory + f'/{stat2}_{stat1}_Folded.sac_1_DISP.1'

    return None

def getDist(stat1,stat2,componentDirectory):
    """Given two stats, finds dist from original SAC file"""
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
    fig, ax = plt.subplots() #pylint: disable=unused-variable

    issues = list(issue_dict.keys())
    count = list(issue_dict.values())

    ax.bar(issues, count)

    plt.xticks(rotation=65,ha='right')
    plt.title(label)
    plt.show()

def setupFTANDirectory(FMSTDirectory,period,projectCode,component,_overwrite=False):
    """
    Creates a new FMST run directory inside the base FMST directory.

    FMST requires a list of inputs that are called by several different scripts,
    but all of them can be called sequentially using ttomoss. The travel time
    file will be unique to that period, so every period we want to perform an
    inversion for needs its own directory.

    This uses a master template file called {projectCode}_Master inside the FMST
    directory. It contains all the input files that should be standard across
    all runs. Those files are copied into the new directory, which will be called
    {project_code}_{period}s_{component}.

    Example, project Rainier will look for its template in Rainier_Master, then
    copy its contents into Rainier_5s_ZZ if we give this function a period of 5
    and a component of ZZ.

    Parameters
    ----------
    FMSTDirectory : str
        Full path to the FMST directory..
    period : int
        Period of interest.
    projectCode : str
        Name of project. See full docstring for details.
    component : str
        Component. Ex. "ZZ","NE".
    _overwrite : bool, optional
        If set True, will overwrite the existing directory if it already exists
        This is useful if you're playing around with the inversion parameters and
        want a clean slate. The default is False.

    Returns
    -------
    fmstPath : str
        Path to the newly created FMST run directory.
    """
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

def getAvgVelocity(period,df,component):
    """Uses the phase velocity DataFrame to find the average velocity for that period"""
    if df.empty is True:
        raise ValueError('Empty DataFrame given to getAvgVelocity')
    try:
        vels = df[f'{float(period)}s'].to_list()
    except KeyError:
        print('Key error. Almost certainly moves this period is not in the CSV')
        return 'no_key_in_csv'

    vels = [x for x in vels if str(x) != 'nan']

    avgVel = np.mean(vels)

    return avgVel

def editBackgroundVel(fmstPath,avgPhvel,colorbar_margin: float=0.2):
    """
    Edits the background velocity in grid2dss.in, the initial grid creation
    input file for FMST, to match the average measured phase velocity from FTAN.

    Parameters
    ----------
    fmstPath : str
        Full path to the FMST base directory.
    avgPhvel : float or int
        Average phase velocity.

    Returns
    -------
    None.

    """
    filepath = f'{fmstPath}/mkmodel/grid2dss.in'
    if not os.path.isfile(filepath):
        raise ValueError('Could not find grid2dss.in')

    with open(filepath,"r",encoding='utf-8') as infile:
        lines = infile.readlines()
        newline =  str(avgPhvel) + lines[14][len(str(avgPhvel)):]
        lines[14] = newline

    with open(filepath,"w",encoding='utf-8') as outfile:
        for line in lines:
            outfile.write(line)

    filepath = f'{fmstPath}/gmtplot/mapping.py'
    print(filepath)

    minvel = avgPhvel - (avgPhvel*colorbar_margin)
    maxvel = avgPhvel + (avgPhvel*colorbar_margin)
    with open(filepath, 'r') as infile:
        lines = infile.readlines()
        print(len(lines))
        for i, line in enumerate(lines):
            if line.strip() == 'plot_phase_vel_modern()':
                newline = f'    plot_phase_vel_modern(minvel={minvel},maxvel={maxvel})'
                lines[i] = newline

    with open(filepath,"w",encoding='utf-8') as outfile:
        for line in lines:
            outfile.write(line)

def moveFMSTInputs(fmstPath,tomoDirectory,_overwrite=False):
    """Moves FMST inputs to their corresponding FMST directory"""
    files = ['sources.dat','receivers.dat','otimes.dat']
    for file in files:
        filepath = fmstPath + f'/{file}'
        if os.path.isfile(filepath) is True:
            if _overwrite is False:
                continue

            os.remove(filepath)

        shutil.copy(f'{tomoDirectory}/{file}',fmstPath)

def findAllFinalTomoImages(fmstPath,projectCode,component,periods):
    """Finds all final tomo images and puts them inside the primary FMST dir"""
    imgdir = fmstPath + '/phvelMaps'
    if not os.path.isdir(imgdir):
        os.mkdir(imgdir)

    for period in periods:
        filename = f'{projectCode}_{period}s_{component}/gmtplot/tmp.png'
        fullpath = f'{fmstPath}/{filename}'
        shutil.copyfile(src=fullpath,
                        dst=f'{imgdir}/phvel_{period}s.png')

"""
<<<<<<< Updated upstream
=======
"""

def create_output_files(fmstDirectory: str,
                        periods: List[float]):
    """Creates the output files for running many combos of inversion parameters"""
    for period in periods:
        file_name = f'{fmstDirectory}/{period}s_Outputs'
        if os.path.exists(file_name):
            os.remove(file_name)

        with open(file_name, 'w',encoding='utf-8') as file:
            file.write('')

def get_output_info(fmstDir: str):
    """
    Parameters
    ----------
    fmstDir : str
        Path to the FMST directory for this exact run (not the overall FMST
        directory, but the local directory for this exact run.

    Returns
    -------
    rms : float
        Root mean square misfit in milliseconds of the last iteration
    variamce : float
        Variance of the last iteration
    """

    with open (f'{fmstDir}/residuals.dat',encoding='utf-8') as infile:
        lines = infile.readlines()
        rms, variance = lines[-1].split('    ')
        variance = variance.rstrip()

    return rms, variance

def edit_inversion_params(fmstPath: str,
                          smoothing: float,
                          damping: float,
                          lon_grids: int,
                          lat_grids: int):
    """
    Parameters
    ----------
    fmstDir : str
        Path to the FMST directory for this exact run (not the overall FMST
        directory, but the local directory for this exact run.
    smoothing : float
        Smoothing parameter.
    damping : float
        Damping parameter.
    lon_grids : int
        Number of grid cells in the longitudinal axis.
    lat_grids : int
        Number of grid cells in the latitudinal axis.

    Returns
    -------
    None.

    """
    # Editing grid size
    filepath = f'{fmstPath}/mkmodel/grid2dss.in'
    if not os.path.isfile(filepath):
        raise ValueError('Could not find grid2dss.in')

    with open(filepath,"r",encoding='utf-8') as infile:
        lines = infile.readlines()
        lines[7] = f'{lat_grids}                    c: Number of grid points in theta (N-S)\n'
        lines[8] = f'{lon_grids}                    c: Number of grid points in theta (N-S)\n'

    with open(filepath,"w",encoding='utf-8') as outfile:
        for line in lines:
            outfile.write(line)

    # Editing smoothing and damping
    filepath = f'{fmstPath}/subinvss.in'
    if not os.path.isfile(filepath):
        raise ValueError('Could not find grid2dss.in')

    with open(filepath,"r",encoding='utf-8') as infile:
        lines = infile.readlines()
        lines[11] = f'{damping}                       c: Damping factor (epsilon)\n'
        lines[14] = f'{smoothing}                      c: Smoothing factor (eta)\n'

    with open(filepath,"w",encoding='utf-8') as outfile:
        for line in lines:
            outfile.write(line)

def make_residuals_dfs(ftanDir: str,
                       fmstDir: str,
                       normalize: bool=True,
                       debug: bool=False):
    """
    Makes a DataFrame

    Parameters
    ----------
    ftanDir : str
        FTAN Directory
    fmstDir : str
        DESCRIPTION.
    debug : bool, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    rtravel_df : TYPE
        DESCRIPTION.
    otimes_df : TYPE
        DESCRIPTION.
    residuals_df : TYPE
        DESCRIPTION.

    """
    rtravel = f'{fmstDir}/rtravel.out' # Source-receiver travel times through solution model
    otimes = f'{fmstDir}/otimes.dat' # Source-receiver times from real data
    distances = f'{ftanDir}/distances.dat'

    with open(rtravel,"r",encoding='utf-8') as rtravelfile:
        travel_times = []
        indices = []
        use_list = []
        for i, line in enumerate(rtravelfile.readlines()):
            travel_times.append(float(line[:-4].split('   ')[-1]))
            use_list.append(int(line[:-4].split('   ')[-2][-1:]))
            indices.append(i)

        df_dict = {'Index': indices, 'Travel_Times': travel_times, 'Use': use_list}
        rtravel_df = pd.DataFrame(df_dict)

    with open(otimes,"r",encoding='utf-8') as otimes_file:
        travel_times = []
        use = []
        indices = []
        for i, line in enumerate(otimes_file.readlines()):
            travel_times.append(float(line.split(' ')[1]))
            use.append(float(line.split(' ')[0]))
            indices.append(i)

        df_dict = {'Index': indices, 'Use': use, 'Travel_Times': travel_times}
        otimes_df = pd.DataFrame(df_dict)

    dists = []
    with open(distances,'r',encoding='utf-8') as distfile:
        lines = distfile.readlines()

        for line in lines:
            entries = line.split(' ')
            dist = entries[0]
            dists.append(float(dist))

    residuals = []
    rtravel_times = rtravel_df['Travel_Times'].to_list()
    use_travel = rtravel_df['Use'].to_list()
    otravel_times = otimes_df['Travel_Times'].to_list()

    """
    # PLACEHOLDER PLACEHOLDER
    dists = [1 for x in rtravel_times]
    # PLACEHOLDER PLACEHOLDER
    """

    if debug:
        print('===============PRINTING LARGE RESIDUAL CALCULATIONS===============')
    for i, otravel_time in enumerate(otravel_times):
        rtravel_time = rtravel_times[i]
        use_here = use_travel[i]
        dist = dists[i]

        # If the predicted travel time isn't actually used, set the residual to 0
        if use_here == 0:
            residual = 0
        else:
            if normalize:
                residual = (rtravel_time - otravel_time) / dist
            else:
                residual = rtravel_time - otravel_time
        if debug:
            if np.abs(residual) > 20:
                print(f'Pre {rtravel_time} - Obs {otravel_time} = {residual} in line {i+1}')

        residuals.append(residual)

    if debug:
        print('==================================================================')

    df_dict = {'Index': indices, 'Residuals': residuals}
    residuals_df = pd.DataFrame(df_dict)

    return rtravel_df, otimes_df, residuals_df

def get_residual_statistics(residuals_df: DataFrame):
    filtered_df = residuals_df[(residuals_df['Residuals'] != 0)]
    max_residuals = max([np.abs(x) for x in filtered_df['Residuals'].to_list()])
    avg_residuals = np.mean(filtered_df['Residuals'].to_list())
    std_residuals = np.std(filtered_df['Residuals'].to_list())
    avg_abs_residuals = np.mean([np.abs(x) for x in filtered_df['Residuals'].to_list()])
    print('==================================================================')
    print(f'Maximum Residual: {max_residuals:.5f}')
    print(f'Mean: {avg_residuals:.5f}, Standard Deviation: {std_residuals:.5f}')
    print(f'Mean Residual Magnitude {avg_abs_residuals:.5f}')

    return std_residuals

def remove_worst_fits(otimes_df: DataFrame,
                    residuals_df: DataFrame,
                    std_residuals: float,
                    stds: float,
                    debug: bool=False):

    # Set travel times to not be used if they come from
    cleaned_otimes_df = otimes_df.copy()
    min_residual = stds*std_residuals * -1
    max_residual = stds*std_residuals

    if debug:
        print(f"Total Residuals: {len(residuals_df[(residuals_df['Residuals'] != 0)])}")
        print(f"Total Otimes: {len(otimes_df[(otimes_df['Use'] != 0)])}")

    outside_range_mask = (residuals_df['Residuals'] < min_residual) | (residuals_df['Residuals'] > max_residual)

    if debug:
        print('BAD RESIDUALS THAT ARE BEING REMOVED')
        print(residuals_df[outside_range_mask]['Residuals'].to_list())

    cleaned_otimes_df.loc[outside_range_mask, 'Use'] = 0.0
    indices = np.where(outside_range_mask)[0].tolist()

    print(f"Number to Remove: {len(residuals_df[outside_range_mask]['Index'].to_list())}")
    print(f"Otimes After Cleaning: {len(cleaned_otimes_df[(cleaned_otimes_df['Use'] != 0)])}")
    print('==================================================================')

    return cleaned_otimes_df, indices

def rewrite_otimes(fmstDir: str,
                   otimes_df: DataFrame):

    use_list = [int(x) for x in otimes_df['Use'].to_list()]
    times_list = otimes_df['Travel_Times'].to_list()
    otimes_file = f'{fmstDir}/otimes.dat' # Source-receiver times from real data

    with open(otimes_file,'w',encoding='utf-8') as otimes:
        for i, use in enumerate(use_list):
            travel_time = times_list[i]
            if float(travel_time) == 0.0:
                travel_time = '0.0000'
            line = f'{use} {travel_time} 1.0\n'
            otimes.write(line)

def plot_residuals(residuals_df: DataFrame,
                   title: str,
                   normalized: bool=True,
                   period: float=None,
                   percent: float=None,
                   binsize: float=5):
    filtered_df = residuals_df[(residuals_df['Residuals'] != 0)]
    residuals = filtered_df['Residuals'].to_list()
    fig, ax = plt.subplots(figsize=(4,8))
    if normalized is True:
        ax.hist(residuals)
    else:
        ax.hist(residuals,bins=range(int(min(residuals)),int(max(residuals))+binsize,binsize))

    if normalized is True:
        ax.set_xlabel('Normalized Travel Time\nResiduals (s/km)')
    else:
        ax.set_xlabel('Travel Time Residuals (s)')
    ax.set_ylabel('Frequency')

    title_str = f'Residuals {title} Removing\nWorst {percent:.4f} Percent'

    if normalized is True:
        title_str = 'Normalized ' + title_str

    if period is not None:
        title_str = title_str + f'\n{int(period)}s Period'

    ax.set_title(title_str)

    return fig, ax

def plot_smoothing_damping(outputs_file: str,
                           smoothing_list: list,
                           damping_list: list,
                           label_cells: bool=False):
    """
    Visually plots a matrix that shows the model RMS for different combinations
    of smoothing and damping

    Parameters
    ----------
    outputs_file : str
        File containing the information outputted by runInversions in main.
    smoothing_list : list
        List of unique smoothing parameters that were used. NOT a list of every
        smoothing value used (5 smoothing and 5 damping results in 25 combos,
        give the list of 5, not 25)
    damping_list : list
        List of unique damping parameters that were used. NOT a list of every
        damping value used (5 smoothing and 5 damping results in 25 combos,
        give the list of 5, not 25)
    label_cells : bool, optional
        If True, will label the combination inside every cell. Illegible and bad
        but I figured I would leave it just in case? The default is False.

    Returns
    -------
    None.

    """
    cols = ['Smoothing','Damping','Lon_Grids','Lat_Grids','RMS','Variance','Variability', 'Variance_KMS','Roughness']
    df = pd.read_csv(outputs_file,names=cols,sep='\s+')

    rms = np.array(df['RMS'].to_list())
    print(rms)
    rms_matrix = rms.reshape(len(smoothing_list),len(damping_list))

    smoothing_labels = [str(x) for x in smoothing_list]
    damping_labels = [str(x) for x in damping_list]

    fig, ax = plt.subplots()
    im = ax.imshow(rms_matrix,interpolation='none',cmap='viridis')

    if label_cells:
        for i in range(len(damping_labels)):
            for j in range(len(smoothing_labels)):
                ax.text(j,i,f'{smoothing_labels[j]},{damping_labels[i]}',ha='center',va='center',color='w')

    ax.set_xticks(np.arange(len(smoothing_labels)))
    ax.set_xticklabels(smoothing_labels)

    ax.set_yticks(np.arange(len(damping_labels)))
    ax.set_yticklabels(damping_labels)

    ax.set_xlabel('Smoothing')
    ax.set_ylabel('Damping')

    cbar = ax.figure.colorbar(im,ax=ax)
    cbar.ax.set_ylabel('Total Misfit (ms)',rotation=-90,va='bottom')

    #plt.colorbar()
    plt.show()

    return fig

def reset_FMST_directory(dataDirectory: str,
                         fmstDirectory: str,
                         ftanDirectory: str,
                         component: str,
                         period: int,
                         projectCode: str,
                         snr: float,
                         min_wavelengths: float,
                         damping: float,
                         smoothing: float,
                         lon_grids: int,
                         lat_grids: int,
                         colorbar_margin: float=0.2,
                         included_networks=None):
    """Sets up a properly formatted FMST directory"""

    if included_networks is not None:
        print('WARNING: included networks is defined, meaning certain networks may be excluded. Set to None to use all')
        if isinstance(included_networks,list) != True:
            raise TypeError('included_networks must be a list of strings, where each string is a station')
        else:
            for _net_ in included_networks:
                if isinstance(_net_, str) != True:
                    raise TypeError('included_networks must be a list of strings, where each string is a station')


    tomoDirectory = getTomoDirectory(dataDirectory,component) + f'/{period}s'
    fmstPath = setupFTANDirectory(FMSTDirectory=fmstDirectory,
                                                period=period,
                                                projectCode=projectCode,
                                                component=component,
                                                _overwrite=True)

    moveFMSTInputs(fmstPath=fmstPath,
                            tomoDirectory=tomoDirectory,
                            _overwrite=False)

    df = pd.read_csv(f'{ftanDirectory}/PhaseVelocities_SNR{snr}_WL{min_wavelengths}.csv')

    # We want to use the new average velocity if we only use certain networks
    if included_networks is not None:
        for _network in included_networks:
            df = df[df.Station1.str.contains('|'.join(_network))]
            df = df[df.Station2.str.contains('|'.join(_network))]

    if df.empty is True:
        raise ValueError('Read in Phase Velocity CSV is empty. Check that FTAN/Folded is correct')

    avgPhvel = getAvgVelocity(period,df,'ZZ')

    # Edit starting model velocity, edit colorbar
    editBackgroundVel(fmstPath,avgPhvel,colorbar_margin=colorbar_margin)
    # I wrote the above function to adjust the colorbar too but I broke that shit because I'm dumb
    edit_inversion_params(fmstPath,smoothing,damping,lon_grids,lat_grids)

def haversine_distance(lat1: float, lon1: float,
                       lat2: float, lon2: float) -> float:
    """
    Calculates the great circle distance between two points on the Earth
    in kilometers.

    Parameters
    ----------
    lat1 : float
        Latitude of point 1.
    lon1 : float
        Longitude of point 1.
    lat2 : float
        Latitude of point 2.
    lon2 : float
        Longitude of point 2.

    Returns
    -------
    distance : float
        Distance in kilometers.

    """

    R = 6371 # Earth's radius in kilometers

    # Convert decimal degrees to radians
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])

    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * atan2(sqrt(a),sqrt(1-a))
    distance = R * c

    return distance

def get_removed_paths(indices: list,
                      ftanDir: str) -> list:
    """
    Parameters
    ----------
    indices : list
        List of the indices for pairs that were thrown out.
    ftanDir : str
        FTAN home directory.

    Returns
    -------
    pairs : list
        List of lists, where each sublist is lat1,lon1,lat2,lon2.

    """
    distfile = f'{ftanDir}/distances.dat'
    pairs = []
    with open(distfile,'r') as dists:
        for i, line in enumerate(dists.readlines()):
            if i in indices:
                val = line.split(' ')
                lat1 = float(val[1])
                lon1 = float(val[2])
                lat2 = float(val[3])
                lon2 = float(str(val[4])[:-2])
                pairs.append([lat1,lon1,lat2,lon2])

    return pairs

def read_velocity_grid(fmstPeriodDir: str,
                       input_grid: bool=False):
    """
    Reads the velocity grid files from FMST. Can pull either a list of all the
    velocities if input_grid is False, or the just the background velocity used
    as the a priori model if input_grid is set True.

    Parameters
    ----------
    fmstPeriodDir : str
        Path to the FMST directory for the specific period of interest
    input_grid : bool, optional
        If True, returns the a priori background vel. If false, returns a list
        of all model parameters (velocities). The default is False.

    Returns
    -------
    float or list[float]
        If input_grid is False, returns a list of all velocities in the grid.
        If input_grid is True, returns the a priori background velocity.
    """
    vels = []

    if input_grid is False:
        velgrid = f'{fmstPeriodDir}/gridc.vtx'

        with open(velgrid,'r',encoding='utf-8') as gridfile:
            lines = gridfile.readlines()
            for i, line in enumerate(lines):
                if i < 3:
                    continue
                vel = float(line.strip())
                vels.append(vel)

        return vels

    else:
        velgrid = f'{fmstPeriodDir}/gridi.vtx'

        with open(velgrid,'r',encoding='utf-8') as gridfile:
            lines = gridfile.readlines()
            for i, line in enumerate(lines):
                if i != 4:
                    continue
                vel = float(line.split(' ')[2])
                return vel

def get_model_variability(fmstPeriodDir: str) -> float:
    """
    This function returns the variability of the model. Each model parameter is
    the grid cell's percentage different from the a priori background velocity.
    Ex. If the background velocity is 3 km/s, and the cell has a velocity of
    3.3km/s, the percent would be 0.10. ALl of these model parameters are squared
    and summed as a way of measuring how big the variations in the model are.

    Parameters
    ----------
    fmstPeriodDir : str
        Path to the FMST directory for the specific period of interest

    Returns
    -------
    float
        Variability of the model.

    """
    model_velocities = read_velocity_grid(fmstPeriodDir,input_grid=False)
    background_vel = read_velocity_grid(fmstPeriodDir,input_grid=True)

    percentages = [((x-background_vel)/background_vel) for x in model_velocities]

    variability = 0
    for percent in percentages:
        variability += (percent**2)

    return variability

def plot_tradeoff_curve(outputs_file):
    cols = ['Smoothing','Damping','Lon_Grids','Lat_Grids','RMS','Variance','Variability', 'Variance_KMS','Roughness']
    df = pd.read_csv(outputs_file,names=cols,sep='\s+')

    variance = df['Variance'].to_list()
    misfit = df['RMS'].to_list()
    damping = df['Damping'].to_list()

    norm_variance = [x / max(variance) for x in variance]
    norm_misfit = [x / max(misfit) for x in misfit]

    fig, ax = plt.subplots()
    ax.plot(norm_variance,norm_misfit,'o')
    ax.set_xlabel('Normalized Variance')
    ax.set_ylabel('Normalized Misfit')
    ax.set_title('Misfit vs. Variance, Damping')

    plt.show()

def get_grid_spacing(bound_box,lon_grids,lat_grids):
    lon_deg = abs(bound_box[3] - bound_box[2]) / lon_grids
    lat_deg = abs(bound_box[1] - bound_box[0]) / lat_grids
    lon_km = lon_deg * 111
    lat_km = lat_deg * 111

    print(f'=====================STATION SPACING=============================')
    print(f'Longitude Spacing | Degrees: {lon_deg}, Kilometers: ~{lon_km}')
    print(f'Latitude Spacing | Degrees: {lat_deg}, Kilometers: ~{lat_km}')








