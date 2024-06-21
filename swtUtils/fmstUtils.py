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

from tqdm import tqdm
import numpy as np
import obspy
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
from deprecated import deprecated

import swtUtils.ftan as ft #pylint: disable=import-error
#pylint: disable=invalid-name

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
    """Returns the tomography directory for the component of interest"""
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
    client = Client("IRIS")
    inventory = client.get_stations(network = network,
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

def getInventoryLength(inventory):
    """Returns the number of stations in your inventory"""
    length = 0
    for network in inventory:
        length += len(network)

    return length

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

def getReferenceVelocity(stationDict,dataDirectory,FTANDirectory,period,component,minSNR,minWavelengths):
    componentDirectory = getComponentDirectory(dataDirectory,component)
    stationList = list(stationDict.keys())
    interpErrorDict = makeInterpErrorDict()

    phvels = []
    for stat1 in tqdm(stationList):
        for stat2 in stationList:
            if stat1 == stat2:
                continue

            filepath = checkIfFTANExists(stat1,stat2,FTANDirectory)
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
    return avgPhvel

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

    phvels = []
    with open(timesFile, 'a',encoding='utf-8') as outfile:
        #outfile.write(f'{len(stationList)**2}\n')
        for stat1 in tqdm(stationList):
            for stat2 in stationList:
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

@deprecated(reason="Instead of using the avg velocity for all data, it's" +
            "better to use the reference curve made in main using getReferenceVelocity")
def getAvgVelocity(dataDirectory,period,component):
    """Reads the previously pickled average velocity for that period"""
    tomoDirectory = getTomoDirectory(dataDirectory,component) + f'/{period}s'
    avgPhvel = loadObj(tomoDirectory + '/avgPhvel.pkl')
    return avgPhvel

def editBackgroundVel(fmstPath,avgPhvel):
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
