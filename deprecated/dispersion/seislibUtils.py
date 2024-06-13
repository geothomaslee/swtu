#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 10:03:38 2024

@author: Thomas Lee
University of New Mexico
Department of Earth and Planetary Science
"""

import os
import pickle
from glob import glob

from seislib.tomography import SeismicTomography
from seislib.plotting import make_colorbar
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import numpy as np
from tqdm import tqdm

from . import ftan as ft

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

def getTomoDirectory(dataDirectory,component,period):
    """Returns tomography working directory"""
    if not isinstance(component, str):
        raise TypeError('Component must be 2 character string. Ex. "ZZ"')
    if len(component) != 2:
        raise ValueError('Component must be 2 character string. Ex. "ZZ"')

    tomoDirectory = os.path.expanduser(dataDirectory + f'/Tomography/{component}/{period}s')

    if not os.path.isdir(dataDirectory + '/Tomography'):
        os.mkdir(dataDirectory + '/Tomography')

    if not os.path.isdir(dataDirectory + f'/Tomography/{component}'):
        os.mkdir(dataDirectory + f'/Tomography/{component}')

    if not os.path.isdir(tomoDirectory):
        os.mkdir(tomoDirectory)

    return tomoDirectory

def saveObj(obj, filename):
    """Quick function for pickling a file"""
    with open(filename, 'wb') as f:
        pickle.dump(obj, f)

def loadObj(filename):
    """Quick function for pickling a file"""
    with open(filename, 'rb') as f:
        return pickle.load(f)

def getSACFileList(componentDirectory):
    """Returns a list of SAC files in the component directory"""
    fileList = glob(componentDirectory + '/*.sac')
    return fileList

def getStationNames(filepath):
    """Pulls the name of the two stations from the name of the file"""
    basename = os.path.basename(filepath)
    name_no_sac = basename[:-4]
    stat1, stat2 = name_no_sac.split('_')
    return stat1, stat2

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
                       'could not find left or right' : 0,
                       'unreasonable phvel from interp' : 0,
                       'unreasonable snr' : 0,
                       'bad left phvel' : 0,
                       'bad right phvel' : 0}
    return interpErrorDict

def makeSeisLibInput(dataDirectory,component,period,minSNR,minWavelengths,veltype='p',detailedError=False,rerunIfExists=False):
    """For making the inputs for seislib"""
    componentDirectory = getComponentDirectory(dataDirectory,component)
    tomoDirectory = getTomoDirectory(dataDirectory,component,period)
    fileList = getSACFileList(componentDirectory)
    if len(fileList) == 0:
        raise ValueError('Could not find any files')

    # Setting up output file
    SeisLibInputFile = tomoDirectory + '/seisLibInput.txt'
    if os.path.isfile(SeisLibInputFile) and rerunIfExists is False:
        return
    open(SeisLibInputFile, 'w',encoding='utf-8').close()

    # Creating issue dicts
    issue_dict = makeIssueDict()
    fpDict = makefpDict()
    interpErrorDict = makeInterpErrorDict()

    with open(SeisLibInputFile, "a",encoding='utf-8') as outfile:
        for file in tqdm(fileList):
            # Remove autocorrelations
            stat1, stat2 = getStationNames(file)
            if stat1 == stat2:
                issue_dict['autocorrelation'] += 1
                continue

            # Perform FTAN. If error handler returns none, then add to counter,
            # but error handler will put error type in fpDict
            fparam_out = pt.FTAN(file,tmin=1,tmax=40,nfin=41)
            fparam_out, fpDict = _fparamErrorCodeHandler(fparam_out,fpDict)
            if fparam_out is None:
                issue_dict['fparam_out_none'] += 1
                continue

            # Pull relevant info
            fparam = fparam_out[0]
            stats = fparam_out[1]
            obper1, gvel1, phvel1, snr = pt.pullRelevantFTANInfo(fparam)

            # Choose group or phase velocity
            if veltype == 'p':
                vel = phvel1
            elif veltype == 'g':
                vel = gvel1
            else:
                return ValueError('veltype must be "g" or "p"')

            # Interpolate the period/velocity curve and put any errors inside
            # the interpErrorDict
            interpOut = interpPeriod(period,obper1,vel,snr)
            interpOut, interpErrorDict = _interpPeriodErrorHandler(interpOut,interpErrorDict)
            if interpOut is None:
                issue_dict['interpOut none'] += 1
                continue

            vel = interpOut[0]
            snr = interpOut[1]

            # Remove unreasonable velocities
            if vel < 1.5 or vel > 5:
                issue_dict['bad phvel'] += 1
                continue

            # Remove low SNR
            if snr < minSNR:
                issue_dict['low snr'] += 1
                continue

            # Remove stations that are too close for the period of interest
            minDist = minWavelengths * getWavelength(period,vel)
            if stats.sac['dist'] < minDist:
                issue_dict['too close'] += 1
                continue

            # Write stations/velocity to the output file
            lat1 = stats.sac['stla']
            lon1 = stats.sac['stlo']
            lat2 = stats.sac['evla']
            lon2 = stats.sac['evlo']
            outfile.write(f'{lat1} {lon1} {lat2} {lon2} {vel}\n')
            issue_dict['good'] += 1

    if detailedError is True:
        saveObj(issue_dict, f'{tomoDirectory}/issueDict.pkl')
        saveObj(fpDict, f'{tomoDirectory}/fpDict.pkl')
        saveObj(interpErrorDict,f'{tomoDirectory}/interpErrorDict.pkl')

    return None

def _fparamErrorCodeHandler(fparam_out,fpDict):
    """Counts the types of errors from fparam"""
    if isinstance(fparam_out,tuple):
        if fparam_out[0] is None:
            fpDict['fparam returning none'] += 1
            return None, fpDict
        return fparam_out, fpDict

    if fparam_out == 1:
        fpDict['infDepen'] += 1
    if fparam_out == 2:
        fpDict['omdom'] += 1
    if fparam_out == 3:
        fpDict['zero'] += 1
    if fparam_out == 4:
        fpDict['slice'] += 1
    if fparam_out == 5:
        fpDict['cannot unpack none'] += 1

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

def getWavelength(per,vel):
    """Returns the wavelength given period and velocity"""
    return per * vel

def runTomo(dataPath, rdamp,hc=25,cs=0.04):
    """Runs seislib tomography"""
    tomo = SeismicTomography(cell_size=cs,
                             regular_grid=False,
                             #latmin=box[0],
                             #latmax=box[1],
                             #lonmin=-box[2],
                             #lonmax=-box[3],
                             verbose=True)
    tomo.add_data(src=dataPath)

    tomo.grid.set_boundaries(latmin=tomo.latmin_data,
                             latmax=tomo.latmax_data,
                             lonmin=tomo.lonmin_data,
                             lonmax=tomo.lonmax_data)
    tomo.compile_coefficients(keep_empty_cells=False)
    tomo.refine_parameterization(hitcounts=hc, keep_empty_cells=True)
    tomo.refine_parameterization(hitcounts=hc, keep_empty_cells=True)
    tomo.refine_parameterization(hitcounts=hc, keep_empty_cells=True)
    c = 1 / tomo.solve(rdamp=rdamp)
    return tomo, c

def plot_tomo(tomo, c, period, box = None,fig_size=(5, 7), cm='jet_r',veltype='g'):
    """Plots tomography outputs"""
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Miller())
    ax.coastlines(resolution='10m', color='k', lw=1, zorder=100)

    # Add gridlines and labels with correct formatters
    gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = True  # Enable top labels for longitude
    gl.right_labels = False  # Disable right labels for latitude
    gl.bottom_labels = False  # Optionally, disable bottom labels if only top labels are desired
    gl.left_labels = True  # Enable left labels for latitude
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12, 'color': 'gray'}
    gl.ylabel_style = {'size': 12, 'color': 'gray'}
    img = tomo.colormesh(mesh=tomo.grid.mesh,
                         c=c,
                         ax=ax,
                         cmap=cm,
                         shading='flat',
                         edgecolors='face')

    if box is None:
        map_boundaries = (tomo.grid.lonmin, tomo.grid.lonmax,
                          tomo.grid.latmin, tomo.grid.latmax)

        print(map_boundaries)
    else:
        map_boundaries = (box[2],box[3],box[0],box[1])

    ax.set_extent(map_boundaries, ccrs.PlateCarree()) # had a backslash at the end
    if veltype == 'g':
        ax.set_title(f'Group Velocity of {period}s Rayleigh Waves')
    elif veltype == 'p':
        ax.set_title(f'Phase Velocity of {period}s Rayleigh Waves')
    else:
        print(veltype)
        raise ValueError('veltype must be "g" or "p"')

    cb = make_colorbar(ax, img, orientation='horizontal')


    transform = ccrs.PlateCarree()._as_mpl_transform(ax) #pylint: disable=protected-access
    ax.plot(-121.7603,46.8523,'ok',transform=ccrs.PlateCarree())
    ax.annotate('Mount Rainier',
                xy=(-121.7603,46.88),
                xycoords=transform,
                fontsize=12)

    if veltype == 'g':
        cb.set_label(label='Group Velocity [km/s]', labelpad=10, fontsize=10)
    if veltype == 'p':
        cb.set_label(label='Phase Velocity [km/s]', labelpad=10, fontsize=10)

    dataDirectory = '/Volumes/NewHDant/RainierAmbient'
    plt.savefig(fname=f'{dataDirectory}/Tomography/Figures/{period}s.png')
    plt.show()

def runEntireWorkflow(period,veltype='g'):
    """Runs the entire tomography workflow"""
    dataDirectory = '/Volumes/NewHDant/RainierAmbient'
    makeSeisLibInput(dataDirectory,
                     component='ZZ',
                     period=period,
                     minSNR=5,
                     minWavelengths=3,
                     veltype=veltype,
                     detailedError=True,
                     rerunIfExists=False)

    bound_box = [46.1,47.3,-122.5,-120.9]
    tomo, c = runTomo(f'{dataDirectory}/Tomography/ZZ/{period}s/SeisLibInput.txt',
                      rdamp=0.003,
                      hc=5,
                      cs=0.05)

    plot_tomo(tomo, c,
              period=period,
              veltype=veltype)
              #box=bound_box)




testPeriod = 8
runEntireWorkflow(testPeriod)

print(loadObj(f'/Volumes/NewHDant/RainierAmbient/Tomography/ZZ/{testPeriod}s/interpErrorDict.pkl'))
print(loadObj(f'/Volumes/NewHDant/RainierAmbient/Tomography/ZZ/{testPeriod}s/issueDict.pkl'))
print(loadObj(f'/Volumes/NewHDant/RainierAmbient/Tomography/ZZ/{testPeriod}s/fpDict.pkl'))
