#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 08:15:34 2024

@author: Thomas Lee
University of New Mexico
Department of Earth and Planetary Science

Citations:
    [1] Bensen, G. D., et al. (2007). Processing seismic ambient noise data to
    obtain reliable broad-band surface wave dispersion measurements.
        Geophysical Journal International, 169(3), 1239-1260.
    [2] Rawlinson, N. and Sambridge M., 2005. "The fast marching method: An
        effective tool for tomographic imaging and tracking multiple phases in
        complex layered media", Explor. Geophys., 36, 341-350.
"""
#pylint: disable=invalid-name
import subprocess
import shutil
import os
import time

import matplotlib.pyplot as plt
import numpy as np

#pylint: disable=import-error
from swtUtils import ftan, fmstUtils


def main(
        foldTraces: bool=True,
        makeReferenceVelocities: bool=True,
        runFTAN: bool=True,
        makeFMSTInputs: bool=True,
        setupFMSTDirectory: bool=True,
        runInversion: bool=True,
        runOnlyGMT: bool=False):
    """
    The main function of swtu. You should read each of the 6 input parameters
    as their own individual steps, and documentation is included progressively
    throughout the code and in the README.md on the GitHub page for this repo.

    Parameters
    ----------
    foldTraces : bool, optional
        Folds the cross-correlations. The default is True.
    runFTAN : bool, optional
        Performs FTAN on the folded traces. The default is True.
        This is performed on a copy of the folded traces that is moved into the
        FTAN directory, not the clean copies saved in the dataDirectory.
    makeReferenceVelocities: bool, optional
        FMST requires a starting velocity that is used for every cell in the
        inversion grid, and the model is iteratively improved from there. It is
        therefore extremely sensitive to the starting velocity, so this function
        makes reference velocities for every period.
    makeFMSTInputs : bool, optional
        Takes the FTAN outputs and makes properly formatted inputs for FMST, on
        a per period basis. The default is True.
    setupFMSTDirectory : bool, optional
        Creates a properly formatted FMST directory for each period.
        The default is True.
    runInversion : bool, optional
        Runs FMST. The default is True.
    runOnlyGMT : bool, optional
        Will run plotgmt6 without re-running the inversion. The default is False.

    """

    dataDirectory = '/Volumes/NewHDant/RainierAmbient2'
    ftanDirectory = '/Users/thomaslee/FTAN'
    fmstDirectory = '/Users/thomaslee/fmst_v1.1'
    component='ZZ'
    projectCode='Rainier'
    periods = list(np.arange(1,16,1))
    network = 'UW,CC,XU,XD,TA,YH,YW'
    channel='BH*,HH*,EH*'
    bound_box = [46.0796,47.8242,-122.8625,-120.25]

    stationList = fmstUtils.getLocalStations(dataDirectory,'ZZ',forceOverwrite=False)
    print('Acquiring list of valid stations from IRIS...')
    stationDict = fmstUtils.getValidStations(network,bound_box,channel,stationList)
    print('Stations acquired')

    if foldTraces is True:
        """
        First step is to fold all of the two-sided traces
        This helps increase signal recovery, especially in places where the cross-correlations
        are heavily one-sided, such as near the ocean
        These are placed inside FTAN/Folded
        """
        print('======================FOLDING TRACES==========================')
        ftan.foldAllTraces(dataDirectory,component)
        print(' ')

    if runFTAN is True:
        """
        Now we perform FTAN analysis using runFTAN.csh
        This script is a hand-me-down script through several generations of students
        in Brandon Schmandt's group, that I updated the documentation for extensively

        October 25, 2024 Update: This code doesn't work, I don't know why, but
        just run the script directly in your terminal to do it manually.
        """

        print('=====RUNNING FTAN ANALYSIS=====')
        print('Copying folded SAC files to FTAN directory...')
        if os.path.isdir(ftanDirectory +'/Folded'):
            shutil.rmtree(f'{ftanDirectory}/Folded')

        shutil.copytree(src=f'{dataDirectory}/Stacks/{component}/Folded',
                        dst=f'{ftanDirectory}/Folded')

        print('Calling runFTAN.csh through subprocess.call()')
        print('Output from the terminal is not displayed until the end')
        print('Usually takes 2ish minutes to run, but your mileage may vary...')

        subprocess.call(f'{ftanDirectory}/runFTAN.csh',cwd=ftanDirectory)
        print(' ')

    if makeReferenceVelocities is True:
        """
        This creates a reference velocity curve that represents an average for
        the entire region, and is also used as a starting velocity for every cell
        in the tomographic region later on.
        """
        print('=======CREATING REFERENCE VELOCITY CURVE=======')
        fmstUtils.makeTomoDirectory(dataDirectory,periods,component)
        refMinSNR=6
        refMinWavelengths=2
        refPhvelDict = {}
        refperiods = np.arange(0,20,0.2)
        for period in refperiods:
            refPhvelDict[period],phvels = fmstUtils.getReferenceVelocity(
                                            stationDict=stationDict,
                                            dataDirectory=dataDirectory,
                                            FTANDirectory=f'{ftanDirectory}/Folded',
                                            period=period,
                                            component=component,
                                            minSNR=refMinSNR,
                                            minWavelengths=refMinWavelengths)
            #fmstUtils.saveObj(phvels,f'{dataDirectory}/Tomography/{component}/{period}s_refPhvels.pkl')

        fmstUtils.saveObj(refPhvelDict,f'{dataDirectory}/Tomography/{component}/refPhvelCurve.pkl')

        refPhvelDict = fmstUtils.loadObj(f'{dataDirectory}/Tomography/{component}/refPhvelCurve.pkl')

        plt.plot(list(refPhvelDict.keys()),list(refPhvelDict.values()))
        plt.title(f'SNR {refMinSNR} at >{refMinWavelengths} wavelengths')
        plt.xlabel('Period')
        plt.ylabel('Phase Velocity (km/s)')
        plt.show()

    if makeFMSTInputs is True:
        """
        This creates properly formatted input files for FMST which are used
        later to set up a full FMST directory.

        """

        fmstUtils.makeTomoDirectory(dataDirectory,periods,component)
        for period in periods:
            print(f'=====Working on FMST Outputs for {period}s...=====')

            phvels = fmstUtils.makeFMSTInputs(stationDict=stationDict,
                                              dataDirectory=dataDirectory,
                                              FTANDirectory=f'{ftanDirectory}/Folded',
                                              period=period,
                                              component=component,
                                              minSNR=3,
                                              minWavelengths=1.5,
                                              detailedError=True)

            fmstUtils.saveObj(phvels,f'{dataDirectory}/Tomography/{component}/{period}s_allPhvel.pkl')

            plt.figure()
            plt.hist(phvels,bins=np.arange(1.5,5,0.2))
            plt.title(f'{period}s Rayleigh Wave Phase Velocities')
            plt.xlabel('Phase Velocity (km/s)')
            plt.savefig(f'{dataDirectory}/Figures/{period}s_PhaseDistribution.png')

            avgPhvel = round(float(np.mean(phvels)),4)
            tomoDirectory = fmstUtils.getTomoDirectory(dataDirectory,component) + f'/{period}s'
            fmstUtils.saveObj(avgPhvel,f'{tomoDirectory}/avgPhvel.pkl')

            print(fmstUtils.loadObj(f'{dataDirectory}/Tomography/ZZ/{period}s/interpErrorDict.pkl'))
            print(fmstUtils.loadObj(f'{dataDirectory}/Tomography/ZZ/{period}s/issueDict.pkl'))
            print(fmstUtils.loadObj(f'{dataDirectory}/Tomography/ZZ/{period}s/fpDict.pkl'))
            print(' ')

    if setupFMSTDirectory is True:
        """
        This uses the previously created travel time inputs as well as a "master"
        FMST folder that contains files used for the inversion for every period,
        to ensure that inversion parameters (except for the actual measurements)
        are used consistently in all the phase inversions.
        """
        for period in periods:
            tomoDirectory = fmstUtils.getTomoDirectory(dataDirectory,component) + f'/{period}s'
            fmstPath = fmstUtils.setupFTANDirectory(FMSTDirectory=fmstDirectory,
                                                    period=period,
                                                    projectCode=projectCode,
                                                    component=component,
                                                    _overwrite=True)

            fmstUtils.moveFMSTInputs(fmstPath=fmstPath,
                                     tomoDirectory=tomoDirectory,
                                     _overwrite=False)

            avgPhvel = fmstUtils.getAvgVelocity(dataDirectory,period,component)
            fmstUtils.editBackgroundVel(fmstPath,avgPhvel)

    if runInversion is True:
        """
        This runs the inversions!
        """
        for i,period in enumerate(periods):
            fmstDir = f'{fmstDirectory}/{projectCode}_{period}s_{component}'
            if runOnlyGMT is False:
                print(f'=====PERFORMING INVERSION FOR {period}s...======')
                inversion_start = time.perf_counter()
                mkmodelDir = fmstDir + '/mkmodel'
                subprocess.run('grid2dss',cwd=mkmodelDir,shell=True,check=False)
                shutil.copy(f'{mkmodelDir}/grid2d.vtx',f'{fmstDir}/gridi.vtx')

                subprocess.run('ttomoss',cwd=fmstDir,check=False)

                print(f'Inversion took {time.perf_counter() - inversion_start} seconds')

            if runOnlyGMT is True:
                if i == 0:
                    print('=====REBUILDING PHASE VELOCITY MAPS WITH MASTER GMT TEMPLATE======')
                print(f'Making map for {period}s...')
                shutil.copy(src=f'{fmstDirectory}/{projectCode}_Master/gmtplot/plotgmt6',
                            dst=f'{fmstDir}/gmtplot')

            gmtplotDir = fmstDir + '/gmtplot'
            subprocess.run('tslicess',cwd=gmtplotDir,check=False)
            subprocess.run(['chmod','+x', './plotgmt6'],cwd=gmtplotDir,check=False)
            subprocess.call(f'{gmtplotDir}/plotgmt6',cwd=gmtplotDir)

    """
    fmstUtils.findAllFinalTomoImages(fmstPath=fmstDirectory,
                                     projectCode=projectCode,
                                     component=component,
                                     periods=periods)
    """

if __name__ == '__main__':
    main(foldTraces=False,
         runFTAN=False,
         makeReferenceVelocities=False,
         makeFMSTInputs=True,
         setupFMSTDirectory=False,
         runInversion=False,
         runOnlyGMT=False)
