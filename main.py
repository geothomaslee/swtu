#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 08:15:34 2024

@author: Thomas Lee
University of New Mexico
Department of Earth and Planetary Science
"""

import subprocess
import shutil
import os

import matplotlib.pyplot as plt
import numpy as np

import swtUtils.ftan as ftan
import swtUtils.fmstUtils as fmstUtils


def main(foldTraces=True,runFTAN=True,makeFMSTInputs=True,setupFMSTDirectory=True):

    dataDirectory = '/Volumes/NewHDant/RainierAmbient'
    ftanDirectory = '/Users/thomaslee/FTAN'
    fmstDirectory = '/Users/thomaslee/fmst_v1.1'
    component='ZZ'
    projectCode='Rainier'
    periods = list(np.arange(1,20,1))

    if foldTraces is True:
        # First step is to fold all of the two-sided traces
        # This helps increase signal recovery, especially in places where the cross-correlations
        # are heavily one-sided, such as near the ocean
        print('====FOLDING TRACES=====')
        ftan.foldAllTraces(dataDirectory,component)
        print(' ')

    if runFTAN is True:
        # Now we perform FTAN analysis using runFTAN.csh
        # This script is a hand-me-down script through several generations of students
        # in Brandon Schmandt's group, that I updated the documentation for extensively
        # Originally this script would throw out low SNR and stations that it estimated
        # to be too close, but I don't like how it does that so I do that myself later

        print('=====RUNNING FTAN ANALYSIS=====')
        print('Copying folded SAC files to FTAN directory...')
        if os.path.isdir(ftanDirectory +'/Folded'):
            shutil.rmtree(f'{ftanDirectory}/Folded')

        shutil.copytree(src=f'{dataDirectory}/Stacks/{component}/Folded',
                        dst=f'{ftanDirectory}/Folded')

        print('Calling runFTAN.csh through subprocess.call()')
        print('Output from the terminal is not displayed until the end')
        print('Usually takes 2ish minutes to run, but your mileage may vary...')

        output, error = subprocess.call(f'{ftanDirectory}/runFTAN.csh',cwd=ftanDirectory)
        print(output.decode())
        print(error.decode())

        print(' ')

    if makeFMSTInputs is True:
        # The next step is to define the periods we want to make phase vel maps for

        fmstUtils.makeTomoDirectory(dataDirectory,periods,component)
        for period in periods:
            print(f'===Working on FMST Outputs for {period}s...')
            network = 'UW,CC,XU,XD,TA,YH'
            channel='BH*,HH*,EH*'
            bound_box = [46.1,47.3,-122.5,-120.9]

            stationList = fmstUtils.getLocalStations(dataDirectory,'ZZ')
            stationDict = fmstUtils.getValidStations(network,bound_box,channel,stationList)

            phvels = fmstUtils.makeFMSTInputs(stationDict=stationDict,
                                              dataDirectory=dataDirectory,
                                              FTANDirectory=f'{ftanDirectory}/Folded',
                                              period=period,
                                              component=component,
                                              minSNR=3,
                                              minWavelengths=1.5,
                                              detailedError=True)

            plt.hist(phvels,bins=np.arange(1.5,5,0.2))
            plt.title(f'{period}s Rayleigh Wave Phase Velocities')
            plt.xlabel('Phase Velocity (km/s)')
            plt.show()

            print(fmstUtils.loadObj(f'/Volumes/NewHDant/RainierAmbient/Tomography/ZZ/{period}s/interpErrorDict.pkl'))
            print(fmstUtils.loadObj(f'/Volumes/NewHDant/RainierAmbient/Tomography/ZZ/{period}s/issueDict.pkl'))
            print(fmstUtils.loadObj(f'/Volumes/NewHDant/RainierAmbient/Tomography/ZZ/{period}s/fpDict.pkl'))
            print(' ')

    if setupFMSTDirectory is True:
        for period in periods:
            ftanPath = fmstUtils.setupFTANDirectory(FMSTDirectory=fmstDirectory,
                                                    period=period,
                                                    projectCode=projectCode,
                                                    component=component,
                                                    _overwrite=False)


if __name__ == '__main__':
    main(foldTraces=False,
         runFTAN=False,
         makeFMSTInputs=False,
         setupFMSTDirectory=True)
