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
from math import ceil

import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
from glob import glob
import pandas as pd
from tqdm import tqdm

#pylint: disable=import-error
from swtUtils import ftan, fmstUtils, mapping


def main(
        foldTraces: bool=True,
        create_station_pairs_df: bool=True,
        runFTAN: bool=True,
        makeFMSTInputs: bool=True,
        setupFMSTDirectory: bool=True,
        runCheckerboardTests: bool=True,
        runInversion: bool=True,
        runOnlyGMT: bool=False,
        plot_tradeoff_curves: bool=True,
        plot_inversion_params=True,
        collect_plots=False):
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
    check_specific_station_pair = True

    component='ZZ' # Component pair
    projectCode='Rainier' # Project code
    min_per = 5 # Minimum period
    max_per = 8 # Maximum period
    dt = 1 # Spacing between period measurements
    whole_periods = np.arange(min_per,max_per+1,1) # Whole number periods, if dt is not 1
    periods = list(np.arange(min_per,max_per+dt,dt)) # Periods with a spacing of dt
    periods = [round(dt*np.floor(round(x / dt,2)),1) for x in periods] # Fix floating point errors
    network = 'UW,CC,XU,XD,TA,YH,YW' # Network list
    channel='BH*,HH*,EH*' # Channel list
    bound_box = [46.0796,47.8242,-122.8625,-120.25] # [min_lat,max_lat,min_lon,max_lon]
    snr = 5 # Signal to noise ratio. Can be list if testing multiple
    min_wavelengths = 2 # Minimum wavelength for measurements to be accepted. Can be list if testing multiple

    stationList = fmstUtils.getLocalStations(dataDirectory,'ZZ',forceOverwrite=False)
    print('Acquiring list of valid stations from IRIS...')
    stationDict = fmstUtils.getValidStations(network,bound_box,channel,stationList)
    print('Stations acquired')

    if makeFMSTInputs is True:
        create_station_pairs_df = True

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

    if create_station_pairs_df:
        all_files = glob(f'{ftanDirectory}/Folded/*2_DISP.1')
        stat_pairs = []

        for file in all_files:
            stat1 = file.split('/')[-1].split('_')[0]
            stat2 = file.split('/')[-1].split('_')[1]
            lat1 = stationDict[stat1][0]
            lon1 = stationDict[stat1][1]
            lat2 = stationDict[stat2][0]
            lon2 = stationDict[stat2][1]
            stat_pairs.append((stat1,lat1,lon1,stat2,lat2,lon2))

        df = pd.DataFrame(stat_pairs,columns=['Station1','Station1_Lat',
                                              'Station1_Lon','Station2',
                                              'Station2_Lat','Station2_Lon'])

        period_columns = [f'{i}s' for i in periods]
        for col in period_columns:
            df[col] = np.nan

    if makeFMSTInputs is True:
        """
        This creates properly formatted input files for FMST which are used
        later to set up a full FMST directory.

        """
        SNRs_To_Test = [5]
        plot_travel_times = True

        for _snr in SNRs_To_Test:
            for i,period in enumerate(periods):
                if i == 0:
                    writeDists=True
                else:
                    writeDists=False
                print(f'=====Working on FMST Outputs for {period}s...=====')
                #print(f'{ftanDirectory}/Folded_SNR{_snr}_WL{min_wavelengths}')
                phvels, dists, stat1_list,stat2_list = fmstUtils.makeFMSTInputs(
                                        stationDict=stationDict,
                                        dataDirectory=dataDirectory,
                                        FTANDirectory=f'{ftanDirectory}/Folded_SNR{_snr}_WL{min_wavelengths}',
                                        period=period,
                                        component=component,
                                        minSNR=_snr,
                                        minWavelengths=min_wavelengths,
                                        writeDists=writeDists,
                                        detailedError=True)

                for i, phvel in enumerate(phvels):
                    stat1 = stat1_list[i]
                    stat2 = stat2_list[i]

                    mask = ((df['Station1'] == stat1) & (df['Station2'] == stat2)) | \
                           ((df['Station1'] == stat2) & (df['Station1'] == stat1))

                    if mask.any():
                        df.loc[mask,f'{period}s'] = phvel

                fmstUtils.saveObj(phvels,f'{dataDirectory}/Tomography/{component}/{period}s_allPhvel.pkl')

                avgvel = np.mean(phvels)

                if plot_travel_times:
                    fig, ax = plt.subplots()
                    travel_times = [a/b for a,b in zip(dists,phvels)]
                    ax.scatter(dists,travel_times)
                    ref_vel = 3.113
                    ref_dists = np.arange(0,ceil(max(dists)),1)
                    ref_times = ref_dists / ref_vel
                    ax.plot(ref_dists,ref_times,'g-')

                    plt.title(f'Measured Travel Times for {period}s')
                    ax.set_xlabel('Distance')
                    ax.set_ylabel('Travel Time (s)')
                    plt.show()

                plt.figure()
                plt.hist(phvels,bins=np.arange(1.5,5,0.2))
                plt.title(f'{period}s Rayleigh Wave Phase Velocities: {avgvel} km/s Average')
                plt.xlabel('Phase Velocity (km/s)')
                plt.savefig(f'{dataDirectory}/Figures/{period}s_PhaseDistribution_SNR_{snr}.png')

                avgPhvel = round(float(np.mean(phvels)),4)
                tomoDirectory = fmstUtils.getTomoDirectory(dataDirectory,component) + f'/{period}s'
                fmstUtils.saveObj(avgPhvel,f'{tomoDirectory}/avgPhvel.pkl')

                #print(fmstUtils.loadObj(f'{dataDirectory}/Tomography/ZZ/{period}s/interpErrorDict.pkl'))
                #print(fmstUtils.loadObj(f'{dataDirectory}/Tomography/ZZ/{period}s/issueDict.pkl'))
                #print(fmstUtils.loadObj(f'{dataDirectory}/Tomography/ZZ/{period}s/fpDict.pkl'))
                print(' ')

        df.round(4)
        df.to_csv(f'{ftanDirectory}/PhaseVelocities_SNR{snr}_WL{min_wavelengths}.csv')

    if check_specific_station_pair is True:

        df = pd.read_csv('PhaseVelocities.csv')

        stat1 = 'TA.C05A'
        stat2 = 'TA.D05A'
        mask = ((df['Station1'] == stat1) & (df['Station2'] == stat2)) | \
               ((df['Station1'] == stat2) & (df['Station1'] == stat1))

    damping_list = [20]
    smoothing_list = [0.5]


    lon_grids = 12
    lat_grids = 8
    plot_histograms = False
    stds = 2 # If the residual is this many standard deviations away, the
    # corresponding observed travel time won't be used in the 2nd inversion step

    colorbar_margin=0.4

    percent_removed = (1 - (norm.cdf(stds) - norm.cdf(-1*stds)))*100

    if setupFMSTDirectory is True:
        fmstUtils.create_output_files(fmstDirectory,periods)

    if runCheckerboardTests is True:
        pass

    for damping in damping_list:
        for smoothing in smoothing_list:
            if setupFMSTDirectory is True:
                """
                This uses the previously created travel time inputs as well as a "master"
                FMST folder that contains files used for the inversion for every period,
                to ensure that inversion parameters (except for the actual measurements)
                are used consistently in all the phase inversions.
                """
                for period in periods:
                    fmstUtils.reset_FMST_directory(dataDirectory=dataDirectory,
                                                   fmstDirectory=fmstDirectory,
                                                   ftanDirectory=ftanDirectory,
                                                   component=component,
                                                   period=period,
                                                   projectCode=projectCode,
                                                   snr=snr,
                                                   min_wavelengths=min_wavelengths,
                                                   damping=damping,
                                                   smoothing=smoothing,
                                                   lon_grids=lon_grids,
                                                   lat_grids=lat_grids)

            if runInversion is True:
                """
                This runs the inversions! See FMST documentation for more.
                """
                for i,period in enumerate(periods):

                    fmstPeriodDir = f'{fmstDirectory}/{projectCode}_{period}s_{component}'
                    #fmstPeriodDir='/Users/thomaslee/fmst_v1.1/5.0_CheckerboardTest'
                    print(f'==========PERFORMING INVERSION FOR {period}s...===========')
                    print('Running first inversion pass...     ')
                    inversion_start = time.perf_counter()
                    mkmodelDir = fmstPeriodDir + '/mkmodel'
                    subprocess.run('grid2dss',cwd=mkmodelDir,shell=True,check=False)
                    shutil.copy(f'{mkmodelDir}/grid2d.vtx',f'{fmstPeriodDir}/gridi.vtx')
                    subprocess.run('ttomoss',cwd=fmstPeriodDir,check=False)
                    # We need to run this to get the bounds file for plotting removed raypaths
                    gmtplotDir = fmstPeriodDir + '/gmtplot'
                    subprocess.run('tslicess',cwd=gmtplotDir,check=False)

                    rerun_inversion = False

                    if rerun_inversion:

                        print('Removing worst fitting measurements     ')
                        # Calculate residuals and pull them into a DataFrame
                        rtravel_df, otimes_df, residuals_df = fmstUtils.make_residuals_dfs(ftanDir=ftanDirectory,
                                                                                           fmstDir=fmstPeriodDir)
                        if plot_histograms:
                            fig, ax = fmstUtils.plot_residuals(residuals_df=residuals_df,
                                                               title='Before',
                                                               percent=percent_removed,
                                                               period=period)
                            plt.show()

                        # Calculate the standard deviation of residuals
                        std_residuals = fmstUtils.get_residual_statistics(residuals_df)


                        rms1, var1 = fmstUtils.get_output_info(fmstPeriodDir)
                        print(f'Model RMS: {rms1} Variance: {var1}')

                        # Remove the worst ones and return a new otimes DataFrame
                        otimes_df, indices = fmstUtils.remove_worst_fits(otimes_df,
                                                                residuals_df,
                                                                std_residuals,
                                                                stds=stds)

                        raypaths = fmstUtils.get_removed_paths(indices,ftanDirectory)

                        fig = mapping.plot_removed_ray_paths(fmstPeriodDir=fmstPeriodDir,
                                                             raypaths=raypaths,
                                                             period=period)
                        fig.show()

                        print('Worst fitting measurements removed...')

                        print('Resetting FMST directory...')
                        fmstUtils.reset_FMST_directory(dataDirectory=dataDirectory,
                                                       fmstDirectory=fmstDirectory,
                                                       ftanDirectory=ftanDirectory,
                                                       component=component,
                                                       period=period,
                                                       projectCode=projectCode,
                                                       snr=snr,
                                                       min_wavelengths=min_wavelengths,
                                                       damping=damping,
                                                       smoothing=smoothing,
                                                       lon_grids=lon_grids,
                                                       lat_grids=lat_grids,
                                                       colorbar_margin=colorbar_margin)

                        fmstUtils.rewrite_otimes(fmstPeriodDir,otimes_df)

                        print('Re-running inversion....      ')
                        print('==================================================')
                        subprocess.run('grid2dss',cwd=mkmodelDir,shell=True,check=False)
                        shutil.copy(f'{mkmodelDir}/grid2d.vtx',f'{fmstPeriodDir}/gridi.vtx')
                        subprocess.run('ttomoss',cwd=fmstPeriodDir,check=False)


                        # Re-create dataframes
                        rtravel_df, otimes_df, residuals_df = fmstUtils.make_residuals_dfs(ftanDir=ftanDirectory,
                                                                                           fmstDir=fmstPeriodDir)
                        if plot_histograms:
                            fig, ax = fmstUtils.plot_residuals(residuals_df=residuals_df,
                                                               title='After',
                                                               percent=percent_removed,
                                                               period=period)
                            plt.show()

                        # Calculate the standard deviation of residuals
                        std_residuals = fmstUtils.get_residual_statistics(residuals_df)
                        print(f'Two-step inversion took {(time.perf_counter() - inversion_start):.4f} seconds')

                        # Run misfitss to get the variance and roughness
                        result = subprocess.run('misfitss',cwd=fmstPeriodDir,check=False,capture_output=True,text=True)
                        output_lines = result.stdout.splitlines()
                        assert len(output_lines) == 2, f'Expected 2 lines, got {len(output_lines)}'
                        print('MISFITSSSS TEST OUTPUT')
                        try:
                            variance_kms = float(output_lines[0].split(' ')[-1])
                            roughness = float(output_lines[1].split(' ')[-1])
                            print('==========================')

                            # Write inversion statistics to output file
                            rms2, var2 = fmstUtils.get_output_info(fmstPeriodDir)
                            print(f'Model RMS: {rms2} Variance: {var2}')
                            print(f'RMS Reduction after Removing Worst: {(float(rms1)-float(rms2)):.4f}')
                            print(f'Variance Reduction After Removing Worst: {(float(var1)-float(var2)):.4f}')
                            variability = fmstUtils.get_model_variability(fmstPeriodDir)
                            write_str = f'{float(smoothing)} {float(damping)} {lon_grids} {lat_grids} {rms2} {var2} {variability} {variance_kms} {roughness}\n'

                            with open (f'{fmstDirectory}/{period}s_Outputs','a') as outfile:
                                outfile.write(write_str)

                        except:
                            pass

                    # Plot the final results
                    gmtplotDir = fmstPeriodDir + '/gmtplot'
                    subprocess.run('tslicess',cwd=gmtplotDir,check=False)
                    os.chmod(f'{gmtplotDir}/mapping.py', 0o755)
                    subprocess.run(['python',f'{gmtplotDir}/mapping.py'],
                                   cwd=gmtplotDir,
                                   check=True)

    if plot_tradeoff_curves is True:
        period = 5.0
        outputs_file = f'{fmstDirectory}/{period}s_Outputs'

        fmstUtils.plot_tradeoff_curve(outputs_file)

    if plot_inversion_params is True:
        fmstPeriodDir = f'{fmstDirectory}/{projectCode}_5.0s_{component}'
        im = fmstUtils.plot_smoothing_damping(f'{fmstDirectory}/5.0s_Outputs',smoothing_list,damping_list)

        #rtravel_df, otimes_df, residuals_df = fmstUtils.make_residuals_dfs(fmstPeriodDir)

    if collect_plots is True:

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

    """
    fmstUtils.findAllFinalTomoImages(fmstPath=fmstDirectory,
                                     projectCode=projectCode,
                                     component=component,
                                     periods=periods)
    """

if __name__ == '__main__':
    main(foldTraces=True,
         runFTAN=False,
         create_station_pairs_df=False,
         makeFMSTInputs=False,
         setupFMSTDirectory=False,
         runCheckerboardTests=False,
         runInversion=False,
         runOnlyGMT=False,
         plot_tradeoff_curves=False,
         plot_inversion_params=False,
         collect_plots=False)

