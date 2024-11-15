# Surface Wave Tomography Utils
#### A collection of scripts for making phase velocity maps from ambient-noise cross-correlations, using FTAN and FMST

Primary scripts used:
Automatic Frequency-Time Analysis (FTAN) from Bensen et. al., 2007
Bensen, G. D., Ritzwoller, M. H., Barmin, M. P., Levshin, A. L., Lin, F., Moschetti, M. P., ... & Yang, Y. (2007). Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements. Geophysical journal international, 169(3), 1239-1260.

Fast Marching Surface Tomography (FMST) from Rawlinson and Sambridge, 2005
Rawlinson, N. and Sambridge M., 2005. "The fast marching method: An effective tool for tomographic imaging and tracking multiple phases in complex layered media", Explor. Geophys., 36, 341-350.

## Setup
First, you'll need to install both FMST and FTAN from their respective original authors. You will also need to compile them, and instructions for doing this can be found from the original authors as well.
It is assumed that these will both be installed in their own directories and compiled there. I installed them in my own home directory. The primary directory for each program will be called its home directory for the rest of the documentation.
- `/Users/{your_username}/FTAN` - example FTAN home directory
- `/Users/{your_username}/fmst_v1.1` - example FMST home directory

<details>
<summary>A Note on Compiling (READ ESPECIALLY IF ON MAC): </summary>
- FMST is done entirely with fortran, and the compileall script it comes with defaults to the ifort compiler which works fine if you have it, but I would assume in this day and age most of us already have and use gfortran.
- The makefile in FTAN/bin defaults to using gfortran for its compiler so that likely does not need to be changed
- FTAN has some of its scripts done in C, and therefore a C compiler is needed, which defaults to gcc, which you should have if you have gfortran already
- A NOTE FOR MAC USERS:
  - MACs ship with the Clang compiler for C, and it force aliases the phrase "gcc" to call Clang.
  - Clang cannot interpret most of the flags for gcc that FTAN's makefile wants
  - If you change the compiler to gcc, that alias will force it to call Clang
    - As far as I know, this cannot be overwritten
  - To get around this, you have to set the compiler to the specific gcc version you have. For me, that is gcc-13.
  - MACs also do not have /usr/local/lib on the Path, and therefore the -l flag for loading libraries won't be able to find fftw3, a required input library
    - You can either add it to the path via ~/.bash_profile OR you can use the -L flag in the makefile in FTAN/bin to tell it to load libraries from that directory
</details>

#### File Structure
This program is meant to work with my package ambient2-tlee-fork (a fork of ambient2 by Ross Maguire).
`dataDirectory` refers to the directory containing ambient2runs and the folder `/Stacks`, 
- `/Stacks` contains `/Stacks/{component}` : ex. `/Stacks/ZZ`
  - Each component directory contains a list of two-sided cross-correlation functions which need 5 fields defined: `stla` `stlo` `evla` `evlo` `dist`
  - It doesn't matter which station is defined as the station or event, these will be folded anyway

## FTAN
#### Fold Traces
A workflow is entirely laid out in `swtu.main()`
First, the traces need to be folded. This is done by setting `foldTraces=True` in `main()`.
Foldtraces takes the two sided traces, folds them, and saves them as a new trace with relevant headers pulled from the unfolded trace. These are saved in `/Stacks/{component}/Folded`

#### Run FTAN
First, the folded cross-correlations are moved into the FTAN home directory into `/Folded`. If the folded directory already exists and contains cross-correlations, it is removed entirely and replaced with the ones saved in `dataDirectory/Stacks/{component}/Folded`. Next the subprocess module is used to call `runFTAN.csh`. NOTE: Recently this script has stopped working when called from within Python, so the best practice is to just manually run that script from your own terminal. That script contains its own documentation, but it should output 3 files for every traces.
1. A .txt containing SNR information, that is used by the script to throw out low SNRs before even performing FTAN.
2. A file ending in `.sac_1_DISP.1` which is the initial results, and should not be used.
3. A file ending in `.sac_2_DISP.1` which contains the final FTAN results and should be used.

#### Create Station Pairs DF
This will create an empty DataFrame that contains a row for every possible station pair, each stations coordinates, and a column for every period of interest that we intend to get FTAN values for. The actual results for this will be stored into the DataFrame in the next step.

## FMST
#### Make FMST Inputs
This is handled in `main()` by `makeFMSTInputs`. The documentation for that function is extensive and should be used for detailed reference. Generally, it makes and formats the `receivers.dat` `sources.dat` and `otimes.dat` files that are used by FMST for the inversion. `otimes.dat` specifically contains the travel times of the surface waves and is calculated using the distance between the stations and the phase velocity of the period of interest. 
This is done for every period given by `periods`, defined at the top of `main()` just below the parameters. These files are saved in a directory called `dataDirectory/Tomography/{component}/{period}s`. All of the phase velocities found by FTAN for each station pair are also written into the DataFrame created in the last step, and the DataFrame is saved into the FTAN directory as a csv for future reference.

#### Setup FMST Directory 
All periods we want to perform the inversion for need their own directory, but they should be using the same inversion parameters (obviously different travel times). You first need to setup an master template for these inside the FMST home directory, called `{projectCode}_Master`. The project code just refers to the general name of the project, ex. 'Rainier' in my case.

The master directory should include the following files, which will likely be obtained from the two examples that come with FMST.
- `{projectCode}_Master`
  - fm2dss.in
  - misfitss.in
  - residualss.in
  - subinvss.in
  - subiter.in
  - ttomoss.in
  - gmtplot
    - plotgmt6
    - studybounds.dat
    - tslicess.in
    - velgradabs.cpt
    - velgradproj.cpt
  - mkmodel
    - grid2dss.in
    - subiter.in
    - residuals.dat
  - mkdata
    - synthtss.in
   
Next, a new directory will be created with the name `{projectCode}_{period}s_{component}` and all files from the master will be copied into it. The background velocity set in `grid2dss.in` will be changed to the average velocity for that period, which was pickled and saved while setting up the input files. Finally, the FMST input files will for each period will be copied into their proper FMST directory.

#### Run Inversion
The inversion will be run according to the steps outlined in the documentation for FSMT. First, an input grid is created by running `grid2dss` and moving the output file into the main run directory and renaming it to `gridi.vtx`. Then the inversion is actually performed using `ttomoss`. `tslicess` is run in `/gmtplot` to create an input grid formatted for GMT, then `plotgmt6` is called. `plotgmt6` is a script written by Wenkai Song and heavily edited by me, that just adjusts the original `plotgmt` file that comes with FMST to work with GMT6.5.0. My version of `plotgmt6` contains all the same things as the original, but with a couple additions of things relevant to me - Ex. Mount Rainier is plotted. If you want to play around with these custom plotting features, always change the master plotgmt6 file so that it is applied to all periods. All the phase velocity maps are moved into `{fmstHomeDirectory}/phvelMaps` and renamed to their corresponding period.

An output file is made in the main FMST directory called `{period}s_outputs` and for every run of the inversion (after the two-step inversion if two-step is running), a list of information about the runs is output to compare different inversion parameters.
The columns in that file correspond to `smoothing, damping, number of longitudinal grid cells, number of latitudinal grid cells, rms, variance`
This inversion can also be set to run multiple times with different combinations of damping and smoothing, to find the ideal parameters.

#### Two-Step Inversion
This is currently hard coded into the previous step, but may become an optional event in the near future. The residuals between the predicted travel time and the observed travel time are calculated. The observed travel times come from `otimes.dat` and the predicted travel times come from `rtravel.dat`. Next, residuals that are more than `stds` standard deviations away from 0 are found. Their corresponding observed travel times are then identified, and set to be unused in `otimes.dat`. Finally, the period specific FMST directory is reset using `fmstUtils.reset_FMST_directory()` and the updated `otimes.dat` is written into the directory. The inversion is then run again with otherwise identical parameters. Essentially, this just removes a small number of measurements that the model is having the hardest time fitting.
