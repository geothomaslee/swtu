# Surface Wave Tomography Utils
#### A collection of scripts for making phase velocity maps from ambient-noise cross-correlations, using FTAN and FMST

Primary scripts used:
Automatic Frequency-Time Analysis (FTAN) from Bensen et. al., 2007
Bensen, G. D., Ritzwoller, M. H., Barmin, M. P., Levshin, A. L., Lin, F., Moschetti, M. P., ... & Yang, Y. (2007). Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements. Geophysical journal international, 169(3), 1239-1260.

Fast Marching Surface Tomography (FMST) from Rawlinson and Sambridge, 2005
Rawlinson, N. and Sambridge M., 2005. "The fast marching method: An effective tool for tomographic imaging and tracking multiple phases in complex layered media", Explor. Geophys., 36, 341-350.


### Setup
First, you'll need to install both FMST and FTAN from their respective original authors. You will also need to compile them, and instructions for doing this can be found from the original authors as well.
It is assumed that these will both be installed in their own directories and compiled there. I installed them in my own home directory. The primary directory for each program will be called its home directory for the rest of the documentation.
`/Users/{your_username}/FTAN` - example FTAN home directory
`Users/{your_username}/fmst_v1.1` - example FMST home directory

<summary>A Note on Compiling: </summary>
<details>
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
