#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:19:10 2025

@author: Thomas Lee
Rice University
"""

import os
import shutil
from glob import glob

def setupCheckerboardTest(fmstDirectory: str,
                          period: float,
                          projectCode: str='',
                          component: str='ZZ',
                          overwriteMode: bool=True):

    chdir = f'{fmstDirectory}/Checkerboard_{projectCode}_{period}s_{component}'
    master_dir = f'{fmstDirectory}/{projectCode}_Master'
    period_dir = f'{fmstDirectory}/{projectCode}_{period}s_{component}'

    if os.path.isdir(chdir):
        if overwriteMode is False:
            if glob(f'{chdir}/**') != []:
                raise ValueError('Checkerboard test directory already exists and contains files, but overwrite mode set to False\n' +
                    'Set overwrite mode to True (default) or manually sort out if you want to overwrite the previous run')
        else:
            shutil.rmtree(chdir)
            shutil.copytree(master_dir,chdir)
    else:
        shutil.copytree(master_dir,chdir)

    shutil.copy(f'{period_dir}/receivers.dat',chdir)
    shutil.copy(f'{period_dir}/sources.dat',chdir)
    shutil.copy(f'{period_dir}/otimes.dat',chdir)

    """
    with open(f'{chdir}/mkmodel/grid2dss.in','r',encoding='utf-8') as gridfile:
        lines = gridfile.readlines()
        #lines[29][0] = '1'
        #lines[31][0] = '1'
        #gridfile.writelines()
    """









fmstDirectory = '/Users/thomaslee/fmst_v1.1'
setupCheckerboardTest(fmstDirectory=fmstDirectory,
                      period=5.0,
                      projectCode='Rainier',
                      component='ZZ')