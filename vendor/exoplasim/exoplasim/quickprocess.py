#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 15:41:42 2023

@author: Mo Cohen
"""

## Little script to just use pyburn to postprocess
## raw ExoPlaSim output into a netCDF4 file

## Needs to be run from the top-level exoplasim directory which contains
## pyburn, gcmt, etc.

## The .py scripts inside the top-level exoplasim directory have been 
## hacked by commenting out the exoplasim.pyburn etc. type imports, which
## throw an error when importing in the Python interpreter inside the
## directory. This way I can access the functions in these modules
## without having to reinstall the package.

#import iris
import exoplasim
import numpy as np
from pyburn import *

#out = 'MOST.00010.npz'
#raw = 'MOST.00010'

raw_directory = '/home/s1144983/Repos/exodev/exoplasim/wolf_1262_5e-07_crashed/'

for file_no in np.arange(0,1):
    print(file_no)
    raw = f'MOST.0000{file_no}'
    out = f'MOST.0000{file_no}.npz'
    print(raw, out)

    postprocess(rawfile=raw_directory+raw, outfile=raw_directory+out, radius=1.66, gravity=12.1,
                gascon=295.37, namelist='/home/s1144983/Repos/exodev/exoplasim/plasim/run/example.nl')

