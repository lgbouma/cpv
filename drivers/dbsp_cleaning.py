"""
Given a night's worth of DBSP data, update the DBSP header frametypes to
"bias", "science", "arc", "flat", and "standard", if the latter exists.
"""
from astropy.io import fits
from glob import glob
import os
from os.path import join

##########################################
# vvv change below here vvv

# NIGHT: 2023/11/11
nightdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP/20231111'
bias_indices = range(1, 11+1)
arc_indices = [13,14,15, 85,86]
domeflat_indices = [20,21,22, 87,88] #range(20, 22+1)
science_indices = range(24, 83+1) # 24-79: LP12-502, 80-83, two others
objectname_indices = range(24, 79+1) # 1 to 11
junk_indices = [12, 19, 23, 84]

# NIGHT: 2023/11/12
nightdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP/20231112'
bias_indices = range(9, 18+1)
arc_indices = [1,2,3, 85,86]
domeflat_indices = [4,5,6,7,8,87,88] #range(20, 22+1)
science_indices = range(19, 84+1) # 24-79: LP12-502, 80-83, two others
objectname_indices = range(22, 77+1)
junk_indices = [81,82]

# NIGHT: 2023/11/12
nightdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP/20231207'
bias_indices = range(9, 18+1)
arc_indices = [1,2,3, 84,85]
domeflat_indices = [4,5,6,7,8, 86,87] #range(20, 22+1)
science_indices = range(19, 83+1) # 24-79: LP12-502, 80-83, two others
objectname_indices = range(23, 69+1)
junk_indices = []

# NIGHT: 2023/12/18
nightdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP/20231218'
bias_indices = range(9, 18+1)
arc_indices = [1,2,3, 49, 50, 51]
domeflat_indices = [4,5,6,7,8, 52,53,54] #range(20, 22+1)
science_indices = range(19, 48+1) # 24-79: LP12-502, 80-83, two others
objectname_indices = range(19, 38+1)
junk_indices = []

# NIGHT: 2023/12/18
nightdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP/20240115'
bias_indices = range(7, 16+1)
arc_indices = [1,2,3, 80,81,82]
domeflat_indices = [4,5,6, 78,79] #range(20, 22+1)
science_indices = range(18, 73) # 18-50: LP12-502, 51-73: TIC141
objectname_indices = range(18, 50+1)
junk_indices = [17, 74, 75,76,77]

# NIGHT: 2024/04/28
nightdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP/20240428'
bias_indices = range(7, 16+1)
arc_indices = [1,2,3, 66,67]
domeflat_indices = [4,5,6] #range(20, 22+1)
science_indices = range(17, 65+1) # 18-50: LP12-502, 51-73: TIC141
objectname_indices = range(17, 65+1)
junk_indices = []

#OBJECTNAME = 'LP_12-502'
OBJECTNAME = 'DG_CVn'

DO_BIAS_UPDATE = 1
DO_ARC_UPDATE = 1
DO_SCIENCE_UPDATE = 1
DO_OBJECTNAME_UPDATE = 1
DO_FLAT_UPDATE = 1
DO_JUNK_UPDATE = 1


# ^^^ change above here ^^^
##########################################

# convert indices to zfill'd file names
bias_indices = [f"{str(ix).zfill(4)}" for ix in bias_indices]
arc_indices = [f"{str(ix).zfill(4)}" for ix in arc_indices]
domeflat_indices = [f"{str(ix).zfill(4)}" for ix in domeflat_indices]
science_indices = [f"{str(ix).zfill(4)}" for ix in science_indices]
objectname_indices = [f"{str(ix).zfill(4)}" for ix in objectname_indices]
junk_indices = [f"{str(ix).zfill(4)}" for ix in junk_indices]

# update "IMGTYPE" to be "bias"...
if DO_BIAS_UPDATE:
    for c in ['red','blue']:
        for bias_index in bias_indices:
            fname = f"{c}{bias_index}.fits"
            fitspath = join(nightdir, fname)
            assert os.path.exists(fitspath)
            fits.setval(fitspath, 'IMGTYPE', value='bias')
            print(f'Updated IMGTYPE to "bias" for {fitspath}')

# update "IMGTYPE" to be "arc"...
if DO_ARC_UPDATE:
    for c in ['red','blue']:
        for arc_index in arc_indices:
            fname = f"{c}{arc_index}.fits"
            fitspath = join(nightdir, fname)
            assert os.path.exists(fitspath)
            fits.setval(fitspath, 'IMGTYPE', value='arcs')
            fits.setval(fitspath, 'OBJECT', value='arcs')
            print(f'Updated IMGTYPE to "arcs" for {fitspath}')

# update "IMGTYPE" to be "science"...
if DO_SCIENCE_UPDATE:
    for c in ['red','blue']:
        for sci_index in science_indices:
            fname = f"{c}{sci_index}.fits"
            fitspath = join(nightdir, fname)
            if os.path.exists(fitspath):
                fits.setval(fitspath, 'IMGTYPE', value='object')
                print(f'Updated IMGTYPE to "object" for {fitspath}')

# update "IMGTYPE" to be "science"...
if DO_FLAT_UPDATE:
    for c in ['red','blue']:
        for flat_index in domeflat_indices:
            fname = f"{c}{flat_index}.fits"
            fitspath = join(nightdir, fname)
            assert os.path.exists(fitspath)
            fits.setval(fitspath, 'IMGTYPE', value='flat')
            print(f'Updated IMGTYPE to "flat" for {fitspath}')

# update "IMGTYPE" to be "science"...
if DO_JUNK_UPDATE:
    for c in ['red','blue']:
        for index in junk_indices:
            fname = f"{c}{index}.fits"
            fitspath = join(nightdir, fname)
            if os.path.exists(fitspath):
                fits.setval(fitspath, 'IMGTYPE', value='junk')
                fits.setval(fitspath, 'OBJECT', value='junk')
                print(f'Updated IMGTYPE to "junk" for {fitspath}')


# update "IMGTYPE" to be "science"...
if DO_OBJECTNAME_UPDATE:
    for c in ['red','blue']:
        for index in objectname_indices:
            fname = f"{c}{index}.fits"
            fitspath = join(nightdir, fname)
            assert os.path.exists(fitspath)
            fits.setval(fitspath, 'OBJECT', value=OBJECTNAME)
            print(f'Updated OBJECT to "{OBJECTNAME}" for {fitspath}')
