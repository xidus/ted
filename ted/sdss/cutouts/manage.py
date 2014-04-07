#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Mon 17 Mar 2014
#   Initial build.
#

"""
Manage cutouts

"""

import os

# ... : ted
# from ... import msg
# from ... import env

# . : ted.sdss.cutouts
from . import load_all_cutout_sequences


def remove_unwanted_data():

    import glob

    # Load all cutout sequences
    css, targets = load_all_cutout_sequences()
    dirs_ok = [os.path.split(cs.path('coord'))[1] for cs in css]

    # /path/to/data/cutouts/101x101
    base = css[0].path('root')
    print 'base:', base

    # Put the log file inthe parent directory
    # /path/to/data/cutouts
    opath = os.path.split(base)[0]
    ofname = os.path.join(opath, 'remove.log')

    with open(ofname, 'a') as fsock:
        fsock.write('\nNext remove log\n')

    # Get pattern for all the files in the directory `base`
    iglob = os.path.join(base, '*')
    print 'iglob', iglob

    # Assuming that only directories are present.
    # If not, these files should not be there either.
    dirs_present = [os.path.split(path)[1] for path in glob.glob(iglob)]

    for dir_present in dirs_present:
        if dir_present not in dirs_ok:
            path = os.path.join(base, dir_present)
            cmd = 'rm -rf {}'.format(path)
            print cmd
            os.system(cmd)
            with open(ofname, 'a') as fsock:
                fsock.write('{}\n'.format(path))


def remove_flags():

    import glob

    # Load all cutout sequences
    css, targets = load_all_cutout_sequences()

    # /path/to/data/cutouts/101x101
    base = css[0].path('root')
    print 'base:', base

    opath = os.path.split(base)[0]
    ofname = os.path.join(opath, 'remove_flagged.log')

    with open(ofname, 'a') as fsock:
        fsock.write('\nNext remove flagged log\n')

    # Get pattern for all the files in the directory `base`
    iglob = os.path.join(base, '*')
    print 'iglob', iglob

    # Create paths to flagged filenames
    fnames = [
        os.path.join(path, 'flagged')
        for path in glob.glob(iglob)
        if os.path.isdir(path)
    ]
    # os.remove raises an exception if files to be removed do not exist.
    fnames = [fname for fname in fnames if os.path.isfile(fname)]

    for fname in fnames:
        cmd = 'rm {}'.format(path)
        print cmd
        os.system(cmd)
        with open(ofname, 'a') as fsock:
            fsock.write('{}\n'.format(fname))


def remove_file(filename=None):

    import glob

    # Load all cutout sequences
    css, targets = load_all_cutout_sequences()

    # /path/to/data/cutouts/101x101
    base = css[0].path('root')
    print 'base:', base

    opath = os.path.split(base)[0]
    oname = os.path.splitext(filename)[0]
    ofname = os.path.join(opath, 'remove_{}.log'.format(oname))

    with open(ofname, 'a') as fsock:
        fsock.write('\nNext remove log\n')

    # Get pattern for all the files in the directory `base`
    iglob = os.path.join(base, '*')
    print 'iglob', iglob

    # Create full paths to filenames
    fnames = [
        os.path.join(path, filename)
        for path in glob.glob(iglob)
        if os.path.isdir(path)
    ]
    # os.remove raises an exception if files to be removed do not exist.
    fnames = [fname for fname in fnames if os.path.isfile(fname)]

    for fname in fnames:
        cmd = 'rm {}'.format(fname)
        print cmd
        print os.system(cmd)
        with open(ofname, 'a') as fsock:
            fsock.write('{}\n'.format(fname))

