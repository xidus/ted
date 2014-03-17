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

# . : cutouts
from . import load_all_cutout_sequences


def remove_unwanted_data():

    import glob

    # Load all cutout sequences
    css, targets = load_all_cutout_sequences()
    dirs_ok = [os.path.split(cs.path('coord'))[1] for cs in css]

    base = css[0].path('root')
    print 'base:', base

    opath = os.path.split(base)[0]
    ofname = os.path.join(opath, 'remove.log')

    with open(ofname, 'a') as fsock:
        fsock.write('\nNext remove log\n')

    # paths have structure: /path/to/data/cutouts/101x101/coord
    # Get the parent directory
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


