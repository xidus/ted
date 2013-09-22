#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sun 1 Sep 2013
#   Initial build.
#

import os

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy import wcs

from .. import paths
from ..converters import *
from ..sdss import stripe82


# http://www.diveintopython.net/unit_testing/stage_1.html
class PipeError(Exception): pass
class OutOfRangeError(PipeError): pass


# Define the needed paths
_path_data_small = paths['data']['small']['path']
_path_data_big = paths['data']['big']['path']


def skyml_get_cutout(ra, dec, cutout_size=(101, 101), cutout_name=None, fn_fields=None):
    """Generate a sequence of image files of given width and height
    that have the given coordinate pair at their centre.

    ALso, get only frames that were captured within the specified
    time limits. If not given, produce the full possible sequence.

    Parameters
    ----------

    ra : degrees, floating-point

    dec : degrees, floating-point

    cutout_size : tuple of two positive integers or a single integer.
        If given as a single integer int_, it will be converted into
        a tuple (int_, int_).

    Steps
    -----

    - Check that (ra, dec) lie within Stripe 82.

    - Validate cutout_size and make sure it is a tuple.

    - Create a directory within

        $DATA_BIG/cutouts

      named by the coordinate pair in the format

        $DATA_BIG/fpC/cutouts/<RAdeg><Decdeg> (?) (to within two decimal places)

      The filenames for the cutouts should be
      the MJD of the first-row readout.

    - Load into memory the .csv file with the unique
      list of all the frames in the data set.

    - Query the list to get only the entries that cover the given coordinate.

    - Build the local path to each file. (should use a sub function)

    - Load into memory each frame and do the following:

        * Map the coordinate into the pixel location of the frame.
        * Check to see if the cutout window is completely covered by the frame.
        * If so, extract the cutout (slicing) and save it
          as a (FITS or PNG?) in the dedicated folder.

    - Produce a movie from the cutouts.


    BEFORE 2013-09-02
    -----------------

    (RA, Dec) should be in degrees.

    Check that they lie within the range of Stripe 82.

    Create a directory within

        $DATA_BIG/cutouts

    named by the coordinate pair in the format

        $DATA_BIG/cutouts/<RAdeg><Decdeg>

    Here, save the response from the SQL query.

    Create directories `fpC` and `output` for the frames (fpC)
    and final cutouts.

    The filenames for the cutouts should be the MJD of the
    first-row readout.

    Download the frames using the custom output directory above.

    Then, go through each frame and check if a cutout can be made,
    i.e. if the frame actually covers the whole pixel area.

    Then do the cutout (slicing) and save the image in the output
    directory.

    ---

    OLD THOUGHTS:

    Convert pixel-distance from the point to the cutout edges
    into World Coordinate System for the FITS files.

        * Assume that this is the same for each frame?
        * What about registration?

    Construct SQL query that returns fields that cover the entire
    cutout area as given by the

    Calibrate the central pixel:

    Determine from FITS-file what pixel coordinate that is closest
    to the given (ra, dec). This pixel point will be the centre of
    the image.

    """

    if not (stripe82.ra_min <= ra <= stripe82.ra_max):
        raise OutOfRangeError, 'Your Right Ascension is outside of Stripe 82'

    if not (stripe82.dec_min <= dec <= stripe82.dec_max):
        raise OutOfRangeError, 'Your Declination is outside of Stripe 82'

    if fn_fields is None:
        # Choose the unique list
        fn_fields = os.path.join(_path_data_small, 'fpC_unique_URIs.txt')

    if cutout_name is None:
        # Use the coordinates
        cutout_name = ('{:0>+7.3f}' * 2).format(ra, dec)

    # Define the directory path for the cutouts for the given coordinate
    opath_cutouts = os.path.join(_path_data_big, 'cutouts', cutout_name)
    print opath_cutouts

    # Load the list
    # How many in total?
    df = pd.read_csv(fn_fields, sep=';')
    print 'Number fo fields in total:', len(df)

    # Get frames that cover the coordinate
    # How many of the total number cover the coordinate?
    ra_ix = (df['raMin'].values <= ra) & (df['raMax'].values >= ra)
    dec_ix = (df['decMin'].values <= dec) & (df['decMax'].values >= dec)
    cover_ix = ra_ix & dec_ix

    df_cover = df[cover_ix]

    # Load each frame and determine whether it allows for a cutout of the desired dimensions.
    # How many of the covering frames are good enough (valid)?
    # Save the number of pixels that they each miss in order to be included.
    # ...OR save the number of pixels that would allow them to BE included in the cutout output.

    # Make a cutout of each valid frame, and save it as gray-scale PNG in the cutoutfolder.

    fields = CAS_query(sql_fields)
    if 'error' in fields.lower():
        sys.exit('ERROR: {}: CAS says something is wrong ...'.format(
            sys._getframe().f_code.co_name)
        )
    with open(os.path.join(opath, '{}.csv'.format(SDSS_id)), 'w+') as fsock:
        fsock.write(fields)


def main():
    ra = np.random.random() * stripe82.stripe_width - stripe82.stripe_width / 2.
    dec = np.random.random() * stripe82.stripe_height - stripe82.stripe_height / 2.
    cutout_name = 'TestSN'
    cutout_size = (101, 101)
    skyml_get_cutout(ra=ra, dec=dec, cutout_size=cutout_size, cutout_name=cutout_name)

if __name__ == '__main__':
    main()

