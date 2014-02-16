#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sat 4 Jan 2014
#   Based on notebook cutouts.ipynb.
#

"""
Create cutouts around a point in the sky.


"""

from __future__ import unicode_literals
# from __future__ import print_function

import os
import datetime as dt

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.io import fits
import scipy.misc as misc

# .. : ted
from .. import env
from .. import TEDError

# from ..time import mjd2date
from ..time import tai2date
from ..image import lingray

# . : sdss
# from . import iscoadded
from .das import frame_path
# from .aux import HTML_create_cutout_display
from .stripe82 import (
    ra_min, ra_max,
    # dec_min, dec_max,
    # stripe_width, stripe_height,
    stripe_extent_ra, stripe_extent_dec, w2h
)


class FITSHeaderError(TEDError): pass


class CutoutSequence(object):
    """
    A cutout sequence for a given coordinate and cutout size.

    To make a cutout, the following is needed:

    *   a coordinate covered by the data, e.g. from

        *   snlist (choose number, when loading)
        *   gxlist

    *   The size of the cutout

    From this, everything else can be built, processed and analysed.

    *   For each coordinate, prepare the following
        *   its root directory from radec coord
        *   subdirectories (should be a standard list)

    *   It should load the available fields from the .csv-file

    *   Selection I
        *   Choose only the frames that meet the following criteria
            *   <del>lower_bound < dRA < upper_bound</del>
                not needed, since this is already performed on the field list.
            *   coordinate covered by frame
        *   Create a plot for visual inspection (report and for stitching).

    *   Check the consistency of each file
        *   Uses astropy to open the file...
        *   Download if necessary (this should be a function in ted.sdss.das?)

    *

    """
    def __init__(self, radec=None, size=(101, 101), is_sn=False):

        # get_crd(ix=self.source_ix, source=self.source_key)
        self.radec = np.array(radec)[None, :]
        self.size = tuple(size)
        # May probably need to be changed to some more general labeling scheme.
        self.is_sn = is_sn

        # Get uniquely-formatted coordinate and size string
        self.crd_str = crd2str(self.radec.ravel())
        self.size_str = '{}x{}'.format(*self.size)

        # Define the pixel distance to the edges of the cutout image
        self.dy = np.floor(self.size[1] / 2)
        self.dx = np.floor(self.size[0] / 2)
        # print 'Distance to image edge pixels (dy, dx)', dy, dx

        # Prepare the coordinate input for the conversion from WCS to pixels
        # one where input is ordered (RA, Dec) and
        # one where input is ordered (Dec, RA)
        self.world = radec[None, :]
        self.world_r = radec[::-1][None, :]
        # print '(RA, Dec) =', world
        # print '(Dec, RA) =', world_r

        # Create empty dictionary for the background models (templates)
        self.bg_models = {}

        # Set modes
        self.initialised = False
        self.loaded = False

    # The time consumer
    # -----------------

    def initialise(self):
        """
        Initialise cutout sequence

        Steps
        -----
        * Create directories for the cutout sequence.
        * Load all fields and select those that cover the given coordinate.
        * Based on the list of covering fields, obtain filenames of valid
          FITS files that can be loaded. (Check file integrity for each
          file by loading it. If it is loadable, it is considered valid.)
        * Open each file in the list of valid FITS files. If there is
          sufficient space available around the pixel corresponding to the
          given coordinate, create a cutout and save it and .PNG and .fit.gz.
        * Gunzip the raw cutouts before registration (WCSREMAP requires this).
        * Perform WCS-registration using the first-recorded frame as template.
        * Create background models.

        The registered cutouts can now be loaded into a datacube for analysis.

        """

        # Define and build directory structure for this coordinate
        # Creates self.cpaths
        self.create_directories()

        # # Selection I
        # Get covering fields
        # self.fields = get_covering_fields(self.radec)
        self.get_covering_fields()

        # Selection II
        # creates self.files,
        # the list of filenames from which to extract raw cutouts.
        self.get_consistent_frames()

        # Selection III
        # Selects frames for which a cutout can be made.
        # Creates files on disk to be processed externally.
        self.create_raw_cutouts()

        # Cutout overview
        # Creates plots and statistics.
        # Where to svae these?
        # self.cutout_overview()

        # Process raw cutouts
        # -------------------

        # Gunzip raw cutouts
        self.gunzip()

        # USe WCSREMAP to astrometrically register cutouts to the oldest cutout
        self.wcsremap() # index=0

        # Set initialised mode
        self.initialised = True

    # The memory hurdle
    # -----------------

    def load(self, pad=0, bg_model='median', force=False):
        """
        Load remapped cutouts and perform preliminary steps for the analysis

        Parameters
        ----------
        pad : int, 0 < min(self.size)
            The registered frames have border effects from the resampling. Use
            pad to focus the analysis only on the central part of the frames.
            Only pixels that are farther than <pad> number of pixels from the
            border of the frames.
        bg_model : str, choice
            Which background model to subtract from the remapped frames to
            obtain the residual images?
        force : bool
            Whether or not to force running `self.initialise()`, if this step
            has not been done prior to running this method.

        Side effects
        ------------
        See respective methods.

        """

        if not self.initialised:
            print 'CutoutSequence not initialised ...'
            if force or 'y' == raw_input('Do you wish to initialise? (y/N)'):
                print 'Initialising CutoutSequence ...'
                self.initialise()
            else:
                return None

        # Load data cube
        # self.load_cutouts(which='residual', pad=pad)
        self.load_registered_cutouts(pad=pad)

        # Create background models
        # mean and median
        self.calculate_background_models()

        # Residuals
        self.calculate_residuals(bg_model=bg_model)

        # The data are now available for analysis.
        self.loaded = True

    # Work horses
    # -----------

    def is_transient(self, experiment='any'):
        """
        Finds out if there is a transient event in the sequence or not.

        Parameters
        ----------
        sigma : float
            Width of the LoG filter
        tau : float, 0 < tau < 1
            The relative intensity threshold above img(max - min) * tau
            of which at least one local maximum has to be in any of the
            images in order for the function to return True.

        Returns
        -------
        has_detected_transient : bool
            Whether or not a single local maximum in any of the cutouts is
            above the threshold defined by tau and at the given blob scale.

        """

        if experiment == 'any':
            """
            Experiment: Just see if there is anything.

            Is there any positive reponse at all and anywhere
            in the sequence, i.e. the cube of images?

            """
            return np.any(self.cube_cut.ravel())

            """
            Experiment: Positive only if consistens consecutive position of a signal.

            """

    # Service methods for self.initialise()
    # -------------------------------------

    def create_directories(self):
        """
        Generate a path structure that should be unique for each coordinate

        Side effects
        ------------
        self.cpaths : dict
            Dictionary of paths relevant for cutouts

        """

        # Initialise path dictionary
        cpaths = {}

        # The root directory for all the output generated for this coordinate
        # A directory for each cutout sequence (defined by its coordinate str)
        # will reside within this directory.
        cpaths['root'] = os.path.join(env.paths.get('cutouts'), self.size_str)

        # Define directory particular for this cutout sequence
        cpaths['coord'] = os.path.join(env.paths.get('cutouts'), self.crd_str)

        # Define subdirectories within the directory
        # particular for this cutout sequence
        subdirs = [
            ('fits', 'fit.gz'),
            ('png', 'png'),
            ('gunzip', 'fit'),
            ('wcsremap', 'wcsremap'),
            # ('residual', 'residual'),
            ('results', 'results'),
        ]
        for key, subdir in subdirs:
            cpaths[key] = os.path.join(cpaths.get('coord'), subdir)

        # Save paths
        self.cpaths = cpaths

        """Create all paths in self.cpaths"""
        for path in self.cpaths.itervalues():
            if not os.path.exists(path):
                os.makedirs(path)

    # @property
    def path(self, key):
        """Shortcut to retrieving a path from self.cpaths."""
        return self.cpaths.get(key)

    def get_covering_fields(self):
        """
        These fields form the basis input for the cutout algorithm which will
        determine the coverage further and chose frames that cover enough of
        the surrounding area that a cutout of the required size can be made.

        """

        self.fields = get_covering_fields(self.radec)

    def get_consistent_frames(self, redownload=True):
        """
        Check the consistency of each file.

        Should preferably have been done before the cutout algorithm, and
        non-repairable and non-downloadable files should be flagged so that
        they can be separated out above when testing the other criteria.

        I did not have time for this, however.

        Side effects
        ------------
        self.files : list
            list of filenames for the covering and existing frames

        """

        nframes = self.fields.shape[0]
        self.files = []
        self.notfiles = []
        self.coadded = []

        for i in range(nframes):
            # print '{: >3d}'.format(i)
            field = self.fields.iloc[i:i + 1]
            filename = frame_path(field)

            if field.run in (106, 206):
                self.coadded.append(i)
                # print 'CA', filename

            elif not os.path.isfile(filename):
                self.notfiles.append(i)
                # print u'!E', filename

            else:
                self.files.append(filename)

            # Some of the files that have been downloaded may be corrupted.
            # Try to load the file and catch any *IOError*s
            try:
                # print 'Checking file integrity ...'
                hdus = fits.open(filename)

            except IOError:
                # print 'IO', filename

                if not redownload:

                    self.files.pop()
                    continue

                else:

                    # Get download function
                    from .das import download_URI

                    # Close file before redownloading
                    if 'hdus' in dir():
                        hdus.close()

                    # print 'Trying to (re-)download ...'
                    URI = frame_path(field, local=False)
                    http_status, URI_, filename_ = download_URI((URI, filename))

                    if http_status is 200:
                        # print 'File was downloaded ...'

                        try:
                            # print 'Checking file integrity ...'
                            hdus = fits.open(filename)

                        except IOError:
                            # print 'Oh no! Something went wrong again. Skipping file for now ...'
                            self.files.pop()
                            continue

                        finally:
                            if 'hdus' in dir():
                                hdus.close()

                    else:
                        # print 'Response status code:', http_status
                        # print 'Skipping file for now ...'
                        self.files.pop()
                        continue

            finally:
                if 'hdus' in dir():
                    hdus.close()

        # self.fpCs = files
        # print 'Done ...'

    def create_raw_cutouts(self):
        """
        Cutout algorithm

        This block of code opens each frame in the list of files that
        made it through the above selections.

        When loaded, its WCS transformation is performed to get the
        (fractional) pixel coordinate of the corresponding input
        sky coordinate.

        If the distance from the pixel coordinate to any edge of the
        image is too small to make a cutout of the desired size around
        the coordinate, the image is skipped.

        Cutouts are made with the rest of the images.

        They are saved as {.png, .fit.gz} in the subfolders {png, fits},
        respectively.

        A small report on the cutout is also made.

        This report is saved as text, but an .HTML document is also
        created in the png directory, so that all the output cutouts
        and the stats can be viewd in one single document.

        """

        # Accounting
        # ----------

        # Which files where used?
        files_indices = []

        # Locations of output files which were written successfully to disk.
        # Also used to generate the .HTML report that will show each image.
        cutout_files = []

        # For the timeline
        cutout_dates = []

        # Filenames of files which the program failed to write to disk.
        failed_writes = []

        # Formats used to read and write formatted dates
        fmts = [
            '%Y-%m-%dT%H:%M:%S',
            '%Y-%m-%d',
            '%y/%m/%d',
        #     '%H:%M:%S',
        ]

        # For finding out what the largest possible cutout size is which
        # allows for a descent number of cutouts.
        # What is the largest possible pixel side length possible for
        # each frame?
        # = min( min(y, x), min(data.shape[0] - y, data.shape[1] - x) ) * 2 - 1
        pxmax = []

        # I want to map out what is going on.
        # Facts: The frames picked out above, clearly cover the candidate.
        #        The pixel coordinates should therefore also at least be within the region covered by the image data.

        # I need to find out how to get the pixels to be within the image data region.

        # Input:
        #   * the transformation in each frame.
        #   * the coordinate that the cutout has to be made around

        # Given the transformation,
        #   * what is the order of the returned pixel coordinates?
        #   * what is the expected order of the world coordinates
        #     that are fed into the conversion function?

        # Pixels
        # ------

        # All the returned pixel coordinates when the order of the input coordinates are consistently (RA, Dec),
        # and the returned pixel coordinates are consistently assumed to be (col_ix, row_ix)
        row_ixes_all = []
        col_ixes_all = []

        # All the returned pixel coordinates when the order of the input coordinates are consistently (RA, Dec),
        # and the returned pixel coordinates are consistently assumed to be (col_ix, row_ix)
        row_ixes_all_r = []
        col_ixes_all_r = []

        # Actual, used row and col indices of the cutouts that can be made
        row_ixes = []
        col_ixes = []

        # Clean up
        # --------

        # astropy.io.fits will not overwrite existing files,
        # so I delete them all before creating the cutouts.
        print 'Can not overwrite existing files, so remove old files ...'
        cmd_fits_clean = 'rm {}'.format(
            os.path.join(
                self.path('fits'), '*.fit.gz'
            )
        )
        print cmd_fits_clean
        print ''
        os.system(cmd_fits_clean)

        print 'Load each image and create cutout when possible (marked by \'*\' at the end)'
        # print 'In DS9 fpCs have East up, North right'
        for i, ifname_cutout in enumerate(self.fpCs):

            # Load the file
            hdulist = fits.open(ifname_cutout)

            # Get the primary HDU
            primary = hdulist[0]

            # Get the header of the primary HDU
            phd = primary.header

            # Get the image data
            image = primary.data

            # Load the coordinate transformation
            w = wcs.WCS(phd)
            # w.printwcs()

            # Are the image-data axis order and the naxisi order the same?
            # No, the image is loaded into a numpy array, and the order of the axes is row-first
            # whereas the order of the axes in the FITS header is column-first

            # i.e. print out:
            # print 'image.shape:      ({}, {})'.format(*image.shape)
            # print '(naxis1, naxis2): ({}, {})'.format(w.naxis1, w.naxis2)

            # Gives:
            # image.shape:      (1489, 2048)
            # (naxis1, naxis2): (2048, 1489)

            # Get the reference world coordinates
            crval1, crval2 = phd.get('CRVAL1'), phd.get('CRVAL2')
            # crpix1, crpix2 = phd.get('CRPIX1'), phd.get('CRPIX2')

            # Use reference world coordinates to retrieve the given reference pixel coordinates
            reference_pixels = w.wcs_world2pix([[crval1, crval2]], 1)

            # What if we switch them around?
            # reference_pixels_r = w.wcs_world2pix([[crval2, crval1]], 1)

            # print 'reference_pixels:  ', reference_pixels
            # print 'reference_pixels_r:', reference_pixels_r

            # reference_pixels:   [[ 1024.5   744.5]]
            # reference_pixels_r: [[ 358421.2155787 -303921.6830727]]

            # Then one gets values like the ones that I am more or less consistently getting.

            # SOLUTIONs:
            #
            # One way of getting the correct pixel coordinates could then be
            #  to just check for outragous pixel-coordinate values, and
            #  switch the world coordinate input order if necessary.

            # First sanity check: compare given and retrieved reference pixels.
            pix1, pix2 = reference_pixels[0, 0], reference_pixels[0, 1]

            # Second sanity check: Transform using both input orderings
            pixels = w.wcs_world2pix(self.world, 1)
            pixels_r = w.wcs_world2pix(self.world_r, 1)

            # print pixels
            # print pixels2
            # print phd.get('CTYPE1')

            # [[ 357706.51209122 -302862.75875357]]
            # [[ 1821.48191956   494.61177651]]
            # DEC--TAN

            # [[ 1947.93272917  1482.63382086]]
            # [[ 357461.49779276 -301603.01175018]]
            # RA---TAN

            # Are naxis1 and naxis2 as we would expect? (col_ix, row_ix)
            if pix1 > pix2:
                # i.e. we compare sth like (1024.5 > 744.5)

                # Compared to the image data array,
                # the returned pixel indices are
                # returned in the order
                #
                #   (column index, row index)
                #
                # I switch the assignment order, so that
                #     * y is the row index, and
                #     * x is the column index.

                if 'RA' in phd.get('CTYPE1'):

                    # The input order of the world coordinates has to be (RA, Dec)

                    # Columns first
                    x, y = pixels[0, 0], pixels[0, 1]

                    # Save all indices for debugging
                    col_ixes_all.append(x)
                    row_ixes_all.append(y)

                    col_ixes_all_r.append(pixels_r[0, 0])
                    row_ixes_all_r.append(pixels_r[0, 1])

                else:
                    # Still columns first
                    x, y = pixels_r[0, 0], pixels_r[0, 1]

                    # Save all indices for debugging
                    # Notice that the post fixes ('', '_r') are switched
                    col_ixes_all.append(x)
                    row_ixes_all.append(y)

                    row_ixes_all_r.append(pixels[0, 0])
                    col_ixes_all_r.append(pixels[0, 1])

            else:

                 # (Assuming that pix1 != pix2)

                # Compared to the image data array,
                # the returned pixel indices are
                # returned in the order
                #
                #   (row index, column index)
                #

                if 'RA' in phd.get('CTYPE1'):
                    y, x = pixels[0, 0], pixels[0, 1]

                    # Save all indices for debugging
                    col_ixes_all.append(x)
                    row_ixes_all.append(y)

                    col_ixes_all_r.append(pixels_r[0, 0])
                    row_ixes_all_r.append(pixels_r[0, 1])

                else:
                    x, y = pixels_r[0, 0], pixels_r[0, 1]

                    # Save all indices for debugging
                    # Notice that the post fixes ('', '_r') are switched
                    col_ixes_all.append(x)
                    row_ixes_all.append(y)

                    row_ixes_all_r.append(pixels[0, 0])
                    col_ixes_all_r.append(pixels[0, 1])

            # Now we have the x (data col) and y (data row) coordinates for the image data (in NumPy coordinates)
            # Calculate the maximum-possible-square-cutout side length
            max_pixel_side_length = min(
                min(y, x),
                min(
                    image.shape[0] - y,
                    image.shape[1] - x)
            ) * 2 - 1 # Has to be within the image data
            pxmax.append(max_pixel_side_length)

            # Print info for all frames
            print '{: >4d} {: >8s} (row, col) = ({: >4.0f}, {: >4.0f})'.format(i, phd.get('CTYPE1'), y, x),

            # Recap: y is the row index
            # Recap: x is the col index
            # print 'ROW {: >4.0f} {: >7.2f} {: >4.0f}'.format(dy, y, image.shape[0] - dy),
            # print 'COL {: >4.0f} {: >7.2f} {: >4.0f}'.format(dx, x, image.shape[1] - dx),

            is_within_ypad = (0 + self.dy <= y <= image.shape[0] - self.dy)
            is_within_xpad = (0 + self.dx <= x <= image.shape[1] - self.dx)

            if is_within_ypad and is_within_xpad:

                # Mark this to the screen
                print '*'

                # Append the pixel coordinates for plotting below
                row_ixes.append(y)
                col_ixes.append(x)

                # Save the index for the included input file
                files_indices.append(i)

                # Get cutout
                cutout_slice = image[
                    y - self.dy:y + self.dy + 1,
                    x - self.dx:x + self.dx + 1
                ]
                # cutout_slice = image

                # Get observation time
                tai = phd.get('TAI', None)

                try:
                    datetime_observed = tai2date(tai)

                except Exception:
                    # Use DATE-OBS instead
                    date_obs = phd.get('DATE-OBS', None)
                    time_obs = phd.get('TAIHMS', None)
                    datetime_obs = date_obs + 'T' + time_obs
                    try:
                        datetime_observed = dt.datetime.strptime(datetime_obs, fmts[0])

                    except Exception:
                        try:
                            datetime_observed = dt.datetime.strptime(date_obs, fmts[1])

                        except Exception:
                            try:
                                datetime_observed = dt.datetime.strptime(date_obs, fmts[2])

                            except Exception:
                                raise FITSHeaderError('The header does not contain any readable observation timestamp.')

                # Save the cutout date for the timeline
                cutout_dates.append(datetime_observed)

                # Format the date for the output filenames
                date_str = dt.datetime.strftime(datetime_observed, fmts[0])
                # date_str = dt.datetime.strftime(datetime_observed, '%s')

                # Save the cutout

                # As .FITS
                # --------
                # ...
                # Create new PrimaryHDU object
                # Attach the changed image to it
                # Add to the hdulist, change reference-CR+PX and so on...
                # Ref: http://prancer.physics.louisville.edu/astrowiki/index.php/Image_processing_with_Python_and_SciPy
                hdulist[0].data = cutout_slice
                hdulist[0].header.set('ORIGINAL', ifname_cutout[ifname_cutout.rfind('/') + 1:])

                # Could be nice to have the world-coordinate extent saved in the cutout, but it is not needed at the moment.
                if 0:
                    hdulist[0].header.set('RA max', )
                    hdulist[0].header.set('RA min', )
                    hdulist[0].header.set('Dec max', )
                    hdulist[0].header.set('Dec min', )

                # hdulist[0].header.set('', '')
                # Do we need all the testing above?
                # astropy.io.fits is loading the FITS image data in a NumPy array which is indexed row-first,
                # whereas astropy.wcs does what it does independently of how the data are being loaded by astropy.io.fits.
                # This means that we should always expect the return order of the pixel coordinates from astropy.wcs to be
                # (column index, row index), and that only the input order needs to be checked for.
                # Assuming that the input order is all that one needs to know in order to correctly
                # set the reference coordinates below.
                hdulist[0].header.set('CRPIX1', phd.get('CRPIX1') - (x - self.dx))
                hdulist[0].header.set('CRPIX2', phd.get('CRPIX2') - (y - self.dy))
                # hdulist[0].header.set('CRPIX1', x)
                # hdulist[0].header.set('CRPIX2', y)
                # if 'RA' in phd.get('CTYPE1'):
                #     # The order is (RA, Dec)
                #     hdulist[0].header.set('CRVAL1', world[0, 0])
                #     hdulist[0].header.set('CRVAL2', world[0, 1])
                # else:
                #     # The order is (Dec, RA)
                #     hdulist[0].header.set('CRVAL1', world_r[0, 0])
                #     hdulist[0].header.set('CRVAL2', world_r[0, 1])
                fname_fits = date_str + '.fit.gz'
                ofname_fits = os.path.join(self.path('fits'), fname_fits)
                hdulist.writeto(ofname_fits)

                # Was the file successfully saved to disk?
                if not os.path.isfile(ofname_fits):
                    failed_writes.append(i)
                else:
                    cutout_files.append(fname_fits)

                # As .PNG
                # -------

                # re-scale intensities (should only be done for the .png images)
                cutout_slice = lingray(cutout_slice)
                # cutout_slice = lingray(np.log10(cutout_slice))

                # Normalise
                cutout_slice = cutout_slice / cutout_slice.sum()

                # Save the image
                fname_png = date_str + '.png'
                ofname_png = os.path.join(self.path('png'), fname_png)
                misc.imsave(ofname_png, cutout_slice)

            # OK, the coordinate was not far enough from the edge of the image to make a cutout
            elif y < 0 or x < 0:

                # Some transformations gave this, print it out
                print '  NEGATIVE pixel coordinate !!!'

            else:

                # Do not mark anything
                print ''

            hdulist.close()

        # Summary
        # -------

        # This is what I need for cutout_overview()
        # self.files_indices = files_indices
        # self.cutout_files = cutout_files
        # self.cutout_dates = cutout_dates
        # self.failed_writes = failed_writes
        # self.pxmax = pxmax
        # self.row_ixes_all = row_ixes_all
        # self.row_ixes_all_r = row_ixes_all_r
        # self.row_ixes = row_ixes
        # self.col_ixes_all = col_ixes_all
        # self.col_ixes_all_r = col_ixes_all_r
        # self.col_ixes = col_ixes

    ### Image processing and analysis on cutouts for the given coordinate

    def gunzip(self):

        # Gunzip

        # Original command
        # From: https://superuser.com/questions/139419/how-do-i-gunzip-to-a-different-destination-directory
        #
        #  for f in *.gz; do
        #  STEM=$(basename $f .gz)
        #  gunzip -c "${f}" > fit/"${STEM}"
        #  done
        #
        # Alternative
        # <http://docs.python.org/2/library/gzip.html> (use Python for all the things!)
        #

        cmd_kw = dict(fits=self.path('fits'), gunzip=self.path('gunzip'))
        cmd_gunzip = '''\
cd {fits}
for f in *.gz; do
STEM=$(basename $f .gz)
gunzip -c "${{f}}" > {gunzip}/"${{STEM}}"
done\
'''.format(**cmd_kw)

        print(cmd_gunzip)
        os.system(cmd_gunzip)

    def wcsremap(self):

        # WCSREMAP
        import glob

        iglob = os.path.join(self.path('gunzip'), '*.fit')
        files = sorted(glob.glob(iglob))

        """FIXME"""
        # Need to be valid file that could have a cutout made out of it.
        template_frame = self.files[0]

        cmd_kw = dict(template=template_frame)
        cmd_wcsremap_fstr = '''\
wcsremap \
-template {template} \
-source {ifname} \
-outIm {ofname}\
'''

        # Do I need noise-output images?
        # This is a measure of the pixel-to-pixel
        # covariance for the resampled frames.

        print 'Processing files:'
        for i, ifname in enumerate(files):
            print '{: >3d}, '.format(i),
            if i and not (i + 1) % 10:
                print ''
            fname = os.path.basename(ifname)
            ofname = os.path.join(self.path('wcsremap'), fname)
            cmd_kw.update(ifname=ifname, ofname=ofname)
            cmd_wcsremap = cmd_wcsremap_fstr.format(**cmd_kw)

            # print(cmd_wcsremap)
            os.system(cmd_wcsremap)

    # Service methods for loading and initialising the data for the analysis
    # ----------------------------------------------------------------------

    # def create_residuals(self, bg_model='median'):
    #     """
    #     Creates residual images from the WCSREMAP'ed frames
    #     and saves them to disk in FITS format.

    #     Side effects
    #     ------------
    #     self.cube_residual : 3D-array
    #     Files on disk located at self.path('residual')

    #     """

    #     # Get self.cube_remap if not already in memory
    #     if not 'cube_remap' in dir(self):
    #         self._load_cube_remap()

    #     # Get self.template_median|mean
    #     self.calculate_background_models()

    #     # Get self.cube_residual
    #     self.calculate_residuals(bg_model=bg_model)

    #     # Get output filenames for residuals
    #     ifiles = self.get_filenames('wcsremap', '*.fit')
    #     ofiles = self.replace_in_string_list(ifiles, 'wcsremap', 'residual')

    #     # Get WCS header from template file.
    #     hdulist = fits.open(ifiles[0])
    #     _header = hdulist[0].header
    #     hdulist.close()

    #     # import copy

    #     for i, ofname in enumerate(ofiles):
    #         print '{: >3d}'.format(i),
    #         if i and not (i + 1) % 10: print ''
    #         # header = copy.deepcopy(_header)
    #         # header.set('', )
    #         fitskw = dict(data=self.cube_residual[:, :, i], header=_header)
    #         hdu = fits.PrimaryHDU(**fitskw)
    #         hdu.writeto(ofname)

    # def replace_in_string_list(self, list_=None, from_=None, to=None):
    #     """Do same replacement for each string in a list."""
    #     return map(lambda s: s.replace(from_, to), list_)

    def get_filenames(self, path=None, match=None):
        """Get all files within directory specified by path."""
        import glob
        iglob = os.path.join(self.path(path), match)
        return sorted(glob.glob(iglob))

    # def load_cutouts(self, which=None, pad=0):
    #     """
    #     Loads a data cube of the chosen kind.

    #     Parameters
    #     ----------
    #     which : str | list of strings
    #         What data to load? Choose among the following:
    #         * remap : the WCSREMAP'ed .fit-files
    #         * residual : the residual images .fit-files

    #     Side effects
    #     ------------
    #     self.cube_<which> : 3D-array
    #         Data cube containing cutout frames of the <which> kind.

    #     Future
    #     ------

    #     """

    #     if which == 'remap':
    #         self._load_cube_remap()

    #     else:
    #         # which := residual
    #         self._load_cube_residual(pad=pad)

    def load_registered_cutouts(self, pad=0):
        """
        Load all remapped images

        Assumptions
        -----------
        Using `files` assumes that no failed writings to disk occurred.
            (see Future below)

        Side effects
        ------------
        self.cube_remap : 3D-array
            Data cube containing the remapped cutout frames.

        Future
        ------
        Check integrity of the fits files found in the directory.
        Try: to load each image and save it into the data cube.
        Except: failed read-ins. Save the filename in list self.wcsremap_failed
        Else: Add the loaded image data to the cube. Save the filename into
            list for future reference.

        """

        files = self.get_filenames('wcsremap', '*.fit')
        cube_shape = self.size + (len(files),)
        self._cube_remap = np.zeros(cube_shape)
        for i, ifname in enumerate(files):
            print '{: >3d}, '.format(i),
            if i and not (i + 1) % 10: print ''
            hdulist = fits.open(ifname)
            self._cube_remap[:, :, i] = hdulist[0].data
            hdulist.close()

        self.cube_remap = self.get_view(self._cube_remap, pad=pad)

    # def _load_cube_residual(self, pad=0):

    #     files = self.get_filenames('residual', '*.fit')
    #     cube_shape = self.size + (len(files),)
    #     self._cube_residual = np.zeros(cube_shape)
    #     for i, ifname in enumerate(files):
    #         print '{: >3d}, '.format(i),
    #         if i and not (i + 1) % 10: print ''
    #         hdulist = fits.open(ifname)
    #         self._cube_residual[:, :, i] = hdulist[0].data
    #         hdulist.close()

    #     # Create view
    #     self.cube_residual = self.get_view(self._cube_residual, pad=pad)

    def get_view(self, data_cube, pad=0):
        """Creates a padded view into the data cube."""
        Ny, Nx = self.size
        return data_cube[pad:Ny - pad, pad:Nx - pad, :]

    def calculate_background_models(self):
        """
        Get the template images to subtract from all the images

        Side effects
        ------------
        self.bg_models[<bg_model>] : dict
            Updated with the possible background models.

        """

        # Median gives sharper edges close to the zero-padding
        self.bg_models['median'] = np.median(self.cube_remap, axis=2)

        # Smoother edges, but results are affected by sharp transients
        # self.template_mean = np.mean(self.cube_remap, axis=2)
        self.bg_models['mean'] = np.mean(self.cube_remap, axis=2)

    def calculate_residuals(self, bg_model=None):
        """
        Calculate the residual images

        Parameters
        ----------
        bg_model : str, choice
            Chosen background model.

        Side effects
        ------------
        self.template : 2D-array
            The chosen background model if available (None otherwise)
        self.cube_residuals : 3D-array

        """

        self.template = self.bg_models.get(bg_model)
        if self.template is None:
            print 'Model not in ({}).'.format(', '.join(self.bg_models.keys()))
            return None

        # Subtract the background template
        self.cube_residual = self.cube_remap - self.template[:, :, None]

    # Service methoids for the work horse method `self.is_transient`
    # --------------------------------------------------------------

    def calculate_LoG(self, **kwargs):
        """
        Calculate the Laplacian-of-Gaussian-filtered images

        kwargs should contain at least on of the following parameters.

        Parameters
        ----------
        radius : float, int
            The pixel radius of the blob. Probably an estimate.
        sigma : float
            The width of the LoG filter. If radius is given, the given value
            of sigma is ignored and instead calculated from radius.

        Side effects
        ------------
        self.cube_LoG : float, 3D-array
            The LoG-filtered residuals.

        """

        if 'radius' in kwargs:
            radius = kwargs.get('radius')
            # np.sqrt(2) * sigma = radius (width) of LoG filter
            # LoG(x, y) = (x ** 2 + y ** 2 - 2 * sigma ** 2) * G(x, y; sigma) / sigma ** 4
            # G(x, y; sigma) = exp( - (x ** 2 + y ** 2) / (2 * sigma ** 2) )
            sigma = radius / np.sqrt(2)

        elif 'sigma' in kwargs:
            # Now we have the width of the filter
            sigma = kwargs.get('sigma')

        else:
            print 'None of the required arguments were passed to the function.'
            return

        from scipy.ndimage import gaussian_laplace

        print 'width = {: >.2f} px =>'.format(radius)
        print 'sigma = {: >.2f} px'.format(sigma)

        # For each image obtain LoG(x, y; sigma)
        self.cube_LoG = np.zeros_like(self.cube_remap)
        for i in range(self.cube_remap.size[2]):
            self.cube_LoG[:, :, i] = gaussian_laplace(
                self.cube_residual[:, :, i], sigma=sigma
            )

    def compare_neighbours(self):
        """
        Eight-neighbour comparison

        For each residual image I_R, create four copies
        that are shifted, respectively,

            N   W
        * [-1, -1] (NW)
        * [-1, +1] (NE)
        * [+1, -1] (SW)
        * [+1, +1] (SE)
            S   E

        The local minima are the True entries in the matrix B

        B = (I_R < NW) & (I_R < NE) & (I_R < SW) & (I_R < SE)

        This binary field can then be overplotted on the original
        (remapped) cutout.

        Side effects
        ------------
        self.cube_minima_locs : float, 3D-array
            The data cube of bit fields for each 8-neighbour comparison.
            True : pixel location marks local minimum in LoG-filtered frames.
                This translates to a local maximum in the remapped frames.
            False : Not a local minimum of the LoG-filtered frames.

        """

        # Parameters for the padding below
        N_y, N_x = self.cube_remap.shape[:2]
        pad_y = 2 * N_y
        pad_x = 2 * N_x

        # Directions for 8-neighbour comparison
        directions = [
            (-1, 0), # N
            (+1, 0), # S
            (0, +1), # E
            (0, -1), # W

            (-1, -1), # NW
            (-1, +1), # NE
            (+1, -1), # SW
            (+1, +1)  # SE
        ]

        self.cube_minima_locs = np.zeros_like(self.cube_remap)
        for i in range(self.cube_remap.size[2]):
            I_res = self.cube_LoG[:, :, i]

            # this may be too expensive for larger cutout sizes,
            # but I did not have time to fiddle around with indices.
            I_res_3x3 = np.pad(I_res, (pad_y, pad_x), 'wrap')

            # Create matrix for each offset direction
            I_cmp = [
                I_res_3x3[
                    N_y + i: 2 * N_y + i,
                    N_x + j: 2 * N_x + j
                ] for (i, j) in directions]

            # Compare
            I_min = np.ones_like(I_res).astype(bool)
            for i in range(len(I_cmp)):
                I_min &= (I_res <= I_cmp[i])

            # Save into data cube of minima locations
            self.cube_minima_locs[:, :, i] = I_min

    def threshold_intensities(self, tau=.5):
        """
        Threshold each image in sequence.

        Side effects
        ------------
        self.cube_cut : float, 3D-array
            Data cube containing a bit field for each frame.
            True : Intensity of the local maximum in the remapped
                image is above threshold, i.e. the intensity is
                larger than the background. (Does it mean this?)
            False : Value was not larger than threshold.

        """

        self.cube_cut = np.zeros_like(self.cube_remap)
        for i in range(self.cube_remap.size[2]):

            # Best way to threshold?
            # Initial idea (KSP agreed after discussing it with him)
            t_cut = tau * (self.cube_remap[:].max() - self.cube_remap[:].min())
            # thresholding i SIP course was done on max(image),
            # i.e. negative intensities do not influence the positive range of
            # values used to determine the intensity value of the threshold.
            # This may introduce a bias for images with negative intensities.
            # t_cut = tau * (self.cube_remap[:].max())

            # Get thresholded binary field
            self.cube_cut[:, :, i] = (
                self.cube_minima_locs * self.cube_remap[:, :, i]
                ) > t_cut

    # def threshold_LoG_values(self, t_cut=-30):

    #     # USE: LoG intensities

    #     # t_cut = -13.
    #     t_cut = -30.

    #     # Get thresholded binary field
    #     I_min_thres = (cube_LoG[:, :, cut_ix] * I_min) < t_cut

    ## END class Cutout definition --------------------------------------------


def DEPRECATED_get_crd(ix=None, src='snlist'):
    """
    Get source coordinate.

    Parameters
    ----------

    Returns
    -------
    radec : np.array()

    """

    # Available sources
    src_lists = ('snlist', 'galaxies')

    if src in src_lists:

        if src == 'snlist':
            df_list = pd.read_csv(env.files[src], sep=';')
            radec = np.array([df_list.Ra[ix], df_list.Dec[ix]])

            # peak_MJD = df_list.Peak_MJD[ix]
            # peak_date = mjd2date(peak_MJD)
            # print peak_date

        elif src == 'galaxies':
            print('Not implemented yet...')
            raise SystemExit

        # print '(RA, Dec) = ({}, {})'.format(*radec)
        return radec


def crd2str(radec):
    """
    Returns a unique string-representation of coordinate

    """
    # print radec.flatten()
    # print '{:f}, {:f}'.format(*radec.flatten())
    return '{:0>10.5f}__{:0>10.5f}'.format(*radec).replace('.', '_')


def get_covering_fields(radec):
    """
    Get all fields that cover the coordinate.

    Parameters
    ----------
    radec : array
        length-2 array containing the ra and dec coordinate.

    Returns
    -------
    fields : pandas.DataFrame
        List of fields that cover the input coordinate.

    """
    radec_ = radec.ravel()
    # print radec_.shape
    # raise SystemExit

    # Raw list of fields
    fields = pd.read_csv(env.files.get('fields'), sep=',')

    # Get indices for fields that cover the chosen coordinate
    RA_geq_ramin = (fields.raMin <= radec_[0])
    RA_leq_ramax = (fields.raMax >= radec_[0])
    Dec_geq_decmin = (fields.decMin <= radec_[1])
    Dec_leq_decmax = (fields.decMax >= radec_[1])

    # The conjunction of these boolean arrays gives
    # the fields that satisfy the above need.
    ix = (RA_geq_ramin & RA_leq_ramax & Dec_geq_decmin & Dec_leq_decmax)

    return fields.loc[ix]

# ---


def plot_covering(radec, fields, opath):
    """
    Visualise selection

    """

    # Get indices for those fields that have RA \in [300; 360]
    subtract_ix = fields.raMax > ra_max

    # Make copies of these arrays for the visualisation
    df_vis = fields.copy()
    radec_vis = radec.copy()

    # Translate the coordinate
    df_vis.raMax[subtract_ix] -= 360
    df_vis.raMin[subtract_ix] -= 360

    # Only if there wore any subtrations to be made,
    # subtract 360 from the chosen coordinate so that
    # its position in the stripe can be visualised.
    if subtract_ix.sum():
        radec_vis[0] -= 360

    # Prepare plot input
    pkw = dict(mec='None', mfc='w', ls='None', marker='o', ms=5)
    rkw = dict(fill=True, fc=mpl.rcParams['axes.color_cycle'][0], ec='None', alpha=.3)

    # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15., 15. * 8 / w2h * 2))
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15., 15. * 8 / w2h * 3))

    # TOP plot
    # For each covering field
    for i in range(df_vis.shape[0]):
        cov_ramin, cov_ramax = df_vis.raMin.iloc[i], df_vis.raMax.iloc[i]
        cov_decmin, cov_decmax = df_vis.decMin.iloc[i], df_vis.decMax.iloc[i]
        dra, ddec = (cov_ramax - cov_ramin), (cov_decmax - cov_decmin)
        ax1.add_patch(mpl.patches.Rectangle((cov_ramin, cov_decmin), dra, ddec, **rkw))
    ax1.plot([radec_vis[0]], [radec_vis[1]], **pkw)
    ax1.plot(stripe_extent_ra, stripe_extent_dec, c=mpl.rcParams['axes.color_cycle'][1])
    ax1.set_xlim(ra_max + 1, ra_min - 1)
    ax1.set_xlabel(r'$\alpha\ /\ \mathrm{deg}$')
    ax1.set_ylabel(r'$\delta\ /\ \mathrm{deg}$')

    # MIDDLE plot
    for i, (rmn, rmx) in enumerate(zip(df_vis.raMin, df_vis.raMax)):
        ax2.plot([rmn, rmx], [i] * 2)
    ax2.axvline(x=radec_vis[0], c='w')
    ax2.set_xlim(*ax2.get_xlim()[::-1])
    ax2.set_xlabel(r'$\alpha\ /\ \mathrm{deg}$')
    ax2.set_ylabel(rmath('Frame index'))
    # ax2.ticklabel_format(scilimits=(-10, 10))
    ax2.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))

    # BOTTOM plot
    for i, (dmn, dmx) in enumerate(zip(fields.decMin, fields.decMax)):
        ax3.plot([i] * 2, [dmn, dmx])
    ax3.axhline(y=radec[1], c='w')
    ax3.set_xlabel(rmath('Frame index'))
    ax3.set_ylabel(r'$\delta\ /\ \mathrm{deg}$')

    fig.tight_layout()

    # Save it
    plt.savefig(os.path.join(opath, 'coverage.pdf'))
    plt.savefig(os.path.join(opath, 'coverage.png'), dpi=72)


def plot_possible_cutouts(pxmax, opath):
    """
    POSIBLE Cutouts
        # cutout_overview()
    """

    # How many cutouts of a given pixel side length can be made
    # out of the covering frames of the given coordinate?

    # NB: These are not all the covering images from the initial selection (Selection I)
    #     They are only those images that allowed for the given cutout size, e.g. 101x101.

    nbins = 100

    # for i, pxsz in enumerate(pxmax):
    #     print i, pxsz

    fig, ax = plt.subplots(1, 1, figsize=(15., 4))

    h = ax.hist(pxmax, bins=np.linspace(0, 1489, nbins + 1), zorder=10)

    Hcum = np.cumsum(h[0][::-1])[::-1]
    Hoff = h[1][:-1]
    ax.bar(left=Hoff, height=Hcum, width=Hoff[1] - Hoff[0], alpha=.5, fc=mpl.rcParams['axes.color_cycle'][1])

    ax.set_ylabel(rmath('No. of frames with cutout size possible'))
    ax.set_xlabel(rmath('Cutout side length in pixels'))

    # ax.set_xlim(.0, 1500.)
    ax.set_xlim(ax.get_xlim()[0], 1500.)

    fig.tight_layout()

    plt.savefig(os.path.join(opath, 'max_cutout_size.pdf'))
    plt.savefig(os.path.join(opath, 'max_cutout_size.png'), dpi=72)


def plot_pixel_indices(co):

    pkw = dict(ls='none', marker='o', ms=12, mec='None', alpha=.8)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15., 15. * 1489. / 2048. / 2))

    ax1.plot(co.col_ixes_all, co.row_ixes_all, **pkw)
    ax1.plot(co.col_ixes, co.row_ixes, **pkw)
    ax1.set_xlim(0, 2048 - 1)
    ax1.set_ylim(1489 - 1, 0)

    ax2.plot(co.col_ixes_all_r, co.row_ixes_all_r, **pkw)

    ax1.set_ylabel(rmath('y px coordinate'))
    for ax in (ax1, ax2):
        ax.set_xlabel(rmath('x px coordinate'))

    fig.tight_layout()

    plt.savefig(os.path.join(co.path('dim'), 'pixel_indices_joint.pdf'))
    plt.savefig(os.path.join(co.path('dim'), 'pixel_indices_joint.png'), dpi=72)


def plot_time_coverage(cutout_dates, opath):

    # Get matplotlib dates for the timeline
    cutout_mdates = mpl.dates.date2num(sorted(cutout_dates))
    cutout_mdates_diff = cutout_mdates[1:] - cutout_mdates[:-1]
    cmdiff_mean = cutout_mdates_diff.mean()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15., 4))

    ax1.axhline(y=1., c=mpl.rcParams['axes.color_cycle'][0])
    for mdate in cutout_mdates:
        ax1.axvline(x=mdate, ymin=.45, ymax=.55, lw=1., c=mpl.rcParams['axes.color_cycle'][1])
    ax1.set_xlim(cutout_mdates.min(), cutout_mdates.max())
    ax1.set_ylim(.0, 2.)
    ax1.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y'))
    ax1.xaxis.set_major_locator(mpl.dates.YearLocator())
    ax1.set_yticklabels([])

    ax2.hist(np.log10(cutout_mdates_diff))
    plt.draw()
    xticklabels = [10 ** float(t.get_text().replace(u'\u2212', u'-')) for t in ax2.get_xticklabels()]
    # print xticklabels
    ax2.set_xticklabels(xticklabels)

    fig.tight_layout()
    # Save it in the dimension subdirectory as above.
    plt.savefig(os.path.join(opath, 'stats.pdf'))
    plt.savefig(os.path.join(opath, 'stats.png'), dpi=72)


def write_cutout_summary(co):

    summary = ''
    summary += 'SUMMARY\n'
    summary += '\nSky coordinate: (RA, Dec) = ({:4f}, {:4f})'.format(*co.radec)
    summary += '\n\nInitial number of covering frames: {}'.format(co.fields.count().max())
    summary += '\nCo-added frames (left out for now): {}'.format(len(co.coadded))
    summary += '\nNon-existing file: {}'.format(len(co.notfiles))
    summary += '\nValid files to use for processing cutouts: {}'.format(len(co.fpCs))
    summary += '\n\nCutout dimension (rows, cols) = ({}, {})'.format(*co.size)
    summary += '\nPossible cutouts with this dimension: {}'.format(len(co.cutout_files))
    summary += '\nFailed writings to disk: {}'.format(len(co.failed_writes))
    # summary += '\n\nMean number of days between cutouts: {:.0f}d {:.3f}h'.format(
    #     cmdiff_mean // 1,
    #     (cmdiff_mean % 1) * 24
    # )

    ofname_summary = os.path.join(co.path('dim'), 'summary.txt')
    with open(ofname_summary, 'w+') as fsock:
        fsock.write(summary)

    print summary


def plot_background_models(self):
    # Set common vmin and vmax
    vmin = 0.
    vmax = np.max([self.template_median.max(), self.template_mean.max()])
    print vmax, np.log10(vmax)

    # Choose overall colour map
    cmap = mpl.cm.cubehelix

    # extent : None | (left, right, bottom, top)
    # default : assigns zero-based row, column indices
    #           to the `x`, `y` centers of the pixels
    extent = [0, cube_shape[0] - 1, 0, cube_shape[1] - 1]

    imkw = dict(cmap=cmap, interpolation='nearest')
    imkw.update(aspect='equal')  # , origin='lower')
    # extent = [0, cutout_size[0] - 1, 0, cutout_size[1] - 1]
    imkw.update(extent=extent)

    cbkw = dict(extend='neither', drawedges=False)
    cbkw.update(pad=.005, orientation='horizontal')

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(15., 5.))

    imkw.update(vmin=vmin, vmax=np.log10(vmax))
    ax1.imshow(np.log10(template_median), **imkw)
    ax2.imshow(np.log10(template_mean), **imkw)
    ax1.set_title(rmath('Log10(Median)'))
    ax2.set_title(rmath('Log10(Mean)'))

    imkw.update(vmin=vmin, vmax=vmax)
    ax3.imshow(template_median, **imkw)
    ax4.imshow(template_mean, **imkw)
    ax3.set_title(rmath('Median'))
    ax4.set_title(rmath('Mean'))

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_axis_off()
        cbar = fig.colorbar(mappable=im, ax=axes[i], **cbkw)

    if 0:
        """
        Maybe use this for the illustration that whos the uncropped frames.
        """
        padx = [pad - 1, extent[1] - pad + 1, extent[1] - pad + 1, pad - 1] + [pad - 1]
        pady = [extent[3] - pad + 1, extent[3] - pad + 1, pad - 1, pad - 1] + [extent[3] - pad + 1]
        print padx
        print pady
        for ax in (ax1, ax2, ax3, ax4):
            ax.plot(padx, pady)

    fig.tight_layout()
    plt.savefig(os.path.join(self.path('dim'), 'background_templates.pdf'))


def plot_residual_sample(self):
    # Show a sample of difference images

    # Offset to start from
    # Use it to see the change between SN present and not
    offset = 13

    # Use mask when plotting residuals
    # cube_res = np.ma.masked_where(cutout_mask_3D, cube_residual)
    cube_res = cube_residual

    vmin = cube_res[:, :, offset:offset + 4].min()
    vmax = cube_res[:, :, offset:offset + 4].max()
    print vmin, vmax

    # cmap = mpl.cm.bone
    cmap = mpl.cm.cubehelix

    imkw = dict(cmap=cmap, interpolation='nearest')
    imkw.update(aspect='equal')  # , origin='lower')
    imkw.update(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(1, 4, figsize=(15., 5.))

    cbkw = dict(extend='neither', drawedges=False)
    cbkw.update(pad=.005, orientation='horizontal')

    for i in range(4):
        index = offset + i
        im = axes[i].imshow(cube_res[:, :, index], **imkw)
        axes[i].set_title('cube_res[:, :, {}]'.format(index))
        axes[i].set_axis_off()
        cbar = fig.colorbar(mappable=im, ax=axes[i], **cbkw)

    fig.tight_layout()

    ofname_check_residual = os.path.join(path_LoG, 'check_residual.pdf')
    plt.savefig(ofname_check_residual)


def plot_LoG_samples(self):
    # Plot a couple of LoG-filtered images to check if filter size is reasonable (clear spot when SN present)

    # Also plot them as 3D surfaces to get a better feel for the toplogy of the the LoG-filtered difference images
    from mpl_toolkits.mplot3d import axes3d

    # Mask reference cube
    #cube_ref = np.ma.masked_where(cutout_mask_3D, cube_LoG)
    cube_ref = cube_LoG

    # Use log-scale?
    if 0:
        cube_ref = np.log10(cube_ref)
        title_fstr = 'Log10(cube_LoG[:, :, {}])'

    else:
        title_fstr = 'cube_LoG[:, :, {}]'

    # Offset to start from
    # Use it to see the change between SN present and not
    offset = 13

    vmin = cube_ref[:, :, offset:offset + 4].min()
    vmax = cube_ref[:, :, offset:offset + 4].max()
    print 'vmin = {:.2f}\nvmax = {:.2f}'.format(vmin, vmax)

    cmap = mpl.cm.cubehelix

    # Imshow kwargs
    imkw = dict(cmap=cmap, interpolation='nearest')
    imkw.update(aspect='equal')  # , origin='lower')
    imkw.update(vmin=vmin, vmax=vmax)

    # Axes3D kwargs
    axkw = dict(cmap=cmap, rstride=5, cstride=5)
    axkw.update(vmin=vmin, vmax=vmax)

    # Colour bar kwargs
    cbkw = dict(extend='neither', drawedges=False)
    cbkw.update(pad=.005, orientation='horizontal')

    # Coordinates
    y = np.arange(cube_shape[0])
    x = np.arange(cube_shape[1])
    X, Y = np.meshgrid(x, y)

    # Plot
    fig = plt.figure(figsize=(15., 8))
    # fig, axes = plt.subplots(1, 4, figsize=(15., 5.))

    for i in range(4):

        index = offset + i
        j = i + 4

        print i, j

        axi = fig.add_subplot(241 + i)
        axj = fig.add_subplot(241 + j, projection='3d')

        # In the top row
        im = axi.imshow(cube_ref[:, :, index], **imkw)
        axi.set_title(title_fstr.format(index))
        axi.set_axis_off()
        cbar = fig.colorbar(mappable=im, ax=axi, **cbkw)

        # In the bottom row
        Z = cube_ref[:, :, index]
        # axj.plot_wireframe(X, Y, Z, **axkw)
        surf = axj.plot_surface(X, Y, Z,  **axkw)
        # cbar = fig.colorbar(mappable=surf, ax=axes[i], **cbkw)
        axj.set_zlim(vmin, vmax)

    fig.tight_layout()

    ofname_check_LoG = os.path.join(path_LoG, 'check_LoG.pdf')
    plt.savefig(ofname_check_LoG)


def plot_detection_sample(self):

    # Plot I_min

    cmap = mpl.cm.binary

    # extent = [0, 80, 0, 80]
    imkw = dict(cmap=cmap, interpolation='nearest')
    imkw.update(aspect='equal')  # , origin='lower')
    imkw.update(extent=extent)

    pkw = dict(ls='None', marker='o', ms=12, mec='None', alpha=.5)
    pkw.update(c=mpl.rcParams.get('axes.color_cycle')[1])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11., 5.))

    ax1.imshow(I_min, **imkw)
    ax1.plot([40], [40], **pkw)
    ax1.set_title(rmath('Minima'))

    ax2.imshow(I_min_thres, **imkw)
    ax2.plot([40], [40], **pkw)
    ax2.set_title(rmath('Minima (intensity > threshold)'))

    # for ax in (ax1, ax2):
    #     ax.plot(padx, pady, c=mpl.rcParams.get('axes.color_cycle')[0])

    fig.tight_layout()

    ofname_check_minima = os.path.join(path_LoG, 'check_minima.pdf')
    plt.savefig(ofname_check_minima)


def plot_intensity_histogram(self):

    # Find a suitable threshold

    # Where are the True entries?
    y_ix, x_ix = I_min.nonzero()

    # The intensities of the LoG-minima
    cube_masked = cube_remap[y_ix, x_ix, cut_ix]

    # Indices of intensities that are larger than zero.
    cix = cube_masked > 0

    # Histogram of pixel values for the original remapped image
    # pkw.update(c=mpl.rcParams.get('axes.color_cycle')[1])
    fig, ax = plt.subplots(1, 1, figsize=(15., 4.))
    h = ax.hist(cube_masked[cix], bins=50)
    # ax.set_yscale('log')
    ax.set_xlabel(rmath('Intensity'))
    ax.set_ylabel(rmath('Count'))
    ax.set_title(rmath('Pixel intensities from remapped frame'))
    fig.tight_layout()

    # print h[0]

    ofname_check_intensities = os.path.join(path_LoG, 'check_intensities_remap.pdf')
    plt.savefig(ofname_check_intensities)


def plot_histogram_LoG(self):

    # How are the LoG pixel values distributed?

    LoG = cube_LoG[:, :, cut_ix]
    # cix = LoG < np.inf
    cix = LoG < -1

    # Histogram of pixel values for the original remapped image
    # pkw.update(c=mpl.rcParams.get('axes.color_cycle')[1])
    fig, ax = plt.subplots(1, 1, figsize=(15., 4.))
    # h = ax.hist(array1D_LoG, bins=50)
    # h = ax.hist(cube_LoG[:, :, cut_ix].flat, bins=50)
    h = ax.hist(LoG[cix], bins=50)
    # ax.set_yscale('log')
    ax.set_title(rmath('Pixel intensities of the LoG-filtered residual'))
    fig.tight_layout()

    # print h[0]

    ofname_check_intensities = os.path.join(path_LoG, 'check_intensities_LoG.pdf')
    plt.savefig(ofname_check_intensities)


def plot_histogram_ranked_LoG(self):

    # How are the relative changes in pixel intensity of the ranked and sorted LoG values distributed?

    LoGs = LoG.copy().flatten()
    # print LoGs.shape
    LoGs.sort()
    # print LoGs[0], LoGs[-1]

    # Relative differences
    N_rel = 10
    LoGr = LoGs[:N_rel] / LoGs[1:N_rel + 1] - 1.

    print '; '.join(['{:.5%}'.format(R) for R in LoGr])


    if 0:
        # Histogram of pixel values for the original remapped image
        # pkw.update(c=mpl.rcParams.get('axes.color_cycle')[1])
        fig, ax = plt.subplots(1, 1, figsize=(15., 4.))
        # h = ax.hist(array1D_LoG, bins=50)
        # h = ax.hist(cube_LoG[:, :, cut_ix].flat, bins=50)
        h = ax.hist(LoG[cix], bins=50)
        # ax.set_yscale('log')
        ax.set_title(rmath('Pixel intensities of the LoG-filtered residual'))
        fig.tight_layout()

        # print h[0]

        ofname_check_intensities = os.path.join(path_LoG, 'check_intensities_LoG.pdf')
        plt.savefig(ofname_check_intensities)


def get_cutout_sequences():
    """
    Returns
    -------
    cutout_sequences : np.array, 1D
        List of CutoutSequence objects.
    targets : np.array, 1D
        List of targets (feature labels) for each cutout.
        For now, just 1 feature : bool is_sn

    """

    tlist = pd.read_csv(env.files.get('tlist'))

    # print tlist.info()
    # print tlist.head()

    params = dict(size=(101, 101))

    cutout_sequences = []
    targets = []
    for i in range(tlist.shape[0]):
        # print '{: >3d}'.format(i)

        # Create cutout object
        params.update(
            radec=np.array([tlist.Ra[i], tlist.Dec[i]]),
            is_sn=tlist.is_sn[i]
        )
        cutout_sequences.append(CutoutSequence(**params))
        targets.append(tlist.is_sn[i])

    return np.array(cutout_sequences), np.array(targets)


def create_cutout_data():
    """
    Initialise directories and create data for analysis.

    """
    css, targets = get_cutout_sequences()
    for i, cs in enumerate(css):
        # print '{: >3d}'.format(i)
        cs.initialise()


def do_grid_search(N_folds=5):

    css, targets = get_cutout_sequences()

    # Create training and test data sets
    N_css = css.size
    N_items = N_css / N_folds

    # 1-fold CV to begin with
    css_test = css[:N_items]
    targets_test = targets[:N_items]

    css_train = css[N_items:]
    targets_train = targets[N_items:]

    # Define hyper-parameter space
    N_sigmas = 3
    N_taus = 3
    radii = np.linspace(5, 7, N_sigmas)
    # sigmas = np.linspace(, N_sigmas)
    sigmas = radii / np.sqrt(2)
    taus = np.linspace(.50, .70, N_taus)

    # Cube of predictions
    # N_sigmas x N_taus
    N_train = N_css - N_items
    cube_predict = np.zeros((N_sigmas, N_taus, N_train)).astype(bool)

    """
    For each cutout sequence in the training set,
    do grid search and report prediction.

    One slice in the third dimension of cube_predict
    contains predictions for one cutout sequence.

    """

    # Load and setup the data for all the cutouts
    for k in range(N_train):

        # Select the training data point (the cutout sequence)
        cs = css_train[k]

        # Load registered frames and get residuals using <bg_model>
        cs.load(pad=10, bg_model='median')

        # Grid search
        # -----------

        # Sigmas
        for i in range(N_sigmas):
            # Blobs
            # No need to perform LoG every time I change between thresholds
            cs.calculate_LoG(sigma=sigmas[i])
            cs.compare_neighbours()

            # Threshold values
            for j in range(N_taus):
                cs.threshold_intensities(tau=taus[j])

                # Save the resulting prediction
                cube_predict[i, j, k] = cs.is_transient()

    # Training accuracy
    # N_sigmas x N_taus
    # ^ : XOR
    accuracies = (cube_predict ^ targets_train[None, None, :]).sum(axis=2) / N_train

    # Get indices of the maximum
    ii, jj = np.unravel_index(accuracies.argmax(), (N_sigmas, N_taus))

    max_entries = (accuracies == accuracies[ii, jj])
    N_max_entries = max_entries.sum()
    if N_max_entries > 1:
        print 'More than one best hyper-parameter combination.'
        # print max_entries.nonzero()

    """
    For the test data, there will be only a vector of predictions,
    i.e. one for each cutout sequence in the test set.

    The accuracy of this set of predictions is the main
    result of the 1-fold grid search.

    """
    sigma_ii = sigmas[ii]
    tau_jj = taus[jj]

    vec_predict = np.zeros(N_items).astype(bool)
    for i in range(N_items):
        cs = css_test[i]
        cs.load(pad=10, bg_model='median')
        cs.calculate_LoG(sigma=sigma_ii)
        cs.compare_neighbours()
        cs.threshold_intensities(tau=tau_jj)
        vec_predict[i] = cs.is_transient()

    # Accuracy is a scalar, since we sum over the only dimension
    accuracy = (vec_predict ^ targets_test).sum() / N_items

    # Save results
    ofname_training = os.path.join(cs.path('root'), 'accuracies.csv')
    with open(ofname_training, 'w+') as fsock:
        for row in accuracies:
            fsock.write(','.join(row) + '\n')

    ofname_test = os.path.join(cs.path('root'), 'accuracy_test')
    with open(ofname_test, 'w+') as fsock:
        fsock.write('{}\n'.format(accuracy))


def main():
    """
    Control script for the functions beneath.

    Managagement script to make only the needed steps.
    Keep a logfile of steps completed?
    Yes.

    Basically a dictionary with key = crd2path

    And then do a complete cleanup, when needed?
    Yes

    In this way, I can also clean up entire root
    directories by running through the registration file.

    Should there be a source in this file? E.g. from gxlist or snlist?
    Yes.

    """

    # from mplconf import mplrc
    # from mplconf import rmath

    # mplrc('publish_digital')


    # plot_covering(co.radec, co.fields_covering, co.path('root'))
    # plot_possible_cutouts(pxmax=co.pxmax, opath=co.path('root'))
    # plot_time_coverage(cutout_dates=co.cutout_dates, opath=co.path('dim'))
    # plot_pixel_indices(co=co)
    # plot_time_coverage(cutout_dates=co.cutout_dates, opath=co.path('dim'))
    # write_cutout_summary(co=co)
    # HTML_create_cutout_display(co.cutout_files, co.path('dim'))


