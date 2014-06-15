#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sat 4 Jan 2014
#   Based on notebook cutouts.ipynb.
#

"""
Create cutouts around a point in the sky.

"""

# from __future__ import unicode_literals # Some numpy functionality crashes,
# from __future__ import print_function

import os
import datetime as dt

import numpy as np
import pandas as pd
from astropy import wcs
from astropy.io import fits
import scipy.misc as misc
import yaml

# ... : ted
from ... import msg
from ... import env
from ... import TEDError

# from ...time import mjd2date
from ...time import tai2date
from ...image import lingray

# . : sdss
# from . import iscoadded
from ..das import frame_path
# from .aux import HTML_create_cutout_display

from .plotting import plot_covering
from .plotting import plot_possible_cutouts
# from .plotting import plot_pixel_indices
from .plotting import plot_time_coverage

# MAinly called from illustration control
from .plotting import plot_background_models
from .plotting import plot_residual_samples
from .plotting import plot_LoG_samples
from .plotting import plot_binary_fields

# from .plotting import plot_histogram_intensity

# from .plotting import plot_histogram_LoG

# from .plotting import plot_histogram_ranked_LoG



class FITSHeaderError(TEDError): pass
class CutoutSequenceError(TEDError): pass


class Log(object):

    def __init__(self, fname, sep=','):
        self.fname = fname
        self.sep = sep

    def _log(self, s):
        with open(self.fname, 'a') as fsock:
            fsock.write('{}\n'.format(s))

    def log(self, s):
        now = dt.datetime.now().isoformat()
        self._log('{now}{sep}{s}'.format(
                now=now, sep=self.sep, s=s))

    def log_empty(self):
        self._log('')


class CutoutSequence(object):
    """
    A cutout sequence for a given coordinate and cutout size.

    To make a cutout, the following is needed:

    *   a coordinate covered by the data, e.g. from

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

    """
    Since there are more than 17000 available coordinates from the galaxy list
    to use, I am setting this pretty high to begin with.
    This threshold is something that should be played with.
    Setting it too high may eliminate events from the event list, since they
    will be replaced by galaxy coordinates, thus reducing the ratio of events
    to the number of non-events.
    """
    MIN_NUMBER_OF_CUTOUTS = 50

    def __init__(self, radec=None, size=(101, 101), is_sn=False):

        # get_crd(ix=self.source_ix, source=self.source_key)
        self.radec = np.array(radec)[None, :]
        self.size = tuple(size)
        # May probably need to be changed to some more general labeling scheme.
        self.is_sn = bool(is_sn)

        # Get uniquely-formatted coordinate and size string
        self.crd_str = crd2str(self.radec.ravel())
        self.size_str = '{}x{}'.format(*self.size)

        # Define the pixel distance to the edges of the cutout image
        self.dy = np.floor(self.size[1] / 2)
        self.dx = np.floor(self.size[0] / 2)
        # print 'Distance to image edge pixels (dy, dx)', dy, dx

        # Prepare the coordinate input for the conversion from WCS to pixels
        # where input is ordered (RA, Dec) and
        # where input is ordered (Dec, RA)
        self.world = radec[None, :]
        self.world_r = radec[::-1][None, :]

        # Define directory paths for this instance
        self.define_directory_paths()

        # Create empty dictionary for the background models (templates)
        self.bg_models = {}

        # Create a dictionary for saved steps
        self._file = {}

    # Helpers
    # -------

    def logged(fun):
        """A decorator for logging"""
        def wrapped(self, *args, **kwargs):
            if not self.step_finished(fun.func_name):
                fun(self, *args, **kwargs)
                self.log(step=fun.func_name, status='finished')
        return wrapped

    def log_init(self):
        """Create log file to save the steps performed so far."""

        # Default settings for the log file
        self._fn_log = os.path.join(self.path('coord'), 'log.yaml')
        self._log_kw = dict(default_flow_style=False, explicit_start=True)
        self._log_kw.update(indent=2)

        # Create the file, if it does not exist
        if os.path.isfile(self._fn_log):
            self.log_load()

        else:
            # self._log = cs.OrderedDict()
            self._log = {}
            self.log_save()

    def log_load(self):
        """Load the logged steps"""
        with open(self._fn_log, 'r') as fsock:
            # self._log = cs.OrderedDict(yaml.load(fsock.read()))
            self._log = yaml.load(fsock.read())

    def log_save(self):
        """Overwrite last saved log file."""
        with open(self._fn_log, 'w+') as fsock:
            yaml.dump(self._log, stream=fsock, **self._log_kw)

    def step_finished(self, step=None):
        """Return status of the step"""
        return self._log.get(step) == 'finished'

    def log(self, step=None, status=None):
        """Log method to be run at the end of a method."""
        self._log.update({step: status})
        self.log_save()

    # It this instance for some reason is not usable, it can be flagged.
    @logged
    def flag(self, reason='unspecified'):
        """Write flag reason to file in the coordinate directory"""
        ofname = os.path.join(self.path('coord'), 'flagged')
        with open(ofname, 'w+') as fsock:
            fsock.write('{}\n'.format(reason))

    @property
    def flagged(self):
        return os.path.isfile(os.path.join(self.path('coord'), 'flagged'))

    # The time consumer
    # -----------------

    def initialise(self):
        """
        Initialise cutout sequence

        Steps
        -----
        * Load all fields and select those that cover the given coordinate.
        * Based on the list of covering fields, obtain filenames of valid
          FITS files that can be loaded. (Check file integrity for each
          file by loading it. If it is loadable, it is considered valid.)
        * Open each file in the list of valid FITS files. If there is
          sufficient space available around the pixel corresponding to the
          given coordinate, create a cutout and save it and .PNG and .fit.gz.
        * Gunzip the raw cutouts before registration (WCSREMAP requires this).
        * Perform WCS-registration using the first-recorded frame as template.

        The registered cutouts can now be loaded into a datacube for analysis.

        """

        # Create directory structure for this coordinate
        msg('Creating directories')
        self.create_directories()

        # Save first entry to log files
        msg('Initialising log file')
        self.log_init()

        # Is this instance flagged?
        if self.flagged:
            # Then stop here
            err = 'CutoutSequence is flagged!'
            raise CutoutSequenceError(err)

        # Create raw cutout sequence
        # --------------------------

        # # Selection I
        # Get covering fields
        # self.fields = get_covering_fields(self.radec)
        msg('Getting covering fields')
        # if not self.step_finished(step='get_covering_fields'):
        self.get_covering_fields()

        # Selection II
        # creates self.files,
        # the list of filenames from which to extract raw cutouts.
        # self.get_consistent_frames()
        # Since it sucks big time that I can not run TSOCKS properly on
        # imagerserver[1-3], I am leaving out the steps that try to download
        # missing or corrupted files.
        msg('Getting consistent frames')
        self.get_consistent_frames()

        # Selection III
        # Selects frames for which a cutout can be made.
        # Creates files on disk to be processed externally.
        msg('Creating raw cutouts')
        self.create_raw_cutouts()

        # Cutout overview
        # Creates plots and statistics.
        # Where to svae these?
        # self.cutout_overview()

        # Process raw cutouts
        # -------------------

        # Gunzip raw cutouts
        msg('GUNZIP')
        self.gunzip()

        # Use WCSREMAP to astrometrically register cutouts to the oldest cutout
        msg('WCSREMAP')
        self.wcsremap() # index=0

    # The memory hurdle
    # -----------------

    def load(self, clip=None, bg_model=None, quality=None):
        """
        Load remapped cutouts and perform preliminary steps for the analysis

        Parameters
        ----------
        clip : int, 0 < min(self.size)
            The registered frames have border effects from the resampling. Use
            clip to focus the analysis only on the central part of the frames.
            Only pixels that are farther than <clip> number of pixels from the
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

        # Set the clip size
        if clip is not None: self.set_clip(clip)

        # Set the background model once for easy reference
        if bg_model is not None: self.set_bg_model(bg_model)

        # Set qualities to load into memory
        if quality is not None: self.set_load_quality(quality)

        # Set the quality to use to the same as the load quality.
        self.set_quality(quality=self.load_quality)

        # Load data into memory
        self.load_registered_cutouts()

        return self

    def set_cs_parameters(self, clip, bg_model, quality):
        """
        Set all CutoutSequence parameters in one go.

        Notes
        -----
        Used by CVHelper() in crossvalidation module.

        """

        # Set the clip size
        self.set_clip(clip)

        # Set the background model once for easy reference
        self.set_bg_model(bg_model)

        # Set qualities to load into memory
        self.set_load_quality(quality)

        # Set the quality to use to the same as the load quality.
        self.set_quality(quality=self.load_quality)


    def set_clip(self, clip=0):
        self.clip = clip

    def set_bg_model(self, bg_model=None):
        self.bg_model = bg_model

    @staticmethod
    def _check_quality(quality):
        if not isinstance(quality, list):
            print '`quality` is not a list instance ...'
            raise SystemExit

    def set_load_quality(self, quality):
        """Set the qualities to load into memory"""
        self._check_quality(quality)
        self.load_quality = quality

    def set_quality(self, quality):
        """Set the quality to use for the analysis"""
        self._check_quality(quality)
        print 'Setting quality to use:', quality
        self.quality = quality

    @property
    def template(self):
        return self.bg_models.get(self.bg_model)

    def cleanup(self):
        """Clears loaded and generated data variables from memory."""
        for key in dir(self):
            if '__' not in key and 'cube_' in key:
                print 'Clearing', key, '...'
                delattr(self, key)

    # Service methods for self.initialise()
    # -------------------------------------

    def define_directory_paths(self):
        """
        Define a path structure that should be unique for each coordinate

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
        cpaths['coord'] = os.path.join(cpaths.get('root'), self.crd_str)

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

    def create_directories(self):
        """Create all paths in self.cpaths"""
        for path in self.cpaths.values():
            if not os.path.exists(path):
                os.makedirs(path)

    # @property
    def path(self, key):
        """Shortcut to retrieving a path from self.cpaths."""
        return self.cpaths.get(key)

    def file(self, key):
        """Shortcut to retrieving a file path from self._file."""
        return self._file.get(key)

    def add_file(self, key=None, value=None):
        """Shortcut to adding a file path to self._file."""
        if key is not None:
            self._file.update({key:value})

    def get_covering_fields(self):
        """
        These fields form the basis input for the cutout algorithm which will
        determine the coverage further and chose frames that cover enough of
        the surrounding area that a cutout of the required size can be made.

        """

        # Define filename for previously/to-be saved output
        fname = os.path.join(self.path('coord'), 'covering_fields.csv')
        self.add_file(key='fields', value=fname)

        # Load a file or calculate again?
        if os.path.isfile(fname):
            self.fields = load_fields()

            # Plot coverage overview for given cutout
            # Un-commenting reason: see last lines of self.create_raw_cutouts()
            plot_covering(self.radec, self.fields, self.path('coord'))

        else:
            self.fields = get_covering_fields(self.radec)
            self.fields.to_csv(fname, index=False, headers=True)

        # Display output to screen
        msg('Found {:d} fields that cover coordinate {}'.format(
            self.fields.shape[0], self.crd_str)
        )

    def get_consistent_frames(self):
        """
        Check the consistency of each file.

        Should preferably have been done before the cutout algorithm, and
        non-repairable and non-downloadable files should be flagged so that
        they can be separated out above when testing the other criteria.

        I did not have time for this, however.

        Raises
        ------
        CutoutSequenceError : if no cutouts could be made at all.

        Side effects
        ------------
        self.frames : str, 1D-array
            Filenames for the existing and non-coadded frames
        self.frames_coadded : str, 1D-array
            Filenames for the existing and coadded frames
        self.frames_invalid_FITS : str, 1D-array
            Filenames for the existing frames woth invalid content
        self.frames_nonexsitent : str, 1D-array
            Filenames for the non-existing frames

        Dependencies
        ------------
        self.get_covering_fields()

        """

        keys = [
            'frames',
            'frames_coadded',
            'frames_invalid_FITS',
            'frames_nonexistent',
        ]

        # File name for saved output
        for key in keys:
            fname = os.path.join(self.path('coord'), key + '.dat')
            self.add_file(key=key, value=fname)

        # Load a file or calculate again?
        if os.path.isfile(self.file('frames')):
            # self.frames = pd.read_csv(self.file('frames'))
            self.frames = np.loadtxt(self.file('frames'), dtype=str)# .tolist()

        else:

            self.frames_nonexistent = np.array([])
            self.frames_invalid_FITS = np.array([])
            self.frames_coadded = np.array([])
            self.frames = np.array([])
            for i in range(self.fields.shape[0]):

                coadded = False

                field = self.fields.iloc[i:i + 1]
                filename = frame_path(field)

                # Skip CO-ADDED frames and non-existent files
                if field.run in (106, 206):
                    self.frames_coadded = np.append(
                        self.frames_coadded, filename)
                    coadded = True

                if not os.path.isfile(filename):

                    self.frames_nonexistent = np.append(
                        self.frames_nonexistent, filename)

                else:

                    try:
                        # Does it open?
                        hdus = fits.open(filename)

                    except IOError:
                        # If not, then it is invalid FITS
                        self.frames_invalid_FITS = np.append(
                            self.frames_invalid_FITS, filename)

                    else:
                        if coadded:
                            self.frames_coadded = np.append(
                                self.frames_coadded, filename)
                        else:
                            # It is valid fits AND not a coadded frame.
                            # Let's use it
                            self.frames = np.append(self.frames, filename)

                    finally:
                        # This check turned out to be necessary, when running
                        # the program on the whole data set. Don'ø't know why.
                        # Error is about the variable `hdus` not esxisting.
                        if 'hdus' in dir():
                            hdus.close()

            # Save the lists
            for key in keys:
                fname = self.file(key=key)
                with open(fname, 'w+') as fsock:
                    fsock.write('\n'.join(getattr(self, key)))

        fstr = 'Of {:d} covering fields, {:d} are valid'
        fstr += ' AND non-coadded FITS files'
        # msg(fstr.format(self.fields.shape[0], len(self.frames)))
        msg(fstr.format(self.fields.shape[0], self.frames.size))

    # @logged
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

        Dependencies
        ------------
        self.get_consistent_frames() needs to have been run before this step.

        """

        ofname_cutouts = os.path.join(self.path('coord'), 'cutouts.csv')
        ofname_pxmax = os.path.join(self.path('coord'), 'pxmax.dat')
        # ofname_pxcrd = os.path.join(self.path('coord'), 'pxcrd.csv')
        # ofname_pxcrd_all = os.path.join(self.path('coord'), 'pxcrd_all.csv')

        if (os.path.isfile(ofname_cutouts)
            and
            # os.path.isfile(ofname_pxmax)
            os.path.isfile(ofname_pxmax)):
            # and
            # os.path.isfile(ofname_pxcrd)
            # and
            # os.path.isfile(ofname_pxcrd_all)):

            def s2d(s):
                for fmt in ['%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S']:
                    try: return dt.datetime.strptime(s, fmt)
                    except: continue

            # Load max cutout px size
            self.pxmax = list(np.loadtxt(ofname_pxmax))

            # Load cutout dates and filenames
            # dtype = [('cutout_dates', dt.datetime), ('cutouts', str)]
            dtype = [('cutout_dates', dt.datetime)]
            lkw = dict(converters={0:s2d}, delimiter=',', dtype=dtype)
            lkw.update(usecols=[0])
            data = np.loadtxt(ofname_cutouts, **lkw)
            [setattr(self, name, data[name]) for name in data.dtype.names]

            # Get the cutout filenames
            # I have not solved the problem of why the filenames can not be
            # loaded, so I am getting them in this way instead.
            self.cutouts = self.get_filenames(path='fits', match='*.fit.gz')

            # # Load row and col pixel indices
            # # cutouts only
            # ds = ['row_ixes', 'col_ixes']
            # ts = [int] * 2
            # dtype = [(d, t) for (d, t) in zip(ds, ts)]
            # lkw = dict(delimiter=',', dtype=dtype)
            # data = np.loadtxt(ofname_pxcrd, **lkw)
            # [setattr(self, name, data[name]) for name in data.dtype.names]
            # # self.row_ixes_all,
            # # self.col_ixes_all = [data[name] for name in data.dtype.names]

            # # Load row and col pixel indices
            # # ALL frames
            # ds = [
            #     'row_ixes_all', 'col_ixes_all',
            #     'row_ixes_all_r', 'col_ixes_all_r'
            # ]
            # ts = [int] * 4
            # dtype = [(d, t) for (d, t) in zip(ds, ts)]
            # lkw = dict(delimiter=',', dtype=dtype)
            # data = np.loadtxt(ofname_pxcrd, **lkw)
            # [setattr(self, name, data[name]) for name in data.dtype.names]
            # # self.row_ixes_all,
            # # self.col_ixes_all,
            # # self.row_ixes_all_r,
            # # self.col_ixes_all_r = [data[name] for name in data.dtype.names]

        else:

            """Create the cutouts"""

            # Accounting
            # ----------

            # Which files where used?
            frames_used_ix = []

            # Locations of output files which were written successfully to disk.
            # Also used to generate the .HTML report that will show each image.
            self.cutouts = []

            # For the timeline
            self.cutout_dates = []

            # Filenames of files which the program failed to write to disk.
            cutouts_failed_writes = []

            # For finding out what the largest possible cutout size is which
            # allows for a descent number of cutouts.
            # What is the largest possible pixel side length possible for
            # each frame?
            # = min( min(y, x), min(data.shape[0] - y, data.shape[1] - x) ) * 2 - 1
            self.pxmax = []

            # I want to map out what is going on.
            # Facts: The frames picked out above, clearly cover the candidate.
            #        The pixel coordinates should therefore also at least be within
            #        the region covered by the image data.

            # I need to find out how to get the pixels to be within the image data
            # region.

            # Input:
            #   * the transformation in each frame.
            #   * the coordinate that the cutout has to be made around

            # Given the transformation,
            #   * what is the order of the returned pixel coordinates?
            #   * what is the expected order of the world coordinates
            #     that are fed into the conversion function?

            # Pixels
            # ------

            # All the returned pixel coordinates when the order of the input coord-
            # inates are consistently (RA, Dec), and the returned pixel coordinates
            # are consistently assumed to be (col_ix, row_ix)
            self.row_ixes_all = []
            self.col_ixes_all = []

            # All the returned pixel coordinates when the order of the input coord-
            # inates are consistently (RA, Dec), and the returned pixel coordinates
            # are consistently assumed to be (col_ix, row_ix)
            self.row_ixes_all_r = []
            self.col_ixes_all_r = []

            # Actual, used row and col indices of the cutouts that can be made
            self.row_ixes = []
            self.col_ixes = []

            # Clean up
            # --------

            # astropy.io.fits will not overwrite existing files,
            # so I delete them all before creating the cutouts.
            print 'astropy.io.fits can not overwrite existing files,'
            print 'so first remove any previously saved files ...'
            cmd_fits_clean = 'rm {}'.format(
                os.path.join(
                    self.path('fits'), '*.fit.gz'
                )
            )
            print cmd_fits_clean
            print ''
            os.system(cmd_fits_clean)

            print 'Load each image and create cutout when possible'
            print '(marked by \'*\' at the end)'
            # print 'In DS9 fpCs have East up, North right'
            for i, ifname_cutout in enumerate(self.frames):

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
                # No, the image is loaded into a numpy array, and the order of the
                # axes is row-first
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

                # Use reference world coordinates to retrieve the given reference
                # pixel coordinates
                reference_pixels = w.wcs_world2pix([[crval1, crval2]], 1)

                # What if we switch them around?
                # reference_pixels_r = w.wcs_world2pix([[crval2, crval1]], 1)

                # print 'reference_pixels:  ', reference_pixels
                # print 'reference_pixels_r:', reference_pixels_r

                # reference_pixels:   [[ 1024.5   744.5]]
                # reference_pixels_r: [[ 358421.2155787 -303921.6830727]]

                # Then one gets values like the ones that I am more or less consis-
                # tently getting.

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

                        # The input order of the world coordinates has to be
                        # (RA, Dec)

                        # Columns first
                        x, y = pixels[0, 0], pixels[0, 1]

                        # Save all indices for debugging
                        self.col_ixes_all.append(x)
                        self.row_ixes_all.append(y)

                        self.col_ixes_all_r.append(pixels_r[0, 0])
                        self.row_ixes_all_r.append(pixels_r[0, 1])

                    else:
                        # Still columns first
                        x, y = pixels_r[0, 0], pixels_r[0, 1]

                        # Save all indices for debugging
                        # Notice that the post fixes ('', '_r') are switched
                        self.col_ixes_all.append(x)
                        self.row_ixes_all.append(y)

                        self.row_ixes_all_r.append(pixels[0, 0])
                        self.col_ixes_all_r.append(pixels[0, 1])

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
                        self.col_ixes_all.append(x)
                        self.row_ixes_all.append(y)

                        self.col_ixes_all_r.append(pixels_r[0, 0])
                        self.row_ixes_all_r.append(pixels_r[0, 1])

                    else:
                        x, y = pixels_r[0, 0], pixels_r[0, 1]

                        # Save all indices for debugging
                        # Notice that the post fixes ('', '_r') are switched
                        self.col_ixes_all.append(x)
                        self.row_ixes_all.append(y)

                        self.row_ixes_all_r.append(pixels[0, 0])
                        self.col_ixes_all_r.append(pixels[0, 1])

                # Now we have the x (data col) and y (data row) coordinates for the
                # image data (in NumPy coordinates)
                # Calculate the maximum-possible-square-cutout side length
                max_pixel_side_length = min(
                    min(y, x),
                    min(
                        image.shape[0] - y,
                        image.shape[1] - x)
                ) * 2 - 1 # Has to be within the image data
                self.pxmax.append(max_pixel_side_length)

                # Print info for all frames
                print '{: >4d} {: >8s} (row, col) = ({: >4.0f}, {: >4.0f})'.format(
                    i, phd.get('CTYPE1'), y, x),

                # Recap: y is the row index
                # Recap: x is the col index
                # print 'ROW {: >4.0f} {: >7.2f} {: >4.0f}'.format(
                #     dy, y, image.shape[0] - dy),
                # print 'COL {: >4.0f} {: >7.2f} {: >4.0f}'.format(
                #     dx, x, image.shape[1] - dx),

                is_within_ypad = (0 + self.dy <= y <= image.shape[0] - self.dy)
                is_within_xpad = (0 + self.dx <= x <= image.shape[1] - self.dx)

                if is_within_ypad and is_within_xpad:

                    # Mark this to the screen
                    print '*'

                    # Append the pixel coordinates for plotting below
                    self.row_ixes.append(y)
                    self.col_ixes.append(x)

                    # Save the index for the included input file
                    frames_used_ix.append(i)

                    # Get cutout
                    cutout_slice = image[
                        y - self.dy:y + self.dy + 1,
                        x - self.dx:x + self.dx + 1
                    ]
                    # cutout_slice = image

                    try:
                        datetime_observed = fits_get_observation_time(
                            header=phd)

                    except FITSHeaderError as err:
                        print '  FITSHeaderError:', err.message
                        hdulist.close()
                        continue

                    # For one frame, `datetime_observed` had the year 1858,
                    # and so it seems that this frame had an offset of zero
                    # seconds.
                    #   It is not clear what causes this.
                    #   It the frame a coadded image?`(Do coadded frames have
                    #  zero offset?)
                    #   Is is because the observation time is simply missing?
                    #   What other reasons could there be?
                    try:
                        # Format the date for the output filenames
                        fmt = '%Y-%m-%dT%H:%M:%S'
                        date_str = dt.datetime.strftime(datetime_observed, fmt)
                        # date_str = dt.datetime.strftime(datetime_observed, '%s')

                    except ValueError as err:
                        print '  ValueError:', err.message,
                        print '; Date:', datetime_observed
                        hdulist.close()
                        continue

                    # Save the cutout date for the timeline
                    self.cutout_dates.append(datetime_observed)

                    # Save the cutout

                    # As .FITS
                    # --------
                    # ...
                    # Create new PrimaryHDU object
                    # Attach the changed image to it
                    # Add to the hdulist, change reference-CR+PX and so on...
                    # Ref: http://prancer.physics.louisville.edu/astrowiki/index.\
                    # php/Image_processing_with_Python_and_SciPy
                    hdulist[0].data = cutout_slice
                    hdulist[0].header.set(
                        'ORIGINAL', ifname_cutout[ifname_cutout.rfind('/') + 1:])

                    # Could be nice to have the world-coordinate extent saved in
                    # the cutout, but it is not needed at the moment.
                    if 0:
                        hdulist[0].header.set('RA max', )
                        hdulist[0].header.set('RA min', )
                        hdulist[0].header.set('Dec max', )
                        hdulist[0].header.set('Dec min', )

                    # hdulist[0].header.set('', '')
                    # Do we need all the testing above?
                    # astropy.io.fits is loading the FITS image data in a NumPy
                    # array which is indexed row-first,
                    # whereas astropy.wcs does what it does independently of how
                    # the data are being loaded by astropy.io.fits.
                    # This means that we should always expect the return order of
                    # the pixel coordinates from astropy.wcs to be
                    # (column index, row index), and that only the input order
                    # needs to be checked for.
                    # Assuming that the input order is all that one needs to know
                    # in order to correctly
                    # set the reference coordinates below.
                    hdulist[0].header.set(
                        'CRPIX1', phd.get('CRPIX1') - (x - self.dx))
                    hdulist[0].header.set(
                        'CRPIX2', phd.get('CRPIX2') - (y - self.dy))
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
                        cutouts_failed_writes.append(i)
                    else:
                        self.cutouts.append(fname_fits)

                    # As .PNG
                    # -------

                    # re-scale intensities (should only be done for the .png
                    # images)
                    cutout_slice = lingray(cutout_slice)
                    # cutout_slice = lingray(np.log10(cutout_slice))

                    # Normalise
                    cutout_slice = cutout_slice / cutout_slice.sum()

                    # Save the image
                    fname_png = date_str + '.png'
                    ofname_png = os.path.join(self.path('png'), fname_png)
                    misc.imsave(ofname_png, cutout_slice)

                # OK, the coordinate was not far enough from the edge of the image
                # to make a cutout
                elif y < 0 or x < 0:

                    # Some transformations gave this, print it out
                    print '  NEGATIVE pixel coordinate !!!'

                else:

                    # Do not mark anything
                    print ''

                hdulist.close()

            # --- END for loop
            msg('Creating plots for visual inspection')

            # Create plots for visual inspection
            plot_time_coverage(self.cutout_dates, opath=self.path('coord'))
            plot_possible_cutouts(self.pxmax, opath=self.path('coord'))
            # plot_pixel_indices(self, opath=self.path('coord'))

        # --- END if cutouts already made

        # Recast variables
        # ----------------

        # Gather the output in a structured array so that they can be
        # sorted together after the observation date.

        # Rationale: When loading the data from disk, the cutout filenames
        # are sorted by filename which is names after the observation date.
        # This approach ensures that entries in either array matches the same
        # cutout.

        if 0:
            ts = [str, dt.datetime]
            ds = ['cutouts', 'cutout_dates']
            dtype = [(d, t) for (d, t) in zip(ds, ts)]
            data = np.zeros(len(self.cutouts), dtype=dtype)
            data_ = [self.cutouts, self.cutout_dates]
            for d in data_:
                print len(d)
            # print data.shape, data.dtype
            # print data
            # for i, field in enumerate(ds):
            #     print i, field
            #     data[field] = data_[i]
            #     print data[field]
            data['cutouts'] = self.cutouts
            data['cutout_dates'] = np.array(self.cutout_dates, dtype=[dtype[1]])

            # Do the sorting
            data.sort(order='cutout_dates')

            # Re-assign the sorted arrays
            [setattr(self, name, data[name]) for name in data.dtype.names]

        elif 0:
            data =[]
            for tuple_ in zip(self.cutout_dates, self.cutouts):
                data.append(tuple_)
            data = sorted(data, key=lambda k: k[0])

            cutout_dates = []
            cutouts = []
            for (d, c) in data:
                cutout_dates.append(d)
                cutouts.append(c)

            self.cutouts = np.array(cutouts)
            self.cutout_dates = np.array(cutout_dates)

        else:
            self.cutouts = np.array(sorted(self.cutouts))
            self.cutout_dates = np.array(sorted(self.cutout_dates))

        # This array is just to make histograms
        # so the sort order does not matter.
        self.pxmax = np.array(self.pxmax)

        # Summary
        # -------

        # Are the number of cutouts made above the desired threshold?
        if self.cutouts.size >= self.MIN_NUMBER_OF_CUTOUTS:
            msg('Succesfully created {:d} cutouts'.format(len(self.cutouts)))

            msg('Saving cutout data')

            # This is what I need for cutout_overview()

            # Save cutout dates and filenames
            with open(ofname_cutouts, 'w+') as fsock:
                fstr = '{:s},{:s}\n'
                for (c, d) in zip(self.cutouts, self.cutout_dates):
                    fsock.write(fstr.format(d.isoformat(), c))

            # Save max cutout px size
            with open(ofname_pxmax, 'w+') as fsock:
                fsock.write('\n'.join(np.array(self.pxmax).astype(str)))

            # # Save row and col pixel indices (cutouts only)
            # data = np.array([self.row_ixes, self.col_ixes]).T
            # np.savetxt(ofname_pxcrd, data, delimiter=',')

            # # Save row and col pixel indices (ALL frames)
            # data = np.array([
            #     self.row_ixes_all, self.col_ixes_all,
            #     self.row_ixes_all_r, self.col_ixes_all_r,
            # ]).T
            # np.savetxt(ofname_pxcrd_all, data, delimiter=',')

        else:
            msg_fstr = 'Not enough cutouts could be made'
            msg_fstr += ' (only {:d} of {:d} needed)'
            msg(msg_fstr.format(self.cutouts.size, self.MIN_NUMBER_OF_CUTOUTS))

            # Take this coordinate off the list.
            self.flag(reason='Not enough cutouts')

            """
            Must stop here and return a value to the program that created
            this instance and tell it to flag this coordinate and pick a
            new (non-flagged) one.

            No need to flag the ones that are inluded in tlist to begin with.

            """
            raise CutoutSequenceError('Not enough cutouts')

    ### Image processing and analysis on cutouts for the given coordinate

    @logged
    def gunzip(self):

        # Gunzip

        # Original command
        # From: https://superuser.com/questions/139419/how-do-i-gunzip-
        # to-a-different-destination-directory
        #
        #  for f in *.gz; do
        #  STEM=$(basename $f .gz)
        #  gunzip -c "${f}" > fit/"${STEM}"
        #  done
        #
        # Alternative
        # <http://docs.python.org/2/library/gzip.html>
        # (use Python for all the things!)
        #

        # abs_glob = os.path.join(self.path('fits'), '*.gz')
        cmd_kw = dict(fits=self.path('fits'), gunzip=self.path('gunzip'))
        cmd_gunzip = '''\
cd {fits}
for f in *.gz; do
STEM=$(basename $f .gz)
/bin/gunzip -c "${{f}}" > {gunzip}/"${{STEM}}"
done\
'''.format(**cmd_kw)

        print(cmd_gunzip)
        os.system(cmd_gunzip)

    @logged
    def wcsremap(self):

        # Dummy code for an idea
        # self.requires(methods=['gunzip'])
        # OR
        # @logged(requires=['gunzip'])
        # OR
        # dependencies(order=['gunzip'])
        # ---

        # Get sorted list of filenames in the gunzip directory
        files = self.get_filenames(path='gunzip', match='*.fit')

        try:
            # Set the template frame to be the first in the time series.
            template_frame = files[0]

        except IndexError:
            # wcsremap() was run before gunzip()
            # or gunzip() did not produce anything.
            print 'IndexError'
            print 'Can not proceed without reference to cutout frames'
            print 'Flagging sequence ...'
            self.flag(reason='wcsremap: no files found')

            if 0:
                # Make gunzip() runable again
                # a status other than 'finished' will force it to run again.
                print 'WCSREMAP: Forcing self.gunzip() to run again'
                self.log(step='gunzip', status='redo')

                msg('GUNZIP (called from CutoutSequence.wcsremap())')
                self.gunzip()

                try:
                    files = self.get_filenames(path='gunzip', match='*.fit')
                    template_frame = files[0]

                except:
                    print 'WCSREMAP: No files were found ...'
                    print 'WCSREMAP: Can not proceed. Flagging sequence ...'
                    self.flag(reason='wcsremap: no files found')

        except:
            print 'Failed for unknown reason ...'
            print 'Can not proceed without reference to cutout frames'
            print 'Flagging sequence ...'
            self.flag(reason='wcsremap: unknown reason')

        # Did instance get flagged above?
        if self.flagged:
            raise CutoutSequenceError('Flagged by CutoutSequence.wcsremap()')

        # Not specifying complete path, since environment should be set
        # on my own machine as well as on the image server.
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
            # print ifname
            # print ofname
            cmd_kw.update(ifname=ifname, ofname=ofname)
            cmd_wcsremap = cmd_wcsremap_fstr.format(**cmd_kw)

            # print(cmd_wcsremap)
            os.system(cmd_wcsremap)
        print '\n'

    # Service methods for loading and initialising the data for the analysis
    # ----------------------------------------------------------------------

    def get_filenames(self, path=None, match=None):
        """Get all files within directory specified by path."""
        import glob
        iglob = os.path.join(self.path(path), match)
        return sorted(glob.glob(iglob))

    def load_cutoutsio_wrapper(self):
        return self.select_usable(
            self.select_by_chosen_quality(
                self.load_cutoutsio()))

    @property
    def fname_cutoutsio(self):
        return os.path.join(self.path('coord'), 'cutoutsio.csv')

    def load_cutoutsio(self, verbose=False):
        # Define filename for previously/to-be saved output
        if verbose:
            print 'Loading', self.fname_cutoutsio, '...'
        return pd.read_csv(self.fname_cutoutsio)

    def select_by_chosen_quality(self, files):
        """Select the cutouts with the specified quality flags"""

        # Start out with none selected:
        qix = np.zeros(files.shape[0]).astype(bool)

        # Add if is of chosen quality:
        for q in self.load_quality:
            qix |= (files.quality == q)

        return files.iloc[qix]

    def select_usable(self, files):
        return files.iloc[(files.use == 1).values]

    def load_registered_cutouts(self):
        """
        Load all remapped images

        Assumptions
        -----------
        Using `files` assumes that no failed writings to disk occurred.
            (see Future below)

        Side effects
        ------------
        self._files : pandas.DataFrame
            Holds the list of filenames that were loaded after excluding
            filenames for which the file could not be loaded or the file
            quality did not match the desired list of qualities.

        self.cube_remap : 3D-array
            Data cube containing the remapped cutout frames.
        self._cube_remap : 3D-array
            Data cube containing the *unclipped* remapped cutout frames.

        """

        # list of files to load
        files = self.load_cutoutsio()
        files = self.select_by_chosen_quality(files)
        files = self.select_usable(files)

        # Save the full images so that view can be changed during runtime.
        self._cube_remap = np.zeros(self.size + (files.shape[0],))
        breakout = False
        for i in range(files.shape[0]):

            print '{: >3d}, '.format(i),
            if i and not (i + 1) % 10: print ''

            ifname, qflag, uflag = files.iloc[i]
            hdulist = fits.open(ifname)

            try:
                self._cube_remap[:, :, i] = hdulist[0].data

            except ValueError as err:
                print '\nValueError:', err.message
                print ifname.replace(env.paths.get('data'), '')
                print 'Flagging cutout and starting over ...'
                files.iloc[i:i + 1]['use'] = 0
                files.to_csv(self.fname_cutoutsio, index=False, header=True)
                breakout = True
                break

            finally:
                hdulist.close()

        print '\n'

        if breakout:
            self.load_registered_cutouts()

        else:

            # Attach `files` to the instance
            self._files = files
            # self._qualities = files.quality.values

            # Calibrate the data used for the subsequent analysis
            self.calibrate()

    def calibrate(self, quality=None, clip=None):
        """
        Calibrate the data used for analysis

        Parameters
        ----------
        quality : int, list, optional
            List of quality flags to use in the analysis.
            If None, the value set by self.quality is used.
        clip : int, scalar, optional
            How many pixels inwards to create the centered NumPy view of
            each frame in the cube of images to use for the analysis.
            If None, the value set by self.clip is used.

        """

        if quality is None:
            quality = self.quality

        msg('C A L I B R A T I O N :: quality == {}'.format(quality))
        # raise SystemExit

        if clip is None:
            clip = self.clip

        # Get indices for the data based on the list of files with load_quality
        qix = np.zeros(self._files.shape[0]).astype(bool)
        for q in self.quality:
            # Add if chosen
            qix |= (self._files.quality == q)

        # Save the stuff that is used in the analysis
        self.cube_remap = self.get_view(self._cube_remap[:, :, qix], clip=clip)
        self.files = self._files.iloc[qix]
        print 'self._files.shape[0] =', self._files.shape[0]
        print 'self.files.shape[0] =', self.files.shape[0]

    def get_view(self, data_cube, clip=0):
        """Creates a clipped view into the data cube."""
        Ny, Nx = self.size
        return data_cube[clip:Ny - clip, clip:Nx - clip, :]

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
        self.bg_models['mean'] = np.mean(self.cube_remap, axis=2)

    def calculate_residuals(self):
        """
        Calculate the residual images

        Side effects
        ------------
        self.template : float, 2D-array
            The chosen background model if available (None otherwise)
        self.cube_residuals : 3D-array

        """
        # Make sure that the template is updated
        self.calculate_background_models()

        # Subtract the background template
        self.cube_residual = self.cube_remap - self.template[:, :, None]

    # Service methods for the work-horse method `self.is_transient`
    # --------------------------------------------------------------

    # Single prediction
    # -----------------

    def predict(self, exp='any', **params):
        """Get prediction from running the given experiment once."""
        return getattr(self, 'experiment_' + exp)(**params)

    def experiment_any(self, sigma, tau):
        """
        Parameters
        ----------
        sigma : float
            scale for the LoG filter
        tau : float
            relative threshold

        Returns
        -------
        prediction : bool, scalar
            prediction

        """
        self.calculate_LoG(sigma=sigma)
        self.compare_neighbours()
        self.threshold_intensities(tau=tau)
        return np.any(self.cube_threshold)

    def experiment_many(self, sigma, tau):
        """
        Parameters
        ----------
        sigma : float
            scale for the LoG filter
        tau : float
            relative threshold

        Returns
        -------
        signals : int, vector
            signal vector of the same length as the cutout sequence
            with each entry counting the number of signals there were
            in the cutout frame at the same position.

        """
        self.calculate_LoG(sigma=sigma)
        self.compare_neighbours()
        self.threshold_intensities(tau=tau)

        # Build vector counting number of signals in each frame
        signals = self.cube_threshold.sum(axis=0).sum(axis=0)

        # Save copy on disk
        self.save_predictions(signals)

        return signals

    # Grid search
    # -----------

    def gridsearch(self, exp='any', **params):
        """Run grid search for given experiment"""
        return getattr(self, 'gridsearch_' + exp)(**params)

    def gridsearch_any(self, sigmas, taus):
        """Run grid search for experiment 'ANY'"""

         # Matrix of predictions (MoP)
        predictions = np.zeros((sigmas.size, taus.size)).astype(bool)

        # Sigmas
        for sigma_ix, sigma in enumerate(sigmas):
            # Set the scale of the objects that are looked for
            # No need to perform LoG every time I change between thresholds
            self.calculate_LoG(sigma=sigma)
            self.compare_neighbours()

            # Threshold values
            for tau_ix, tau in enumerate(taus):
                self.threshold_intensities(tau=tau)
                # Save the resulting prediction (bool)
                predictions[sigma_ix, tau_ix] = np.any(self.cube_threshold)

        # Save copy on disk
        self.save_predictions(predictions)

        return predictions

    def gridsearch_any2(self, sigmas, taus): self.gridsearch_any(sigmas, taus)

    def gridsearch_blr(self, sigmas, taus):
        """Baseline experiment using random predictions"""
        predictions = (np.random.random(size=(sigmas.size, taus.size)) >= .5)
        self.save_predictions(predictions)
        return predictions

    def gridsearch_bla(self, sigmas, taus):
        """Baseline experiment returning True for every parameter pair"""
        predictions = np.ones((sigmas.size, taus.size)).astype(bool)
        self.save_predictions(predictions)
        return predictions

    def gridsearch_bln(self, sigmas, taus):
        """Baseline experiment returning False for every parameter pair"""
        predictions = np.zeros((sigmas.size, taus.size)).astype(bool)
        self.save_predictions(predictions)
        return predictions

    def gridsearch_blr2(self, sigmas, taus): self.gridsearch_blr(sigmas, taus)
    def gridsearch_bla2(self, sigmas, taus): self.gridsearch_bla(sigmas, taus)
    def gridsearch_bln2(self, sigmas, taus): self.gridsearch_bln(sigmas, taus)

    # Grid search I/O
    # ---------------

    def save_predictions(self, predictions):
        """Save matrix of predictions from grid-search run to self._fname_gsp"""
        df = pd.DataFrame(data=predictions)
        df.to_csv(self._fname_gsp, index=False, header=False)

    @property
    def _fname_gsp(self):
        return self._fname_gs_prediction

    def set_fname_gsp(self, fname):
        self._fname_gs_prediction = os.path.join(self.path('results'), fname)
        return self._fname_gs_prediction

    @property
    def has_gs_prediction(self):
        return os.path.isfile(self._fname_gsp)

    @property
    def gs_prediction(self):
        return pd.read_csv(self._fname_gsp, header=None).values

    @property
    def gs_prediction_time(self):
        return os.path.getmtime(self._fname_gsp)

    # Scientific methods for the experiments
    # --------------------------------------

    def calculate_LoG(self, sigma):
        """
        Calculate the Laplacian-of-Gaussian-filtered images

        kwargs should contain at least on of the following parameters.

        Parameters
        ----------
        sigma : float
            The width of the LoG filter.

        Side effects
        ------------
        self.cube_LoG : float, 3D-array
            The LoG-filtered residuals.

        """
        msg('sigma = {: >5.2f} px'.format(sigma))
        self.cube_LoG = LoG2Dstack(self.cube_residual, sigma=sigma)

    def compare_neighbours(self):
        """
        Side effects
        ------------
        self.cube_minima_locs : float, 3D-array
            The data cube of bit fields for each 8-neighbour comparison.
            True : pixel location marks local minimum in LoG-filtered frames.
                This translates to a local maximum in the remapped frames.
            False : Not a local minimum of the LoG-filtered frames.

        """
        self.cube_minima_locs = minloc2Dstack(self.cube_LoG)

    def threshold_intensities(self, tau):
        """
        Threshold each image in sequence.

        Taking the logical conjunction with this binary fields marking the
        local minima of the LoG-filtered image and the original remapped image
        will result in an image where only the local minima have an intensity.
        These intensities will then serve in the (non-)detection of a signal
        given that they are above some defined threshold or not.

        Uses
        ----
        self.cube_remap
        self.cube_minima_locs


        Side effects
        ------------
        self.cube_threshold : float, 3D-array
            Data cube containing a bit field for each frame.
            True : Intensity of the local maximum in the remapped
                image is above threshold, i.e. the intensity is
                larger than the background. (Does it mean this?)
            False : Value was not larger than threshold.

            In other words:
                A boolean matrix where True means that a signal was found.

        """
        self.cube_threshold = threshold2Dstack(
            stack=self.cube_remap, stack_locs=self.cube_minima_locs, tau=tau)

    # END Scientific methods for the experiments

    @property
    def sigix(self):
        """
        Signal indices

        Returns
        -------
        None : NoneType
            If `cube_threshold` has not been created, return None instead.
        V : int, 1D-array
            Cutout indices for cutouts that have a signal in `cube_threshold`

        """
        if 'cube_threshold' not in dir(self):
            return None
        # Sum cube_threshold along both image axes and return indices of
        # entries in resulting vector with non-zero sum.
        return (self.cube_threshold.sum(axis=0).sum(axis=0) > 0).nonzero()[0]

    def report_signals(self):
        """
        Create a report of the cutouts in the given instance.

        """
        # plot_background_models(self)
        for i in self.sigix:
            plot_LoG_samples(self, offset=i)
            plot_residual_samples(self, offset=i)
            plot_binary_fields(self, offset=i)

    def __len__(self):
        if 'files' in dir(self):
            files = self.files
        else:
            files = self.load_cutoutsio_wrapper()
        return files.shape[0]

    def __getitem__(self, key):
        """
        Returns
        -------
        FD : 2D-array
            [
                [date_0, file_0],
                [date_1, file_1],
                [date_2, file_2],
                ...
                [date_n, file_n],

                dtype=object
            ]


        Python documentation
        --------------------
        Called to implement evaluation of self[key]. For sequence types,
        the accepted keys should be integers and slice objects. Note
        that the special interpretation of negative indexes (if the
        class wishes to emulate a sequence type) is up to the
        __getitem__() method. If key is of an inappropriate type,
        TypeError may be raised; if of a value outside the set of
        indexes for the sequence (after any special interpretation of
        negative values), IndexError should be raised. For mapping
        types, if key is missing (not in the container), KeyError should
        be raised.

        Note: for loops expect that an IndexError will be raised for
        illegal indexes to allow proper detection of the end of the
        sequence.

        """
        return self.cube_remap[key]

    def qix(self, q, which='remap'):
        _qix = (self.files.quality.values == q)
        return getattr(self, 'cube_' + which)[:, :, _qix]

    ## END class CutoutSequence definition ------------------------------------

# class Cutout(object):
#     def __init__(self, slicable):
#         self._slicable = slicable

#     def __getitem__(self, key):
#         return self._slicable[key]

#     def __setitem__(self, key, value):
#         self._slicable[key] = value

#     ## END class Cutout definition --------------------------------------------

def crd2str(radec):
    """Returns a unique string-representation of coordinate"""
    # return '{:010.5f}__{:010.5f}'.format(*radec).replace('.', '_')
    return '{:011.6f}__{:011.6f}'.format(*radec).replace('.', '_')


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
    fields = load_fields()

    # Get indices for fields that cover the chosen coordinate
    RA_geq_ramin = (fields.raMin <= radec_[0])
    RA_leq_ramax = (fields.raMax >= radec_[0])
    Dec_geq_decmin = (fields.decMin <= radec_[1])
    Dec_leq_decmax = (fields.decMax >= radec_[1])

    # The conjunction of these boolean arrays gives
    # the fields that satisfy the above need.
    ix = (RA_geq_ramin & RA_leq_ramax & Dec_geq_decmin & Dec_leq_decmax)

    return fields.loc[ix]


def fits_get_observation_time(header=None):
    """
    Extract observation time from given FITS header

    """

    if header is None:
        raise FITSHeaderError('Given header is None')

    # Formats used to read and write formatted dates
    fmts = [
        '%Y-%m-%dT%H:%M:%S',
        '%Y-%m-%d',
        '%y/%m/%d',
    #     '%H:%M:%S',
    ]
    tai = header.get('TAI', None)

    try:
        datetime_observed = tai2date(tai)

    except Exception:
        # Use DATE-OBS instead
        date_obs = header.get('DATE-OBS', None)
        time_obs = header.get('TAIHMS', None)
        datetime_obs = date_obs + 'T' + time_obs
        try:
            datetime_observed = dt.datetime.strptime(
                datetime_obs, fmts[0])

        except Exception:
            try:
                datetime_observed = dt.datetime.strptime(
                    date_obs, fmts[1])

            except Exception:
                try:
                    datetime_observed = dt.datetime.strptime(
                        date_obs, fmts[2])

                except Exception:
                    err = 'The header does not contain any'
                    err += ' readable observation timestamp.'
                    raise FITSHeaderError(err)

    return datetime_observed


def LoG2Dstack(stack, sigma):

    from scipy.ndimage import gaussian_laplace as gl

    # For each image in the stack obtain LoG(x, y; sigma)
    stack_LoG = np.zeros_like(stack)
    for i in range(stack.shape[2]):
        stack_LoG[:, :, i] = gl(stack[:, :, i], sigma=sigma)

    return stack_LoG


def minloc2Dstack(stack):
    """
    Find the local minima in each image of the stack of images (stack) and
    using eight-neighbour comparison.

    For each image in the stack, find local minima by comparing each pixel
    intensity with the values of all eight of its neighbours.

    The outer-most 1px-border of the image can not be a local minima (or
    maxima), since these pixels only have their inward neighbours. To be
    sure that these are not picked as local minima, the comparison is only
    made for the part of the image excluding the border pixels.

    The local minima are the True entries in the matrix `stack_minlocs`

    Parameters
    ----------
    stack : int or float, 3D-array
        A stack of 2D images.
        Stacking order is along the last dimension.

    Returns
    -------
    stack_minlocs : bool, 3D-array
        A stack of 2D binary fields, True where a local minimum was found.
        Stacking order is along the last dimension.

    """
    # Directions for 8-neighbour comparison
    directions = [
        ((None, -2), (None, -2), (None, None)),   # NW
        ((None, -2), (1, -1), (None, None)),      # N
        ((None, -2), (2, None), (None, None)),    # NE

        ((1, -1), (None, -2), (None, None)),      # W
        ((1, -1), (2, None), (None, None)),       # E

        ((2, None), (None, -2), (None, None)),   # SW
        ((2, None), (1, -1), (None, None)),      # S
        ((2, None), (2, None), (None, None)),    # SE

    ]
    slices = [
        (slice(*rows), slice(*cols), slice(*cutouts))
        for (rows, cols, cutouts)
        in directions
    ]

    # Slice indices for the centre of the image clipped 1px from each side.
    self_centre = (slice(1, -1), slice(1, -1), slice(None, None))

    # Indices for the 1px- border around the image
    borders = [
        # Top side
        (0, slice(None, None), slice(None, None)),
        # Right side
        (slice(None, None), -1, slice(None, None)),
        # Bottom side
        (-1, slice(None, None), slice(None, None)),
        # Left side
        (slice(None, None), 0, slice(None, None)),
    ]

    # Get the clipped images
    stack_centre = stack[self_centre]

    # Compare neighbours
    # ------------------

    # Alle entries are set to True to begin with.
    stack_minlocs = np.ones_like(stack).astype(bool)

    # Clear the border of all the images.
    for border in borders:
        stack_minlocs[border] = False

    # For each direction, take the logical conjunction between
    # the central rows and columns of the cube (the centre of
    # each image in the stack) and the truth values of the
    # comparison between values in the image centres and their
    # neighbours in the direction specified by the slice pairs.
    #   Having passed all 8 comparisons, the true local minima
    # are now located at the indices of the True entries.
    for slice_pair in slices:
        stack_minlocs[self_centre] &= (stack_centre <= stack[slice_pair])
    #     print slice_pair
    #     # for i in range(stack.shape[2]):
    #     i = 0
    #     print stack[:, :, i][slice_pair[:2]]
    #     # print stack_centre[:, :, i]
    #     print stack_centre[:, :, i] <= stack[:, :, i][slice_pair[:2]]
    #     print ''
    # raise SystemExit

    return stack_minlocs

def threshold2Dstack(stack, stack_locs, tau):
    """
    Threshold image stack.

    Method
    ------
    * For each image in the stack, get its maximum and minimum.
    * For each maximum and minimum, use relative threshold `tau` to calculate
      the absolute threshold based on the entire range of intensities.

        absolute_threshold = relative_threshold * (max - min)

    * Threshold intensities at given index locations.
    * Return stack of binary images, where True entries mark the indices for
      which the original intensity of the image was larger than the absolute
      threshold for that image.

    Notes
    -----
    Best way to threshold?
    Initial idea (KSP agreed after discussing it with him)
    t_cut = tau * (stack[:].max() - stack[:].min())
    A better way

    """

    # Keep only entries in the stack images where a location was specified.
    stack_maxima_only = (stack_locs * stack)

    if 1:

        """Base maxima and minima values on the whole image."""

        # Generate <N stack depth> long vector
        # of minima and maxima for each image.
        minima = stack.min(axis=0).min(axis=0)
        maxima = stack.max(axis=0).max(axis=0)

        # Now t_cut is a vector as well
        t_cut = tau * (maxima - minima)

    else:

        """Base maxima and minima values on specified locations only."""

        # Generate <N stack depth> long vector
        # of minima and maxima for each image.
        minima = stack_maxima_only.min(axis=0).min(axis=0)
        maxima = stack_maxima_only.max(axis=0).max(axis=0)

        # Now t_cut is a vector as well
        t_cut = tau * (maxima - minima)

    # Subtract the minimum value of each image from these entries,
    # so that their values lie within the range within which tau
    # was used to set the cut value for t_cut (<- unclear).
    stack_maxima_only -= minima[None, None, :]

    # Return boolean matrix where True means that a signal was found.
    return (stack_maxima_only > t_cut[None, None, :])



# ---

def write_cutout_sequence_summary(cs):

    summary = ''
    summary += 'SUMMARY\n'
    summary += '\nSky coordinate: (RA, Dec) = ({:4f}, {:4f})'.format(*cs.radec)
    summary += '\n\nInitial number of covering frames: {}'.format(cs.fields.shape[0])
    summary += '\nCo-added frames (left out for now): {}'.format(len(cs.frames_coadded))
    summary += '\nNon-existing file: {}'.format(len(cs.frames_nonexistent))
    summary += '\nValid files to use for processing cutouts: {}'.format(len(cs.files))
    summary += '\n\nCutout dimension (rows, cols) = ({}, {})'.format(*cs.size)
    summary += '\nPossible cutouts with this dimension: {}'.format(len(cs.cutout_files))
    summary += '\nFailed writings to disk: {}'.format(len(cs.failed_writes))
    # summary += '\n\nMean number of days between cutouts: {:.0f}d {:.3f}h'.format(
    #     cmdiff_mean // 1,
    #     (cmdiff_mean % 1) * 24
    # )

    ofname_summary = os.path.join(cs.path('coord'), 'summary.txt')
    with open(ofname_summary, 'w+') as fsock:
        fsock.write(summary)

    print summary

# ---

# Loaders
# -------

def load_fields():
    return pd.read_csv(env.files.get('fields'), sep=',')


def load_snlist():
    return pd.read_csv(env.files.get('snlist'), sep=';')


def load_gxlist():
    return pd.read_csv(env.files.get('gxlist'))


def load_tlist():
    return pd.read_csv(env.files.get('tlist'))


def load_fp2q():
    ifname = env.files.get('fp2q')
    try:
        return pd.read_csv(ifname, index_col=[0])
    except:
        print 'There was an error when trying to load', ifname, '...'
        return None


def load_all_cutout_sequences():
    """
    Returns
    -------
    cutout_sequences : np.array, 1D
        List of CutoutSequence objects.
    targets : np.array, 1D
        List of targets (feature labels) for each cutout.
        For now, just 1 feature : bool is_sn

    """

    snlist = load_snlist()
    gxlist = load_gxlist()

    params = dict(size=(101, 101))

    cutout_sequences = []
    targets = []

    # First take the supernovae
    params.update(is_sn=True)
    for i in range(snlist.shape[0]):
        params.update(radec=np.array([snlist.Ra[i], snlist.Dec[i]]))
        cutout_sequences.append(CutoutSequence(**params))

    # Then the galaxy coordinates
    params.update(is_sn=False)
    for i in range(gxlist.shape[0]):
        params.update(radec=np.array([gxlist.Ra[i], gxlist.Dec[i]]))
        cutout_sequences.append(CutoutSequence(**params))

    # Add all the targets
    targets.extend([True] * snlist.shape[0])
    targets.extend([False] * gxlist.shape[0])

    return np.array(cutout_sequences), np.array(targets)


def load_cutout_sequences():
    """
    Returns
    -------
    cutout_sequences : np.array, 1D
        List of CutoutSequence objects.
    targets : np.array, 1D
        List of targets (feature labels) for each cutout.
        For now, just 1 feature : bool is_sn

    """

    tlist = load_tlist()

    # print tlist.info()
    # print tlist.head()

    params = dict(size=(101, 101))

    cutout_sequences = []
    targets = []
    # coords = []
    for i in range(tlist.shape[0]):
        # print '{: >3d}'.format(i)

        # Create cutout object
        params.update(
            radec=np.array([tlist.Ra[i], tlist.Dec[i]]),
            is_sn=tlist.is_sn[i]
        )
        cutout_sequences.append(CutoutSequence(**params))
        targets.append(tlist.is_sn[i])
        # if not targets[-1]:
        #     coords.append(cutout_sequences[-1].crd_str)

    # print len(coords), np.unique(coords).size

    return np.array(cutout_sequences), np.array(targets).astype(bool)


# Cutout-sequence creation
# ------------------------

def update_tlist(css, targets):
    """Overwrite the existing tlist (having made the backup)"""

    import shutil

    # First save a copy of the original tlist, if not done already
    fname_tlist_backup = env.files.get('tlist') \
                            + '.backup.' + dt.datetime.now().isoformat()
    shutil.copyfile(env.files.get('tlist'), fname_tlist_backup)

    # Now that it is saved, the possibly changed list can be saved as tlist.
    # For now it seems easier to just change the original dataframe and save.
    tlist = load_tlist()
    for i, (cs, target) in enumerate(zip(css, targets)):
        Ra, Dec = cs.radec.flatten()
        tlist.iloc[i]['Ra'] = Ra
        tlist.iloc[i]['Dec'] = Dec
        tlist.iloc[i]['is_sn'] = int(target)

    # Write to disk
    tlist.to_csv(env.files.get('tlist'), index=False, headers=True)


def log_tlist_initialise():
    """
    Prepare a fresh log file for the latest run of create_cutout_data().

    """

    if os.path.isfile(env.files.get('log_tlist')):

        # Save a copy of the original tlist, if not done already

        import shutil

        fname_backup = env.files.get('log_tlist') \
                                + '.backup.' + dt.datetime.now().isoformat()
        shutil.copyfile(env.files.get('log_tlist'), fname_backup)

    columns = [
        'Ra', 'Dec', 'is_sn', 'N_fields', 'N_frames', 'N_cutouts',
    ]

    # Save the column names
    with open(env.files.get('log_tlist'), 'w+') as fsock:
        fsock.write('{}\n'.format(','.join(columns)))


def log_tlist_stats(cs=None):
    if cs is not None:

        if not os.path.isfile(env.files.get('log_tlist')):
            log_tlist_initialise()

        # Get the information
        Ra, Dec = cs.radec.flatten()
        is_sn = int(cs.is_sn)
        N_fields = cs.fields.shape[0]
        N_frames = len(cs.frames)
        N_cutouts = len(cs.cutouts)
        line = [
            str(d)
            for d in (Ra, Dec, is_sn, N_fields, N_frames, N_cutouts)
        ]

        with open(env.files.get('log_tlist'), 'a') as fsock:
            fsock.write('{}\n'.format(','.join(line)))


def create_cutout_data():
    """
    Initialise directories and create data for analysis.

    """

    # Create a new log file
    log_tlist_initialise()
    log = Log(env.files.get('log_create_cutout_data'))
    log.log_empty()
    log_while = ' while too_few_cutouts:'
    log_except = ' except:'
    log_while2 = ' while flagged and in_use:'
    log_passed = ' passed !'

    # Load the data
    css, targets = load_cutout_sequences()

    # There will be flagged cutouts sequences, so load these data sets
    # separately instead of every time a flag is encountered.
    tlist = load_tlist()
    gxlist = load_gxlist()
    # Replacement coordinates will only be drawn from gxlist
    # where all coordinates are non-events per definition
    params = dict(size=(101, 101), is_sn=False)

    # Run through each sequence
    for i, cs in enumerate(css):

        msg('', char=' ')
        msg('CutoutSequence[{:d}].init() - crd_str: {}'.format(i, cs.crd_str))

        log_base = 'css[{: >4d}]: for:'.format(i)
        log.log(log_base)

        too_few_cutouts = True
        while too_few_cutouts:

            log.log(log_base + log_while)

            try:
                cs.initialise()

            except CutoutSequenceError as err:
                msg('CutoutSequenceError: ' + err.message.upper())

                log.log(log_base + log_while + log_except)

                """
                cs is now flagged, and the flagging has been saved to the
                logfile for that coordinate.
                """

                # Begin log entry
                # Entry index of tlist.csv.backup
                line = ['{:d}'.format(i)]
                line += cs.radec.flatten().astype(str).tolist()
                line += [str(int(cs.is_sn))]
                line += [dt.datetime.now().isoformat()]
                with open(env.files.get('log_cut'), 'a') as fsock:
                    fsock.write('{}\n'.format(','.join(line)))

                # Now keep looking for a replacement for this entry
                # until a sequence can be made with the given coordinate.

                # Load both datasets so can check coordinate usage.
                N_gx = gxlist.shape[0]
                flagged = in_use = True
                while flagged and in_use:

                    log_str = log_base + log_while + log_except + log_while2
                    log.log(log_str)

                    RA, Dec = gxlist.iloc[np.random.randint(0, N_gx)]

                    # Check that the newly drawn coord is NOT in tlist already.
                    tix = (tlist.Ra.values == RA) & (tlist.Dec.values == Dec)
                    if tix.sum() > 0:
                        log_str += ' in_use ...'
                        log.log(log_str)
                        continue
                    else:
                        in_use = False

                    # Otherwise create instance
                    params.update(radec=np.array([RA, Dec]))
                    # Updating the variable cs
                    cs = CutoutSequence(**params)

                    # Was the instance previously flagged?
                    flagged = cs.flagged
                    log_str += ' flagged ...' if flagged else ' NOT flagged'
                    log.log(log_str)

            else:

                log.log(log_base + log_while + log_passed)
                too_few_cutouts = False

        # Outside of outer while loop, but inside for loop.
        # change the entry css
        css[i] = cs
        targets[i] = int(cs.is_sn)

        # Save the stats for this instance
        log_tlist_stats(cs)

    # Outside the outer loop
    # Save the possibly updated tlist to disk
    update_tlist(css, targets)

"""
Cutout sequences were made before I employed the quality flag in fields.csv.

The following two functions make it possible to map the quality of a field
and the cutout frame that was made from the data for this field.
"""

def create_cutout_original_to_field_quality_dict():
    """
    Build dictionary (YAML) for all the fields where the key is the
    basename of the filename for every possible frame (including invalid
    and non-exisiting), and the value is the quality of the field that
    corresponds to this frame.

    Side effects
    ------------
    framepath2quality.yaml : file
        The resulting dictionary created from the above operations and
        saved as a YAML document that can be used to generate list of
        cutout 2 quality mapping for each cutout sequence.

    """

    fields = load_fields()

    if 0:
        fp2q = {}
        for i in range(fields.shape[0]):

            print '{: >6d}'.format(i),
            if i and not ((i + 1) % 10):
                print ''

            field = fields.iloc[i:i + 1]
            filename = frame_path(field)
            basename = os.path.basename(filename)
            fp2q[basename] = int(fields.quality[i])

        # Make sure that the next output starts on a new line
        print '\n'

        # Save the data
        ykw = dict(default_flow_style=False, explicit_start=True, indent=2)
        with open(env.files.get('fp2q'), 'w+') as fsock:
            yaml.dump(fp2q, stream=fsock, **ykw)

    else:

        fnames = []
        qflags = []
        for i in range(fields.shape[0]):

            print '{: >6d}'.format(i),
            if i and not ((i + 1) % 10):
                print ''

            field = fields.iloc[i:i + 1]
            filename = frame_path(field)
            basename = os.path.basename(filename)
            fnames.append(basename)
            qflags.append(int(fields.quality[i]))
        print '\n'

        df = pd.DataFrame(data=qflags, columns=['quality'], index=fnames)
        df.to_csv(env.files.get('fp2q'), index=True, header=True)


def create_cutout2quality_mapping():
    """
    Run through every coordinate on tlist and load every raw cutout that
    was made. From the header of each of these files, obtain the quality
    flag by using the ORIGINAL header keyword as key in the dictionary
    above.

    Side effects
    ------------
    qualityflags

    """

    # tlist = load_tlist()
    css, targets = load_cutout_sequences()

    # Load framepath2quality dictionary
    fp2q = load_fp2q()

    if fp2q is None:
        print 'fp2q dictionary could not be loaded ...'
        raise SystemExit

    for i, cs in enumerate(css):

        print '{: >4d}'.format(i),
        if i and not ((i + 1) % 10):
            print ''

        # Define output filename for the records
        ofname = os.path.join(cs.path('coord'), 'cutoutsio.csv')

        # Get the filenames (to become one column in the DataFrame)
        filenames = cs.get_filenames(path='wcsremap', match='*.fit')

        # Get the corresponding quality flag from the previously created fp2q.
        qflags = []
        # In case a key error is caught, only save the filenames with key.
        fnames = []
        for fname in filenames:
            h = fits.getheader(fname)
            original = h['ORIGINAL']
            try:
                quality = int(fp2q.ix[original])
            except KeyError:
                print 'E!    ',
                continue
            qflags.append(quality)
            fnames.append(fname)
        print '\n'

        # Save the data in `ofname`
        uflags = np.ones_like(qflags).tolist()
        data = np.array([fnames, qflags, uflags]).T
        columns = ['fname', 'quality', 'use']
        df = pd.DataFrame(data=data, columns=columns)
        df.to_csv(ofname, index=False, header=True)

