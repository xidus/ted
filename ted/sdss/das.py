#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sat 27 Apr 2013
#   Initial build.
#

"""
Scripts for dealing with the SDSS Data Archive Server (DAS)

"""

import os
import sys
import time

# if 'threading' in sys.modules:
#     raise Exception('threading module loaded before patching!')
#import gevent.monkey
#gevent.monkey.patch_all()
from gevent.pool import Pool

import requests
import numpy as np
import pandas as pd

from .. import env
from . import iscoadded
# from . import load_SNe

_path_data = env.paths.get('data')
_path_das = env.paths.get('das')
_path_cas = env.paths.get('cas')

_formatstrings = env.formatstrings

_proxies = env.proxies

_sdss_das = 'http://das.sdss.org/'

def fits_path_formatter(fmt='fpC'):
    """Returns correctly formatted URI"""
    pass


def create_download_URIs():
    """Returns a list/np.array of desired URIs."""
    pass


class Field(object):

    _filters = 'ugriz'

    def __init__(self, run, rerun, camcol, field):

        self.run = run
        self.rerun = rerun
        self.camcol = camcol
        self.field = field

    # As a property, I can use Field.field_parameters
    # instead of calling it explicitly as a method
    @property
    def field_parameters(self):
        """Return a dictionary with the field parameters"""
        return dict(
            run=self.run, rerun=self.rerun,
            camcol=self.camcol, field=self.field
        )

    def frame_parameters(self, filt=None):
        """
        Return a dictionary with the parameters
        used in the frame format string.

        """
        if filt not in self._filters:
            err_fstr = 'Argument `filt` must have one of the following values {} ...'
            raise FrameFormatError(err_fstr.format(','.join(list(self._filters))))
        frame_parameters = self.field_parameters
        frame_parameters.update(filt=filt)
        return frame_parameters

    def relative_path(self, frame_type='fpC'):
        """Returns correctly formatted URI"""
        return _formatstrings.get(frame_type).format(*self.parameters)

    def download(self, frame_type='fpC', only_if_not_exists=True):
        # Get online and local paths
        # Download and check file
        # Save frame information to database.
        pass

    # def local_file_checks_out(self, frame_type='fpC'):
    # def local_file_is_incorrupt(self, frame_type='fpC'):
    def local_file_is_readable(self, frame_type='fpC'):
        # Requires astropy.io.fits...
        # ... or should I just use zlib to try to unzip the file?
        import gzip
        try:
            ifname = path_local(self.relative_path(frame_type=frame_type))
            with gzip.open(ifname, 'rb') as fsock:
                pass
        except:
            return False
        else:
            return True

    def local_file_exists(self, frame_type='fpC'):
        pass

    def online_file_exists(self, frame_type='fpC'):
        pass

    def to_db(self):
        pass


def path_local(relative_URI, subdirs=[]):
    return os.path.join(_path_das, *subdirs + [relative_URI])


def path_online(relative_URI, subdirs=[]):
    return os.path.join(_sdss_das, *subdirs + [relative_URI])


def register_fields(frame_type='fpC', filt='r'):
    """

    """
    pass


# def register_frame(ifname):
def frame2db(run, rerun, camcol, field, filt, frame_type='fpC'):
    """
    Add the field to the
    """
    pass


def download_frames_by_radec(radec, frame_type='fpC'):
    """
    Requires passing queries to the CAS.
    """
    pass
    cas_root = os.path.join(_path_cas, 'radec')
    das_root = os.path.join(_path_das, frame_type)


###############################################################################
def download_frames_by_complete_result_list(frame_type='fpC', filt='r', pool_size=10):
    """
    Parameters
    ----------
    pool_size : int > 0
        the number of threads to use

    """

    URIs, ofnames = get_field_locations(
        ifname=env.files.get('fields'),
        frame_type=frame_type,
        filt=filt
    )

    # Sort and take only unique entries
    # But how? --> Same problem in download_frames_by_sn()
    # print np.unique(URIs_all).size, np.unique(ofnames_all).size
    # raise SystemExit

    # Download in parallel
    download_frames(URIs=URIs, ofnames=ofnames, pool_size=pool_size)


def download_frames_by_sn(bix=None, eix=None, frame_type='fpC', filt='r', pool_size=10):
    """
    Parameters
    ----------
    bix : int
        Start index for the array of query-result filenames
    eix : int
        End index for the array of query-result filenames
    pool_size : int > 0
        the number of threads to use

    """

    import glob

    # filters = 'ugriz'
    iglob = os.path.join(_path_cas, 'fields', '*.csv')
    fnames = sorted(glob.glob(iglob))[bix:eix]

    if fnames:
        URIs_all = []
        ofnames_all = []
        for ifname in fnames:

            URIs, ofnames = get_field_locations(
                ifname=ifname,
                frame_type=frame_type,
                filt=filt
            )
            URIs_all += URIs
            ofnames_all += ofnames

        # Sort and take only unique entries
        # But how? --> Same problem in download_frames_by_complete_result_list()
        # print np.unique(URIs_all).size, np.unique(ofnames_all).size
        # raise SystemExit

        # Download in parallel
        download_frames(URIs=URIs, ofnames=ofnames, pool_size=pool_size)


def get_field_locations(ifname=None, frame_type='fpC', filt='r'):
    """
    Loads query results from a .csv file and
    builds lists of local and online addresses.

    Parameters
    ----------
    ifname : string, required
        the location of the file containing the query results in .csv format.
    frame_type : string
        type of frame to download
        The URI format depends on it.
        Locally, fields will be saved in *_path_cas*/*frame_type*
    filt : str, single character
        Frame filter.

    """
    filters = 'ugriz'
    _sdss_das = 'http://das.sdss.org'
    _frame_path = os.path.join(_path_das, frame_type)
    _frame_fstr = _formatstrings.get(frame_type)
    ofnames = []
    URIs = []

    # Load data
    df = pd.read_csv(ifname, sep=',')
    dflen = df.shape[0]

    # Correct runs
    coadd_ix = (df.run.values == 106) + (df.run.values == 206)
    df.run.ix[coadd_ix] = (df.run.ix[coadd_ix] - 6) * 1000 + 6

    # Get local and online URIs for each result
    for i in range(dflen):
        parameters = dict(
            run=df.run.ix[i],
            rerun=df.rerun.ix[i],
            camcol=df.camcol.ix[i],
            field=df.field.ix[i],
            filt=filt
        )
        URIs.append(os.path.join(_sdss_das, _frame_fstr.format(**parameters)))
        ofnames.append(os.path.join(_frame_path, _frame_fstr.format(**parameters)))
    return URIs, ofnames


def download_frames(URIs=None, ofnames=None, pool_size=10):
    pool = Pool(pool_size)
    time_total_beg = time.time()
    for URI, ofname in zip(URIs, ofnames):
        print 'Downloading: Image {} ...'.format(URI)
        pool.spawn(download_URI, (URI, ofname))
    time_total_end = time.time()
    print 'Finished downloading {} frames in {:.0f} seconds'.format(
        len(URIs),
        time_total_end - time_total_beg
    )


def download_URI((URI, ofname)):

    response = requests.get(URI)

    if not os.path.isfile(ofname):

        opath = os.path.dirname(ofname)
        if not os.path.exists(opath):
            os.makedirs(opath)

        with open(ofname, 'wb+') as fsock:
            fsock.write(response.content)

    cmd = 'echo "{},{},{}" >> {}'.format(
        response.status_code, URI, ofname, env.files.get('log_das')
    )
    os.system(cmd)

    return response.status_code, URI, ofname


def frame_path(df_entry, frame_type='fpC', local=True, filt_ix=2):
    """
    Returns the path of the given image.

    Parameters
    ----------
    df_entry : single-entry panda.DataFrame

    """
    keys = ('run', 'rerun', 'camcol', 'field')
    filters = 'ugriz'

    fstr_frame = env.formatstrings.get(frame_type)

    if local:
        base_path = env.paths.get(frame_type)
    else:
        base_path = 'http://das.sdss.org'

    string_params = {key: df_entry[key].values[0] for key in keys}

    if iscoadded(string_params.get('run')):
        # Fix the value for the URI
        string_params['run'] = (string_params['run'] - 6) * 1000 + 6

    # Manually
    string_params['filt'] = filters[filt_ix]

    # Join base path and filename
    return os.path.join(base_path, fstr_frame.format(**string_params))


def check_field_list(do_get_frames=True):
    """
    Check the consistency of each file.

    Should preferably have been done before the cutout algorithm, and
    non-repairable and non-downloadable files should be flagged so that
    they can be separated out above when testing the other criteria.

    Side effects
    ------------
    files : list
        list of filenames for the covering and existing frames

    Steps
    -----

    Go through each offline frame location and
    verify the downloaded files' integrity by

        1.  Verifying that the file exists.
        2.  Opening it with astropy.io.fits
        3.  If file exists, but is corrupted
        4.  Try to download it once more
        5.  Repeat step 2. through 3.
        6.  If file is corrupted, exclude the corresponding
            field from the the fields.csv file.

    """

    import shutil

    fname = env.files.get('fields')
    fn_nonfpCs = os.path.splitext(fname)[0] + '_invalid_nonfpCs.csv'

    fields = pd.read_csv(fname)

    # Selection indices
    six = np.array([]).astype(bool)

    for i in range(fields.shape[0]):

        print '{: >7,d}'.format(i),

        # Note the small difference between the variable names here.
        do_get_frame = do_get_frames

        # field = fields.iloc[i]
        field = fields[i:i + 1]
        filename = frame_path(field)

        # If frames is coadded, leave it out of the final list for now
        if field.run in (106, 206):
            print u'CA', filename, '... Skipping ...'

            six = np.append(six, False)
            continue

        # If it does not exist, download
        if not os.path.isfile(filename):
            print u'!E', filename,

            if do_get_frame:
                print '... Trying to download ...',

                URI = frame_path(field, local=False)
                http_status, URI_, filename_ = download_URI((URI, filename))

                # Do we have it?
                if http_status is 200:
                    print '... Downloaded ...',

                    # No need to get it again, if integrity test fails.
                    do_get_frame = False

                else:
                    print '... Download failed', http_status, '... Skipping ...'

                    six = np.append(six, False)
                    continue

            else:
                print '... Not downloading ... Field kept ...'
                six = np.append(six, True)
                continue

        # At this point, it can be assumed that the file is on disk
        if do_get_frame:
            print u' E', filename,

        print '... Checking file integrity ...',

        if not is_FITS_file(filename):

            if not do_get_frame:
                print '... No more download tries ... Skipping'
                six = np.append(six, False)
                continue

            else:
                # Try to download one more time
                URI = frame_path(field, local=False)
                http_status, URI_, filename_ = download_URI((URI, filename))

                # Do we have it?
                if http_status is 200:
                    print '... Downloaded ...',

                else:
                    print '... Download failed', http_status, '... Skipping ...'
                    six = np.append(six, False)
                    continue

            # At this point, we have the file, but have not tested its integrity
            if is_FITS_file(filename):
                print 'VALID ...'
                six = np.append(six, True)
            else:
                print 'NOT VALID ...'
                six = np.append(six, False)

    # How many where excluded?

    # Now separate the invalid/non-existing files from the final list.
    fields_fpCs = fields.iloc[six]
    fields_nonfpCs = fields.iloc[~six]

    fields_fpCs.to_csv(fname, index=False, header=True)
    fields_nonfpCs.to_csv(fn_nonfpCs, index=False, header=True)


def is_FITS_file(ifname):
    """
    Check to see if `ifname` is a valid, i.e. readable, FITS file.

    """

    from astropy.io import fits

    try:
        hdus = fits.open(ifname)

    except IOError:
        retval = False

    else:
        retval = True

    finally:

        if 'hdus' in dir():
            hdus.close()

        return retval



###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
def DEPRECATED_download_fields_from_list(pool_size=10):
    """
    Downloads FITS files based on URIs in the file

        fpC_unique_URIs.txt

    """

    pool = Pool(pool_size)

    opath = os.path.join(_path_das, 'fpC')
    if not os.path.exists(opath):
        os.makedirs(opath)

    # NOTE THE DOUBLE PARENTHESIS !!!
    with open('fpC_unique_URIs.txt', 'r') as fsock:
        URIs = fsock.read().split('\n')

    N_images = len(URIs)

    time_total_beg = time.time()
    for i, URI in enumerate(URIs):
        fname = os.path.basename(URI)
        ofname = os.path.join(opath_full, fname)
        # print ofname
        # s = 'Downloading: Image {: >4d} out of {: >4d} ...\r'.format((i + 1), N_images)
        s = 'Downloading: Image {} ...\r'.format(URI)
        pool.spawn(download_URI, (URI, ofname))
        # download_URI(URI=URI, ofname=ofname)
        sys.stdout.write(s)
        sys.stdout.flush()
    sys.stdout.write('\n')
    sys.stdout.flush()

    time_total_end = time.time()
    print 'Done downloading field catalogue ... in {:.0f} seconds'.format(
        time_total_end - time_total_beg
    )


###############################################################################
def DEPRECATED_get_FITS(parameter_dict=None, url_fmt='fpC', subdir=None):

    if parameter_dict is None:
        print 'ERROR: No parameter dict'
        raise SystemExit

    # from pprint import pprint as pp
    # pp(parameter_dict)

    url_formats = dict(
        fpC='http://das.sdss.org/imaging/{run:d}/{rerun:d}/corr/{camcol:d}/fpC-{run:06d}-{filt:s}{camcol:d}-{field:04d}.fit.gz'
    )

    if url_formats.get(url_fmt) is None:
        print 'ERROR: Unknown DAS-URL format'
        raise SystemExit

    url_DAS = url_formats[url_fmt].format(**parameter_dict)
    fname = os.path.basename(url_DAS)

    opath = os.path.join(*['.', url_fmt] + [subdir] if subdir is not None else [])
    ofname = os.path.join(opath, fname)

    if not os.path.exists(opath):
        print 'Creating directory {}'.format(os.path.realpath(opath))
        os.makedirs(opath)

    print 'Requesting URI {} ...'.format(url_DAS),
    response = requests.get(url_DAS)
    print response.status_code

    print 'Saving file {} ...'.format(ofname)
    with open(ofname, 'wb') as fsock:
        fsock.write(response.content)


###############################################################################
def DEPRECATED_export_fpC_URIs():
    """
    Looks through cas directory and loads all .csv files and for
    each image frame in each .csv file use its parameters to build an URI
    for the corresponding fpC frame to be downloaded from the DAS, and,
    finally, append this URI to a list of URIs

    When done appending URI strings, count the number of frames before
    and after removing possible duplicates.

    Save the unique-entries-only list in a file and upload it to the server
    where the files are to be downloaded to.

    On the server, the file list will be read and each FITS file will
    be downloaded.

    Estimate the time for downloading the files in serial.

    """

    import glob

    filters = 'ugriz'
    fstr_fpC = 'http://das.sdss.org/imaging/{run:d}/{rerun:d}/corr/{camcol:d}/fpC-{run:06d}-{filt:s}{camcol:d}-{field:04d}.fit.gz'
    URI_list = []

    ipath = os.path.join(_path_cas, 'fields')
    iglob = os.path.join(ipath, '*.csv')

    ofname = os.path.join(_path_cas, 'fpC_unique_URIs.txt')

    fnames = sorted(glob.glob(iglob))

    if fnames:

        count_total = len(fnames)
        for i, ifname in enumerate(fnames):

            # fname = os.path.splitext(os.path.basename(ifname))[0]
            # print 'Reading results for {} ...'.format(ifname)
            sys.stdout.write('Reading file {: >4d} out of {: >4d} ...\r'.format(i+1, count_total))
            sys.stdout.flush()

            df = pd.read_csv(ifname, sep=',')

            # Begin quick fix
            # Efficient?
            coadd_ix = (df.run.values == 106) + (df.run.values == 206)
            df.run.values[coadd_ix] -= 6
            df.run.values[coadd_ix] *= 1000
            df.run.values[coadd_ix] += 6
            # End quick fix

            for i in range(len(df)):
                parameter_dict = dict(
                    run=df.run.values[i],
                    rerun=df.rerun.values[i],
                    camcol=df.camcol.values[i],
                    field=df.field.values[i],
                    filt=filters[2]  # r
                )
                URI_list.append(fstr_fpC.format(**parameter_dict))

        # Done
        sys.stdout.write('\n')
        sys.stdout.flush()

        # Prepare screen output and content for saving.
        URI_list_size_before = len(URI_list)
        URI_list = sorted(list(set(URI_list)))
        URI_list_size_after = len(URI_list)

        print 'Total unique count:', URI_list_size_after
        print 'Duplicates:', URI_list_size_before - URI_list_size_after

        with open(ofname, 'w+') as fsock:
            fsock.write('\n'.join(URI_list))

