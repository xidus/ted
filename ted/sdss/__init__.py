#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sat 27 Apr 2013
#   Initial build.
#

"""
The purpose of this module is to provide helper functions for
downloading and handling SDSS imaging data for the Stripe 82 survey.

TODO
----

*   Make a csv-file with all unique frames. (bash script does this, but header has to be moved from the bottom to the top afterwards.)

*   Automate making a list of the relative paths to the downloaded images

*   Make a list of the number of recorded frames that covers each SN. (How did I make it before?)

*   Make CAS_get_fields() check how far the download is,
    by looking at what has already been downloaded.

*   <del>Make CAS_get_field_single() check if the file already exists</del>

"""

import os
import sys
import time

if 'threading' in sys.modules:
    raise Exception('threading module loaded before patching!')
    # print '\030[1;31mthreading module loaded before patching!\030[0m'
import gevent.monkey
gevent.monkey.patch_all()
from gevent.pool import Pool

import requests
import numpy as np
import pandas as pd

# from IPython import embed

from .. import env
# from ..parse import (
#     dec2deg,
#     ra2deg,
# )

_path_data = env.paths['data']


def URI_exists(uri):
    # import requests
    response = requests.head(uri)
    return response.status_code in (200, 302)


def DAS_create_download_URIs():
    """Returns a list/np.array of desired URIs."""
    pass


def DAS_export_fpC_URIs():
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

    # global (
    #     _path_data
    #     _path_data
    # )

    import glob

    filters = 'ugriz'
    fstr_fpC = 'http://das.sdss.org/imaging/{run:d}/{rerun:d}/corr/{camcol:d}/fpC-{run:06d}-{filt:s}{camcol:d}-{field:04d}.fit.gz'
    URI_list = []

    ipath = os.path.join(_path_data, 'cas', 'fields')
    iglob = os.path.join(ipath, '*.csv')

    ofname = os.path.join(_path_data, 'generated', 'fpC_unique_URIs.txt')

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


def DAS_download_fields_from_list(pool_size=10):
    """
    Downloads FITS files based on URIs in the file

        fpC_unique_URIs.txt

    """

    pool = Pool(pool_size)

    # embed()
    # sys.exit()

    opath = os.path.join(_path_data, 'das', 'fpC')
    if not os.path.exists(opath):
        os.makedirs(opath)

    # NOTE THE DOUBLE PARENTHESIS !!!
    def DAS_get_URI((URI, ofname)):
        if URI is None or ofname is None:
            return

        if os.path.exists(ofname):
            return

        response = requests.get(URI)
        with open(ofname, 'wb+') as fsock:
            fsock.write(response.content)

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
        pool.spawn(DAS_get_URI, (URI, ofname))
        # DAS_get_URI(URI=URI, ofname=ofname)
        sys.stdout.write(s)
        sys.stdout.flush()
    sys.stdout.write('\n')
    sys.stdout.flush()

    time_total_end = time.time()
    print 'Done downloading field catalogue ... in {:.0f} seconds'.format(
        time_total_end - time_total_beg
    )


def DAS_download_fields(in_parallel=False, pool_size=10):

    import glob

    filters = 'ugriz'

    ipath = 'cas'
    fnames = sorted(glob.glob(os.path.join(ipath, '*.csv')))

    if fnames:
        # print 'Using first result from file {} ...'.format(fnames[0])
        # fieldID,run,rerun,camcol,field,raMin,raMax,decMin,decMax
        for ifname in fnames[:1]:  # For now only take the first SN

            fname = os.path.splitext(os.path.basename(ifname))[0]
            print 'Reading results for {} ...'.format(ifname)
            df = pd.read_csv(ifname, sep=',')

            for i in range(len(df)):
                parameter_dict = dict(
                    run=df.run.values[i],
                    rerun=df.rerun.values[i],
                    camcol=df.camcol.values[i],
                    field=df.field.values[i],
                    filt=filters[2]  # r
                )
                DAS_get_FITS(parameter_dict=parameter_dict, subdir=fname)


def DAS_get_FITS(parameter_dict=None, url_fmt='fpC', subdir=None):

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


def CAS_get_galaxies():

    opath = os.path.join(_path_data, 'cas', 'galaxies')
    ofname = os.path.join(opath, 'galaxies.csv')

    sql_galaxies = """\
SELECT
    objID,
    fieldID,
    -- psfMag_r,
    -- specobjid,
    run, rerun, camcol, field
    -- ra, dec, b, l
    -- raErr, decErr, raDecCorr,

FROM
    PhotoObjAll

WHERE
    run in (106, 206)
  AND
    psfMag_r < 20.0
  AND
    type = 3
"""
    time_total_beg = time.time()
    galaxies = CAS_query(sql_galaxies)
    print 'Downloading galaxy objects from PhotoObjAll'
    time_total_end = time.time()
    print 'Query executed in {:.0f} seconds'.format(
        time_total_end - time_total_beg
    )
    if 'error' in galaxies.lower():
        print('ERROR: {}: CAS says something is wrong ...'.format(
            sys._getframe().f_code.co_name)
        )
    with open(ofname, 'w+') as fsock:
        fsock.write(galaxies)


def CAS_get_fields():
    """For each coordinate in the SNe file, get the frames
    from the Field table that cover that coordinate.

    Saves each result in a seperate file named <SDSS_ID>.csv .

    """

    # Clean up local cas directory
    # CAS_field_clean_local_dir()
    # Do it manually forom the main program instead
    # This way, if the program is interrupted, it does not need to begin all
    # over, since the single-field-getter skips requests that have already been made.

    # with open('template_fields_that_cover_SN.sql', 'r') as fsock:
    #     sql_fields_fstr = fsock.read()

    sql_fields_fstr = """\
SELECT
    fieldID,
    -- skyVersion,
    run, rerun, camcol, field,
    -- nObjects,
    -- numStars_r,
    -- nCR_r,
    -- nBrightObj_r,
    -- nFaintObj_r,
    quality,        -- Quality of field in terms of acceptance
    mjd_r,          -- Julian Date when row 0 was read
    -- Do I need Astrometric transformation constants?
    -- Do I need Airmass measurements?
    raMin, raMax, decMin, decMax

FROM
    -- TARGStripe82..Field
    Field

WHERE
    raMin < {obj_ra} AND {obj_ra} < raMax
  AND
    decMin < {obj_dec} AND {obj_dec} < decMax

ORDER BY
    run ASC,
    rerun ASC,
    camcol ASC,
    field ASC,
    raMin ASC,
    decMin ASC
"""

    df = load_SNe()
    SNe_len = df.count().max()

    print 'Beginning search through all confirmed SNe ...'

    time_total_beg = time.time()

    for i, (ra, dec, SDSS_id) in enumerate(zip(df.Ra, df.Dec, df.SDSS_id)):
        sql_fields = sql_fields_fstr.format(obj_ra=ra, obj_dec=dec)
        CAS_get_field_single((sql_fields, SDSS_id))
        s = 'Downloading: SN {: >4d} out of {: >4d} ...\r'.format((i + 1), SNe_len)
        sys.stdout.write(s)
        sys.stdout.flush()

    sys.stdout.write('\n')
    sys.stdout.flush()

    time_total_end = time.time()
    minutes = (time_total_end - time_total_beg) / 60.

    print 'Done downloading field catalogue ... in {:.0f} minutes'.format(minutes)
    # print 'Downloaded field catalogues for {} SNe'.format(i+1)


def CAS_get_field_single((sql_fields, SDSS_id)):
    """Get query results from a single request for fields."""

    # Define output structure
    opath = os.path.join(_path_data, 'cas', 'fields')
    ofname = os.path.join(opath, '{}.csv'.format(SDSS_id))

    # Don't bother requesting anything if the file already exists
    if os.path.isfile(ofname):
        return

    # If the file does not exists, make sure that the path does
    if not os.path.exists(opath):
        os.makedirs(opath)

    # Request the data
    fields = CAS_query(sql_fields)
    if 'error' in fields.lower():
        sys.exit('ERROR: {}: CAS says something is wrong ...'.format(
            sys._getframe().f_code.co_name)
        )

    # And save it
    with open(ofname, 'w+') as fsock:
        fsock.write(fields)


def CAS_query(sql_raw):
    """Sends SQL query to the CAS and returns the raw text of the response."""

    import mechanize

    form_url = 'http://cas.sdss.org/stripe82/en/tools/search/sql.asp'
    return_format = 'csv'

    sql_filtered = ''
    for line in sql_raw.split('\n'):
        sql_filtered += line.split('--')[0] + ' ' + os.linesep

    # Browser
    br = mechanize.Browser()

    # User-Agent
    br.addheaders = [
        (
            'User-agent',
            'Mozilla/5.0 (X11; Ubuntu; Linux i686; rv:21.0) Gecko/20100101 Firefox/21.0'
        )
    ]

    br.open(form_url)
    br.select_form(nr=0)

    # User credentials
    br.form['cmd'] = sql_filtered
    br.form['format'] = [return_format]

    # Search and return the content of the resulting CSV-file
    return br.submit().read()


def CAS_field_clean_local_dir():

    opath = os.path.join(_path_data, 'cas', 'fields')
    if not os.path.exists(opath):
        return
    oglob = os.path.join(opath, '*')
    cmd = 'rm {}'.format(oglob)
    print cmd
    if not oglob == '/':
        os.system(cmd)
    else:
        print 'Clean-up command not executed ...'


def load_SNe():
    """
    Loads the SNe list from the CSV-file.

    Column titles:
    --------------

    * SDSS SN Id : string
    * Type : string
    * IAUC Name : string
    * Right Ascension (hh:mm:ss) : float
    * Declination(dd:mm:ss) : float
    * Redshift : float
    * Peak MJD (approx) : float

    Returns
    -------

    pd.DataFrame object with the SNe

    """

    ifname = env.files['snlist']

    if not os.paths.exists(ifname):
        SDSS_get_snlist()

    return pd.read_csv(ifname, sep=';')


def SDSS_get_snlist():
    """NOT IMPLEMENTED YET."""

    return

    # The best way to do this seems to be to grab the HTML table and manipulate it manually.

    # Download the lists as HTML tables
    # Edit their source codes to get .csv format.
    # Use the merge function to merge the lists into the big master list.
    # The name of the list is currently given by

    if 0:
        # List of
        URI_snlist = 'http://www.sdss.org/supernova/snlist.dat'
        fname = os.path.basename(URI_snlist)
        # File-name path to local copy of snlist.dat
        ofname = os.path.join(_path_data, fname)
        # File-name path to local copy of snlist.csv (to be loaded from Pandas)
        ifname = os.path.join(_path_data, os.path.splitext(fname)[0] + '.csv')

        if 1 or not os.path.exists(ifname):

            if not os.path.exists(ofname):
                print 'Downloading {} ... '.format(URI_snlist)
                response = requests.get(URI_snlist)
                print 'HTTP response: {} ...'.format(response.status_code)
                print 'Saving file {} ...'.format(ofname)
                with open(ofname, 'w+') as fsock:
                    fsock.write(response.content)

            else:
                print 'Downloaded file {} found ...'.format(fname)

            print 'Loading {} to create csv-file ...'.format(fname)
            with open(ofname, 'r') as fsock:
                sncontent = fsock.read()
            print 'Replacing spaces entries with semicolons ...'
            sncontent = sncontent.replace('   ', ';').replace('---', '')
            print 'Modifying column names ...'
            snlist = sncontent.split('\n')
            snlist[0] = 'SDSS_id;SN_type;IAUC_id;Ra;Dec;redshift;Peak_MJD'

            if 0:
                # Debug
                ncols = np.zeros(len(snlist))
                for i, line in enumerate(snlist):
                    ncols[i] = len(line.split(';'))
                    print i + 1, ncols[i]

                for i, cnt in enumerate(np.unique(ncols)):
                    print cnt, (ncols == cnt).sum()

            print 'Saving file {} ...'.format(ifname)
            with open(ifname, 'w+') as fsock:
                fsock.write('\n'.join(snlist))

        # return pd.read_csv(
        #     # ofname, sep=';', skiprows=2, converters=dict(Ra=ra2deg, Dec=dec2deg)
        #     ifname, sep=';', converters=dict(Ra=ra2deg, Dec=dec2deg)
        # )


