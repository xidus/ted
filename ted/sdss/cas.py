#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sat 27 Apr 2013
#   Initial build.
#

"""
Scripts for dealing with the SDSS Catalogue Archive Server (CAS)

"""

import os
import sys
import time

"""
Will mechanize work with proxies?

I should create a wrapper or use socksipy to have all connections use the proxy server.
"""
import mechanize

from .. import env
from . import load_SNe_candidate_list

_path_data = env.paths.get('data')
_path_das = env.paths.get('das')
_path_cas = env.paths.get('cas')

_proxies = env.proxies


def get_galaxies():

    opath = os.path.join(_path_cas, 'galaxies')
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
    galaxies = query(sql_galaxies)
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


def get_fields():
    """For each coordinate in the SNe file, get the frames
    from the Field table that cover that coordinate.

    Saves each result in a seperate file named <SDSS_ID>.csv .

    """

    # Clean up local cas directory
    # field_clean_local_dir()
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

    df = load_SNe_candidate_list()
    SNe_len = df.count().max()

    print 'Beginning search through all confirmed SNe ...'

    time_total_beg = time.time()

    for i, (ra, dec, SDSS_id) in enumerate(zip(df.Ra, df.Dec, df.SDSS_id)):
        sql_fields = sql_fields_fstr.format(obj_ra=ra, obj_dec=dec)
        get_field_result(sql_fields, SDSS_id)
        s = 'Downloading: SN {: >4d} out of {: >4d} ...\r'.format((i + 1), SNe_len)
        sys.stdout.write(s)
        sys.stdout.flush()

    sys.stdout.write('\n')
    sys.stdout.flush()

    time_total_end = time.time()
    minutes = (time_total_end - time_total_beg) / 60.

    print 'Done downloading field catalogue ... in {:.0f} minutes'.format(minutes)
    # print 'Downloaded field catalogues for {} SNe'.format(i+1)


def get_field_result(sql_fields, SDSS_id):
    """Get query results from a single request for fields."""

    # Define output structure
    opath = os.path.join(_path_cas, 'fields')
    ofname = os.path.join(opath, '{}.csv'.format(SDSS_id))

    # Don't bother requesting anything if the file already exists
    if os.path.isfile(ofname):
        return

    # If the file does not exists, make sure that the path does
    if not os.path.exists(opath):
        os.makedirs(opath)

    # Request the data
    fields = query(sql_fields)
    if 'error' in fields.lower():
        sys.exit('ERROR: {}: CAS says something is wrong ...'.format(
            sys._getframe().f_code.co_name)
        )

    # And save it
    with open(ofname, 'w+') as fsock:
        fsock.write(fields)


def create_unique_field_list():
    """
    Creates a file containing all results from *get_fields()*
    and saves it in the CAS root directory.

    Looks through folder `fields` in the CAS root directory and loads
    all .csv files one at the time and adds the lines in each file to
    a list. This list is then sorted and only unique entries are kept
    before the list is saved in the CAS root directory.

    """

    ipath = os.path.join(_path_cas, 'fields')
    iglob = os.path.join(ipath, '*.csv')
    ofname = os.path.join(_path_cas, 'fields.csv')
    tfname = os.path.join(_path_cas, 'fields.tmp')

    # Clean up first, since the file is only appended to in the following
    if os.path.isfile(ofname):
        os.remove(ofname)

    commands = [

        # Build one big file with all the results
        'cat {iglob} >> {t}'.format(iglob=iglob, t=tfname),

        # Sort and remove duplicates
        'cat {t} | sort | uniq > {o}'.format(t=tfname, o=ofname),

        # Remove the temporayr file
        'rm {t}'.format(t=tfname),
    ]

    for cmd in commands:
        print cmd
        os.system(cmd)

    # Move last line (with the CSV headers) to the top
    with open(ofname, 'r') as fsock:
        lines = fsock.readlines()

    lines = [lines[-1]] + lines[:-1]

    with open(ofname, 'w') as fsock:
        fsock.write(''.join(lines))


def count_field_records():
    """
    Count the number for field records obtained for
    each SN, and save those numbers for later plotting.

    """

    ipath = os.path.join(_path_cas, 'fields')
    iglob = os.path.join(ipath, '*.csv')
    nfname = os.path.join(_path_cas, 'nfieldrecords.dat')

    commands = [

        # Count the number of records found for each coordinate
        # Including the header which is subtracted when showing th stats.
        'wc -l {iglob} > {n}'.format(iglob=iglob, n=nfname),

        # # Extract the number of results. Save output to temporary file.
        # 'sed "$d" < {n} > tempfile'.format(n=nfname),

        # # Move the files
        # 'mv {n} {n}.old'.format(n=nfname),
        # 'mv tempfile {n}'.format(n=nfname),
    ]

    for cmd in commands:
        print cmd
        os.system(cmd)

    # This is instead of running the last three commands in the commands list.
    # ... since they did not work, when run from Python.

    # Doing it in Python instead.
    import numpy as np
    nfieldrecords = np.loadtxt(nfname, usecols=[0], dtype=[('', int)])
    with open(nfname, 'w+') as fsock:
        fsock.write('\n'.join(nfieldrecords[:-1].astype(str)))


def query(sql_raw):
    """Sends SQL query to the CAS and returns the raw text of the response."""

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


def field_clean_local_dir():
    """
    Todo
    ----
    Make it possible to supply paths to folders that should have stuff removed.

    """

    opath = os.path.join(_path_cas, 'fields')
    if not os.path.exists(opath):
        return
    oglob = os.path.join(opath, '*')
    cmd = 'rm {}'.format(oglob)
    print cmd
    if not oglob == '/':
        os.system(cmd)
    else:
        print 'Clean-up command not executed ...'

