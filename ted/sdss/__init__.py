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

X   Make a csv-file with all unique frames. (bash script does this, but header has to be moved from the bottom to the top afterwards.)
X   Make a list of the number of recorded frames that covers each SN. (How did I make it before?)
X   Make CAS_get_fields() check how far the download is,
    by looking at what has already been downloaded.

*   Automate making a list of the relative paths to the downloaded images

*   <del>Make CAS_get_field_single() check if the file already exists</del>

"""

import os

import requests
import numpy as np
import pandas as pd

from .. import env
from ..parse import ra2deg, dec2deg

_path_data = env.paths.get('data')
_path_das = env.paths.get('das')
_path_cas = env.paths.get('cas')
_path_sql = env.paths.get('sql')

_proxies = env.proxies


def iscoadded(run):
    """Test a single run number."""
    return run in (106, 206)


def URI_exists(uri):
    # import requests
    response = requests.head(uri)
    return response.status_code in (200, 302)


def load_SNe_candidate_list():
    """
    Loads the SNe list from the CSV-file.

    Column titles:
    --------------
    SDSS SN Id : string
    Type : string
    IAUC Name : string
    Right Ascension (hh:mm:ss) : float
    Declination(dd:mm:ss) : float
    Redshift : float
    Peak MJD (approx) : float

    Returns
    -------
    pd.DataFrame object with the SNe

    """
    ifname = env.files.get('snlist')
    if not os.path.exists(ifname):
        # SDSS_get_snlist()
        merge_sne_lists()
    return pd.read_csv(ifname, sep=';')


def merge_sne_lists():
    """
    Reads in the two lists of SNe candidates and merges them into one
    bigger list (it turns out). This list is then having an extra column
    added with a special flag, and the final list is saved as .csv

    """

    fn_snlist_sdss_org = env.files.get('snlist_902')
    fn_snlist_SNII_upd = env.files.get('snlist_1030')

    if not (
        os.path.isfile(fn_snlist_sdss_org) and
        os.path.isfile(fn_snlist_SNII_upd)):
            raise SystemExit('''\
Source files do not exist!\
Please add them manually first.'''
            )

    ofn_snlists_merged = env.files.get('snlist')

    def col_count(df):
        # names = ['SDSS_id', 'SN_type', 'IAUC_id', 'redshift', 'Peak_MJD']
        # tests = ['', '', '', np.nan, np.nan]
        names = ['SN_type', 'IAUC_id', 'redshift', 'Peak_MJD']
        tests = ['', '', np.nan, np.nan]
        counts = [(df_right[name].values != test).sum() for name, test in zip(names, tests)]
        print ('{: >10s} ' * len(names)).format(*names)
        print ('{: >10d} ' * len(names)).format(*counts)
        print ''

    def snid_sortable(SDSS_id):
        _id = int(SDSS_id[2:])
        # print SDSS_id, _id
        return 'SN{:0>5d}'.format(_id)

    # String to value or NaN
    s2valornan = lambda s: s or np.nan
    conv = dict(SDSS_id=snid_sortable, Ra=ra2deg, Dec=dec2deg,
        redshift=s2valornan, Peak_MJD=s2valornan)
    df_left = pd.read_csv(fn_snlist_sdss_org, sep=';', converters=conv)
    df_right = pd.read_csv(fn_snlist_SNII_upd, sep=';', converters=conv)
    # df_left = pd.read_csv(fn_snlist_sdss_org, sep=';', converters=conv, index_col='SDSS_id')
    # df_right = pd.read_csv(fn_snlist_SNII_upd, sep=';', converters=conv, index_col='SDSS_id')

    SDSS_id_left = list(df_left.SDSS_id)
    SDSS_id_right = list(df_right.SDSS_id)
    SDSS_ids = np.array(SDSS_id_left + SDSS_id_right)
    SNids = np.unique(SDSS_ids)
    SNids.sort()

    columns = list(df_left.columns)
    # print columns
    # raise SystemExit
    ncols = len(columns) - 1
    df = pd.DataFrame(columns=columns)
    print df.describe()
    tests = [''] * 2 + [np.nan] * 4
    # raise SystemExit

    for SNid in SNids:
        # print SNid,
        lix = (df_left.SDSS_id == SNid)
        rix = (df_right.SDSS_id == SNid)
        ls = lix.sum()
        rs = rix.sum()

        # First, do both sides have the SDSS_id?
        # If not, just add the given one to the DataFrame.
        if ls and rs:
            # print 'is a duplicate ...'
            row_left = df_left[lix]
            row_right = df_right[rix]
            values = row_left.copy()
            # print values.values[0, 0], values.shape
            # raise SystemExit
            test_left = [(row_left[col].values[0] != test) for col, test in zip(columns, tests)]
            test_right = [(row_right[col].values[0] != test) for col, test in zip(columns, tests)]
            for i in range(1, ncols):

                col = columns[i]

                if test_left[i] and test_right[i]:
                    # Both have valid values
                    # Choose the value from the newest list
                    values.values[0, i] = row_right.values[0, i]

                elif test_right[i]:
                    values.values[0, i] = row_right.values[0, i]

                else:
                    values.values[0, i] = row_left.values[0, i]

            df = df.append(values)

        elif rs:
            # print 'is unique (right) ...'
            df = df.append(df_right[rix])

        elif ls:
            # print 'is unique (left) ...'
            df = df.append(df_left[lix])

    df.sort_index(by='SDSS_id', ascending=True)
    # print df.head(10)
    # raise SystemExit

    # ---

    # print df_left.head(10)
    # raise SystemExit

    # df = pd.merge(df_left, df_right, on=['SDSS_id'], sort=True)
    # col_count(df_right)
    # df = pd.merge(df_left, df_right, how='outer', sort=True)
    # col_count(df)

    # from IPython import embed
    # embed()

    # df = df.sort_index(by='SDSS_id', ascending=True)
    # df = pd.merge(df_left, df_right, how='outer')
    # df = pd.concat([df_left, df_right], join='outer', join_axes=[''], ignore_index=True, verify_integrity=True)
    # df = df_left.combine_first(df_right)
    # df_left.update(df_right)
    # df = df_left

    # col_count(df_right)
    # df_right.update(df_left)
    # col_count(df_right)
    # df = df_right

    # df = df.sort_index(by='SDSS_id', ascending=True)
    # df = pd.merge(df_left, df_right, on=df_left.columns.values, sort=True)
    # df = pd.concat([df_left, df_right], join='outer', ignore_index=True, verify_integrity=True)
    # df = df.sort_index(by='SDSS_id', ascending=True)

    # print df.head(10)
    # raise SystemExit

    # Check if confirmed by IAUC, i.e. that it has an ID.
    confirmed = (df['IAUC_id'].values != np.nan)
    flags = np.zeros_like(df['IAUC_id'].values).astype(str)
    flags[confirmed] = 'C'
    flags[~confirmed] = ''
    df['Flag'] = flags

    print 'Confirmed SNe (has IAUC ID):', confirmed.sum()
    print 'Internally confirmed SNe:', (~confirmed).sum()

    df.to_csv(ofn_snlists_merged, sep=';', index=False, header=True)


def sql_fill_table_SNe():
    """
    Loads merged list of SNe candidates and
    inserts the records into an SQLite3 db.
    """

    import sqlite3

    ofn_sqlite_skyml = env.files.get('skymldb')
    fn_create_table_sne = os.path.join(_path_sql, 'skyml_create_table_SNe.sql')

    df = load_SNe_candidate_list()

    # # Read SQL
    with open(fn_create_table_sne, 'r') as fsock:
        sql_create_table_sne = fsock.read()

    # Connect to database file
    with sqlite3.connect(ofn_sqlite_skyml) as con:
        cur = con.cursor()
        cur.execute('DROP TABLE IF EXISTS Supernovae')
        cur.execute(sql_create_table_sne)
        # from IPython import embed
        # embed()
        sql_insert = '''\
INSERT INTO Supernovae
    (SDSS_id, SN_type, IAUC_id, Ra, Dec, redshift, Peak_MJD, Flag)
    VALUES
    (?, ?, ?, ?, ?, ?, ?, ?)
'''
        # cur.executemany(sql_insert, df.values)
        for i, row in enumerate(df.values):
            # print i
            cur.execute(sql_insert, row)
        # con.commit()
        # df.to_sql(name='Supernovae', con=con, if_exists='replace')


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


from . import env


def wcs_world_order(w):
    """
    Prints out the return-order of world coordinates,
    when using the wcs.WCS.wcs_pix2world() method.


    Parameters
    ----------
    w : wcs.WCS object, required
        the wcs.WCS object created from a given header
        e.g.

            hdulist = fits.open('image.fits')
            w = wcs.WCS(hdulist[0].header)


    Returns
    -------
    None


    References
    ----------
    * astropy.wcs.WCS.wcs_pix2world.__doc__

    """

    lat_order = w.wcs.lat
    lng_order = w.wcs.lng
    lat_type = w.wcs.lattyp
    lng_type = w.wcs.lngtyp
    world_order = {}
    if lat_order < lng_order:
        world_order['first'] = ['lat', lat_type]
        world_order['last'] = ['lng', lng_type]
    else:
        world_order['first'] = ['lng', lng_type]
        world_order['last'] = ['lat', lat_type]

    for order, coord_info in world_order.iteritems():
        print '{: <5s}: {} / {}'.format(order, *coord_info)


