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


"""
Content
-------

__Functions:__

* iscoadded
* URL_exists
* load_SNe_candidate_list
* merge_sne_lists
* sql_fill_table_SNe

"""


###############################################################################
def iscoadded(run):
    """Test a single run number."""
    return run in (106, 206)


###############################################################################
def URI_exists(uri):
    # import requests
    response = requests.head(uri)
    return response.status_code in (200, 302)


###############################################################################
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
        # get_snlist()
        merge_sne_lists()
    return pd.read_csv(ifname, sep=';')


###############################################################################
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

    # String to value or NaN (hence the name valornan)
    snid_sortable = lambda SDSS_id: 'SN{:0>5d}'.format(int(SDSS_id[2:]))
    s2valornan = lambda s: s or np.nan
    conv = dict(SDSS_id=snid_sortable, Ra=ra2deg, Dec=dec2deg,
        redshift=s2valornan, Peak_MJD=s2valornan)
    df_left = pd.read_csv(fn_snlist_sdss_org, sep=';', converters=conv)
    df_right = pd.read_csv(fn_snlist_SNII_upd, sep=';', converters=conv)

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

    """Check for duplicate coordinates"""

    report = ''

    # Check for duplicate coordinate pairs
    coords = np.array([])
    for i in range(df.shape[0]):
        coords = np.append(coords, '{:014.9f}_{:014.9f}'.format(
            df.Ra.values[i], df.Dec.values[i])
        )
    ucoords, indices = np.unique(coords, return_inverse=True)
    report += 'Number of list entries:     {: >4d}\n'.format(df.shape[0])
    report += 'Number of unique entries:   {: >4d}\n'.format(ucoords.size)
    # print 'Number of unique entry IDs: {: >4d}'.format(np.unique(snlist.SDSS_id.values).size)

    # `indices` is an array of the same length as `coords`
    # An entry in `indices` is itself an index for an entry in `ucoords`,
    # i.e. `coords[i]` is `ucoords[indices[i]]`, `coords.size == indices.size`
    # Thus repeated entries in `indices` means that there were more than one
    # entry in `coords` which held this value; in this case a coordinate pair.
    duplicates = []
    for ix in np.unique(indices):
        if (indices == ix).sum() > 1:
            # There is a duplicate of this coordinate
            # Save the index for `ucoords`
            duplicates.append(ix)

    # Now we have the indices for the entries in `ucoords` whose values in
    # `coords` appear more than once. Let's retrieve the corresponding indices
    # for these values in `coords`.
    coord_indices = []
    for ix in duplicates:
        report += '\n'
        for i, uc in enumerate(ucoords[indices[indices == ix]]):
            if i == 0:
                # We only need to search for the indices of the duplicates
                # from one of the duplicates.
                coord_indices.append((coords == uc).nonzero()[0])
            # Report the actual coordinate strings for evry duplicate so that a
            # visual inspection can also verify that they are in fact congruent
            report += '{}'.format(uc)

    report += '\nIndices of `ucoords`: {}'.format(duplicates)
    report += '\nIndices of `coords`: {}'.format(repr(coord_indices))
    report += '\n'

    report += 'Entries from snlist:'
    for cices in coord_indices:
        report += '\n'
        report += '{}'.format(df.iloc[cices])

    # Selection of the entry in the list which gets to stay.
    # I choose the first occurrence, since the list is at this point already
    # sorted after SDSS_id which increases with time.
    # For this selection to be fair, the duplicate coordinates have to also
    # refer to the same object.
    # How do I know this?
    # I looked at the output of the previous for-loop, and each pair of entries
    # of the three duplicate coordinates had an estimated Modified Julian Date
    # Peak time which were separated in time by an interval that is of the
    # same order as the timespan in which a supernova is typically visible.
    # It is interesting that the survey found at least two peak dates for what
    # I, at least for now, assume is the same object. or two very closely
    # located objects. I do not know what the odds of spotting two different
    # events along the same line of sight and coming from two separate galaxies
    # within this relatively short amount of time; and considering three such
    # events within such a short list seems even less plausible.

    # POP the last of the two (or in principle more) entries from the list,
    # before saving it.
    # OR simply let Pandas do ALL the above work and also remove the duplicates.
    # report is still useful for the report, i.e. with the argument above.
    df.drop_duplicates(cols=['Ra', 'Dec'], inplace=True)

    with open(env.files.get('log_snlist'), 'w+') as fsock:
        fsock.write(report)

    if 0:
        """Add flag"""

        # ADDING A FLAG
        # Check if confirmed by IAUC, i.e. that it has an ID.
        confirmed = (df['IAUC_id'].values != np.nan)
        flags = np.zeros_like(df['IAUC_id'].values).astype(str)
        flags[confirmed] = 'C'
        flags[~confirmed] = ''
        df['Flag'] = flags

        print 'Confirmed SNe (has IAUC ID):', confirmed.sum()
        print 'Internally confirmed SNe:', (~confirmed).sum()

    df.to_csv(ofn_snlists_merged, sep=';', index=False, header=True)


###############################################################################
def sql_fill_table_SNe():
    """
    Loads merged list of SNe candidates and
    inserts the records into an SQLite3 db.
    """

    import sqlite3

    ofn_sqlite_db = env.files.get('db')
    fn_create_table_sne = os.path.join(_path_sql, 'skyml_create_table_SNe.sql')

    df = load_SNe_candidate_list()

    # # Read SQL
    with open(fn_create_table_sne, 'r') as fsock:
        sql_create_table_sne = fsock.read()

    # Connect to database file
    with sqlite3.connect(ofn_sqlite_db) as con:
        cur = con.cursor()
        cur.execute('DROP TABLE IF EXISTS Supernovae')
        cur.execute(sql_create_table_sne)

        sql_insert = '''\
INSERT INTO Supernovae
    (SDSS_id, SN_type, IAUC_id, Ra, Dec, redshift, Peak_MJD)
    VALUES
    (?, ?, ?, ?, ?, ?, ?)
'''
        # cur.executemany(sql_insert, df.values)
        for i, row in enumerate(df.values):
            cur.execute(sql_insert, row)

        con.commit()

        # df.to_sql(name='Supernovae', con=con, if_exists='replace')


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
def get_snlist():
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


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
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


