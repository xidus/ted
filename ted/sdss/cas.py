#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sat 27 Apr 2013
#   Initial build.
#

"""
Scripts for dealing with the SDSS
Catalogue Archive Server (CAS) and
the results obtained herefrom.

"""

import os
import sys
import time

# Will mechanize work with proxies?
# Maybe I should create a wrapper or use socksipy to have all connections use the proxy server.
import mechanize

# .. : ted
from .. import env
from ..time import format_HHMMSS_diff
from . import load_SNe_candidate_list

_path_data = env.paths.get('data')
_path_das = env.paths.get('das')
_path_cas = env.paths.get('cas')

_proxies = env.proxies

"""
Content
-------

__Functions:__

* get_galaxies
* create_galaxy_list
* plot_gxlist
* build_tlist

* get_fields
* get_field_result
* create_unique_field_list
* filter_invalid_from_unique_field_list
* count_field_records

* query

* field_clean_local_dir

"""


def get_galaxies():
    """
    This query takes too long?

    """

    ipath = env.paths.get('sql')
    fn_sql_galaxies = os.path.join(ipath, 'cas', 'stripe82_galaxies.sql')
    with open(fn_sql_galaxies, 'r') as fsock:
        sql_galaxies = fsock.read()
    print sql_galaxies
    print ''

    print 'Downloading galaxy objects from PhotoObjAll'
    time_total_beg = time.time()
    galaxies = query(sql_galaxies)
    time_total_end = time.time()
    print 'Query executed in {:.0f} seconds'.format(
        time_total_end - time_total_beg
    )
    if 'error' in galaxies.lower():
        print('ERROR: {}: CAS says something is wrong ...'.format(
            sys._getframe().f_code.co_name)
        )
        print '\nContent of returned result:\n\n{}'.format(galaxies)
    with open(env.files.get('galaxies'), 'w+') as fsock:
        fsock.write(galaxies)


def create_galaxy_list():
    """
    Create a list of galaxy coordinates `gxlist.csv` by
    filtering response from SDSS Skyserver in the following way:

        0. Exclude duplicate coordinates.
        1. Exclude coordinates that are not within covered stripe area.
        2. Exclude coordinates that are too close to any supernova coordinate.

    Returns
    -------
    None

    Side effects
    ------------
    File `gxlist.csv` in env.paths.get('cas')

    """

    from .cutouts import get_covering_fields

    import numpy as np
    import pandas as pd

    # Manipulation and analysis of geometric objects in the Cartesian plane.
    # from shapely.geometry import Polygon, Point
    # from shapely.ops import cascaded_union

    from IPython import embed

    print 'Loading data ...'

    cols = dict(usecols=['Ra', 'Dec'])
    fcols = dict(usecols=['raMin', 'raMax', 'decMin', 'decMax'])
    galaxies = pd.read_csv(env.files.get('galaxies'), sep=',', **cols)
    snlist = pd.read_csv(env.files.get('snlist'), sep=';', **cols)
    fields = pd.read_csv(env.files.get('fields'), sep=',', **fcols)

    # Step 0
    # ------

    print 'Step 0 : Exclude duplicate coordinates ...'

    gra = galaxies.Ra.values
    gdec = galaxies.Dec.values

    # coords = np.zeros_like(gra).astype(str)
    coords = []
    for i in np.arange(galaxies.shape[0]):
        # coords[i] = '{},{}'.format(gra[i], gdec[i])
        coords.append('{},{}'.format(gra[i], gdec[i]))
    u_coords, uix = np.unique(coords, return_index=True)

    print 'Reducing data for the next step ...'

    u_galaxies = galaxies.iloc[uix]
    u_gra = gra[uix]
    u_gdec = gdec[uix]

    # Step 1
    # ------

    print 'Step 1 : Exclude coordinates that are not within covered stripe area ...'

    ipath = os.path.dirname(env.files.get('gxlist'))
    ifname = os.path.join(ipath, 'cu_gxlist.csv')
    if not os.path.isfile(ifname):

        """This step takes ~40 minutes on my laptop."""

        ramin, ramax = fields.raMin.values, fields.raMax.values
        decmin, decmax = fields.decMin.values, fields.decMax.values

        N_gx = u_galaxies.shape[0]
        N_fields = fields.shape[0]

        # How many times can I have an array of <N_rows> rows by
        N_rows = 100
        N_gx_eee = N_gx // N_rows
        N_gx_rem = N_gx % N_rows

        cix = np.array([]).astype(bool)

        # Create vectors and matrices that are repeatedly used
        N_fields_ZEROS = np.zeros(N_fields)[None, :]

        RAMIN = np.zeros(N_rows)[:, None] + ramin[None, :]
        RAMAX = np.zeros(N_rows)[:, None] + ramax[None, :]

        DECMIN = np.zeros(N_rows)[:, None] + decmin[None, :]
        DECMAX = np.zeros(N_rows)[:, None] + decmax[None, :]

        for n in range(N_gx_eee):

            # How far are we?

            beg = n * N_rows
            end = (n + 1) * N_rows

            print 'n = {: >4d}; {: >6d}; {: >6d};'.format(n, beg, end)

            # Create matrices
            GRA = u_gra[beg:end, None] + N_fields_ZEROS
            GDEC = u_gdec[beg:end, None] + N_fields_ZEROS

            CMP = np.ones((N_rows, N_fields)).astype(bool)
            CMP &= (GRA > RAMIN)
            CMP &= (GRA < RAMAX)
            CMP &= (GDEC > DECMIN)
            CMP &= (GDEC < DECMAX)

            # Append the booleans to my master index file
            cix = np.append(cix, np.any(CMP, axis=1))

            # Clean up
            del GRA, GDEC, CMP

        if N_gx_rem > 0:

            # Finally, the remaining less than <N_rows> coordinates
            beg = (n + 1) * N_rows
            end = beg + N_gx_rem

            # Create matrices
            GRA = u_gra[beg:end, None] + N_fields_ZEROS
            GDEC = u_gdec[beg:end, None] + N_fields_ZEROS

            RAMIN = np.zeros(N_gx_rem)[:, None] + ramin[None, :]
            RAMAX = np.zeros(N_gx_rem)[:, None] + ramax[None, :]

            DECMIN = np.zeros(N_gx_rem)[:, None] + decmin[None, :]
            DECMAX = np.zeros(N_gx_rem)[:, None] + decmax[None, :]

            CMP = np.ones((N_gx_rem, N_fields)).astype(bool)
            CMP &= (GRA > RAMIN)
            CMP &= (GRA < RAMAX)
            CMP &= (GDEC > DECMIN)
            CMP &= (GDEC < DECMAX)

            # Append the booleans to my master index file
            cix = np.append(cix, np.any(CMP, axis=1))

            # Check
            print ''
            print 'N_gx =', N_gx
            print 'cix.size =', cix.size
            print 'cix.dtype =', cix.dtype

        # Embed so that I do not need to re-do this step again...
        embed()

        print 'Reducing data for the next step ...'

        cu_galaxies = u_galaxies.iloc[cix]
        cu_gra = u_gra[cix]
        cu_gdec = u_gdec[cix]

    else:

        print 'Step 1. already performed. Loading result ...'

        cu_galaxies = pd.read_csv(ifname, sep=',')
        cu_gra = cu_galaxies.Ra.values
        cu_gdec = cu_galaxies.Dec.values


    # Step 2
    # ------

    print 'Step 2 : Exclude coordinates that are too close to any supernova coordinate ...'

    # Criteria?
    # Unknown, but for now it should just lie
    # outside the range of the cutout extent, i.e. more than 101 px away.
    # 101 px in the SDSS frames correspond to about 10 ** -2 degrees.
    criteria_distance = .001  # [deg]

    # Count how many rows that are left at this step
    N_gx = cu_galaxies.shape[0]
    N_sn = snlist.shape[0]

    # How many times can I have an array of <N_rows> rows by
    N_rows = 10000
    N_gx_eee = N_gx // N_rows
    N_gx_rem = N_gx % N_rows

    dix = np.array([]).astype(bool)

    # Create repeatedly used vectors
    N_sn_ZEROS = np.zeros(N_sn)[None, :]

    RA_sn = np.zeros(N_rows)[:, None] + snlist.Ra.values[None, :]
    DEC_sn = np.zeros(N_rows)[:, None] + snlist.Dec.values[None, :]

    print 'Creating matrices that can calculate all distances simultaneously ...'

    # Loop for as many times as needed
    for n in range(N_gx_eee):

        beg = n * N_rows
        end = (n + 1) * N_rows

        # How far are we?
        print 'n = {: >4d}; {: >6d}; {: >6d};'.format(n, beg, end)

        # Create matrices
        # Broadcast shapes to get a N_gx-by-N_sn
        RA_gx = cu_gra[beg:end, None] + N_sn_ZEROS
        DEC_gx = cu_gdec[beg:end, None] + N_sn_ZEROS

        # print 'Calculating differences for each coordinate type ...'

        # Differences
        dRA = RA_gx - RA_sn
        dDEC = DEC_gx - DEC_sn

        # print 'Calculating the distances between every possible set of coordinates ...'

        # Distances from each coordinate to each supernova
        dS = np.sqrt(dRA ** 2 + dDEC ** 2)

        # print 'Creating boolean vector for each coordinate ...'

        # Are there any SNe too close for a given coordinate?
        # Check along the columns, i.e .return boolean vector of rows (galaxies)
        # that met the criteria of being far enough away that it should be outside
        # a cutout which also has a known SDSS supernova candidate within it.
        # distance indices
        dix = np.append(dix, np.any(dS > criteria_distance, axis=1))

    if N_gx_rem > 0:

        # Finally, the remaining less than <N_rows> coordinates
        beg = (n + 1) * N_rows
        end = beg + N_gx_rem

        # Create matrices

        RA_gx = cu_gra[beg:end, None] + N_sn_ZEROS
        DEC_gx = cu_gdec[beg:end, None] + N_sn_ZEROS

        RA_sn = np.zeros(N_gx_rem)[:, None] + snlist.Ra.values[None, :]
        DEC_sn = np.zeros(N_gx_rem)[:, None] + snlist.Dec.values[None, :]

        # print 'Calculating differences for each coordinate type ...'

        # Differences
        dRA = RA_gx - RA_sn
        dDEC = DEC_gx - DEC_sn

        # print 'Calculating the distances between every possible set of coordinates ...'

        # Distances from each coordinate to each supernova
        dS = np.sqrt(dRA ** 2 + dDEC ** 2)

        # print 'Creating boolean vector for each coordinate ...'

        # Append the booleans to my master index file
        dix = np.append(dix, np.any(dS > criteria_distance, axis=1))

    # Check
    print 'N_gx =', N_gx
    print 'dix.size =', dix.size

    print 'Reducing data for the next step ...'

    dcu_galaxies = cu_galaxies.iloc[dix]

    # Step 3:
    # -------
    print 'Step 3 : Exclude coordinates that are too close to other non-events ...'

    """
    Rationale
    ---------
    When I discovered that I had not checked that the distance between the
    galaxy coordinates themselves were not too close to essentially stay out
    of each other's cutouts, I realised that I never made sure that I could
    order them by observation date and then use only the first-observed object
    (or it could be the same, but it would make no difference for my algorithm
    to begin with), so that I could argue that I should keep the first one.
        As it is (2014-02-26), I do not have any means to use observation
    dates, since I do not have the luxury of starting over and begin with
    Step 0. Also, More importantly, I do not have time to spend on finding out
    how to obtain observation dates from the SDSS database.

    As this is basically just a proof of concept, where I just need enough
    coordinates that as far as the merged snlists are concerned do not have any
    known transient events happening in them that are classified as SNe. There
    does not seem to be any reason why I should not just choose arbitrarily
    between two coordinates whose cutouts (101x101) will overlap.

    Algorithm
    ---------
    1. I go through the list from the top.
    2. For each coordinate I calculate the distance to all other coordinates.
    3. I get an array of the same length, containing the distances.
    4. Entries with coordinates too close to be outside the given coordinate's
       cutout `view` will be removed from the list.
    5. The now potentially reduced list is used in the next step. Since we
       start from the top, the previous coordinate entries will all be there in
       the reduced list.
    6. Increase the entry index and repeat from step 3. until the end of the
       final entry is reached (no more remaining possibly too close coordinates
       left).

    """

    i = 0
    ddcu_galaxies = dcu_galaxies.copy()
    while i < ddcu_galaxies.shape[0]:
        if i and not i % 1000:
            print 'i = {: >5d}'.format(i)
        templist = ddcu_galaxies.copy()
        entry = templist[i:i + 1]
        dRa = entry.Ra.values[0] - templist.Ra.values
        dDec = entry.Dec.values[0] - templist.Dec.values
        dS = np.sqrt(dRa ** 2 + dDec ** 2)
        dix = (dS > criteria_distance)
        dix[i] = True
        ddcu_galaxies = templist.iloc[dix].copy()
        del templist
        # Yikes :O !!! This turned out to be important :D
        i += 1

    # print 'Final size of gxlist: {:,.0f}'.format(ddcu_galaxies.shape[0])

    # Step 4:
    # -------
    print 'Step 4 : Exclude coordinates that are covered by too few fields ...'

    """
    This was determined from the coordinates that were in the first tlist
    made before this step, where I discovered that some of the choden coordintes
    had as little as 1 field covering it. With fewer fields covering, the chance
    of getting a number of cutouts in a sequence that matches the average number of cutouts
    for the coordinates in snlist is smaller. There could be a bias here.
    """
    # Lowest frame count 69.0
    # Highest frame count 162.0
    MIN_NUMBER_OF_FIELDS = 69
    MAX_NUMBER_OF_FIELDS = 162
    # Check if enough covering fields
    # If not enough, exclude the coordinate.
    n_fields = np.array([])
    for i in range(ddcu_galaxies.shape[0]):
        print '{: >5d}'.format(i),
        Ra, Dec = ddcu_galaxies.iloc[i]
        N = get_covering_fields(np.array([Ra, Dec])[None, :]).shape[0]
        n_fields = np.append(n_fields, N)
        print '- Fields: {: >3d}'.format(N)

    fix = (
        (n_fields >= MIN_NUMBER_OF_FIELDS)
        &
        (n_fields <= MAX_NUMBER_OF_FIELDS)
    )

    fddcu_galaxies = ddcu_galaxies.iloc[fix]

    print 'Final size of gxlist: {:,.0f}'.format(fddcu_galaxies.shape[0])

    # Finalise
    # --------

    print 'Step finalise : save the resulting list to disk.'

    fddcu_galaxies.to_csv(env.files.get('gxlist'), index=False, header=True)


def plot_gxlist():
    """
    Generate figure showing the galaxy coordinates within the stripe
    plotted over the regions covered by the fields that are available.
    """
    pass


###############################################################################
def build_tlist():
    """
    Build the event/non-event data set and save it as a file.

    """
    import numpy as np
    import pandas as pd

    from ..parse import ra2deg, dec2deg

    if 1:
        snlist = pd.read_csv(env.files.get('snlist'), sep=';')
    else:
        snid_sortable = lambda SDSS_id: 'SN{:0>5d}'.format(int(SDSS_id[2:]))
        s2valornan = lambda s: s or np.nan
        conv = dict(SDSS_id=snid_sortable, Ra=ra2deg, Dec=dec2deg,
            redshift=s2valornan, Peak_MJD=s2valornan)
        lkw = dict(sep=';', converters=conv)
        snlist = pd.read_csv(env.files.get('snlist_1030'), **lkw)
    gxlist = pd.read_csv(env.files.get('gxlist'), sep=',')

    # print gxlist.info()
    # print gxlist.head(10)
    print snlist.info()
    print snlist.head(10)

    # How many needed in total
    N_sne = snlist.shape[0]
    N_gx = gxlist.shape[0]

    N_needed = np.round(N_sne, decimals=-2) * 2
    N_gx_c = N_needed - N_sne

    gx_ix = np.unique(np.random.randint(0, N_gx, size=N_gx_c))
    while gx_ix.size != N_gx_c:
        N_cur = gx_ix.size
        N_left = N_gx_c - N_cur
        gx_ix = np.unique(
            np.append(
                gx_ix,
                np.random.randint(0, N_gx, size=N_left)
            )
        )

    gx_chosen = gxlist.iloc[gx_ix]
    # print gx_chosen.info()
    # print gx_chosen.head(10)
    # raise SystemExit

    # Build data set
    ra = np.append(snlist.Ra.values, gx_chosen.Ra.values)
    dec = np.append(snlist.Dec.values, gx_chosen.Dec.values)
    is_sn = np.append(np.ones(N_sne), np.zeros(N_gx_c)).astype(bool)

    # Collect and shuffle the lines, so that I only need to split
    # the data set N-fold, when using the data.
    dataset = np.array([ra,
                        dec,
                        is_sn]).T
    # Do in-place shuffle
    """
    This is an in-place operation on a view of the original array.
    It does not create a new, shuffled array, so there's no need
    to transpose the result.
    REF: https://stackoverflow.com/questions/20546419/shuffle-columns-of-an-array-with-numpy
    """
    np.random.shuffle(dataset)

    if 0:
        coords = np.array([])
        # for i in range(dataset.shape[0]):
        #     coords = np.append(coords, '{:014.9f}_{:014.9f}'.format(
        #         dataset[i, 0], dataset[i, 1])
        #     )
        for i in range(snlist.shape[0]):
            coords = np.append(coords, '{:014.9f}_{:014.9f}'.format(
                snlist.Ra.values[i], snlist.Dec.values[i])
            )
        ucoords, indices = np.unique(coords, return_inverse=True)
        print ucoords.size, np.unique(snlist.SDSS_id.values).size, np.unique(snlist.Peak_MJD.values).size
        raise SystemExit

    # tlist = pd.DataFrame(data=dict(Ra=ra, Dec=dec, is_sn=is_sn))
    tlist = pd.DataFrame(
        data=dataset,
        columns=['Ra', 'Dec', 'is_sn']
    )
    print tlist.info()
    print tlist.head(10)
    tlist.to_csv(env.files.get('tlist'), index=False, header=True)


###############################################################################
def build_tlist_sample(N=5):
    """
    Build *SAMPLE* event data set and save it as a file.

    """
    import numpy as np
    import pandas as pd

    snlist = pd.read_csv(env.files.get('snlist'), sep=';')

    ra = snlist.Ra.values[:N]
    dec = snlist.Dec.values[:N]
    is_sn = np.ones(N).astype(bool)

    dataset = np.array([ra,
                        dec,
                        is_sn]).T

    tlist = pd.DataFrame(
        data=dataset,
        columns=['Ra', 'Dec', 'is_sn']
    )
    print tlist.info()
    print tlist.head(N)
    tlist.to_csv(env.files.get('tlist'), index=False, header=True)


###############################################################################
def check_snlist():
    """
    Check for duplicates in snlist_1030

    """
    import numpy as np
    import pandas as pd

    from ..parse import ra2deg, dec2deg

    # Converters
    snid_sortable = lambda SDSS_id: 'SN{:0>5d}'.format(int(SDSS_id[2:]))
    s2valornan = lambda s: s or np.nan

    if 0:
        ifname = env.files.get('snlist_1030')
        conv = dict(SDSS_id=snid_sortable, Ra=ra2deg, Dec=dec2deg,
            redshift=s2valornan, Peak_MJD=s2valornan)

    else:
        ifname = env.files.get('snlist')
        conv = dict(SDSS_id=snid_sortable, redshift=s2valornan,
                    Peak_MJD=s2valornan)

    lkw = dict(sep=';', converters=conv)
    snlist = pd.read_csv(ifname, **lkw)

    # Check for duplicate coordinate pairs
    coords = np.array([])
    for i in range(snlist.shape[0]):
        coords = np.append(coords, '{:014.9f}_{:014.9f}'.format(
            snlist.Ra.values[i], snlist.Dec.values[i])
        )
    ucoords, indices = np.unique(coords, return_inverse=True)
    print 'Number of list entries:     {: >4d}'.format(snlist.shape[0])
    print 'Number of unique entries:   {: >4d}'.format(ucoords.size)
    # print 'Number of unique entry IDs: {: >4d}'.format(np.unique(snlist.SDSS_id.values).size)

    duplicates = []
    for ix in np.unique(indices):
        if (indices == ix).sum() > 1:
            duplicates.append(ix)

    if duplicates:

        coord_indices = []
        for ix in duplicates:
            print ''
            for i, uc in enumerate(ucoords[indices[indices == ix]]):
                if i == 0:
                    coord_indices.append((coords == uc).nonzero()[0])
                print uc

        print '\nIndices of the original list:', duplicates
        print coord_indices
        print ''

        print 'Entries from snlist:'
        for cices in coord_indices:
            print
            print snlist.iloc[cices]

    else:
        print 'No duplicates found ...'

###############################################################################
def check_tlist():
    """
    For each entry in `tlist.csv`, find out how many fields cover it.

    """
    from .cutouts import get_covering_fields

    import numpy as np
    import pandas as pd

    ifname = env.files.get('tlist')
    tlist = pd.read_csv(ifname)

    sn_fields = []
    gx_fields = []
    for i in range(tlist.shape[0]):
        print '{: >4d} -'.format(i),
        Ra, Dec, is_sn = tlist.iloc[i]
        n_fields = get_covering_fields(np.array([Ra, Dec])[None, :]).shape[0]
        if is_sn:
            sn_fields.append(n_fields)
            print 'SN',
        else:
            gx_fields.append(n_fields)
            print 'GX',
        print '- Fields: {: >3d}'.format(n_fields)

    for data, name in zip([sn_fields, gx_fields], ['SN', 'GX']):
        ofname = os.path.join(
            env.paths.get('data'),
            'nfieldrecords_{}.csv'.format(name)
        )
        with open(ofname, 'w+') as fsock:
            fsock.write('\n'.join(np.array(data).astype(str)))

    if 0:
        import matplotlib as mpl
        mpl.use('pdf')
        import matplotlib.pyplot as plt

        from mplconf import mplrc
        from mplconf import rmath

        mplrc('publish_digital')
        fig, ax = plt.subplots(figsize=(12.5, 4))
        ax.hist(sn_fields, bins=100)
        ax.hist(gx_fields, bins=100)
        ax.set_xlabel(rmath('Number of fields for a coordinate'))
        ax.set_ylabel(rmath('Counts'))
        fig.tight_layout()
        ofname_fig = os.path.join(
            env.paths.get('data'),
            'tlist_nfieldrecords.pdf'
        )
        plt.savefig(ofname_fig)

###############################################################################
def get_fields(skip_if_exists=True):
    """For each coordinate in the SNe file, get the frames
    from the Field table that cover that coordinate.

    Saves each result in a seperate file named <SDSS_ID>.csv .

    """

    # Clean up local cas directory
    # field_clean_local_dir()
    # Do it manually from the main program instead
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
    -- Stripe82..Field
    Field

WHERE
    raMin < {obj_ra} AND {obj_ra} < raMax
  AND
    decMin < {obj_dec} AND {obj_dec} < decMax
  AND
    (raMax - raMin) > {dra_min}
  AND
    (raMax - raMin) < {dra_max}

ORDER BY
    run ASC,
    rerun ASC,
    camcol ASC,
    field ASC,
    raMin ASC,
    decMin ASC
"""

    opath = os.path.join(_path_cas, 'fields')
    if not os.path.exists(opath):
        os.makedirs(opath)

    df = load_SNe_candidate_list()
    SNe_len = df.count().max()

    print 'Beginning search through all confirmed SNe ...'

    time_total_beg = time.time()

    for i, (ra, dec, SDSS_id) in enumerate(zip(df.Ra, df.Dec, df.SDSS_id)):
        sql_fields = sql_fields_fstr.format(obj_ra=ra, obj_dec=dec, dra_min=.1, dra_max=1.)
        # get_field_result(sql_fields, SDSS_id)
        # Define output structure
        ofname = os.path.join(opath, '{}.csv'.format(SDSS_id))

        # Don't bother requesting anything if the file already exists
        if os.path.isfile(ofname) and skip_if_exists:
            return

        # Update progress output
        s = 'Downloading: SN {: >4d} out of {: >4d}, {} ...\r'.format(
            (i + 1), SNe_len, format_HHMMSS_diff(time_total_beg, time.time())
        )
        sys.stdout.write(s)
        sys.stdout.flush()

        # Request the data
        fields = query(sql_fields)
        if 'error' in fields.lower():
            sys.exit('ERROR: {}: CAS says something is wrong ...'.format(
                sys._getframe().f_code.co_name)
            )

        # And save it
        with open(ofname, 'w+') as fsock:
            fsock.write(fields)

    sys.stdout.write('\n')
    sys.stdout.flush()

    time_total_end = time.time()
    time_total = format_HHMMSS_diff(time_total_beg, time_total_end)

    print 'Done downloading field catalogue ... in {}'.format(time_total)
    # print 'Downloaded field catalogues for {} SNe'.format(i+1)


def get_field_result(sql_fields, SDSS_id):
    """Get query results from a *single* request for fields."""

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
    ofname = env.files.get('fields')
    tfname = os.path.join(_path_cas, 'fields.tmp')

    # Clean up first, since the file is only appended to in the following
    if os.path.isfile(ofname):
        os.remove(ofname)

    commands = [

        # Build one big file with all the results
        'cat {iglob} >> {t}'.format(iglob=iglob, t=tfname),

        # Sort and remove duplicates
        'cat {t} | sort | uniq > {o}'.format(t=tfname, o=ofname),

        # Remove the temporary file
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


def filter_invalid_from_unique_field_list(dra_min=.1, dra_max=1.):
    """
    Remove invalid entries from the unique field
    list created by *create_unique_field_list()*.

    Parameters
    ----------
    dra_min : float
        the minimum allowable angle separating the RA start and
        end coordinate for a given field in the unique field list.

    dra_max : float
        the maximum allowable angle separating the RA start and
        end coordinate for a given field in the unique field list.

    Side effects
    ------------
    <del>Creates a backup of the original field list.</del>
    Saves the results back into the original destination.

    """
    import pandas as pd
    import shutil

    fname = env.files.get('fields')
    # fn_valid = os.path.splitext(fname)[0] + '_valid.csv'
    fn_invalid = os.path.splitext(fname)[0] + '_invalid.csv'

    shutil.copyfile(fname, '{}.orig'.format(fname))

    df = pd.read_csv(fname, sep=',')

    """
    2014-07-28
    ----------
    There is an error in this procedure, but the result is the intended:
    to remove fields for which the physical extent---as given by the
    coordinates raMax, raMin, decMax, decMin---gives too small or even
    negative side lengths. This problem is only observed in the RA
    coordinates.

    The are five observed cases for where the RA coordinate extents land.

    Case 1: raMax \in [  0;  60] and raMin \in [  0;  60]
    Case 2: raMax \in [  0;  60] and raMin \in [300; 360]
    Case 3: raMax \in [300; 360] and raMin \in [300; 360]
    Case 4: raMax \in [  0;  60] and raMin \in [300; 360]
    Case 5: raMax > 360          and raMin \in [  0;  60]

    Case 4 should not occur, since this means that the field is obtained
    end-first and beginning-last. These fields are considered invalid.

    Case 5 is OK, if subtracting 360 deg from raMax, the value is larger
    than raMin. Otherwise, the coordinate difference produces a negative
    side length again.
      It turns out that this is the case, so these fields are also invalid.

    The procedure below removes all fields for which the coordinate difference
    raMax - raMin < dra_min = .1 (default). Since this also removes all Case-4
    and Case-5 records above, the result is what is intended; but with wrong
    assumptions.

    Since this is the way that my data were processed, I leave it at that,
    if I need to reproduce the field list again.

    """

    dra = df.raMax - df.raMin

    dra_too_small_ix = (dra < dra_min)
    dra_too_large_ix = (dra > dra_max)

    df_invalid = df.loc[dra_too_small_ix | dra_too_large_ix]
    df_valid = df.loc[(~dra_too_small_ix) & (~dra_too_large_ix)]

    df_valid.to_csv(fname, index=False, header=True)
    df_invalid.to_csv(fn_invalid, index=False, header=True)


def count_field_records():
    """
    Count the number for field records obtained for
    each SN, and save those numbers for later plotting.

    """

    import glob
    import pandas as pd

    iglob = os.path.join(_path_cas, 'fields', '*.csv')
    filenames = sorted(glob.glob(iglob))

    df_fields = pd.read_csv(env.files.get('fields'), sep=',')

    counts = []
    beg = time.time()
    for i, ifname in enumerate(filenames):
        df_results = pd.read_csv(ifname, sep=',')
        count = 0
        for j in range(df_results.shape[0]):
            count += (df_fields['fieldID'] == df_results.iloc[j]['fieldID']).sum()
        step = time.time()
        dt_str = format_HHMMSS_diff(beg, step)
        print '{: >4d}, {: >3d}, {}'.format(i, count, dt_str)
        counts.append(str(count))

    # print len(filenames), len(counts)

    with open(env.files.get('nrecords'), 'w+') as fsock:
        fsock.write('\n'.join(counts))


def count_field_records_by_quality():
    """
    Count the number for field records obtained for
    each SN, and save those numbers for later plotting.

    """

    import glob
    import pandas as pd

    iglob = os.path.join(_path_cas, 'fields', '*.csv')
    filenames = sorted(glob.glob(iglob))

    df_fields = pd.read_csv(env.files.get('fields'), sep=',')
    print 'Number of fields:', df_fields.shape[0]

    # FieldQuality Data values
    # name          value   description
    # BAD           1       Not acceptable for the survey
    # ACCEPTABLE    2       Barely acceptable for the survey
    # GOOD          3       Fully acceptable -- no desire for better data
    # MISSING       4       No objects in the field, because data is missing.
    #                       We accept the field into the survey as a HOLE
    # HOLE          5       Data in this field is not acceptable, but we will
    #                       put the field into the survey as a HOLE, meaning
    #                       none of the objects in the field are part of the
    #                       survey.
    # See: http://cas.sdss.org/stripe82/en/help/browser/enum.asp?n=FieldQuality
    qualities = range(1, 4)
    qices = [(df_fields.quality.values == q) for q in qualities]
    fdict = {q: df_fields.iloc[qix] for (q, qix) in zip(qualities, qices)}
    cdict = {q: [] for q in qualities}

    for i, qix in enumerate(qices):
        print 'Number of fields with quality {:d}: {: >3d}'.format(
            i + 1, qix.sum())
    print 'For qualities 4 and 5, there were not present in my filtered dataset'

    beg = time.time()
    for i, ifname in enumerate(filenames):
        df_results = pd.read_csv(ifname, sep=',')
        counts = [0] * len(qualities)
        for j in range(df_results.shape[0]):
            for k, q in enumerate(qualities):
                ices = (fdict[q]['fieldID'] == df_results.iloc[j]['fieldID'])
                counts[k] += ices.sum()
        step = time.time()
        Dt = format_HHMMSS_diff(beg, step)
        print '{: >4d}, {}: {: >3d}, {: >3d}, {: >3d}'.format(i, Dt, *counts)
        for k, q in enumerate(qualities):
            cdict[q].append(counts[k])

    list_of_lists = [cdict[q] for q in qualities]
    with open(env.files.get('nrecords_q'), 'w+') as fsock:
        for row in zip(*list_of_lists):
            fsock.write('{},{},{}\n'.format(*row))


def query(sql_raw):
    """Sends SQL query to the CAS and returns the raw text of the response."""

    # form_data_receive_only_url = 'http://cas.sdss.org/astro/en/tools/search/x_sql.asp'
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

