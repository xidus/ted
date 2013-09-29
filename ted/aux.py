#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sun 29 Sep 2013
#   Initial build.
#

import os

from . import Environment

env = Environment()


def iscoadded(run):
    return run in (106, 206)


def frame_path(df_entry, frame_fmt='fpC', which='local'):
    """Returns the path of the given image."""
    keys = ('run', 'rerun', 'camcol', 'field')
    filters = 'ugriz'
    fstr_frame = env.fstrings.get('fpC', None)

    if which == 'local':
        base_prefix = os.path.join(env.paths['data'], 'das', frame_fmt)

    else:
        base_prefix = 'http://das.sdss.org'

    string_params = {key: df_entry[key].values[0] for key in keys}

    if iscoadded(string_params['run']):
        # For now, pass over the coadded frames (since these do not have an observation date).
        # return None
        # or return anyway and parse somewhere else
        string_params['run'] = (string_params['run'] - 6) * 1000 + 6

    string_params['filt'] = filters[2]  # r

    # for key, val in string_params.iteritems():
    #     print key, val, type(val)

    fname = fstr_frame.format(**string_params)
    return os.path.join(base_prefix, fname)


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


