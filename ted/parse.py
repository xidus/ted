#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Tue 17 Sep 2013
#   Initial build.
#


def dec2deg(dec):
    """Converts Declination to decimal notation."""

    # It is not enough to just check the sign of the first part.
    # if this is -0, then np.sign() returns a zero,
    # but this is not what we want when either of the minutes-
    # and seconds-part of the coordinate are non-zero.

    if '-' == dec.strip()[0]:
        sgn = -1
    else:
        sgn = +1

    deg, minutes, seconds = [float(s) for s in dec.split(':')]
    return deg + sgn * (minutes / 60. + seconds / 3600.)


# ra_counts = 0

def ra2deg(ra):
    """Converts Rect Ascension decimal notation."""

    # global ra_counts
    # ra_counts += 1
    # print ra_counts,
    # print len(ra.strip()), ra, ra.strip()

    # It is not enough to just check the sign of the first part.
    # if this is -0, then np.sign() returns a zero,
    # but this is not what we want when either of the minutes-
    # and seconds-part of the coordinate are non-zero.

    if '-' == ra.strip()[0]:
        sgn = -1
    else:
        sgn = +1

    hours, minutes, seconds = [float(s) for s in ra.split(':')]
    return hours * 15. + sgn * (minutes / 4. + seconds / 240.)


def mjd2date(mjd):
    """Converts *Modified* Julian Date (MJD) to Common-Era (CE) Python datetime object."""
    import datetime as dt
    # Smithsonian Date
    mjd_epoch = dt.datetime(1858, 11, 17, 0, 0, 0)
    if not isinstance(mjd, float):
        # Assume mjd is iterable
        return [mjd_epoch + dt.timedelta(days=nfracdays) for nfracdays in mjd]
    else:
        # Return a single datetime object
        return mjd_epoch + dt.timedelta(days=mjd)

