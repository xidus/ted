#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sun 29 Sep 2013
#   Initial build.
#

def format_HHMMSS_diff(beg, end):
        diff = end - beg
        HH = diff // 60 // 60
        MM = diff // 60 - HH * 60
        SS = diff % 60
        return '{:.0f}:{:0>2.0f}:{:0>6.3f}'.format(HH, MM, SS)



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


# def seconds_since_1858_11_17_to_datetime(seconds=None):
def tai2date(seconds=None):
    """
    Converts *International Atomic Time (TAI)* to Common-Era (CE) Python datetime object.

    TAI: Number of econds since 1858-11-17.

    Parameters
    ----------
    seconds : positive int or float, required
        The number of seconds since 17 November 1858

    Returns
    -------
    Python datetime.datetime

    References
    ----------
    Survey Interface File for Corrected Frame
        http://www.sdss.org/DR7/dm/flatFiles/fpC.html

    """
    import datetime as dt

    if seconds is None:
        raise ValueError('Please provide number of seconds')

    try:
        seconds = float(seconds)

    except Exception:
        raise ValueError('Seconds must be of type int or float')

    return dt.datetime(1858, 11, 17) + dt.timedelta(seconds=seconds)

