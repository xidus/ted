#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Tue 17 Sep 2013
#   Initial build.
#

"""
Parsers and converters for astronomical and geospatial conversions.

Content
-------
dec2deg : Declination in sexagesimal system to floating-point degrees
ra2deg : Right Ascension in sexagesimal system to floating-point degrees
deg2ra : Floating-point degrees to Right Ascension in sexagesimal system

sex2deg : Degrees in sexagesimal system to floating-point degrees
deg2sex : Floating-point degrees to degrees in sexagesimal system

deg2lon : Floating-point degrees to longitude in sexagesimal system

"""


def dec2deg(dec):
    """
    Converts Declination to degrees in decimal notation.

    References
    ----------
    [1]: https://en.wikipedia.org/wiki/Celestial_coordinates#Notes_on_conversion

    """

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


def ra2deg(ra):
    """
    Converts Right Ascension to degrees in decimal notation.

    References
    ----------
    [1]: https://en.wikipedia.org/wiki/Right_ascension

    """

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


def deg2ra(deg):
    """
    Converts degrees to Right Ascension.

    Parameters
    ----------
    deg : float
        the degree measure which is to be converted.

    Returns
    -------
    ra : string, i.e hh:mm:ss.sss
        the hour measure in sexagesimal system

    """

    print deg

    if deg >= 360:
        deg %= 360

    elif deg < 0:
        # Parentheses included for clarity.
        deg = (deg % -60) + 360

    print deg

    # Now degrees are within the desired range.

    # Get integral number of hours and the remainder
    # that is to be resolved into minutes and seconds.

    # There are 360 deg / 24 hours = 15 degrees per hour
    hh = deg // 15
    hh_remainder = deg % 15  # [deg]

    # There are 360 deg / (24 * 60 minutes) = 1/4 degrees per minute
    minutes = hh_remainder * 4.  # [deg] * [min / deg] = [min]
    mm = minutes // 1
    mm_remainder = minutes % 1

    # Seconds
    ss = mm_remainder * 60

    # Create the output string
    fstr = '{hh:.0f}:{mm:0>2.0f}:{ss:0>6.3f}'
    return fstr.format(hh=hh, mm=mm, ss=ss)


def sex2deg(sex):
    """
    Converts sexagesimal to degrees.

    Parameters
    ----------
    sex : string, i.e deg:mm:ss.sss
        the degree measure in sexagesimal system

    Returns
    -------
    deg : float
        the degree measure which is to be converted.

    Note
    ----
    Does not take number of significant figures into account.

    """

    sgn = -1 if '-' == sex.strip()[0] else +1

    deg, minutes, seconds = [float(s) for s in sex.split(':')]
    return deg + sgn * (minutes / 60. + seconds / 3600.)


def deg2sex(deg, limits=None):
    """
    Converts degrees to sexagesimal system.

    Parameters
    ----------
    deg : float
        the degree measure which is to be converted.
    limits : [degree_min, degree_max]
        TO BE IMPLEMENTED
        Choose the 360 degree interval, e.g. [-180; +180[ .

    Returns
    -------
    sex : string, i.e deg:mm:ss.sss
        the degree measure in sexagesimal system

            0 <= deg < 360
            0 <= mm < 60
            0 <= ss.sss < 60

    Note
    ----
    Does not take number of significant figures into account.

    """

    # If deg is not (0 <= deg < 360), convert it to this range

    """
    There are different ways of thinking about remainders when you deal
    with negative numbers, and he is probably confusing two of them. The
    mod function is defined as the amount by which a number exceeds the
    largest integer multiple of the divisor that is not greater than that
    number. In this case, -340 lies between -360 and -300, so -360 is the
    greatest multiple LESS than -340; we subtract 60 * -6 = -360 from -340
    and get 20:

    -420 -360 -300 -240 -180 -120  -60   0   60   120  180  240  300  360
    --+----+----+----+----+----+----+----+----+----+----+----+----+----+--
           | |                                                    |  |
       -360| |-340                                             300|  |340
           |=|                                                    |==|
            20                                                     40

    Working with a positive number like 340, the multiple we subtract is
    smaller in absolute value, giving us 40; but with negative numbers, we
    subtract a number with a LARGER absolute value, so that the mod
    function returns a positive value. This is not always what people
    expect, but it is consistent.

    -- [Mod Function and Negative Numbers][1]

    References
    ----------
    [1]: http://mathforum.org/library/drmath/view/52343.html

    """

    if deg >= 360:
        deg %= 360

    elif deg < 0:
        # Parentheses included for clarity.
        deg = (deg % -60) + 360

    # Now degrees are within the desired range.

    # Get integral number of degrees and the remainder
    # that is to be resolved into minutes and seconds.
    dd = deg // 1
    dd_remainder = deg % 1

    # Minutes
    minutes = dd_remainder * 60
    mm = minutes // 1
    mm_remainder = minutes % 1

    # Seconds
    ss = mm_remainder * 60

    # Create the output string
    fstr = '{dd:.0f}:{mm:0>2.0f}:{ss:0>6.3f}'
    return fstr.format(dd=dd, mm=mm, ss=ss)


def deg2lon(deg):
    """
    Converts degrees to corresponding longitude on the globe.

    Parameters
    ----------
    deg : float
        the degree measure which is to be converted.

    Returns
    -------
    lon : string, i.e deg:mm:ss.sss
        the degree measure in lonagesimal format

            -180 <= deg < +180
            0 <= mm < 60
            0 <= ss.sss < 60

    Note
    ----
    Does not take number of significant figures into account.

    """

    if not (-180 <= deg < 180):
        print('Degrees not within allowed range')
        raise SystemExit

    sgn = -1 if deg < 0 else +1

    # Make sure degrees are within desired range
    # if deg < 0, the modulus function has to use -180.
    # deg %= sgn * 180
    # deg = sgn * (deg % (sgn * 180))

    # Get the integral number of degrees
    # and the remainder used to get the minutes and seconds.

    # We want to keep the integral degrees negative
    dd = sgn * (deg // sgn)

    # Whereas for the rest of the calculations we only,
    # need the positive numerical value of the remainder
    dd_remainder = deg % sgn

    # Convert remaining fractional degrees to minutes
    minutes = dd_remainder * 60
    mm = minutes // sgn
    mm_remainder = minutes % 1

    # Seconds
    ss = mm_remainder * 60

    # Create the output string
    fstr = '{dd:.0f}:{mm:0>2.0f}:{ss:0>6.3f}'
    return fstr.format(dd=dd, mm=mm, ss=ss)

