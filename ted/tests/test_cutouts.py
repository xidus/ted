#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Tue 13 Aug 2013
#   Initial build.
#

"""
Test ted.sdss.cutouts
"""

def test_minlocs2Dstack():

    import numpy as np

    from ted.sdss.cutouts import minloc2Dstack

    # Local minima is in (1, 1)
    mat0 = [
        [2, 1, 2, 0],
        [2, 0, 1, 1],
        [2, 1, 1, 2],
        [2, 1, 1, 2],
    ]
    # Local minima is in (2, 1)
    mat1 = [
        [5, 4, 5, 3],
        [5, 4, 4, 5],
        [5, 3, 4, 4],
        [5, 4, 4, 5],
    ]
    # Local minima is in (2, 2)
    mat2 = [
        [8, 7, 8, 6],
        [8, 7, 7, 8],
        [8, 7, 6, 7],
        [8, 7, 7, 8],
    ]

    # Stack arrays depth wise
    cube = np.dstack((mat0, mat1, mat2))

    # Define answers
    answers = [
        (np.array([1]), np.array([1])),
        (np.array([2]), np.array([1])),
        (np.array([2]), np.array([2])),
    ]

    # Get the results
    minima = minloc2Dstack(cube)

    # Compare with expected
    for i in range(cube.shape[2]):
        assert answers[i] == minima[:, :, i].nonzero()


def test_threshold2Dstack():

    import numpy as np

    from ted.sdss.cutouts import threshold2Dstack

    """stack"""

    # Local minima is in (1, 1)
    s0 = [
        [100,  50, 100],
        [100,   0,  50],
        [100,  50,  50],
    ]
    # Local minima is in (2, 1)
    s1 = [
        [200, 100, 200],
        [200,  50, 100],
        [200, 100, 100],
    ]
    # Local minima is in (2, 2)
    s2 = [
        [300, 200, 300],
        [300, 100, 200],
        [300, 200, 200],
    ]

    # Stack arrays depth wise
    stack = np.dstack((s0, s1, s2))

    """stack_locs"""

    # Local minima is in (1, 1)
    sl0 = [
        [ True, False, False],
        [False,  True, False],
        [False, False,  True],
    ]
    # Local minima is in (2, 1)
    sl1 = [
        [ True, False, False],
        [False,  True, False],
        [False,  True, False],
    ]
    # Local minima is in (2, 2)
    sl2 = [
        [False, False, False],
        [False,  True,  True],
        [False, False, False],
    ]

    # Stack arrays depth wise
    stack_locs = np.dstack((sl0, sl1, sl2))

    """answers"""

    # Local minima is in (1, 1)
    a0 = [
        [ True, False, False],
        [False, False, False],
        [False, False, False],
    ]
    # Local minima is in (2, 1)
    a1 = [
        [ True, False, False],
        [False, False, False],
        [False, False, False],
    ]
    # Local minima is in (2, 2)
    a2 = [
        [False, False, False],
        [False, False,  True],
        [False, False, False],
    ]

    # Define answers
    answer = np.dstack((a0, a1, a2))

    # Get the results
    result = threshold2Dstack(stack, stack_locs, tau=.5)

    # The test
    assert np.all(answer == result)

if __name__ == '__main__':
    test_threshold2Dstack()
