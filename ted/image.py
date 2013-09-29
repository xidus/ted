#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sun 29 Sep 2013
#   Initial build.
#

import numpy as np


# From: http://prancer.physics.louisville.edu/astrowiki/index.php/Image_processing_with_Python_and_SciPy
# Define a function for making a linear gray scale
def lingray(x, a=None, b=None):
    """
    Auxiliary function that specifies the linear gray scale.
    a and b are the cutoffs : if not specified, min and max are used
    """
    if a == None:
            a = np.min(x)

    if b == None:
            b = np.max(x)

    return 255.0 * (x - float(a)) / (b - a)

