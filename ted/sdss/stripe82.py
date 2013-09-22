#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Tue 17 Sep 2013
#   Initial build.
#

# Stripe 82 coverage

# Stripe width defined along the parallels (Right Ascension)
ra_min = -60.
ra_max = 60.
stripe_width = ra_max - ra_min

# Stripe height defined along the meridians (Declination)
stripe_height = 2.5
dec_min = -stripe_height / 2.
dec_max = dec_min + stripe_height

# Define point pairs for the corners of the stripe extend
# Vertices define a rectangle moving from the lower right of the plot,
# to the left, then up, then right, and back to the lower right.
rect_coords_ra = [ra_min, ra_max, ra_max, ra_min, ra_min]
rect_coords_dec = [dec_min, dec_min, dec_max, dec_max, dec_min]

# Width to height
w2h = stripe_width / stripe_height
