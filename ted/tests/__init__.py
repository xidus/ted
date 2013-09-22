#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Tue 13 Aug 2013
#   Initial build.
#

"""
Test module.

"""

import unittest


class EnvironmentSetupTest(unittest.TestCase):
    """"""

    def test_path(self):
        assert 'path' in locals(), 'Environment settings not found.'


if __name__ == '__main__':
    unittest.main()
