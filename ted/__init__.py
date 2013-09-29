#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Tue 13 Aug 2013
#   Initial build.
#

__author__ = 'Joachim Mortensen <mojo@fys.ku.dk>'
__version__ = '0.1'

# __all__ = ['pipe', 'sdss']


import os

import yaml

# from . import (
#     sdss,
#     pipe,
# )


class TEDError(Exception): pass
class EnvironmentLoadError(TEDError): pass
class NoneExistingFileError(TEDError): pass

# _pkg_home_dir = __path__[0] if __path__ else os.path.dirname(os.path.realpath(__file__))

# __path__ is not automatically defined, when simply running this file.
_pkg_home_dir = os.path.dirname(os.path.realpath(__file__))


class Environment(object):
    """Looks for an environment-configuration file and loads the paths."""

    paths = None
    files = None

    def __init__(self, ifname=None):

        self._env_file = os.path.join(_pkg_home_dir, 'env.yaml')

        if ifname is not None:

            if os.path.isfile(ifname):
                self._env_file = ifname

            else:
                raise NoneExistingFileError('File does not exists ...')

        try:
            with open(self._env_file, 'r') as fsock:
                doc = yaml.load(fsock.read())

        except:
            raise EnvironmentLoadError('Path-configuration file could not load ...')

        self.paths = doc.get('paths', None)
        self.fstrings = doc.get('fstrings', None)

        self.files = {}
        for file_key, path_list in doc.get('files').iteritems():
            self.files[file_key] = os.path.join(*path_list)


if __name__ == '__main__':
    from pprint import pprint

    # Initialise
    env = Environment()

    pprint(env.paths)
    print ''
    pprint(env.files)
