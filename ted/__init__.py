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

    # Class attributes
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

        # This is the important part
        for env_type in ('files', 'paths'):
            setattr(self, env_type, {})
            for env_type_key, env_type_list in doc.get(env_type).iteritems():
                getattr(self, env_type)[env_type_key] = os.path.join(*env_type_list)

        # Add the formatstrings...
        self.formatstrings = doc.get('formatstrings')

        # Add the proxies that requests will use
        self.proxies = doc.get('proxies')

# Initialise
env = Environment()

if __name__ == '__main__':
    from pprint import pprint
    pprint(env.paths)
    print ''
    pprint(env.files)
