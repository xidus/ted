#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Tue 13 Aug 2013
#   Initial build.
#

__author__ = 'Joachim Mortensen <mojo@fys.ku.dk>'
__version__ = '0.1'

__all__ = []  # ['pipe', 'sdss']


# import sys
import os

import yaml


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

    def __init__(self, *args, **kwargs):

        self._paths_file = os.path.join(_pkg_home_dir, 'env.yaml')
        # print os.path.isfile(self._paths_file)

        try:
            with open(self._paths_file, 'r') as fsock:
                doc = yaml.load(fsock.read())

        except:
            raise EnvironmentLoadError('Path-configuration file could not load ...')

        self.paths = doc['paths']
        self.files = {}
        for file_key in doc['files']:
            file_info = doc['files'][file_key]
            base_key = file_info.get('base')
            base = self.paths[base_key] if base_key is not None else ''
            subdirs = file_info['subdirs'] if file_info.get('subdirs') is not None else ['']
            fname = file_info['fname']
            self.files[file_key] = os.path.join(*([base] + subdirs + [fname]))

            # This does not work for file paths that are defined for files that have yet to be created.
            # # Check if file exists
            # if not os.path.isfile(self.files[file_key]):
            #     raise NoneExistingFileError('File does not exists ({})...'.format(self.files[file_key]))

        # Add the attributes
        # self.dict2attr(doc)

    # def dict2attr(self):
    #     pass

# Initialise
env = Environment()


if __name__ == '__main__':
    from pprint import pprint
    pprint(env.paths)
    print ''
    pprint(env.files)
