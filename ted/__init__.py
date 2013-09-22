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
        # for file_key, file_info in doc['files'].iteritems():
        #     fname = file_info.get('fname')
        #     if fname is None:
        #         continue
        #     base = file_info.get('base', '')
        #     subdirs = file_info.get('subdirs', [''])
        #     path_parts = [base] + subdirs + [fname]
        #     self.files[file_key] = os.path.join(*path_parts)
        for file_key, path_list in doc['files'].iteritems():
            self.files[file_key] = os.path.join(*path_list)

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
