#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Tue 13 Aug 2013
#   Initial build.
#

"""
The Transient-Event Detector (TED)

"""

__author__ = 'Joachim Mortensen <mojo@fys.ku.dk>'
__version__ = '0.1'

# __all__ = ['pipe', 'sdss']

"""
Notes to self
-------------

When `ted/sdss/__init__.py` contains

```python
from .. import env
```

Then loading both modules

In [1]: import ted
In [2]: import ted.sdss

Does not recreate the variables within the scope of `ted`.
They are simply referred to by variables within the scope of whatever
module which imported from `..` and into its own namespace.

In [3]: ted.sdss.env
Out[3]: <ted.Environment at 0xa00be2c>

In [4]: ted.env
Out[4]: <ted.Environment at 0xa00be2c>

---

"""

import os

import yaml

# from . import (
#     sdss,
#     pipe,
# )


# _pkg_home_dir = __path__[0] if __path__ else os.path.dirname(os.path.realpath(__file__))
# __path__ is not automatically defined, when simply running this file.
_pkg_home_dir = os.path.dirname(os.path.realpath(__file__))


class TEDError(Exception): pass
class EnvironmentLoadError(TEDError): pass
class NoneExistingFileError(TEDError): pass

class Namespace(object): pass

# class TEDRelPath(yaml.YAMLObject, list):

#     yaml_tag = u'!TEDRelPath'

#     def __init__(self, rel_path):

#         self.rel_path = rel_path
#         self.abs_path = os.path.join(_pkg_home_dir, self.rel_path)

#     def __repr__(self):
#         return [self.abs_path]

#     def __str__(self):
#         return [self.abs_path]


# __ted = TEDRelPath('sql')
# print __ted
# print '\n'.join(*__ted)
# print repr(env.paths.get('sql'))


def msg(istr, width=78, char='-'):
    if len(char) > 1: char = char[:1]
    width_left = max(0, width - len(istr))
    padding = char * max(0, width_left // 2 - 1)
    ostr = '{padding} {istr} {padding}'.format(padding=padding, istr=istr)
    while len(ostr) != width:
        ostr += char[:1]
    print '\n' + ostr + '\n'


class Environment(object):
    """Looks for an environment-configuration file and loads the paths."""

    _ifname_base = os.path.join(_pkg_home_dir, 'base.yaml')
    _ifname_base_default = os.path.join(_pkg_home_dir, 'base.default.yaml')
    _ifname_setup = os.path.join(_pkg_home_dir, 'setup.yaml')
    _ifname_env = os.path.join(_pkg_home_dir, 'env.yaml')

    def __init__(self):

        # self._check_source_files()
        self.update_base()
        self.load()

    def load(self):

        try:
            with open(self._ifname_env, 'r') as fsock:
                doc = yaml.load(fsock.read())

        except:
            err = 'Path-configuration file could not load ...'
            raise EnvironmentLoadError(err)

        # This is the important part
        for env_type in ('files', 'paths'):
            # Re-initialise
            setattr(self, env_type, {})
            for env_type_key, env_type_list in doc.get(env_type).items():
                getattr(self, env_type)[env_type_key] = os.path.join(
                    *env_type_list)

        # Add file params.yaml in the same directory as the environment files.
        getattr(self, 'files')['params'] = os.path.join(
            _pkg_home_dir, 'parameters.yaml')

        # Add the sql-script directory
        getattr(self, 'paths')['sql'] = os.path.join(_pkg_home_dir, 'sql')

        # Add the formatstrings...
        self.formatstrings = doc.get('formatstrings')

        # Add the proxies that requests will use
        self.proxies = doc.get('proxies')

    def _check_source_files(self):

        if not os.path.isfile(self._ifname_env):
            self.update_base()

    def update_base(self):

        if not os.path.isfile(self._ifname_base):
            ifname_base = self._ifname_base_default

        else:
            ifname_base = self._ifname_base

        with open(ifname_base, 'r') as fsock:
            base = fsock.read()

        with open(self._ifname_setup, 'r') as fsock:
            setup = fsock.read()

        with open(self._ifname_env, 'w+') as fsock:
            fsock.write(base)
            fsock.write(setup)

# Initialise
env = Environment()


if __name__ == '__main__':
    from pprint import pprint
    pprint(env.paths)
    print ''
    pprint(env.files)
