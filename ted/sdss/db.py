#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Sat 1 Feb 2014
#


import sqlite3

import pandas as pd

from .. import env

_db = env.files('db')


def get_table(name=None):
    with sqlite3.connect(_db) as con:
        cur = con.cursor()
        cur.execute('SELECT * from ?', name)
        rows = cur.fetchall()
    return rows


def get_table_as_df(name=None):
    with sqlite3.connect(_db) as con:
        df = pd.read_sql('SELECT * from {}'.format(name), con, con=con)
    return df


def create_table_snlist():
    pass


def create_table_fields():
    pass


def create_table_frames():
    pass


def create_table_galaxies():
    pass


def create_table_gxlist():
    pass


def create_table_tlist():
    pass


