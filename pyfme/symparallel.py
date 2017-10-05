#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2017 Stephane Caron <stephane.caron@normalesup.org>
#
# This file is part of PyFME <https://github.com/stephane-caron/PyFME>.
#
# PyFME is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PyFME is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PyFME. If not, see <http://www.gnu.org/licenses/>.

import itertools
import logging
import multiprocessing
import sympy
import time

from sign import eval_sign


def slice_and_pool(input_list, prepare_args, map_fun, callback, msg,
                   nb_proc=8, slice_len=400):
    input_len = len(input_list)
    if input_len < 1:
        return
    for slice_index in xrange(1 + input_len / slice_len):
        slice_start = slice_index * slice_len
        slice_end = min(input_len, slice_start + slice_len)
        cur_slice = itertools.islice(input_list, slice_start, slice_end)
        logging.debug("%s slice [%d, %d]" % (msg, slice_start, slice_end - 1))
        args = prepare_args(cur_slice)
        pool = multiprocessing.Pool(nb_proc)
        res_it = pool.imap_unordered(map_fun, args)
        for x in res_it:
            callback(x)
        pool.close()
        pool.join()


def simplify_rows_map(args):
    i, j, r_ij = args
    return (i, j, r_ij.simplify())


def simplify_rows(rows):
    """Simplify a matrix given as list of rows."""
    nb_rows = len(rows)
    nb_cols = len(rows[0].v)
    row_indexes = xrange(nb_rows)
    msg = "simplify_rows(%d rows)" % nb_rows

    def prepare_args(row_indexes_slice):
        return ((i, j, rows[i].v[j]) for i in row_indexes_slice for j in
                xrange(nb_cols))

    def callback(map_result):
        i, j, r_ij = map_result
        rows[i].v[j] = r_ij

    slice_and_pool(row_indexes, prepare_args, simplify_rows_map, callback, msg)


def simplify_rows_gcd_map(args):
    """
    Note
    ----
    As of sympy 0.7.5, there may be some issues if some x \in Wi are not
    simplified but simplify to an integer. See
    <https://github.com/sympy/sympy/issues/8384>.
    """
    i, Wi, pos_polynoms = args
    gcd_i = sympy.gcd([x for x in Wi if x != 0])
    abs_gcd_i = eval_sign(gcd_i, pos_polynoms) * gcd_i
    return (i, abs_gcd_i)


def simplify_rows_gcd(rows, pos_polynoms=None):
    assert type(rows) is list and type(rows[0].v) is sympy.Matrix
    nb_rows = len(rows)
    row_indexes = xrange(nb_rows)
    msg = "simplify_rows_gcd(%d rows)" % nb_rows

    def prepare_args(row_indexes_slice):
        return ((i, rows[i].v, pos_polynoms) for i in row_indexes_slice)

    def callback(map_result):
        i, abs_gcd_i = map_result
        rows[i].v /= abs_gcd_i

    slice_and_pool(
        row_indexes, prepare_args, simplify_rows_gcd_map, callback, msg)


def simplify_rows_lcm_map(args):
    i, Wi, pos_polynoms = args
    lcm_i = sympy.lcm([sympy.denom(x.cancel()) for x in Wi if x != 0])
    abs_lcm_i = eval_sign(lcm_i, pos_polynoms) * lcm_i
    return (i, abs_lcm_i)


def simplify_rows_lcm(rows, pos_polynoms=None):
    assert type(rows) is list and type(rows[0].v) is sympy.Matrix
    nb_rows = len(rows)
    row_indexes = xrange(nb_rows)
    msg = "simplify_rows_lcm(%d rows)" % nb_rows

    def prepare_args(row_indexes_slice):
        return ((i, rows[i].v, pos_polynoms) for i in row_indexes_slice)

    def callback(x):
        i, abs_lcm_i = x
        rows[i].v *= abs_lcm_i

    slice_and_pool(
        row_indexes, prepare_args, simplify_rows_lcm_map, callback, msg)


def remove_duplicates_map(args):
    j, rj, ri = args
    row_diff = (rj - ri)
    assert row_diff.shape[0] == 1
    assert row_diff.shape[1] > 1
    for i in xrange(row_diff.shape[1]):
        if row_diff[0, i].simplify() != 0:  # lazy simplify
            return None
    return j


def group_rows(rows, hash_fun):
    d = {}
    for i, row in enumerate(rows):
        key = hash_fun(row)
        if key in d:
            d[key].add(i)
        else:
            d[key] = set([i])
    return d


def remove_duplicates(rows):
    assert type(rows) is list and type(rows[0].v) is sympy.Matrix
    del_rows = set([])

    def hash_fun(row):
        return row.v[0].expand()

    rows_table = group_rows(rows, hash_fun)
    start_time = time.time()

    for i, row in enumerate(rows):
        if i in del_rows:
            continue
        key = hash_fun(row)
        next_rows = rows_table[key] - set([i])

        def prepare_args(next_rows_slice):
            return ((j, rows[j].v, rows[i].v) for j in next_rows_slice)

        def callback(map_result):
            if map_result is not None:
                del_rows.add(map_result)

        msg = "remove_duplicates(%d rows) row %d" % (len(rows), i)
        slice_and_pool(next_rows, prepare_args, remove_duplicates_map,
                       callback, msg)

    rows_in_removal_order = sorted(set(del_rows), reverse=True)
    for i in rows_in_removal_order:
        del rows[i]
    logging.debug("%d duplicate rows removed in %.1f s" % (
        len(del_rows), time.time() - start_time))
