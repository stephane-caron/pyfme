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

import logging
import sympy

import symparallel
from sign import eval_sign as virgin_eval_sign
from sign import has_unknown_sign as virgin_has_unknown_sign
from sign import UnknownSign
from symparallel import slice_and_pool


def right_inverse(M):
    return M.T * (M * M.T).inv()


class ConeRow(object):

    """
    Row of a polyhedral convex cone.

    Parameters
    ----------
    v : sympy.Matrix
        Row vector of shape (1, m).
    H : set
        Historical subset.
    E : set
        Set of effectively eliminated variables.
    I : set
        Set of implicitely eliminated variables.
    """

    def __init__(self, v, H=None, E=None, I=None):
        assert type(v) is sympy.Matrix and v.shape[0] == 1
        self.E = set() if E is None else E
        self.I = set() if I is None else I
        self.H = set() if H is None else H
        self.v = v

    def copy(self, v=None):
        v = v if v is not None else self.v
        copy_row = ConeRow(v)
        copy_row.E = set(self.E)
        copy_row.I = set(self.I)
        copy_row.H = set(self.H)
        return copy_row

    def without_col(self, col):
        v2 = sympy.Matrix(self.v)
        v2.col_del(col)
        return self.copy(v2)

    def __str__(self):
        s = str(self.v)
        s = s.replace("Matrix([[", "ConeRow([")
        s = s.replace("]])", "])")
        return s


class PivotSet(object):

    def __init__(self, cone, pivot_col):
        self.cone = cone
        self.lbounds = []
        self.kohler_matrix = cone.kohler_matrix
        self.pos_polynoms = cone.pos_polynoms
        self.gen_pos_polynoms = cone.gen_pos_polynoms
        self.inf_pos_polynoms = cone.inf_pos_polynoms
        self.rows = []
        self.ubounds = []

        self.pivot_col = pivot_col
        self.pivot_input_id = cone.input_ids[pivot_col]
        self.input_ids = list(cone.input_ids)
        self.input_ids.remove(self.pivot_input_id)

        self.officially_elim = set(cone.officially_elim)
        self.officially_elim.add(self.pivot_input_id)

        for row in cone.rows:
            fact = row.v[pivot_col]
            fact_sign = cone.eval_sign(fact)
            if fact_sign == 0:
                self.rows.append(row.without_col(self.pivot_col))
            elif fact_sign < 0:
                self.lbounds.append((row, fact))
            else:  # fact_sign > 0
                self.ubounds.append((row, fact))

    @property
    def worst_case_size(self):
        return len(self.lbounds) * len(self.ubounds) + len(self.rows)

    def compute_rows(self):
        print "PivotSet: %d lower bounds" % len(self.lbounds),
        print " and %d upper bounds" % len(self.ubounds)
        init_rows = len(self.rows)

        for (li, lbound) in enumerate(self.lbounds):
            msg = "PivotSet: combinations for lower bound %d / %d" % (
                li, len(self.lbounds))
            input_list = self.ubounds

            def prepare_args(cur_ubounds):
                return ((self.cone, self.officially_elim, self.pivot_col,
                         self.pivot_input_id, lbound, ubound)
                        for ubound in cur_ubounds)

            def callback(new_row):
                if new_row is not None:
                    self.rows.append(new_row)

            slice_and_pool(input_list, prepare_args, pivot_set_combine_bounds,
                           callback, msg)

        print "PivotSet => %d new rows" % (len(self.rows) - init_rows)

    def compute_cone(self):
        self.compute_rows()
        symparallel.simplify_rows_gcd(self.rows, self.cone.pos_polynoms)
        symparallel.simplify_rows(self.rows)
        symparallel.remove_duplicates(self.rows)
        return Cone(origin=self)


def pivot_set_combine_bounds(args):
    cone, O, pivot_col, pivot_input_id, lbound, ubound = args
    (lrow, lfact), (urow, ufact) = lbound, ubound

    H = set.union(lrow.H, urow.H)
    E = set.union(lrow.E, urow.E)
    E.add(pivot_input_id)

    v = lrow.v * ufact - urow.v * lfact
    assert type(v) is sympy.Matrix
    assert v[pivot_col] == 0
    v.simplify()  # needed for detection of implicit eliminations

    I = set.union(lrow.I, urow.I) - E
    for i in xrange(cone.nb_inputs):
        input_id = cone.input_ids[i]
        if input_id in I or input_id in E:
            continue
        elif v[i] == 0 and lrow.v[i] != 0 and urow.v[i] != 0:
            I.add(input_id)

    imbert_lb = 1 + len(E)
    imbert_ub = 1 + len(set.union(E, set.intersection(I, O)))

    if len(H) > imbert_ub:  # Imbert's first acceleration theorem
        logging.debug("PivotSet: used Imbert's #1")
        return None

    if len(H) > imbert_lb:  # Imbert's second acceleration theorem
        kohler_rank = cone.kohler_matrix[list(H), :].rank()
        if kohler_rank != len(H) - 1:
            logging.debug("PivotSet: used Kohler's criterion")
            return None
    else:  # len(H) == imbert_lb
        logging.debug("PivotSet: used Imbert's #2")
        pass

    v.col_del(pivot_col)
    return ConeRow(v, H, E, I)


class Cone(object):

    """
    Polyhedral convex cone described by:

    .. math::

        T x \leq 0

    The goal of the elimination is to convert `x` variables into their
    projections :math:`y = M x`.

    Parameters
    ----------
    T : sympy.Matrix
        Constraint matrix on input variables.
    M : sympy.Matrix
        Linear mapping from input to output variables.
    origin : Cone or PivotSet, optional
        Used for initialization from another object.
    """

    def __init__(self, T=None, M=None, origin=None):
        self.gen_pos_polynoms = []
        self.inf_pos_polynoms = []
        self.pos_polynoms = []
        self.unknown_sign_counts = {}
        if T is not None and M is not None:
            return self.init_from_matrices(T, M)
        elif origin is not None:
            return self.init_from_origin(origin)
        raise Exception("no empty constructor for Cone")

    def init_from_matrices(self, T, M):
        assert T.shape[1] == M.shape[1]
        if M.shape[1] < M.shape[0]:
            raise ValueError("FME not needed as M has more rows than columns")
        Mrinv = right_inverse(M)
        Mrinv.simplify()

        ker_basis = M.nullspace()
        ker_dim = len(ker_basis[0])
        K = sympy.Matrix(ker_dim, len(ker_basis), lambda i, j: ker_basis[j][i])
        assert (M * K).norm().simplify() == 0

        U, V = T * K, T * Mrinv   # Constraint:  U * z + V * w  <= 0
        W = U.row_join(V)         # Constraint:  W * [ z  w ] <= 0
        nb_rows = W.shape[0]

        self.input_ids = range(len(ker_basis))
        self.kohler_matrix = U
        self.officially_elim = set([])
        self.rows = [ConeRow(W[i, :], H=set([i])) for i in xrange(nb_rows)]

        symparallel.simplify_rows_lcm(self.rows, self.pos_polynoms)
        symparallel.simplify_rows(self.rows)
        symparallel.simplify_rows_gcd(self.rows, self.pos_polynoms)
        symparallel.simplify_rows(self.rows)

    def init_from_origin(self, origin):
        self.inf_pos_polynoms = list(origin.inf_pos_polynoms)
        self.input_ids = list(origin.input_ids)
        self.gen_pos_polynoms = list(origin.gen_pos_polynoms)
        self.kohler_matrix = sympy.Matrix(origin.kohler_matrix)
        self.officially_elim = set(origin.officially_elim)
        self.pos_polynoms = list(origin.pos_polynoms)
        self.rows = list(origin.rows)
        if 'unknown_sign_counts' not in origin.__dict__:
            origin.unknown_sign_counts = {}
        self.unknown_sign_counts = dict(origin.unknown_sign_counts)

    def copy(self):
        return Cone(origin=self)

    @property
    def nb_inputs(self):
        return len(self.input_ids)

    def __str__(self):
        s = "Cone with %d hyperplanes " % len(self.rows)
        s += "and %d inputs" % self.nb_inputs
        return s

    def add_pos_assumption(self, poly):
        self.gen_pos_polynoms.append(poly)
        self.pos_polynoms.append(poly)

    def infer_pos_polynom(self, comb):
        assert len(comb) == len(self.gen_pos_polynoms)
        for x in comb:
            assert self.eval_sign(sympy.S(x)) >= 0, str(x)
        gen_poly = self.gen_pos_polynoms
        new_poly = sum([v_i * gen_poly[i] for (i, v_i) in enumerate(comb)])
        self.inf_pos_polynoms.append(new_poly)
        self.pos_polynoms.append(new_poly)

    def eval_sign(self, expr):
        return virgin_eval_sign(expr, self.pos_polynoms)

    def has_unknown_sign(self, expr):
        return virgin_has_unknown_sign(expr, self.pos_polynoms)

    def get_matrix(self):
        return sympy.Matrix([r.v for r in self.rows])

    def get_unknown_signs(self, col=None, input_id=None):
        assert col is not None or input_id is not None
        if col is None:
            col = self.input_ids.index(input_id)
        unknown_signs = set()
        for r in self.rows:
            try:
                self.eval_sign(r.v[col])
            except UnknownSign as e:
                unknown_signs.add(e.expr.expand())
        return unknown_signs

    def count_signs(self):
        if 'sign_count' not in self.__dict__:
            self.sign_count = {}
        if self.sign_count:
            return self.sign_count
        print "count_signs(%s)" % str(self)
        for pivot_col, input_id in enumerate(self.input_ids):
            self.sign_count[input_id] = {'+': 0, '-': 0, '?': 0, '0': 0}
            usign_set = set()
            for r in self.rows:
                try:
                    expr = r.v[pivot_col]
                    sign = self.eval_sign(expr)
                    c = '+' if sign > 0 else '-' if sign < 0 else '0'
                    self.sign_count[input_id][c] += 1
                except UnknownSign as e:
                    self.sign_count[input_id]['?'] += 1
                    usign_set.add(e.expr.expand())
            self.sign_count[input_id]['?u'] = len(usign_set)
        return self.sign_count
