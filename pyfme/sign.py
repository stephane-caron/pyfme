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

import sympy


class UnknownSign(Exception):

    def __init__(self, expr=None):
        # default expr=None because of http://bugs.python.org/issue1692335
        self.expr = expr if expr else "?"

    def __str__(self):
        return str(self.expr)

    def __repr__(self):
        return "UnknownSign(%s)" % str(self.expr)


def get_sign_expr(sign):
    if sign.func == sympy.functions.elementary.complexes.sign:
        assert len(sign.args) == 1
        return sign.args[0]
    minus_sign = -sign
    assert minus_sign.func == sympy.functions.elementary.complexes.sign
    assert len(minus_sign.args) == 1
    return -minus_sign.args[0]


def eval_sign(expr, pos_polynoms=None):
    sign = sympy.sign(expr.factor())
    if sign in [0, -1, +1]:
        return sign
    if pos_polynoms is None:
        pos_polynoms = []
    expr = get_sign_expr(sign)
    for pos_poly in pos_polynoms:
        quot, rem = sympy.div(expr, pos_poly)
        if quot == 0:
            continue
        try:
            squot = eval_sign(quot, pos_polynoms)
            srem = eval_sign(rem, pos_polynoms)
            if squot * srem >= 0:
                return squot
        except UnknownSign:
            continue
    raise UnknownSign(expr)


def eval_sign_or_none(expr, pos_polynoms=None):
    try:
        return eval_sign(expr, pos_polynoms)
    except UnknownSign:
        return None


def has_unknown_sign(expr, pos_polynoms=None):
    try:
        eval_sign(expr, pos_polynoms)
        return False
    except UnknownSign:
        return True
