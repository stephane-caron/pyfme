#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2017 Stephane Caron <stephane.caron@normalesup.org>
#
# This file is part of PyFME <https://github.com/stephane-caron/PyFME>.
#
# PyFME is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# PyFME is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with PyFME. If not, see <http://www.gnu.org/licenses/>.

import IPython
import os
import sympy

try:
    import pyfme
except ImportError:
    import sys
    script_path = os.path.realpath(__file__)
    sys.path.append(os.path.dirname(script_path) + '/../')
    import pyfme

"""
This example derives automatically the calculation of the single-contact wrench
cone reported in <https://scaron.info/research/icra-2015.html>.
"""

X = sympy.Symbol('X', real=True, positive=True)
Y = sympy.Symbol('Y', real=True, positive=True)
mu = sympy.Symbol('mu', real=True, positive=True)

F = sympy.Matrix([
    # fx  fy   fz
    [-1,  0, -mu],
    [+1,  0, -mu],
    [0,  +1, -mu],
    [0,  -1, -mu]])

T = sympy.diag(F, F, F, F)

M = sympy.Matrix([
    # f1x f1y f1z f2x f2y f2z f3x f3y f3z f4x f4y f4z
    [1,    0,  0,  1,  0,  0,  1,  0,  0,  1,  0,  0],
    [0,    1,  0,  0,  1,  0,  0,  1,  0,  0,  1,  0],
    [0,    0,  1,  0,  0,  1,  0,  0,  1,  0,  0,  1],
    [0,    0, +Y,  0,  0, -Y,  0,  0, -Y,  0,  0, +Y],
    [0,    0, -X,  0,  0, -X,  0,  0, +X,  0,  0, +X],
    [-Y,  +X,  0, +Y, +X,  0, +Y, -X,  0, -Y, -X, +0]])

if __name__ == "__main__":
    interactive = '-ni' not in sys.argv
    figpath = "wrench_cone.pdf"
    cone = pyfme.Cone(T, M)
    dag = pyfme.ReductionDAG(cone)
    dag.savefig(figpath)
    if interactive:
        print "First, open the file '%s' in your favorite PDF viewer." % figpath
        raw_input("Then, press [Enter] ")

    pivot_seq = [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)]
    for (node_id, pivot_col) in pivot_seq:
        if interactive:
            print "\nMake sure you refresh '%s' in your PDF viewer." % figpath
            raw_input("Press [Enter] to call pivot(%d, %d) " % (
                node_id, pivot_col))
        dag.pivot(node_id, pivot_col)
        dag.savefig(figpath)

    T = dag.nodes[-1].cone.get_matrix()
    assert T.shape[1] == 6  # output is a 6D wrench

    print "\nCalculations complete!"
    print "The resulting wrench cone has %d inequalities." % T.shape[0]
    print "Its formula is given by:"
    print repr(T)

    if IPython.get_ipython() is None:
        IPython.embed()
