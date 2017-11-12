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

import copy
import itertools
import logging
import multiprocessing
import os
import sympy

from cone import PivotSet, Cone
from sign import UnknownSign


def elim_cost_to_dhm(score):
    mn = score / 200
    dhm = (mn / 60 / 24, mn / 60 % 24, mn % 60)
    return \
        ("%d d " % dhm[0] if dhm[0] > 0 else "") + \
        ("%d h " % dhm[1] if dhm[1] > 0 else "") + \
        "%d min" % dhm[2]


class SignSet(object):

    def __init__(self, signs=[]):
        self.signs = set()
        self.extend(signs)

    def __contains__(self, sign):
        e = sign.expand()
        return e in self.signs or -e in self.signs

    def add(self, sign):
        e = sign.expand()
        if e not in self:
            self.signs.add(e)

    def extend(self, signs):
        for sign in signs:
            self.add(sign)

    def __iter__(self):
        return self.signs.__iter__()

    def next(self):
        return self.signs.next()


class UnknownSigns(object):

    def __init__(self, signs):
        self.signs = SignSet(signs)


def reduction_node_recompute_pivots_map(args):
    try:
        cone, pivot_col = args
        pivot_set = PivotSet(cone, pivot_col)
        return pivot_col, pivot_set
    except UnknownSign:
        signs = cone.get_unknown_signs(col=pivot_col)
        return pivot_col, UnknownSigns(signs)


class ReductionNode(object):

    def __init__(self, node_id, cone, parent, parent_pivot_col=None,
                 parent_sign_expr=None, imbert_effect=None, elim_order=None,
                 compute_pivots=True):
        self.cone = cone
        self.dotted = False
        self.elim_order = [] if elim_order is None else elim_order
        self.elim_cost = {}
        self.elim_pivots = {}
        self.imbert_effect = imbert_effect
        self.is_test = False
        self.node_id = node_id
        self.not_pivoted = {}
        self.parent = parent
        self.parent_pivot_col = parent_pivot_col
        self.parent_sign_expr = parent_sign_expr
        self.sons = {}
        if cone is not None and compute_pivots:
            self.recompute_pivots()

    def copy(self):
        copycat = ReductionNode(
            self.node_id, self.cone, self.parent, self.parent_pivot_col,
            self.parent_sign_expr, self.imbert_effect, self.elim_order,
            compute_pivots=False)
        copycat.not_pivoted = copy.copy(self.not_pivoted)
        copycat.sons = copy.copy(self.sons)
        return copycat

    @property
    def depth(self):
        return 1 + self.parent.depth if self.parent else 0

    @property
    def is_leaf(self):
        return len(self.sons) < 1

    def get_var_elim_order(self):
        if 'elim_order' in self.__dict__:
            pass
        else:  # TMP
            print "Reconstructing node(%d).elim_order..." % self.node_id
            if self.parent is None:  # root node
                self.elim_order = []
            elif self.parent_pivot_col is None:  # assumption node
                self.elim_order = self.parent.get_var_elim_order()
            else:
                var_id = self.parent.cone.input_ids[self.parent_pivot_col]
                self.elim_order = self.parent.get_var_elim_order() + [var_id]
        return self.elim_order

    @property
    def elim_set(self):
        return set(self.elim_order)

    def count_signs(self):
        return self.cone.count_signs()

    @property
    def available_pivots(self):
        return [(pivot_col, pivot_set)
                for pivot_col, pivot_set in self.not_pivoted.iteritems()
                if type(pivot_set) is PivotSet]

    def get_pivot_cols(self):
        return [pivot_col
                for pivot_col, pivot_set in self.not_pivoted.iteritems()
                if type(pivot_set) is PivotSet]

    def get_usign_cols(self):
        return [(pivot_col, x.signs)
                for pivot_col, x in self.not_pivoted.iteritems()
                if type(x) is UnknownSigns]

    def get_all_unknown_signs(self):
        signs = SignSet()
        for pivot_col, pivot_signs in self.get_usign_cols():
            signs.extend(pivot_signs)
        return signs

    def nb_pivots_needing(self, sign_expr):
        return len([1 for pivot_col, pivot_signs in self.get_usign_cols()
                    if sign_expr in pivot_signs])

    def find_best_test(self, pivot_col, base_symbols):
        from sign import eval_sign
        pivot_signs = self.not_pivoted[pivot_col].signs
        # pivot_signs = self.get_all_unknown_signs()

        def sign_score(sign_expr):
            score = [0, 0]
            pos_case = self.cone.pos_polynoms + [+sign_expr]
            neg_case = self.cone.pos_polynoms + [-sign_expr]
            for expr in pivot_signs:
                reward = self.nb_pivots_needing(expr) ** 2
                try:
                    eval_sign(expr, pos_case)
                    score[0] += reward
                except UnknownSign:
                    pass
                try:
                    eval_sign(expr, neg_case)
                    score[1] += reward
                except UnknownSign:
                    pass
            return score

        candidate_signs = SignSet(pivot_signs)
        for expr in pivot_signs:
            for Z in base_symbols:
                q, r = sympy.div(expr, Z, Z)
                if q not in candidate_signs and self.cone.has_unknown_sign(q):
                    candidate_signs.add(q)
                if r not in candidate_signs and self.cone.has_unknown_sign(r):
                    candidate_signs.add(r)

        l = sorted([(cs, sign_score(cs)) for cs in candidate_signs], key=lambda
                   x: -(x[1][0] * x[1][1]))
        for (cs, score) in l:
            print "(%2d, %2d): %s" % (score[0], score[1], cs)
        return l[0][0]

    def expansion_if_elim(self, var_id):
        for pivot_col, x in self.not_pivoted.iteritems():
            if self.cone.input_ids[pivot_col] == var_id:
                if type(x) is not PivotSet:
                    raise x
                return x.worst_case_size
        return 0

    def has_pivots(self):
        for pivot_col, pivot_set in self.not_pivoted.iteritems():
            if type(pivot_set) is PivotSet:
                return True
        for son in self.sons.itervalues():
            if son.has_pivots():
                return True
        return False

    def recompute_pivots(self):
        logging.info("recompute_pivots() for node %d" % self.node_id)
        self.not_pivoted = {}
        nb_proc = self.cone.nb_inputs
        pool = multiprocessing.Pool(nb_proc + 1)
        args = ((self.cone, pcol) for pcol in xrange(self.cone.nb_inputs))
        res_it = pool.imap_unordered(reduction_node_recompute_pivots_map, args)
        for x in res_it:
            pivot_col, pivot_set = x
            self.not_pivoted[pivot_col] = pivot_set
        pool.close()
        pool.join()

    def pivot(self, pivot_col, next_id):
        var_id = self.cone.input_ids[pivot_col]
        pivot_set = self.not_pivoted.pop(pivot_col)
        worst_case_size = pivot_set.worst_case_size
        son_cone = pivot_set.compute_cone()
        node = ReductionNode(
            next_id, son_cone, parent=self, parent_pivot_col=pivot_col,
            imbert_effect=(worst_case_size, len(son_cone.rows)),
            elim_order=self.elim_order + [var_id])
        self.sons[pivot_col] = node
        self.remove_unknown_pivots()
        return node

    def remove_unknown_pivots(self):
        del_keys = [pivot_col for pivot_col, x in self.not_pivoted.iteritems()
                    if type(x) is not PivotSet]
        for k in del_keys:
            del self.not_pivoted[k]

    def make_test(self, sign_expr, next_id1, next_id2):
        def append_son(sign_expr, next_id):
            cone = self.cone.copy()
            cone.add_pos_assumption(sign_expr)
            node = ReductionNode(
                next_id, cone, self, parent_sign_expr=sign_expr,
                elim_order=self.elim_order)
            self.sons["%s >= 0" % sign_expr] = node
            return node
        assert not self.sons, "Node %d has descendants already" % self.node_id
        node1 = append_son(+sign_expr, next_id1)
        node2 = append_son(-sign_expr, next_id2)
        self.not_pivoted = {}
        self.is_test = True
        return node1, node2

    def compute_elimination_cost(self, var_id):
        if var_id in self.elim_cost:
            return self.elim_cost[var_id]
        tmp_dag = ReductionDAG(self.cone)
        tmp_dag.recursive_make_test_leaves(var_id)
        plist = tmp_dag.get_son_leaf_pivots(var_id)
        psets = (tmp_dag.nodes[nid].not_pivoted[pcol] for nid, pcol in plist)
        score = sum(pset.worst_case_size for pset in psets)
        self.elim_cost[var_id] = score
        self.elim_pivots[var_id] = len(plist)

    def compute_elimination_costs(self):
        for var_id in self.cone.input_ids:
            if var_id not in self.elim_cost:
                self.compute_elimination_cost(var_id)
        return self.elim_cost

    def to_dot(self, next_id=10000, leaves_only=False):
        dot_str = ""
        if self.dotted:
            return dot_str, next_id
        write_self = self.is_leaf or not leaves_only
        if write_self:
            next_str, next_id = self.pivots_to_dot(next_id)
            dot_str = self.dot_str + next_str
        for pivot_col, son in self.sons.iteritems():
            if son.parent_pivot_col is not None:
                var_id = self.cone.input_ids[pivot_col]
                edge_str = "  %d -> %d [label=\"z%d\"];\n" % (
                    self.node_id, son.node_id, var_id)
                next_str, next_id = son.to_dot(
                    next_id, leaves_only=leaves_only)
                dot_str += (edge_str if write_self else "") + next_str
            elif son.parent_sign_expr:
                sign_expr = son.parent_sign_expr
                edge_str = "  %d -> %d [" % (self.node_id, son.node_id)
                edge_str += "color=blue,fontcolor=blue,"
                edge_str += "label=\"(%s >= 0)\"];\n" % str(sign_expr)
                next_str, next_id = son.to_dot(
                    next_id, leaves_only=leaves_only)
                dot_str += (edge_str if write_self else "") + next_str
        self.dotted = True
        return dot_str, next_id

    @property
    def dot_label(self):
        label = "Node %d\\n" % self.node_id
        label += "%d hyperplanes\\n" % len(self.cone.rows)
        elim_label = ", ".join(["z%d" % i for i in self.elim_order])
        if elim_label:
            label += "Eliminated %s\\n" % elim_label
        if 'sign_count' in self.cone.__dict__:  # kron
            label += "\\nSigns:\\n"
            for input_id, d in self.cone.sign_count.iteritems():
                label += "z%d: %3d+ / %3d- / %3dz / %3d? / %3d?u\\n" % (
                    input_id, d['+'], d['-'], d['0'], d['?'], d['?u'])
        if 'elim_cost' not in self.__dict__:
            self.elim_cost = {}
        if self.elim_cost:
            label += "\\nElim. costs:\\n"
            label += "\\n".join(["z%d: +%d (%d p / %s)" % (
                k, v, self.elim_pivots[k], elim_cost_to_dhm(v))
                for k, v in self.elim_cost.iteritems()])
            label += "\\n"
        if self.cone.gen_pos_polynoms:
            label += "\\nFree cone:\\n"
            for (i, pos_poly) in enumerate(self.cone.gen_pos_polynoms):
                label += "%s >= 0  (%d)\\n" % (str(pos_poly), i)
        if self.cone.inf_pos_polynoms:
            label += "\\nInferred:\\n"
            for pos_poly in self.cone.inf_pos_polynoms:
                label += "%s >= 0\\n" % str(pos_poly)
        return label

    @property
    def dot_str(self):
        dot_str = "  %d [style=filled,shape=box," % self.node_id
        color_str = \
            "#FFBB66" if self.parent_pivot_col is not None else \
            "#9999FF" if self.is_test else \
            "#CCCCCC"
        dot_str += "color=\"%s\",fontsize=18," % color_str
        dot_str += "label=\"%s\"];\n" % self.dot_label
        return dot_str

    def pivots_to_dot(self, next_id):
        dot_str = ""
        for pivot_col, x in self.not_pivoted.iteritems():
            var_id = self.cone.input_ids[pivot_col]
            if type(x) is PivotSet and x.worst_case_size >= 0:
                label = "pivot(%d, %d)\\n" % (self.node_id, pivot_col)
                label += "+%d\\n" % x.worst_case_size
                label += "(%s)" % elim_cost_to_dhm(x.worst_case_size)
                dot_str += "  %d -> %d [" % (self.node_id, next_id)
                dot_str += "color=yellow,fontcolor=\"#FFBB66\",weight=2,"
                dot_str += "label=\"z%d\"];\n" % var_id
                dot_str += "  %d [style=filled,color=" % next_id
                dot_str += "yellow,label=\"%s\"]" % label
            elif type(x) is UnknownSigns:
                label = "\\n".join([str(e) for e in x.signs])
                dot_str += "  %d -> %d [color=red," % (self.node_id, next_id)
                dot_str += "fontcolor=red,weight=2,label=\"z%d\"];\n" % var_id
                dot_str += "  %d [style=filled,color=" % next_id
                dot_str += "\"#FF9999\",label=\"%s\"];\n" % label
            else:
                assert type(x) is UnknownSign
                print "Node %d needs to compute unknown signs (pivot_col=%d)" \
                    % (self.node_id, pivot_col)
            next_id += 1
        return dot_str, next_id

    def __str__(self):
        return "Node %d (%s)" % (self.node_id, str(self.cone).lower())

    def __repr__(self):
        return str(self)


class ReductionDAG(object):

    def __init__(self, cone=None, node=None, base_symbols=None):
        assert cone or node
        if cone is not None:
            assert type(cone) is Cone
            self.root = ReductionNode(0, cone, parent=None)
        elif node is not None:
            assert type(node) is ReductionNode
            self.root = ReductionNode(0, node.cone, parent=None)
            self.root.elim_order = node.elim_order
        self.nodes = [self.root]
        self.base_symbols = []
        self.elim_cost = {}
        self.elim_pivots = {}
        if base_symbols is not None:
            self.base_symbols = base_symbols

    def copy(self):
        copycat = ReductionDAG(cone=None)
        copycat.root = self.root.copy()
        copycat.nodes = [node.copy() for node in self.nodes]
        return copycat

    @property
    def internal_nodes(self):
        for node in self.nodes:
            if not node.is_leaf:
                yield node

    @property
    def leaves(self):
        for node in self.nodes:
            if node.is_leaf:
                yield node

    @property
    def nb_leaves(self):
        return len(list(self.leaves))

    def merge_equiv_nodes(self):
        def get_node_pools(node):
            other_pools = []
            my_pool = [] if node.is_test else [node]
            for son in node.sons.itervalues():
                son_pool, other_pools2 = get_node_pools(son)
                if node.is_test:
                    other_pools.append(son_pool)
                else:
                    my_pool.extend(son_pool)
                other_pools.extend(other_pools2)
            return my_pool, other_pools

        root_pool, other_pools = get_node_pools(self.root)
        node_pools = filter(None, [root_pool] + other_pools)
        for pool in node_pools:
            node_pairs = ((n1, n2) for n1 in pool for n2 in pool)
            for (n1, n2) in node_pairs:
                if n1.elim_set == n2.elim_set \
                        and len(n1.cone.rows) < len(n2.cone.rows):
                    print "Node %d <- Node %d" % (n2.node_id, n1.node_id)
                    pivot_col = n2.parent_pivot_col
                    n2.parent.sons[pivot_col] = n1

    def pivot(self, node_id, pivot_col):
        logging.info("dag.pivot(%d, %d)" % (node_id, pivot_col))
        next_node_id = len(self.nodes)
        new_node = self.nodes[node_id].pivot(pivot_col, next_node_id)
        if new_node:
            self.nodes.append(new_node)

    def pivot_all(self, root_id=0, var_id=None, max_depth=None):
        """
        Minimize estimated number of hyperplanes.

        Parameters
        ----------
        root_id : int
            Index of node to start pivoting from.
        var_id : int, optional
            Index of variable to pivot.
        max_depth : int, optional
            Maximum recursion depth.
        """
        nodes = [self.root]
        while nodes:
            node = nodes.pop(0)
            if max_depth and node.depth >= max_depth:
                continue
            for pivot_col in node.get_pivot_cols():
                if var_id is None or node.cone.input_ids[pivot_col] == var_id:
                    print "dag.pivot(%d, %d)" % (node.node_id, pivot_col)
                    self.pivot(node.node_id, pivot_col)
            nodes.extend(node.sons.values())

    def get_son_leaf_pivots(self, var_id, node_id=None):
        node_id = node_id if node_id is not None else 0
        node = self.nodes[node_id]
        input_ids = node.cone.input_ids
        node_pivots = node.get_pivot_cols()
        if not node.is_leaf:
            son_pivots = [self.get_son_leaf_pivots(var_id, son.node_id)
                          for son in node.sons.itervalues()]
            return list(itertools.chain(*son_pivots))
        return [(node_id, pivot_col) for pivot_col in node_pivots
                if var_id is None or input_ids[pivot_col] == var_id]

    def pivot_leaves(self, var_id, root_id=None, ask_confirmation=True):
        root_id = 0 if root_id is None else root_id
        queue = self.get_son_leaf_pivots(var_id, root_id)
        total_size = 0
        print "%d pivots to be performed:" % len(queue)
        for node_id, pivot_col in queue:
            node = self.nodes[node_id]
            pivot_size = node.not_pivoted[pivot_col].worst_case_size
            print "- pivot(%d, %d): +%d" % (node_id, pivot_col, pivot_size)
            total_size += pivot_size
        print "Estimated time: %s" % elim_cost_to_dhm(total_size)
        if not ask_confirmation or raw_input("OK? [y/N] ") in ['y', 'Y']:
            for node_id, pivot_col in queue:
                print "dag.pivot(%d, %d)" % (node_id, pivot_col)
                self.pivot(node_id, pivot_col)

    def make_test(self, node_id, sign_expr):
        logging.info("dag.make_test(%d, %s)" % (node_id, sign_expr))
        node = self.nodes[node_id]
        next_id1 = len(self.nodes)
        next_id2 = len(self.nodes) + 1
        node1, node2 = node.make_test(sign_expr, next_id1, next_id2)
        self.nodes.append(node1)
        self.nodes.append(node2)

    def recursive_make_test_leaves(self, var_id, node_id=None):
        node_id = node_id if node_id is not None else 0
        node = self.nodes[node_id]
        input_ids = node.cone.input_ids
        if not node.is_leaf:
            for son in node.sons.itervalues():
                self.recursive_make_test_leaves(var_id, son.node_id)
        else:
            for pivot_col, pivot_usigns in node.get_usign_cols():
                if var_id is None or input_ids[pivot_col] == var_id:
                    usign_expr = node.find_best_test(
                        pivot_col, self.base_symbols)
                    print "Testing %s at node %d" % (usign_expr, node_id)
                    self.make_test(node_id, usign_expr)
                    return self.recursive_make_test_leaves(var_id, node_id)

    def eliminate_variable(self, var_id):
        self.recursive_make_test_leaves(var_id)
        self.pivot_leaves(var_id, ask_confirmation=True)

    def count_signs(self):
        for node in self.nodes:
            node.count_signs()

    def compute_elimination_cost(self, var_id):
        elim_cost = 0
        elim_pivots = 0
        for leaf in self.leaves:
            for pcol, x in leaf.not_pivoted.iteritems():
                if leaf.cone.input_ids[pcol] == var_id:
                    assert type(x) is PivotSet, "make all tests first"
                    elim_cost += x.worst_case_size
                    elim_pivots += 1
        return elim_cost, elim_pivots

    def compute_elimination_costs(self):
        elim_cost = {var_id: 0 for var_id in self.root.cone.input_ids}
        elim_pivots = {var_id: 0 for var_id in self.root.cone.input_ids}
        for leaf in self.leaves:
            leaf.compute_elimination_costs()
            for var_id in leaf.elim_cost:
                elim_cost[var_id] += leaf.elim_cost[var_id]
                elim_pivots[var_id] += leaf.elim_pivots[var_id]
        self.elim_cost = elim_cost
        self.elim_pivots = elim_pivots
        return elim_cost

    def savefig(self, fname, root_id=None, leaves_only=False):
        """
        Save a drawing of the DAG as a PDF figure.

        Parameters
        ----------
        fname : string
            File name (must end with ".pdf").
        root_id : int, optional
            Root node index.
        leaves_only : bool, optional
            If True, only leaf nodes are drawn.
        """
        assert fname.endswith('.pdf')
        for node in self.nodes:
            node.dotted = False
        root_id = root_id if root_id is not None else 0
        if 'elim_cost' not in self.__dict__:
            self.elim_cost = {}
        if 'elim_pivots' not in self.__dict__:
            self.elim_pivots = {i: -1 for i in self.root.cone.input_ids}
        dot_str = ""
        if self.elim_cost:
            label = "DAG\\n\\n"
            elim = [var_id
                    for var_id, v in self.elim_cost.iteritems()
                    if v == 0]
            label += "\\n".join(["Eliminated z%d" % var_id for var_id in elim])
            label += "\\n\\n" if elim else "\\n"
            label += "\\n".join([
                "z%d: +%6d / %d leaves / %s" % (
                    var_id, v, self.elim_pivots[var_id], elim_cost_to_dhm(v))
                for var_id, v in self.elim_cost.iteritems()
                if v > 0])
            label += "\\n"
            dot_str += "  DAG [style=filled,shape=box,"
            dot_str += "color=\"#AACCFF\",fontsize=24,"
            dot_str += "label=\"%s\"];\n" % label
        dot_str += self.nodes[root_id].to_dot(leaves_only=leaves_only)[0]
        with open("graphviz.dot", "w") as f:
            f.write("digraph G {\n%s\n}\n" % dot_str)
        os.system('dot -Tpdf graphviz.dot -o %s' % fname)
        os.unlink('graphviz.dot')
