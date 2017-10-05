# SymPyFME

This library can be used to [project
polytopes](https://scaron.info/teaching/projecting-polytopes.html) using
[Fourier-Motzkin
elimination](https://en.wikipedia.org/wiki/Fourier–Motzkin_elimination) with
the two Imbert acceleration theorems. It is implemented in Python using
[SymPy](http://www.sympy.org/en/index.html) for symbolic computations and
multiprocessing to leverage the high degree of parallelization achievable with
this method.

## Trying out variable eliminations

SymPyFME was mainly designed as a *graphical tool* to help a human user explore
various variable-elimination schemes to project a polytope given in symbolic
form. The user can explore the [directed acyclic
graph](https://en.wikipedia.org/wiki/Directed_acyclic_graph) (DAG) of
variable-elimination sequences, adding if-then-else nodes to the graph when
desired. A calculation DAG can be plotted at any time using
[GraphViz](http://www.graphviz.org/) by calling:

```python
dag.savefig("elim-dag.pdf")
```

It is also possible to tell the library to explore all possible computations
(at the cost of a tremendous computation time) by:

```python
dag.pivot_all()
```

The result of a full DAG compilation is either a single inequality system, or a
tree with if-then-else tests as internal nodes and inequality systems as
external nodes.

## Related libraries

For better performance on numerical systems, you can check out one of these
libraries:

- [PANDA](http://comopt.ifi.uni-heidelberg.de/software/PANDA/)
- [cdd](https://www.inf.ethz.ch/personal/fukudak/cdd_home/index.html) or [pycddlib](https://github.com/mcmtroffaes/pycddlib))
- [lrs](http://cgm.cs.mcgill.ca/~avis/C/lrs.html)
