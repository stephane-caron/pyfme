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

from distutils.core import setup

setup(
    name='PyFME',
    version='1.0.0',
    description="Fourier-Motzkin elimination with Imbert accelerations",
    url="https://github.com/stephane-caron/PyFME",
    author="Stéphane Caron",
    author_email="stephane.caron@normalesup.org",
    license="GPL",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Programming Language :: Python :: 2.7'],
    packages=['pyfme']
)
