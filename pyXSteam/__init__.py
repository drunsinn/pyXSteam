#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. codeauthor::  drunsinn <dr.unsinn@googlemail.com>
.. Licensed under the GNU General Public License v2 (GPLv2)
.. Source at https://github.com/drunsinn/pyXSteam

this is the pyXSteam package"""

import logging

from .XSteam import XSteam
from .XSteam_HW import XSteam_HW
from .UnitConverter import UnitConverter
from .Constants import UnitSystem, IceType

logging.getLogger(__name__).addHandler(logging.NullHandler())

__version__ = "1.0.0"
__author__ = "drunsinn"
__license__ = "MIT"
__email__ = "dr.unsinn@googlemail.com"
