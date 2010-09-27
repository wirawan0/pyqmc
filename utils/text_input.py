#!/usr/bin/python
# $Id: text_input.py,v 1.4 2010-09-27 19:58:29 wirawan Exp $
#
# pyqmc.utils.text_input module
# Quick-n-dirty text input utilities
#
# Wirawan Purwanto
# Created: 20090601
#
# Routines put here are commonly used in my own scripts.
# They are not necessarily suitable for general-purpose uses; evaluate
# your needs and see if they can them as well.
#
# TODO
# - book-keep the line number. Also note superfile must have its own line
#   number keeping.
#

"""
Obsoleted. Now loads wpylib.iofmt.text_input .
"""

import sys
sys.stderr.write("pyqmc.utils.text_input is obsoleted. Please use wpylib.iofmt.text_input instead.\n")
sys.stderr.flush()

from wpylib.iofmt.text_input import *
