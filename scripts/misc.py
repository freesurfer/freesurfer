# Original author - Krish Subramaniam
# $Id: misc.py,v 1.1 2009/02/05 21:49:02 krish Exp $

"""
This module essentially contains useful classes and functions that don't fit
elsewhere. 
"""
from __future__ import print_function
__all__ = ['callback_var']

# optparse can't handle variable number of arguments for an option.
# this callback allows that.
def callback_var(option, opt_str, value, parser):
    value = []
    rargs = parser.rargs
    while rargs:
        arg = rargs[0]
        if ((arg[:2] == '--' and len(arg) > 2) or
            (arg[:1] == '-' and len(arg) > 1 and arg[1] != '-')):
            break
        else:
            value.append(arg)
            del rargs[0]
    setattr(parser.values, option.dest, value)
