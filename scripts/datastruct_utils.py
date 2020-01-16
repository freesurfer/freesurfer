# Original author - Krish Subramaniam
# $Id: datastruct_utils.py,v 1.6 2014/11/21 23:07:41 greve Exp $
from __future__ import generators, print_function

__all__ = ['Ddict', 'TableWriter', 'StableDict', 'unique_union',
           'intersect_order']

class Ddict(dict):
    """
    This datastructure is used to store 2d Table 
    Mainly for a*stats2table 
    Usage:
    >>> tb = Ddict( dict )
    >>> tb['bert']['putamen'] = .05
    >>> tb['bert']['caudate'] = 1.6
    >>> tb['fsaverage']['putamen'] = 2.2
    >>> car_details
    {'fsaverage': {'putamen': 2.2}, 'bert': {'putamen': 0.05, 'caudate':
        1.6}}
    """
    def __init__(self, default=None):
        self.default = default

    def __getitem__(self, key):
        if key not in self:
            self[key] = self.default()
        return dict.__getitem__(self, key)

class TableWriter:
    """
    This class writes a 2d Table of type Ddict(dict) to a file. Some essential
    things needs to be set for this class 
    rows - a sequence of text which go in the first column
    columns - a sequence of text which go in the first row
    table - a Ddict(dict) object which has *exactly* len(columns) x len(rows) elements
    row1col1 - the text which goes in 1st row and 1st column
    delimiter - what separates each element ( default is a tab )
    filename - the filename to write to.
    """
    def __init__(self, r, c, table):
        self.rows = r
        self.columns = c
        self.table = table
        self.pretext = ''
        self.posttext = ''

    def assign_attributes(self, filename='stats.table', row1col1='',
                          delimiter='\t'):
        self.filename = filename
        self.row1col1 = row1col1
        self.delimiter = delimiter

    def decorate_col_titles(self, pretext, posttext):
        self.pretext = pretext
        self.posttext = posttext

    def write(self):
        fp = open(self.filename, 'w')
        fp.write(self.row1col1)
        for c in self.columns:
            if((c == 'eTIV' or c == 'BrainSegVolNotVent') and (self.pretext == 'lh_' or self.pretext == 'rh_')):
                # For eTIV in aparc stats file
                fp.write(self.delimiter + c);
            else:
                fp.write(self.delimiter + self.pretext + c + self.posttext)
        fp.write('\n')
        
        for r in self.rows:
            fp.write(r)
            for c in self.columns:
                fp.write(self.delimiter + '%s' %self.table[r][c])
            fp.write('\n')
        fp.close()    

    def write_transpose(self):
        fp = open(self.filename, 'w')
        fp.write(self.row1col1)
        for r in self.rows:
            fp.write(self.delimiter + r)
        fp.write('\n')

        for c in self.columns:
            if((c == 'eTIV' or c == 'BrainSegVolNotVent') and (self.pretext == 'lh_' or self.pretext == 'rh_')):
                # For eTIV in aparc stats file
                fp.write(c)
            else:
                fp.write(self.pretext + c + self.posttext)
            for r in self.rows:
                fp.write(self.delimiter + '%g' %self.table[r][c])
            fp.write('\n')
        fp.close()    
        

# stabledict.py,v 1.13 2007/08/28 07:53:44 martin Exp
# This module has been tested with CPython 2.2.3, 2.3.5, 2.4.4, 2.5 and  2.5.1.
# It won't compile on Python 2.1.X or lower because of missing features.
# It won't work with Python 3.X because of the changed dict protocol. (PEP3106)
# From package StableDict, PSF licensed by Martin Kammerhofer
# http://pypi.python.org/pypi/StableDict/0.2

"""
A dictionary class remembering insertion order.

Order (i.e. the sequence) of insertions is remembered (internally
stored in a hidden list attribute) and replayed when iterating. A
StableDict does NOT sort or organize the keys in any other way.
"""

from warnings import warn as _warn

# Helper metaclass-function.  Not exported by default but accessible
# as StableDict.__metaclass__.
def copy_baseclass_docs(classname, bases, dict, metaclass=type):
    """Copy docstrings from baseclass.

    When overriding methods in a derived class the docstrings can
    frequently be copied from the base class unmodified.  According to
    the DRY principle (Don't Repeat Yourself) this should be
    automated. Putting a reference to this function into the
    __metaclass__ slot of a derived class will automatically copy
    docstrings from the base classes for all doc-less members of the
    derived class.
    """
    for (name, member) in dict.iteritems():
        if getattr(member, "__doc__", None):
            continue
        for base in bases: # look only in direct ancestors
            basemember = getattr(base, name, None)
            if not basemember:
                continue
            basememberdoc = getattr(basemember, "__doc__", None)
            if basememberdoc:
                member.__doc__ = basememberdoc
    return metaclass(classname, bases, dict)


# String constants for Exceptions / Warnings:
_ERRsizeChanged = "StableDict changed size during iteration!"
_WRNnoOrderArg  = "StableDict created/updated from unordered mapping object"
_WRNnoOrderKW   = \
              "StableDict created/updated with (unordered!) keyword arguments"


# Note: This class won't work with Python 3000 because the dict
#       protocol will change according to PEP3106. (However porting it
#       to Python 3.X will not be much of an effort.)
class StableDict(dict):
    """Dictionary remembering insertion order

    Order of item assignment is preserved and repeated when iterating
    over an instance.

    CAVEAT: When handing an unordered dict to either the constructor
    or the update() method the resulting order is obviously
    undefined. The same applies when initializing or updating with
    keyword arguments; i.e. keyword argument order is not preserved. A
    runtime warning will be issued in these cases via the
    warnings.warn function."""

    __metaclass__ = copy_baseclass_docs # copy docstrings from base class

    # Python 2.2 does not mangle __* inside __slots__
    __slots__ = ("_StableDict__ksl",) # key sequence list aka __ksl

    # @staticmethod
    def is_ordered(dictInstance):
        """Returns true if argument is known to be ordered."""
        if isinstance(dictInstance, StableDict):
            return True
        try: # len() may raise an exception.
            if len(dictInstance) <= 1:
                return True # A length <= 1 implies ordering.
        except:
            pass
        return False # Assume dictInstance.keys() is _not_ ordered.
    is_ordered = staticmethod(is_ordered)

    def __init__(self, *arg, **kw):
        if arg:
            if len(arg) > 1:
                raise TypeError("at most one argument permitted")
            arg = arg[0]
            if hasattr(arg, "keys"):
                if not self.is_ordered(arg):
                    _warn(_WRNnoOrderArg, RuntimeWarning, stacklevel=2)
                super(StableDict, self).__init__(arg, **kw)
                self.__ksl = list(arg.keys())
            else: # Must be a sequence of 2-tuples.
                super(StableDict, self).__init__(**kw)
                self.__ksl = []
                for pair in arg:
                    if len(pair) != 2:
                        raise ValueError("not a 2-tuple", pair)
                    self.__setitem__(pair[0], pair[1])
                if kw:
                    ksl = self.__ksl
                    for k in super(StableDict, self).iterkeys():
                        if k not in ksl:
                            ksl.append(k)
                    self.__ksl = ksl
        else: # No positional argument given.
            super(StableDict, self).__init__(**kw)
            self.__ksl = list(super(StableDict, self).keys())
        if len(kw) > 1:
            # There have been additionial keyword arguments.
            # Since Python passes them in an (unordered) dict
            # we cannot possibly preserve their order (without
            # inspecting the source or byte code of the call).
            _warn(_WRNnoOrderKW, RuntimeWarning, stacklevel=2)

    def update(self, *arg, **kw):
        if arg:
            if len(arg) > 1:
                raise TypeError("at most one non-keyword argument permitted")
            arg = arg[0]
            if hasattr(arg, "keys"):
                if not self.is_ordered(arg):
                    _warn(_WRNnoOrderArg, RuntimeWarning, stacklevel=2)
                super(StableDict, self).update(arg)
                ksl = self.__ksl
                for k in arg.keys():
                    if k not in ksl:
                        ksl.append(k)
                self.__ksl = ksl
            else: # Must be a sequence of 2-tuples.
                 for pair in arg:
                    if len(pair) != 2:
                        raise ValueError("not a 2-tuple", pair)
                    self.__setitem__(pair[0], pair[1])
        if kw:
            # There have been additionial keyword arguments.
            # Since Python passes them in an (unordered) dict
            # we cannot possibly preserve their order (without
            # inspecting the source or byte code of the call).
            if len(kw) > 1:
                _warn(_WRNnoOrderKW, RuntimeWarning, stacklevel=2)
            super(StableDict, self).update(kw)
            ksl = self.__ksl
            for k in kw.iterkeys():
                if k not in ksl:
                    ksl.append(k)
            self.__ksl = ksl

    def __str__(self):
        def _repr(x):
            if x is self:
                return "StableDict({...})" # Avoid unbounded recursion.
            return repr(x)
        return ( "StableDict({" + ", ".join([
                 "%r: %s" % (k, _repr(v)) for k, v in self.iteritems()])
                 + "})" )

    # Try to achieve: self == eval(repr(self))
    def __repr__(self):
        def _repr(x):
            if x is self:
                return "StableDict({...})" # Avoid unbounded recursion.
            return repr(x)
        return ( "StableDict([" + ", ".join([
                 "(%r, %s)" % (k, _repr(v)) for k, v in self.iteritems()])
                 + "])" )

    def __setitem__(self, key, value):
        super(StableDict, self).__setitem__(key, value)
        if key not in self.__ksl:
            self.__ksl.append(key)

    def __delitem__(self, key):
        if key in self.__ksl:
            self.__ksl.remove(key)
        super(StableDict, self).__delitem__(key)

    def __iter__(self):
        length = len(self)
        for key in self.__ksl[:]:
            yield key
        if length != len(self):
            raise RuntimeError(_ERRsizeChanged)

    def keys(self):
        return self.__ksl[:]

    def iterkeys(self):
        return self.__iter__()

    def values(self):
        return [ self[k] for k in self.__ksl ]

    def itervalues(self):
        length = len(self)
        for key in self.__ksl[:]:
            yield self[key]
        if length != len(self):
            raise RuntimeError(_ERRsizeChanged)

    def items(self):
        return [ (k, self[k]) for k in self.__ksl ]

    def iteritems(self):
        length = len(self)
        for key in self.__ksl[:]:
            yield ( key, self[key] )
        if length != len(self):
            raise RuntimeError(_ERRsizeChanged)

    def clear(self):
        super(StableDict, self).clear()
        self.__ksl = []

    def copy(self):
        return StableDict(self)

    def pop(self, k, *default):
        if k in self.__ksl:
            self.__ksl.remove(k)
        return super(StableDict, self).pop(k, *default)

    def popitem(self):
        item = super(StableDict, self).popitem()
        try:
            self.__ksl.remove(item[0])
        except:
            raise ValueError("cannot remove", item, self.__ksl, self)
        return item

# Our metaclass function became a method.  Make it a function again.
StableDict.__metaclass__ = staticmethod(copy_baseclass_docs)

# Given a sequence, return a sequence with unique items with order intact
def unique_union(seq):
    seen = {}
    result = []
    for item in seq:
        if item in seen: continue
        seen[item] = 1
        result.append(item)
    return result

# Given 2 sequences return the intersection with order intact as much as
# possible 
def intersect_order(s1, s2):
    seen = {}
    result = []
    for item in s1:
        if item in seen: continue
        seen[item] = 1
    for item in s2:
        if item not in seen: continue
        result.append(item)
    return result
