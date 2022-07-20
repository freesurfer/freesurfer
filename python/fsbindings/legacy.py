# Long ago, before the days of freesurfer 7, python files in the scripts
# directory imported some utilities from a fsutils.py file. Now that we
# have a legitimate python library, all the tools from fsutils have been
# copied here. Those particular python scripts, like aparcstats2table, will
# now import 'fsuilts' via:
#
# import fsbindings.legacy as fsutils
#
# Ideally, everything in this file should either be phased out or integrated
# cleanly with the rest of the freesurfer package, but for now this is an easy
# way to get the old scripts running with fspython without having to modify much.

import os
import sys
import logging
import itertools


# -----------------------------------------------------------------------------
#                                   MISC
# -----------------------------------------------------------------------------


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


def check_subjdirs():
    """
    Quit if SUBJECTS_DIR is not defined as an environment variable. This is not
    a function which returns a boolean. Execution is stopped if not found.
    If found, returns the SUBJECTS_DIR
    """
    if 'SUBJECTS_DIR' not in os.environ:
        print('ERROR: SUBJECTS_DIR environment variable not defined!')
        sys.exit(1)
    return os.environ['SUBJECTS_DIR']


# -----------------------------------------------------------------------------
#                              DATA STRUCT UTILS
# -----------------------------------------------------------------------------


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
_WRNnoOrderKW   = "StableDict created/updated with (unordered!) keyword arguments"


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


# -----------------------------------------------------------------------------
#                             LOGGING AND PARSING
# -----------------------------------------------------------------------------


ch = logging.StreamHandler()
#create logger
aseglogger = logging.getLogger("asegstats2table")
aseglogger.setLevel(logging.INFO)
aseglogger.addHandler(ch)

aparclogger = logging.getLogger("aparcstats2table")
aparclogger.setLevel(logging.INFO)
aparclogger.addHandler(ch)

tractlogger = logging.getLogger("tractstats2table")
tractlogger.setLevel(logging.INFO)
tractlogger.addHandler(ch)

overlaplogger = logging.getLogger("overlap2table")
overlaplogger.setLevel(logging.INFO)
overlaplogger.addHandler(ch)

class BadFileError(Exception):
    def __init__(self, filename):
        self.filename = filename
       
    def __str__(self):
        return self.filename

"""
This is the base class which parses the .stats files. 
"""
class StatsParser:

    measure_column_map = {} #to be defined in subclass
    fp = None # filepointer
    include_structlist = None #parse only these structs
    exclude_structlist = None #exclude these structs

    structlist = None #list of structures
    measurelist = None #list of measures corresponding to the structures
    # len ( structlist ) must be equal to len (measurelist )

    # constructor just needs a .stats filename
    # load it and report error
    def __init__(self, filename):
        self.filename = filename
        # raise exception if file doesn't exist or
        # is too small to be a valid stats file
        if not os.path.exists(filename):
            raise BadFileError(filename)
        if os.path.getsize(filename) < 10:
            raise BadFileError(filename)
        self.fp = open(filename, 'r')

        self.include_structlist = [] 
        self.exclude_structlist = [] 
        self.structlist = []
        self.measurelist = []
        self.include_vol_extras = 1;

    # parse only the following structures
    def parse_only(self, structlist):
        # this is simply an Ordered Set
        # we need this because if inputs repeat, this will take care
        self.include_structlist = StableDict()
        for struct in structlist:
            self.include_structlist[struct] = 1
        self.structlist = []
        self.measurelist = []
        
    # exclude the following structures
    def exclude_structs(self, structlist):
        # this is simply an Ordered Set
        # we need this because if inputs repeat, this will take care
        self.exclude_structlist = StableDict()
        for struct in structlist:
            self.exclude_structlist[struct] = 1
        self.structlist = []
        self.measurelist = []

        
    # actual parsing will be done by subclass
    def parse(self):
        pass
    
"""
aseg.stats parser ( or parser for similarly formatted .stats files )
Derived from StatsParser
"""
class AsegStatsParser(StatsParser):

    measure_column_map = {'nvoxels':2,'nvertices':2,'volume':3,'Area_mm2':3, 'mean':5, 'std':6, 'max':8, 'snr':10}
    maxsegno = None
    id_name_map = StableDict()

    # segs > max_seg_no are discarded
    def set_maxsegno(self, maxsegno):
        self.maxsegno = int(maxsegno)

    # we take in the measure we need..
    def parse(self, measure):
        self.id_name_map = StableDict()
        for line in self.fp:
            # valid line
            if line.rfind('#') == -1:
                strlst = line.split()
                seg = strlst[1] # segmentation id is in the 2nd field
                if self.maxsegno and int(seg) > self.maxsegno:
                    continue
                if self.exclude_structlist and seg in self.exclude_structlist.keys():
                    continue
                self.id_name_map[seg] = strlst[4] # Col 4 is Seg Name
                # Float of the measure we are interested in 
                self.measurelist.append( float(strlst[self.measure_column_map[measure]] ))

        # if we have a spec which instructs the table to have only specified segs,
        # we need to make sure the order also has to be maintained
        # hence we create temp dicts and copy them to the return dicts
        # If the segmentation name is not present in the current statsfile,
        # it is simply named as 'Placeholder_Segmentation' in the hope that we'd bump
        # into this in the future 
        if self.include_structlist: 
            tmp_id_name_map = StableDict()
            tmp_ml = []
            for oseg in self.include_structlist.keys():
                idlist = self.id_name_map.keys()
                if oseg in idlist:
                    corres_meas_index = idlist.index( oseg )
                    tmp_id_name_map[oseg] = self.id_name_map[oseg] 
                    tmp_ml.append( self.measurelist[corres_meas_index])
                else:
                    tmp_id_name_map[oseg] ='Placeholder_Segmentation'
                    tmp_ml.append( 0.0 )
            self.id_name_map = tmp_id_name_map
            self.measurelist = tmp_ml[:]

        # edge case in case of measure = volume
        self.fp.seek(0)
        for line in self.fp:
            # volume measures.. which are found at the beginning of files and the 
            # values are in the 4th column ( comma separated )
            # segidlist is not useful here, all the values of segidlist will be 
            # 'Placeholder_ID'
            # this is because len of all 3 lists -- structlist, measurelist and segidlists should be same
            if measure == 'volume' and self.include_vol_extras:
                beg_struct_tuple = (
                        ('# Measure lhCortex, lhCortexVol,', 'lhCortexVol'), 
                        ('# Measure rhCortex, rhCortexVol,', 'rhCortexVol'), 
                        ('# Measure Cortex, CortexVol,', 'CortexVol'),
                        ('# Measure lhCorticalWhiteMatter, lhCorticalWhiteMatterVol,','lhCorticalWhiteMatterVol'),
                        ('# Measure rhCorticalWhiteMatter, rhCorticalWhiteMatterVol,','rhCorticalWhiteMatterVol'),
                        ('# Measure CorticalWhiteMatter, CorticalWhiteMatterVol,','CorticalWhiteMatterVol'),
                        ('# Measure lhCerebralWhiteMatter, lhCerebralWhiteMatterVol,','lhCerebralWhiteMatterVol'),
                        ('# Measure rhCerebralWhiteMatter, rhCerebralWhiteMatterVol,','rhCerebralWhiteMatterVol'),
                        ('# Measure CerebralWhiteMatter, CerebralWhiteMatterVol,','CerebralWhiteMatterVol'),
                        ('# Measure SubCortGray, SubCortGrayVol,','SubCortGrayVol'),
                        ('# Measure TotalGray, TotalGrayVol,','TotalGrayVol'),
                        ('# Measure SuperTentorial, SuperTentorialVol,','SuperTentorialVol'),
                        ('# Measure SupraTentorial, SupraTentorialVol,','SupraTentorialVol'),
                        ('# Measure SupraTentorialNotVent, SupraTentorialVolNotVent,','SupraTentorialVolNotVent'),
                        ('# Measure SupraTentorialNotVentVox, SupraTentorialVolNotVentVox,','SupraTentorialVolNotVentVox'),
                        ('# Measure IntraCranialVol, ICV,','IntraCranialVol'),
                        ('# Measure EstimatedTotalIntraCranialVol, eTIV,','EstimatedTotalIntraCranialVol'),
                        ('# Measure Mask, MaskVol,','MaskVol'),
                        ('# Measure BrainVol-to-eTIV, BrainVol-to-eTIV,','BrainVol-to-eTIV'),
                        ('# Measure BrainSegVol-to-eTIV, BrainSegVol-to-eTIV,','BrainSegVol-to-eTIV'),
                        ('# Measure MaskVol-to-eTIV, MaskVol-to-eTIV,','MaskVol-to-eTIV'),
                        ('# Measure lhSurfaceHoles, lhSurfaceHoles,','lhSurfaceHoles'),
                        ('# Measure rhSurfaceHoles, rhSurfaceHoles,','rhSurfaceHoles'),
                        ('# Measure SurfaceHoles, SurfaceHoles,','SurfaceHoles'),
                        ('# Measure BrainSeg, BrainSegVol,','BrainSegVol'),
                        ('# Measure BrainSegNotVent, BrainSegVolNotVent,','BrainSegVolNotVent'),
                        ('# Measure BrainSegNotVentSurf, BrainSegVolNotVentSurf,','BrainSegVolNotVentSurf'),)
                c = 0
                for start, structn in beg_struct_tuple:
                    c = c + 1
                    if line.startswith(start) :
                        strlst = line.split(',')
                        self.id_name_map['Placeholder_ID'+str(c)] = structn 
                        self.measurelist.append( float(strlst[3]))

        # check sanity
        assert( len(self.id_name_map) == len(self.measurelist))
        # return
        return self.id_name_map, self.measurelist

                  
"""
?h.aparc*.stats parser ( or parser for similarly formatted .stats files )
Derived from StatsParser
"""
class AparcStatsParser(StatsParser):

    # this is a map of measure requested and its corresponding column# in ?h.aparc*.stats
    measure_column_map = {'area':2, 'volume':3, 'thickness':4, 'thickness.T1':4, 'thicknessstd':5, 'meancurv':6, 'gauscurv':7, 'foldind':8, 'curvind':9 }
    parc_measure_map = StableDict()

    # we take in the measure we need..
    def parse(self, measure):
        self.parc_measure_map = StableDict()
        for line in self.fp:
            # a valid line is a line without a '#'
            if line.rfind('#') == -1:
                strlist = line.split()
                # for every parcellation
                parcid = strlist[0]
                val = float(strlist[self.measure_column_map[measure]])
                self.parc_measure_map[parcid] = val

        # if we have a spec which instructs the table to have only specified parcs,
        # we need to make sure the order has to be maintained
        if self.include_structlist: 
            tmp_parc_measure_map = StableDict()
            for oparc in self.include_structlist.keys():
                parclist = self.parc_measure_map.keys()
                if oparc in parclist:
                    tmp_parc_measure_map[oparc] = self.parc_measure_map[oparc] 
                else:
                    tmp_parc_measure_map[oparc] = 0.0
            self.parc_measure_map = tmp_parc_measure_map

        # measures which are found at the beginning of files. 
        self.fp.seek(0)
        for line in self.fp:
            beg_struct_tuple = (('# Measure EstimatedTotalIntraCranialVol, eTIV', 'eTIV'),)
            for start, structn in beg_struct_tuple:
                if line.startswith(start):
                    strlst = line.split(',')
                    self.parc_measure_map[structn] = float( strlst[3])
            beg_struct_tuple = (('# Measure BrainSegNotVent, BrainSegVolNotVent', 'BrainSegVolNotVent'),)
            for start, structn in beg_struct_tuple:
                if line.startswith(start):
                    strlst = line.split(',')
                    self.parc_measure_map[structn] = float( strlst[3])
            if measure == 'area':
                beg_struct_tuple = (
                    ('# Measure Cortex, WhiteSurfArea,', 'WhiteSurfArea'),
                )
                for start, structn in beg_struct_tuple:
                    if line.startswith(start):
                        strlst = line.split(',')
                        self.parc_measure_map[structn] = float( strlst[3])
            if measure == 'thickness':
                beg_struct_tuple = (
                    ('# Measure Cortex, MeanThickness,', 'MeanThickness'),
                )
                for start, structn in beg_struct_tuple:
                    if line.startswith(start):
                        strlst = line.split(',')
                        self.parc_measure_map[structn] = float( strlst[3])

        return self.parc_measure_map


"""
tracula overall.stats file parser ( or parser for similarly formatted .stats files )
Derived from StatsParser
"""
class TractOverallStatsParser(StatsParser):
    # this is a map of measure and its corresponding value in tract.overall.stats
    measure_value_map = StableDict()

    # we take in the measure we need..
    def parse(self):
        self.measure_value_map = StableDict()
        pre_tuple = ( ('# pathwayname','pathway'),
                      ('# subjectname','subject'),
                    )
        pthwy_subj = StableDict()
        for line in self.fp:
            # a valid line is a line without a '#'
            if line.rfind('#') == -1:
                strlist = line.split()
                # for every measure assign value 
                meas = strlist[0]
                # if measure is in excluded measure list, we are not interested
                if self.exclude_structlist and meas in self.exclude_structlist.keys():
                    continue
                val = float(strlist[1])
                self.measure_value_map[meas] = val
            else:
                # lookfor subjectname and pathwayname
                for start, what in pre_tuple:
                    if line.startswith(start):
                        strlst = line.split(' ')
                        name = strlst[2].strip()
                        pthwy_subj[what] = name

        # if we have a spec which instructs the table to have only specified measures,
        # we need to make sure the order has to be maintained
        if self.include_structlist: 
            tmp_measure_value_map = StableDict()
            for omeasure in self.include_structlist.keys():
                measurelist = self.measure_value_map.keys()
                if omeasure in measurelist:
                    tmp_measure_value_map[omeasure] = self.measure_value_map[omeasure] 
                else:
                    print('ERROR: The measure you requested: '+ omeasure + ' is not valid')
                    print('It is not found one of the overall statistics file')
                    sys.exit(1)
            self.measure_value_map = tmp_measure_value_map
        
        return pthwy_subj, self.measure_value_map

 
"""
tracula byvoxel.stats file parser ( or parser for similarly formatted .stats files )
Derived from StatsParser
"""
class TractByvoxelStatsParser(StatsParser):
    # this is a map of measure and its corresponding value in tract.byvoxel.stats
    measure_value_map = StableDict()

    # map of measures needed and their corresponding column in byvoxel.stats file
    measure_column_map = {'AD':3, 'RD':4, 'MD':5, 'FA':6 }

    # we take in the measure we need..
    def parse(self, measure):
        self.measure_value_map = StableDict()
        pre_tuple = ( ('# pathwayname','pathway'),
                      ('# subjectname','subject'),
                    )
        pthwy_subj = StableDict()

        # scroll to the start of data. find subjectname
        # and pathwayname meanwhile
        line = self.fp.readline()
        while not line.startswith('# pathway start'):
            # lookfor subjectname and pathwayname
            for start, what in pre_tuple:
                if line.startswith(start):
                    strlst = line.split(' ')
                    name = strlst[2].strip()
                    pthwy_subj[what] = name
            line = self.fp.readline()

        line = self.fp.readline() # skip to the 'x y z AD RD MD FA line
        line = self.fp.readline() # skip to the actual data
        
        # scroll until the end of data collecting the spec-ed measure
        count = 0
        while not line.startswith('# pathway end'):
            strlst = line.split()
            val = float(strlst[ self.measure_column_map[measure] ])
            count = count + 1
            self.measure_value_map[str(count)] = val
            line = self.fp.readline()

        return pthwy_subj, self.measure_value_map

"""
Overlap stats file parser ( for files created by mri_compute_overlap )
Derived from StatsParser
"""
class OverlapStatsParser(StatsParser):
    # this is a map of measure and its corresponding value in ...
    measure_value_map = StableDict()
    # map of measures needed and their corresponding column in byvoxel.stats file
    measure_column_map = {'voldiff':1, 'dice':2, 'jacc':3 }

    # we take in the measure we need..
    def parse(self, measure):
        self.measure_value_map = StableDict()
        for line in self.fp:
            # a valid line is a line without a '#'
            if line.rfind('#') == -1:
                strlist = line.split()
                # for every label
                labelid = strlist[0]
                val = float(strlist[self.measure_column_map[measure]])
                self.measure_value_map[labelid] = val
        return self.measure_value_map


# -----------------------------------------------------------------------------
#                                    LONG
# -----------------------------------------------------------------------------


"""
This is the class representing the long qdectables files. 
"""
class LongQdecTable:

    subjects_tp_map = StableDict()
    variables = []
    subjectsdir = ""
    filename = ""
    commonval = ""
    cross = True
    
    # constructor
    def __init__(self, stpmap=None, vari=None, sdir=None,cval=None,cross=True):
        if stpmap == None:
           return
        if isinstance(stpmap,str):
           self.parse(stpmap)
        else:
           self.subjects_tp_map = stpmap
           self.variables = vari
           self.subjectsdir = sdir
           self.commonval = cval
           self.cross = cross
    
    # append a new subject (base)
    def append(self,bid,alltpdata,varlist):
        if self.variables != varlist:
            print('\nERROR: append: variables do not agree\n')
            sys.exit(1)        
        if bid in self.subjects_tp_map:
            print('\nERROR: append: subject '+bid+' seems to exists already?\n')
            sys.exit(1)        
        self.subjects_tp_map[bid] = alltpdata
        #later maybe if base exists, append tpdata?: self.subjects_tp_map[bid].append( alltpdata  )
        
                        
    # we read in the file
    def parse(self,filename,warncross=True):
        self.filename = filename
        # raise exception if file doesn't exist
        if not os.path.exists(filename):
            raise BadFileError(filename)
        fp = open(filename, 'r')
    
        self.subjects_tp_map = StableDict()
        self.variables = []
        self.subjectsdir = ""
        self.cross = True
        
        gotheaders = False
        for line in fp:
            # remove WS from begining and end:
            line = line.strip()
            # skip empty line
            if len(line) == 0:
                continue
            # skip comment
            if line[0] == '#':
                continue

            strlst = line.split()
            #print(strlst[0].upper())
            if strlst[0].upper() == 'SUBJECTS_DIR':
                self.subjectsdir = strlst[1]
                
            elif strlst[0].upper() == 'FSID':
                gotheaders = True
                if not strlst[1].upper() == 'FSID-BASE':
                    if warncross:
                        print('\nWarning: second column is not \'fsid-base\' assuming cross sectional qdec table\n')
                    #print('\nERROR: Make sure second column is \'fsid-base\' to specify the subject tempate (base)\n')
                    #sys.exit(1)
                    self.variables= strlst[1:]  # 0 is tpid, 1 is templateid
                else:
                    self.cross = False
                    self.variables= strlst[2:]  # 0 is tpid, 1 is templateid
            elif strlst[0].startswith('Measure:') or strlst[0].startswith('h.aparc',1,):
                print('Input is probably stats table, reading it as cross sectional...\n')
                self.cross = True
                self.variables = strlst[1:] # 0 is subject id
                gotheaders = True
                
            else:
                if not gotheaders:
                    print('\nERROR: qdec table missing correct column headers?')
                    print('       Make sure first column is labeled \'fsid\' for the time point and')
                    print('       second column is \'fsid-base\' to specify the subject tempate (base), e.g.:\n')
                    print(' fsid    fsid-base   age ')
                    print(' me1     me          22.3 ')
                    print(' me2     me          23.2 ')
                    print(' you1    you         21.6 ')
                    print(' you2    you         22.5\n'                )
                    sys.exit(1)
                    
                # check if time point already exists:
                #base is in second column:
                key = strlst[1]
                # if cross format, use fsid for base:
                if self.cross:
                    key = strlst[0]
                # fsid is the time point:    
                tp  = strlst[0]
                # if base exists
                if key in self.subjects_tp_map:
                    # if cross, make sure fsid does not have duplicates
                    if self.cross:
                        print('\nERROR: no fsid-base in header, but fsid '+key+' seems to exists multiple times?\n')
                        sys.exit(1)
                    # check if tp is already in this base
                    for tpdata in self.subjects_tp_map[key]:
                        if tpdata[0] == tp:
                            print('ERROR: Multiple occurence of time point (fsid) \''+tp+'\' in (fsid-base) '+key+'!')
                            sys.exit(1)
                else:
                    self.subjects_tp_map[key] = []
                # append tp and data to this subject
                if self.cross:
                    self.subjects_tp_map[key].append( [tp] + strlst[1:]  )
                else:
                    self.subjects_tp_map[key].append( [tp] + strlst[2:]  )
        fp.close()
        if len(self.subjects_tp_map) == 0:
            raise BadFileError(filename)        
        return self.subjects_tp_map, self.variables, self.subjectsdir, self.cross


    def split(self,col='fsid-base'):
        # Will split the table based on the colum (eg. gender into male and female)
        # Default is the subject name (fsid-base)
        # NOTE: col values need to be identical across all time points for a given subject!
        # Returns a list of LongQdecTables
        
        alltables = []
        
        if col == 'fsid':
            print('ERROR: cannot split fsid (one timepoint per file)?')
            sys.exit(1)        
        
        if col == 'fsid-base':
            for key,value in self.subjects_tp_map.items():
                stpmap = StableDict()
                stpmap[key] = value
                alltables.append(LongQdecTable(stpmap,self.variables,"",key,self.cross))        

        elif col in self.variables:
            allasdict = StableDict()
            poscols = [i for i,x in enumerate(self.variables) if x == col]
            if len(poscols) != 1:
                print('ERROR: did not find '+col+' or found it in several columns!')
                sys.exit(1)
            colnum = poscols[0] + 1
              
            for bid,value in self.subjects_tp_map.items():
                key = value[0][colnum]
                print('Key: '+str(key)+'\n')
                for tpdata in value:
                    if tpdata[colnum] != key:
                        print('ERROR: split: '+col+' value needs to be the same within each subject ('+bid+')!')
                        sys.exit(1)
                if key not in allasdict:
                    stpmap = StableDict()
                    stpmap[bid] = value
                    allasdict[key] = LongQdecTable(stpmap,self.variables,"",key,self.cross)
                else:
                    allasdict[key].append(bid,value,self.variables)
            
            for key in allasdict:
               alltables.append(allasdict[key]) 
            
        else:
            print('ERROR: column "'+col+'" unknown!\n')
            sys.exit(1)
        
        return alltables


    def make_cross(self):
        # collapsing long table to cross
        # numerical values will be averaged
        # for other values we take the first occurence (hopefully basline if the table was sorted)
        #subjects_new_map = StableDict()        
        for key,value in self.subjects_tp_map.items():
            subjentry = [ key ]
            for i in range(1,len(value[0])):
                try:
                    col = [float(entry[i]) for entry in value]
                    subjentry.append(str(sum(col)/len(col)))
                except ValueError:
                    subjentry.append(str(value[0][i]))
            #subjects_new_map[key] = [ subjentry ] 
            self.subjects_tp_map[key] = [ subjentry ] 
        self.cross = True       


    def sort(self,col):
        # sort table within each subject according to col (usually time) variable
        
        #first find column number:
        colnum = ''
        if col=='tpid':
            colnum = 0
        else:
            poscols = [i for i,x in enumerate(self.variables) if x == col]
            if len(poscols) != 1:
               print('ERROR: did not find '+col+' or found it in several columns!')
               sys.exit(1)
            colnum = poscols[0]   
        
        for key,value in self.subjects_tp_map.items():
            #print('Key before: '+key+'  ->  '+str(value)+'\n')
            #a = sorted(value, key=lambda tpdata: tpdata[colnum])
            #print('Key after : '+key+'  ->  '+str(a)+'\n')
           
            self.subjects_tp_map[key] = sorted(value, key=lambda tpdata: tpdata[colnum])
    
    def append_table(self,filename):
        if self.cross:
            print('ERROR: append_table not supported for type cross!')
            sys.exit(1)
        
        # append columns from another table (read it from disk) to this
        # it is assumed that a row exists for each subject.tp in this table
        print('Parsing the qdec table: '+filename)
        statstable = LongQdecTable()
        statstable.parse(filename,False) #don't warn about being cross sectional table
        #print(statstable.variables)
        
        self.variables = list(itertools.chain(*[self.variables, statstable.variables]))
        first = True
        crossnames = True
        #iterate through current table
        for subjectid, tplist in self.subjects_tp_map.items():
            
        
            if not statstable.cross:
                print('statstable is not corss (= long)\n')
                # table to append is in long format
                #  check if subject is here
                if subjectid not in statstable.subjects_tp_map:
                    print('ERROR: did not find '+subjectid+' in table '+filename+'!')
                    sys.exit(1)
                    
                # get that data
                addtplist = statstable.subjects_tp_map[subjectid]
                
                # check if all time points are in same order
                for i,tpdata,addtpdata in itertools.izip(itertools.count(),tplist,addtplist):
                    if tpdata[0] != addtpdata[0]:
                        print('ERROR: time point id'+tpdata[0]+' not found in other table!')
                        sys.exit(1)
                    # append all columns (except the id)
                    self.subjects_tp_map[subjectid][i] = list(itertools.chain(*[self.subjects_tp_map[subjectid][i], addtplist[i][1:]]))
            else:
                # if saved in cross format
                for i,tpdata in itertools.izip(itertools.count(),tplist):
                    #determin if fsid is cross or long (only once)
                    if first:
                        first=False
                        tpid = tpdata[0]
                        if tpid in statstable.subjects_tp_map:
                            crossnames = True
                        elif tpid+'.long.'+subjectid in statstable.subjects_tp_map:
                            crossnames = False
                        else:
                            print('ERROR: time point id'+tpid+' not found in other table!')
                            sys.exit(1)
                    # get the name
                    tpid = tpdata[0]
                    if not crossnames:
                        tpid = tpdata[0]+'.long.'+subjectid
                    # get the data
                    addtplist = statstable.subjects_tp_map[tpid]
                    # append all columns (except the id)
                    self.subjects_tp_map[subjectid][i] = list(itertools.chain(*[self.subjects_tp_map[subjectid][i], addtplist[0][1:]]))

    def write(self,filename):
        self.filename = filename
        fp = open(filename, 'w')
        if self.subjectsdir != "":
            fp.write('subjects_dir '+self.subjectsdir+'\n')
        if self.cross:
            fp.write('fsid '+" ".join(self.variables)+'\n')
        else:
            fp.write('fsid fsid-base '+" ".join(self.variables)+'\n')
        for key,value in self.subjects_tp_map.items():
            for tpdata in value:
                if self.cross:
                    fp.write(tpdata[0]+' '+ ' '.join(tpdata[1:])+'\n')
                else:
                    fp.write(tpdata[0]+' '+key+' '+ ' '.join(tpdata[1:])+'\n')
        fp.close()
