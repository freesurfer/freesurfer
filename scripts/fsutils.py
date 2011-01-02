# Original author - Krish Subramaniam
# $Id: fsutils.py,v 1.3 2011/01/02 05:24:00 krish Exp $
import os
import logging
from misc import *
from subject_info import *
from datastruct_utils import *

# logging for both stats2table

ch = logging.StreamHandler()
#create logger
aseglogger = logging.getLogger("asegstats2table")
aseglogger.setLevel(logging.INFO)
aseglogger.addHandler(ch)

aparclogger = logging.getLogger("aparcstats2table")
aparclogger.setLevel(logging.INFO)
aparclogger.addHandler(ch)

class BadFileError(Exception):
    def __init__(self, filename):
        self.filename = filename
       
    def __str__(self):
        return self.filename

"""
This is the base class which parses the .stats files. 
Both aparcstats and asegstats parsers should derive from this
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

    measure_column_map = {'volume':3, 'mean':5, 'std':6}
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
            if measure == 'volume':
                beg_struct_tuple = (
                        ('# Measure lhCortex, lhCortexVol,', 'lhCortexVol'), 
                        ('# Measure rhCortex, rhCortexVol,', 'rhCortexVol'), 
                        ('# Measure Cortex, CortexVol,', 'CortexVol'),
                        ('# Measure lhCorticalWhiteMatter, lhCorticalWhiteMatterVol,','lhCorticalWhiteMatterVol'),
                        ('# Measure rhCorticalWhiteMatter, rhCorticalWhiteMatterVol,','rhCorticalWhiteMatterVol'),
                        ('# Measure SubCortGray, SubCortGrayVol,','SubCortGrayVol'),
                        ('# Measure TotalGray, TotalGrayVol,','TotalGrayVol'),
                        ('# Measure SuperTentorial, SuperTentorialVol,','SuperTentorialVol'),
                        ('# Measure IntraCranialVol, ICV,','IntraCranialVol'),
                        ('# Measure BrainSeg, BrainSegVol,','BrainSegVol'),)
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
    measure_column_map = {'area':2, 'volume':3, 'thickness':4, 'thicknessstd':5, 'meancurv':6 }
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

        return self.parc_measure_map
