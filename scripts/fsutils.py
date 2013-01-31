# Original author - Krish Subramaniam
# $Id: fsutils.py,v 1.15 2013/01/31 19:22:45 greve Exp $
import os
import logging
import sys
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

tractlogger = logging.getLogger("tractstats2table")
tractlogger.setLevel(logging.INFO)
tractlogger.addHandler(ch)

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
            if measure == 'volume':
                beg_struct_tuple = (
                        ('# Measure lhCortex, lhCortexVol,', 'lhCortexVol'), 
                        ('# Measure rhCortex, rhCortexVol,', 'rhCortexVol'), 
                        ('# Measure Cortex, CortexVol,', 'CortexVol'),
                        ('# Measure lhCorticalWhiteMatter, lhCorticalWhiteMatterVol,','lhCorticalWhiteMatterVol'),
                        ('# Measure rhCorticalWhiteMatter, rhCorticalWhiteMatterVol,','rhCorticalWhiteMatterVol'),
                        ('# Measure CorticalWhiteMatter, CorticalWhiteMatterVol,','CorticalWhiteMatterVol'),
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
    measure_column_map = {'area':2, 'volume':3, 'thickness':4, 'thicknessstd':5, 'meancurv':6, 'gauscurv':7, 'foldind':8, 'curvind':9 }
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
                    print 'ERROR: The measure you requested: '+ omeasure + ' is not valid'
                    print 'It is not found one of the overall statistics file'
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
