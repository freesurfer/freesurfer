# Original author - Martin Reuter
# $Id: qdectable_parser.py,v 1.1 2011/03/01 00:11:25 mreuter Exp $
import os
import logging
import sys
from datastruct_utils import *

class BadFileError(Exception):
    def __init__(self, filename):
        self.filename = filename
       
    def __str__(self):
        return self.filename

"""
This is the class which parses the qdectables files. 
"""
class QdecTableParser:

    fp = None # filepointer
    subjects_tp_map = StableDict()
    variables = []
    subjectsdir = ""

    # constructor just needs a filename
    # load it and report error
    def __init__(self, filename):
        self.filename = filename
        # raise exception if file doesn't exist
        if not os.path.exists(filename):
            raise BadFileError(filename)
        self.fp = open(filename, 'r')
                
    # we read in the file
    def parse(self):
        self.subjects_tp_map = StableDict()
        for line in self.fp:
            # valid line

            if line.rfind('#') == -1:
                strlst = line.split()
                #print strlst[0].upper()
                if strlst[0].upper() == 'SUBJECTS_DIR':
                    self.subjectsdir = strlst[1]
                    
                elif strlst[0].upper() == 'FSID':
                    if not strlst[1].upper() == 'FSID-BASE':
                        print 'ERROR: Make sure second column is \'fsid-base\' to specify the subject tempate (base)'
                        sys.exit(1)
                    self.variables= strlst[2:]  # 0 is tpid, 1 is templateid
                    
                else:
                    # check if time point already exists:
                    key = strlst[1]
                    tp  = strlst[0]
                    if key in self.subjects_tp_map:
                        for tpdata in self.subjects_tp_map[key]:
                            if tpdata[0] == tp:
                                print 'ERROR: Multiple occurence of \''+tp+'\' in '+key+'!'
                                sys.exit(1)
                    else:
                        self.subjects_tp_map[key] = []
                    # append tp and data to this subject
                    self.subjects_tp_map[key].append( [tp] + strlst[2:]  )
                
                

        # return
        return self.subjects_tp_map, self.variables, self.subjectsdir

