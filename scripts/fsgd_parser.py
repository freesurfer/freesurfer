# Original author - Martin Reuter
# $Id: fsgd_parser.py,v 1.2 2011/02/17 23:10:39 mreuter Exp $
import os
import logging
import sys
#from misc import *
#from subject_info import *
from datastruct_utils import *

class BadFileError(Exception):
    def __init__(self, filename):
        self.filename = filename
       
    def __str__(self):
        return self.filename

"""
This is the class which parses the fsgd files. 
"""
class FsgdParser:

    fp = None # filepointer
    subjects_tp_map = StableDict()
    variables = []
    defaultvar = ""

    # constructor just needs a fsgd filename
    # load it and report error
    def __init__(self, filename):
        self.filename = filename
        # raise exception if file doesn't exist
        if not os.path.exists(filename):
            raise BadFileError(filename)
        self.fp = open(filename, 'r')
                
    # we read in the fsgd
    def parse(self):
        self.subjects_tp_map = StableDict()
        for line in self.fp:
            # valid line

            if line.rfind('#') == -1:
                strlst = line.split()
                #print strlst[0].upper()
                if strlst[0].upper() == 'CLASS':
                    self.subjects_tp_map[strlst[1]] = []

                if strlst[0].upper() == 'VARIABLES':
                    self.variables= strlst[1:]
                   
                if strlst[0].upper() == 'INPUT':
                    # check if class exists
                    key = strlst[2]
                    tp  = strlst[1]
                    if not key in self.subjects_tp_map:
                        print 'ERROR: Subject \''+key+'\' was not defined as "Class"!'
                        sys.exit(1)
                    # check if time point already exists:
                    for tpdata in self.subjects_tp_map[key]:
                        if tpdata[0] == tp:
                            print 'ERROR: Multiple occurence of \''+tp+'\' in '+key+'!'
                            sys.exit(1)
                    # append tp and data to this subject
                    self.subjects_tp_map[key].append([ tp ] + strlst[3:]  )
                    
                if strlst[0].upper() == 'DEFAULTVARIABLE':
                    self.defaultvar = strlst[1]
                    
        if not self.defaultvar == '':
            index=-1
            for index in (i for i in xrange(len(self.variables)) if self.variables[i].upper()==self.defaultvar.upper()):
                break
            if index == -1 or not self.variables[index].upper()==self.defaultvar.upper():
                print 'ERROR: DefaultVariable \''+str(self.defaultvar)+'\' not found in Variables: '+str(self.variables)
                sys.exit(1)



        # return
        return self.subjects_tp_map, self.variables, self.defaultvar

