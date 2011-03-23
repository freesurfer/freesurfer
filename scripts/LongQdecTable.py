# Original author - Martin Reuter
# $Id: LongQdecTable.py,v 1.1.2.2 2011/03/23 16:03:47 mreuter Exp $
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
This is the class representing the long qdectables files. 
"""
class LongQdecTable:

    subjects_tp_map = StableDict()
    variables = []
    subjectsdir = ""
    filename = ""
    commonval = ""
    
    # constructor
    def __init__(self, stpmap, vari=None, sdir=None,cval=None):
        if isinstance(stpmap,str):
           self.parse(stpmap)
        else:
           self.subjects_tp_map = stpmap
           self.variables = vari
           self.subjectsdir = sdir
           self.commonval = cval
    
                        
    # we read in the file
    def parse(self,filename):
        self.filename = filename
        # raise exception if file doesn't exist
        if not os.path.exists(filename):
            raise BadFileError(filename)
        fp = open(filename, 'r')
    
        self.subjects_tp_map = StableDict()
        self.variables = []
        self.subjectsdir = ""
        
        for line in fp:
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

        fp.close()
        return self.subjects_tp_map, self.variables, self.subjectsdir


    def split(self,col='fsid-base'):
        # Will split the table based on the colum (eg. gender into male and female)
        # Default is the subject name (fsid-base)
        # NOTE: col values need to be identical across all time points for a given subject!
        # Returns a list of LongQdecTables
        
        alltables = []
        
        if col == 'fsid':
            print 'ERROR: cannot split fsid (one timepoint per file)?'
            sys.exit(1)        
        
        if col == 'fsid-base':
            for key,value in self.subjects_tp_map.items():
                stpmap = StableDict()
                stpmap[key] = value
                alltables.append(LongQdecTable(stpmap,self.variables,"",key))        

        elif col in self.variables:
#            for key,value in self.subject_tp_map:
            print 'Sorry, not implemented yet!\n'
            sys.exit(1)
            
        else:
            print 'ERROR: column "'+col+'" unknown!\n'
            sys.exit(1)
        
        return alltables


    def cross(self):
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
                    


    def sort(self,col):
        # sort table within each subject according to col (usually time) variable
        
        #first find column number:
        colnum = ''
        if col=='tpid':
            colnum = 0
        else:
            poscols = [i for i,x in enumerate(self.variables) if x == col]
            if len(poscols) != 1:
               print 'ERROR: did not find '+col+' or found it in several columns!'
               sys.exit(1)
            colnum = poscols[0]   
        
        for key,value in self.subjects_tp_map.items():
           #print 'Key before: '+key+'  ->  '+str(value)+'\n'
           #a = sorted(value, key=lambda tpdata: tpdata[colnum])
           #print 'Key after : '+key+'  ->  '+str(a)+'\n'
           
           self.subjects_tp_map[key] = sorted(value, key=lambda tpdata: tpdata[colnum])
           

    def write(self,filename):
        self.filename = filename
        fp = open(filename, 'w')
        if self.subjectsdir != "":
            fp.write('subjects_dir '+self.subjectsdir+'\n')
        fp.write('fsid fsid-base '+" ".join(self.variables)+'\n')
        for key,value in self.subjects_tp_map.items():
            for tpdata in value:
                fp.write(tpdata[0]+' '+key+' '+ ' '.join(tpdata[1:])+'\n')
        fp.close()
