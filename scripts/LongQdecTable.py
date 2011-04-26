# Original author - Martin Reuter
# $Id: LongQdecTable.py,v 1.2 2011/04/26 17:00:22 mreuter Exp $
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
    cross = True
    
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
        self.cross = True
        
        gotheaders = False
        for line in fp:
            # valid line

            if line.rfind('#') == -1:
                strlst = line.split()
                #print strlst[0].upper()
                if strlst[0].upper() == 'SUBJECTS_DIR':
                    self.subjectsdir = strlst[1]
                    
                elif strlst[0].upper() == 'FSID':
                    gotheaders = True
                    if not strlst[1].upper() == 'FSID-BASE':
                        print '\nWarning: second column is not \'fsid-base\' assuming cross sectional qdec table\n'
                        #print '\nERROR: Make sure second column is \'fsid-base\' to specify the subject tempate (base)\n'
                        #sys.exit(1)
                    else:
                        self.cross = False
                    self.variables= strlst[2:]  # 0 is tpid, 1 is templateid
                    
                else:
                    if not gotheaders:
                        print '\nERROR: qdec table missing correct column headers?'
                        print '       Make sure first column is labeled \'fsid\' for the time point and'
                        print '       second column is \'fsid-base\' to specify the subject tempate (base), e.g.:\n'
                        print ' fsid    fsid-base   age '
                        print ' me1     me          22.3 '
                        print ' me2     me          23.2 '
                        print ' you1    you         21.6 '
                        print ' you2    you         22.5\n'                
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
                            print '\nERROR: no fsid-base in header, but fsid '+key+' seems to exists multiple times?\n'
                            sys.exit(1)
                        # append this time point to this base
                        for tpdata in self.subjects_tp_map[key]:
                            if tpdata[0] == tp:
                                print 'ERROR: Multiple occurence of time point (fsid) \''+tp+'\' in (fsid-base) '+key+'!'
                                sys.exit(1)
                    else:
                        self.subjects_tp_map[key] = []
                    # append tp and data to this subject
                    if self.cross:
                        self.subjects_tp_map[key].append( [tp] + strlst[1:]  )
                    else:
                        self.subjects_tp_map[key].append( [tp] + strlst[2:]  )

        fp.close()
        return self.subjects_tp_map, self.variables, self.subjectsdir, self.cross


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
