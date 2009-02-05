# Original author - Krish Subramaniam
# $Id: datastruct_utils.py,v 1.1 2009/02/05 21:49:02 krish Exp $

__all__ = ['Ddict', 'TableWriter']

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
        if not self.has_key(key):
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
            fp.write(self.pretext + c + self.posttext)
            for r in self.rows:
                fp.write(self.delimiter + '%g' %self.table[r][c])
            fp.write('\n')
        fp.close()    
        
