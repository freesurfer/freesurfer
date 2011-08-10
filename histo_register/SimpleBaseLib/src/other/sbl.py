import sblc
import os


# set working dir to SBL's dataPath
# fix(later): should do this directly from C++ code
os.chdir( sblc.dataPath() )


# convert a python value to a string value suitable to for SBL commands/configs
def strProc( val ):
    valStr = str( val );
    if valStr == "False":
        valStr = "0"
    if valStr == "True":
        valStr = "1"
    return valStr


# represents an entry in a configuration file
class ConfigEntry:
   
    # basic constructor
    def __init__( self, _name, _value, _comment ):
        self.name = _name
        self.value = _value
        self.comment = _comment


# represents a configuration file
class Config:

    # basic constructor
    def __init__( self ):
#        print( "adding entries" )
#        self._entries = []
        self.__dict__[ "_entries" ] = []

    # add a config entry
    def __setattr__( self, name, value ):
        if not name.startswith( "_" ):
            found = False
            for e in self._entries: # fix(later): use dict (though want to maintain order)
                if e.name == name:
                    e.value = value 
                    found = True
            if not found:
                self._entries.append( ConfigEntry( name, value, "" ) )
#        elif name == "_entries":
#            print( "adding entries again" );
#            self.__dict__[ "_entries" ] = []
        

    # read a config entry
    def __getattr__( self, name ):
        if not name.startswith( "_" ):
            for e in self._entries: # fix(later): use dict (though want to maintain order)
                if e.name == name:
                    return e.value
        raise AttributeError

    # create a string version suitable for passing to an SBL command
    def __str__( self ):
        s = ""
        for e in self._entries:
            if e.name:
                s += e.name + "=" + strProc( e.value ) + " "
        return s

    # load a configuration file (in SBL format)
    def load( self, fileName ):
        f = open( fileName, "r" )
        if f:
            for line in f:
                line = line.strip()

                # get comments/meta-data
                preComment = line
                comment = ""
                if '[' in line:
                    split = line.split( '[', 1 )
                    preComment = split[ 0 ]
                    comment = "[" + split[ 1 ]
                elif '#' in line:
                    split = line.split( '#', 1 )
                    preComment = split[ 0 ]
                    comment = "#" + split[ 1 ]

                # get name and value (if any)
                name = ""
                value = ""
                split = preComment.split()
                if len( split ) >= 2:
                    name = split[ 0 ]
                    value = split[ 1 ]

                # append an entry (even for blank lines)
                self._entries.append( ConfigEntry( name, value, comment ) )

    # save this configuration file (in SBL format)
    def save( self, fileName ):
        f = open( fileName, "w" )
        if f:
            for e in self._entries:
                if e.name:
                    f.write( e.name )
                    f.write( " " )
                    f.write( strProc( e.value ) )
                    if e.comment:
                        f.write( " " )
                if e.comment:
                    f.write( e.comment )
                f.write( "\n" )
            

# provides a simple interface to SBL commands
class CommandRouter:

    # return true if user has requested that the current command stop running
    def checkCommandCancel( self ):
        return sblc.checkCommandEvents()

    # display a message
    def disp( self, indent, message ):
        sblc.disp( 0, indent, message )

    # display a warning
    def warning( self, message ):
        sblc.disp( 1, 0, message )

    # display a fatal error (will terminate program)
    def fatalError( self, message ):
        sblc.disp( 2, 0, message )

    # assume all other method calls are commands; send to SBL C++ command system
    def __getattr__( self, name ):
        if not name.startswith( "_" ):
            def runCommand( *args, **keywords ):
                cmdStr = name + " " + " ".join( [strProc( a ) for a in args] )
                sblc.execCommand( cmdStr )
            return runCommand
        else:
            raise AttributeError
