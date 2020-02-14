from _operator import itemgetter
from collections import namedtuple


GMMparameter = namedtuple('GMMparameter', 'mergedName numberOfComponents searchStrings')


def kvlReadSharedGMMParameters(fileName):
    # Read text file where each line is formatted as:
    #   mergedName numberOfComponents searchString(s)
    sharedGMMParameters = []
    with open(fileName) as fid:
        for textLine in fid.readlines():
            # Remove leading and trailing blanks
            components = textLine.strip().split()
            if len(components) > 0 and components[0] != '#':
                if len(components) < 2:
                    raise ValueError( 'Can''t parse line: {0}'.format(textLine))
                mergedName = components[0]
                numberOfComponents = int(components[1])
                searchStrings = components[2:]
                # Add to sharedGMMParameters structure array
                sharedGMMParameters.append(GMMparameter(mergedName, numberOfComponents, searchStrings))
    return sharedGMMParameters


def kvlWriteSharedGMMParameters( sharedGMMParameters, fileName ):
    #
    with open( fileName, 'w' ) as fid:
        print( '# The format is: mergedName numberOfComponents searchStrings\n\n', file=fid ) 
        numberOfStructures = len( sharedGMMParameters )
        
        for GMMparameter in sharedGMMParameters:
            print( GMMparameter.mergedName, GMMparameter.numberOfComponents, 
                  *GMMparameter.searchStrings, file=fid )


def kvlReadCompressionLookupTable(fileName):
    # Format is "FreeSurferLabel compressedLabel name R G B A"
    table = []
    with open(fileName) as fid:
        for line in fid.readlines():
            FreeSurferLabel, compressedLabel, name, R, G, B, A = [
                data_type(value) for data_type, value in zip(
                    (int, int, str, int, int, int, int),
                    line.split())]
            # Add contents to output matrices
            table.append({
                'FreeSurferLabel': FreeSurferLabel,
                'compressedLabel': compressedLabel,
                'name': name,
                'color': [R, G, B, A],
            })
    # Sort output according to compressedLabel
    table = sorted(table, key=itemgetter('compressedLabel'))
    FreeSurferLabels, names, colors = [[entry[key] for entry in table] for key in ['FreeSurferLabel', 'name', 'color']]
    return FreeSurferLabels, names, colors


def kvlWriteCompressionLookupTable( fileName, FreeSurferLabels, names, colors ):
    #
    # Format is "FreeSurferLabel compressedLabel name R G B A"
    #
    with open( fileName, 'w' ) as fid:
        #
        for i in range( 0, len( names ) ):
            #
            print( FreeSurferLabels[ i ], i, names[ i ], *colors[ i ], file=fid )


