from collections import namedtuple
GMMparameter = namedtuple('GMMparameter', 'mergedName numberOfComponents searchStrings')
def kvlReadSharedGMMParameters( fileName ):
    #
    # Read text file where each line is formatted as:
    #
    #  mergedName numberOfComponents searchString(s)
    #
    print( 'Reading contexts of file {0}'.format(fileName))

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
