import logging
import math
import numpy as np
import colorsys


def meshValidityTest(alphas, name):
    probability_discrepancy = np.max(np.abs(np.sum(alphas, axis=1) - 1))
    if probability_discrepancy > 1e-5:
        message = '%s invalid: class probabilities in the mesh should sum to one in all nodes' % name
        raise ValueError(message)


def kvlGetMergingFractionsTable( names, mergeOptions ):
    '''Computes a numerOfClasses x numberOfStructures matrix where each column indicates the fractions 
    of the various classes (super-structures) in the corresponding structure. So each column sums to 1.'''

    #
    numberOfClasses = len( mergeOptions )

    fractionsTable = np.zeros( ( numberOfClasses, len( names ) ) )
    mergedNames = [];
    for classNumber, mergeOption in enumerate( mergeOptions ):
        mergedNames.append( mergeOption.mergedName.strip() )
        for searchString in mergeOption.searchStrings:
            for structureNumber, name in enumerate( names ):
                if searchString in name:
                    fractionsTable[ classNumber, structureNumber ] = 1.0

    if not fractionsTable.any( axis=0 ).all():
        raise ValueError( 'some structures are not associated with any super-structures' )

    fractionsTable = fractionsTable / np.sum(fractionsTable, 0)

    # Print out merge info
    for classNumber in range( numberOfClasses ):
        print( mergedNames[ classNumber ] )
        for structureNumber in range( len( names ) ):
            percentage = int( fractionsTable[ classNumber, structureNumber ] * 100 )
            if percentage > 0:
                print( '    %s (%d%%)' % ( names[ structureNumber ].ljust( len( max( names, key=len ) ) ), percentage ) )

    return fractionsTable, mergedNames


def kvlMergeAlphas( alphas, fractionsTable ):
    '''Creates a 'mergedAlphas' matrix where one or more columns of the 'alphas'
    matrix have been "merged" (i.e., added together).'''

    # Make sure we're dealing with a valid mesh
    meshValidityTest( alphas, 'alphas' )

    # Do the actual merging
    mergedAlphas = np.dot( alphas, fractionsTable.T )

    # Make sure we're dealing with a valid mesh
    meshValidityTest( mergedAlphas, 'mergedAlphas' )

    return mergedAlphas
