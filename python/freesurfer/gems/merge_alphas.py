import logging
import math
import numpy as np
import colorsys


def meshValidityTest(alphas, name):
    probability_discrepancy = np.max(np.abs(np.sum(alphas, axis=1) - 1))
    if probability_discrepancy > 1e-5:
        message = '%s invalid: class probabilities in the mesh should sum to one in all nodes' % name
        raise ValueError(message)


def kvlMergeAlphas(alphas, names, mergeOptions, FreeSurferLabels=None, colors=None):
    '''Creates a 'mergedAlphas' matrix where one or more columns of the 'alphas'
    matrix have been "merged" (i.e., added together).'''

    alpha_count, label_count = alphas.shape
    if FreeSurferLabels is None:
        FreeSurferLabels = [math.nan] * label_count
    if colors is None:
        colors = [[math.nan] * 4] * label_count

    # Make sure we're dealing with a valid mesh
    meshValidityTest(alphas, 'alphas')

    numberOfClasses = len(mergeOptions)

    translationTable = np.zeros((numberOfClasses, len(names)))
    mergedNames = [];
    for classNumber, mergeOption in enumerate(mergeOptions):
        mergedNames.append(mergeOption.mergedName.strip())
        for searchString in mergeOption.searchStrings:
            for structureNumber, name in enumerate(names):
                if searchString in name:
                    translationTable[classNumber, structureNumber] = 1.0

    if not translationTable.any():
        raise ValueError('some structures are not associated with any super-structures')

    translationTable = translationTable / np.sum(translationTable, 0)

    mergedAlphas = np.dot(alphas, translationTable.T)
    meshValidityTest(mergedAlphas, 'mergedAlphas')

    mergedFreeSurferLabels = -(np.arange(numberOfClasses) + 1)

    mergedColors = []
    for n in range(numberOfClasses):
        hue = n / numberOfClasses
        mergedColors.append([255 * component for component in colorsys.hsv_to_rgb(hue, 1, 1)] + [255])

    # Print out merge info
    for classNumber in range(numberOfClasses):
        print(mergedNames[classNumber])
        for structureNumber in range(len(names)):
            percentage = int(translationTable[classNumber, structureNumber] * 100)
            if percentage > 0:
                print('    %s (%d%%)' % (names[structureNumber].ljust(len(max(names, key=len))), percentage))

    return mergedAlphas, mergedNames, mergedFreeSurferLabels, mergedColors, translationTable
