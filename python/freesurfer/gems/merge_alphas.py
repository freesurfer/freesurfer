import logging
import math
import numpy as np
import colorsys


logger = logging.getLogger(__name__)


def meshValidityTest(alphas, name):
    probability_discrepancy = np.max(np.abs(np.sum(alphas, axis=1) - 1))
    if probability_discrepancy > 1e-5:
        message = '%s invalid: class probabilities in the mesh should sum to one in all nodes' % name
        raise ValueError(message)


def kvlMergeAlphas(alphas, names, mergeOptions, FreeSurferLabels=None, colors=None):
    #
    # Creates a 'mergedAlphas' matrix where one or more columns of the 'alphas'
    # matrix have been "merged" (i.e., added together)
    #

    alpha_count, label_count = alphas.shape
    if FreeSurferLabels is None:
        FreeSurferLabels = [math.nan] * label_count
    if colors is None:
        colors = [[math.nan] * 4] * label_count

    # Make sure we're dealing with a valid mesh
    meshValidityTest(alphas, 'alphas')


    # ath ---------------------------------------------------------------

    # to-port: numberOfStructures = size(alphas, 2);  # Number of freesurfer labels being segmented/outputed
    # to-port: translationTable = zeros( numberOfClasses, numberOfStructures );
    mergedNames = [];
    for classNumber, mergeOption in enumerate(mergeOptions):
        mergedNames.append(mergeOption.mergedName.strip())
        for searchString in mergeOption.searchStrings:
            # to-port: some search stuff
            translationTable[classNumber, structureNumbers] = 1
    # to-port: error check
    # to-port: translationTable = translationTable ./ sum( translationTable );
    mergedAlphas = alphas * translationTable.T

    # ath ---------------------------------------------------------------



    # Start by copying everything
    mergedAlphas = np.copy(alphas)
    mergedNames = [str(name).strip() for name in names]
    mergedFreeSurferLabels = [label for label in FreeSurferLabels]
    mergedColors = [color for color in colors]

    # Also keep track of what class ends up being merged into what
    mergingTable = [set([index]) for index in range(mergedAlphas.shape[1])]

    zeroes_for_appending = np.zeros((alpha_count, 1), dtype=np.double)
    for mergeOption in mergeOptions:
        logger.debug('starting shape = %s', str(mergedAlphas.shape))
        mergedName = mergeOption.mergedName.strip()
        searchStrings = mergeOption.searchStrings
        if isinstance(searchStrings, str):
            searchStrings = [searchStrings]
        searchStrings = {str(value).strip() for value in searchStrings}
        def is_target_string(candidate):
            for pattern in searchStrings:
                if pattern in candidate:
                    return True
            return False
            logger.debug('searching for %s', searchStrings)
        # Retrieve the class number corresponding to the merged label. If it doesn't exist yet,
        # create a new class
        matching_name_indices = [index for index, name in enumerate(mergedNames) if name == mergedName]
        if len(matching_name_indices) == 1:
            classNumberOfMerged = matching_name_indices[0]
            logger.debug('keeping class number %d', classNumberOfMerged)
        else:
            classNumberOfMerged = mergedAlphas.shape[1]
            logger.debug('adding class number %d', classNumberOfMerged)
            mergedAlphas = np.append(mergedAlphas, zeroes_for_appending, axis=1)
            mergedNames.append(mergedName)
            mergedFreeSurferLabels.append(None)
            mergedColors.append(None)
            mergingTable.append(set())
        # Get class numbers of all structures to be merged into classNumberOfMerged
        mergingClassNumbers = [index for index, name in enumerate(mergedNames) if is_target_string(name)]
        logger.debug('merging from %s', str(mergingClassNumbers))
        # Now merge (i.e., add the probabilities) of those class numbers to that of classNumberOfMerged
        alphaValuesToMerge = np.zeros((alpha_count, ), dtype=np.double)
        for mergeIndex in mergingClassNumbers:
            alphaValuesToMerge += mergedAlphas[:, mergeIndex]
            mergingTable[classNumberOfMerged] |= mergingTable[mergeIndex]
        mergedAlphas[:, classNumberOfMerged] = alphaValuesToMerge
        # Remove the corresponding rows/columns in the resulting variables
        disappearingClassNumbers = [number for number in mergingClassNumbers if number != classNumberOfMerged]
        disappearingClassNumbers.sort()
        disappearingClassNumbers.reverse()
        logger.debug('removing these %s', str(disappearingClassNumbers))
        for disappearingClassNumber in disappearingClassNumbers:
            mergedAlphas = np.delete(mergedAlphas, disappearingClassNumber, 1)
        mergedFreeSurferLabels = [label for index, label in enumerate(mergedFreeSurferLabels) if not index in disappearingClassNumbers]
        mergedNames = [name for index, name in enumerate(mergedNames) if not index in disappearingClassNumbers]
        mergedColors = [color for index, color in enumerate(mergedColors) if not index in disappearingClassNumbers]
        mergingTable = [sources for index, sources in enumerate(mergingTable) if not index in disappearingClassNumbers]
        logger.debug('ending shape = %s', str(mergedAlphas.shape))
    # Make sure we still have a valid mesh
    meshValidityTest(mergedAlphas, 'mergedAlphas')
    # Take care of NaN's we temporily put in the mergedFreeSurferLabels and mergedColors
    newlyInventedClasses = [index for index, label in enumerate(mergedFreeSurferLabels) if label is None]
    numberOfNewlyInvitedClasses = len(newlyInventedClasses)
    saturation = 1.0
    color_value = 1.0
    for index, newIndex in enumerate(newlyInventedClasses):
        mergedFreeSurferLabels[newIndex] = -(1 + index)
        hue = index / numberOfNewlyInvitedClasses
        mergedColors[newIndex]= [255 * component for component in colorsys.hsv_to_rgb(hue, saturation, color_value)] + [255]
    # Also compute lookup table indicating for each original class number (column number in alphas) the final
    # class number (column number) in mergedAlphas
    translationTable = [0] * label_count
    for mergedClassNumber, mergingTableValueSet in enumerate(mergingTable):
        for mergingTableValue in mergingTableValueSet:
            translationTable[mergingTableValue] = 1 + mergedClassNumber # 1 based index
    return mergedAlphas, mergedNames, mergedFreeSurferLabels, mergedColors, translationTable
