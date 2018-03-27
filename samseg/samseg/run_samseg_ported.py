import logging
import os
import numpy as np

from samseg.command_arguments import parse_args
from samseg.kvl_read_compression_lookup_table import kvlReadCompressionLookupTable
from samseg.kvl_read_shared_gmm_parameters import kvlReadSharedGMMParameters
from samseg.process_timer import ProcessTimer
from samseg.register_atlas_ported import samseg_registerAtlas
from samseg.run_utilities import find_or_create_save_path, specify_model, determine_optimization_options, \
    find_samseg_data_dir
from gems2python import GEMS2Python

from samseg.samseg_ported import samsegment

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)  # TODO: configurable logging

def force_fortran_order(funcname):
    old_func = getattr(np, funcname)
    def new_func(*args, **kwargs):
        kwargs['order'] = 'F'
        return old_func(*args, **kwargs)
    setattr(np, funcname, new_func)

for funcname in ('array', 'zeros', 'empty', 'zeros_like', 'empty_like'):
    force_fortran_order(funcname)

def run_samseg_from_cmdargs(cmdargs):
    # Run with no arguments to get help
    verbose = cmdargs.verbose
    atlas_only = cmdargs.atlas_only
    savePath = cmdargs.output
    numberOfThreads = cmdargs.threads
    showFigures = cmdargs.showfigs
    noBrainMasking = cmdargs.nobrainmask
    useDiagonalCovarianceMatrices = cmdargs.diagcovs
    RegMatFile = cmdargs.regmat
    InitLTAFile = cmdargs.InitLTAFile
    imageFileNames = cmdargs.image_file_names
    # Display input
    display_cmdargs(cmdargs)

    return run_samseg(
            imageFileNames,
            savePath,
            showFigures,
            noBrainMasking,
            useDiagonalCovarianceMatrices,
            verbose,
            numberOfThreads,
            atlas_only=atlas_only,
            InitLTAFile=InitLTAFile
    )


def run_samseg(
    imageFileNames,
    savePath,
    showFigures=False,
    noBrainMasking=False,
    useDiagonalCovarianceMatrices=False,
    verbose=False,
    numberOfThreads=None,
    InitLTAFile=None,
    atlas_only=False,
    checkpoint_manager=None
):
    # Create the output folder
    savePath = find_or_create_save_path(savePath)

    # set SAMSEG_DATA_DIR as an environment variable, eg,
    # setenv SAMSEG_DATA_DIR /autofs/cluster/koen/koen/GEMSapplications/wholeBrain
    samsegDataDir = find_samseg_data_dir()

    process_timer = ProcessTimer('samseg begin')
    # Specify the maximum number of threads the C++ stuff will use. The more threads you can use
    # the faster, provided you have a matching amount of cores on your system - up to
    # the point where memory becomes the bottle neck.
    # If the following command is not provided, the number of cores on your system will be used
    if numberOfThreads is not None:
        GEMS2Python.setGlobalDefaultNumberOfThreads(numberOfThreads)

    # Affine registration
    templateFileName = os.path.join(samsegDataDir, 'template.nii')
    affineRegistrationMeshCollectionFileName = os.path.join(samsegDataDir, 'atlasForAffineRegistration.txt.gz')
    worldToWorldTransformMatrix = None
    # TODO: load RegMatFile
    # if ( ~isempty( RegMatFile ) )
    #   load( RegMatFile, 'worldToWorldTransformMatrix' );
    # end
    #
    worldToWorldTransformMatrix, transformedTemplateFileName = samseg_registerAtlas(
        imageFileNames[0],
        affineRegistrationMeshCollectionFileName,
        templateFileName,
        savePath,
        showFigures,
        worldToWorldTransformMatrix,
        InitLTAFile,
        checkpoint_manager
    )
    process_timer.mark_time('registration done')
    if atlas_only:
        print('Registration-only requested, so quiting now')
        return
    # FreeSurfer (http://surfer.nmr.mgh.harvard.edu) has a standardized way of representation segmentations,
    # both manual and automated, as images in which certain intensity levels correspond to well-defined
    # anatomical structures - for instance an intensity value 17 always corresponds to the left hippocampus.
    # The text file FreeSurferColorLUT.txt distributed with FreeSurfer contains all these definitions, as well
    # as a color for each structure in
    # RGBA (Red-Green-Blue and Alpha (opaqueness)) format with which FreeSurfer will always represent segmented
    # structures in its visualization tools
    #
    # Let's read the contents of the "compressionLookupTable.txt" file, and show the names of the structures
    # being considered. The results are automatically sorted according to their "compressed label", i.e., the
    # first result corresponds to the first entry in the vector of probabilities associated with each node in
    # our atlas mesh.
    compressionLookupTableFileName = os.path.join(samsegDataDir, 'compressionLookupTable.txt')
    FreeSurferLabels, names, colors = kvlReadCompressionLookupTable(compressionLookupTableFileName)

    # Because we have many labels to segment, and each of these labels has its own Gaussian mixture model
    # whose parameters (mean, variance, mixture weight) we have to estimate from the data, it may makes sense to restrict
    # the degrees of freedom in the model somewhat by specifying that some of these labels have the same parameters
    # governing their Gaussian mixture model. For example, we'd expect no intensity differences between the left and right
    # part of each structure.
    # The way we implement this is by defining "super-structures" (i.e., a global white matter tissue class), and therefore
    # work with a simplied ("reduced") model during the entire parameter estimation phase. At the same time we also build
    # an inverse lookup table (mapping from original class number onto a reduced class number (super-structure)) that we
    # will need to compute the final segmentation.
    sharedGMMParametersFileName = os.path.join(samsegDataDir, 'sharedGMMParameters.txt')
    sharedGMMParameters = kvlReadSharedGMMParameters(sharedGMMParametersFileName)

    # Set various model specifications
    modelSpecifications = specify_model(FreeSurferLabels, noBrainMasking, useDiagonalCovarianceMatrices, sharedGMMParameters, names,
                                        colors, samsegDataDir)

    # Set various optimization options
    optimizationOptions = determine_optimization_options(verbose, samsegDataDir)

    FreeSurferLabels, names, volumesInCubicMm = samsegment(imageFileNames, transformedTemplateFileName,
                                                           modelSpecifications, optimizationOptions,
                                                           savePath, showFigures,
                                                           checkpoint_manager)
    print('names', names)
    print('volumesInCubicMm', volumesInCubicMm)
    process_timer.mark_time('samseg done')

def display_cmdargs(cmdargs):
    log_image_file_names(cmdargs.image_file_names)
    logger.info("output to %s", cmdargs.output)
    logger.info("threads=%d", cmdargs.threads)
    logger.info("init lta is %s", cmdargs.InitLTAFile)
    log_mode('verbose', cmdargs.verbose)


def log_image_file_names(image_file_names, title='image file names'):
    logger.info('%s:', title)
    for image_file_name in image_file_names:
        logger.info('    %s', image_file_name)


def log_mode(name, is_on):
    value = 'on' if is_on else 'off'
    logger.info('%s is %s', name, value)


if __name__ == '__main__':
    run_samseg_from_cmdargs(parse_args())
