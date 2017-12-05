# function retval = run_samseg(varargin)
# % Run with no arguments to get help
import logging

from as_python.samseg.command_arguments import parse_args
from as_python.samseg.process_timer import ProcessTimer
from as_python.samseg.run_utilities import update_recipe_with_calculated_paths, determine_transformed_template_filename, \
    exvivo_shared_gmm_parameters, standard_shared_gmm_parameters
from as_python.samseg.samsegment import samsegment

logger = logging.getLogger(__name__)

BAD_RESULT = 1  # retval = 1;
GOOD_RESULT = 0


def run_samseg(recipe):
    process_timer = ProcessTimer()
    update_recipe_with_calculated_paths(recipe)
    shared_gmm_parameters = run_or_retrieve_registration_process(recipe)
    process_timer.mark_time('registration done')
    run_segmentation_process(recipe, shared_gmm_parameters)
    process_timer.mark_time('samseg done')


def run_or_retrieve_registration_process(recipe):
    display_registration_input(recipe)
    if recipe.exvivo:
        shared_gmm_parameters = exvivo_shared_gmm_parameters()
        recipe.affine_file_names = build_tailored_affine_registration_atlas(recipe)
    else:
        shared_gmm_parameters = standard_shared_gmm_parameters()
        recipe.affine_file_names = build_standard_affine_registration_atlas(recipe)

    if recipe.regmat:
        retrieve_registration_process(recipe)
    else:
        run_registration_process(recipe)

    create_and_write_transformations(recipe)
    return shared_gmm_parameters


def display_registration_input(recipe):
    logger.info("inputs:")
    for image_file_name in recipe.image_file_names:
        logger.info('    %s', image_file_name)
    logger.info("output to %s", recipe.save_path)
    logger.info("threads=%d", recipe.threads)
    if recipe.exvivo:
        logger.info("exvivo mode is on")
    else:
        logger.info("exvivo mode is off")
    if recipe.missing_structures:
        logger.info("missing structures:")
        for missing_structure in recipe.missing_structures:
            logger.info("    %s", missing_structure)
    else:
        logger.info("no missing structures")


def build_tailored_affine_registration_atlas(recipe):
    #   % Create a tailor-made atlas for affine registration purposes
    #
    #   % Read mesh
    #   meshCollection = kvlReadMeshCollection( meshCollectionFileName );
    #   mesh = kvlGetMesh( meshCollection, -1 );
    #   [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
    #
    #   % Get a Matlab matrix containing a copy of the probability vectors in each mesh node (size numberOfNodes x
    #   % numberOfLabels ).
    #   alphas = kvlGetAlphasInMeshNodes( mesh );
    #
    #   % Remove non-existing structures
    #   mergeOptions = struct;
    #   mergeOptions( 1 ).mergedName = 'Unknown';
    #   mergeOptions( 1 ).searchStrings = missingStructureSearchStrings;
    #   [ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, mergeOptions, FreeSurferLabels, colors );
    #
    #   % Get global tissue types
    #   [ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, sharedGMMParameters, FreeSurferLabels, colors );
    #
    #   % Additionally move some deep gray matter structures into global GM
    #   mergeOptions = struct;
    #   mergeOptions( 1 ).mergedName = 'Global GM';
    #   mergeOptions( 1 ).searchStrings = { 'Thalamus', 'Pallidum', 'Putamen' };
    #   [ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, mergeOptions, FreeSurferLabels, colors );
    #
    #   % Create tailored atlas
    #   kvlSetAlphasInMeshNodes( mesh, alphas );
    #   [ template, transform ] = kvlReadImage( templateFileName );
    #   templateImageBuffer = kvlGetImageBuffer( template );
    #   priors = kvlRasterizeAtlasMesh( mesh, size( templateImageBuffer ) );
    #   transformMatrix = kvlGetTransformMatrix( transform );
    #   [ affineRegistrationMeshCollectionFileName, affineRegistrationTemplateFileName ] = ...
    #                                                 createAtlasMeshForAffineRegistration( priors, transformMatrix, savePath );
    pass


def build_standard_affine_registration_atlas(recipe):
    #   affineRegistrationMeshCollectionFileName = sprintf( '%s/SPM12_6classes_30x30x30_meshCollection.txt.gz', AvgDataDir );
    #   affineRegistrationTemplateFileName = sprintf( '%s/SPM12_6classes_30x30x30_template.nii', AvgDataDir );
    pass


def retrieve_registration_process():
    #   fprintf('Not performing registration:\n');
    #   fprintf('  Loading reg file %s\n',RegMatFile);
    #   load(RegMatFile);
    #   fname = sprintf('%s/SPM12_6classes_30x30x30_template_coregistrationMatrices.mat',savePath);
    #   save(fname,'worldToWorldTransformMatrix','imageToImageTransformMatrix');
    pass


def run_registration_process(
        affine_registration_mesh_collection_file_name,
        affine_registration_template_file_name,
        image_file_names,
        save_path,
        show_figures
):
    #   fprintf('entering registerAtlas\n');
    #   showFigures = false;
    #   worldToWorldTransformMatrix = samseg_registerAtlas( imageFileNames{ 1 }, ...
    #                                                       affineRegistrationMeshCollectionFileName, ...
    #                                                       affineRegistrationTemplateFileName, ...
    #                                                       savePath, ...
    #                                                       showFigures );
    pass


def create_and_write_transformations(kvl, template_file_name):
    # % For historical reasons the samsegment script figures out the affine transformation from
    # % a transformed MNI template (where before transformation this template defines the segmentation
    # % mesh atlas domain). This is a bit silly really, but for now let's just play along and make
    # % sure we generate it
    # [ origTemplate, origTemplateTransform ] = kvlReadImage( templateFileName );
    # kvlWriteImage( origTemplate, transformedTemplateFileName, ...
    #                kvlCreateTransform( double( worldToWorldTransformMatrix * kvlGetTransformMatrix( origTemplateTransform ) ) ) );
    pass


def run_segmentation_process(recipe, shared_gmm_parameters):
    #         compression_lookup_table_file_name,
    #         exvivo,
    #         kvl,
    #         image_file_names,
    #         mesh_collection_file_name,
    #         missing_structure_search_strings,
    #         save_path,
    #         shared_gmm_parameters,
    #         show_figures
    # ):
    transformed_template_filename = determine_transformed_template_filename(recipe.save_path)

    model_specifications = specify_model(recipe.exvivo, recipe.missing_structures, shared_gmm_parameters)
    optimization_options = determine_optimization_options()

    # [ FreeSurferLabels, names, volumesInCubicMm ] = samsegment( imageFileNames, transformedTemplateFileName, meshCollectionFileName, ...
    #                                                             compressionLookupTableFileName, modelSpecifications, ...
    #                                                             optimizationOptions, savePath, showFigures );
    [FreeSurferLabels, names, volumesInCubicMm] = samsegment(
        recipe,
        transformed_template_filename,
        model_specifications,
        optimization_options,
    )

    show_segmentation_results(names, volumesInCubicMm)


def specify_model(exvivo, missing_structure_search_strings, shared_gmm_parameters):
    # exvivo = cmdargs.exvivo;
    # % Set various model specifications
    # modelSpecifications = struct;
    # modelSpecifications.missingStructureSearchStrings = missingStructureSearchStrings;
    # modelSpecifications.sharedGMMParameters = sharedGMMParameters;
    # modelSpecifications.useDiagonalCovarianceMatrices = false;
    # modelSpecifications.brainMaskingSmoothingSigma = 3; % sqrt of the variance of a Gaussian blurring kernel
    # modelSpecifications.brainMaskingThreshold = 0.01;
    # modelSpecifications.K = 0.1; % Stiffness of the mesh
    # modelSpecifications.biasFieldSmoothingKernelSize = 50.0;  % Distance in mm of sinc function center to first zero crossing
    # if exvivo
    #   modelSpecifications.brainMaskingThreshold = -Inf; % Disable brain masking
    #   modelSpecifications.useDiagonalCovarianceMatrices = true;
    # end
    pass


def determine_optimization_options():
    # % Set various optimization options
    # optimizationOptions = struct;
    # optimizationOptions.multiResolutionSpecification = struct;
    # optimizationOptions.multiResolutionSpecification( 1 ).meshSmoothingSigma = 2.0; % In mm
    # optimizationOptions.multiResolutionSpecification( 1 ).targetDownsampledVoxelSpacing = 2.0; % In mm
    # optimizationOptions.multiResolutionSpecification( 1 ).maximumNumberOfIterations = 100;
    # optimizationOptions.multiResolutionSpecification( 1 ).estimateBiasField = true;
    # optimizationOptions.multiResolutionSpecification( 2 ).meshSmoothingSigma = 0.0; % In mm
    # optimizationOptions.multiResolutionSpecification( 2 ).targetDownsampledVoxelSpacing = 1.0; % In mm
    # optimizationOptions.multiResolutionSpecification( 2 ).maximumNumberOfIterations = 100;
    # optimizationOptions.multiResolutionSpecification( 2 ).estimateBiasField = true; % Switching this off will use the bias field estimated
    #                                                                                 % at lower resolution(s)
    # optimizationOptions.maximumNumberOfDeformationIterations = 20;
    # optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion = 1e-4;
    # optimizationOptions.verbose = 0;
    # optimizationOptions.maximalDeformationStopCriterion = 0.001; % Measured in pixels
    # optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion = optimizationOptions.maximalDeformationStopCriterion; % Idem
    # % optimizationOptions.relativeCostDecreaseStopCriterion = 1e-6;
    # optimizationOptions.maximalDeformationAppliedStopCriterion = 0.0;
    # optimizationOptions.BFGSMaximumMemoryLength = 12;
    pass


def show_segmentation_results(names, volumesInCubicMm):
    # names
    # volumesInCubicMm
    pass


if __name__ == '__main__':
    run_samseg(parse_args())
