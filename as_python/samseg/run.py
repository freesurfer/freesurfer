# function retval = run_samseg(varargin)
# % Run with no arguments to get help
import logging

from easydict import EasyDict

from as_python.samseg.command_arguments import parse_args
from as_python.samseg.process_timer import ProcessTimer
from as_python.samseg.register_atlas import samseg_register_atlas
from as_python.samseg.run_utilities import update_recipe_with_calculated_paths, determine_transformed_template_filename, \
    determine_optimization_options, specify_model, \
    determine_shared_gmm_parameters, use_standard_affine_registration_atlas
from as_python.samseg.samsegment import samsegment

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)  # TODO: configurable logging


def run_samseg(recipe):
    process_timer = ProcessTimer('samseg begin')
    recipe = update_recipe_with_calculated_paths(recipe)
    display_recipe(recipe)
    shared_gmm_parameters = run_or_retrieve_registration_process(recipe)
    process_timer.mark_time('registration done')
    run_segmentation_process(recipe, shared_gmm_parameters)
    process_timer.mark_time('samseg done')


def display_recipe(recipe):
    log_image_file_names(recipe.image_file_names)
    logger.info("output to %s", recipe.save_path)
    logger.info("threads=%d", recipe.threads)
    log_mode('exvivo', recipe.exvivo)
    log_mode('verbose', recipe.verbose)
    log_missing_structures(recipe.missing_structures)
    logger.info('atlas_directory is %s', recipe.avg_data_dir)


def log_image_file_names(image_file_names, title='image file names'):
    logger.info('%s:', title)
    for image_file_name in image_file_names:
        logger.info('    %s', image_file_name)


def log_missing_structures(missing_structures):
    if missing_structures:
        logger.info("missing structures:")
        for missing_structure in missing_structures:
            logger.info("    %s", missing_structure)
    else:
        logger.info("no missing structures")


def log_mode(name, is_on):
    value = 'on' if is_on else 'off'
    logger.info('%s is %s', name, value)


def run_or_retrieve_registration_process(recipe):
    shared_gmm_parameters = determine_shared_gmm_parameters(recipe.exvivo)
    if recipe.exvivo:
        recipe.affine_file_names = create_tailored_affine_registration_atlas(recipe)
    else:
        recipe.affine_file_names = use_standard_affine_registration_atlas(recipe.avg_data_dir)

    if recipe.regmat:
        world_to_world_transform_matrix = retrieve_registration_process(recipe)
    else:
        world_to_world_transform_matrix = run_registration_process(recipe)

    create_and_write_transformations(recipe, world_to_world_transform_matrix)
    return shared_gmm_parameters


def create_tailored_affine_registration_atlas(recipe):
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


def retrieve_registration_process():
    #   fprintf('Not performing registration:\n');
    #   fprintf('  Loading reg file %s\n',RegMatFile);
    #   load(RegMatFile);
    #   fname = sprintf('%s/SPM12_6classes_30x30x30_template_coregistrationMatrices.mat',savePath);
    #   save(fname,'worldToWorldTransformMatrix','imageToImageTransformMatrix');
    pass


def run_registration_process(recipe):
    registration_recipe = EasyDict()
    registration_recipe.verbose = recipe.verbose
    registration_recipe.image_file_name = recipe.image_file_names[0]
    registration_recipe.mesh_collection_file_name = \
        recipe.affine_file_names.mesh_collection_file_name
    registration_recipe.template_file_name = \
        recipe.affine_file_names.template_file_name
    registration_recipe.save_path = recipe.save_path
    registration_recipe.show_figures = recipe.show_registration_figures

    world_to_world_transform_matrix = samseg_register_atlas(registration_recipe)

    #         affine_registration_mesh_collection_file_name,
    #         affine_registration_template_file_name,
    #         image_file_names,
    #         save_path,
    #         show_figures
    # ):
    #   fprintf('entering registerAtlas\n');
    #   showFigures = false;
    #   worldToWorldTransformMatrix = samseg_registerAtlas( imageFileNames{ 1 }, ...
    #                                                       affineRegistrationMeshCollectionFileName, ...
    #                                                       affineRegistrationTemplateFileName, ...
    #                                                       savePath, ...
    #                                                       showFigures );
    return world_to_world_transform_matrix


def create_and_write_transformations(recipe, world_to_world_transform_matrix):
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
    optimization_options = determine_optimization_options(recipe.verbose)

    # [ FreeSurferLabels, names, volumesInCubicMm ] = samsegment( imageFileNames, transformedTemplateFileName, meshCollectionFileName, ...
    #                                                             compressionLookupTableFileName, modelSpecifications, ...
    #                                                             optimizationOptions, savePath, showFigures );
    [free_surfer_labels, names, volumes_in_cubic_mm] = samsegment(
        recipe,
        transformed_template_filename,
        model_specifications,
        optimization_options,
    )

    show_segmentation_results(names, volumes_in_cubic_mm)


def show_segmentation_results(names, volumes_in_cubic_mm):
    logger.info('volumes in cubic meters:')
    for name, volume in zip(names, volumes_in_cubic_mm):
        logger.info('   %s=%f', name, volume)


if __name__ == '__main__':
    run_samseg(parse_args())
