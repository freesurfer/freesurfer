# function retval = run_samseg(varargin)
# % Run with no arguments to get help
from samseg.command_arguments import create_cmdargs
from samseg.kvl import KVL
from samseg.process_timer import ProcessTimer
from samseg.samsegment import samsegment

BAD_RESULT = 1  # retval = 1;
GOOD_RESULT = 0


def run_samseg(varargin):
    process_timer = ProcessTimer()

    cmdargs = create_cmdargs(varargin)
    if (cmdargs):
        kvl = create_kvl(cmdargs)
        exvivo = is_exvivo(cmdargs)
        save_path = find_or_create_save_path(cmdargs)
        image_file_names = determine_image_file_names(cmdargs)
        missing_structure_search_strings = determine_missing_structure_search_strings(cmdargs)
        avg_data_dir = find_avg_data_dir()
        mesh_collection_file_name = determine_mesh_collection_file_name(avg_data_dir)
        compression_lookup_table_file_name = determine_compression_lookup_table_file_name(avg_data_dir)

        shared_gmm_parameters = run_or_retrieve_registration_process(
            cmdargs,
            compression_lookup_table_file_name,
            exvivo,
            kvl,
            image_file_names,
            mesh_collection_file_name,
            missing_structure_search_strings,
            save_path,
        )

        run_segmentation_process(
            compression_lookup_table_file_name,
            kvl,
            image_file_names,
            exvivo,
            mesh_collection_file_name,
            missing_structure_search_strings,
            save_path,
            shared_gmm_parameters,
            show_figures=should_show_segmentation_figures(cmdargs)
        )
        process_timer.mark_time('samseg done')
        return GOOD_RESULT
    else:
        return BAD_RESULT


def create_kvl(cmdargs):
    number_of_threads = determine_number_of_threads(cmdargs)
    return KVL(number_of_threads)


def find_or_create_save_path():
    # savePath = cmdargs.outdir;
    # % Create the output folder
    # if ~exist(savePath, 'dir')
    #   mkdir(savePath);
    # end
    pass


def run_or_retrieve_registration_process(
        cmdargs,
        compression_lookup_table_file_name,
        exvivo,
        kvl,
        image_file_names,
        mesh_collection_file_name,
        missing_structure_search_strings,
        save_path,
):
    template_file_name = determine_template_file_name(cmdargs)
    number_of_threads = kvl.number_of_threads

    display_registration_input(
        image_file_names,
        save_path,
        number_of_threads,
        exvivo,
        missing_structure_search_strings
    )
    if exvivo:
        shared_gmm_parameters = exvivo_shared_gmm_parameters()
        affine_file_names = build_tailored_affine_registration_atlas(
            compression_lookup_table_file_name,
            kvl,
            mesh_collection_file_name,
            template_file_name
        )
    else:
        shared_gmm_parameters = standard_shared_gmm_parameters()
        affine_file_names = build_standard_affine_registration_atlas(kvl)
    affine_registration_mesh_collection_file_name, affine_registration_template_file_name = affine_file_names

    reg_mat_file = find_reg_mat_file(cmdargs)
    # if(isempty(RegMatFile))
    if (reg_mat_file):
        run_registration_process(
            affine_registration_mesh_collection_file_name,
            affine_registration_template_file_name,
            image_file_names,
            save_path,
            show_figures=should_show_registration_figures(cmdargs)

        )
    else:
        reg_mat_file = retrieve_registration_process()

    create_and_write_transformations(kvl, template_file_name)
    return shared_gmm_parameters


def determine_template_file_name(cmdargs):
    # % set SAMSEG_DATA_DIR as an environment variable, eg,
    # % setenv SAMSEG_DATA_DIR /autofs/cluster/koen/koen/GEMSapplications/wholeBrain
    # AvgDataDir = getenv( 'SAMSEG_DATA_DIR' );
    # templateFileName = sprintf('%s/mni305_masked_autoCropped.mgz',AvgDataDir);
    pass


def determine_number_of_threads(cmdargs):
    # numberOfThreads = cmdargs.nthreads;
    pass


def is_exvivo(cmdargs):
    # exvivo = cmdargs.exvivo;
    pass


def determine_missing_structure_search_strings(cmdargs):
    # numMissing = size(cmdargs.missingStructures,1);
    # missingStructureSearchStrings = cell(0,0);
    # for n = 1:numMissing
    #   missingStructureSearchStrings{n} = deblank(cmdargs.missingStructures(n,:));
    # end
    pass


def display_registration_input(
        image_file_names,
        save_path,
        number_of_threads,
        exvivo,
        missing_structure_search_strings
):
    # % Display input
    # imageFileNames
    # savePath
    # numberOfThreads
    # exvivo
    # missingStructureSearchStrings
    pass


def find_reg_mat_file(cmdargs):
    # RegMatFile = cmdargs.regmatfile;
    pass


def exvivo_shared_gmm_parameters():
    #   % Specify which classes share the same intensity Gaussian mixture model (and the number of components within each model)
    #   sharedGMMParameters = struct;
    #   sharedGMMParameters( 1 ).mergedName = 'Unknown'; % Background and what is normally CSF
    #   sharedGMMParameters( 1 ).searchStrings = { 'Unknown', 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus' };
    #   sharedGMMParameters( 1 ).numberOfComponents = 1;
    #   sharedGMMParameters( 2 ).mergedName = 'Global WM'; % WM
    #   sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
    #   sharedGMMParameters( 2 ).numberOfComponents = 1;
    #   sharedGMMParameters( 3 ).mergedName = 'Global GM'; % GM
    #   sharedGMMParameters( 3 ).searchStrings = { 'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities', 'Putamen' };
    #   sharedGMMParameters( 3 ).numberOfComponents = 1;
    #   sharedGMMParameters( 4 ).mergedName = 'Thalamus'; % Thalamus
    #   sharedGMMParameters( 4 ).searchStrings = { 'Thalamus' };
    #   sharedGMMParameters( 4 ).numberOfComponents = 1;
    #   sharedGMMParameters( 5 ).mergedName = 'Pallidum'; % Pallidum
    #   sharedGMMParameters( 5 ).searchStrings = { 'Pallidum' };
    #   sharedGMMParameters( 5 ).numberOfComponents = 1;
    pass


def standard_shared_gmm_parameters():
    #   % Specify which classes share the same intensity Gaussian mixture model (and the number of components within each model)
    #   sharedGMMParameters = struct;
    #   sharedGMMParameters( 1 ).mergedName = 'Unknown'; % Background
    #   sharedGMMParameters( 1 ).searchStrings = { 'Unknown'};
    #   sharedGMMParameters( 1 ).numberOfComponents = 3;
    #   sharedGMMParameters( 2 ).mergedName = 'Global WM'; % WM
    #   sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
    #   sharedGMMParameters( 2 ).numberOfComponents = 2;
    #   sharedGMMParameters( 3 ).mergedName = 'Global GM'; % GM
    #   sharedGMMParameters( 3 ).searchStrings = { 'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities' };
    #   sharedGMMParameters( 3 ).numberOfComponents = 3;
    #   sharedGMMParameters( 4 ).mergedName = 'Global CSF'; % CSF
    #   sharedGMMParameters( 4 ).searchStrings = { 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus' };
    #   sharedGMMParameters( 4 ).numberOfComponents = 3;
    #   sharedGMMParameters( 5 ).mergedName = 'Thalamus'; % Thalamus
    #   sharedGMMParameters( 5 ).searchStrings = { 'Thalamus' };
    #   sharedGMMParameters( 5 ).numberOfComponents = 2;
    #   sharedGMMParameters( 6 ).mergedName = 'Pallidum'; % Pallidum
    #   sharedGMMParameters( 6 ).searchStrings = { 'Pallidum' };
    #   sharedGMMParameters( 6 ).numberOfComponents = 2;
    #   sharedGMMParameters( 7 ).mergedName = 'Putamen'; % Putamen
    #   sharedGMMParameters( 7 ).searchStrings = { 'Putamen' };
    #   sharedGMMParameters( 7 ).numberOfComponents = 2;
    pass


def build_tailored_affine_registration_atlas(
        compression_lookup_table_file_name,
        kvl,
        mesh_collection_file_name,
        template_file_name
):
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


def build_standard_affine_registration_atlas(kvl):
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


def run_segmentation_process(
        compression_lookup_table_file_name,
        exvivo,
        kvl,
        image_file_names,
        mesh_collection_file_name,
        missing_structure_search_strings,
        save_path,
        shared_gmm_parameters,
        show_figures
):
    transformed_template_filename = determine_transformed_template_filename(save_path)

    model_specifications = specify_model(exvivo, missing_structure_search_strings, shared_gmm_parameters)
    optimization_options = determine_optimization_options()

    # [ FreeSurferLabels, names, volumesInCubicMm ] = samsegment( imageFileNames, transformedTemplateFileName, meshCollectionFileName, ...
    #                                                             compressionLookupTableFileName, modelSpecifications, ...
    #                                                             optimizationOptions, savePath, showFigures );
    [FreeSurferLabels, names, volumesInCubicMm] = samsegment(
        kvl,
        image_file_names,
        transformed_template_filename,
        mesh_collection_file_name,
        compression_lookup_table_file_name,
        model_specifications,
        optimization_options,
        save_path,
        show_figures
    )

    show_segmentation_results(names, volumesInCubicMm)


def determine_image_file_names(cmdargs):
    # ninvolfiles = size(cmdargs.involfiles,1);
    # imageFileNames = cell(0,0);
    # for n = 1:ninvolfiles
    #   imageFileNames{n} = deblank(cmdargs.involfiles(n,:));
    # end
    pass


def determine_transformed_template_filename(save_path):
    # transformedTemplateFileName = sprintf( '%s/mni305_masked_autoCropped_coregistered.mgz', savePath );
    pass


def find_avg_data_dir():
    # AvgDataDir = getenv( 'SAMSEG_DATA_DIR' );
    pass


def determine_mesh_collection_file_name(avg_data_dir):
    # meshCollectionFileName = sprintf('%s/CurrentMeshCollection30New.txt.gz',AvgDataDir);
    pass


def determine_compression_lookup_table_file_name(avg_data_dir):
    # compressionLookupTableFileName = sprintf('%s/namedCompressionLookupTable.txt',AvgDataDir);
    pass


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


def should_show_registration_figures(cmdargs):
    return False


def should_show_segmentation_figures(cmdargs):
    # exvivo = cmdargs.exvivo;
    # showFigures = false; % Set this to true if you want to see some figures during the run.
    # if exvivo
    #   showFigures = true;
    # end
    return is_exvivo(cmdargs)


def show_segmentation_results(names, volumesInCubicMm):
    # names
    # volumesInCubicMm
    pass
