import logging
import os

import numpy as np

logger = logging.getLogger(__name__)


def find_samseg_data_dir(samseg_data_dir=None):
    if samseg_data_dir is None:
        return os.environ.get('SAMSEG_DATA_DIR')
    return samseg_data_dir


def determine_compression_lookup_table_file_name(avg_data_dir):
    return '{0}/compressionLookupTable.txt'.format(avg_data_dir)


def find_or_create_save_path(recipe, makedirs=os.makedirs):
    save_path = recipe.output
    makedirs(save_path, exist_ok=True)
    return save_path


class Specification:
    def __init__(self, params):
        for key, value in params.items():
            self.__setattr__(key, value)


class ResolutionSpecification(Specification):
    pass


class OptimizationOptions(Specification):
    pass


class ModelSpecification(Specification):
    pass


def determine_optimization_options(verbose=False, avg_data_dir=None):
    avg_data_dir = find_samseg_data_dir(avg_data_dir)
    return OptimizationOptions({
        'multiResolutionSpecification': [
            # % Set various optimization options
            # optimizationOptions = struct;
            # optimizationOptions.multiResolutionSpecification = struct;
            # optimizationOptions.multiResolutionSpecification(1).atlasFileName = fullfile(samsegDataDir, 'atlas_level1.txt.gz');
            # optimizationOptions.multiResolutionSpecification( 1 ).targetDownsampledVoxelSpacing = 2.0; % In mm
            # optimizationOptions.multiResolutionSpecification( 1 ).maximumNumberOfIterations = 100;
            # optimizationOptions.multiResolutionSpecification( 1 ).estimateBiasField = true;
            ResolutionSpecification({
                'atlasFileName': os.path.join(avg_data_dir, 'atlas_level1.txt.gz'),
                'targetDownsampledVoxelSpacing': 2.0,
                'maximumNumberOfIterations': 100,
                'estimateBiasField': True,
            }),
            # optimizationOptions.multiResolutionSpecification( 2 ).atlasFileName = fullfile( samsegDataDir, 'atlas_level2.txt.gz' );
            # optimizationOptions.multiResolutionSpecification( 2 ).targetDownsampledVoxelSpacing = 1.0; % In mm
            # optimizationOptions.multiResolutionSpecification( 2 ).maximumNumberOfIterations = 100;
            # optimizationOptions.multiResolutionSpecification( 2 ).estimateBiasField = true; % Switching this off will use the bias field estimated
            #                                                                                 % at lower resolution(s)
            ResolutionSpecification({
                'atlasFileName': os.path.join(avg_data_dir, 'atlas_level2.txt.gz'),
                'targetDownsampledVoxelSpacing': 1.0,
                'maximumNumberOfIterations': 100,
                'estimateBiasField': True,
            }),
        ],
        # optimizationOptions.maximumNumberOfDeformationIterations = 20;
        'maximumNumberOfDeformationIterations': 20,
        # optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion = 1e-4;
        'absoluteCostPerVoxelDecreaseStopCriterion': 1e-4,
        # optimizationOptions.verbose = 0;
        'verbose': verbose,
        # optimizationOptions.maximalDeformationStopCriterion = 0.001; % Measured in pixels
        'maximalDeformationStopCriterion': 0.001,
        # optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion = optimizationOptions.maximalDeformationStopCriterion; % Idem
        'lineSearchMaximalDeformationIntervalStopCriterion': 0.001,
        # % optimizationOptions.relativeCostDecreaseStopCriterion = 1e-6;
        # optimizationOptions.maximalDeformationAppliedStopCriterion = 0.0;
        'maximalDeformationAppliedStopCriterion': 0.0,
        # optimizationOptions.BFGSMaximumMemoryLength = 12;
        'BFGSMaximumMemoryLength': 12,
    })


def specify_model(FreeSurferLabels, noBrainMasking, useDiagonalCovarianceMatrices, shared_gmm_parameters, names, colors,
                  avg_data_dir=None):
    avg_data_dir = find_samseg_data_dir(avg_data_dir)
    return ModelSpecification({
        'FreeSurferLabels': FreeSurferLabels,
        # modelSpecifications.atlasFileName = fullfile( samsegDataDir, 'atlas_level2.txt.gz' );
        'atlasFileName': os.path.join(avg_data_dir, 'atlas_level2.txt.gz'),
        # modelSpecifications.FreeSurferLabels = FreeSurferLabels;
        # modelSpecifications.names = names;
        'names': names,
        # modelSpecifications.colors = colors;
        'colors': colors,
        # modelSpecifications.sharedGMMParameters = sharedGMMParameters;
        'sharedGMMParameters': shared_gmm_parameters,
        # modelSpecifications.useDiagonalCovarianceMatrices = useDiagonalCovarianceMatrices;
        'useDiagonalCovarianceMatrices': useDiagonalCovarianceMatrices,
        # modelSpecifications.brainMaskingSmoothingSigma = 3; % sqrt of the variance of a Gaussian blurring kernel
        'brainMaskingSmoothingSigma': 3.0,
        # if noBrainMasking
        #   modelSpecifications.brainMaskingThreshold = -Inf;
        # else
        #   modelSpecifications.brainMaskingThreshold = 0.01;
        'brainMaskingThreshold': -np.inf if noBrainMasking else 0.01,
        # end
        # modelSpecifications.K = 0.1; % Stiffness of the mesh
        'K': 0.1,
        # modelSpecifications.biasFieldSmoothingKernelSize = 50.0;  % Distance in mm of sinc function center to first zero crossing
        'biasFieldSmoothingKernelSize': 50,
    })
