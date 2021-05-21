import numpy as np
import os
from freesurfer.samseg.Samseg import Samseg
from freesurfer.samseg.utilities import Specification
from freesurfer.samseg.SamsegUtility import *
from freesurfer.samseg import gems
import nibabel as nib
from os.path import join
import logging
from freesurfer.samseg.merge_alphas import kvlMergeAlphas, kvlGetMergingFractionsTable
import pickle
from skimage.transform import resize
from freesurfer import samseg
from scipy.io import loadmat
from scipy.ndimage.interpolation import affine_transform
import itertools
from scipy.ndimage.morphology import binary_fill_holes
from sklearn.metrics import log_loss
import pandas as pd
from freesurfer.samseg.io import kvlReadCompressionLookupTable, kvlReadSharedGMMParameters, GMMparameter
from scipy.ndimage import map_coordinates
import scipy
from shutil import copyfile
import tensorflow as tf
import tensorflow_probability as tfp


eps = np.finfo(float).eps


class SamsegTumor(Samseg):
    def __init__(self, imageFileNames, atlasDir, savePath, userModelSpecifications={}, userOptimizationOptions={},
                 imageToImageTransformMatrix=None, visualizer=None, saveHistory=None, savePosteriors=None,
                 saveWarp=None, saveMesh=None, threshold=None, thresholdSearchString=None,
                 targetIntensity=None, targetSearchStrings=None, modeNames=None, pallidumAsWM=True,
                 saveModelProbabilities=False, nTumorComp=2, betaMRF=1.5,useMeanConstr=True,useMRF=True,
                 saveDeformed=False,saveCorrectedImages=False, numberOfSamplingSteps=50, numberOfBurnInSteps=20,useShapeModel=True,
                 groundTruthFileName=None
                 ):
        Samseg.__init__(self, imageFileNames, atlasDir, savePath, userModelSpecifications, userOptimizationOptions,
                 imageToImageTransformMatrix, visualizer, saveHistory, savePosteriors,
                 saveWarp, saveMesh, threshold, thresholdSearchString,
                 targetIntensity, targetSearchStrings, modeNames, pallidumAsWM=pallidumAsWM,
                 saveModelProbabilities=saveModelProbabilities)

        self.numberOfSamplingSteps = numberOfSamplingSteps
        self.numberOfBurnInSteps = numberOfBurnInSteps
        self.nTumorComp = nTumorComp
        self.betaMRF = betaMRF
        self.useMRF = useMRF
        self.useMeanConstr = useMeanConstr
        self.saveDeformed = saveDeformed
        self.saveCorrectedImages = saveCorrectedImages
        self.tumorPriorMRF = None
        self.modelSpecifications = setTumorComp(self.nTumorComp, self.atlasDir, self.modelSpecifications)
        self.useShapeModel = useShapeModel
        self.groundTruthFileName = groundTruthFileName

        if len(self.imageFileNames)>1:
            self.pallidumAsWM = False

    def getClassNumber(self, structureSearchString):
        #
        if structureSearchString is None:
            return None

        #
        structureClassNumber = None
        for classNumber, mergeOption in enumerate(self.modelSpecifications.sharedGMMParameters):
            for searchString in mergeOption.searchStrings:
                if structureSearchString in searchString:
                    structureClassNumber = classNumber

        if structureClassNumber is None:
            raise RuntimeError('Could not find "%s" in model. Make sure you are using the correct atlas' % structureSearchString)

        return structureClassNumber

    def checkConditions(self, searchString, checkStructureOwnClass=True):

        # The implementation here only works for the special scenario where
        #   (a) The structure of interest has its own class (mixture model) not shared with any other structure
        #       Not checked if checkStructureOwnClass=False
        #   (b) This class (mixture model) has only a single component
        #   (c) The structure of interest is not a mixture of two or more classes (mixture models)
        # Let's test for that here

        # Get class number
        classNumber = self.getClassNumber(searchString)

        # Get class fractions
        numberOfGaussiansPerClass = [param.numberOfComponents for param in self.modelSpecifications.sharedGMMParameters]

        classFractions, _ = kvlGetMergingFractionsTable(self.modelSpecifications.names,
                                                        self.modelSpecifications.sharedGMMParameters)

        structureNumbers = np.flatnonzero(classFractions[classNumber, :] == 1)
        gaussianNumbers = [sum(numberOfGaussiansPerClass[0: classNumber])]
        if checkStructureOwnClass and structureNumbers.size is not 1:
            raise Exception('Structure of interest should correspond to exactly one class (mixture model) and vice versa')
        if len(gaussianNumbers) is not 1:
            raise Exception('Structure of interest should have a mixture model with only a single component')

        return structureNumbers[0], classNumber, gaussianNumbers[0]

    def initializeGMM(self):
        Samseg.initializeGMM(self)
        data = self.imageBuffers[self.mask, :]
        mesh = self.probabilisticAtlas.getMesh(join(self.atlasDir, "atlas_level2.txt.gz"))
        priors = mesh.rasterize_2(self.imageBuffers.shape[:3])
        maxVal = priors.max()
        priors = priors[self.mask, :] / maxVal
        priors = priors @ self.classFractions.T
        self.gmm.initializeGMMParameters(data, priors)
        tumorMeans = np.zeros_like(self.gmm.means[-self.nTumorComp:])
        dataMeans = data.mean(axis=0)
        dataStd = data.std(axis=0)
        order = self.modeNames
        for i in range(len(order)):  # table 4 mikaels paper
            if order[i] == "t1":
                tumorMeans[0, i] = dataMeans[i] + 0.2 * dataStd[i]  # edema
                tumorMeans[1, i] = dataMeans[i] + 0.2 * dataStd[i]  # core
            if order[i] == "t2":
                tumorMeans[0, i] = dataMeans[i] + 0.7 * dataStd[i]
                tumorMeans[1, i] = dataMeans[i] + 0.7 * dataStd[i]
            if order[i] == "flair":
                tumorMeans[0, i] = dataMeans[i] + 1.0 * dataStd[i]
                tumorMeans[1, i] = dataMeans[i] + 1.0 * dataStd[i]
            if order[i] in ["t1c", "t1ce", "t1gd"]:
                tumorMeans[0, i] = dataMeans[i] + 0.2 * dataStd[i]
                tumorMeans[1, i] = dataMeans[i] + 1.5 * dataStd[i]
        self.gmm.means[-self.nTumorComp:] = tumorMeans


    def initializeMeanConstraints(self):
        if not self.useMeanConstr:
            self.meanConstraints = [None,None]
            return
        self.wm_comp_number = self.getClassNumber("White")
        self.gm_comp_number = self.getClassNumber("Cortex")
        self.ed_comp_number = self.getClassNumber("Tumor")
        self.enh_comp_number = self.ed_comp_number + 1
        n_contrasts = self.imageBuffers.shape[-1]
        n_constraints = 6
        A = np.zeros((n_constraints, n_contrasts * self.gmm.numberOfGaussians))
        b = np.zeros(n_constraints)
        modeNames = [i.lower() for i in self.modeNames]
        if "flair" in modeNames:
            self.flair_idx = modeNames.index("flair")
            A[0, self.ed_comp_number * n_contrasts + self.flair_idx] = -1
            A[0, self.wm_comp_number * n_contrasts + self.flair_idx] = 1
            b[0] = 1.15
            A[1, self.ed_comp_number * n_contrasts + self.flair_idx] = -1
            A[1, self.gm_comp_number * n_contrasts + self.flair_idx] = 1
            b[1] = 1.15
            A[2, self.enh_comp_number * n_contrasts + self.flair_idx] = -1
            A[2, self.wm_comp_number * n_contrasts + self.flair_idx] = 1
            b[2] = 1.0
            A[3, self.enh_comp_number * n_contrasts + self.flair_idx] = -1
            A[3, self.gm_comp_number * n_contrasts + self.flair_idx] = 1
            b[3] = 1.0
        else:
            self.flair_idx = None
            A = A[4:]
            b = b[4:]
        if "t1ce" in modeNames or "t1gd" in modeNames or "t1c" in modeNames:
            self.t1ce_idx = modeNames.index([i for i in modeNames if i in ["t1ce","t1gd","t1c"]][0])
            A[4, self.enh_comp_number * n_contrasts + self.t1ce_idx] = -1
            A[4, self.wm_comp_number * n_contrasts + self.t1ce_idx] = 1
            b[4] = 1.1
            A[5, self.enh_comp_number * n_contrasts + self.t1ce_idx] = -1
            A[5, self.gm_comp_number * n_contrasts + self.t1ce_idx] = 1
            b[5] = 1.1
        else:
            self.t1ce_idx = None
            if not "flair" in modeNames:
                raise ValueError("Cannot initialize constraints on means. Neither flair or t1ce modes found.")
            A = A[:4]
            b = b[:4]
        b = -np.log(b)
        self.meanConstraints = [A,b]


    def estimateModelParameters( self, initialBiasFieldCoefficients=None, initialDeformation=None,
                                 initialDeformationAtlasFileName=None,
                                 skipGMMParameterEstimationInFirstIteration=False,
                                 skipBiasFieldParameterEstimationInFirstIteration=True):

        logger = logging.getLogger(__name__)
        self.optimizationHistory = []
        self.optimizationSummary = []
        if self.useMeanConstr:
            self.initializeMeanConstraints()

        # Convert optimizationOptions from dictionary into something more convenient to access
        optimizationOptions = Specification(self.optimizationOptions)
        source = optimizationOptions.multiResolutionSpecification
        optimizationOptions.multiResolutionSpecification = []
        for levelNumber in range(len(source)):
            optimizationOptions.multiResolutionSpecification.append(Specification(source[levelNumber]))
        print('====================')
        print(optimizationOptions)
        print('====================')

        self.deformation, self.deformationAtlasFileName = initialDeformation, initialDeformationAtlasFileName
        self.biasField.setBiasFieldCoefficients(initialBiasFieldCoefficients)

        # Loop over resolution levels
        numberOfMultiResolutionLevels = len(optimizationOptions.multiResolutionSpecification)
        numberOfGaussiansPerClass = [param.numberOfComponents for param in self.modelSpecifications.sharedGMMParameters]

        for multiResolutionLevel in range(numberOfMultiResolutionLevels):

            logger.debug('multiResolutionLevel=%d', multiResolutionLevel)
            self.visualizer.start_movie(window_id='Mesh deformation (level ' + str(multiResolutionLevel) + ')',
                                        title='Mesh Deformation - the movie (level ' + str(multiResolutionLevel) + ')')

            maximumNumberOfIterations = optimizationOptions.multiResolutionSpecification[
                multiResolutionLevel].maximumNumberOfIterations
            estimateBiasField = optimizationOptions.multiResolutionSpecification[multiResolutionLevel].estimateBiasField
            historyOfCost = [1 / eps]
            logger.debug('maximumNumberOfIterations: %d', maximumNumberOfIterations)

            # Downsample the images, the mask, the mesh, and the bias field basis functions (integer)
            logger.debug('Setting up downsampled model')
            downSamplingFactors = np.uint32(np.round(optimizationOptions.multiResolutionSpecification
                                                     [
                                                         multiResolutionLevel].targetDownsampledVoxelSpacing / self.voxelSpacing))
            downSamplingFactors[downSamplingFactors < 1] = 1
            downSampledImageBuffers, downSampledMask, downSampledMesh, downSampledInitialDeformationApplied, \
            downSampledTransform = self.getDownSampledModel(
                optimizationOptions.multiResolutionSpecification[multiResolutionLevel].atlasFileName,
                downSamplingFactors)

            self.biasField.downSampleBasisFunctions(downSamplingFactors)

            # Also downsample the strength of the hyperprior, if any
            self.gmm.downsampledHyperparameters(downSamplingFactors)

            # Save initial position at the start of this multi-resolution level
            initialNodePositions = downSampledMesh.points
            initialNodePositionsInTemplateSpace = self.probabilisticAtlas.mapPositionsFromSubjectToTemplateSpace(
                initialNodePositions, downSampledTransform)

            # Set priors in mesh to the merged (super-structure) ones
            mergedAlphas = kvlMergeAlphas(downSampledMesh.alphas, self.classFractions)
            downSampledMesh.alphas = mergedAlphas

            #
            cropping = self.cropping
            mask = self.mask
            orig_shape = np.array(nib.load(self.imageFileNames[0]).shape).astype(np.int)
            uncropped_mask = np.zeros(orig_shape)
            uncropped_mask[cropping[0].start:cropping[0].stop, cropping[1].start:cropping[1].stop,
            cropping[2].start:cropping[2].stop] = mask
            uncropped_mask = uncropped_mask.astype(np.bool)
            uncropped_mask_ds = uncropped_mask[
                                cropping[0].start % downSamplingFactors[0]::downSamplingFactors[0],
                                cropping[1].start % downSamplingFactors[1]::downSamplingFactors[1],
                                cropping[2].start % downSamplingFactors[2]::downSamplingFactors[2]]
            orig_shape_ds = uncropped_mask_ds.shape
            self.visualizer.show(mesh=downSampledMesh, images=downSampledImageBuffers,
                                 window_id='Mesh deformation (level ' + str(multiResolutionLevel) + ')',
                                 title='Mesh Deformation (level ' + str(multiResolutionLevel) + ')')
            if self.saveHistory:
                levelHistory = {'historyWithinEachIteration': []}

            if self.useMRF:
                """Pre-compute some arrays that will be used for fast indexing for MRF"""
                qshape = orig_shape_ds
                neighbors = getNeighborIndices(qshape)
                XX, YY, ZZ = np.meshgrid(np.arange(qshape[0]), np.arange(qshape[1]), np.arange(qshape[2]), indexing='ij')
                x_parts = [XX[::2, ::2, ::2], XX[::2, ::2, 1::2], XX[::2, 1::2, ::2], XX[::2, 1::2, 1::2],
                           XX[1::2, ::2, ::2], XX[1::2, ::2, 1::2], XX[1::2, 1::2, ::2], XX[1::2, 1::2, 1::2]]
                y_parts = [YY[::2, ::2, ::2], YY[::2, ::2, 1::2], YY[::2, 1::2, ::2], YY[::2, 1::2, 1::2],
                           YY[1::2, ::2, ::2], YY[1::2, ::2, 1::2], YY[1::2, 1::2, ::2], YY[1::2, 1::2, 1::2]]
                z_parts = [ZZ[::2, ::2, ::2], ZZ[::2, ::2, 1::2], ZZ[::2, 1::2, ::2], ZZ[::2, 1::2, 1::2],
                           ZZ[1::2, ::2, ::2], ZZ[1::2, ::2, 1::2], ZZ[1::2, 1::2, ::2], ZZ[1::2, 1::2, 1::2]]

                neighbor_parts = [neighbors[::2, ::2, ::2], neighbors[::2, ::2, 1::2], neighbors[::2, 1::2, ::2],
                                  neighbors[::2, 1::2, 1::2],
                                  neighbors[1::2, ::2, ::2], neighbors[1::2, ::2, 1::2], neighbors[1::2, 1::2, ::2],
                                  neighbors[1::2, 1::2, 1::2]]
            imlist = []
            for iterationNumber in range(maximumNumberOfIterations):
                logger.debug('iterationNumber=%d', iterationNumber)

                # Part I: estimate Gaussian mixture model parameters, as well as bias field parameters using EM.

                # Get the priors at the current mesh position
                tmp = downSampledMesh.rasterize_2(downSampledMask.shape, -1)
                downSampledClassPriors = tmp[downSampledMask] / 65535

                # Initialize the model parameters if needed
                if self.gmm.means is None:
                    self.gmm.initializeGMMParameters(downSampledImageBuffers[downSampledMask, :],
                                                     downSampledClassPriors)

                if self.biasField.coefficients is None:
                    numberOfBasisFunctions = [functions.shape[1] for functions in self.biasField.basisFunctions]
                    numberOfContrasts = downSampledImageBuffers.shape[-1]
                    initialBiasFieldCoefficients = np.zeros((np.prod(numberOfBasisFunctions), numberOfContrasts))
                    self.biasField.setBiasFieldCoefficients(initialBiasFieldCoefficients)

                # Start EM iterations
                historyOfEMCost = [1 / eps]
                EMIterationNumber = 0

                while True:
                    logger.debug('EMIterationNumber=%d', EMIterationNumber)

                    # Precompute intensities after bias field correction for later use (really only caching something that
                    # doesn't really figure in the model
                    downSampledBiasFields = self.biasField.getBiasFields(downSampledMask)
                    downSampledData = downSampledImageBuffers[downSampledMask, :] - downSampledBiasFields[
                                                                                    downSampledMask, :]
                    self.visualizer.show(image_list=[downSampledBiasFields[..., i]
                                                     for i in range(downSampledBiasFields.shape[-1])],
                                         auto_scale=True, window_id='bias field', title='Bias Fields')

                    # E-step: compute the downSampledGaussianPosteriors based on the current parameters
                    if self.useMRF and (EMIterationNumber > 0 or iterationNumber > 0):
                        likelihoods = self.gmm.getLikelihoods(downSampledData, self.classFractions)
                        likelihoods = np.dot(likelihoods, self.classFractions.T)
                        mergetable = downSampledClassPriors[:, :-1] / np.expand_dims(
                            downSampledClassPriors[:, :-1].sum(axis=-1), 1)
                        likelihood_0 = np.sum(likelihoods[:, :-1] * mergetable, axis=-1)
                        likelihood_1 = likelihoods[:, -1]
                        likelihood_ = np.zeros([2] + list(uncropped_mask_ds.shape))
                        likelihood_[0][uncropped_mask_ds] = likelihood_0
                        likelihood_[1][uncropped_mask_ds] = likelihood_1
                        likelihood_[0][~uncropped_mask_ds] = 1

                        q = np.zeros_like(uncropped_mask_ds).astype(np.float)  # q is probability of tumor class
                        q1 = q.copy()
                        q0 = q.copy()
                        q[uncropped_mask_ds] = downSampledGaussianPosteriors[:, -self.nTumorComp:].sum(axis=1)

                        beta = self.betaMRF
                        if multiResolutionLevel==0 and numberOfMultiResolutionLevels>1:
                            beta /= 8
                        tumor_prior_prob = 0.5
                        prior_ = [1 - tumor_prior_prob, tumor_prior_prob]  # uniform prior for tumor class
                        n_iter = 10
                        for iteration in range(n_iter):
                            for j in range(8):
                                q_flat = q.flatten()
                                neighbors_q = q_flat[neighbor_parts[j]]
                                gamma_1 = np.exp(-beta * np.sum(1 - neighbors_q, axis=-1)) * prior_[1]
                                gamma_0 = np.exp(-beta * np.sum(neighbors_q, axis=-1)) * prior_[0]
                                gamma_den = (gamma_1 + gamma_0)

                                gamma_1 /= gamma_den
                                gamma_0 /= gamma_den

                                new_q1 = likelihood_[1][x_parts[j], y_parts[j], z_parts[j]] * gamma_1 + eps
                                new_q0 = likelihood_[0][x_parts[j], y_parts[j], z_parts[j]] * gamma_0 + eps
                                new_q_den = new_q1 + new_q0
                                new_q = new_q1 / (new_q_den)
                                q0[x_parts[j], y_parts[j], z_parts[j]] = new_q0
                                q1[x_parts[j], y_parts[j], z_parts[j]] = new_q1
                                q[x_parts[j], y_parts[j], z_parts[j]] = new_q

                        llh_mrf = np.log(q1[uncropped_mask_ds] + q0[uncropped_mask_ds]).sum() #TODO: Make sure this is right

                        mrf_class_priors = np.concatenate(
                            [(1 - q[uncropped_mask_ds][:, np.newaxis]) * downSampledClassPriors[:, :-1] + 1e-16,
                             np.expand_dims(q[uncropped_mask_ds] * downSampledClassPriors[:, -1] + 1e-16, 1)],
                            axis=1)
                        mrf_class_priors /= np.expand_dims(mrf_class_priors.sum(axis=1), 1)
                        downSampledGaussianPosteriors, minLogLikelihood = self.gmm.getGaussianPosteriors(
                            downSampledData,
                            mrf_class_priors)
                        minLogLikelihood = -llh_mrf
                        self.tumorq = q[uncropped_mask_ds]
                    else:
                        downSampledGaussianPosteriors, minLogLikelihood = self.gmm.getGaussianPosteriors(
                            downSampledData,
                            downSampledClassPriors)

                    # Compute the log-posterior of the model parameters, and check for convergence
                    minLogGMMParametersPrior = self.gmm.evaluateMinLogPriorOfGMMParameters()

                    historyOfEMCost.append(minLogLikelihood + minLogGMMParametersPrior)

                    self.visualizer.plot(historyOfEMCost[1:], window_id='history of EM cost',
                                         title='History of EM Cost (level: ' + str(multiResolutionLevel) +
                                               ' iteration: ' + str(iterationNumber) + ')')
                    EMIterationNumber += 1
                    changeCostEMPerVoxel = (historyOfEMCost[-2] - historyOfEMCost[-1]) / downSampledData.shape[0]
                    changeCostEMPerVoxelThreshold = optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion

                    tumor_seg = np.zeros(downSampledMask.shape)
                    is_tumor = downSampledGaussianPosteriors[:, -self.nTumorComp:].sum(axis=1) > 0.5
                    tumor_mask = np.zeros(downSampledMask.shape)
                    tumor_mask[downSampledMask] = is_tumor
                    tumor_mask = tumor_mask.astype(np.bool)
                    tumor_seg[tumor_mask] = np.argmax(downSampledGaussianPosteriors[is_tumor][:, -self.nTumorComp:],
                                                      axis=1) + 1

                    if (EMIterationNumber == 100) or (changeCostEMPerVoxel < changeCostEMPerVoxelThreshold):
                        # Converged
                        print('EM converged!')
                        break

                    # M-step: update the model parameters based on the current posterior
                    #
                    # First the mixture model parameters
                    if not ((iterationNumber == 0) and skipGMMParameterEstimationInFirstIteration):
                        if self.useMeanConstr:
                            self.gmm.fitGMMParametersWithConstraints(downSampledData, downSampledGaussianPosteriors,self.meanConstraints[0],self.meanConstraints[1])
                        else:
                            self.gmm.fitGMMParameters(downSampledData, downSampledGaussianPosteriors)

                    # Now update the parameters of the bias field model.
                    if (estimateBiasField and not ((iterationNumber == 0)
                                                   and skipBiasFieldParameterEstimationInFirstIteration)):
                        self.biasField.fitBiasFieldParameters(downSampledImageBuffers, downSampledGaussianPosteriors,
                                                              self.gmm.means, self.gmm.variances, downSampledMask)
                    # End test if bias field update

                # End loop over EM iterations
                historyOfEMCost = historyOfEMCost[1:]

                # Visualize the posteriors
                if hasattr(self.visualizer, 'show_flag'):
                    tmp = np.zeros(downSampledMask.shape + (downSampledGaussianPosteriors.shape[-1],))
                    tmp[downSampledMask, :] = downSampledGaussianPosteriors
                    self.visualizer.show(probabilities=tmp, images=downSampledImageBuffers,
                                         window_id='EM Gaussian posteriors',
                                         title='EM Gaussian posteriors (level: ' + str(multiResolutionLevel) +
                                               ' iteration: ' + str(iterationNumber) + ')')

                # Part II: update the position of the mesh nodes for the current mixture model and bias field parameter estimates
                optimizationParameters = {
                    'Verbose': optimizationOptions.verbose,
                    'MaximalDeformationStopCriterion': optimizationOptions.maximalDeformationStopCriterion,
                    'LineSearchMaximalDeformationIntervalStopCriterion': optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion,
                    'MaximumNumberOfIterations': optimizationOptions.maximumNumberOfDeformationIterations,
                    'BFGS-MaximumMemoryLength': optimizationOptions.BFGSMaximumMemoryLength
                }
                historyOfDeformationCost, historyOfMaximalDeformation, maximalDeformationApplied, minLogLikelihoodTimesDeformationPrior = \
                    self.probabilisticAtlas.deformMesh(downSampledMesh, downSampledTransform, downSampledData,
                                                       downSampledMask,
                                                       self.gmm.means, self.gmm.variances, self.gmm.mixtureWeights,
                                                       self.gmm.numberOfGaussiansPerClass, optimizationParameters)

                # print summary of iteration
                print('iterationNumber: %d' % iterationNumber)
                print('maximalDeformationApplied: %.4f' % maximalDeformationApplied)
                print('=======================================================')
                self.visualizer.show(mesh=downSampledMesh, images=downSampledImageBuffers,
                                     window_id='Mesh deformation (level ' + str(multiResolutionLevel) + ')',
                                     title='Mesh Deformation (level ' + str(multiResolutionLevel) + ')')

                # Save history of the estimation
                if self.saveHistory:
                    levelHistory['historyWithinEachIteration'].append({
                        'historyOfEMCost': historyOfEMCost,
                        'mixtureWeights': self.gmm.mixtureWeights,
                        'means': self.gmm.means,
                        'variances': self.gmm.variances,
                        'biasFieldCoefficients': self.biasField.coefficients,
                        'historyOfDeformationCost': historyOfDeformationCost,
                        'historyOfMaximalDeformation': historyOfMaximalDeformation,
                        'maximalDeformationApplied': maximalDeformationApplied
                    })

                # Check for convergence
                historyOfCost.append(minLogLikelihoodTimesDeformationPrior + minLogGMMParametersPrior)
                self.visualizer.plot(historyOfCost[1:],
                                     window_id='history of cost (level ' + str(multiResolutionLevel) + ')',
                                     title='History of Cost (level ' + str(multiResolutionLevel) + ')')
                previousCost = historyOfCost[-2]
                currentCost = historyOfCost[-1]
                costChange = previousCost - currentCost
                perVoxelDecrease = costChange / np.count_nonzero(downSampledMask)
                perVoxelDecreaseThreshold = optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion
                if perVoxelDecrease < perVoxelDecreaseThreshold:
                    break

            if self.useMRF:
                downSamplingFactors = np.uint32(np.round(optimizationOptions.multiResolutionSpecification
                                                         [-1].targetDownsampledVoxelSpacing / self.voxelSpacing))
                downSamplingFactors[downSamplingFactors < 1] = 1
                if (downSamplingFactors>1).any():
                    q = resize(q, orig_shape, order=3)
                    self.tumorPriorMRF=q[uncropped_mask]

            #TODO upsample q, set self.tumorq=q
            # End loop over coordinate descent optimization (intensity model parameters vs. atlas deformation)

            # Visualize the mesh deformation across iterations
            self.visualizer.show_movie(window_id='Mesh deformation (level ' + str(multiResolutionLevel) + ')')

            # Log the final per-voxel cost
            self.optimizationSummary.append({'numberOfIterations': iterationNumber + 1,
                                             'perVoxelCost': currentCost / np.count_nonzero(downSampledMask)})

            # Get the final node positions
            finalNodePositions = downSampledMesh.points

            # Transform back in template space (i.e., undoing the affine registration that we applied), and save for later usage
            finalNodePositionsInTemplateSpace = self.probabilisticAtlas.mapPositionsFromSubjectToTemplateSpace(
                finalNodePositions, downSampledTransform)

            # Record deformation delta here in lieu of maintaining history
            self.deformation = finalNodePositionsInTemplateSpace - initialNodePositionsInTemplateSpace + downSampledInitialDeformationApplied
            self.deformationAtlasFileName = optimizationOptions.multiResolutionSpecification[
                multiResolutionLevel].atlasFileName

            # Save history of the estimation
            if self.saveHistory:
                levelHistory['downSamplingFactors'] = downSamplingFactors
                levelHistory['downSampledImageBuffers'] = downSampledImageBuffers
                levelHistory['downSampledMask'] = downSampledMask
                levelHistory['downSampledTransformMatrix'] = downSampledTransform.as_numpy_array
                levelHistory['initialNodePositions'] = initialNodePositions
                levelHistory['finalNodePositions'] = finalNodePositions
                levelHistory['initialNodePositionsInTemplateSpace'] = initialNodePositionsInTemplateSpace
                levelHistory['finalNodePositionsInTemplateSpace'] = finalNodePositionsInTemplateSpace
                levelHistory['historyOfCost'] = historyOfCost
                levelHistory['priorsAtEnd'] = downSampledClassPriors
                levelHistory['posteriorsAtEnd'] = downSampledGaussianPosteriors
                self.optimizationHistory.append(levelHistory)



    def postProcess(self):
        # =======================================================================================
        #
        # Segment the data using the estimate model parameters, and write results out
        #
        # =======================================================================================

        # OK, now that all the parameters have been estimated, try to segment the original, full resolution image
        # with all the original labels (instead of the reduced "super"-structure labels we created)
        if self.useShapeModel:
            posteriors, biasFields, nodePositions,data, priors, posteriorsTumor = self.computeFinalSegmentationShapeModel()
        else:
            posteriors, biasFields, nodePositions, data, priors, atlasPriors, posteriorsTumor = self.computeFinalSegmentation()

        # Write out segmentation and bias field corrected volumes
        volumesInCubicMm = self.writeResults(biasFields, posteriors, posteriorsTumor)

        # Save the template warp
        if self.saveWarp:
            print('Saving the template warp')
            self.saveWarpField(os.path.join(self.savePath, 'template.m3z'))

        # Save the final mesh collection
        if self.saveMesh:
            print('Saving the final mesh in template space')
            deformedAtlasFileName = os.path.join(self.savePath, 'mesh.txt')
            self.probabilisticAtlas.saveDeformedAtlas(self.modelSpecifications.atlasFileName, deformedAtlasFileName, nodePositions)

        # Save the history of the parameter estimation process
        if self.saveHistory:
            history = {'input': {
                'imageFileNames': self.imageFileNames,
                'imageToImageTransformMatrix': self.imageToImageTransformMatrix,
                'modelSpecifications': self.modelSpecifications,
                'optimizationOptions': self.optimizationOptions,
                'savePath': self.savePath
            }, 'imageBuffers': self.imageBuffers, 'mask': self.mask,
                'historyWithinEachMultiResolutionLevel': self.optimizationHistory,
                "labels": self.modelSpecifications.FreeSurferLabels, "names": self.modelSpecifications.names,
                "volumesInCubicMm": volumesInCubicMm, "optimizationSummary": self.optimizationSummary}
            with open(os.path.join(self.savePath, 'history.p'), 'wb') as file:
                pickle.dump(history, file, protocol=pickle.HIGHEST_PROTOCOL)

        return self.modelSpecifications.FreeSurferLabels, self.modelSpecifications.names, volumesInCubicMm, self.optimizationSummary


    def writeResults(self, biasFields, posteriors, posteriorsTumor):
        # Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention
        names = self.modelSpecifications.names
        structureNumbers = np.array(np.argmax(posteriors, 1), dtype=np.uint32)

        segmentation = np.zeros(self.imageBuffers.shape[:3], dtype=np.int32)
        fslabels = np.array(self.modelSpecifications.FreeSurferLabels, dtype=np.int32)
        if self.useShapeModel:
            fslabels = np.array(list(fslabels[:-1]) + [89,90,91])
            names = names[:-1] + ["Core","Edema","EnhCore"]
        segmentation[self.mask] = fslabels[structureNumbers]
        scalingFactors = scaleBiasFields(biasFields, self.imageBuffers, self.mask, posteriors, self.targetIntensity, self.targetSearchStrings, names)
        # Get corrected intensities and bias field images in the non-log transformed domain
        expImageBuffers, expBiasFields = undoLogTransformAndBiasField(self.imageBuffers, biasFields, self.mask)

        if self.saveDeformed:
            """Transform the posteriors to the template space and save results"""
            deformedPath = os.path.join(self.savePath,"templateSpaceResults")
            os.makedirs(deformedPath,exist_ok=True)
            referenceMeshFileName = join(self.atlasDir,"atlas_level2.txt.gz")
            atlas = samseg.ProbabilisticAtlas()
            referenceMesh = atlas.getMesh(referenceMeshFileName)
            posteriorsTmp = np.zeros(list(self.imageBuffers.shape[:3])+[posteriors.shape[-1]],dtype=np.float32)
            posteriorsTmp[self.mask] = posteriors
            posteriorsTemplateSpace = imagesInAtlasSpace(posteriorsTmp, self.deformation, self.transform,
                                                           referenceMesh)
            posteriorsTemplateSpaceHard = fslabels[np.argmax(posteriorsTemplateSpace,axis=-1).ravel().astype(np.int32)].astype(np.int32)
            posteriorsTemplateSpaceHard = posteriorsTemplateSpaceHard.reshape(posteriorsTemplateSpace.shape[:3])
            self.writeDeformedImage(posteriorsTemplateSpaceHard, os.path.join(deformedPath, 'seg.mgz'), saveLabels=True)
            posteriorsTumorTmp = np.zeros(list(self.imageBuffers.shape[:3])+[posteriorsTumor.shape[-1]],dtype=np.float32)
            posteriorsTumorTmp[self.mask] = posteriorsTumor
            posteriorsTumorTemplateSpace = imagesInAtlasSpace(posteriorsTumorTmp, self.deformation, self.transform,
                                                           referenceMesh)
        # Write out various images - segmentation first
        self.writeImage(segmentation, os.path.join(self.savePath, 'seg.mgz'), saveLabels=True)
        if self.saveCorrectedImages:
            for contrastNumber, imageFileName in enumerate(self.imageFileNames):
                # Contrast-specific filename prefix
                contastPrefix = os.path.join(self.savePath, self.modeNames[contrastNumber])

                # Write bias field and bias-corrected image
                self.writeImage(expBiasFields[..., contrastNumber],   contastPrefix + '_bias_field.mgz')
                self.writeImage(expImageBuffers[..., contrastNumber], contastPrefix + '_bias_corrected.mgz')

                # Save a note indicating the scaling factor
                with open(contastPrefix + '_scaling.txt', 'w') as fid:
                    print(scalingFactors[contrastNumber], file=fid)
        posteriorPath = os.path.join(self.savePath, 'posteriors')
        os.makedirs(posteriorPath, exist_ok=True)
        if self.saveDeformed:
            os.makedirs(join(deformedPath, "posteriors"), exist_ok=True)
        if self.savePosteriors:
            for structureNumber, name in enumerate(names):
                # Write the posteriors to seperate volume files
                for searchString in self.savePosteriors:
                    if searchString in name:
                        posteriorVol = np.zeros(self.imageBuffers.shape[:3], dtype=np.float32)
                        posteriorVol[self.mask] = posteriors[:, structureNumber]
                        self.writeImage(posteriorVol, os.path.join(posteriorPath, name + '.mgz'))
                        if self.saveDeformed:
                            self.writeDeformedImage(posteriorsTemplateSpace[:,:,:,structureNumber], os.path.join(join(deformedPath,"posteriors"), name + '.mgz'))

        for tumorComp in range(posteriorsTumor.shape[-1]):
            name = "TumorComp%d"%tumorComp
            posteriorVol = np.zeros(self.imageBuffers.shape[:3], dtype=np.float32)
            posteriorVol[self.mask] = posteriorsTumor[:, tumorComp].astype(np.float32)
            self.writeImage(posteriorVol, os.path.join(posteriorPath, name + '.mgz'))
            if self.saveDeformed:
                self.writeDeformedImage(posteriorsTumorTemplateSpace[:, :, :, tumorComp],
                                        os.path.join(join(deformedPath, "posteriors"), name + '.mgz'))

        # Compute volumes in mm^3
        # TODO: cache the source geometry in __init__, as this is also loaded by writeImage
        exampleImage = gems.KvlImage(self.imageFileNames[0])
        volumeOfOneVoxel = np.abs(np.linalg.det(exampleImage.transform_matrix.as_numpy_array[:3, :3]))
        volumesInCubicMm = np.sum(posteriors, axis=0) * volumeOfOneVoxel

        # Write structural volumes
        with open(os.path.join(self.savePath, 'samseg.stats'), 'w') as fid:
            for volume, name in zip(volumesInCubicMm, names):
                fid.write('# Measure %s, %.6f, mm^3\n' % (name, volume))

        # Write intracranial volume
        sbtiv = icv(zip(*[names, volumesInCubicMm]))
        with open(os.path.join(self.savePath, 'sbtiv.stats'), 'w') as fid:
            fid.write('# Measure Intra-Cranial, %.6f, mm^3\n' % sbtiv)

        return volumesInCubicMm

    def writeDeformedImage(self, data, path, saveLabels=False):
        # Read source geometry
        geom = fs.Volume.read(os.path.join(self.atlasDir, 'template.nii')).geometry()
        # Uncrop image
        volume = fs.Volume(data, affine=geom.affine, voxsize=geom.voxsize)
        if saveLabels:
            volume.lut = fs.LookupTable.read_default()
        volume.write(path)

    def computeFinalSegmentation(self):
        # Get the final mesh
        mesh = self.probabilisticAtlas.getMesh(self.modelSpecifications.atlasFileName, self.transform,
                                               initialDeformation=self.deformation,
                                               initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)

        # Get the priors as dictated by the current mesh position
        priors = mesh.rasterize(self.imageBuffers.shape[0:3], -1)
        priors = priors/priors.max()
        priors = priors[self.mask, :]

        # Get bias field corrected data
        # Make sure that the bias field basis function are not downsampled
        # (this might happens if the parameters estimation is made only with one downsampled resolution)
        self.biasField.downSampleBasisFunctions([1, 1, 1])
        biasFields = self.biasField.getBiasFields()
        data = self.imageBuffers[self.mask, :] - biasFields[self.mask, :]
        # Compute the posterior distribution of the various structures
        if self.useMRF:
            priorsMRF = np.concatenate([
                (1 - self.tumorq)[:, np.newaxis] * priors[:, :-1] + 1e-16,
                np.expand_dims(self.tumorq * priors[:, -1] + 1e-16, 1)
            ], axis=1)
            priorsMRF /= np.expand_dims(priorsMRF.sum(axis=1), 1)
            atlasPriors = priors
            priors = priorsMRF
            posteriors = self.gmm.getPosteriors(data, priorsMRF, self.classFractions)
            posteriorsGaussian,_ = self.gmm.getGaussianPosteriors(data, priorsMRF@self.classFractions.T)
        else:
            posteriors = self.gmm.getPosteriors(data, priors, self.classFractions)
            posteriorsGaussian,_ = self.gmm.getGaussianPosteriors(data, priors@self.classFractions.T)
            atlasPriors = priors

        estimatedNodePositions = self.probabilisticAtlas.mapPositionsFromSubjectToTemplateSpace(mesh.points,
                                                                                                self.transform)
        return posteriors, biasFields, estimatedNodePositions, data, priors, atlasPriors, posteriorsGaussian[:,-self.nTumorComp:]


    def computeFinalSegmentationShapeModel(self):
        posteriors, biasFields, nodePositions,data, priors, atlasPriors, posteriorsTumor = self.computeFinalSegmentation()
        imageSize = self.mask.shape
        tumorClassNumber = self.getClassNumber("Tumor")
        tumorStructureNumber = np.flatnonzero(self.classFractions[tumorClassNumber, :] == 1)[0]
        trainToTemplateTransformFileName = os.path.join( self.atlasDir, 'brats_template_transforms.mat' ) #because the VAE is trained on brats data
        trainToTemplateMat = np.linalg.inv(loadmat(trainToTemplateTransformFileName)['imageToImageTransformMatrix'])
        trainToSubjectMat = self.transform.as_numpy_array @ trainToTemplateMat
        subjectToTrainMat = np.linalg.inv(trainToSubjectMat)
        trainToSubjectMat[-1] = np.array([0, 0, 0, 1])
        subjectToTrainMat[-1] = np.array([0, 0, 0, 1])

        likelihoods = self.gmm.getLikelihoods(data, self.classFractions)
        posteriors = likelihoods * priors
        posteriors /= np.expand_dims(np.sum(posteriors, axis=1) + eps, 1)

        vaeShape = np.array([240, 240, 155])
        vae = CatVAE(vaeShape[0], vaeShape[1], vaeShape[2], 4, 64, 1)
        vae.load_weights( os.path.join( self.atlasDir, 'model.h5' ) )
        if self.groundTruthFileName is not None:
            self.groundTruth = nib.load(self.groundTruthFileName).get_fdata()
            self.groundTruth[self.groundTruth==4]=3
        else:
            self.groundTruth = None

        """Next we introduce a new tumor component and fit all 3 of them to data where the currently the label is tumor"""
        isTumor = np.array(np.argmax(posteriors, 1), dtype=np.uint32) == tumorStructureNumber
        isTumorMask = np.zeros(self.mask.shape).astype(np.bool)
        isTumorMask[self.mask] = isTumor
        isTumorMask = binary_fill_holes(isTumorMask) #TODO: Decide if this is ok
        isTumorMask[~self.mask] = False
        isTumor = isTumorMask[self.mask]
        changeCostEMPerVoxelThreshold = self.optimizationOptions["absoluteCostPerVoxelDecreaseStopCriterion"]
        tumorGmm = samseg.GMM([3], numberOfContrasts=self.imageBuffers.shape[-1],
                              useDiagonalCovarianceMatrices=self.modelSpecifications.useDiagonalCovarianceMatrices)
        tumorGmm.fullHyperMeansNumberOfMeasurements = np.zeros((3))
        tumorGmm.fullHyperVariancesNumberOfMeasurements = np.ones(3) * 3
        tmpPrior = np.ones([isTumor.sum(), 3]) / 3
        tumorGmm.initializeGMMParameters(data[isTumor, :], tmpPrior)
        print("tumor size:", isTumor.sum())
        nll_hist = [1 / eps]
        for j in range(150):
            tumorPosterior, nll = tumorGmm.getGaussianPosteriors(data[isTumor], tmpPrior)
            tumorPosterior /= np.expand_dims(tumorPosterior.sum(axis=1), 1)
            tumorGmm.fitGMMParameters(data[isTumor], tumorPosterior)
            nll_hist.append(nll)
            if nll_hist[-2] - nll_hist[-1] < changeCostEMPerVoxelThreshold * isTumor.sum() / 100:
                print("Initalized new tumor component: ran %d iterations" % (j + 1))
                break

        """Now we use the VAE to select the most likely configuration of the tumor components"""
        configurations = list(itertools.permutations([1, 2, 3]))
        tumor = np.zeros(list(imageSize) + [self.nTumorComp + 2])
        for tumorCompNumber in range(self.nTumorComp + 1):
            tumor[isTumorMask, 1 + tumorCompNumber] = np.argmax(tumorPosterior, axis=-1) == tumorCompNumber
        tumor[~isTumorMask, 0] = 1  # order: Background, tumorComp1, 2 ,3, (the initial configuration)
        tumorTrainSpace = np.zeros(np.concatenate([vaeShape, [self.nTumorComp + 2]]))
        for tumorCompNumber in range(self.nTumorComp + 2):
            tumorTrainSpace[:, :, :, tumorCompNumber] = affine_transform(tumor[:, :, :, tumorCompNumber],
                                                                         trainToSubjectMat, output_shape=vaeShape,
                                                                         order=1)
        configLoss = []
        for config in configurations:
            tumorTmp = tumorTrainSpace[:, :, :, np.array([0] + list(config))]
            mean, logvar = vae.encode(np.expand_dims(tumorTmp, 0))
            z = vae.reparameterize(mean, logvar)
            tumorVAE = vae.decode(z, apply_softmax=True).numpy()[0]
            vaeLoss = log_loss(np.argmax(tumorTmp.reshape(-1, tumor.shape[-1]), axis=1),
                               tumorVAE.reshape(-1, tumorVAE.shape[-1]))
            configLoss.append(vaeLoss)

        bestConfig = np.array(configurations[np.argmin(configLoss)])
        tumor = tumor[:, :, :, np.array([0] + list(
            bestConfig))]  # now the order is assumed to be the same as in the VAE: Background, core, edema, enhancing core
        tumorGmm.means = tumorGmm.means[bestConfig - 1]
        tumorGmm.variances = tumorGmm.variances[bestConfig - 1]
        tumorGmm.mixtureWeights = tumorGmm.mixtureWeights[bestConfig - 1]
        # create a new GMM with 1 class per tumor comp -- In the end, this new gmm, replaces self.gmm
        gmm = samseg.GMM(self.gmm.numberOfGaussiansPerClass[:-1] + list(np.ones(self.nTumorComp + 1).astype(np.int)),
                         numberOfContrasts=self.imageBuffers.shape[-1],
                         useDiagonalCovarianceMatrices=self.modelSpecifications.useDiagonalCovarianceMatrices)
        gmm.initializeGMMParameters(data, np.ones([data.shape[0], gmm.numberOfGaussians]) / gmm.numberOfGaussians)
        gmm.means[:self.gmm.numberOfGaussians - self.nTumorComp] = self.gmm.means[:-self.nTumorComp]
        gmm.variances[:self.gmm.numberOfGaussians - self.nTumorComp] = self.gmm.variances[:-self.nTumorComp]
        gmm.means[self.gmm.numberOfGaussians - self.nTumorComp:] = tumorGmm.means
        gmm.variances[self.gmm.numberOfGaussians - self.nTumorComp:] = tumorGmm.variances
        # not the mixture weights because they are just all ones
        """Create new extended arrays for the priors, likelihoods and posteriors, with 2 extra tumor-classes"""
        priorsNew = np.zeros(priors.shape + np.array([0, 2]))
        tumorIdx = np.zeros(priors.shape[1]).astype(np.bool)
        tumorIdx[tumorStructureNumber] = True
        priorsNew[:, :priors.shape[1] - 1] = atlasPriors[:, ~tumorIdx]
        priorsNew[:, priors.shape[1] - 1:] = atlasPriors[:, tumorIdx] * tumorGmm.mixtureWeights
        likelihoodsNew = np.zeros(priorsNew.shape)
        likelihoodsNew[:, :likelihoods.shape[1] - 1] = likelihoods[:, ~tumorIdx]
        likelihoodsNew[:, likelihoods.shape[1] - 1:] = np.array(
            [tumorGmm.getGaussianLikelihoods(data, np.expand_dims(tumorGmm.means[i], 1), tumorGmm.variances[i]) for i in
             range(self.nTumorComp + 1)]).T
        posteriors = likelihoodsNew * priorsNew
        posteriors /= np.expand_dims(np.sum(posteriors, axis=1) + eps, 1)
        likelihoods = likelihoodsNew
        priors = priorsNew
        classFractionsNew = np.zeros([gmm.numberOfClasses, likelihoods.shape[1]])
        tumorIdx = np.zeros(self.classFractions.shape[0]).astype(np.bool)
        tumorIdx[tumorClassNumber] = True
        classFractionsNew[:self.classFractions.shape[0] - 1, :self.classFractions.shape[1]] = self.classFractions[
            ~tumorIdx]
        classFractionsNew[self.classFractions.shape[0] - 1, :self.classFractions.shape[1]] = self.classFractions[
            tumorIdx]
        classFractionsNew[self.classFractions.shape[0]:, self.classFractions.shape[1]:] = np.eye(2)
        self.classFractions = classFractionsNew
        averagePosteriors = np.zeros(priors.shape)
        averageTumorPosteriors = np.zeros([data.shape[0], self.nTumorComp + 1])
        if not self.groundTruth is None:
            dice_score_hist,dice_score_hist_avg = [],[]
            scores = np.array(
                [dice_score(self.groundTruth[self.cropping][self.mask].flatten() == k,
                            np.argmax(tumor[self.mask], axis=-1).flatten() == k) for
                 k in range(4)])
            dice_score_hist.append(scores)
            dice_score_hist_avg.append(scores)
        configurations_hist = [] #TODO: Implement sampling of configurations
        currentConfiguration = np.array([1,2,3])
        for sweepNumber in range(self.numberOfBurnInSteps + self.numberOfSamplingSteps):
            tumorTrainSpace = np.zeros(np.concatenate([vaeShape, [self.nTumorComp + 2]]))
            for tumorCompNumber in range(self.nTumorComp + 2):
                tumorTrainSpace[:, :, :, tumorCompNumber] = affine_transform(tumor[:, :, :, tumorCompNumber],
                                                                             trainToSubjectMat, output_shape=vaeShape,
                                                                             order=1)
            if sweepNumber < self.numberOfBurnInSteps and False: #can use this to switch configurations, leave it out for now
                configLoss = []
                configs = configurations
                for config in configs:
                    tumorTmp = tumorTrainSpace[:, :, :, np.array([0] + list(config))]
                    mean, logvar = vae.encode(np.expand_dims(tumorTmp, 0))
                    z = vae.reparameterize(mean, logvar)
                    tumorVAE = vae.decode(z, apply_softmax=True).numpy()[0]
                    vaeLoss = log_loss(np.argmax(tumorTmp.reshape(-1, tumor.shape[-1]), axis=1),
                                       tumorVAE.reshape(-1, tumorVAE.shape[-1]))
                    configLoss.append(vaeLoss)
                currentConfiguration = np.array(configurations[np.argmin(configLoss)])
                if not (currentConfiguration==np.array([1,2,3])).all():
                    print("Changed configuration")
            else:
                currentConfiguration = np.array([1,2,3])

            mean, logvar = vae.encode(
                np.expand_dims(tumorTrainSpace[:, :, :, np.array([0] + list(currentConfiguration))], 0))
            z = vae.reparameterize(mean, logvar)
            tumorVAETrainSpace = vae.decode(z, apply_softmax=True).numpy()[0]
            tumorPriorSubjSpace = np.zeros(tumor.shape)
            for tumorCompNumber in range(self.nTumorComp + 2):
                tumorPriorSubjSpace[:, :, :, tumorCompNumber] = affine_transform(
                    tumorVAETrainSpace[:, :, :, tumorCompNumber], subjectToTrainMat, output_shape=self.mask.shape,
                    order=1)

            tumorPriorSubjSpace = tumorPriorSubjSpace[self.mask]
            newMeans, newVariances = [], []
            if sweepNumber > 0: #skip sampling parameters in the first iteration
                for j,compNumber in enumerate(np.arange(self.nTumorComp + 1)[currentConfiguration-1]):
                    constraints = None
                    if self.useMeanConstr:
                        if j == 1 and self.flair_idx is not None:
                            constraint = (self.flair_idx, [max(self.gmm.means[self.wm_comp_number, self.flair_idx],
                                                               self.gmm.means[
                                                                   self.gm_comp_number, self.flair_idx]) + np.log(1.15),
                                                           np.inf])
                            constraints = [constraint]
                        if j == 2:
                            constraints = []
                            if self.flair_idx is not None:
                                constraints.append(
                                    (self.flair_idx, [max(self.gmm.means[self.wm_comp_number, self.flair_idx],
                                                          self.gmm.means[self.gm_comp_number, self.flair_idx]),
                                                      np.inf]))
                            if self.t1ce_idx is not None:
                                constraints.append(
                                    (self.t1ce_idx, [max(self.gmm.means[self.wm_comp_number, self.t1ce_idx],
                                                         self.gmm.means[self.gm_comp_number, self.t1ce_idx]) + np.log(
                                        1.1), np.inf]))
                    newMean, newVariance = gmm.sampleMeansAndVariancesConditioned(data, np.expand_dims(
                        posteriors[:, posteriors.shape[1] - 3 + compNumber], 1), sum(self.gmm.numberOfGaussiansPerClass[:-1]) + compNumber,
                                                                                  constraints=constraints)
                    newMeans.append(newMean)
                    newVariances.append(newVariance)
                # newMixtureWeights = dirichlet.rvs(np.ones(self.nTumorComp)+posteriorsTumor.sum(axis=0))[0]
                # gmm.mixtureWeights[tumorClassNumber:tumorClassNumber+self.nTumorComp] = newMixtureWeights
                gmm.means[-(self.nTumorComp+1):] = np.array(newMeans).squeeze()
                gmm.variances[-(self.nTumorComp+1):] = np.array(newVariances)
            effectivePriors = priors.copy()
            effectivePriors[:, -(self.nTumorComp + 1):] *= tumorPriorSubjSpace[:, 1:]
            effectivePriors[:, :-(self.nTumorComp + 1)] *= tumorPriorSubjSpace[:, 0, np.newaxis]
            effectivePriors /= np.expand_dims(np.sum(effectivePriors, axis=1) + eps, 1)

            likelihoods = gmm.getLikelihoods(data, self.classFractions)
            posteriors = effectivePriors * likelihoods
            posteriors /= np.expand_dims(np.sum(posteriors, axis=1) + eps, 1)
            posteriorsTumor = posteriors[:, -(self.nTumorComp + 1):]

            tumorp = np.concatenate([np.expand_dims(1 - posteriorsTumor.sum(axis=1), 1), posteriorsTumor], axis=1)
            tumorp[tumorp<eps] = 0
            tumorp /= np.expand_dims(tumorp.sum(axis=1), 1)
            tumor[self.mask] = multinomial_rvs(1, tumorp)
            if self.groundTruth is not None: #TODO: Remove this after development
                scores = np.array(
                    [dice_score(self.groundTruth[self.cropping][self.mask].flatten() == k, np.argmax(tumorp, axis=-1).flatten() == k) for
                     k in range(4)])
                scores_avg = np.array(
                    [dice_score(self.groundTruth[self.cropping][self.mask].flatten() == k, np.argmax(averageTumorPosteriors, axis=-1).flatten() == k) for
                     k in range(4)])
                dice_score_hist.append(scores)
                dice_score_hist_avg.append(scores_avg)
                print("sample dice: ",scores)
                print("avg dice: ",scores_avg)
            # Collect data after burn in steps
            if sweepNumber >= self.numberOfBurnInSteps:
                print('Sample ' + str(sweepNumber + 1 - self.numberOfBurnInSteps) + ' times')
                averagePosteriors += posteriors / self.numberOfSamplingSteps
                averageTumorPosteriors += posteriorsTumor / self.numberOfSamplingSteps
            else:
                print('Burn-in ' + str(sweepNumber + 1) + ' times')
        if not self.groundTruth is None:
            scoreDf = pd.DataFrame()
            dice_score_hist = np.array(dice_score_hist)
            dice_score_hist_avg = np.array(dice_score_hist_avg)
            scoreDf["BG"] = dice_score_hist[:,0]
            scoreDf["NE"] = dice_score_hist[:,1]
            scoreDf["ED"] = dice_score_hist[:,2]
            scoreDf["EC"] = dice_score_hist[:,3]
            scoreDf["BGavg"] = dice_score_hist_avg[:,0]
            scoreDf["NEavg"] = dice_score_hist_avg[:,1]
            scoreDf["EDavg"] = dice_score_hist_avg[:,2]
            scoreDf["ECavg"] = dice_score_hist_avg[:,3]
            scoreDf.to_csv(os.path.join(self.savePath,"diceHistory.csv"))
        self.gmm = gmm
        return averagePosteriors, biasFields, nodePositions, data, priors, averageTumorPosteriors





def imagesInAtlasSpace(images,deformation,transform,referenceMesh):
    """
    images: array of n volumes with shape (x,y,z,n) (in subject space)
    deformation: the estimated deformation in reference space
    transform: the transform matrix from reference to subject space
    referenceMesh: original reference mesh
    """
    imageSize = list((referenceMesh.points.max(axis=0)+1).astype(np.int))
    denseDeformation = np.zeros(imageSize + [3], dtype=np.double)
    mesh = referenceMesh #just to be clear that it's no longer the "reference" mesh
    mesh.points = mesh.points+deformation
    if np.max(np.absolute(deformation)) > 0:
        maxDeformation = np.max(deformation)
        minDeformation = np.min(deformation)
        deltaDeformation = maxDeformation - minDeformation
        mesh.alphas = (deformation - minDeformation) / deltaDeformation
        tmp = mesh.rasterize_warp(imageSize, -1)
        # Unvisited voxels are marked by zeroes in all three coordinates - except possibly for the origin
        # which has all three coordinates zero as its natural state
        validMask = np.absolute(tmp)
        validMask = np.sum(validMask, axis=3)  ## only zero where all three coordinates are zero
        denseDeformation = (tmp) / (2 ** 16 - 1) * (maxDeformation - minDeformation) + minDeformation

        # Due to tetrahedral inside/outside checking performed in the rasterizer, some voxels on the boundary of the
        # image grid are never visited. There also seems to be non-visited voxels inside the image grid (bug in rasterizer).
        # We can resolve both issues by giving to each non-visited voxel the deformation field value of the closest voxel
        # that has been visited, augmented by knowledge that certain deformation components of voxels on the image grid
        # boundaries are zero by theory (sliding boundary conditions)
        closestIndices = scipy.ndimage.morphology.distance_transform_edt(
            validMask == 0, return_indices=True, return_distances=False)
        denseDeformation = denseDeformation[closestIndices[0], closestIndices[1], closestIndices[2]]

        # Boundaries are known a priori
        denseDeformation[0,:,:,0] = 0
        denseDeformation[-1,:,:,0] = 0
        denseDeformation[:,0,:,1] = 0
        denseDeformation[:,-1,:,1] = 0
        denseDeformation[:,:,0,2] = 0
        denseDeformation[:,:,-1,2] = 0

    linspaces = [np.linspace(0,i-1,i) for i in imageSize]
    X,Y,Z = np.meshgrid(linspaces[0],linspaces[1],linspaces[2],indexing='ij')
    densePositions = np.array([X,Y,Z]).transpose(1,2,3,0)
    densePositions = densePositions + denseDeformation
    densePositions = densePositions.reshape(-1,densePositions.shape[-1])
    atlas = samseg.ProbabilisticAtlas()
    densePositionsInSubjectSpace = atlas.mapPositionsFromTemplateToSubjectSpace(densePositions,transform)
    imagesOut = np.zeros(imageSize + [images.shape[-1]])
    for i in range(images.shape[-1]):
        denseIntensitiesInTemplateSpace = map_coordinates(images[:,:,:,i],densePositionsInSubjectSpace.T,order=1,mode='nearest')
        imagesOut[:,:,:,i] = denseIntensitiesInTemplateSpace.reshape(imageSize)
    return imagesOut


def setTumorComp(nTumorComp, atlasDir, modelSpecifications):

    # Create default model specifications as a dictionary
    FreeSurferLabels, names, colors = kvlReadCompressionLookupTable(os.path.join(atlasDir, 'compressionLookupTable.txt'))
    sharedGMMParameters = kvlReadSharedGMMParameters(os.path.join(atlasDir, 'sharedGMMParameters.txt'))
    newGMMParameter = GMMparameter("Tumor",nTumorComp,["Tumor"])
    TumorGMMNumber = None
    for classNumber, mergeOption in enumerate(sharedGMMParameters):
        if 'Tumor' == mergeOption.mergedName:
            TumorGMMNumber = classNumber
    sharedGMMParameters.pop(TumorGMMNumber)
    sharedGMMParameters.append(newGMMParameter)
    modelSpecifications.sharedGMMParameters = sharedGMMParameters
    return modelSpecifications

def makeTumorAtlas(atlasDirOrig,atlasDirNew,tumorPrior=0.2,skull=True):
    os.makedirs(atlasDirNew,exist_ok=True)
    levelFiles = [join(atlasDirOrig, i) for i in os.listdir(atlasDirOrig) if "level" in i.lower()]
    for meshFile in levelFiles:
        meshFileOut = meshFile.replace(atlasDirOrig, atlasDirNew)
        meshcoll = gems.KvlMeshCollection()
        meshcoll.read(meshFile)
        mesh = meshcoll.reference_mesh
        newAlphas = np.zeros((mesh.alphas.shape[0], mesh.alphas.shape[1] + 1))
        newAlphas[:, :-1] = mesh.alphas.copy()
        if skull:
            tumorPossible = mesh.alphas[:,:3].sum(axis=-1)<0.5
        else:
            tumorPossible = mesh.alphas[:, 0] < 0.5  # can only have tumor inside the brain
        newAlphas[tumorPossible, -1] = tumorPrior
        newAlphas[tumorPossible, :-1] *= (1-tumorPrior)
        sumAlphas = newAlphas[tumorPossible, :].sum(axis=1)
        newAlphas[tumorPossible, :] /= (sumAlphas[:, np.newaxis])
        mesh.alphas = newAlphas
        meshcoll.reference_mesh.alphas = mesh.alphas
        meshcoll.write(meshFileOut.replace(".gz", ""))

    sharedParamFile = join(atlasDirOrig, "sharedGMMParameters.txt")
    sharedParamLines = open(sharedParamFile).readlines()
    sharedParamLines.append("Tumor 1 Tumor\n")
    with open(join(atlasDirNew, "sharedGMMParameters.txt"), "w") as f:
        for line in sharedParamLines:
            f.write(line)
    LUTFile = join(atlasDirOrig, "compressionLookupTable.txt")
    LUTLines = open(LUTFile).readlines()
    LUTLines.append("99 %d Tumor                        255    0    0  255\n" % len(LUTLines))
    with open(join(atlasDirNew, "compressionLookupTable.txt"), "w") as f:
        for line in LUTLines:
            f.write(line)
    otherFiles = [join(atlasDirOrig, i) for i in os.listdir(atlasDirOrig) if
                  not "level" in i.lower() and not "parameters.txt" in i.lower() and not "compressionlookuptable.txt" in i.lower()]
    _ = [copyfile(i, i.replace(atlasDirOrig, atlasDirNew)) for i in otherFiles]


def getNeighborIndices(shape):
    """
    Given a 3D volume shape, compute for each voxel, the index of all neighboring voxels
    output is a 4D array: (shape[0],shape[1],shape[2],26)
    """
    ix = np.arange(np.prod(shape)).reshape(shape)
    diffs = []
    for j in range(27):
        diff = np.array([j//9,(j//3)%3,j%3])-1
        if (diff == 0).all():
            continue
        diffs.append(diff)
    diffs = np.array(diffs)
    XX, YY, ZZ = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), np.arange(shape[2]), indexing='ij')
    XYZ = np.array([XX,YY,ZZ]).transpose(1,2,3,0)
    index_arr = (XYZ[:,:,:,np.newaxis,:]-diffs[np.newaxis,np.newaxis,np.newaxis,:])
    index_arr = np.clip(index_arr,np.array([0,0,0]),np.array(shape)-1)
    neighbors_x = index_arr[:,:,:,:,0].reshape(-1,26)
    neighbors_y = index_arr[:,:,:,:,1].reshape(-1,26)
    neighbors_z = index_arr[:,:,:,:,2].reshape(-1,26)
    index_arr = ix[neighbors_x,neighbors_y,neighbors_z].reshape(list(shape)+[26])
    return index_arr

def multinomial_rvs(n, p):
    """
    Sample from the multinomial distribution with multiple p vectors.

    * n must be a scalar.
    * p must an n-dimensional numpy array, n >= 1.  The last axis of p
      holds the sequence of probabilities for a multinomial distribution.

    The return value has the same shape as p.
    """
    count = np.full(p.shape[:-1], n)
    out = np.zeros(p.shape, dtype=int)
    ps = p.cumsum(axis=-1)
    # Conditional probabilities
    with np.errstate(divide='ignore', invalid='ignore'):
        condp = p / ps
    condp[np.isnan(condp)] = 0.0
    for i in range(p.shape[-1]-1, 0, -1):
        binsample = np.random.binomial(count, condp[..., i])
        out[..., i] = binsample
        count -= binsample
    out[..., 0] = count
    return out

def dice_score(gt,seg):
    try:
        dice = np.sum(seg[gt==1])*2.0 / (np.sum(seg) + np.sum(gt))
    except:
        return np.nan
    return dice
  
  
class CatVAE(tf.keras.Model):
    def __init__(self, width, height, depth, numClasses, latent_dim, temperature):
        super(CatVAE, self).__init__()

        self.latent_dim = latent_dim
        self.width = width
        self.height = height
        self.depth = depth
        self.numClasses = numClasses
        self.temperature = temperature

        self.inference_net = tf.keras.models.Sequential(
            [
                tf.keras.layers.InputLayer(input_shape=(self.width, self.height, self.depth, self.numClasses)),

                # First layer
                tf.keras.layers.Conv3D(filters=32, kernel_size=[3, 3, 3], strides=[2, 2, 2],
                                       padding='SAME', activation=tf.nn.leaky_relu),
                tf.keras.layers.Conv3D(filters=32, kernel_size=[3, 3, 3], strides=[2, 2, 2],
                                       padding='SAME', activation=tf.nn.leaky_relu),
                tf.keras.layers.MaxPool3D(pool_size=[2, 2, 2], strides=[2, 2, 2], padding='SAME'),

                # Second layer
                tf.keras.layers.Conv3D(filters=32, kernel_size=[3, 3, 3], strides=[2, 2, 2],
                                       padding='SAME', activation=tf.nn.leaky_relu),
                tf.keras.layers.Conv3D(filters=32, kernel_size=[3, 3, 3], strides=[1, 1, 1],
                                       padding='SAME', activation=tf.nn.leaky_relu),
                tf.keras.layers.MaxPool3D(pool_size=[2, 2, 2], strides=[2, 2, 2], padding='SAME'),

                # Third layer
                tf.keras.layers.Conv3D(filters=32, kernel_size=[3, 3, 3], padding='SAME', activation=tf.nn.leaky_relu),
                tf.keras.layers.Conv3D(filters=32, kernel_size=[3, 3, 3], padding='SAME', activation=tf.nn.leaky_relu),
                tf.keras.layers.MaxPool3D(pool_size=[2, 2, 2], strides=[2, 2, 2], padding='SAME'),

                # Flatten and use fully connected for mean and variance
                tf.keras.layers.Flatten(),

                tf.keras.layers.Dense(units=256, activation=tf.nn.leaky_relu),

                tf.keras.layers.Dense(self.latent_dim + self.latent_dim)
            ]
        )

        self.generative_net = tf.keras.models.Sequential(
            [
                tf.keras.layers.InputLayer(input_shape=(self.latent_dim,)),

                tf.keras.layers.Dense(units=2 * 2 * 2 * 32, activation=tf.nn.leaky_relu),
                tf.keras.layers.Dense(units=4 * 4 * 3 * 32, activation=tf.nn.leaky_relu),
                tf.keras.layers.Reshape(target_shape=(4, 4, 3, 32)),

                # First deconv layer
                tf.keras.layers.Conv3DTranspose(filters=32, kernel_size=[3, 3, 3], strides=[1, 1, 1],
                                                padding='SAME', activation=tf.nn.leaky_relu),
                tf.keras.layers.Conv3DTranspose(filters=32, kernel_size=[3, 3, 3], strides=[1, 1, 1],
                                                padding='SAME', activation=tf.nn.leaky_relu),

                # Second deconv layer
                tf.keras.layers.Conv3DTranspose(filters=32, kernel_size=[3, 3, 3], strides=[2, 2, 2],
                                                padding='SAME', activation=tf.nn.leaky_relu),
                tf.keras.layers.Conv3DTranspose(filters=32, kernel_size=[3, 3, 3], strides=[2, 2, 2],
                                                padding='SAME', activation=tf.nn.leaky_relu),

                # Third deconv layer
                tf.keras.layers.Conv3DTranspose(filters=32, kernel_size=[3, 3, 3], strides=[2, 2, 2],
                                                padding='SAME', activation=tf.nn.leaky_relu),
                tf.keras.layers.Conv3DTranspose(filters=32, kernel_size=[3, 3, 3], strides=[2, 2, 2],
                                                padding='SAME', activation=tf.nn.leaky_relu),

                # Fourth deconv layer
                tf.keras.layers.Conv3DTranspose(filters=32, kernel_size=[3, 3, 3], strides=[2, 2, 2],
                                                padding='SAME', activation=tf.nn.leaky_relu),
                tf.keras.layers.Conv3DTranspose(filters=self.numClasses, kernel_size=[3, 3, 3],
                                                strides=[2, 2, 2], padding='SAME'),
            ]
        )

    @tf.function
    def sample(self, eps=None, samples=1, sample_decode=False, apply_softmax=False):
        if eps is None:
            eps = tf.random.normal(shape=(samples, self.latent_dim))
        return self.decode(eps, sample_decode=sample_decode, apply_softmax=apply_softmax)

    def encode(self, x):
        mean, logvar = tf.split(self.inference_net(x), num_or_size_splits=2, axis=1)
        return mean, logvar

    def reparameterize(self, mean, logvar):
        eps = tf.random.normal(shape=mean.shape)
        return eps * tf.exp(logvar * .5) + mean

    def decode(self, z, sample_decode=False, apply_softmax=False):
        logits = self.generative_net(z)
        # scale logits to desired shape
        logits = self.padToProperSize(logits)
        if sample_decode:
            sample = tfp.distributions.RelaxedOneHotCategorical(temperature=self.temperature, logits=logits).sample()
            return sample
        if apply_softmax:
            return tf.nn.softmax(logits)
        return logits

    def printSummary(self):
        self.inference_net.summary()
        self.generative_net.summary()

    def padToProperSize(self, logits):
        #paddings = tf.constant([[0, 0], [0, 0], [0, 0], [5, 5], [0, 0]])
        #return tf.pad(logits, paddings, "CONSTANT")
        shape = logits.shape
        return logits[:, 8:shape[1]-8, 8:shape[2]-8, 19:shape[3]-18, :]
  
