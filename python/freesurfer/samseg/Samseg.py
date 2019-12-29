import os
import numpy as np
import logging
import pickle
import scipy.io
import freesurfer as fs

from .figures import initVisualizer
from .utilities import requireNumpyArray, Specification
from .BiasField import BiasField
from .ProbabilisticAtlas import ProbabilisticAtlas
from .GMM import GMM
from .Affine import Affine
from .SamsegUtility import *
from .merge_alphas import kvlMergeAlphas, kvlGetMergingFractionsTable
from . import gems

eps = np.finfo(float).eps


class Samseg:
    def __init__(self, imageFileNames, atlasDir, savePath, userModelSpecifications=None, userOptimizationOptions=None,
                 transformedTemplateFileName=None, visualizer=None, saveHistory=None, savePosteriors=None,
                 saveWarp=None, saveMesh=None, threshold=None, thresholdSearchString=None,
                 targetIntensity=None, targetSearchStrings=None):

        # Store input parameters as class variables
        self.imageFileNames = imageFileNames
        self.savePath = savePath
        self.atlasDir = atlasDir
        self.threshold = threshold
        self.thresholdSearchString = thresholdSearchString
        self.targetIntensity = targetIntensity
        self.targetSearchStrings = targetSearchStrings

        # Initialize some objects
        self.affine = Affine()
        self.probabilisticAtlas = ProbabilisticAtlas()

        # Get full model specifications and optimization options (using default unless overridden by user)
        self.modelSpecifications = getModelSpecifications(atlasDir, userModelSpecifications)
        self.optimizationOptions = getOptimizationOptions(atlasDir, userOptimizationOptions)

        # Get transformed template, if any
        self.transformedTemplateFileName = transformedTemplateFileName

        # Print specifications
        print('##----------------------------------------------')
        print('              Samsegment Options')
        print('##----------------------------------------------')
        print('output directory:', savePath)
        print('input images:', ', '.join([imageFileName for imageFileName in imageFileNames]))
        if self.transformedTemplateFileName is not None:
            print('transformed template:', self.transformedTemplateFileName)
        print('modelSpecifications:', self.modelSpecifications)
        print('optimizationOptions:', self.optimizationOptions)

        # Convert modelSpecifications from dictionary into something more convenient to access
        self.modelSpecifications = Specification(self.modelSpecifications)

        # Setup a null visualizer if necessary
        if visualizer is None:
            self.visualizer = initVisualizer(False, False)
        else:
            self.visualizer = visualizer

        self.saveHistory = saveHistory
        self.savePosteriors = savePosteriors
        self.saveWarp = saveWarp
        self.saveMesh = saveMesh

        # Make sure we can write in the target/results directory
        os.makedirs(savePath, exist_ok=True)

        # Class variables that will be used later
        self.biasField = None
        self.gmm = None
        self.imageBuffers = None
        self.mask = None
        self.classFractions = None
        self.cropping = None
        self.transform = None
        self.voxelSpacing = None
        self.optimizationSummary = None
        self.optimizationHistory = None
        self.deformation = None
        self.deformationAtlasFileName = None

    def register(self):
        # =======================================================================================
        #
        # Perform affine registration if needed
        #
        # =======================================================================================
        if self.transformedTemplateFileName is None:
            templateFileName = os.path.join(self.atlasDir, 'template.nii')
            affineRegistrationMeshCollectionFileName = os.path.join(self.atlasDir, 'atlasForAffineRegistration.txt.gz')
            _, self.transformedTemplateFileName, _ = self.affine.registerAtlas(self.imageFileNames[0],
                                                                               affineRegistrationMeshCollectionFileName,
                                                                               templateFileName,
                                                                               self.savePath,
                                                                               visualizer=self.visualizer)

    def preProcess(self):
        # =======================================================================================
        #
        # Preprocessing (reading and masking of data)
        #
        # =======================================================================================

        # Read the image data from disk. At the same time, construct a 3-D affine transformation (i.e.,
        # translation, rotation, scaling, and skewing) as well - this transformation will later be used
        # to initially transform the location of the atlas mesh's nodes into the coordinate system of the image.
        self.imageBuffers, self.transform, self.voxelSpacing, self.cropping = readCroppedImages(
            self.imageFileNames,
            self.transformedTemplateFileName)

        # Background masking: simply setting intensity values outside of a very rough brain mask to zero
        # ensures that they'll be skipped in all subsequent computations
        self.imageBuffers, self.mask = maskOutBackground(self.imageBuffers,
                                                                      self.modelSpecifications.atlasFileName,
                                                                      self.transform,
                                                                      self.modelSpecifications.brainMaskingSmoothingSigma,
                                                                      self.modelSpecifications.brainMaskingThreshold,
                                                                      self.probabilisticAtlas)

        # Let's prepare for the bias field correction that is part of the imaging model. It assumes
        # an additive effect, whereas the MR physics indicate it's a multiplicative one - so we log
        # transform the data first.
        self.imageBuffers = logTransform(self.imageBuffers, self.mask)

        # Visualize some stuff
        if hasattr(self.visualizer, 'show_flag'):
            self.visualizer.show(
                mesh=self.probabilisticAtlas.getMesh(self.modelSpecifications.atlasFileName, self.transform),
                shape=self.imageBuffers.shape,
                window_id='samsegment mesh', title='Mesh',
                names=self.modelSpecifications.names, legend_width=350)
            self.visualizer.show(images=self.imageBuffers, window_id='samsegment images',
                                 title='Samsegment Masked and Log-Transformed Contrasts')

    def process(self):
        # =======================================================================================
        #
        # Parameter estimation
        #
        # =======================================================================================
        self.initializeBiasField()
        self.initializeGMM()
        self.estimateModelParameters()

    def postProcess(self):
        # =======================================================================================
        #
        # Segment the data using the estimate model parameters, and write results out
        #
        # =======================================================================================

        # OK, now that all the parameters have been estimated, try to segment the original, full resolution image
        # with all the original labels (instead of the reduced "super"-structure labels we created)
        posteriors, biasFields, nodePositions, _, _ = self.segment()

        # Write out segmentation and bias field corrected volumes
        volumesInCubicMm = writeResults(self.imageFileNames, self.savePath, self.imageBuffers, self.mask,
                                                     biasFields,
                                                     posteriors, self.modelSpecifications.FreeSurferLabels,
                                                     self.cropping,
                                                     self.targetIntensity, self.targetSearchStrings,
                                                     self.modelSpecifications.names,
                                                     self.threshold, self.thresholdSearchString,
                                                     savePosteriors=self.savePosteriors)

        # Save the template warp
        if self.saveWarp:
            print('Saving the template warp')
            self.saveWarpField(os.path.join(self.savePath, 'template.m3z'))

        # Save the final mesh collection
        if self.saveMesh:
            print('Saving the final mesh in template space')
            image_base_path, _ = os.path.splitext(self.imageFileNames[0])
            _, scanName = os.path.split(image_base_path)
            deformedAtlasFileName = os.path.join(self.savePath, scanName + '_meshCollection.txt')
            self.probabilisticAtlas.saveDeformedAtlas(self.modelSpecifications.atlasFileName, deformedAtlasFileName,
                                                      nodePositions)

        # Save the history of the parameter estimation process
        if self.saveHistory:
            history = {'input': {
                'imageFileNames': self.imageFileNames,
                'transformedTemplateFileName': self.transformedTemplateFileName,
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

    def saveWarpField(self, filename):
        # extract node positions in image space
        nodePositions = self.probabilisticAtlas.getMesh(
            self.modelSpecifications.atlasFileName,
            self.transform,
            initialDeformation=self.deformation,
            initialDeformationMeshCollectionFileName=self.deformationAtlasFileName
        ).points

        # extract geometries
        imageGeom = fs.Volume.read(self.imageFileNames[0]).geometry()
        templateGeom = fs.Volume.read(os.path.join(self.atlasDir, 'template.nii')).geometry()

        # extract image-to-image template transform
        matricesFileName = os.path.join(self.savePath, 'template_coregistrationMatrices.mat')
        matrix = scipy.io.loadmat(matricesFileName)['imageToImageTransformMatrix']

        # rasterize the final node coordinates (in image space) using the initial template mesh
        mesh = self.probabilisticAtlas.getMesh(self.modelSpecifications.atlasFileName)
        coordmap = mesh.rasterize_values(templateGeom.shape, nodePositions).astype('float')

        # adjust for the offset introduced by volume cropping
        coordmap[..., :] += [slc.start for slc in self.cropping]

        # the rasterization is a bit buggy and some voxels are not filled - mark these as invalid
        invalid = np.any(coordmap == 0, axis=-1)
        coordmap[invalid, :] = -1

        # write the warp file
        fs.Warp(coordmap, source=imageGeom, target=templateGeom, affine=matrix).write(filename)

    def getDownSampledModel(self, atlasFileName, downSamplingFactors):

        # Downsample the images and basis functions
        numberOfContrasts = self.imageBuffers.shape[-1]
        downSampledMask = self.mask[::downSamplingFactors[0], ::downSamplingFactors[1], ::downSamplingFactors[2]]
        downSampledImageBuffers = np.zeros(downSampledMask.shape + (numberOfContrasts,), order='F')
        for contrastNumber in range(numberOfContrasts):
            # logger.debug('first time contrastNumber=%d', contrastNumber)
            downSampledImageBuffers[:, :, :, contrastNumber] = self.imageBuffers[::downSamplingFactors[0],
                                                               ::downSamplingFactors[1],
                                                               ::downSamplingFactors[2],
                                                               contrastNumber]

        # Compute the resulting transform, taking into account the downsampling
        downSamplingTransformMatrix = np.diag(1. / downSamplingFactors)
        downSamplingTransformMatrix = np.pad(downSamplingTransformMatrix, (0, 1), mode='constant', constant_values=0)
        downSamplingTransformMatrix[3][3] = 1
        downSampledTransform = gems.KvlTransform(
            requireNumpyArray(downSamplingTransformMatrix @ self.transform.as_numpy_array))

        # Get the mesh
        downSampledMesh, downSampledInitialDeformationApplied = self.probabilisticAtlas.getMesh(atlasFileName,
                                                                                                downSampledTransform,
                                                                                                self.modelSpecifications.K,
                                                                                                self.deformation,
                                                                                                self.deformationAtlasFileName,
                                                                                                returnInitialDeformationApplied=True)

        return downSampledImageBuffers, downSampledMask, downSampledMesh, downSampledInitialDeformationApplied, \
               downSampledTransform,

    def initializeBiasField(self):

        # Our bias model is a linear combination of a set of basis functions. We are using so-called "DCT-II" basis functions,
        # i.e., the lowest few frequency components of the Discrete Cosine Transform.
        self.biasField = BiasField(self.imageBuffers.shape[0:3],
                                   self.modelSpecifications.biasFieldSmoothingKernelSize / self.voxelSpacing)

        # Visualize some stuff
        if hasattr(self.visualizer, 'show_flag'):
            import matplotlib.pyplot as plt  # avoid importing matplotlib by default
            plt.ion()
            f = plt.figure('Bias field basis functions')
            for dimensionNumber in range(3):
                plt.subplot(2, 2, dimensionNumber + 1)
                plt.plot(self.biasField.basisFunctions[dimensionNumber])
            plt.draw()

    def initializeGMM(self):
        # The fact that we consider neuro-anatomical structures as mixtures of "super"-structures for the purpose of model
        # parameter estimation, but at the same time represent each of these super-structures with a mixture of Gaussians,
        # creates something of a messy situation when implementing this stuff. To avoid confusion, let's define a few
        # conventions that we'll closely follow in the code as follows:
        #
        #   - classNumber = 1 ... numberOfClasses  -> indexes a specific super-structure (there are numberOfClasses superstructures)
        #   - numberOfGaussiansPerClass            -> a numberOfClasses-dimensional vector that indicates the number of components
        #                                             in the Gaussian mixture model associated with each class
        #   - gaussianNumber = 1 .... numberOfGaussians  -> indexes a specific Gaussian distribution; there are
        #                                                   numberOfGaussians = sum( numberOfGaussiansPerClass ) of those in total
        #   - classFractions -> a numberOfClasses x numberOfStructures table indicating in each column the mixing weights of the
        #                       various classes in the corresponding structure
        numberOfGaussiansPerClass = [param.numberOfComponents for param in self.modelSpecifications.sharedGMMParameters]
        self.classFractions, _ = kvlGetMergingFractionsTable(self.modelSpecifications.names,
                                                             self.modelSpecifications.sharedGMMParameters)

        # Parameter initialization.
        self.gmm = GMM(numberOfGaussiansPerClass, numberOfContrasts=self.imageBuffers.shape[-1],
                       useDiagonalCovarianceMatrices=self.modelSpecifications.useDiagonalCovarianceMatrices)

    def estimateModelParameters(self, initialBiasFieldCoefficients=None, initialDeformation=None,
                                initialDeformationAtlasFileName=None,
                                skipGMMParameterEstimationInFirstIteration=False,
                                skipBiasFieldParameterEstimationInFirstIteration=True):

        #
        logger = logging.getLogger(__name__)
        self.optimizationHistory = []
        self.optimizationSummary = []

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
            self.visualizer.show(mesh=downSampledMesh, images=downSampledImageBuffers,
                                 window_id='Mesh deformation (level ' + str(multiResolutionLevel) + ')',
                                 title='Mesh Deformation (level ' + str(multiResolutionLevel) + ')')

            if self.saveHistory:
                levelHistory = {'historyWithinEachIteration': []}

            # Main iteration loop over both EM and deformation
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
                    downSampledGaussianPosteriors, minLogLikelihood = self.gmm.getGaussianPosteriors(downSampledData,
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
                    if (EMIterationNumber == 100) or (changeCostEMPerVoxel < changeCostEMPerVoxelThreshold):
                        # Converged
                        print('EM converged!')
                        break

                    # M-step: update the model parameters based on the current posterior
                    #
                    # First the mixture model parameters
                    if not ((iterationNumber == 0) and skipGMMParameterEstimationInFirstIteration):
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
                    tmp = np.zeros(downSampledMask.shape + (downSampledGaussianPosteriors.shape[-1], ))
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

        # End resolution level loop

    def segment(self):
        # Get the final mesh
        mesh = self.probabilisticAtlas.getMesh(self.modelSpecifications.atlasFileName, self.transform,
                                               initialDeformation=self.deformation,
                                               initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)

        # Get the priors as dictated by the current mesh position
        priors = mesh.rasterize(self.imageBuffers.shape[0:3], -1)
        priors = priors[self.mask, :]

        # Get bias field corrected data
        # Make sure that the bias field basis function are not downsampled
        # (this might happens if the parameters estimation is made only with one downsampled resolution)
        self.biasField.downSampleBasisFunctions([1, 1, 1])
        biasFields = self.biasField.getBiasFields()
        data = self.imageBuffers[self.mask, :] - biasFields[self.mask, :]

        # Compute the posterior distribution of the various structures
        posteriors = self.gmm.getPosteriors(data, priors, self.classFractions)

        #
        estimatedNodePositions = self.probabilisticAtlas.mapPositionsFromSubjectToTemplateSpace(mesh.points,
                                                                                                self.transform)

        #
        return posteriors, biasFields, estimatedNodePositions, data, priors
