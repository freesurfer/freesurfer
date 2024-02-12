import os
import numpy as np
import pickle

import gems
from gems.SamsegUtility import *
from gems.utilities import requireNumpyArray
from gems.figures import initVisualizer
from gems.Affine import Affine
from gems.ProbabilisticAtlas import ProbabilisticAtlas
from gems.Samseg import Samseg


"""
Longitudinal version of samsegment
The idea is based on the generative model in the paper

  Iglesias, Juan Eugenio, et al.
  Bayesian longitudinal segmentation of hippocampal substructures in brain MRI using subject-specific atlases.
  Neuroimage 141 (2016): 542-555,
in which a subject-specific atlas is obtained by generating a random warp from the usual population atlas, and
subsequently each time point is again randomly warped from this subject-specific atlas. The intermediate
subject-specific atlas is effectively a latent variable in the model, and it's function is to encourage the
various time points to have atlas warps that are similar between themselves, without having to define a priori
what these warps should be similar to. In the implementation provided here, the stiffness of the first warp
(denoted by K0, as in the paper) is taken as the stiffness used in ordinary samseg, and the stiffness of the
second warp (denoted by K1) is controlled by the setting

  strengthOfLatentDeformationHyperprior

so that K1 = strengthOfLatentDeformationHyperprior * K0. In the Iglesias paper, the setting

  strengthOfLatentDeformationHyperprior = 1.0

was used.

The overall idea is extended here by adding also latent variables encouraging corresponding Gaussian mixture
models across time points to be similar across time -- again without having to define a priori what exactly
they should look like. For given values of these latent variables, they effectively act has hyperparameters
on the mixture model parameters, the strength of which is controlled through the setting

  strengthOfLatentGMMHyperprior

which weighs the relative strength of this hyperprior relative to the relevant (expected) data term in the
estimation procedures of the mixture model parameters at each time point. This aspect can be switched off by
setting

  strengthOfLatentGMMHyperprior = 0.0

NOTE: The general longitudinal pipeline in FreeSurfer 6.0 is described in the paper

  Reuter, Martin, et al.
  Within-subject template estimation for unbiased longitudinal image analysis.
  Neuroimage 61.4 (2012): 1402-1418,

which is based on the idea of retaining some temporal consistency by simply initializing the model fitting
procedure across all time points in exactly the same way. This is achieved by first creating a subject-specific
template that is subsequently analyzed and the result of which is then used as a (very good) initialization
in each time point. This behavior can be mimicked by setting

  initializeLatentDeformationToZero = True
  numberOfIterations = 1
  strengthOfLatentGMMHyperprior = 0.0
  strengthOfLatentDeformationHyperprior = 1.0

in the implementation provided here.
"""

eps = np.finfo( float ).eps


class SamsegLongitudinal:
    def __init__(self,
        imageFileNamesList,
        atlasDir,
        savePath,
        userModelSpecifications={},
        userOptimizationOptions={},
        visualizer=None, 
        saveHistory=False,
        saveMesh=None,
        targetIntensity=None,
        targetSearchStrings=None,
        numberOfIterations=5,
        strengthOfLatentGMMHyperprior=0.5,
        strengthOfLatentDeformationHyperprior=20.0,
        saveSSTResults=True,
        updateLatentMeans=True,
        updateLatentVariances=True,
        updateLatentMixtureWeights=True,
        updateLatentDeformation=True,
        initializeLatentDeformationToZero=False,
        threshold=None,
        thresholdSearchString=None,
        modeNames=None,
        pallidumAsWM=True,
        savePosteriors=False,
        saveModelProbabilities=False,
        tpToBaseTransforms=None,
        ):

        # Store input parameters as class variables
        self.imageFileNamesList = imageFileNamesList
        self.numberOfTimepoints = len(self.imageFileNamesList)
        self.savePath = savePath
        self.atlasDir = atlasDir
        self.threshold = threshold
        self.thresholdSearchString = thresholdSearchString
        self.targetIntensity = targetIntensity
        self.targetSearchStrings = targetSearchStrings
        self.numberOfIterations = numberOfIterations
        self.strengthOfLatentGMMHyperprior = strengthOfLatentGMMHyperprior
        self.strengthOfLatentDeformationHyperprior = strengthOfLatentDeformationHyperprior
        self.modeNames = modeNames
        self.pallidumAsWM = pallidumAsWM
        self.savePosteriors = savePosteriors
        self.tpToBaseTransforms = tpToBaseTransforms
        self.saveModelProbabilities = saveModelProbabilities

        # Check if all time point to base transforms are identity matrices.
        # If so, we can derive a combined 4D mask during preprocessing
        self.allIdentityTransforms = True
        if tpToBaseTransforms is not None:
            for tp, transform in enumerate(self.tpToBaseTransforms):
                if not np.allclose(transform.matrix, np.eye(4)):
                    self.allIdentityTransforms = False


        # Initialize some objects
        self.probabilisticAtlas = ProbabilisticAtlas()

        # Get full model specifications and optimization options (using default unless overridden by user)
        self.userModelSpecifications = userModelSpecifications
        self.userOptimizationOptions = userOptimizationOptions

        # Setup a null visualizer if necessary
        if visualizer is None:
            self.visualizer = initVisualizer(False, False)
        else:
            self.visualizer = visualizer

        self.saveHistory = saveHistory
        self.saveMesh = saveMesh
        self.saveSSTResults = saveSSTResults
        self.updateLatentMeans = updateLatentMeans
        self.updateLatentVariances = updateLatentVariances
        self.updateLatentMixtureWeights = updateLatentMixtureWeights
        self.updateLatentDeformation = updateLatentDeformation
        self.initializeLatentDeformationToZero = initializeLatentDeformationToZero

        # Make sure we can write in the target/results directory
        os.makedirs(savePath, exist_ok=True)

        # Here some class variables that will be defined later
        self.sstModel = None
        self.timepointModels = None
        self.imageBuffersList = None
        self.sstFileNames = None
        self.combinedImageBuffers = None
        self.latentDeformation = None
        self.latentDeformationAtlasFileName = None
        self.latentMeans = None
        self.latentVariances = None
        self.latentMixtureWeights = None
        self.latentMeansNumberOfMeasurements = None
        self.latentVariancesNumberOfMeasurements = None
        self.latentMixtureWeightsNumberOfMeasurements = None
        self.timepointVolumesInCubicMm = []
        self.optimizationSummary = None
        self.history = None
        self.historyOfTotalCost = None
        self.historyOfTotalTimepointCost = None
        self.historyOfLatentAtlasCost = None

    def segment(self, saveWarp=False,initTransformFile=None):
        # =======================================================================================
        #
        # Main function that runs the whole longitudinal segmentation pipeline
        #
        # =======================================================================================
        self.constructAndRegisterSubjectSpecificTemplate(initTransformFile)
        self.preProcess()
        self.fitModel()
        return self.postProcess(saveWarp=saveWarp)

    def constructAndRegisterSubjectSpecificTemplate(self,initTransformFile=None):
        # =======================================================================================
        #
        # Construction and affine registration of subject-specific template (sst)
        #
        # =======================================================================================

        # Initialization transform for registration
        initTransform = None
        if initTransformFile:
            trg = self.validateTransform(sf.load_affine(initTransformFile))
            initTransform = convertRASTransformToLPS(trg.convert(space='world').matrix)

        # Generate the subject specific template (sst)
        self.sstFileNames = self.generateSubjectSpecificTemplate()
        sstDir, _ = os.path.split(self.sstFileNames[0])

        # Affine atlas registration to sst
        templateFileName = os.path.join(self.atlasDir, 'template.nii')
        affineRegistrationMeshCollectionFileName = os.path.join(self.atlasDir, 'atlasForAffineRegistration.txt.gz')

        affine = Affine(imageFileName=self.sstFileNames[0],
                         meshCollectionFileName=affineRegistrationMeshCollectionFileName,
                         templateFileName=templateFileName)
        self.imageToImageTransformMatrix, _ = affine.registerAtlas(savePath=sstDir, visualizer=self.visualizer,initTransform=initTransform)


    def preProcess(self):

        # construct sstModel
        self.constructSstModel()

        # =======================================================================================
        #
        # Preprocessing (reading and masking of data)
        #
        # =======================================================================================

        templateFileName = os.path.join(self.atlasDir, 'template.nii')
        self.sstModel.imageBuffers, self.sstModel.transform, self.sstModel.voxelSpacing, self.sstModel.cropping = readCroppedImages(self.sstFileNames, templateFileName, self.imageToImageTransformMatrix)

        self.imageBuffersList = []
        self.voxelSpacings = []
        self.transforms = []
        self.masks = []
        self.croppings = []
        if self.allIdentityTransforms:

            self.imageBuffersList = []
            for imageFileNames in self.imageFileNamesList:
                imageBuffers, _, _, _ = readCroppedImages(imageFileNames, templateFileName,
                                                          self.imageToImageTransformMatrix)
                self.imageBuffersList.append(imageBuffers)

            # Put everything in a big 4-D matrix to derive one consistent mask across all time points
            imageSize = self.sstModel.imageBuffers.shape[:3]
            numberOfContrasts = self.sstModel.imageBuffers.shape[-1]
            self.combinedImageBuffers = np.zeros(imageSize + (numberOfContrasts * (1 + self.numberOfTimepoints),))
            self.combinedImageBuffers[..., 0:numberOfContrasts] = self.sstModel.imageBuffers
            for timepointNumber in range(self.numberOfTimepoints):
                self.combinedImageBuffers[..., (timepointNumber + 1) * numberOfContrasts:
                                               (timepointNumber + 2) * numberOfContrasts] = self.imageBuffersList[
                    timepointNumber]

            self.combinedImageBuffers, self.sstModel.mask = maskOutBackground(self.combinedImageBuffers,
                                                                              self.sstModel.modelSpecifications.atlasFileName,
                                                                              self.sstModel.transform,
                                                                              self.sstModel.modelSpecifications.maskingProbabilityThreshold,
                                                                              self.sstModel.modelSpecifications.maskingDistance,
                                                                              self.probabilisticAtlas,
                                                                              self.sstModel.voxelSpacing)
            combinedImageBuffers = logTransform(self.combinedImageBuffers, self.sstModel.mask)

            # Retrieve the masked sst and time points
            self.sstModel.imageBuffers = combinedImageBuffers[..., 0:numberOfContrasts]
            for timepointNumber in range(self.numberOfTimepoints):
                self.imageBuffersList[timepointNumber] = combinedImageBuffers[...,
                                                         (timepointNumber + 1) * numberOfContrasts:
                                                         (timepointNumber + 2) * numberOfContrasts]
                self.masks.append(self.sstModel.mask)
                self.croppings.append(self.sstModel.cropping)
                self.voxelSpacings.append(self.sstModel.voxelSpacing)
                self.transforms.append(self.sstModel.transform)

        else:

            for timepointNumber, imageFileNames in enumerate(self.imageFileNamesList):

                # Compute transformation from population atlas (p) to time point (tp), passing through the template space (s)
                # The transformation needs to be in vox to vox space as self.imageToImageTransformMatrix
                # We need to concatenate the following transformations
                # self.imageToImageTransformMatrix -> population to template space - vox to vox transform
                # tmp_s.geom.vox2world -> template space - vox to world transform
                # self.tpToBaseTransforms[timepointNumber].inv() -> template to time point space - world to world transform
                # tmp_tp.geom.world2vox -> time point space - world to vox transform
                tmpTp = sf.load_volume(imageFileNames[0])
                tmpS = sf.load_volume(os.path.join(self.savePath, "base", "template_coregistered.mgz"))
                pToTpTransform = tmpTp.geom.world2vox @ self.tpToBaseTransforms[timepointNumber].inv() @ tmpS.geom.vox2world @ self.imageToImageTransformMatrix

                imageBuffers, transform, voxelSpacing, cropping = readCroppedImages(imageFileNames, templateFileName, pToTpTransform.matrix)

                #
                self.imageBuffersList.append(imageBuffers)
                self.voxelSpacings.append(voxelSpacing)
                self.transforms.append(transform)
                self.croppings.append(cropping)

            # Derive mask for sst model
            imageBuffer, self.sstModel.mask = maskOutBackground(self.sstModel.imageBuffers, self.sstModel.modelSpecifications.atlasFileName,
                                                                self.sstModel.transform,
                                                                self.sstModel.modelSpecifications.maskingProbabilityThreshold,
                                                                self.sstModel.modelSpecifications.maskingDistance,
                                                                self.probabilisticAtlas,
                                                                self.sstModel.voxelSpacing)
            self.sstModel.imageBuffers = logTransform(imageBuffer, self.sstModel.mask)

            # Derive one mask for each time point model
            for timepointNumber in range(self.numberOfTimepoints):
                imageBuffer, timepointMask = maskOutBackground(self.imageBuffersList[timepointNumber],
                                                               self.sstModel.modelSpecifications.atlasFileName,
                                                               self.transforms[timepointNumber],
                                                               self.sstModel.modelSpecifications.maskingProbabilityThreshold,
                                                               self.sstModel.modelSpecifications.maskingDistance,
                                                               self.probabilisticAtlas,
                                                               self.voxelSpacings[timepointNumber])
                imageBuffer = logTransform(imageBuffer, timepointMask)
                self.imageBuffersList[timepointNumber] = imageBuffer
                self.masks.append(timepointMask)

        # construct timepoint models
        self.constructTimepointModels()

        self.visualizer.show(images=self.sstModel.imageBuffers, title='sst')
        for timepointNumber in range(self.numberOfTimepoints):
            self.visualizer.show(images=self.imageBuffersList[timepointNumber], title='time point ' + str(timepointNumber))

    def fitModel(self):

        # =======================================================================================
        #
        # Parameter estimation for SST
        #
        # =======================================================================================

        self.sstModel.fitModel()

        if hasattr(self.visualizer, 'show_flag'):
            import matplotlib.pyplot as plt  # avoid importing matplotlib by default
            plt.ion()
            self.sstModel.biasField.downSampleBasisFunctions([1, 1, 1])
            sstBiasFields = self.sstModel.biasField.getBiasFields( self.sstModel.mask)
            sstData = self.sstModel.imageBuffers[self.sstModel.mask, :] - sstBiasFields[self.sstModel.mask, :]
            axsList = []
            for contrastNumber in range(self.sstModel.gmm.numberOfContrasts):
                f = plt.figure()
                numberOfAxes = 2 + self.numberOfTimepoints
                numberOfRows = np.int(np.ceil(np.sqrt(numberOfAxes)))
                numberOfColumns = np.int(np.ceil(numberOfAxes / numberOfRows))
                axs = f.subplots(numberOfRows, numberOfColumns, sharex=True)
                ax = axs.ravel()[0]
                _, bins, _ = ax.hist(self.sstModel.imageBuffers[self.sstModel.mask, contrastNumber], 100)
                ax.grid()
                ax.set_title('sst before bias field correction')
                ax = axs.ravel()[1]
                ax.hist(sstData[:, contrastNumber], bins)
                ax.grid()
                ax.set_title('sst after bias field correction')
                for timepointNumber in range(self.numberOfTimepoints):
                    ax = axs.ravel()[2 + timepointNumber]
                    ax.hist(self.imageBuffersList[timepointNumber][self.masks[timepointNumber], contrastNumber], bins)
                    ax.grid()
                    ax.set_title('time point ' + str(timepointNumber))
                axsList.append(axs)
            plt.draw()

        if self.saveHistory:
            self.history = {
                "sstMeans": self.sstModel.gmm.means,
                "sstVariances": self.sstModel.gmm.variances,
                "sstMixtureWeights": self.sstModel.gmm.mixtureWeights,
                "sstBiasFieldCoefficients": self.sstModel.biasField.coefficients,
                "sstDeformation": self.sstModel.deformation,
                "sstDeformationAtlasFileName": self.sstModel.deformationAtlasFileName,
                "sstOptimizationSummary": self.sstModel.optimizationSummary,
                "sstOptimizationHistory": self.sstModel.optimizationHistory
            }

        if self.saveSSTResults:
            sstlabels, sstnames, sstVolumesInCubicMm, sstoptimizationSummary = self.sstModel.postProcess()
            if self.saveHistory:
                self.history["sstVolumesInCubicMm"] = sstVolumesInCubicMm

        # =======================================================================================
        #
        # Iterative parameter vs. latent variables estimation, using SST result for initialization
        # and/or anchoring of hyperprior strength
        #
        # =======================================================================================

        # Initialization of the time-specific model parameters
        for timepointNumber in range(self.numberOfTimepoints):
            self.timepointModels[timepointNumber].initializeGMM()
            self.timepointModels[timepointNumber].gmm.means = self.sstModel.gmm.means.copy()
            self.timepointModels[timepointNumber].gmm.variances = self.sstModel.gmm.variances.copy()
            self.timepointModels[timepointNumber].gmm.mixtureWeights = self.sstModel.gmm.mixtureWeights.copy()
            self.timepointModels[timepointNumber].initializeBiasField()

        # Initialization of the latent variables, acting as hyperparameters when viewed from the model parameters' perspective
        self.latentDeformation = self.sstModel.deformation.copy()
        self.latentDeformationAtlasFileName = self.sstModel.deformationAtlasFileName
        self.latentMeans = self.sstModel.gmm.means.copy()
        self.latentVariances = self.sstModel.gmm.variances.copy()
        self.latentMixtureWeights = self.sstModel.gmm.mixtureWeights.copy()

        if self.initializeLatentDeformationToZero:
            for timepointNumber in range(self.numberOfTimepoints):
                self.timepointModels[timepointNumber].deformation = self.latentDeformation.copy()
                self.timepointModels[timepointNumber].latentDeformationAtlasFileName = self.latentDeformationAtlasFileName
                self.latentDeformation[:] = 0

        # Strength of the hyperprior (i.e., how much the latent variables control the conditional posterior of the parameters)
        # is user-controlled.
        #
        # For the GMM part, I'm using the *average* number of voxels assigned to the components in each mixture (class) of the
        # SST segmentation, so that all the components in each mixture are well-regularized (and tiny components don't get to do
        # whatever they want)
        #
        # Note that we need to take into account the possible resolution difference between SST and each time point.
        # Here we assume that these time points have similar resolution, otherwise the mean might not be the best choice
        # Scale latent number of measurements by the voxel spacing ratio between the subject-specific template and the time point mean resolution
        meanTimePointResolution = 0
        for t in range(self.numberOfTimepoints):
            meanTimePointResolution += np.prod(self.timepointModels[t].voxelSpacing)
        meanTimePointResolution /= self.numberOfTimepoints
        voxelSpacingRatio = np.prod(self.sstModel.voxelSpacing) / meanTimePointResolution
        print("Voxel spacing ratio: " + str(voxelSpacingRatio))

        K0 = self.sstModel.modelSpecifications.K  # Stiffness population -> latent position
        K1 = self.strengthOfLatentDeformationHyperprior * K0  # Stiffness latent position -> each time point
        sstEstimatedNumberOfVoxelsPerGaussian = np.sum(self.sstModel.optimizationHistory[-1]['posteriorsAtEnd'], axis=0) * \
                                                np.prod(self.sstModel.optimizationHistory[-1]['downSamplingFactors'])
        numberOfClasses = len(self.sstModel.gmm.numberOfGaussiansPerClass)
        numberOfGaussians = sum(self.sstModel.gmm.numberOfGaussiansPerClass)
        self.latentMeansNumberOfMeasurements = np.zeros(numberOfGaussians)
        self.latentVariancesNumberOfMeasurements = np.zeros(numberOfGaussians)
        self.latentMixtureWeightsNumberOfMeasurements = np.zeros(numberOfClasses)
        for classNumber in range(numberOfClasses):
            #
            numberOfComponents = self.sstModel.gmm.numberOfGaussiansPerClass[classNumber]
            gaussianNumbers = np.array(np.sum(self.sstModel.gmm.numberOfGaussiansPerClass[:classNumber]) +
                                       np.array(range(numberOfComponents)), dtype=np.uint32)
            sstEstimatedNumberOfVoxelsInClass = np.sum(sstEstimatedNumberOfVoxelsPerGaussian[gaussianNumbers])

            self.latentMixtureWeightsNumberOfMeasurements[
                classNumber] = self.strengthOfLatentGMMHyperprior * sstEstimatedNumberOfVoxelsInClass * voxelSpacingRatio

            averageSizeOfComponents = sstEstimatedNumberOfVoxelsInClass / numberOfComponents
            self.latentMeansNumberOfMeasurements[gaussianNumbers] = self.strengthOfLatentGMMHyperprior * averageSizeOfComponents * voxelSpacingRatio
            self.latentVariancesNumberOfMeasurements[gaussianNumbers] = self.strengthOfLatentGMMHyperprior * averageSizeOfComponents * voxelSpacingRatio

        # Estimating the mode of the latentVariance posterior distribution (which is Wishart) requires a stringent condition
        # on latentVariancesNumberOfMeasurements so that the mode is actually defined
        threshold = (self.sstModel.gmm.numberOfContrasts + 2) + 1e-6
        self.latentVariancesNumberOfMeasurements[self.latentVariancesNumberOfMeasurements < threshold] = threshold

        # No point in updating latent GMM parameters if the GMM hyperprior has zero weight. The latent variances are also
        # a bit tricky, as they're technically driven to zero in that scenario -- let's try not to go there...
        if self.strengthOfLatentGMMHyperprior == 0:
            self.updateLatentMeans, self.updateLatentVariances, self.updateLatentMixtureWeights = False, False, False

        # Loop over all iterations
        self.historyOfTotalCost, self.historyOfTotalTimepointCost, self.historyOfLatentAtlasCost = [], [], []
        progressPlot = None
        iterationNumber = 0
        if self.saveHistory:
            self.history = {**self.history,
                       **{
                           "timepointMeansEvolution": [],
                           "timepointVariancesEvolution": [],
                           "timepointMixtureWeightsEvolution": [],
                           "timepointBiasFieldCoefficientsEvolution": [],
                           "timepointDeformationsEvolution": [],
                           "timepointDeformationAtlasFileNamesEvolution": [],
                           "latentMeansEvolution": [],
                           "latentVariancesEvolution": [],
                           "latentMixtureWeightsEvolution": [],
                           "latentDeformationEvolution": [],
                           "latentDeformationAtlasFileNameEvolution": []
                       }
                       }

        # Make latent atlas directory
        latentAtlasDirectory = os.path.join(self.savePath, 'latentAtlases')
        os.makedirs(latentAtlasDirectory, exist_ok=True)

        while True:

            # =======================================================================================
            #
            # Update parameters for each time point using the current latent variable estimates
            #
            # =======================================================================================

            # Create a new atlas that will be the basis to deform the individual time points from
            latentAtlasFileName = os.path.join(latentAtlasDirectory, 'latentAtlas_iteration_%02d.mgz' % (iterationNumber + 1))
            self.probabilisticAtlas.saveDeformedAtlas(self.latentDeformationAtlasFileName, latentAtlasFileName, self.latentDeformation, True)

            # Only use the last resolution level, and with the newly created atlas as atlas
            for timepointNumber in range(self.numberOfTimepoints):
                self.timepointModels[timepointNumber].optimizationOptions = self.sstModel.optimizationOptions
                self.timepointModels[timepointNumber].optimizationOptions['multiResolutionSpecification'] = [
                    self.timepointModels[timepointNumber].optimizationOptions['multiResolutionSpecification'][-1]]
                self.timepointModels[timepointNumber].optimizationOptions['multiResolutionSpecification'][0]['atlasFileName'] = latentAtlasFileName
                print(self.timepointModels[timepointNumber].optimizationOptions)

            # Loop over all time points
            totalTimepointCost = 0
            for timepointNumber in range(self.numberOfTimepoints):
                #
                self.timepointModels[timepointNumber].modelSpecifications.K = K1
                self.timepointModels[timepointNumber].gmm.hyperMeans = self.latentMeans
                self.timepointModels[timepointNumber].gmm.hyperVariances = self.latentVariances
                self.timepointModels[timepointNumber].gmm.hyperMixtureWeights = self.latentMixtureWeights
                self.timepointModels[timepointNumber].gmm.fullHyperMeansNumberOfMeasurements = self.latentMeansNumberOfMeasurements.copy()
                self.timepointModels[timepointNumber].gmm.fullHyperVariancesNumberOfMeasurements = self.latentVariancesNumberOfMeasurements.copy()
                self.timepointModels[timepointNumber].gmm.fullHyperMixtureWeightsNumberOfMeasurements = self.latentMixtureWeightsNumberOfMeasurements.copy()

                self.timepointModels[timepointNumber].estimateModelParameters(
                                            initialBiasFieldCoefficients=self.timepointModels[timepointNumber].biasField.coefficients,
                                            initialDeformation=self.timepointModels[timepointNumber].deformation,
                                            initialDeformationAtlasFileName=self.timepointModels[timepointNumber].deformationAtlasFileName,
                                            skipBiasFieldParameterEstimationInFirstIteration=False,
                                            skipGMMParameterEstimationInFirstIteration=(iterationNumber == 0)
                                            )

                totalTimepointCost += self.timepointModels[timepointNumber].optimizationHistory[-1]['historyOfCost'][-1]

                print('=================================')
                print('\n')
                print('timepointNumber: ', timepointNumber)
                print('perVoxelCost: ', self.timepointModels[timepointNumber].optimizationSummary[-1]['perVoxelCost'])
                print('\n')
                print('=================================')
                if hasattr(self.visualizer, 'show_flag'):
                    import matplotlib.pyplot as plt  # avoid importing matplotlib by default
                    plt.ion()
                    self.timepointModels[timepointNumber].biasField.downSampleBasisFunctions([1, 1, 1])
                    timepointBiasFields = self.timepointModels[timepointNumber].biasField.getBiasFields(self.masks[timepointNumber])
                    timepointData = self.imageBuffersList[timepointNumber][self.masks[timepointNumber], :] - timepointBiasFields[self.masks[timepointNumber], :]
                    for contrastNumber in range(self.sstModel.gmm.numberOfContrasts):
                        axs = axsList[contrastNumber]
                        ax = axs.ravel()[2 + timepointNumber]
                        ax.clear()
                        ax.hist(timepointData[:, contrastNumber], bins)
                        ax.grid()
                        ax.set_title('time point ' + str(timepointNumber))
                    plt.draw()

                # End loop over time points

            # =======================================================================================
            #
            # Check for convergence.
            # =======================================================================================

            # In order to also measure the deformation from the population atlas -> latent position,
            # create:
            #   (1) a mesh collection with as reference position the population reference position, and as positions
            #       the currently estimated time point positions.
            #   (2) a mesh with the current latent position
            # Note that in (1) we don't need those time positions now, but these will come in handy very soon to
            # optimize the latent position
            #
            # The parameter estimation happens in a (potentially) downsampled image grid, so it's import to work in the same space
            # when measuring and updating the latentDeformation

            transformUsedForEstimation = gems.KvlTransform(
                requireNumpyArray(self.sstModel.optimizationHistory[-1]['downSampledTransformMatrix']))
            mesh_collection = gems.KvlMeshCollection()
            mesh_collection.read(self.latentDeformationAtlasFileName)
            mesh_collection.transform(transformUsedForEstimation)
            referencePosition = mesh_collection.reference_position
            timepointPositions = []
            for timepointNumber in range(self.numberOfTimepoints):
                positionInTemplateSpace = self.probabilisticAtlas.mapPositionsFromSubjectToTemplateSpace(referencePosition,
                                                                                 transformUsedForEstimation) + \
                                          self.latentDeformation + self.timepointModels[timepointNumber].deformation
                timepointPositions.append(
                    self.probabilisticAtlas.mapPositionsFromTemplateToSubjectSpace(positionInTemplateSpace, transformUsedForEstimation))
            mesh_collection.set_positions(referencePosition, timepointPositions)


            # Read mesh in sst warp
            mesh = self.probabilisticAtlas.getMesh(latentAtlasFileName, transformUsedForEstimation)

            #
            calculator = gems.KvlCostAndGradientCalculator(mesh_collection, K0, 0.0, transformUsedForEstimation)
            latentAtlasCost, _ = calculator.evaluate_mesh_position(mesh)

            #
            totalCost = totalTimepointCost + latentAtlasCost
            print('*' * 100 + '\n')
            print('iterationNumber: ', iterationNumber)
            print('totalCost: ', totalCost)
            print('   latentAtlasCost: ', latentAtlasCost)
            print('   totalTimepointCost: ', totalTimepointCost)
            print('*' * 100 + '\n')
            self.historyOfTotalCost.append(totalCost)
            self.historyOfTotalTimepointCost.append(totalTimepointCost)
            self.historyOfLatentAtlasCost.append(latentAtlasCost)

            if hasattr(self.visualizer, 'show_flag'):
                import matplotlib.pyplot as plt  # avoid importing matplotlib by default
                plt.ion()
                if progressPlot is None:
                    plt.figure()
                    progressPlot = plt.subplot()
                progressPlot.clear()
                progressPlot.plot(self.historyOfTotalCost, color='k')
                progressPlot.plot(self.historyOfTotalTimepointCost, linestyle='-.', color='b')
                progressPlot.plot(self.historyOfLatentAtlasCost, linestyle='-.', color='r')
                progressPlot.grid()
                progressPlot.legend(['total', 'timepoints', 'latent atlas deformation'])
                plt.draw()

            if self.saveHistory:
                self.history["timepointMeansEvolution"].append(self.timepointModels[timepointNumber].gmm.means.copy())
                self.history["timepointVariancesEvolution"].append(self.timepointModels[timepointNumber].gmm.variances.copy())
                self.history["timepointMixtureWeightsEvolution"].append(self.timepointModels[timepointNumber].gmm.mixtureWeights.copy())
                self.history["timepointBiasFieldCoefficientsEvolution"].append(self.timepointModels[timepointNumber].biasField.coefficients.copy())
                self.history["timepointDeformationsEvolution"].append(self.timepointModels[timepointNumber].deformation)
                self.history["timepointDeformationAtlasFileNamesEvolution"].append(self.timepointModels[timepointNumber].deformationAtlasFileName)
                self.history["latentMeansEvolution"].append(self.latentMeans.copy())
                self.history["latentVariancesEvolution"].append(self.latentVariances.copy())
                self.history["latentMixtureWeightsEvolution"].append(self.latentMixtureWeights.copy())
                self.history["latentDeformationEvolution"].append(self.latentDeformation.copy())
                self.history["latentDeformationAtlasFileNameEvolution"].append(self.latentDeformationAtlasFileName)

            if iterationNumber >= (self.numberOfIterations - 1):
                print('Stopping')
                break

            # =======================================================================================
            #
            # Update the latent variables based on the current parameter estimates
            #
            # =======================================================================================
            self.updateLatentDeformationAtlas(mesh_collection, mesh, K0, K1, transformUsedForEstimation)
            self.updateGMMLatentVariables()

            iterationNumber += 1
            # End loop over parameter and latent variable estimation iterations

    def updateLatentDeformationAtlas(self, mesh_collection, mesh, K0, K1, transformUsedForEstimation):
        # Update the latentDeformation
        if self.updateLatentDeformation:
            # Set up calculator
            calculator = gems.KvlCostAndGradientCalculator(mesh_collection, K0, K1, transformUsedForEstimation)

            # Get optimizer and plug calculator in it
            optimizerType = 'L-BFGS'
            optimizationParameters = {
                'Verbose': self.sstModel.optimizationOptions['verbose'],
                'MaximalDeformationStopCriterion': self.sstModel.optimizationOptions['maximalDeformationStopCriterion'],
                'LineSearchMaximalDeformationIntervalStopCriterion': self.sstModel.optimizationOptions[
                    'lineSearchMaximalDeformationIntervalStopCriterion'],
                'MaximumNumberOfIterations': self.sstModel.optimizationOptions['maximumNumberOfDeformationIterations'],
                'BFGS-MaximumMemoryLength': self.sstModel.optimizationOptions['BFGSMaximumMemoryLength']
            }
            optimizer = gems.KvlOptimizer(optimizerType, mesh, calculator, optimizationParameters)

            # Run deformation optimization
            historyOfDeformationCost = []
            historyOfMaximalDeformation = []
            nodePositionsBeforeDeformation = mesh.points
            while True:
                minLogLikelihoodTimesDeformationPrior, maximalDeformation = optimizer.step_optimizer_samseg()
                print("maximalDeformation=%.4f minLogLikelihood=%.4f" % (
                    maximalDeformation, minLogLikelihoodTimesDeformationPrior))
                historyOfDeformationCost.append(minLogLikelihoodTimesDeformationPrior)
                historyOfMaximalDeformation.append(maximalDeformation)
                if maximalDeformation == 0:
                    break

            nodePositionsAfterDeformation = mesh.points
            maximalDeformationApplied = np.sqrt(
                np.max(np.sum((nodePositionsAfterDeformation - nodePositionsBeforeDeformation) ** 2, 1)))

            #
            nodePositionsBeforeDeformation = self.probabilisticAtlas.mapPositionsFromSubjectToTemplateSpace(
                nodePositionsBeforeDeformation,
                transformUsedForEstimation)
            nodePositionsAfterDeformation = self.probabilisticAtlas.mapPositionsFromSubjectToTemplateSpace(
                nodePositionsAfterDeformation,
                transformUsedForEstimation)
            estimatedUpdate = nodePositionsAfterDeformation - nodePositionsBeforeDeformation
            self.latentDeformation += estimatedUpdate

    def updateGMMLatentVariables(self):

        # Update latentMeans
        if self.updateLatentMeans:
            numberOfGaussians = np.sum(self.sstModel.gmm.numberOfGaussiansPerClass)
            numberOfContrasts = self.latentMeans.shape[-1]
            for gaussianNumber in range(numberOfGaussians):
                # Set up linear system
                lhs = np.zeros((numberOfContrasts, numberOfContrasts))
                rhs = np.zeros((numberOfContrasts, 1))
                for timepointNumber in range(self.numberOfTimepoints):
                    mean = np.expand_dims(self.timepointModels[timepointNumber].gmm.means[gaussianNumber], 1)
                    variance = self.timepointModels[timepointNumber].gmm.variances[gaussianNumber]

                    lhs += np.linalg.inv(variance)
                    rhs += np.linalg.solve(variance, mean)

                # Solve linear system
                latentMean = np.linalg.solve(lhs, rhs)
                self.latentMeans[gaussianNumber, :] = latentMean.T

        # Update latentVariances
        if self.updateLatentVariances:
            numberOfGaussians = np.sum(self.sstModel.gmm.numberOfGaussiansPerClass)
            numberOfContrasts = self.latentMeans.shape[-1]
            for gaussianNumber in range(numberOfGaussians):
                # Precision is essentially averaged
                averagePrecision = np.zeros((numberOfContrasts, numberOfContrasts))
                for timepointNumber in range(self.numberOfTimepoints):
                    variance = self.timepointModels[timepointNumber].gmm.variances[gaussianNumber]

                    averagePrecision += np.linalg.inv(variance)
                averagePrecision /= self.numberOfTimepoints

                latentVarianceNumberOfMeasurements = self.latentVariancesNumberOfMeasurements[gaussianNumber]
                latentVariance = np.linalg.inv(averagePrecision) * \
                                 (latentVarianceNumberOfMeasurements - numberOfContrasts - 2) / latentVarianceNumberOfMeasurements
                self.latentVariances[gaussianNumber] = latentVariance

        # Update latentMixtureWeights
        if self.updateLatentMixtureWeights:
            numberOfClasses = len(self.sstModel.gmm.numberOfGaussiansPerClass)
            for classNumber in range(numberOfClasses):
                numberOfComponents = self.sstModel.gmm.numberOfGaussiansPerClass[classNumber]
                averageInLogDomain = np.zeros(numberOfComponents)
                for componentNumber in range(numberOfComponents):
                    gaussianNumber = sum(self.sstModel.gmm.numberOfGaussiansPerClass[:classNumber]) + componentNumber
                    for timepointNumber in range(self.numberOfTimepoints):
                        mixtureWeight = self.timepointModels[timepointNumber].gmm.mixtureWeights[gaussianNumber]
                        averageInLogDomain[componentNumber] += np.log(mixtureWeight + eps)
                    averageInLogDomain[componentNumber] /= self.numberOfTimepoints

                # Solution is normalized version
                solution = np.exp(averageInLogDomain)
                solution /= np.sum(solution + eps)

                #
                for componentNumber in range(numberOfComponents):
                    gaussianNumber = sum(self.sstModel.gmm.numberOfGaussiansPerClass[:classNumber]) + componentNumber
                    self.latentMixtureWeights[gaussianNumber] = solution[componentNumber]

    def postProcess(self, saveWarp=False):

        # =======================================================================================
        #
        # Using estimated parameters, segment and write out results for each time point
        #
        # =======================================================================================
        #

        sstDir = os.path.join(self.savePath, 'base')
        os.makedirs(sstDir, exist_ok=True)
        baseModel = self.sstModel;
        # Save the final mesh collection
        if self.saveModelProbabilities:
            print('Saving base model probs')
            baseModel.saveGaussianProbabilities(os.path.join(sstDir, 'probabilities') )
        if saveWarp:
            baseModel.saveWarpField(os.path.join(sstDir, 'template.m3z'))

        self.timepointVolumesInCubicMm = []
        for timepointNumber in range(self.numberOfTimepoints):
            timepointModel = self.timepointModels[timepointNumber]
            timepointModel.deformation = self.latentDeformation + timepointModel.deformation
            timepointModel.deformationAtlasFileName = self.latentDeformationAtlasFileName
            posteriors, biasFields, nodePositions, _, _ = timepointModel.computeFinalSegmentation()
            timepointDir = os.path.join(self.savePath, 'tp%03d' % (timepointNumber + 1))
            os.makedirs(timepointDir, exist_ok=True)
            timepointModel.savePath = timepointDir
            volumesInCubicMm = timepointModel.writeResults(biasFields, posteriors)

            # Save the timepoint->template warp
            if saveWarp:
                timepointModel.saveWarpField(os.path.join(timepointDir, 'template.m3z'))

            # Save the final mesh collection
            if self.saveMesh:
                print('Saving the final mesh in template space')
                deformedAtlasFileName = os.path.join(timepointModel.savePath, 'mesh.txt')
                timepointModel.probabilisticAtlas.saveDeformedAtlas(timepointModel.modelSpecifications.atlasFileName,
                                                                    deformedAtlasFileName, nodePositions)

            if self.saveModelProbabilities:
                print('Saving model probs')
                timepointModel.saveGaussianProbabilities( os.path.join(timepointModel.savePath, 'probabilities') )

            # Save the history of the parameter estimation process
            if self.saveHistory:
                history = {'input': {
                    'imageFileNames': timepointModel.imageFileNames,
                    'imageToImageTransformMatrix': timepointModel.imageToImageTransformMatrix,
                    'modelSpecifications': timepointModel.modelSpecifications,
                    'optimizationOptions': timepointModel.optimizationOptions,
                    'savePath': timepointModel.savePath
                }, 'imageBuffers': timepointModel.imageBuffers, 'mask': timepointModel.mask,
                    'cropping': timepointModel.cropping,
                    'transform': timepointModel.transform.as_numpy_array,
                    'historyWithinEachMultiResolutionLevel': timepointModel.optimizationHistory,
                    "labels": timepointModel.modelSpecifications.FreeSurferLabels, "names": timepointModel.modelSpecifications.names,
                    "volumesInCubicMm": volumesInCubicMm, "optimizationSummary": timepointModel.optimizationSummary}
                with open(os.path.join(timepointModel.savePath, 'history.p'), 'wb') as file:
                    pickle.dump(history, file, protocol=pickle.HIGHEST_PROTOCOL)
            

            self.timepointVolumesInCubicMm.append(volumesInCubicMm)

        #
        self.optimizationSummary = {
            "historyOfTotalCost": self.historyOfTotalCost,
            "historyOfTotalTimepointCost": self.historyOfTotalTimepointCost,
            "historyOfLatentAtlasCost": self.historyOfLatentAtlasCost
        }

        if self.saveHistory:
            self.history["labels"] = self.sstModel.modelSpecifications.FreeSurferLabels
            self.history["names"] = self.sstModel.modelSpecifications.names
            self.history["timepointVolumesInCubicMm"] = self.timepointVolumesInCubicMm
            self.history["optimizationSummary"] = self.optimizationSummary
            with open(os.path.join(self.savePath, 'history.p'), 'wb') as file:
                pickle.dump(self.history, file, protocol=pickle.HIGHEST_PROTOCOL)

    def generateSubjectSpecificTemplate(self,saveWarp=False):
        sstDir = os.path.join(self.savePath, 'base')
        os.makedirs(sstDir, exist_ok=True)

        sstFileNames = []
        for contrastNumber, contrastImageFileNames in enumerate(zip(*self.imageFileNamesList)):

            # Read in the various time point images, and compute the average
            numberOfTimepoints = len(contrastImageFileNames)
            image0 = sf.load_volume(contrastImageFileNames[0])
            imageBuffer = image0.transform(affine=self.tpToBaseTransforms[0])
            # Make sure that we are averaging only non zero voxels
            count = np.zeros(imageBuffer.shape)
            count[imageBuffer > 0] += 1
            for timepointNumber in range(1, numberOfTimepoints):
                tmp = sf.load_volume(contrastImageFileNames[timepointNumber]).transform(affine=self.tpToBaseTransforms[timepointNumber]).data
                imageBuffer += tmp
                count[tmp > 0] += 1
            # Make sure that we are not dividing by zero for, e.g., background voxels
            imageBuffer[count > 0] /= count[count > 0]

            # Write image to disk
            sstFilename = os.path.join(sstDir, 'mode%02d_average.mgz' % (contrastNumber + 1))
            imageBuffer.save(sstFilename)

            sstFileNames.append(sstFilename)

        return sstFileNames

    def constructSstModel(self):

        sstDir, _ = os.path.split(self.sstFileNames[0])

        self.sstModel = Samseg(
            imageFileNames=self.sstFileNames,
            atlasDir=self.atlasDir,
            savePath=sstDir,
            imageToImageTransformMatrix=self.imageToImageTransformMatrix,
            userModelSpecifications=self.userModelSpecifications,
            userOptimizationOptions=self.userOptimizationOptions,
            visualizer=self.visualizer,
            saveHistory=True,
            savePosteriors=self.savePosteriors,
            targetIntensity=self.targetIntensity,
            targetSearchStrings=self.targetSearchStrings,
            modeNames=self.modeNames,
            pallidumAsWM=self.pallidumAsWM
        )

    def constructTimepointModels(self):

        self.timepointModels = []

        # Construction of the cross sectional model for each time point
        for timepointNumber in range(self.numberOfTimepoints):
            self.timepointModels.append(Samseg(
                imageFileNames=self.imageFileNamesList[timepointNumber],
                atlasDir=self.atlasDir,
                savePath=self.savePath,
                imageToImageTransformMatrix=self.imageToImageTransformMatrix,
                userModelSpecifications=self.userModelSpecifications,
                userOptimizationOptions=self.userOptimizationOptions,
                visualizer=self.visualizer,
                saveHistory=True,
                targetIntensity=self.targetIntensity,
                targetSearchStrings=self.targetSearchStrings,
                modeNames=self.modeNames,
                pallidumAsWM=self.pallidumAsWM,
                savePosteriors=self.savePosteriors,
            ))

            self.timepointModels[timepointNumber].mask = self.masks[timepointNumber]
            self.timepointModels[timepointNumber].imageBuffers = self.imageBuffersList[timepointNumber]
            self.timepointModels[timepointNumber].voxelSpacing = self.voxelSpacings[timepointNumber]
            self.timepointModels[timepointNumber].transform = self.transforms[timepointNumber]
            self.timepointModels[timepointNumber].cropping = self.croppings[timepointNumber]
