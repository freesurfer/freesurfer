import numpy as np
import os

from freesurfer.samseg.Samseg import Samseg
from freesurfer.samseg.utilities import Specification
from freesurfer.samseg.SamsegUtility import *
from freesurfer.samseg.io import kvlReadSharedGMMParameters
from freesurfer.samseg.VAE import VAE
from freesurfer.samseg.merge_alphas import kvlGetMergingFractionsTable
from . import gems

eps = np.finfo(float).eps


class SamsegLesion(Samseg):
    def __init__(self, imageFileNames, atlasDir, savePath, userModelSpecifications={}, userOptimizationOptions={},
                 imageToImageTransformMatrix=None, visualizer=None, saveHistory=None, savePosteriors=None,
                 saveWarp=None, saveMesh=None, threshold=0.3, thresholdSearchString='Lesion',
                 targetIntensity=None, targetSearchStrings=None, modeNames=None, pallidumAsWM=True,
                 saveModelProbabilities=False,
                 numberOfSamplingSteps=50, numberOfBurnInSteps=50,
                 numberOfPseudoSamplesMean=500, numberOfPseudoSamplesVariance=500, rho=50,
                 intensityMaskingPattern=None, intensityMaskingSearchString='Cortex', gmmFileName=None, sampler=True,
                 ignoreUnknownPriors=False,
                 ):
        Samseg.__init__(self, imageFileNames, atlasDir, savePath, userModelSpecifications, userOptimizationOptions,
                 imageToImageTransformMatrix, visualizer, saveHistory, savePosteriors,
                 saveWarp, saveMesh, threshold, thresholdSearchString,
                 targetIntensity, targetSearchStrings, modeNames, pallidumAsWM=pallidumAsWM,
                 saveModelProbabilities=saveModelProbabilities, gmmFileName=gmmFileName,
                 ignoreUnknownPriors=ignoreUnknownPriors)
        self.numberOfSamplingSteps = numberOfSamplingSteps
        self.numberOfBurnInSteps = numberOfBurnInSteps
        self.numberOfPseudoSamplesMean = numberOfPseudoSamplesMean
        self.numberOfPseudoSamplesVariance = numberOfPseudoSamplesVariance
        self.rho = rho
        self.intensityMaskingClassNumber = self.getClassNumber(intensityMaskingSearchString)
        self.sampler = sampler

        if intensityMaskingPattern is None:
            raise ValueError('Intensity mask pattern must be set')
        if len(intensityMaskingPattern) != len(imageFileNames):
            raise ValueError('Number of lesion mask patterns does not match the number of input images.')
        if not(all(pattern in (0, 1, -1) for pattern in intensityMaskingPattern)):
            raise ValueError('Lesion mask pattern values can be only 0, 1 or -1')

        self.intensityMaskingPattern = intensityMaskingPattern

        # Check conditions on white matter and lesion gaussian/structure and
        # get their structure numbers, class number as well as the gaussian number
        wmSearchString = 'White'
        lesionSearchString = self.thresholdSearchString
        self.lesionStructureNumber, self.lesionClassNumber, self.lesionGaussianNumber = self.checkConditions(lesionSearchString)
        _, _, self.wmGaussianNumber = self.checkConditions(wmSearchString, checkStructureOwnClass=False)

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
        # Here we tied the lesion gaussian to the wm one
        self.gmm.fullHyperMeansNumberOfMeasurements[self.lesionGaussianNumber] = self.numberOfPseudoSamplesMean
        self.gmm.fullHyperVariancesNumberOfMeasurements[self.lesionGaussianNumber] = self.numberOfPseudoSamplesVariance
        self.gmm.tiedGaussiansInit(self.wmGaussianNumber, self.lesionGaussianNumber, self.rho)

    def computeFinalSegmentation(self):

        posteriors, biasFields, nodePositions, data, priors = Samseg.computeFinalSegmentation(self)

        # If no sampler return the segmentation computed by Samseg, so that the VAE is not used.
        if not self.sampler:
            return posteriors, biasFields, nodePositions, data, priors

        #
        numberOfVoxels = data.shape[0]
        imageSize = self.mask.shape

        # Since we sample in subject space, the number of pseudo samples in the gmm needs to be updated accordingly
        self.gmm.downsampledHyperparameters(self.voxelSpacing)

        # Create intensity-based lesion mask
        if self.intensityMaskingClassNumber is not None:
            # We have -1 mask below mean, +1 above, 0 nothing
            intensityMaskingGaussianNumbers = [sum(self.gmm.numberOfGaussiansPerClass[0: self.intensityMaskingClassNumber])]
            intensityMaskingMean = np.sum(self.gmm.means[intensityMaskingGaussianNumbers] *
                                          self.gmm.mixtureWeights[intensityMaskingGaussianNumbers].reshape(-1, 1), axis=0)
            tmp = np.ones(numberOfVoxels, dtype=bool)
            for contrastNumber in range(self.gmm.numberOfContrasts):
                direction = self.intensityMaskingPattern[contrastNumber]
                if direction == -1:
                    tmp = np.logical_and(tmp, data[:, contrastNumber] < intensityMaskingMean[contrastNumber])
                elif direction == 1:
                    tmp = np.logical_and(tmp, data[:, contrastNumber] > intensityMaskingMean[contrastNumber])
            intensityMask = np.zeros(imageSize, dtype=bool)
            intensityMask[self.mask] = tmp
        else:
            intensityMask = np.zeros(imageSize, dtype=bool)
            intensityMask[self.mask] = True

        self.visualizer.show(image_list=[intensityMask.astype(float)], title="Intensity mask")

        # Initialize the structure likelihoods from the initial parameter values.
        # Since only the parameters of a single structure will be altered,
        # only one column in the likelihoods will need to be updated during sampling
        likelihoods = self.gmm.getLikelihoods(data, self.classFractions)

        # Initialize the sampler with a majority-vote lesion segmentation, masked with intensity-based mask
        posteriors = likelihoods * priors
        posteriors /= np.expand_dims(np.sum(posteriors, axis=1) + eps, 1)
        lesion = np.zeros(imageSize)
        lesion[self.mask] = (np.array(np.argmax(posteriors, 1), dtype=np.uint32) == self.lesionStructureNumber)
        lesion *= intensityMask

        self.visualizer.show(image_list=[lesion], title="Initial lesion segmentation")

        # Initialize the VAE tensorflow model and its various settings.
        # Restore from checkpoint the VAE
        vae = VAE(self.atlasDir, self.transform, imageSize)

        # Do the actual sampling of lesion, latent variables of the VAE model, and mean/variance of the lesion intensity model.
        averagePosteriors = np.zeros_like(likelihoods)
        self.visualizer.start_movie(window_id="Lesion prior using VAE only", title="Lesion prior using VAE only -- the movie")
        self.visualizer.start_movie(window_id="Lesion sample", title="Lesion sample -- the movie")
        for sweepNumber in range(self.numberOfBurnInSteps + self.numberOfSamplingSteps):

            # Sample from the VAE latent variables, conditioned on the current lesion segmentation.
            # Implementation-wise we don't store the latent variables, but rather the factorized
            # prior in the visible units (voxels) that they encode.
            lesionPriorVAE = (vae.sample(lesion) * intensityMask)[self.mask]

            if hasattr(self.visualizer, 'show_flag'):
                tmp = np.zeros(imageSize)
                tmp[self.mask] = lesionPriorVAE
                self.visualizer.show(probabilities=tmp, title="Lesion prior using VAE only",
                                 window_id="Lesion prior using VAE only")

            # Sample from the mean and variance, conditioned on the data and the lesion segmentation
            mean, variance = self.gmm.sampleMeansAndVariancesConditioned(data, lesion[self.mask].reshape(-1, 1),
                                                                         self.lesionGaussianNumber)

            # Sample from the lesion segmentation, conditioned on the data and the VAE latent variables
            # (Implementation-wise the latter is encoded in the VAE prior). At the same time we also
            # compute the full posterior of each structure, which is at the end the thing we're averaging
            # over (i.e., the reason why we're sampling)
            effectivePriors = np.array(priors / 65535, dtype=np.float32)
            effectivePriors[:, self.lesionStructureNumber] *= lesionPriorVAE
            if hasattr(self.visualizer, 'show_flag'):
                tmp = np.zeros(imageSize)
                tmp[self.mask] = effectivePriors[:, self.lesionStructureNumber]
                self.visualizer.show(probabilities=tmp, title="Lesion prior using VAE and atlas together",
                                     window_id="Lesion prior using VAE and atlas together")

            # Generative model where the atlas generates *candidate* lesions, and the VAE prior is sampled
            # from *only within the candidates*.
            numberOfStructures = priors.shape[-1]
            otherStructureNumbers = [i for i in range(numberOfStructures) if i != self.lesionStructureNumber]
            otherStructurePriors = effectivePriors[:, otherStructureNumbers]
            otherStructurePriors /= (np.sum(otherStructurePriors, axis=1).reshape(-1, 1) + eps)
            effectivePriors[:, otherStructureNumbers] = otherStructurePriors * \
                                                        (1 - effectivePriors[:, self.lesionStructureNumber].reshape(-1, 1))
            likelihoods[:, self.lesionStructureNumber] = self.gmm.getGaussianLikelihoods(data, mean, variance)
            posteriors = effectivePriors * likelihoods
            posteriors /= np.expand_dims(np.sum(posteriors, axis=1) + eps, 1)
            sample = np.random.rand(numberOfVoxels) <= posteriors[:, self.lesionStructureNumber]
            lesion = np.zeros(imageSize)
            lesion[self.mask] = sample

            self.visualizer.show(image_list=[lesion], title="Lesion sample", window_id="Lesion sample")

            # Collect data after burn in steps
            if sweepNumber >= self.numberOfBurnInSteps:
                print('Sample ' + str(sweepNumber + 1 - self.numberOfBurnInSteps) + ' times')
                averagePosteriors += posteriors / self.numberOfSamplingSteps
            else:
                print('Burn-in ' + str(sweepNumber + 1) + ' times')

        #
        self.visualizer.show_movie(window_id="Lesion prior using VAE only")
        self.visualizer.show_movie(window_id="Lesion sample")

        # Return
        return averagePosteriors, biasFields, nodePositions, data, priors
