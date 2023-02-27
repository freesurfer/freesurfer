import os
from .SamsegLesion import SamsegLesion
from .SamsegLongitudinal import SamsegLongitudinal
from .SamsegUtility import *


class SamsegLongitudinalLesion(SamsegLongitudinal):
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
        numberOfSamplingSteps=50,
        numberOfBurnInSteps=50,
        numberOfPseudoSamplesMean=500,
        numberOfPseudoSamplesVariance=500,
        rho=50,
        intensityMaskingPattern=None,
        intensityMaskingSearchString='Cortex',
        tpToBaseTransforms=None,
                 ):
        SamsegLongitudinal.__init__(self,
        imageFileNamesList=imageFileNamesList,
        atlasDir=atlasDir,
        savePath=savePath,
        userModelSpecifications=userModelSpecifications,
        userOptimizationOptions=userOptimizationOptions,
        visualizer=visualizer,
        saveHistory=saveHistory,
        saveMesh=saveMesh,
        targetIntensity=targetIntensity,
        targetSearchStrings=targetSearchStrings,
        numberOfIterations=numberOfIterations,
        strengthOfLatentGMMHyperprior=strengthOfLatentGMMHyperprior,
        strengthOfLatentDeformationHyperprior=strengthOfLatentDeformationHyperprior,
        saveSSTResults=saveSSTResults,
        updateLatentMeans=updateLatentMeans,
        updateLatentVariances=updateLatentVariances,
        updateLatentMixtureWeights=updateLatentMixtureWeights,
        updateLatentDeformation=updateLatentDeformation,
        initializeLatentDeformationToZero=initializeLatentDeformationToZero,
        threshold=threshold,
        thresholdSearchString=thresholdSearchString,
        modeNames=modeNames,
        pallidumAsWM=pallidumAsWM,
        savePosteriors=savePosteriors,
        tpToBaseTransforms=tpToBaseTransforms
        )

        self.numberOfSamplingSteps = numberOfSamplingSteps
        self.numberOfBurnInSteps = numberOfBurnInSteps
        self.numberOfPseudoSamplesMean = numberOfPseudoSamplesMean
        self.numberOfPseudoSamplesVariance = numberOfPseudoSamplesVariance
        self.rho = rho
        self.intensityMaskingSearchString = intensityMaskingSearchString
        self.intensityMaskingPattern = intensityMaskingPattern

    def constructSstModel(self):

        sstDir, _ = os.path.split(self.sstFileNames[0])

        self.sstModel = SamsegLesion(
            imageFileNames=self.sstFileNames,
            atlasDir=self.atlasDir,
            savePath=sstDir,
            imageToImageTransformMatrix=self.imageToImageTransformMatrix,
            userModelSpecifications=self.userModelSpecifications,
            userOptimizationOptions=self.userOptimizationOptions,
            visualizer=self.visualizer,
            saveHistory=True,
            targetIntensity=self.targetIntensity,
            targetSearchStrings=self.targetSearchStrings,
            modeNames=self.modeNames,
            pallidumAsWM=self.pallidumAsWM,
            numberOfSamplingSteps=self.numberOfSamplingSteps,
            numberOfBurnInSteps=self.numberOfBurnInSteps,
            numberOfPseudoSamplesMean=self.numberOfPseudoSamplesMean,
            numberOfPseudoSamplesVariance=self.numberOfPseudoSamplesVariance,
            rho=self.rho,
            intensityMaskingPattern=self.intensityMaskingPattern,
            intensityMaskingSearchString=self.intensityMaskingSearchString,
            sampler=False
        )

    def constructTimepointModels(self):

        self.timepointModels = []

        # Construction of the cross sectional model for each timepoint
        for timepointNumber in range(self.numberOfTimepoints):
            self.timepointModels.append(SamsegLesion(
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
                numberOfSamplingSteps=self.numberOfSamplingSteps,
                numberOfBurnInSteps=self.numberOfBurnInSteps,
                numberOfPseudoSamplesMean=self.numberOfPseudoSamplesMean,
                numberOfPseudoSamplesVariance=self.numberOfPseudoSamplesVariance,
                rho=self.rho,
                intensityMaskingPattern=self.intensityMaskingPattern,
                intensityMaskingSearchString=self.intensityMaskingSearchString
            ))
            self.timepointModels[timepointNumber].mask = self.sstModel.mask
            self.timepointModels[timepointNumber].imageBuffers = self.imageBuffersList[timepointNumber]
            self.timepointModels[timepointNumber].voxelSpacing = self.sstModel.voxelSpacing
            self.timepointModels[timepointNumber].transform = self.sstModel.transform
            self.timepointModels[timepointNumber].cropping = self.sstModel.cropping

    def initializeLatentVariables(self):

        # First call parent function to initialize all the latent variables
        K0, K1 = SamsegLongitudinal.initializeLatentVariables(self)

        # Now override the lesion latent variables to the WM ones
        self.setLesionLatentVariables()

        return K0, K1

    def updateGMMLatentVariables(self):
        # First call parent function to initialize all the latent variables
        SamsegLongitudinal.updateGMMLatentVariables(self)

        # Now override the lesion latent variables to the WM ones
        self.setLesionLatentVariables()

    def setLesionLatentVariables(self):
        self.latentMeans[self.sstModel.lesionGaussianNumber] = self.sstModel.gmm.means[self.sstModel.wmGaussianNumber]
        self.latentVariances[self.sstModel.lesionGaussianNumber] = self.rho * self.sstModel.gmm.variances[self.sstModel.wmGaussianNumber]
        self.latentMeansNumberOfMeasurements[self.sstModel.lesionGaussianNumber] = self.numberOfPseudoSamplesMean
        self.latentVariancesNumberOfMeasurements[self.sstModel.lesionGaussianNumber] = self.numberOfPseudoSamplesVariance
