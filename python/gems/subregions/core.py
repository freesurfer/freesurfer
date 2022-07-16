import os
import shutil
import tempfile
import numpy as np
import surfa as sf
import scipy.ndimage

import gems
from gems.subregions import utils


class MeshModel:

    def __init__(
        self,
        atlasDir,
        outDir,
        inputImageFileNames,
        inputSegFileName,
        meshStiffness=0.05,
        optimizerType='L-BFGS',
        bbregisterMode=None,
        resolution=0.5,
        useTwoComponents=False,
        tempDir=None,
        fileSuffix='',
        debug=False,
        ):
        """
        MeshModel is a generic base class to facilitate GEMS mesh deformation for given ROIs.
        To implement a mesh model for a particular set of structures, this class must be subclassed
        and the following functions MUST be implemented:

        self.preprocess_images()         : Precompute the mask, segmentation, and image volumes
        self.get_cheating_label_groups() : Return the reduced labels group used to fit the mesh to the initial segmentation
        self.get_cheating_gaussians()    : Return the Gaussian means and variances used to fit the mesh to the initial segmentation
        self.get_label_groups()          : Return the reduced labels group used to fit the mesh to the image
        self.get_gaussian_hyps()         : Return the hyperparameters used to estimate the Gaussian parameters during image-fitting
        self.postprocess_segmentation()  : Update and write the label volumes and discrete segmentation(s)

        Further information is documented in each function definition.

        This is a framework that was meant to facilitate an (almost prefect) port of the subfield matlab code.
        """

        # Set some paths
        self.outDir = outDir
        self.atlasDir = atlasDir
        self.atlasMeshFileName = os.path.join(atlasDir, 'AtlasMesh.gz')
        self.atlasDumpFileName = os.path.join(atlasDir, 'AtlasDump.mgz')
        self.compressionLookupTableFileName = os.path.join(atlasDir, 'compressionLookupTable.txt')
        self.inputImageFileNames = inputImageFileNames
        self.inputSegFileName = inputSegFileName

        # Some settings
        self.meshStiffness = meshStiffness
        self.optimizerType = optimizerType
        self.bbregisterMode = bbregisterMode
        self.resolution = resolution
        self.useTwoComponents = useTwoComponents
        self.tempDir = tempDir
        self.fileSuffix = fileSuffix
        self.debug = debug

        # Some optimization defaults that should be overwritten by each subclass
        self.cheatingMeshSmoothingSigmas = [3.0, 2.0]
        self.cheatingMaxIterations = [300, 150]

        self.meshSmoothingSigmas = [1.5, 0.75, 0]
        self.imageSmoothingSigmas = [0, 0, 0]
        self.maxIterations = [7, 5, 3]

        self.isLong = False
        self.longMeshSmoothingSigmas = [[1.5, 0.75], [0.75, 0]]
        self.longImageSmoothingSigmas = [[0, 0], [0, 0]]
        self.longMaxIterations = [[6, 3], [2, 1]]
        self.maxGlobalLongIterations = 2
        self.longMask = None

        # Here are some options that control how to much to dilate masks throughout different
        # stages. Might be necessary to tune depending on the geometry of the ROI (like brainstem).
        self.atlasTargetSmoothing = 'forward'
        self.cheatingAlphaMaskStrel = 3
        self.alphaMaskStrel = 5

    def cleanup(self):
        """
        Essentially the only thing to do during cleanup is (potentially)
        remove the temporary directory.
        """
        if not self.debug:
            shutil.rmtree(self.tempDir)
        else:
            print(f'Not removing temporary directory: {self.tempDir}')

    def initialize(self):
        """
        Initialize the mesh model by running sanity checks on the input options,
        loading input volumes, creating the temporary directory, and doing image preprocessing.
        """

        # First thing: set up the temporary directory. IO to this space should be relatively
        # limited unless debug mode is enabled
        if self.tempDir is None:
            self.tempDir = tempfile.mkdtemp()
        else:
            os.makedirs(self.tempDir, exist_ok=True)

        # Make sure the output directory exists as well
        os.makedirs(self.outDir, exist_ok=True)

        # Sanity check on the optimizer type
        optimizerTypes = ['FixedStepGradientDescent', 'GradientDescent', 'ConjugateGradient', 'L-BFGS']
        if self.optimizerType not in optimizerTypes:
            sf.system.fatal('Optimizer type must be one of: ' + ', '.join(optimizerTypes))

        # Sanity check on the registration mode for alternative images
        bbregisterModes = [None, 't1', 't2']
        if self.bbregisterMode not in bbregisterModes:
            sf.system.fatal('BBregister mode must be one of: ' + ', '.join(bbregisterModes))

        # Make sure all the atlas files are there
        if not os.path.isfile(self.atlasMeshFileName):
            sf.system.fatal(f'Provided atlas mesh file `{self.atlasMeshFileName}` does not exist.')
        if not os.path.isfile(self.atlasDumpFileName):
            sf.system.fatal(f'Provided atlas image `{self.atlasDumpFileName}` does not exist.')
        if not os.path.isfile(self.compressionLookupTableFileName):
            sf.system.fatal(f'Provided compression LUT `{self.compressionLookupTableFileName}` does not exist.')

        # Load compressed and FreeSurfer label mapping information
        self.labelMapping, self.names, self.FreeSurferLabels = utils.read_compression_lookup_table(self.compressionLookupTableFileName)

        # Set the target mesh file paths
        self.warpedMeshFileName = os.path.join(self.tempDir, 'warpedOriginalMesh.txt')
        self.warpedMeshNoAffineFileName = os.path.join(self.tempDir, 'warpedOriginalMeshNoAffine.txt')

        # Read the input volumes (images and reference segmentation)
        self.inputSeg = sf.load_volume(self.inputSegFileName)
        self.inputImages = [sf.load_volume(path) for path in self.inputImageFileNames]
        self.correctedImages = [img.copy() for img in self.inputImages]
        self.highResImage = np.mean(self.inputImages[0].geom.voxsize) < 0.99

        # Now we define a set of volume members that must be properly computed during
        # the `preprocess_images` stage of all MeshModel subclasses. Further documentation below.
        self.preprocess_images()

    def preprocess_images(self):
        """
        Preprocess the input images for later processing. This function must be redefined in a subclass,
        and the following volumes (at the minimum) must be set during this stage:

            1. self.atlasAlignmentTarget : A binary tissue mask that acts as the target for the initial
                                          affine atlas registration.
            2. self.synthImage : A synthetic image generated from the input segmentation, used for the initial
                                fitting of the mesh to the subject (the `cheating` step).
            3. self.processedImage: An image (or set of images represented by each frame) used for the primary
                                    mesh fitting. It is expected that this image has been properly resampled
                                    to the working target resolution.

        It is expected that these are surfa.Volume objects with proper geometry information.
        """
        raise NotImplementedError('All subclasses of MeshModel must implement the preprocess_images() function!')

    def label_group_names_to_indices(self, labelNames):
        """
        Clean and convert a group of label names (list of lists) to a grouping of label indices.
        """
        labelIndices = [[self.labelMapping.search(name, exact=True) for name in group] for group in labelNames]
        labelIndices = [[i for i in group if i is not None] for group in labelIndices]
        labelIndices = [g for g in labelIndices if g]
        return labelIndices

    def reduce_alphas(self, sameGaussianParameters, alphas=None):
        """
        Compute a set of reduced alpha values given groups of labels. Will use the original
        alpha values if alphas is None.
        """
        if alphas is None:
            alphas = self.originalAlphas

        numberOfReducedLabels = len(sameGaussianParameters)
        # ATH: are alphas always 32-bit floats?
        reducedAlphas = np.zeros((alphas.shape[0], numberOfReducedLabels), dtype='float32')
        reducingLookupTable = np.zeros(alphas.shape[1], dtype='int32')

        # Convert to list so we can use index
        fslabels = list(self.FreeSurferLabels)

        # Reduce the labels
        for reducedLabel in range(numberOfReducedLabels):
            sameGaussians = sameGaussianParameters[reducedLabel]
            for label in sameGaussians:
                compressedLabel = fslabels.index(label)
                reducedAlphas[:, reducedLabel] += alphas[:, compressedLabel]
                reducingLookupTable[compressedLabel] = reducedLabel

        # Make sure classes sum to one
        if np.max(np.abs(np.sum(reducedAlphas, -1) - 1)) > 1e-5:
            sf.system.fatal('The vector of prior probabilities in the mesh nodes must always sum to one over all classes')

        return (reducedAlphas, reducingLookupTable)

    def crop_image_by_atlas(self, image):
        """
        Crop image to the aligned atlas image. Also construct a 3-D affine transformation that will later be used
        to transform the location of the atlas mesh's nodes into the coordinate system of the image.
        """
        trf = image.geom.world2vox @ self.alignedAtlas.geom.vox2world
        template_corners = np.mgrid[:2, :2, :2].T.reshape(-1, 3) * (np.array(self.alignedAtlas.shape[:3]) - 1)
        corners = trf.transform(template_corners)
        lower = corners.min(0).astype(int)
        upper = (corners.max(0) + 1).astype(int)

        image_limit = np.array(image.shape[:3]) - 1
        lower = np.clip(lower, (0, 0, 0), image_limit)
        upper = np.clip(upper, (0, 0, 0), image_limit)

        trf.matrix[:3, -1] -= lower
        cropping = tuple([slice(l, u + 1) for l, u in zip(lower, upper)])

        transform = gems.KvlTransform(np.asfortranarray(trf.matrix))
        return (image[cropping].copy(), transform)

    def align_atlas_to_seg(self):
        """
        The initial stage before mesh fitting involves aligning the atlas coordinates to
        the image coordinates by registering an atlas image to the subject's segmentation.
        This step requires that the atlasAlignmentTarget has been properly configured
        during preprocessing.
        """

        # Make sure the subclass has computed the target mask
        mask = self.atlasAlignmentTarget.copy()
        if mask is None:
            sf.system.fatal('All MeshModel subclasses must compute atlasAlignmentTarget during preprocessing!')

        # No need for a high-resolution alignment here
        mask.data = mask.data > 0
        if np.mean(mask.geom.voxsize) < 0.99:
            mask = mask.resize(1, interp='nearest')

        # Crop mask to the label bounding box
        mask = mask.crop_to_bbox(margin=6)

        # Let's smooth the mask a bit (maybe) with one dilation and erosion pass
        if self.atlasTargetSmoothing == 'forward':
            strel = utils.spherical_strel(1)
            mask.data = scipy.ndimage.morphology.binary_dilation(mask.data, strel)
            mask.data = scipy.ndimage.morphology.binary_erosion(mask.data, strel, border_value=1)
        elif self.atlasTargetSmoothing == 'backward':
            strel = utils.spherical_strel(1)
            mask.data = scipy.ndimage.morphology.binary_erosion(mask.data, strel, border_value=1)
            mask.data = scipy.ndimage.morphology.binary_dilation(mask.data, strel)
        elif self.atlasTargetSmoothing is not None:
            sf.system.fatal(f'Unknown atlasTargetSmoothing option `{self.atlasTargetSmoothing}`.')

        # We're going to use mri_robust_register for this registration, so let's ensure the mask
        # value is 255 and we'll write to disk
        mask.data = mask.data.astype('float32') * 255
        targetMaskFile = os.path.join(self.tempDir, 'targetMask.mgz')
        mask.save(targetMaskFile)

        # Write the atlas as well
        alignedAtlasFile = os.path.join(self.tempDir, 'alignedAtlasImage.mgz')
        # ATH skipping this for now... we'll just copy instead
        # self.atlasImage.save(alignedAtlasFile)
        shutil.copyfile(self.atlasDumpFileName, alignedAtlasFile)

        # Run the actual registration and load the result
        utils.run(f'mri_robust_register --mov {alignedAtlasFile} --dst {targetMaskFile} --lta {self.tempDir}/trash.lta --mapmovhdr {alignedAtlasFile} --sat 50 -verbose 0')
        utils.run(f'mri_robust_register --mov {alignedAtlasFile} --dst {targetMaskFile} --lta {self.tempDir}/trash.lta --mapmovhdr {alignedAtlasFile} --affine --sat 50 -verbose 0')
        self.alignedAtlas = sf.load_volume(alignedAtlasFile)

    def prepare_for_seg_fitting(self):
        """
        Prepare the mesh collection, preprocessed image, reduced alphas, and Gaussians parameters.
        """

        # Make sure the subclass has computed the synthed target
        if self.synthImage is None:
            sf.system.fatal('All MeshModel subclasses must compute synthImage during preprocessing!')

        # Crop the synthesized image by the aligned atlas and compute the new mesh alignment
        self.workingImage, self.transform = self.crop_image_by_atlas(self.synthImage)
        self.workingImageShape = self.workingImage.shape[:3]

        # Read in collection, set stiffness, and apply transform
        self.meshCollection = gems.KvlMeshCollection()
        self.meshCollection.read(self.atlasMeshFileName)
        self.meshCollection.transform(self.transform)
        self.meshCollection.k = self.meshStiffness

        # Retrieve the reference mesh, i.e. the mesh representing the average shape
        self.mesh = self.meshCollection.reference_mesh
        self.originalNodePositions = self.mesh.points.copy(order='K')
        self.originalAlphas = self.mesh.alphas.copy(order='K')

        # Compute the cheating Gaussian label groups
        labelGroups = self.get_cheating_label_groups()
        self.sameGaussianParameters = self.label_group_names_to_indices(labelGroups)

        # Compute the reduced alphas - those referring to the super-structures
        self.reducedAlphas, _ = self.reduce_alphas(self.sameGaussianParameters)
        self.mesh.alphas = self.reducedAlphas 
        mask = (self.mesh.rasterize(self.workingImageShape).sum(-1) / 65535) > 0.99
        if self.cheatingAlphaMaskStrel > 0:
            mask = scipy.ndimage.morphology.binary_erosion(mask, utils.spherical_strel(self.cheatingAlphaMaskStrel), border_value=1)
        self.workingImage.data[mask == 0] = 0

        # Get the inital Gaussian parameters
        self.means, self.variances = self.get_cheating_gaussians(self.sameGaussianParameters)

        # Write the inital and cropped/masked images for debugging purposes
        if self.debug:
            self.synthImage.save(os.path.join(self.tempDir, 'synthImage.mgz'))
            self.workingImage.save(os.path.join(self.tempDir, 'synthImageMasked.mgz'))

    def fit_mesh_to_seg(self):
        """
        The second processing step involves deforming the roughly-aligned mesh to the subject segmentation.
        This is the initial 'cheating' step and requires that the synthImage volume has been properly configured
        during preprocessing.
        """

        # Just get the image buffer (array) and convert to a Kvl image object
        imageBuffer = self.workingImage.data.copy(order='K')
        image = gems.KvlImage(gems.utilities.requireNumpyArray(imageBuffer))

        # Use a multi-resolution approach
        for multiResolutionLevel, meshSmoothingSigma in enumerate(self.cheatingMeshSmoothingSigmas):

            # Set mesh alphas
            self.mesh.alphas = self.reducedAlphas

            # It's good to smooth the mesh, otherwise we get weird compressions of the mesh along the boundaries
            if meshSmoothingSigma > 0:
                print(f'Smoothing mesh collection with kernel size {meshSmoothingSigma}')
                self.meshCollection.smooth(meshSmoothingSigma)

            # Note that it uses variances instead of precisions
            calculator = gems.KvlCostAndGradientCalculator(
                typeName='AtlasMeshToIntensityImage',
                images=[image],
                boundaryCondition='Sliding',
                transform=self.transform,
                means=self.means.reshape((-1, 1)),
                variances=self.variances.reshape((-1, 1, 1)),
                mixtureWeights=np.ones(len(self.means), dtype='float32'),
                numberOfGaussiansPerClass=np.ones(len(self.means), dtype='int32'))

            # Step some optimizer stop criteria
            maximalDeformationStopCriterion = 1e-10
            relativeChangeInCostStopCriterion = 1e-10

            # Get optimizer and plug calculator into it
            optimizationParams = {
                'Verbose': False,
                'MaximalDeformationStopCriterion': maximalDeformationStopCriterion,
                'LineSearchMaximalDeformationIntervalStopCriterion': 1e-10,
                'MaximumNumberOfIterations': 1000,
                'BFGS-MaximumMemoryLength': 12
            }
            optimizer = gems.KvlOptimizer(self.optimizerType, self.mesh, calculator, optimizationParams)

            # Run the optimizations
            history = []
            for iteration in range(self.cheatingMaxIterations[multiResolutionLevel]):

                # Step optimizer
                minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer_samseg()

                # Log step information
                iterationInfo = [
                    f'Res: {multiResolutionLevel + 1:03d}',
                    f'Iter: {iteration + 1:03d}',
                    f'MaxDef: {maximalDeformation:.4f}',
                    f'MinLLxP: {minLogLikelihoodTimesPrior:.4f}',
                ]
                print('  '.join(iterationInfo))

                # Track optimization history
                previous = history[-1] if history else np.finfo(np.float32).max
                history.append(minLogLikelihoodTimesPrior)

                # Check for stop criteria
                relativeChange = np.abs((previous - minLogLikelihoodTimesPrior) / minLogLikelihoodTimesPrior)
                if maximalDeformation <= maximalDeformationStopCriterion or relativeChange < relativeChangeInCostStopCriterion:
                    break

        # OK, we're done. Let's modify the mesh atlas in such a way that our computed mesh node positions are
        # assigned to what was originally the mesh warp corresponding to the first training subject
        self.mesh.alphas = self.originalAlphas 
        self.meshCollection.set_positions(self.originalNodePositions, [self.mesh.points])

        # Write the resulting atlas mesh to file in native atlas space.
        # This is nice because all we need to do is to modify imageDump_coregistered
        # with the T1-to-T2 transform to have the warped mesh in T2 space
        inverseTransform = gems.KvlTransform(np.asfortranarray(np.linalg.inv(self.transform.as_numpy_array)))
        self.meshCollection.transform(inverseTransform)
        self.meshCollection.write(self.warpedMeshFileName)

    def prepare_for_image_fitting(self, compute_hyps=True):
        """
        Prepare the mesh collection, preprocessed image, reduced alphas, and estimated hyperparameters.
        """

        # Make sure the subclass has computed the target image
        if self.processedImage is None:
            sf.system.fatal('All MeshModel subclasses must compute processedImage during preprocessing!')

        # Crop the image by the aligned atlas and compute the new mesh alignment
        self.workingImage, self.transform = self.crop_image_by_atlas(self.processedImage)

        # ATH for now let's squeeze the data, but will need to adapt something better
        # for multi-image cases down the road
        self.workingImage.data = np.asfortranarray(self.workingImage.data.squeeze())
        self.workingImageShape = self.workingImage.data.shape[:3]

        # Read the atlas mesh from file, and apply the previously determined transform to the location of its nodes
        # ATH does this have to be re-read?
        self.meshCollection = gems.KvlMeshCollection()
        self.meshCollection.read(self.warpedMeshFileName)
        self.meshCollection.transform(self.transform)
        self.meshCollection.k = self.meshStiffness

        # Retrieve the correct mesh to use from the meshCollection
        self.mesh = self.meshCollection.get_mesh(0)

        # We're not interested in image areas that fall outside our cuboid ROI where our atlas is defined. Therefore,
        # generate a mask of what's inside the ROI. Also, by convention we're skipping all voxels with zero intensity.
        mask = (self.mesh.rasterize(self.workingImageShape).sum(-1) / 65535) > 0.99
        if self.alphaMaskStrel > 0:
            mask = scipy.ndimage.morphology.binary_erosion(mask, utils.spherical_strel(self.alphaMaskStrel), border_value=1)
        mask = np.asfortranarray(mask & (self.workingImage.data > 0))

        # Apply the mask to the image we're analyzing by setting the intensity of all voxels not belonging
        # to the brain mask to zero. This will automatically discard those voxels in subsequent C++ routines, as
        # voxels with intensity zero are simply skipped in the computations.
        self.workingMask = self.workingImage.new(mask)
        self.workingImage.data[mask == 0] = 0
        # Let's do this to make results more similar to the matlab version
        self.maskIndices = np.unravel_index(np.where(mask.flatten(order='F')), self.workingImageShape, order='F')

        # Write the initial and cropped/masked images for debugging purposes
        if self.debug:
            self.processedImage.save(os.path.join(self.tempDir, 'processedImage.mgz'))
            self.workingImage.save(os.path.join(self.tempDir, 'processedImageMasked.mgz'))
            self.workingMask.save(os.path.join(self.tempDir, 'processedImageMask.mgz'))

        # Compute the Gaussian label groups
        labelGroups = self.get_label_groups()
        self.sameGaussianParameters = self.label_group_names_to_indices(labelGroups)

        # Compute the reduced alphas
        self.reducedAlphas, self.reducingLookupTable = self.reduce_alphas(self.sameGaussianParameters)
        self.mesh.alphas = self.reducedAlphas

        if compute_hyps:
            # Compute the hyperparameters
            self.meanHyper, self.nHyper = self.get_gaussian_hyps(self.sameGaussianParameters, self.mesh)

        # Init empty means and variances
        self.means = None
        self.variances = None

    def fit_mesh_to_image(self):
        """
        Fit mesh to the image data.
        """

        # Just get the original image buffer (array) and convert to a Kvl image object
        imageBuffer = self.workingImage.data.copy(order='K')
        image = gems.KvlImage(gems.utilities.requireNumpyArray(imageBuffer))

        # Useful to have cached
        numMaskIndices = self.maskIndices[0].shape[-1]
        numberOfClasses = len(self.sameGaussianParameters)

        # Multi-resolution loop
        numberOfMultiResolutionLevels = len(self.meshSmoothingSigmas)
        for multiResolutionLevel in range(numberOfMultiResolutionLevels):

            if self.isLong:
                self.mesh = self.meshCollection.get_mesh(0)

            # Special case when we want to recompute reduced alphas for a second-component
            # Note: how should we deal with more than one component during longitudinal global iterations?
            if self.useTwoComponents and multiResolutionLevel == 1:
                # Get second component label groups
                labelGroups = self.get_second_label_groups()
                self.sameGaussianParameters = self.label_group_names_to_indices(labelGroups)
                numberOfClasses = len(self.sameGaussianParameters)
                self.reducedAlphas, self.reducingLookupTable = self.reduce_alphas(self.sameGaussianParameters)
                self.mesh.alphas = self.reducedAlphas
                # Compute new Gaussian hyperparameters
                self.meanHyper, self.nHyper = self.get_second_gaussian_hyps(self.sameGaussianParameters, self.meanHyper, self.nHyper)
                # Reset means and variances to be computed
                self.means = None
                self.variances = None

            # Set the mesh alphas back
            # ATH is this necessary though?
            self.mesh.alphas = self.reducedAlphas

            # Smooth the mesh using a Gaussian kernel
            meshSmoothingSigma = self.meshSmoothingSigmas[multiResolutionLevel]
            if meshSmoothingSigma > 0:
                print(f'Smoothing mesh collection with kernel size {meshSmoothingSigma:.4f}')
                self.meshCollection.smooth(meshSmoothingSigma)

            # Smooth the image using a Gaussian kernel
            imageSigma = self.imageSmoothingSigmas[multiResolutionLevel]
            if imageSigma > 0:
                raise NotImplementedError('Image smoothing not implemented yet!')

            # ATH this is in case the above smoothing only sets the buffer, but this should be removed
            # really since it's not necessary if things are correctly implemented
            image = gems.KvlImage(gems.utilities.requireNumpyArray(imageBuffer))
            
            # Now with this smoothed atlas, we're ready for the real work. There are essentially two sets of parameters
            # to estimate in our generative model: (1) the mesh node locations (parameters of the prior), and (2) the
            # means and variances of the Gaussian intensity models (parameters of the
            # likelihood function, which is really a hugely simplistic model of the MR imaging process). Let's optimize
            # these two sets alternately until convergence.

            positionUpdatingMaximumNumberOfIterations = 30
            maximumNumberOfIterations = self.maxIterations[multiResolutionLevel]
            
            historyOfCost = []
            for iterationNumber in range(maximumNumberOfIterations):
                print(f'Iteration {iterationNumber + 1} of {maximumNumberOfIterations}')

                # Part I: estimate Gaussian mean and variances using EM

                # Get the priors as dictated by the current mesh position as well as the image intensities
                data = imageBuffer[self.maskIndices]

                # Avoid spike in memory during the posterior computation
                priors = np.zeros((numMaskIndices, numberOfClasses), dtype='uint16')
                for l in range(numberOfClasses):
                    priors[:, l] = self.mesh.rasterize(self.workingImageShape, l)[self.maskIndices]
                posteriors = priors / 65535

                # Start EM iterations. Initialize the parameters if this is the first time ever you run this
                if (self.means is None) or (self.variances is None):

                    self.means = np.zeros(numberOfClasses)
                    self.variances = np.zeros(numberOfClasses)

                    thresh = 1e-2
                    for classNumber in range(numberOfClasses):
                        posterior = posteriors[:, classNumber]
                        if np.sum(posterior) > thresh:
                            mu = (self.meanHyper[classNumber] * self.nHyper[classNumber] + data @ posterior) / (self.nHyper[classNumber] + np.sum(posterior) + thresh)
                            variance = (((data - mu) ** 2) @ posterior + self.nHyper[classNumber] * (mu - self.meanHyper[classNumber]) ** 2) / (np.sum(posterior) + thresh)
                            self.means[classNumber] = mu
                            self.variances[classNumber] = variance + thresh
                        else:
                            self.means[classNumber] = self.meanHyper[classNumber]
                            self.variances[classNumber] = 100

                    # Prevents NaNs during the optimization
                    self.variances[self.variances == 0] = 100

                stopCriterionEM = 1e-5
                historyOfEMCost = []
                for EMIterationNumber in range(100):

                    # E-step: compute the posteriors based on the current parameters

                    minLogLikelihood = 0
                    for classNumber in range(numberOfClasses):
                        mu = self.means[classNumber]
                        variance = self.variances[classNumber]
                        prior = priors[:, classNumber] / 65535
                        posteriors[:, classNumber] = (np.exp(-(data - mu) ** 2 / 2 / variance) * prior) / np.sqrt(2 * np.pi * variance)
                        minLogLikelihood = minLogLikelihood + 0.5 * np.log(2 * np.pi * variance) - 0.5 * np.log(self.nHyper[classNumber]) + \
                                                              0.5 * self.nHyper[classNumber] / variance * (mu - self.meanHyper[classNumber]) ** 2

                    normalizer = np.sum(posteriors, -1) + np.finfo(np.float32).eps
                    posteriors /= normalizer[..., np.newaxis]
                    minLogLikelihood = minLogLikelihood - np.sum(np.log(normalizer))  # This is what we're optimizing with EM
                    if np.isnan(minLogLikelihood):
                        sf.system.fatal('minLogLikelihood is NaN')

                    # Log some iteration information
                    iterationInfo = [
                        f'Res: {multiResolutionLevel + 1:03d}',
                        f'Iter: {iterationNumber + 1:03d} | {EMIterationNumber + 1:03d}',
                        f'MinLL: {minLogLikelihood:.4f}',
                    ]
                    print('  '.join(iterationInfo))

                    # Track EM history
                    previous = historyOfEMCost[-1] if historyOfEMCost else np.finfo(np.float32).max
                    historyOfEMCost.append(minLogLikelihood)

                    # Check for convergence
                    relativeChangeCost = (previous - minLogLikelihood) / minLogLikelihood
                    if relativeChangeCost < stopCriterionEM:
                        print('EM converged!')
                        break

                    # M-step: derive parameters from the posteriors

                    # Update parameters of Gaussian mixture model
                    thresh = 1e-2
                    for classNumber in range(numberOfClasses):
                        posterior = posteriors[:, classNumber]
                        if np.sum(posterior) > thresh:
                            mu = (self.meanHyper[classNumber] * self.nHyper[classNumber] + data @ posterior) / (self.nHyper[classNumber] + np.sum(posterior) + thresh)
                            variance = (((data - mu) ** 2) @ posterior + self.nHyper[classNumber] * (mu - self.meanHyper[classNumber]) ** 2) / (np.sum(posterior) + thresh)
                            self.means[classNumber] = mu
                            self.variances[classNumber] = variance + thresh
                        else:
                            self.means[classNumber] = self.meanHyper[classNumber]
                            self.variances[classNumber] = 100

                    # Prevents NaNs during the optimization
                    self.variances[self.variances == 0] = 100

                # Part II: update the position of the mesh nodes for the current set of Gaussian parameters

                if self.isLong:
                    self.mesh = self.meshCollection.get_mesh(0)

                # Keep track if the mesh has moved or not
                haveMoved = False

                # Note that it uses variances instead of precisions
                calculator = gems.KvlCostAndGradientCalculator(
                    typeName='AtlasMeshToIntensityImage',
                    images=[image],
                    boundaryCondition='Sliding',
                    transform=self.transform,
                    means=self.means[..., np.newaxis],
                    variances=self.variances[..., np.newaxis, np.newaxis],
                    mixtureWeights=np.ones(len(self.means), dtype='float32'),
                    numberOfGaussiansPerClass=np.ones(len(self.means), dtype='int32'))

                # Get optimizer and plug calculator into it
                maximalDeformationStopCriterion = 1e-10
                optimizationParameters = {
                    'Verbose': 0,
                    'MaximalDeformationStopCriterion': maximalDeformationStopCriterion,
                    'LineSearchMaximalDeformationIntervalStopCriterion': 1e-10,
                    'MaximumNumberOfIterations': 1000,
                    'BFGS-MaximumMemoryLength': 12
                }
                optimizer = gems.KvlOptimizer(self.optimizerType, self.mesh, calculator, optimizationParameters)

                for positionUpdatingIterationNumber in range(positionUpdatingMaximumNumberOfIterations):

                    # Calculate a good step
                    minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer_samseg()

                    # Log optimization information
                    iterationInfo = [
                        f'Res: {multiResolutionLevel + 1:03d}',
                        f'Iter: {iterationNumber + 1:03d} | {positionUpdatingIterationNumber + 1:03d}',
                        f'MaxDef: {maximalDeformation:.4f}',
                        f'MinLLxP: {minLogLikelihoodTimesPrior:.4f}',
                    ]
                    print('  '.join(iterationInfo))

                    if np.isnan(minLogLikelihoodTimesPrior):
                        print('error: minLogLikelihoodTimesPrior is NaN')

                    if maximalDeformation > 0:
                        haveMoved = True

                    # Check if we need to stop
                    if maximalDeformation <= maximalDeformationStopCriterion:
                        print('maximalDeformation is too small - stopping')
                        break

                # Keep track of the cost function we're optimizing
                previous = historyOfCost[-1] if historyOfCost else np.finfo(np.float32).max
                historyOfCost.append(minLogLikelihoodTimesPrior)
                
                # Determine if we should stop the overall iterations over the two set of parameters
                if not haveMoved or (((previous - minLogLikelihoodTimesPrior) / minLogLikelihoodTimesPrior) < 1e-6):
                    break

    def extract_segmentation(self):
        """
        Extract discrete labels and volumes from the fit mesh.
        """

        # First, undo the collapsing of several structures into super-structures
        self.mesh.alphas = self.originalAlphas
        numberOfClasses = self.originalAlphas.shape[-1]
        numMaskIndices = self.maskIndices[0].shape[-1]

        # Compute normalized posteriors
        imgdata = self.workingImage.data[self.maskIndices]
        posteriors = np.zeros((numMaskIndices, numberOfClasses), dtype='float32')
        for classNumber in range(numberOfClasses):
            prior = self.mesh.rasterize(self.workingImageShape, classNumber)
            mu = self.means[self.reducingLookupTable[classNumber]]
            variance = self.variances[self.reducingLookupTable[classNumber]]
            posteriors[:, classNumber] = (np.exp(-(imgdata - mu) ** 2 / 2 / variance) * (prior[self.maskIndices] / 65535)) / np.sqrt(2 * np.pi * variance)
        normalizer = np.sum(posteriors, -1) + np.finfo(np.float32).eps
        posteriors /= normalizer[..., np.newaxis]
        posteriors = np.round(posteriors * 65535).astype('uint16')

        if self.debug:
            # Write the resulting atlas mesh to file for future reference
            self.meshCollection.write(os.path.join(self.tempDir, 'finalWarpedMesh.txt'))
            # Also write the warped mesh in atlas space
            inverseTransform = gems.KvlTransform(np.asfortranarray(np.linalg.inv(self.transform.as_numpy_array)))
            self.meshCollection.transform(inverseTransform)
            self.meshCollection.write(os.path.join(self.tempDir, 'finalWarpedMeshNoAffine.txt'))

        # Here we do a memory efficient computation of discrete labels and volumes
        self.volumes = {}
        inds = np.zeros(self.workingImageShape, dtype='int32')
        for i in range(numberOfClasses):

            if i == 0:
                sillyAlphas = np.zeros((len(self.originalAlphas), 2), dtype='float32')
                sillyAlphas[:, 0] = self.originalAlphas[:, 0]
                sillyAlphas[:, 1] = 1 - sillyAlphas[:, 0]
                self.mesh.alphas = sillyAlphas
                post = self.mesh.rasterize(self.workingImageShape)[..., 0]
                self.mesh.alphas = self.originalAlphas
            else:
                post = self.mesh.rasterize(self.workingImageShape, i)

            post[self.maskIndices] = posteriors[:, i]

            if i == 0:
                max_post = post
            else:
                M = post > max_post
                inds[M] = i
                max_post[M] = post[M]
            
            # Compute volume
            self.volumes[self.names[i]] = (self.resolution ** 3) * (post.sum() / 65535)

        # Compute all discrete labels and mask
        self.discreteLabels = self.workingImage.new(self.FreeSurferLabels[inds])
        self.discreteLabels[self.workingMask == 0] = 0
        lut_filename = os.path.join(os.environ.get('FREESURFER_HOME'), 'FreeSurferColorLUT.txt')
        self.discreteLabels.labels = sf.load_label_lookup(lut_filename)

        if self.debug:
            self.discreteLabels.save(os.path.join(self.tempDir, 'discreteLabelsAll.mgz'))

    def write_volumes(self, filename, volumes=None):
        """
        Write the cached volume dictionary to a text file.
        """
        if volumes is None:
            volumes = self.volumes 
        with open(filename, 'w') as file:
            for name, volume in volumes.items():
                file.write(f'{name} {volume:.6f}\n')

    def postprocess_segmentation(self):
        """
        This function should perform any necessary modifications to and write the discreteLabels segmentation and labelVolumes.
        """
        raise NotImplementedError('A MeshModel subclass must implement the postprocess_segmentation() function!')

    def get_cheating_label_groups(self):
        """
        This function should return a group (list of lists) of label names that determine the class
        reductions for the initial segmentation-fitting stage.
        """
        raise NotImplementedError('A MeshModel subclass must implement the get_cheating_label_groups() function!')

    def get_cheating_gaussians(self, sameGaussianParameters):
        """
        This function should return a tuple of (means, variances) for the initial segmentation-fitting stage.
        """
        raise NotImplementedError('A MeshModel subclass must implement the get_cheating_gaussians() function!')

    def get_label_groups(self):
        """
        This function should return a group (list of lists) of label names that determine the class
        reductions for the primary image-fitting stage.
        """
        raise NotImplementedError('A MeshModel subclass must implement the get_label_groups() function!')

    def get_gaussian_hyps(self, sameGaussianParameters, mesh):
        """
        This function should return a tuple of (meanHyps, nHyps) for Gaussian parameter estimation.
        """
        raise NotImplementedError('A MeshModel subclass must implement the get_gaussian_hyps() function!')

    def get_second_label_groups(self):
        """
        This optional function should return a group (list of lists) of label names that determine the class
        reductions for the second-component of the primary image-fitting stage.
        """
        raise NotImplementedError('A two-component MeshModel must implement the get_second_label_groups() function!')

    def get_second_gaussian_hyps(self, sameGaussianParameters, meanHyper, nHyper):
        """
        This optional function should return a tuple of (meanHyps, nHyps) for Gaussian parameter estimation
        in the second-component of the primary image-fitting stage.
        """
        raise NotImplementedError('A two-component MeshModel must implement the get_second_gaussian_hyps() function!')
