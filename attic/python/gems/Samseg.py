import os
import sys
import logging
import pickle
import scipy.io
import surfa as sf
from scipy.ndimage.morphology import binary_dilation as dilation

import gems
from gems.utilities import Specification
from gems.BiasField import BiasField
from gems.ProbabilisticAtlas import ProbabilisticAtlas
from gems.GMM import GMM
from gems.Affine import Affine
from gems.SamsegUtility import *
from gems.merge_alphas import kvlMergeAlphas, kvlGetMergingFractionsTable


eps = np.finfo(float).eps


class Samseg:
    def __init__(self,
        imageFileNames,
        atlasDir,
        savePath,
        userModelSpecifications={},
        userOptimizationOptions={},
        imageToImageTransformMatrix=None,
        visualizer=None,
        saveHistory=None,
        savePosteriors=False,
        saveWarp=None,
        saveMesh=None,
        threshold=None,
        thresholdSearchString=None,
        targetIntensity=None,
        targetSearchStrings=None,
        modeNames=None,
        pallidumAsWM=True,
        saveModelProbabilities=False,
        gmmFileName=None,
        ignoreUnknownPriors=False,
        dissectionPhoto=None,
        nthreads=1,
        ):

        # Store input parameters as class variables
        self.imageFileNames = imageFileNames
        self.originalImageFileNames = imageFileNames  # Keep a copy since photo version modifies self.imageFileNames
        self.photo_mask = None  # Useful when working with photos
        self.savePath = savePath
        self.atlasDir = atlasDir
        self.threshold = threshold
        self.thresholdSearchString = thresholdSearchString
        self.targetIntensity = targetIntensity
        self.targetSearchStrings = targetSearchStrings

        # Use defualt if mode names aren't provided
        if not modeNames:
            modeNames = ['mode%02d' % (n + 1) for n in range(len(imageFileNames))]
        elif len(modeNames) != len(imageFileNames):
            raise ValueError('number of mode names does not match number of input images')
        self.modeNames = modeNames

        # Eugenio: there's a bug in ITK that will cause kvlImage to fail if it contqins the string "recon" ...
        # If this problem is not exclusive to the photo mode (RGB), we should move this chunk of code outside the if
        # While at it, we also create a grayscale version and a version with a bit of noise around the cerebrum (so the
        # background class doesn't go bananas)
        # Dissection photos are converted to grayscale
        if dissectionPhoto:

            if len(self.imageFileNames) > 1:
                raise ValueError('In photo mode, you cannot provide more than one input image volume')

            input_vol = sf.load_volume(self.imageFileNames[0])
            input_vol.save(self.savePath + '/original_input.mgz')
            if input_vol.nframes > 1:
                input_vol = input_vol.mean(frames=True)
            input_vol.save(self.savePath + '/grayscale_input.mgz')

            # We also a small band of noise around the mask; otherwise the background/skull/etc may fit the cortex
            self.photo_mask = input_vol.data > 0
            mask_dilated = dilation(self.photo_mask, iterations=5)
            ring = (mask_dilated == True) & (self.photo_mask == False)
            max_noise = input_vol.max() / 50.0
            rng = np.random.default_rng(2021)
            input_vol.data[ring] = max_noise * rng.random(input_vol.data[ring].shape[0])
            self.imageFileNames = []
            self.imageFileNames.append(self.savePath + '/grayscale_input_with_ring.mgz')
            input_vol.save(self.imageFileNames[0])


        # Initialize some objects
        self.affine = Affine( imageFileName=self.imageFileNames[0],
                              meshCollectionFileName=os.path.join(self.atlasDir, 'atlasForAffineRegistration.txt.gz'),
                              templateFileName=os.path.join(self.atlasDir, 'template.nii' ) )
        self.probabilisticAtlas = ProbabilisticAtlas()

        # Get full model specifications and optimization options (using default unless overridden by user)
        # Note that, when processing photos, we point to a different GMM file by default!
        self.optimizationOptions = getOptimizationOptions(atlasDir, userOptimizationOptions)
        if dissectionPhoto and (gmmFileName is None):
            if dissectionPhoto == 'left':
                gmmFileName = self.atlasDir + '/photo.lh.sharedGMMParameters.txt'
            elif dissectionPhoto == 'right':
                gmmFileName = self.atlasDir + '/photo.rh.sharedGMMParameters.txt'
            elif dissectionPhoto == 'both':
                gmmFileName = self.atlasDir + '/photo.both.sharedGMMParameters.txt'
            else:
                sf.system.fatal('dissection photo mode must be left, right, or both')
        self.modelSpecifications = getModelSpecifications(
            atlasDir,
            userModelSpecifications,
            pallidumAsWM=pallidumAsWM,
            gmmFileName=gmmFileName
        )

        # if dissectionPhoto is not None:
        #     self.modelSpecifications['K'] = 0.01
        # TODO: switch this on, and maybe make K dependent on the pixel size?

        # Set image-to-image matrix if provided
        self.imageToImageTransformMatrix = imageToImageTransformMatrix

        # Print specifications
        print('##----------------------------------------------')
        print('              Samsegment Options')
        print('##----------------------------------------------')
        print('output directory:', savePath)
        print('input images:', ', '.join([imageFileName for imageFileName in imageFileNames]))
        print('modelSpecifications:', self.modelSpecifications)
        print('optimizationOptions:', self.optimizationOptions)

        # Convert modelSpecifications from dictionary into something more convenient to access
        self.modelSpecifications = Specification(self.modelSpecifications)

        # Setup a null visualizer if necessary
        if visualizer is None:
            self.visualizer = initVisualizer(False, False)
        else:
            self.visualizer = visualizer
        
        # Configure posterior saving option
        if savePosteriors is True:
            self.savePosteriors = self.modelSpecifications.names
        elif isinstance(savePosteriors, list):
            self.savePosteriors = savePosteriors
        else:
            self.savePosteriors = None

        self.saveModelProbabilities = saveModelProbabilities
        self.saveHistory = saveHistory
        self.saveWarp = saveWarp
        self.saveMesh = saveMesh
        self.ignoreUnknownPriors = ignoreUnknownPriors
        self.dissectionPhoto = dissectionPhoto

        # Make sure we can write in the target/results directory
        os.makedirs(savePath, exist_ok=True)

        # Class variables that will be used later
        self.transform = None
        self.biasField = None
        self.gmm = None
        self.imageBuffers = None
        self.mask = None
        self.classFractions = None
        self.cropping = None
        self.voxelSpacing = None
        self.optimizationSummary = None
        self.optimizationHistory = None
        self.deformation = None
        self.deformationAtlasFileName = None
        self.nthreads = nthreads

    def validateTransform(self, trf):
        # =======================================================================================
        #
        # Internal utility to ensure that a provided affine transform file matches the template and
        # input volume geometries. If the transform is determined to be [image->template], it
        # will be inverted for convenience sake.
        #
        # =======================================================================================

        # Load src (template) and trg (input) geometries (TODO these should really be cached)
        src = sf.load_volume(self.affine.templateFileName).geom
        trg = sf.load_volume(self.imageFileNames[0]).geom

        # Just assume things are okay if no geometries are provided
        if trf.source is None and trf.target is None:
            return trf

        equal = lambda a, b: sf.transform.image_geometry_equal(a, b, tol=1e-2)

        # Make sure at least the source or target geometries match
        if (trf.source is None or equal(trf.source, src)) and (trf.target is None or equal(trf.target, trg)):
            return trf

        # The only remaining possibility is that the transform is inverted (input->template)
        if (trf.source is None or equal(trf.source, trg)) and (trf.target is None or equal(trf.target, src)):
            return trf.inv()

        sf.system.fatal('provided transform does not match input or template geometries')

    def segment(self, costfile=None, timer=None, reg_only=False, transformFile=None, initTransformFile=None):
        # =======================================================================================
        #
        # Main function that runs the whole segmentation pipeline
        #
        # =======================================================================================

        # Initialization transform for registration
        initTransform = None
        if initTransformFile:
            trg = self.validateTransform(sf.load_affine(initTransformFile))
            initTransform = convertRASTransformToLPS(trg.convert(space='world').matrix)

        # Affine transform used to skip registration
        worldToWorldTransformMatrix = None
        if transformFile:
            if transformFile.endswith('.mat'):
                worldToWorldTransformMatrix = scipy.io.loadmat(transformFile).get('worldToWorldTransformMatrix')
            else:
                trf = self.validateTransform(sf.load_affine(transformFile))
                worldToWorldTransformMatrix = convertRASTransformToLPS(trf.convert(space='world').matrix)

        # Register to template, either with SAMSEG code, or externally with FreeSurfer tools (for photos)
        if self.imageToImageTransformMatrix is None:

            if self.dissectionPhoto is not None:
                reference = self.savePath + '/grayscale_input.mgz'
                if self.dissectionPhoto=='left':
                    moving = self.atlasDir + '/exvivo.template.lh.suptent.nii'
                elif self.dissectionPhoto=='right':
                    moving = self.atlasDir + '/exvivo.template.rh.suptent.nii'
                elif self.dissectionPhoto=='both':
                    moving = self.atlasDir + '/exvivo.template.suptent.nii'
                else:
                    sf.system.fatal('dissection photo mode must be left, right, or both')
                transformFile = self.savePath  + '/atlas2image.lta'
                cmd = 'mri_coreg  --seed 2021 --mov ' + moving + ' --ref ' + reference + ' --reg ' + transformFile + \
                      ' --dof 12 --threads ' + str(self.nthreads)
                os.system(cmd)
                trf_val = self.validateTransform(sf.load_affine(transformFile))
                worldToWorldTransformMatrix = convertRASTransformToLPS(trf_val.convert(space='world').matrix)

            self.register(
                costfile=costfile,
                timer=timer,
                reg_only=reg_only,
                worldToWorldTransformMatrix=worldToWorldTransformMatrix,
                initTransform=initTransform
            )

        self.preProcess()
        self.fitModel()
        return self.postProcess()

    def register(self, costfile=None, timer=None, reg_only=False, worldToWorldTransformMatrix=None, initTransform=None):
        # =======================================================================================
        #
        # Perform affine registration if needed
        #
        # =======================================================================================

        # Perform registration on first input image
        self.imageToImageTransformMatrix, optimizationSummary = self.affine.registerAtlas(
            savePath=self.savePath,
            visualizer=self.visualizer,
            worldToWorldTransformMatrix=worldToWorldTransformMatrix,
            initTransform=initTransform
        )

        # Save a summary of the optimization process
        if optimizationSummary and costfile is not None:
            with open(costfile, "a") as file:
                file.write('templateRegistration %d %f\n' % (optimizationSummary['numberOfIterations'], optimizationSummary['cost']))

        # Mark registration time
        if timer is not None:
            timer.mark('atlas registration complete')

        # Exit if specified
        if reg_only:
            print('registration-only requested, so quiting now')
            sys.exit()

    def getMesh(self, *args, **kwargs):
        """
        Load the atlas mesh and perform optional smoothing of the cortex and WM priors.
        """
        if self.modelSpecifications.whiteMatterAndCortexSmoothingSigma > 0:
            competingNames = [['Left-Cerebral-White-Matter',  'Left-Cerebral-Cortex'],
                              ['Right-Cerebral-White-Matter', 'Right-Cerebral-Cortex']]
            competingStructures = [[self.modelSpecifications.names.index(n) for n in names] for names in competingNames]
            kwargs['competingStructures'] = competingStructures
            kwargs['smoothingSigma'] = self.modelSpecifications.whiteMatterAndCortexSmoothingSigma
        return self.probabilisticAtlas.getMesh(*args, **kwargs)

    def preProcess(self):
        # =======================================================================================
        #
        # Preprocessing (reading and masking of data)
        #
        # =======================================================================================

        # Read the image data from disk and crop

        # Historically, the template was resaved with a transformed header and the image-to-image transform
        # was later extracted by comparing the input image and coregistered template transforms, but shear
        # cannot be saved through an ITK image, so a better method is to pass the image-to-image transform matrix
        # directly to samseg. If the SAMSEG_LEGACY_REGISTRATION env var is set, this old method is enabled
        if os.environ.get('SAMSEG_LEGACY_REGISTRATION') is not None:
            print('INFO: using legacy (broken) registration option')
            transformedTemplateFileName = os.path.join(self.savePath, 'template_coregistered_legacy.nii')
            self.imageBuffers, self.transform, self.voxelSpacing, self.cropping = readCroppedImagesLegacy(self.imageFileNames, transformedTemplateFileName)
        else:
            self.imageBuffers, self.transform, self.voxelSpacing, self.cropping = readCroppedImages(
                self.imageFileNames,
                os.path.join(self.atlasDir, 'template.nii'),
                self.imageToImageTransformMatrix
            )

        # Background masking: simply setting intensity values outside of a very rough brain mask to zero
        # ensures that they'll be skipped in all subsequent computations
        self.imageBuffers, self.mask = maskOutBackground(
            self.imageBuffers,
            self.modelSpecifications.atlasFileName,
            self.transform,
            self.modelSpecifications.maskingProbabilityThreshold,
            self.modelSpecifications.maskingDistance,
            self.probabilisticAtlas,
            self.voxelSpacing
        )

        # Let's prepare for the bias field correction that is part of the imaging model. It assumes
        # an additive effect, whereas the MR physics indicate it's a multiplicative one - so we log
        # transform the data first.
        self.imageBuffers = logTransform(self.imageBuffers, self.mask)

        # Visualize some stuff
        if hasattr(self.visualizer, 'show_flag'):
            self.visualizer.show(
                mesh=self.getMesh(self.modelSpecifications.atlasFileName, self.transform),
                shape=self.imageBuffers.shape,
                window_id='samsegment mesh', title='Mesh',
                names=self.modelSpecifications.names, legend_width=350)
            self.visualizer.show(images=self.imageBuffers, window_id='samsegment images',
                                 title='Samsegment Masked and Log-Transformed Contrasts')

    def fitModel(self):
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
        posteriors, biasFields, nodePositions, _, _ = self.computeFinalSegmentation()

        # Write out segmentation and bias field corrected volumes
        volumesInCubicMm = self.writeResults(biasFields, posteriors)

        # Save the template warp
        if self.saveWarp:
            print('Saving the template warp')
            self.saveWarpField(os.path.join(self.savePath, 'template.m3z'))

        # Save the final mesh collection
        if self.saveMesh:
            print('Saving the final mesh in template space')
            deformedAtlasFileName = os.path.join(self.savePath, 'mesh.txt')
            self.probabilisticAtlas.saveDeformedAtlas(self.modelSpecifications.atlasFileName, deformedAtlasFileName, nodePositions)
            
        # Save the Gaussian priors, (normalized) likelihoods and posteriors
        if self.saveModelProbabilities:
            self.saveGaussianProbabilities( os.path.join(self.savePath, 'probabilities') )
 
        # Save the class means and standard deviations
        self.saveGaussianStats()

        # Save the history of the parameter estimation process
        if self.saveHistory:
            history = {'input': {
                'imageFileNames': self.imageFileNames,
                'imageToImageTransformMatrix': self.imageToImageTransformMatrix,
                'modelSpecifications': self.modelSpecifications,
                'optimizationOptions': self.optimizationOptions,
                'savePath': self.savePath
            }, 'imageBuffers': self.imageBuffers, 'mask': self.mask,
                'cropping': self.cropping,
                'transform': self.transform.as_numpy_array,
                'historyWithinEachMultiResolutionLevel': self.optimizationHistory,
                "labels": self.modelSpecifications.FreeSurferLabels, "names": self.modelSpecifications.names,
                "volumesInCubicMm": volumesInCubicMm, "optimizationSummary": self.optimizationSummary}
            with open(os.path.join(self.savePath, 'history.p'), 'wb') as file:
                pickle.dump(history, file, protocol=pickle.HIGHEST_PROTOCOL)
                                

        return self.modelSpecifications.FreeSurferLabels, self.modelSpecifications.names, volumesInCubicMm, self.optimizationSummary

    def writeImage(self, data, path, saveLabels=False):
        # Read source geometry
        target = sf.load_volume(self.imageFileNames[0])

        # Account for multi-frame volumes
        frames = data.shape[-1] if data.ndim == 4 else 1

        # Uncrop image
        uncropped = target.zeros(frames=frames, dtype=data.dtype, order='F')
        uncropped[self.cropping] = data
        if saveLabels:
            uncropped.labels = sf.load_label_lookup(os.path.join(self.atlasDir, 'modifiedFreeSurferColorLUT.txt'))
        uncropped.save(path)

    def writeResults(self, biasFields, posteriors):

        # Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention
        names = self.modelSpecifications.names
        if self.threshold is not None:
            # Figure out the structure number of the special snowflake structure
            for structureNumber, name in enumerate(names):
                if self.thresholdSearchString in name:
                    thresholdStructureNumber = structureNumber
                    break

            # Threshold
            print('thresholding posterior of ', names[thresholdStructureNumber], 'with threshold:', self.threshold)
            tmp = posteriors[:, thresholdStructureNumber].copy()
            posteriors[:, thresholdStructureNumber] = posteriors[:, thresholdStructureNumber] > self.threshold

            # Majority voting
            structureNumbers = np.array(np.argmax(posteriors, 1), dtype=np.uint32)

            # Undo thresholding in posteriors
            posteriors[:, thresholdStructureNumber] = tmp

        else:
            # Majority voting
            structureNumbers = np.array(np.argmax(posteriors, 1), dtype=np.uint32)

        segmentation = np.zeros(self.imageBuffers.shape[:3], dtype=np.int32)
        fslabels = np.array(self.modelSpecifications.FreeSurferLabels, dtype=np.int32)
        segmentation[self.mask] = fslabels[structureNumbers]

        #
        self.scalingFactors = scaleBiasFields(biasFields, self.imageBuffers, self.mask, posteriors, self.targetIntensity, self.targetSearchStrings, names)

        # Get corrected intensities and bias field images in the non-log transformed domain
        expImageBuffers, expBiasFields = undoLogTransformAndBiasField(self.imageBuffers, biasFields, self.mask)

        # Write out various images - segmentation first
        self.writeImage(segmentation, os.path.join(self.savePath, 'seg.mgz'), saveLabels=True)

        # Bias corrected images: depends on whether we're dealing with MRIs or 3D photo reconstructions
        if self.dissectionPhoto is None:  # MRI
            for contrastNumber, imageFileName in enumerate(self.imageFileNames):
                # Contrast-specific filename prefix
                contastPrefix = os.path.join(self.savePath, self.modeNames[contrastNumber])

                # Write bias field and bias-corrected image
                self.writeImage(expBiasFields[..., contrastNumber],   contastPrefix + '_bias_field.mgz')
                self.writeImage(expImageBuffers[..., contrastNumber], contastPrefix + '_bias_corrected.mgz')

                # Save a note indicating the scaling factor
                with open(contastPrefix + '_scaling.txt', 'w') as fid:
                    print(self.scalingFactors[contrastNumber], file=fid)

        else: # photos
            self.writeImage(expBiasFields[..., 0], self.savePath + '/illlumination_field.mgz')
            original_vol = sf.load_volume(self.originalImageFileNames[0])
            bias_native = sf.load_volume(self.savePath + '/illlumination_field.mgz')
            original_vol = original_vol / (1e-6 + bias_native)
            original_vol.save(self.savePath + '/illlumination_corrected.mgz')

        if self.savePosteriors:
            posteriorPath = os.path.join(self.savePath, 'posteriors')
            os.makedirs(posteriorPath, exist_ok=True)
            for structureNumber, name in enumerate(names):
                # Write the posteriors to seperate volume files
                for searchString in self.savePosteriors:
                    if searchString in name:
                        posteriorVol = np.zeros(self.imageBuffers.shape[:3], dtype=np.float32)
                        posteriorVol[self.mask] = posteriors[:, structureNumber]
                        self.writeImage(posteriorVol, os.path.join(posteriorPath, name + '.mgz'))

        # Compute volumes in mm^3
        # TODO: cache the source geometry in __init__, as this is also loaded by writeImage
        exampleImage = gems.KvlImage(self.imageFileNames[0])
        volumeOfOneVoxel = np.abs(np.linalg.det(exampleImage.transform_matrix.as_numpy_array[:3, :3]))
        volumesInCubicMm = np.sum(posteriors, axis=0) * volumeOfOneVoxel

        # Write intracranial volume
        sbtiv = icv(zip(*[names, volumesInCubicMm]))
        with open(os.path.join(self.savePath, 'sbtiv.stats'), 'w') as fid:
            fid.write('# Measure Intra-Cranial, %.6f, mm^3\n' % sbtiv)

        # Write structural volumes
        with open(os.path.join(self.savePath, 'samseg.stats'), 'w') as fid:
            fid.write('# Measure %s, %.6f, mm^3\n' % ('Intra-Cranial', sbtiv))
            for volume, name in zip(volumesInCubicMm, names):
                fid.write('# Measure %s, %.6f, mm^3\n' % (name, volume))

        # Write structural volumes in a csv
        with open(os.path.join(self.savePath, 'samseg.csv'), 'w') as fid:
            fid.write('ROI,volume_mm3,volume_ICV_x1000\n');
            fid.write('%s,%.6f,1000\n' % ('Intra-Cranial', sbtiv))
            for volume, name in zip(volumesInCubicMm, names):
                fid.write('%s,%.6f,%.6f\n' % (name, volume,1000*volume/sbtiv))

        # Write out a freesurfer-style stats file (good for use with asegstats2table)
        with open(os.path.join(self.savePath, 'samseg.fs.stats'), 'w') as fid:
            fid.write('# Measure EstimatedTotalIntraCranialVol, eTIV, Estimated Total Intracranial Volume, %0.6f, mm^3\n'%(sbtiv));
            # Could add other measures here like total brain volume
            fid.write('# ColHeaders  Index SegId NVoxels Volume_mm3 StructName\n');
            k = 0;
            voxsize = self.voxelSpacing[0]*self.voxelSpacing[1]*self.voxelSpacing[2];
            # Sort them by seg index like with the aseg.stats. Exclude Unknown, but keep WM and Cortex
            seglist = [];
            k = 0;
            for volume, name in zip(volumesInCubicMm, names):
                if(name == 'Unknown'): 
                    k = k+1;
                    continue
                idx = fslabels[k];
                seglist.append([idx,volume,name]);
                k = k+1;
            seglist.sort()
            k = 0;
            for seg in seglist:
                idx = seg[0];
                volume = seg[1];
                name = seg[2];
                nvox = round(volume/voxsize,2);
                fid.write('%3d %4d %7d %9.1f %s\n' % (k+1,idx,nvox,volume,name));
                k = k+1;

        return volumesInCubicMm

    def saveWarpField(self, filename):
        # extract node positions in image space
        nodePositions = self.getMesh(
            self.modelSpecifications.atlasFileName,
            self.transform,
            initialDeformation=self.deformation,
            initialDeformationMeshCollectionFileName=self.deformationAtlasFileName
        ).points

        # extract geometries
        source = sf.load_volume(self.imageFileNames[0]).geom
        target = sf.load_volume(os.path.join(self.atlasDir, 'template.nii')).geom

        # extract vox-to-vox template transform
        # TODO: Grabbing the transform from the saved .mat file in either the cross or base
        # directory is pretty messy. Ideally the affine matrix should be stored in this class
        # for both cross-sectional and longitudinal models. Also, it's important to note that
        # longitudinal timepoints might only be aligned in RAS space, not voxel space, so
        # the cached vox->vox transform computed from the base image should be converted for
        # the appropriate image geometries
        matricesFileName = os.path.join(self.savePath, 'template_transforms.mat')
        if not os.path.isfile(matricesFileName):
            matricesFileName = os.path.join(os.path.dirname(self.savePath), 'base', 'template_transforms.mat')
        matrix = scipy.io.loadmat(matricesFileName)['imageToImageTransformMatrix']

        # rasterize the final node coordinates (in image space) using the initial template mesh
        mesh = self.getMesh(self.modelSpecifications.atlasFileName)
        coordmap = mesh.rasterize_values(target.shape, nodePositions)

        # the rasterization is a bit buggy and some voxels are not filled - mark these as invalid
        invalid = np.any(coordmap == 0, axis=-1)
        coordmap[invalid, :] = -1

        # adjust for the offset introduced by volume cropping
        coordmap[~invalid, :] += [slc.start for slc in self.cropping]
        
        # writing a GCA morph requires the fsbindings to be built - this is usually only the case
        # for complete development builds or in distributed freesurfer releases
        try:
            import fsbindings
        except ImportError:
            raise ImportError('the fsbindings package is required to save a GCA morph, but it is not available')

        # write the morph
        fsbindings.write_gca_morph(coordmap, matrix, source, target, filename)

    def saveGaussianProbabilities(self, probabilitiesPath):
        # Make output directory
        os.makedirs(probabilitiesPath, exist_ok=True)
          
        # Get the class priors as dictated by the current mesh position
        mesh = self.getMesh(self.modelSpecifications.atlasFileName, self.transform,
                            initialDeformation=self.deformation,
                            initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)
        mergedAlphas = kvlMergeAlphas( mesh.alphas, self.classFractions )
        mesh.alphas = mergedAlphas
        classPriors = mesh.rasterize(self.imageBuffers.shape[0:3], -1)
        classPriors = classPriors[self.mask, :] / ( 2**16 - 1 )

        # Get bias field corrected data
        self.biasField.downSampleBasisFunctions([1, 1, 1])
        biasFields = self.biasField.getBiasFields()
        data = self.imageBuffers[self.mask, :] - biasFields[self.mask, :]

        # Compute gaussian priors, (normalized) likelihoods and posteriors
        gaussianPosteriors, _ = self.gmm.getGaussianPosteriors( data, classPriors )
        gaussianLikelihoods, _ = self.gmm.getGaussianPosteriors( data, classPriors, priorWeight=0 )
        gaussianPriors, _ = self.gmm.getGaussianPosteriors( data, classPriors, dataWeight=0 )

        # Cycle through gaussians and write volumes
        classNames = [param.mergedName for param in self.modelSpecifications.sharedGMMParameters]
        numberOfGaussiansPerClass = [param.numberOfComponents for param in self.modelSpecifications.sharedGMMParameters]
        for classNumber, className in enumerate(classNames):
            numComponents = numberOfGaussiansPerClass[classNumber]
            for componentNumber in range(numComponents):
                # 
                gaussianNumber = int(np.sum(numberOfGaussiansPerClass[:classNumber])) + componentNumber

                # Write volume
                probabilities = np.zeros( ( *( self.imageBuffers.shape[:3] ), 3 ), dtype=np.float32 )
                probabilities[ self.mask, 0 ] = gaussianPosteriors[ :, gaussianNumber ]
                probabilities[ self.mask, 1 ] = gaussianPriors[ :, gaussianNumber ]
                probabilities[ self.mask, 2 ] = gaussianLikelihoods[ :, gaussianNumber ]
                basename = className
                if numComponents > 1:
                    basename += '-%d' % (componentNumber + 1)
                self.writeImage(probabilities, os.path.join(probabilitiesPath, basename + '.mgz'))

    def saveGaussianStats(self):

        # Determine gaussians classes
        classNames = [param.mergedName for param in self.modelSpecifications.sharedGMMParameters]
        numberOfGaussiansPerClass = [param.numberOfComponents for param in self.modelSpecifications.sharedGMMParameters]
        maxNameSize = len(max(classNames, key=len)) + 2
        space = 20

        with open(os.path.join(self.savePath, 'gaussians.rescaled.txt'), 'w') as fid:

            # Write header
            modes = ''.join([mode.rjust(space) for mode in self.modeNames])
            fid.write(f"{'# class':<{maxNameSize}}{modes}\n")

            # Cycle through gaussians
            for classNumber, className in enumerate(classNames):
                numComponents = numberOfGaussiansPerClass[classNumber]
                for componentNumber in range(numComponents):
                    gaussianNumber = int(np.sum(numberOfGaussiansPerClass[:classNumber])) + componentNumber
                    basename = f'{className}-{componentNumber + 1}' if numComponents > 1 else className

                    stats = []
                    for contrastNumber in range(len((self.modeNames))):

                        # Get gaussian information
                        log_mean = self.gmm.means[gaussianNumber, contrastNumber]
                        log_variance = self.gmm.variances[gaussianNumber, contrastNumber, contrastNumber]

                        # Convert from log-space
                        scaling = self.scalingFactors[contrastNumber]
                        mean = np.exp(log_mean) * scaling
                        std = np.sqrt(log_variance) * mean

                        stats.append(f'{mean:.2f} +/- {std:.2f}'.rjust(space))

                    # Write class
                    fid.write(basename.ljust(maxNameSize) + ''.join(stats) + '\n')

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
        downSampledTransform = gems.KvlTransform(requireNumpyArray(downSamplingTransformMatrix @ self.transform.as_numpy_array))

        # Get the mesh
        downSampledMesh, downSampledInitialDeformationApplied = self.getMesh(atlasFileName,
                                                                             downSampledTransform,
                                                                             self.modelSpecifications.K,
                                                                             self.deformation,
                                                                             self.deformationAtlasFileName,
                                                                             returnInitialDeformationApplied=True)

        return downSampledImageBuffers, downSampledMask, downSampledMesh, downSampledInitialDeformationApplied, downSampledTransform

    def initializeBiasField(self):

        # Our bias model is a linear combination of a set of basis functions. We are using so-called "DCT-II" basis functions,
        # i.e., the lowest few frequency components of the Discrete Cosine Transform.
        self.biasField = BiasField(self.imageBuffers.shape[0:3],
                                   self.modelSpecifications.biasFieldSmoothingKernelSize / self.voxelSpacing, photo_mode=(self.dissectionPhoto is not None))

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
            # When working with 3D reconstructed photos, we don't downsample in z
            if self.dissectionPhoto is not None:
                downSamplingFactors[2] = 1

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
                print('iterationNumber: %d' % iterationNumber)

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
                                                              self.gmm.means, self.gmm.variances, downSampledMask,
                                                              photo_mode=(self.dissectionPhoto is not None))
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
                levelHistory['deformation'] = self.deformation.copy()
                levelHistory['deformationAtlasFileName'] = self.deformationAtlasFileName
                levelHistory['historyOfCost'] = historyOfCost
                levelHistory['priorsAtEnd'] = downSampledClassPriors
                levelHistory['posteriorsAtEnd'] = downSampledGaussianPosteriors
                self.optimizationHistory.append(levelHistory)
                

        # End resolution level loop

    def computeFinalSegmentation(self):
        # Get the final mesh
        mesh = self.getMesh(self.modelSpecifications.atlasFileName, self.transform,
                            initialDeformation=self.deformation,
                            initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)

        # Get the priors as dictated by the current mesh position
        priors = mesh.rasterize(self.imageBuffers.shape[0:3], -1)
        priors = priors[self.mask, :]

        if self.ignoreUnknownPriors:
            unknown_search_strings = next(s for s in self.modelSpecifications.sharedGMMParameters if s.mergedName == 'Unknown').searchStrings
            unknown_label = 0  # TODO: should we assume that 'Unknown' is always zero?
            for label, name in enumerate(self.modelSpecifications.names):
                for search_string in unknown_search_strings:
                    if search_string in name:
                        priors[:, unknown_label] += priors[:, label]
                        priors[:, label] = 0

        # In dissection photos, we merge the choroid with the lateral ventricle
        if self.dissectionPhoto is not None:
            for n in range(len(self.modelSpecifications.names)):
                if self.modelSpecifications.names[n]=='Left-Lateral-Ventricle':
                    llv = n
                elif self.modelSpecifications.names[n]=='Left-choroid-plexus':
                    lcp = n
                elif self.modelSpecifications.names[n]=='Right-Lateral-Ventricle':
                    rlv = n
                elif self.modelSpecifications.names[n]=='Right-choroid-plexus':
                    rcp = n
            if self.dissectionPhoto=='left' or self.dissectionPhoto=='both':
                priors[:, llv] +=  priors[:, lcp]
                priors[:, lcp] = 0
            if self.dissectionPhoto=='right' or self.dissectionPhoto=='both':
                priors[:, rlv] +=  priors[:, rcp]
                priors[:, rcp] = 0

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
