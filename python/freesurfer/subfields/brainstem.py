import os
import shutil
import numpy as np
import scipy.ndimage
import freesurfer as fs

from freesurfer import samseg
from freesurfer.subfields import utils
from freesurfer.subfields.base import MeshModel


class BrainstemSubstructures(MeshModel):

    def __init__(self, **kwargs):

        atlasDir = os.path.join(fs.fshome(), 'average/BrainstemSS/atlas')
        super().__init__(atlasDir=atlasDir, **kwargs)

        # Segmentation mesh-fitting parameters
        self.cheatingMeshSmoothingSigmas = [3.0]
        self.cheatingMaxIterations = [300]

        # Image mesh-fitting parameters
        self.meshSmoothingSigmas = [2, 1, 0]
        self.imageSmoothingSigmas = [0, 0, 0]
        self.maxIterations = [7, 5, 3]

        # Let's not smooth the target mask at all
        self.atlasTargetSmoothing = None
        self.cheatingAlphaMaskStrel = 5
        self.alphaMaskStrel = 0

    def preprocess_images(self):
        """
        Preprocess the input seg and images
        """

        # Define a few hardcoded label constants
        self.BRAINSTEM = 16
        self.DElabelLeft = 28
        self.DElabelRight = 60

        # Atlas alignment target is a masked segmentation
        mask = (self.inputSeg.data == self.BRAINSTEM).astype('float32') * 255
        self.atlasAlignmentTarget = self.inputSeg.copy(mask)

        # This will be the synthetic image used for initial mesh fitting
        self.synthImage = self.inputSeg.copy()
        labels = [self.BRAINSTEM, 7, 8, 15, 28, 46, 47, 60]
        mask = np.isin(self.synthImage.data, labels)
        self.synthImage.data[mask] = 255
        self.synthImage.write(os.path.join(self.tempDir, 'synthImage.mgz'))

        # And also used for image cropping around the brainstem
        brainstemMask = (self.inputSeg.data == self.BRAINSTEM) | (self.inputSeg.data == self.DElabelLeft) | (self.inputSeg.data == self.DElabelRight)
        fixedMargin = int(np.round(15 / np.mean(self.inputSeg.voxsize)))
        imageCropping = self.inputSeg.copy(brainstemMask).bbox(margin=fixedMargin)

        # dilate mask substantially
        radius = int(np.round(5 / self.resolution))
        mask = scipy.ndimage.morphology.binary_dilation(self.inputSeg.data > 0, utils.spherical_strel(radius), border_value=1)

        # Mask and convert to the target resolution
        images = []
        imageMask = self.synthImage.copy(mask)
        for i, image in enumerate(self.inputImages):

            # FS python library does not have cubic interpolation yet, so we'll use mri_convert
            tempFile = os.path.join(self.tempDir, 'tempImage.mgz')
            image[imageCropping].write(tempFile)
            utils.run(f'mri_convert {tempFile} {tempFile} -odt float -rt cubic -vs {self.resolution} {self.resolution} {self.resolution}')
            image = fs.Volume.read(tempFile)
            
            # Resample and apply the image mask in high-resolution target space
            imageMask = imageMask.resample_like(image, interp_method='nearest')
            image.data[imageMask.data == 0] = 0            
            images.append(image.data)

        # Define the pre-processed target image
        self.processedImage = image.copy(np.stack(images, axis=-1))

    def postprocess_segmentation(self):
        """
        Post-process the segmentation and computed volumes.
        """
        segFilePrefix = os.path.join(self.outDir, f'brainstemSsLabels{self.fileSuffix}')
        
        # Recode segmentation
        A = self.discreteLabels.copy()
        A.data[A.data < 170] = 0
        mask = utils.get_largest_cc(A.data > 0)
        A.data[mask == 0] = 0
        A.write(segFilePrefix + '.mgz')
        A.resample_like(self.inputSeg, interp_method='nearest').write(segFilePrefix + '.FSvoxelSpace.mgz')

        # Prune the volumes to what we care about
        validLabels = ['MAC_Medulla', 'MAC_Pons', 'MAC_Midbrain', 'MAC_Sup_Cerebellum_Ped', 'Medulla', 'Pons', 'Midbrain', 'SCP']
        self.volumes = {name: vol for name, vol in self.volumes.items() if name in validLabels}

        # Sum up the total volume
        self.volumes['Whole_brainstem'] = np.sum(list(self.volumes.values()))

        # Write the volumes
        self.write_volumes(segFilePrefix + '.volumes.txt')

    def get_cheating_label_groups(self):
        """
        Return a group (list of lists) of label names that determine the
        class reductions for the initial segmentation-fitting stage.
        """
        labelGroups = [
            # Label group 1
            ['MAC_Medulla', 'MAC_Pons', 'MAC_Midbrain', 'Left-VentralDC', '4th-Ventricle',
             'Left-Cerebellum-White-Matter', 'Left-Cerebellum-Cortex', 'MAC_Sup_Cerebellum_Ped',
             'Medulla', 'Pons', 'SCP', 'Midbrain'],
            # Label group 2
            ['Left-Caudate', 'Left-Accumbens-area', 'Left-Pallidum', '3rd-Ventricle', 'Left-Putamen',
             'Left-Thalamus-Proper', 'Left-Amygdala', 'Left-Lateral-Ventricle', 'Left-choroid-plexus', 'Left-Hippocampus',
             'Left-Cerebral-White-Matter', 'Left-Cerebral-Cortex', 'Background-tissue', 'Background-CSF', 'Background'],
        ]
        return labelGroups

    def get_cheating_gaussians(self, sameGaussianParameters):
        """
        Return a tuple of (means, variances) for the initial segmentation-fitting stage.
        """
        means = np.array([255.0, 1.0])
        variances = np.array([1.0, 1.0])
        return (means, variances)

    def get_label_groups(self):
        """
        Return a group (list of lists) of label names that determine the class reductions for
        the primary image-fitting stage.
        """
        labelGroups = [
            ['MAC_Medulla', 'MAC_Pons', 'MAC_Midbrain', 'MAC_Sup_Cerebellum_Ped', 'Left-VentralDC', 'Medulla', 'Pons', 'SCP', 'Midbrain'],  # Brainstem structures
            ['3rd-Ventricle', 'Left-Lateral-Ventricle', 'Background-CSF', '4th-Ventricle'],  # CSF structures
            ['Left-Amygdala', 'Left-Cerebral-Cortex', 'Left-Hippocampus'],  # Gray matter structures
            ['Left-Caudate'],  # Caudate
            ['Left-Accumbens-area'],  # Accumbens area
            ['Left-Pallidum'], # Pallidum
            ['Left-Putamen'],  # Putamen
            ['Left-Thalamus-Proper'],  # Thalamus
            ['Left-choroid-plexus'],  # Choroid plexus
            ['Left-Cerebral-White-Matter'],  # Cerebral white matter
            ['Background-tissue', 'Background'],  # Background: misc tissue
            ['Left-Cerebellum-White-Matter'],  # cerebellum white matter
            ['Left-Cerebellum-Cortex'],  # cerebellum cortex
        ]
        return labelGroups

    def get_gaussian_hyps(self, sameGaussianParameters, mesh):
        """
        Return a tuple of (meanHyps, nHyps) for Gaussian parameter estimation.
        """

        # TODO this needs to be adapted for multi-image cases as well as masking
        DATA = self.inputImages[0]

        nHyper = np.zeros(len(sameGaussianParameters))
        meanHyper = np.zeros(len(sameGaussianParameters))
        for g in range(len(sameGaussianParameters)):
            
            labels = np.array(sameGaussianParameters[g])

            if any((labels == 3) | (labels == 17) | (labels == 18)):  # gray matter
                listMask = [3, 42, 17, 53, 18, 54]
            elif any(labels == 2):  # white matter
                listMask = [2, 41]
            elif any((labels == 178) | (labels == 34458) | (labels == 28)):  # brainstem + diencephalon
                listMask = [16, 28, 60]
            elif any(labels == 4):  # CSF
                listMask = [4, 43, 14, 15]
            elif any(labels == 11):  # caudate
                listMask = [11, 50]
            elif any(labels == 26):  # accumbens
                listMask = [26, 58]
            elif any(labels == 13):  # pallidum
                listMask = [13, 52]
            elif any(labels == 12):  # putamen
                listMask = [12, 51]
            elif any(labels == 10):  # thalamus
                listMask = [10, 49]
            elif any(labels == 31):  # choroid
                listMask = [31, 63]
            elif any(labels == 0):  # background
                listMask = [0]
            elif any(labels == 7):  # cerebellum WM
                listMask = [7, 46]
            elif any(labels == 8):  # cerebellum CT
                listMask = [8, 47]
            else:
                listMask = []

            if len(listMask) > 0:
                MASK = np.zeros(DATA.data.shape, dtype='bool')
                for l in range(len(listMask)):
                    MASK = MASK | (self.inputSeg.data == listMask[l])  # TODO I think this uses an unmodified seg!
                MASK = scipy.ndimage.morphology.binary_erosion(MASK, utils.spherical_strel(1), border_value=1)
                total_mask = MASK & (DATA.data > 0)
                data = DATA.data[total_mask]
                meanHyper[g] = np.median(data)
                nHyper[g] = 10 + 0.1 * len(data) / np.prod(DATA.voxsize)

        # If any NaN, replace by background
        # ATH: I don't there would ever be NaNs here?
        nans = np.isnan(meanHyper)
        meanHyper[nans] = 55
        nHyper[nans] = 10

        return (meanHyper, nHyper)
