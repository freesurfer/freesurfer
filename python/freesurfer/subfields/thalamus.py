import os
import shutil
import numpy as np
import scipy.ndimage
import freesurfer as fs

from freesurfer import samseg
from freesurfer.subfields import utils
from freesurfer.subfields.core import MeshModel


class ThalamicNuclei(MeshModel):

    def __init__(self, **kwargs):
        atlasDir = os.path.join(fs.fshome(), 'average/ThalamicNuclei/atlas')
        super().__init__(atlasDir=atlasDir, **kwargs)

        # Model thalamus with two components
        self.useTwoComponents = True

        # Segmentation mesh-fitting parameters
        self.cheatingMeshSmoothingSigmas = [3.0, 2.0]
        self.cheatingMaxIterations = [300, 150]

        # Image mesh-fitting parameters
        self.meshSmoothingSigmas = [1.5, 1.125, 0.75, 0]
        self.imageSmoothingSigmas = [0, 0, 0, 0]
        self.maxIterations = [7, 5, 5, 3]

        # Longitudinal mesh-fitting parameters
        self.longMeshSmoothingSigmas = [[1.5, 1.125, 0.75], [1.125, 0.75, 0]]
        self.longImageSmoothingSigmas = [[0, 0, 0], [0, 0, 0]]
        self.longMaxIterations = [[7, 5, 3], [3, 2, 1]]

        # When creating the smooth atlas alignment target, dilate before eroding
        self.atlasTargetSmoothing = 'forward'

    def preprocess_images(self):
        """
        Preprocess the input seg and images
        """

        # Define a few hardcoded label constants
        self.THlabelLeft = 10
        self.THlabelRight = 49
        self.DElabelLeft = 28
        self.DElabelRight = 60

        # Atlas alignment target is a masked segmentation
        match_labels = [self.THlabelLeft, self.THlabelRight, self.DElabelLeft, self.DElabelRight]
        mask = np.isin(self.inputSeg.data, match_labels).astype('float32') * 255
        self.atlasAlignmentTarget = self.inputSeg.copy(mask)

        # Now, the idea is to refine the transform based on the thalamus + ventral DE
        # First, we prepare a modifided SEG that we'll segment
        data = self.inputSeg.data

        # There's a bunch of labels in the SEG that we don't have in our atlas
        # So we'll have to get rid of those
        data[data == 5]  = 4   # left-inf-lat-vent -> left-lat-vent
        data[data == 44] = 4   # right-inf-lat-vent -> left-lat-vent
        data[data == 14] = 4   # 3rd vent -> left-lat-vent
        data[data == 15] = 4   # 4th vent -> LV (we're killing brainstem anyway)
        data[data == 17] = 3   # left HP -> left cortex
        data[data == 53] = 3   # right HP -> left cortex
        data[data == 18] = 3   # left amygdala -> left cortex
        data[data == 54] = 3   # right amygdala -> left cortex
        data[data == 24] = 4   # CSF -> left-lat-vent
        data[data == 30] = 2   # left-vessel -> left WM
        data[data == 62] = 2   # right-vessel -> left WM
        data[data == 72] = 4   # 5th ventricle -> left-lat-vent
        data[data == 77] = 2   # WM hippoint -> left WM
        data[data == 80] = 0   # non-WM hippo -> background
        data[data == 85] = 0   # optic chiasm -> background
        data[data > 250] = 2   # CC labels -> left WM

        # Next we want to remove hemi-specific lables, so we convert right labels to left
        data[data == 41] = 2   # WM
        data[data == 42] = 3   # CT
        data[data == 43] = 4   # LV
        data[data == 46] = 7   # cerebellum WM
        data[data == 47] = 8   # cerebellum CT
        data[data == 50] = 11  # CA
        data[data == 51] = 12  # PU
        data[data == 52] = 13  # PA
        data[data == 58] = 26  # AA
        data[data == 63] = 31  # CP

        # Remove a few remainders
        removal_mask = np.isin(data, [44, 62, 63, 41, 42, 43, 50, 51, 52, 53, 54, 58])
        data[removal_mask] = 0

        # And convert background to 1
        data[data == 0] = 1

        # Now, create a mask with DE merged into thalamus. This will be the
        # synthetic image used for initial mesh fitting
        segMerged = self.inputSeg.copy()
        segMerged.data[segMerged.data == self.DElabelLeft] = self.THlabelLeft
        segMerged.data[segMerged.data == self.DElabelRight] = self.THlabelRight
        self.synthImage = segMerged

        # And also used for image cropping around the thalamus
        thalamicMask = (segMerged.data == self.THlabelLeft) | (segMerged.data == self.THlabelRight)
        fixedMargin = int(np.round(15 / np.mean(self.inputSeg.voxsize)))
        imageCropping = segMerged.copy(thalamicMask).bbox(margin=fixedMargin)

        # Lastly, use it to make the image mask
        struct = np.ones((3, 3, 3))
        mask = scipy.ndimage.morphology.binary_dilation(self.synthImage.data > 1, structure=struct, iterations=2)
        imageMask = self.synthImage.copy(mask)

        # Mask and convert to the target resolution
        images = []
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

        # Recode segmentation
        A = self.discreteLabels.copy()
        A.data[(A.data < 100) & (A.data != 10) & (A.data != 49) ] = 0

        # Kill reticular labels
        leftReticular = self.labelMapping.search('Left-R', exact=True)
        rightReticular = self.labelMapping.search('Right-R', exact=True)
        A.data[A.data == leftReticular] = 0
        A.data[A.data == rightReticular] = 0

        # Get only connected components (sometimes the two thalami are not connected)
        left = utils.get_largest_cc((A.data < 8200) & ((A.data > 100) | (A.data == self.THlabelLeft)))
        right = utils.get_largest_cc((A.data > 8200) | (A.data == self.THlabelRight))
        cc_mask = left | right
        A.data[cc_mask == 0] = 0

        segFilePrefix = os.path.join(self.outDir, f'ThalamicNuclei{self.fileSuffix}')
        A.write(segFilePrefix + '.mgz')
        A.resample_like(self.inputSeg, interp_method='nearest').write(segFilePrefix + '.FSvoxelSpace.mgz')

        # Prune the volumes to what we care about (also let's leave reticular 'R' out)
        validLabels = ['L-Sg', 'LGN', 'MGN', 'PuI', 'PuM', 'H', 'PuL',
                       'VPI', 'PuA', 'MV(Re)', 'Pf', 'CM', 'LP', 'VLa', 'VPL', 'VLp',
                       'MDm', 'VM', 'CeM', 'MDl', 'Pc', 'MDv', 'Pv', 'CL', 'VA', 'VPM',
                       'AV', 'VAmc', 'Pt', 'AD', 'LD']
        isValid = lambda name: (name.replace('Left-', '') in validLabels) or (name.replace('Right-', '') in validLabels)
        self.volumes = {name: vol for name, vol in self.volumes.items() if isValid(name)}

        # Sum up the total volumes per hemisphere 
        self.volumes['Left-Whole_thalamus'] = np.sum([vol for name, vol in self.volumes.items() if name.startswith('Left')])
        self.volumes['Right-Whole_thalamus'] = np.sum([vol for name, vol in self.volumes.items() if name.startswith('Right')])

        # Write the volumes
        self.write_volumes(segFilePrefix + '.volumes.txt')

    def get_cheating_label_groups(self):
        """
        Return a group (list of lists) of label names that determine the
        class reductions for the initial segmentation-fitting stage.
        """
        labelGroups = [
            ['Unknown'],
            ['Left-Cerebral-White-Matter'],
            ['Left-Cerebral-Cortex'],
            ['Left-Cerebellum-Cortex'],
            ['Left-Cerebellum-White-Matter'],
            ['Brain-Stem'],
            ['Left-Lateral-Ventricle'],
            ['Left-choroid-plexus'],
            ['Left-Putamen'],
            ['Left-Pallidum'],
            ['Left-Accumbens-area'],
            ['Left-Caudate'],
        ]
        thalamicLabels = [
            'L-Sg', 'LGN', 'MGN', 'PuI', 'PuM', 'H', 'PuL',
            'VPI', 'PuA', 'R', 'MV(Re)', 'Pf', 'CM', 'LP', 'VLa',
            'VPL', 'VLp', 'MDm', 'VM', 'CeM', 'MDl', 'Pc', 'MDv', 'Pv',
            'CL', 'VA', 'VPM', 'AV', 'VAmc', 'Pt', 'AD', 'LD', 'VentralDC'
        ]
        labelGroups.append(['Left-' + label for label in thalamicLabels])
        labelGroups.append(['Right-' + label for label in thalamicLabels])
        return labelGroups

    def get_cheating_gaussians(self, sameGaussianParameters):
        """
        Return a tuple of (means, variances) for the initial segmentation-fitting stage.
        """
        means = np.zeros(len(sameGaussianParameters))
        variances = 0.01 * np.ones(len(sameGaussianParameters))
        for i in range(len(sameGaussianParameters)):
            label = sameGaussianParameters[i][0]
            if label >= 8100 and label < 8200:
                means[i] = self.THlabelLeft  # left thalamic nuclei + DE -> left TH
            elif label >= 8200:
                means[i] = self.THlabelRight  # right thalamic nuclei + DE -> left TH
            elif label == 0:
                means[i] = 1  # background is 1 instead of 0
            else:
                means[i] = label
        return (means, variances)

    def get_label_groups(self):
        """
        Return a group (list of lists) of label names that determine the class reductions for
        the primary image-fitting stage.
        """
        labelGroups = [
            ['Unknown'],
            ['Left-Cerebral-White-Matter', 'Left-R', 'Right-R'],
            ['Left-Cerebral-Cortex'],
            ['Left-Cerebellum-Cortex'],
            ['Left-Cerebellum-White-Matter'],
            ['Brain-Stem'],
            ['Left-Lateral-Ventricle'],
            ['Left-choroid-plexus'],
            ['Left-Putamen'],
            ['Left-Pallidum'],
            ['Left-Accumbens-area'],
            ['Left-Caudate'],
            ['Left-VentralDC', 'Right-VentralDC'],
        ]

        # Configure left/right thalamic labels
        thalamicLabels = [
            'L-Sg', 'LGN', 'MGN', 'PuI', 'PuM', 'H', 'PuL', 'VPI', 'PuA', 'MV(Re)', 'Pf',
            'CM', 'LP', 'VLa', 'VPL', 'VLp', 'MDm', 'VM', 'CeM', 'MDl', 'Pc', 'MDv', 'Pv',
            'CL', 'VA', 'VPM', 'AV', 'VAmc', 'Pt', 'AD', 'LD',
        ]
        labelGroups.append([f'{side}-{label}' for side in ('Left', 'Right') for label in thalamicLabels])
        return labelGroups

    def get_gaussian_hyps(self, sameGaussianParameters, mesh):
        """
        Return a tuple of (meanHyps, nHyps) for Gaussian parameter estimation.
        """
        nHyper = np.zeros(len(sameGaussianParameters))
        meanHyper = np.zeros(len(sameGaussianParameters))

        # TODO this needs to be adapted for multi-image cases (with masking)
        DATA = self.inputImages[0]

        for g in range(len(sameGaussianParameters)):
            
            labels = np.array(sameGaussianParameters[g])

            if any(labels > 8225):  # thalamus
                listMask = [10, 49]
            elif any(labels == 28): # VDE
                listMask = [28, 60]
            elif any(labels == 0):  # background
                listMask = [1]
            else:
                listMask = labels
            
            if len(listMask) > 0:
                MASK = np.zeros(DATA.data.shape, dtype='bool')
                for l in range(len(listMask)):
                    # Ensure that this uses a modified segmentation
                    MASK = MASK | (self.inputSeg.data == listMask[l])
                radius = np.round(1 / np.mean(DATA.voxsize))
                MASK = scipy.ndimage.morphology.binary_erosion(MASK, utils.spherical_strel(radius), border_value=1)
                total_mask = MASK & (DATA.data > 0)
                data = DATA.data[total_mask]
                meanHyper[g] = np.median(data)
                if any(labels == 28):
                    # Special case: VDE is kind of bimodal in FreeSurfer
                    nHyper[g] = 10
                else:
                    nHyper[g] = 10 + len(data) * np.prod(DATA.voxsize) / (self.resolution ** 3)

        # If any NaN, replace by background
        # ATH: I don't there would ever be NaNs here?
        nans = np.isnan(meanHyper)
        meanHyper[nans] = 55
        nHyper[nans] = 10

        return (meanHyper, nHyper)

    def get_second_label_groups(self):
        """
        Return a group (list of lists) of label names that determine the class reductions for the
        second-component of the primary image-fitting stage.
        """
        labelGroups = [
            ['Unknown'],
            ['Left-Cerebral-White-Matter', 'Left-R', 'Right-R'],
            ['Left-Cerebral-Cortex'],
            ['Left-Cerebellum-Cortex'],
            ['Left-Cerebellum-White-Matter'],
            ['Brain-Stem'],
            ['Left-Lateral-Ventricle'],
            ['Left-choroid-plexus'],
            ['Left-Putamen'],
            ['Left-Pallidum'],
            ['Left-Accumbens-area'],
            ['Left-Caudate'],
            ['Left-VentralDC', 'Right-VentralDC'],
            ['Left-L-Sg', 'Left-LGN', 'Left-MGN', 'Left-H',
            'Left-VPI', 'Left-MV(Re)', 'Left-Pf', 'Left-CM', 'Left-LP', 'Left-VLa', 'Left-VPL', 'Left-VLp',
            'Left-VM', 'Left-CeM', 'Left-Pc', 'Left-MDv', 'Left-Pv', 'Left-CL', 'Left-VA', 'Left-VPM',
            'Left-AV', 'Left-VAmc', 'Left-Pt', 'Left-AD', 'Left-LD', 'Right-L-Sg', 'Right-LGN', 'Right-MGN', 'Right-H',
            'Right-VPI', 'Right-MV(Re)', 'Right-Pf', 'Right-CM', 'Right-LP', 'Right-VLa', 'Right-VPL', 'Right-VLp',
            'Right-VM', 'Right-CeM', 'Right-Pc', 'Right-MDv', 'Right-Pv', 'Right-CL', 'Right-VA', 'Right-VPM',
            'Right-AV', 'Right-VAmc', 'Right-Pt', 'Right-AD', 'Right-LD'],
            ['Left-PuA', 'Left-PuI', 'Left-PuL', 'Left-PuM', 'Left-MDl', 'Left-MDm',
            'Right-PuA', 'Right-PuI', 'Right-PuL', 'Right-PuM', 'Right-MDl', 'Right-MDm']
        ]
        return labelGroups

    def get_second_gaussian_hyps(self, sameGaussianParameters, meanHyper, nHyper):
        """
        Return a tuple of (meanHyps, nHyps) for Gaussian parameter estimation in the second-component
        of the primary image-fitting stage.
        """
        WMind = 1
        GMind = 2
        ThInt = meanHyper[-1]

        # TODO this needs to be enabled with non-T1s are used
        if True:
            # Lateral, brighter
            nHyper[-1] = 25
            meanHyper[-1] = ThInt + 5
            # Medial, darker
            nHyper = np.append(nHyper, 25)
            meanHyper = np.append(meanHyper, ThInt - 5)
        else:
            nHyper[-1] = 25
            nHyper = np.append(nHyper, 25)
            # Lateral, more WM-ish (e.g., darker, in FGATIR)
            meanHyper[-1] = ThInt * (0.95 + 0.1 * (meanHyper[WMind] >= meanHyper[GMind]))
            # Medial, more GM-ish (e.g., brighter, in FGATIR)
            meanHyper = np.append(meanHyper, ThInt * (0.95 + 0.1 * (meanHyper[WMind] < meanHyper[GMind])))
        return (meanHyper, nHyper)
