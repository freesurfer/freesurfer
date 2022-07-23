import os
import shutil
import numpy as np
import scipy.ndimage
import scipy.stats
import surfa as sf

from freesurfer import samseg
from freesurfer.subregions import utils
from freesurfer.subregions.core import MeshModel


class HippoAmygdalaSubfields(MeshModel):

    def __init__(self, side, wmParcFileName, resolution=0.33333, **kwargs):

        atlasDir = os.path.join(os.environ.get('FREESURFER_HOME'), 'average/HippoSF/atlas')
        super().__init__(atlasDir=atlasDir, resolution=resolution, **kwargs)

        # This is a hippocampus-specific setting to specify
        # which hemisphere to segment
        self.side = side

        # Segmentation mesh-fitting parameters
        self.cheatingMeshSmoothingSigmas = [3.0, 2.0]
        self.cheatingMaxIterations = [300, 150]

        # Image mesh-fitting parameters for cross-sectional processing
        self.meshSmoothingSigmas = [1.5, 0.75, 0]
        self.imageSmoothingSigmas = [0, 0, 0]
        self.maxIterations = [7, 5, 3]

        # Image mesh-fitting parameters for longitudinal processing
        self.longMeshSmoothingSigmas = [[1.5, 0.75], [0.75, 0]]
        self.longImageSmoothingSigmas = [[0, 0], [0, 0]]
        self.longMaxIterations = [[6, 3], [2, 1]]

        # Cache some useful info 
        self.wmparc = sf.load_volume(wmParcFileName)

        # When creating the smooth atlas alignment target, erode before dilating
        self.atlasTargetSmoothing = 'backward'

    def preprocess_images(self):

        # Define a few hardcoded label constants
        self.HippoLabelLeft = 17
        self.HippoLabelRight = 53

        # Check the hemi
        if self.side not in ('left', 'right'):
            sf.system.fatal(f'Hemisphere must be either `left` or `right`, but got `{self.side}`')

        # Flip hemi for alignment
        if self.side == 'right':
            atlasImage = sf.load_volume(self.atlasDumpFileName)
            affine = atlasImage.geom.vox2world.matrix.copy()
            affine[0, :] *= -1
            atlasImage.geom.vox2world = affine
            self.atlasDumpFileName = os.path.join(self.tempDir, 'flippedAtlasDump.mgz')
            atlasImage.save(self.atlasDumpFileName)

        # Atlas alignment target is a masked segmentation
        sideHippoLabel = self.HippoLabelLeft if self.side == 'left' else self.HippoLabelRight
        match_labels = [sideHippoLabel, sideHippoLabel + 1]
        mask = np.isin(self.inputSeg.data, match_labels).astype('float32') * 255
        self.atlasAlignmentTarget = self.inputSeg.new(mask)

        # Now, the idea is to refine the transform based on the thalamus + ventral DE
        # First, we prepare a modifided ASEG that we'll segment
        data = self.inputSeg.data.copy(order='K')

        # There's a bunch of labels in the SEG that we don't have in our atlas
        # So we'll have to get rid of those
        data[data == 15] = 0  # 4th vent -> background (we're killing brainstem anyway...)
        data[data == 16] = 0  # get rid of brainstem
        data[data == 7] = 0   # get rid of left cerebellum WM ...
        data[data == 8] = 0   # ... and of left cerebellum CT
        data[data == 46] = 0  # get rid of right cerebellum WM ...
        data[data == 47] = 0  # ... and of right cerebellum CT
        data[data == 80] = 0  # non-WM hippo -> background
        data[data == 85] = 0  # optic chiasm -> background
        data[data == 72] = 4  # 5th ventricle -> left-lat-vent

        if self.side == 'left':
            data[data == 5] = 4   # left-inf-lat-vent -> left-lat-vent
            data[data == 30] = 2  # left-vessel -> left  WM
            data[data == 14] = 4  # 3rd vent -> left-lat-vent
            data[data == 24] = 4  # CSF -> left-lat-vent
            data[data == 77] = 2  # WM hippoint -> left WM
            data[data > 250] = 2  # CC labels -> left WM

            removal_mask = np.isin(data, [44, 62, 63, 41, 42, 43, 49, 50, 51, 52, 53, 54, 58, 60])
            data[removal_mask] = 0

        else:
            bu = data.copy(order='K')
            data.fill(0)
            data[bu == 44] = 4   # right-inf-lat-vent -> left-lat-vent
            data[bu == 62] = 2   # right-vessel -> left  WM
            data[bu == 14] = 4   # 3rd vent -> left-lat-vent
            data[bu == 24] = 4   # CSF -> left-lat-vent
            data[bu == 77] = 2   # WM hippoint -> left WM
            data[bu > 250] = 2   # CC labels -> left WM
            # left to right
            data[bu == 41] = 2   # WM
            data[bu == 42] = 3   # CT
            data[bu == 43] = 4   # LV
            data[bu == 49] = 10  # TH
            data[bu == 50] = 11  # CA
            data[bu == 51] = 12  # PU
            data[bu == 52] = 13  # PA
            data[bu == 53] = 17  # HP
            data[bu == 54] = 18  # AM
            data[bu == 58] = 26  # AA
            data[bu == 60] = 28  # DC
            data[bu == 63] = 31  # CP

        # And convert background to 1
        data[data == 0] = 1

        segMerged = self.inputSeg.new(data)

        # We now merge hippo, amygdala, and cortex. This will be the
        # synthetic image used for initial mesh fitting
        self.synthImage = segMerged.copy()
        self.synthImage[self.synthImage == 17] = 3
        self.synthImage[self.synthImage == 18] = 3

        # And also used for image cropping around the thalamus
        fixedMargin = int(np.round(15 / np.mean(self.inputSeg.geom.voxsize)))
        imageCropping = segMerged.new(self.inputSeg.data == sideHippoLabel).bbox(margin=fixedMargin)

        # Let's dilate this mask (ATH not totally sure why there are two masks here)
        mask = scipy.ndimage.morphology.binary_dilation(segMerged > 1, structure=np.ones((3, 3, 3)), iterations=2)
        mergedMaskDilated = segMerged.new(mask)

        # Lastly, use it to make the image mask
        mask = (segMerged > 16) & (segMerged < 19)
        imageMask = self.synthImage.new(mask)

        # Dilate the mask
        dilatedMask = scipy.ndimage.morphology.binary_dilation(mask, structure=np.ones((3, 3, 3)), iterations=5)
        self.maskDilated5mm = self.synthImage.new(dilatedMask)

        # Mask and convert to the target resolution
        images = []
        for i, image in enumerate(self.inputImages):

            # FS python library does not have cubic interpolation yet, so we'll use mri_convert
            tempFile = os.path.join(self.tempDir, 'tempImage.mgz')
            image[imageCropping].save(tempFile)
            utils.run(f'mri_convert {tempFile} {tempFile} -odt float -rt cubic -vs {self.resolution} {self.resolution} {self.resolution}')
            image = sf.load_volume(tempFile)

            # Resample and apply the first mask in high-resolution target space
            maskTempFile = os.path.join(self.tempDir, 'asegModBinDilatedResampled.mgz')
            mergedMaskDilated.save(maskTempFile)
            utils.run(f'mri_convert {maskTempFile} {maskTempFile} -odt float -rt nearest -rl {tempFile}')
            mask = sf.load_volume(maskTempFile)

            image[mask == 0] = 0

            # Resample and apply the second mask in high-resolution target space
            maskTempFile = os.path.join(self.tempDir, 'hippoMaskResampled.mgz')
            imageMask.save(maskTempFile)
            utils.run(f'mri_convert {maskTempFile} {maskTempFile} -odt float -rt interpolate -rl {tempFile}')
            mask = sf.load_volume(maskTempFile) >= 0.5
            mask = scipy.ndimage.morphology.binary_dilation(mask, structure=np.ones((3, 3, 3)), iterations=int(np.round(3 / self.resolution)))
            image[mask == 0] = 0
            self.longMask = mask

            images.append(image.data)

        # Define the pre-processed target image
        self.processedImage = image.new(np.stack(images, axis=-1))

    def postprocess_segmentation(self):
        """
        Post-process the segmentation and computed volumes.
        """
        segFilePrefix = os.path.join(self.outDir, f'{self.side[0]}h.hippoAmygLabels{self.fileSuffix}')

        A = self.discreteLabels.copy()
        A[A < 200] = 0
        A[(A > 246) & (A < 7000)] = 0
        A[A == 201] = 0
        mask = utils.get_largest_cc(A > 0)
        A[mask == 0] = 0
        A.save(segFilePrefix + '.mgz')
        A.resample_like(self.inputSeg, method='nearest').save(segFilePrefix + '.FSvoxelSpace.mgz')

        # Write merged versions to disk as well

        # First: tail, body, head
        HippoBodyLabel = 231
        HippoHeadLabel = 232

        HPbodyList = ['subiculum-body', 'CA1-body', 'presubiculum-body', 'molecular_layer_HP-body',
                      'CA3-body', 'GC-ML-DG-body', 'CA4-body', 'fimbria']
        HPheadList = ['subiculum-head', 'presubiculum-head', 'CA1-head', 'parasubiculum',
                      'molecular_layer_HP-head', 'GC-ML-DG-head', 'CA4-head', 'CA3-head', 'HATA']

        HPbodyList = [name.lower() for name in HPbodyList]
        HPheadList = [name.lower() for name in HPheadList]

        B = A.copy()
        for c, name in enumerate(self.names):
            name = name.lower().replace(' ', '')

            if name in HPbodyList:
                B[B == self.FreeSurferLabels[c]] = HippoBodyLabel

            if name in HPheadList:
                B[B == self.FreeSurferLabels[c]] = HippoHeadLabel

        # Kill the fissure
        B[B == 215] = 0
        B.save(segFilePrefix + '.HBT.mgz')
        B.resample_like(self.inputSeg, method='nearest').save(segFilePrefix + '.HBT.FSvoxelSpace.mgz')

        # Second: head and body of each subfield
        C = A.copy()
        C[(A == 233) | (A == 234)] = 204  # presubiculum
        C[(A == 235) | (A == 236)] = 205  # subiculum
        C[(A == 237) | (A == 238)] = 206  # CA1
        C[(A == 239) | (A == 240)] = 208  # CA3
        C[(A == 241) | (A == 242)] = 209  # CA4
        C[(A == 243) | (A == 244)] = 210  # GC-DG
        C[(A == 245) | (A == 246)] = 214  # ML
        C.save(segFilePrefix + '.FS60.mgz')
        C.resample_like(self.inputSeg, method='nearest').save(segFilePrefix + '.FS60.FSvoxelSpace.mgz')

        # Third: same as above, but we get rid of internal labels
        D = C.copy()
        D[D == 210] = 209  # GC-DG -> CA4

        # Molecular layer: replace by nearest label that is not background or fissure
        cropping = D.new(D == 214).bbox(margin=2)
        V = D[cropping]
        labels = [l for l in np.unique(V) if l not in (0, 214, 215)]
        mask = V == 214
        for i, label in enumerate(labels):
            dmap = scipy.ndimage.distance_transform_edt(V != label)
            if i == 0:
                mini = dmap[mask]
                seg = label * np.ones(mini.shape)
            else:
                dist = dmap[mask]
                m = dist < mini
                mini[m] = dist[m]
                seg[m] = label
        V[mask] = seg
        D[cropping] = V
        D.save(segFilePrefix + '.CA.mgz')
        D.resample_like(self.inputSeg, method='nearest').save(segFilePrefix + '.CA.FSvoxelSpace.mgz')

        # Extract hippocampal volumes
        validLabels = ['subiculum-body', 'subiculum-head', 'Hippocampal_tail', 'molecular_layer_HP-body', 'molecular_layer_HP-head', 'hippocampal-fissure',
                       'GC-ML-DG-body', 'GC-ML-DG-head', 'CA4-body', 'CA4-head', 'presubiculum-body', 'presubiculum-head', 'CA1-body', 'CA1-head',
                       'parasubiculum', 'fimbria', 'CA3-body', 'CA3-head', 'HATA']
        hippoVolumes = {name: vol for name, vol in self.volumes.items() if name in validLabels}

        # Compute total hippocampal volume (ignore fissure)
        validLabels = [name for name in validLabels if name != 'hippocampal-fissure']
        hippoVolumes['Whole_hippocampus'] = np.sum([vol for name, vol in hippoVolumes.items() if name in validLabels])

        # Compute total hippocampal body volume
        validLabels = ['subiculum-body', 'CA1-body', 'presubiculum-body', 'molecular_layer_HP-body', 'CA3-body', 'GC-ML-DG-body', 'CA4-body', 'fimbria']
        hippoVolumes['Whole_hippocampal_body'] = np.sum([vol for name, vol in hippoVolumes.items() if name in validLabels])

        # Compute total hippocampal head volume
        validLabels = ['subiculum-head', 'presubiculum-head', 'CA1-head', 'parasubiculum', 'molecular_layer_HP-head', 'GC-ML-DG-head', 'CA4-head', 'CA3-head', 'HATA']
        hippoVolumes['Whole_hippocampal_head'] = np.sum([vol for name, vol in hippoVolumes.items() if name in validLabels])

        # Write hippo volumes
        self.write_volumes(os.path.join(self.outDir, f'{self.side[0]}h.hippoSfVolumes{self.fileSuffix}.txt'), hippoVolumes)

        # Extract amygdala volumes
        validLabels = ['Left-Amygdala', 'Lateral-nucleus', 'Paralaminar-nucleus', 'Basal-nucleus', 'Hippocampal-amygdala-transition-HATA',
                       'Accessory-Basal-nucleus', 'Amygdala-background', 'Corticoamygdaloid-transitio', 'Central-nucleus',
                       'Cortical-nucleus', 'Medial-nucleus', 'Anterior-amygdaloid-area-AAA']
        amygdalaVolumes = {name: vol for name, vol in self.volumes.items() if name in validLabels}

        # Compute total amygdala volume
        amygdalaVolumes['Whole_amygdala'] = np.sum(list(amygdalaVolumes.values()))

        # Write amygdala volumes
        self.write_volumes(os.path.join(self.outDir, f'{self.side[0]}h.amygNucVolumes{self.fileSuffix}.txt'), amygdalaVolumes)

    def get_cheating_label_groups(self):
        """
        Return a group (list of lists) of label names that determine the
        class reductions for the initial segmentation-fitting stage.
        """
        labelGroups = [
            ['Left-Cerebral-Cortex', 'Left-Hippocampus', 'alveus', 'subiculum-body', 'subiculum-head', 'Hippocampal_tail' ,
            'molecular_layer_HP-body', 'molecular_layer_HP-head', 'GC-ML-DG-body', 'GC-ML-DG-head',
            'CA4-body', 'CA4-head', 'CA1-body', 'CA1-head', 'CA3-body', 'CA3-head', 'HATA', 'fimbria',
            'presubiculum-body', 'presubiculum-head', 'parasubiculum', 'Left-hippocampus-intensity-abnormality',
            'Left-Amygdala', 'Lateral-nucleus', 'Paralaminar-nucleus', 'Basal-nucleus',
            'Hippocampal-amygdala-transition-HATA', 'Accessory-Basal-nucleus', 'Amygdala-background',
            'Corticoamygdaloid-transitio', 'Central-nucleus', 'Cortical-nucleus', 'Medial-nucleus',
            'Anterior-amygdaloid-area-AAA'],
            ['Left-Cerebral-White-Matter'],
            ['Left-Lateral-Ventricle'],
            ['Left-choroid-plexus'],
            ['Background', 'hippocampal-fissure', 'Background-CSF', 'Background-vessels', 'Background-tissue', 'Unknown'],
            ['Left-VentralDC'],
            ['Left-Putamen'],
            ['Left-Pallidum'],
            ['Left-Thalamus-Proper'],
            ['Left-Accumbens-area'],
            ['Left-Caudate'],
            ['SUSPICIOUS']
        ]
        return labelGroups

    def get_cheating_gaussians(self, sameGaussianParameters):
        """
        Return a tuple of (means, variances) for the initial segmentation-fitting stage.
        """
        means = np.zeros(len(sameGaussianParameters))
        variances = 0.01 * np.ones(len(sameGaussianParameters))
        for l in range(len(sameGaussianParameters)):
            labels = np.array(sameGaussianParameters[l])

            if any((labels >= 200) & (labels <= 226) & (labels != 215)):
                means[l] = 3  # Hippo SF > Hippo
            elif any((labels >= 7000)):
                means[l] = 3  # Amygdala Subnuclei -> Amygdala
            elif any(labels == 0):
                means[l] = 1  # Background is 1 instead of 0
            elif any(labels == 999):
                means[l] = 55
                variances[l] = 55 ** 2  # This is the generic `suspicious` label we use for cysts
            else:
                means[l] = labels[0]
        return (means, variances)

    def get_label_groups(self):
        """
        Return a group (list of lists) of label names that determine the class reductions for
        the primary image-fitting stage.
        """
        if not self.highResImage:
            labelGroups = [
                ['Left-Cerebral-Cortex', 'Left-Hippocampus', 'Left-Amygdala', 'subiculum-head', 'subiculum-body',
                'Hippocampal_tail', 'GC-ML-DG-head', 'GC-ML-DG-body', 'CA4-head', 'CA4-body', 'presubiculum-head', 'presubiculum-body',
                'CA1-head', 'CA1-body', 'parasubiculum', 'CA3-head', 'CA3-body', 'HATA', 'Lateral-nucleus', 'Paralaminar-nucleus',
                'Basal-nucleus', 'Hippocampal-amygdala-transition-HATA', 'Accessory-Basal-nucleus', 'Amygdala-background',
                'Corticoamygdaloid-transitio', 'Central-nucleus', 'Cortical-nucleus', 'Medial-nucleus',
                'Anterior-amygdaloid-area-AAA', 'molecular_layer_HP-body', 'molecular_layer_HP-head']]
        else:
            labelGroups = [
                ['Left-Cerebral-Cortex', 'Left-Hippocampus', 'Left-Amygdala', 'subiculum-head', 'subiculum-body',
                'Hippocampal_tail', 'GC-ML-DG-head', 'GC-ML-DG-body', 'CA4-head', 'CA4-body', 'presubiculum-head', 'presubiculum-body',
                'CA1-head', 'CA1-body', 'parasubiculum', 'CA3-head', 'CA3-body', 'HATA', 'Lateral-nucleus', 'Paralaminar-nucleus',
                'Basal-nucleus', 'Hippocampal-amygdala-transition-HATA', 'Accessory-Basal-nucleus', 'Amygdala-background',
                'Corticoamygdaloid-transitio', 'Central-nucleus', 'Cortical-nucleus', 'Medial-nucleus',
                'Anterior-amygdaloid-area-AAA'],
                ['molecular_layer_HP-body', 'molecular_layer_HP-head']]

        labelGroups.append(['Left-Cerebral-White-Matter', 'fimbria'])
        labelGroups.append(['alveus'])
        labelGroups.append(['Left-Lateral-Ventricle', 'Background-CSF', 'SUSPICIOUS', 'Left-hippocampus-intensity-abnormality'])
        labelGroups.append(['hippocampal-fissure'])
        labelGroups.append(['Left-Pallidum'])
        labelGroups.append(['Left-Putamen'])
        labelGroups.append(['Left-Caudate'])
        labelGroups.append(['Left-Thalamus-Proper'])
        labelGroups.append(['Left-choroid-plexus'])
        labelGroups.append(['Left-VentralDC'])
        labelGroups.append(['Left-Accumbens-area'])
        labelGroups.append(['Unknown', 'Background-tissue'])
        return labelGroups

    def get_gaussian_hyps(self, sameGaussianParameters, mesh):
        """
        Return a tuple of (meanHyps, nHyps) for Gaussian parameter estimation.
        """
        DATA = self.inputImages[0]
        WMPARC = self.wmparc
        mask = (WMPARC == 0) & (self.maskDilated5mm == 0)
        WMPARC[mask] = -1

        nHyper = np.zeros(len(sameGaussianParameters))
        meanHyper = np.zeros(len(sameGaussianParameters))
        for g in range(len(sameGaussianParameters)):
            labels = np.array(sameGaussianParameters[g])

            if any((labels == 3) | (labels == 17) | (labels == 18) | (labels > 7000) | (labels == 226)):
                listMask = 17 if self.side == 'left' else 53
            elif any(labels == 2):
                listMask = [3006, 3007, 3016] if self.side == 'left' else [4006, 4007, 4016]
            elif any(labels == 26):
                listMask = 26 if self.side == 'left' else 58
            elif any(labels == 4):
                listMask = 4 if self.side == 'left' else 43
            elif any(labels == 0):
                listMask = [0]
            elif any(labels == 13):
                listMask = 13 if self.side == 'left' else 52
            elif any(labels == 12):
                listMask = 12 if self.side == 'left' else 51
            elif any(labels == 11):
                listMask = 11 if self.side == 'left' else 50
            elif any(labels == 10):
                listMask = 10 if self.side == 'left' else 49
            elif any(labels == 31):
                listMask = 31 if self.side == 'left' else 63
            elif any(labels == 28):
                listMask = 28 if self.side == 'left' else 60
            else:
                listMask = None

            if listMask is not None:
                if isinstance(listMask, int):
                    listMask = [listMask]

                MASK = np.zeros(DATA.shape, dtype='bool')
                for l in range(len(listMask)):
                    MASK = MASK | (WMPARC == listMask[l])

                radius = np.round(1 / np.mean(DATA.geom.voxsize))
                MASK = scipy.ndimage.morphology.binary_erosion(MASK, utils.spherical_strel(radius), border_value=1)
                total_mask = MASK & (DATA > 0)
                data = DATA[total_mask]

                meanHyper[g] = np.median(data)
                nHyper[g] = 10 + len(data) * np.prod(DATA.geom.voxsize) / (self.resolution ** 3)

        # If any NaN, replace by background
        # ATH: I don't there would ever be NaNs here?
        nans = np.isnan(meanHyper)
        meanHyper[nans] = 55
        nHyper[nans] = 10

        # Here's the part where we simulate partial voluming!
        print('Estimating typical intensities of alveus')

        WMind = None
        GMind = None
        ALind = None
        MLind = None
        FISSind = None
        CSFind = None
        for g in range(len(sameGaussianParameters)):
            labels = np.array(sameGaussianParameters[g])
            if any(labels == 2):
                WMind = g
            if any(labels == 3):
                GMind = g
            if any(labels == 201):
                ALind = g
            if any(labels == 245) and self.highResImage:
                MLind = g
            if any(labels == 215):
                FISSind = g
            if any(labels == 4):
                CSFind = g

        imageShape = (mesh.points.max(axis=0) + 1.5).astype(int)
        priors = mesh.rasterize(imageShape)

        L = np.argmax(priors, axis=-1)
        maskPriors = (priors.sum(-1) / 65535) > 0.97

        I = np.zeros(imageShape)
        for l in range(len(sameGaussianParameters)):
            if l == ALind or l == MLind:
                I[L == l] = meanHyper[WMind]
            elif l == FISSind:
                I[L == l] = meanHyper[CSFind]
            else:
                I[L == l] = meanHyper[l]

        I[maskPriors == 0] = 0
        sigma = np.mean(DATA.geom.voxsize) / (2.355 * self.resolution)
        I_PV = scipy.ndimage.gaussian_filter(I, sigma)

        if ALind is not None:
            data = I_PV[L == ALind]
            # It's multimodal, so regular median won't cut it
            kde = scipy.stats.gaussian_kde(data.flatten())
            v = np.linspace(data.min(), data.max(), 1000)
            meanHyper[ALind] = np.median(v[np.argmax(kde(v))])  # median of argmax??
            nHyper[ALind] = (nHyper[GMind] + nHyper[WMind]) / 2

        if self.highResImage:
            data = I_PV[L == MLind]
            meanHyper[MLind] = np.median(I_PV[L == MLind])
            nHyper[MLind] = (nHyper[WMind] + nHyper[GMind]) / 2

        if FISSind is not None:
            meanHyper[FISSind] = np.median(I_PV[L == FISSind])
            nHyper[FISSind] = (nHyper[CSFind] + nHyper[GMind]) / 2

        return (meanHyper, nHyper)
