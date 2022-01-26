import os
import collections
import numpy as np

from . import fshome


class LookupTable(collections.OrderedDict):
    """
    TODOC
    """

    class Element:
        def __init__(self, name, color):
            self.name = name
            if color is None:
                color = np.random.randint(0, 256, size=(3))
            if len(color) == 3:
                color = np.append(color, 255)
            elif len(color) != 4:
                raise ValueError('Color must be a 4-element RGBA uchar array')
            self.color = np.array(color, dtype='uint8')

    def __str__(self):
        col1 = len(str(max(self.keys()))) + 1
        col2 = max([len(elt.name) for elt in self.values()]) + 2
        lines = []
        for idx, elt in self.items():
            colorstr = '(' + ' '.join([str(c).ljust(3) for c in elt.color]) + ')'
            lines.append(str(idx).ljust(col1) + elt.name.ljust(col2) + colorstr)
        return '\n'.join(lines)

    def add(self, index, name, color=None):
        self[index] = LookupTable.Element(name, color)

    def search(self, name, exact=False):
        if exact:
            return next((idx for idx, elt in self.items() if name == elt.name), None)
        else:
            allcaps = name.upper()
            return [idx for idx, elt in self.items() if allcaps in elt.name.upper()]

    @classmethod
    def read(cls, filename):
        lut = cls()
        with open(filename, 'r') as file:
            lines = file.readlines()
        for line in lines:
            split = line.lstrip().split()
            if split and not split[0].startswith('#'):
                idx, name = split[:2]
                if len(split) >= 5:
                    color = list(map(int, split[2:6]))
                    color[3] = 255 - color[3]  # invert alpha value
                else:
                    color = None
                lut.add(int(idx), name, color)
        return lut

    @classmethod
    def read_default(cls):
        return cls.read(os.path.join(fshome(), 'FreeSurferColorLUT.txt'))

    def write(self, filename):
        col1 = len(str(max(self.keys()))) + 1  # find largest index
        col2 = max([len(elt.name) for elt in self.values()]) + 2  # find longest name
        with open(filename, 'w') as file:
            file.write('#'.ljust(col1) + 'Label Name'.ljust(col2) + 'R   G   B   A\n\n')
            for idx, elt in self.items():
                color = elt.color
                color[3] = 255 - color[3]  # invert alpha value
                colorstr = ' '.join([str(c).ljust(3) for c in color])
                file.write(str(idx).ljust(col1) + elt.name.ljust(col2) + colorstr + '\n')

    def extract(self, labels):
        """
        Extract a new LookupTable from a list of label indices.
        """
        lut = LookupTable()
        for label in labels:
            elt = self.get(label)
            if elt is None:
                raise ValueError(f'Index {label} does not exist in the LookupTable.')
            lut.add(label, elt.name, elt.color)
        return lut

    def copy_colors(self, source_lut):
        """
        Copies colors of matching label indices from a source LookupTable.
        """
        for label in self.keys():
            elt = source_lut.get(label)
            if elt is not None and elt.color is not None:
                self[label].color = elt.color

    def copy_names(self, source_lut):
        """
        Copies names of matching label indices from a source LookupTable.
        """
        for label in self.keys():
            elt = source_lut.get(label)
            if elt is not None:
                self[label].name = elt.name


class RecodingLookupTable(dict):
    """
    TODOC
    """

    def __init__(self):
        self.mapping = {}
        self.target_lut = LookupTable()
        super().__init__()

    def add_target(self, index, name=None, color=None):
        self.target_lut.add(index, name, color)


def default():
    """
    Returns the default freesurfer label lookup table.

    Returns:
        LookupTable: Default lookup table.
    """

    lut = LookupTable()

    # standard brain labels
    lut.add(0,     'Unknown',                                         [  0,   0,   0])
    lut.add(1,     'Left-Cerebral-Exterior',                          [ 70, 130, 180])
    lut.add(2,     'Left-Cerebral-White-Matter',                      [245, 245, 245])
    lut.add(3,     'Left-Cerebral-Cortex',                            [205,  62,  78])
    lut.add(4,     'Left-Lateral-Ventricle',                          [120,  18, 134])
    lut.add(5,     'Left-Inf-Lat-Vent',                               [196,  58, 250])
    lut.add(6,     'Left-Cerebellum-Exterior',                        [  0, 148,   0])
    lut.add(7,     'Left-Cerebellum-White-Matter',                    [220, 248, 164])
    lut.add(8,     'Left-Cerebellum-Cortex',                          [230, 148,  34])
    lut.add(9,     'Left-Thalamus-Unused',                            [  0, 118,  14])
    lut.add(10,    'Left-Thalamus',                                   [  0, 118,  14])
    lut.add(11,    'Left-Caudate',                                    [122, 186, 220])
    lut.add(12,    'Left-Putamen',                                    [236,  13, 176])
    lut.add(13,    'Left-Pallidum',                                   [ 12,  48, 255])
    lut.add(14,    '3rd-Ventricle',                                   [204, 182, 142])
    lut.add(15,    '4th-Ventricle',                                   [ 42, 204, 164])
    lut.add(16,    'Brain-Stem',                                      [119, 159, 176])
    lut.add(17,    'Left-Hippocampus',                                [220, 216,  20])
    lut.add(18,    'Left-Amygdala',                                   [103, 255, 255])
    lut.add(19,    'Left-Insula',                                     [ 80, 196,  98])
    lut.add(20,    'Left-Operculum',                                  [ 60,  58, 210])
    lut.add(21,    'Line-1',                                          [ 60,  58, 210])
    lut.add(22,    'Line-2',                                          [ 60,  58, 210])
    lut.add(23,    'Line-3',                                          [ 60,  58, 210])
    lut.add(24,    'CSF',                                             [ 60,  60,  60])
    lut.add(25,    'Left-Lesion',                                     [255, 165,   0])
    lut.add(26,    'Left-Accumbens-area',                             [255, 165,   0])
    lut.add(27,    'Left-Substancia-Nigra',                           [  0, 255, 127])
    lut.add(28,    'Left-VentralDC',                                  [165,  42,  42])
    lut.add(29,    'Left-Undetermined',                               [135, 206, 235])
    lut.add(30,    'Left-Vessel',                                     [160,  32, 240])
    lut.add(31,    'Left-Choroid-Plexus',                             [  0, 200, 200])
    lut.add(32,    'Left-F3orb',                                      [100,  50, 100])
    lut.add(33,    'Left-lOg',                                        [135,  50,  74])
    lut.add(34,    'Left-aOg',                                        [122, 135,  50])
    lut.add(35,    'Left-mOg',                                        [ 51,  50, 135])
    lut.add(36,    'Left-pOg',                                        [ 74, 155,  60])
    lut.add(37,    'Left-Stellate',                                   [120,  62,  43])
    lut.add(38,    'Left-Porg',                                       [ 74, 155,  60])
    lut.add(39,    'Left-Aorg',                                       [122, 135,  50])
    lut.add(40,    'Right-Cerebral-Exterior',                         [ 70, 130, 180])
    lut.add(41,    'Right-Cerebral-White-Matter',                     [245, 245, 245])
    lut.add(42,    'Right-Cerebral-Cortex',                           [205,  62,  78])
    lut.add(43,    'Right-Lateral-Ventricle',                         [120,  18, 134])
    lut.add(44,    'Right-Inf-Lat-Vent',                              [196,  58, 250])
    lut.add(45,    'Right-Cerebellum-Exterior',                       [  0, 148,   0])
    lut.add(46,    'Right-Cerebellum-White-Matter',                   [220, 248, 164])
    lut.add(47,    'Right-Cerebellum-Cortex',                         [230, 148,  34])
    lut.add(48,    'Right-Thalamus-Unused',                           [  0, 118,  14])
    lut.add(49,    'Right-Thalamus',                                  [  0, 118,  14])
    lut.add(50,    'Right-Caudate',                                   [122, 186, 220])
    lut.add(51,    'Right-Putamen',                                   [236,  13, 176])
    lut.add(52,    'Right-Pallidum',                                  [ 13,  48, 255])
    lut.add(53,    'Right-Hippocampus',                               [220, 216,  20])
    lut.add(54,    'Right-Amygdala',                                  [103, 255, 255])
    lut.add(55,    'Right-Insula',                                    [ 80, 196,  98])
    lut.add(56,    'Right-Operculum',                                 [ 60,  58, 210])
    lut.add(57,    'Right-Lesion',                                    [255, 165,   0])
    lut.add(58,    'Right-Accumbens-area',                            [255, 165,   0])
    lut.add(59,    'Right-Substancia-Nigra',                          [  0, 255, 127])
    lut.add(60,    'Right-VentralDC',                                 [165,  42,  42])
    lut.add(61,    'Right-Undetermined',                              [135, 206, 235])
    lut.add(62,    'Right-Vessel',                                    [160,  32, 240])
    lut.add(63,    'Right-Choroid-Plexus',                            [  0, 200, 221])
    lut.add(64,    'Right-F3orb',                                     [100,  50, 100])
    lut.add(65,    'Right-lOg',                                       [135,  50,  74])
    lut.add(66,    'Right-aOg',                                       [122, 135,  50])
    lut.add(67,    'Right-mOg',                                       [ 51,  50, 135])
    lut.add(68,    'Right-pOg',                                       [ 74, 155,  60])
    lut.add(69,    'Right-Stellate',                                  [120,  62,  43])
    lut.add(70,    'Right-Porg',                                      [ 74, 155,  60])
    lut.add(71,    'Right-Aorg',                                      [122, 135,  50])
    lut.add(72,    '5th-Ventricle',                                   [120, 190, 150])
    lut.add(73,    'Left-Interior',                                   [122, 135,  50])
    lut.add(74,    'Right-Interior',                                  [122, 135,  50])
    lut.add(77,    'WM-Hypointensities',                              [200,  70, 255])
    lut.add(78,    'Left-WM-Hypointensities',                         [255, 148,  10])
    lut.add(79,    'Right-WM-Hypointensities',                        [255, 148,  10])
    lut.add(80,    'Non-WM-Hypointensities',                          [164, 108, 226])
    lut.add(81,    'Left-Non-WM-Hypointensities',                     [164, 108, 226])
    lut.add(82,    'Right-Non-WM-Hypointensities',                    [164, 108, 226])
    lut.add(83,    'Left-F1',                                         [255, 218, 185])
    lut.add(84,    'Right-F1',                                        [255, 218, 185])
    lut.add(85,    'Optic-Chiasm',                                    [234, 169,  30])
    lut.add(86,    'Left_future_WMSA',                                [200, 120, 255])
    lut.add(87,    'Right_future_WMSA',                               [200, 121, 255])
    lut.add(88,    'future_WMSA',                                     [200, 122, 255])
    lut.add(96,    'Left-Amygdala-Anterior',                          [205,  10, 125])
    lut.add(97,    'Right-Amygdala-Anterior',                         [205,  10, 125])
    lut.add(98,    'Dura',                                            [160,  32, 240])
    lut.add(99,    'Lesion',                                          [255, 165,   0])
    lut.add(100,   'Left-wm-intensity-abnormality',                   [124, 140, 178])
    lut.add(101,   'Left-caudate-intensity-abnormality',              [125, 140, 178])
    lut.add(102,   'Left-putamen-intensity-abnormality',              [126, 140, 178])
    lut.add(103,   'Left-accumbens-intensity-abnormality',            [127, 140, 178])
    lut.add(104,   'Left-pallidum-intensity-abnormality',             [124, 141, 178])
    lut.add(105,   'Left-amygdala-intensity-abnormality',             [124, 142, 178])
    lut.add(106,   'Left-hippocampus-intensity-abnormality',          [124, 143, 178])
    lut.add(107,   'Left-thalamus-intensity-abnormality',             [124, 144, 178])
    lut.add(108,   'Left-VDC-intensity-abnormality',                  [124, 140, 179])
    lut.add(109,   'Right-wm-intensity-abnormality',                  [124, 140, 178])
    lut.add(110,   'Right-caudate-intensity-abnormality',             [125, 140, 178])
    lut.add(111,   'Right-putamen-intensity-abnormality',             [126, 140, 178])
    lut.add(112,   'Right-accumbens-intensity-abnormality',           [127, 140, 178])
    lut.add(113,   'Right-pallidum-intensity-abnormality',            [124, 141, 178])
    lut.add(114,   'Right-amygdala-intensity-abnormality',            [124, 142, 178])
    lut.add(115,   'Right-hippocampus-intensity-abnormality',         [124, 143, 178])
    lut.add(116,   'Right-thalamus-intensity-abnormality',            [124, 144, 178])
    lut.add(117,   'Right-VDC-intensity-abnormality',                 [124, 140, 179])
    lut.add(118,   'Epidermis',                                       [255,  20, 147])
    lut.add(119,   'Conn-Tissue',                                     [205, 179, 139])
    lut.add(120,   'SC-Fat-Muscle',                                   [238, 238, 209])
    lut.add(121,   'Cranium',                                         [200, 200, 200])
    lut.add(122,   'CSF-SA',                                          [ 74, 255,  74])
    lut.add(123,   'Muscle',                                          [238,   0,   0])
    lut.add(124,   'Ear',                                             [  0,   0, 139])
    lut.add(125,   'Adipose',                                         [173, 255,  47])
    lut.add(126,   'Spinal-Cord',                                     [133, 203, 229])
    lut.add(127,   'Soft-Tissue',                                     [ 26, 237,  57])
    lut.add(128,   'Nerve',                                           [ 34, 139,  34])
    lut.add(129,   'Bone',                                            [ 30, 144, 255])
    lut.add(130,   'AirCavity',                                       [196, 160, 128])
    lut.add(131,   'Orbital-Fat',                                     [238,  59,  59])
    lut.add(132,   'Tongue',                                          [221,  39, 200])
    lut.add(133,   'Nasal-Structures',                                [238, 174, 238])
    lut.add(134,   'Globe',                                           [255,   0,   0])
    lut.add(135,   'Teeth',                                           [ 72,  61, 139])
    lut.add(136,   'Left-Caudate-Putamen',                            [ 21,  39, 132])
    lut.add(137,   'Right-Caudate-Putamen',                           [ 21,  39, 132])
    lut.add(138,   'Left-Claustrum',                                  [ 65, 135,  20])
    lut.add(139,   'Right-Claustrum',                                 [ 65, 135,  20])
    lut.add(140,   'Cornea',                                          [134,   4, 160])
    lut.add(142,   'Diploe',                                          [221, 226,  68])
    lut.add(143,   'Vitreous-Humor',                                  [255, 255, 254])
    lut.add(144,   'Lens',                                            [ 52, 209, 226])
    lut.add(145,   'Aqueous-Humor',                                   [239, 160, 223])
    lut.add(146,   'Outer-Table',                                     [ 70, 130, 180])
    lut.add(147,   'Inner-Table',                                     [ 70, 130, 181])
    lut.add(148,   'Periosteum',                                      [139, 121,  94])
    lut.add(149,   'Endosteum',                                       [224, 224, 224])
    lut.add(150,   'R-C-S',                                           [255,   0,   0])
    lut.add(151,   'Iris',                                            [205, 205,   0])
    lut.add(152,   'SC-Adipose-Muscle',                               [238, 238, 209])
    lut.add(153,   'SC-Tissue',                                       [139, 121,  94])
    lut.add(154,   'Orbital-Adipose',                                 [238,  59,  59])
    lut.add(155,   'Left-IntCapsule-Ant',                             [238,  59,  59])
    lut.add(156,   'Right-IntCapsule-Ant',                            [238,  59,  59])
    lut.add(157,   'Left-IntCapsule-Pos',                             [ 62,  10, 205])
    lut.add(158,   'Right-IntCapsule-Pos',                            [ 62,  10, 205])

    # infant labels
    lut.add(159,   'Left-Cerebral-WM-unmyelinated',                   [  0, 118,  14])
    lut.add(160,   'Right-Cerebral-WM-unmyelinated',                  [  0, 118,  14])
    lut.add(161,   'Left-Cerebral-WM-myelinated',                     [220, 216,  21])
    lut.add(162,   'Right-Cerebral-WM-myelinated',                    [220, 216,  21])
    lut.add(163,   'Left-Subcortical-Gray-Matter',                    [122, 186, 220])
    lut.add(164,   'Right-Subcortical-Gray-Matter',                   [122, 186, 220])
    lut.add(165,   'Skull',                                           [120, 120, 120])
    lut.add(166,   'Posterior-fossa',                                 [ 14,  48, 255])
    lut.add(167,   'Scalp',                                           [166,  42,  42])
    lut.add(168,   'Hematoma',                                        [121,  18, 134])
    lut.add(169,   'Left-Basal-Ganglia',                              [236,  13, 127])

    # brainstem
    lut.add(170,   'brainstem',                                       [119, 159, 176])
    lut.add(171,   'DCG',                                             [119,   0, 176])
    lut.add(172,   'Vermis',                                          [119, 100, 176])
    lut.add(173,   'Midbrain',                                        [242, 104,  76])
    lut.add(174,   'Pons',                                            [206, 195,  58])
    lut.add(175,   'Medulla',                                         [119, 159, 176])
    lut.add(176,   'Right-Basal-Ganglia',                             [236,  13, 126])
    lut.add(177,   'Vermis-White-Matter',                             [119,  50, 176])
    lut.add(178,   'SCP',                                             [142, 182,   0])
    lut.add(179,   'Floculus',                                        [ 19, 100, 176])

    lut.add(180,   'Left-Cortical-Dysplasia',                         [ 73,  61, 139])
    lut.add(181,   'Right-Cortical-Dysplasia',                        [ 73,  62, 139])
    lut.add(182,   'CblumNodulus',                                    [ 10, 100, 176])
    lut.add(192,   'Corpus_Callosum',                                 [250, 255,  50])
    lut.add(193,   'Left-hippocampal_fissure',                        [  0, 196, 255])
    lut.add(194,   'Left-CADG-head',                                  [255, 164, 164])
    lut.add(195,   'Left-subiculum',                                  [196, 196,   0])
    lut.add(196,   'Left-fimbria',                                    [  0, 100, 255])
    lut.add(197,   'Right-hippocampal_fissure',                       [128, 196, 164])
    lut.add(198,   'Right-CADG-head',                                 [  0, 126,  75])
    lut.add(199,   'Right-subiculum',                                 [128,  96,  64])
    lut.add(200,   'Right-fimbria',                                   [  0,  50, 128])
    lut.add(201,   'alveus',                                          [255, 204, 153])
    lut.add(202,   'perforant_pathway',                               [255, 128, 128])
    lut.add(203,   'parasubiculum',                                   [175, 175,  75])
    lut.add(204,   'presubiculum',                                    [ 64,   0,  64])
    lut.add(205,   'subiculum',                                       [  0,   0, 255])
    lut.add(206,   'CA1',                                             [255,   0,   0])
    lut.add(207,   'CA2',                                             [128, 128, 255])
    lut.add(208,   'CA3',                                             [  0, 128,   0])
    lut.add(209,   'CA4',                                             [196, 160, 128])
    lut.add(210,   'GC-DG',                                           [ 32, 200, 255])
    lut.add(211,   'HATA',                                            [128, 255, 128])
    lut.add(212,   'fimbria',                                         [204, 153, 204])
    lut.add(213,   'lateral_ventricle',                               [121,  17, 136])
    lut.add(214,   'molecular_layer_HP',                              [128,   0,   0])
    lut.add(215,   'hippocampal_fissure',                             [128,  32, 255])
    lut.add(216,   'entorhinal_cortex',                               [255, 204, 102])
    lut.add(217,   'molecular_layer_subiculum',                       [128, 128, 128])
    lut.add(218,   'Amygdala',                                        [104, 255, 255])
    lut.add(219,   'Cerebral_White_Matter',                           [  0, 226,   0])
    lut.add(220,   'Cerebral_Cortex',                                 [205,  63,  78])
    lut.add(221,   'Inf_Lat_Vent',                                    [197,  58, 250])
    lut.add(222,   'Perirhinal',                                      [ 33, 150, 250])
    lut.add(223,   'Cerebral_White_Matter_Edge',                      [226,   0,   0])
    lut.add(224,   'Background',                                      [100, 100, 100])
    lut.add(225,   'Ectorhinal',                                      [197, 150, 250])
    lut.add(226,   'HP_tail',                                         [170, 170, 255])
    lut.add(227,   'Polymorphic-Layer',                               [128, 255, 128])
    lut.add(228,   'Intracellular-Space',                             [204, 153, 204])
    lut.add(231,   'HP_body',                                         [  0, 255,   0])
    lut.add(232,   'HP_head',                                         [255,   0,   0])
    lut.add(233,   'presubiculum-head',                               [ 32,   0,  32])
    lut.add(234,   'presubiculum-body',                               [ 64,   0,  64])
    lut.add(235,   'subiculum-head',                                  [  0,   0, 175])
    lut.add(236,   'subiculum-body',                                  [  0,   0, 255])
    lut.add(237,   'CA1-head',                                        [175,  75,  75])
    lut.add(238,   'CA1-body',                                        [255,   0,   0])
    lut.add(239,   'CA3-head',                                        [  0,  80,   0])
    lut.add(240,   'CA3-body',                                        [  0, 128,   0])
    lut.add(241,   'CA4-head',                                        [120,  90,  50])
    lut.add(242,   'CA4-body',                                        [196, 160, 128])
    lut.add(243,   'GC-ML-DG-head',                                   [ 75, 125, 175])
    lut.add(244,   'GC-ML-DG-body',                                   [ 32, 200, 255])
    lut.add(245,   'molecular_layer_HP-head',                         [100,  25,  25])
    lut.add(246,   'molecular_layer_HP-body',                         [128,   0,   0])
    lut.add(247,   'FreezeSurface',                                   [ 10, 100, 100])
    lut.add(250,   'Fornix',                                          [255,   0,   0])
    lut.add(251,   'CC_Posterior',                                    [  0,   0,  64])
    lut.add(252,   'CC_Mid_Posterior',                                [  0,   0, 112])
    lut.add(253,   'CC_Central',                                      [  0,   0, 160])
    lut.add(254,   'CC_Mid_Anterior',                                 [  0,   0, 208])
    lut.add(255,   'CC_Anterior',                                     [  0,   0, 255])
    lut.add(256,   'Voxel-Unchanged',                                 [  0,   0,   0])
    lut.add(257,   'CSF-ExtraCerebral',                               [ 60,  60,  60])
    lut.add(258,   'Head-ExtraCerebral',                              [150, 150, 200])
    lut.add(259,   'Eye-Fluid',                                       [ 60,  60,  60])
    lut.add(260,   'BoneOrAir',                                       [119, 159, 176])
    lut.add(261,   'PossibleFluid',                                   [120,  18, 134])
    lut.add(262,   'Sinus',                                           [196, 160, 128])
    lut.add(263,   'Left-Eustachian',                                 [119, 159, 176])
    lut.add(264,   'Right-Eustachian',                                [119, 159, 176])
    lut.add(265,   'Left-Eyeball',                                    [ 60,  60,  60])
    lut.add(266,   'Right-Eyeball',                                   [ 60,  60,  60])

    # lymph node and vascular labels
    lut.add(331,   'Aorta',                                           [255,   0,   0])
    lut.add(332,   'Left-Common-IliacA',                              [255,  80,   0])
    lut.add(333,   'Right-Common-IliacA',                             [255, 160,   0])
    lut.add(334,   'Left-External-IliacA',                            [255, 255,   0])
    lut.add(335,   'Right-External-IliacA',                           [  0, 255,   0])
    lut.add(336,   'Left-Internal-IliacA',                            [255,   0, 160])
    lut.add(337,   'Right-Internal-IliacA',                           [255,   0, 255])
    lut.add(338,   'Left-Lateral-SacralA',                            [255,  50,  80])
    lut.add(339,   'Right-Lateral-SacralA',                           [ 80, 255,  50])
    lut.add(340,   'Left-ObturatorA',                                 [160, 255,  50])
    lut.add(341,   'Right-ObturatorA',                                [160, 200, 255])
    lut.add(342,   'Left-Internal-PudendalA',                         [  0, 255, 160])
    lut.add(343,   'Right-Internal-PudendalA',                        [  0,   0, 255])
    lut.add(344,   'Left-UmbilicalA',                                 [ 80,  50, 255])
    lut.add(345,   'Right-UmbilicalA',                                [160,   0, 255])
    lut.add(346,   'Left-Inf-RectalA',                                [255, 210,   0])
    lut.add(347,   'Right-Inf-RectalA',                               [  0, 160, 255])
    lut.add(348,   'Left-Common-IliacV',                              [255, 200,  80])
    lut.add(349,   'Right-Common-IliacV',                             [255, 200, 160])
    lut.add(350,   'Left-External-IliacV',                            [255,  80, 200])
    lut.add(351,   'Right-External-IliacV',                           [255, 160, 200])
    lut.add(352,   'Left-Internal-IliacV',                            [ 30, 255,  80])
    lut.add(353,   'Right-Internal-IliacV',                           [ 80, 200, 255])
    lut.add(354,   'Left-ObturatorV',                                 [ 80, 255, 200])
    lut.add(355,   'Right-ObturatorV',                                [195, 255, 200])
    lut.add(356,   'Left-Internal-PudendalV',                         [120, 200,  20])
    lut.add(357,   'Right-Internal-PudendalV',                        [170,  10, 200])
    lut.add(358,   'Pos-Lymph',                                       [ 20, 130, 180])
    lut.add(359,   'Neg-Lymph',                                       [ 20, 180, 130])

    lut.add(400,   'V1',                                              [206,  62,  78])
    lut.add(401,   'V2',                                              [121,  18, 134])
    lut.add(402,   'BA44',                                            [199,  58, 250])
    lut.add(403,   'BA45',                                            [  1, 148,   0])
    lut.add(404,   'BA4a',                                            [221, 248, 164])
    lut.add(405,   'BA4p',                                            [231, 148,  34])
    lut.add(406,   'BA6',                                             [  1, 118,  14])
    lut.add(407,   'BA2',                                             [120, 118,  14])
    lut.add(408,   'BA1_old',                                         [123, 186, 221])
    lut.add(409,   'BAun2',                                           [238,  13, 177])
    lut.add(410,   'BA1',                                             [123, 186, 220])
    lut.add(411,   'BA2b',                                            [138,  13, 206])
    lut.add(412,   'BA3a',                                            [238, 130, 176])
    lut.add(413,   'BA3b',                                            [218, 230,  76])
    lut.add(414,   'MT',                                              [ 38, 213, 176])
    lut.add(415,   'AIPS_AIP_l',                                      [  1, 225, 176])
    lut.add(416,   'AIPS_AIP_r',                                      [  1, 225, 176])
    lut.add(417,   'AIPS_VIP_l',                                      [200,   2, 100])
    lut.add(418,   'AIPS_VIP_r',                                      [200,   2, 100])
    lut.add(419,   'IPL_PFcm_l',                                      [  5, 200,  90])
    lut.add(420,   'IPL_PFcm_r',                                      [  5, 200,  90])
    lut.add(421,   'IPL_PF_l',                                        [100,   5, 200])
    lut.add(422,   'IPL_PFm_l',                                       [ 25, 255, 100])
    lut.add(423,   'IPL_PFm_r',                                       [ 25, 255, 100])
    lut.add(424,   'IPL_PFop_l',                                      [230,   7, 100])
    lut.add(425,   'IPL_PFop_r',                                      [230,   7, 100])
    lut.add(426,   'IPL_PF_r',                                        [100,   5, 200])
    lut.add(427,   'IPL_PFt_l',                                       [150,  10, 200])
    lut.add(428,   'IPL_PFt_r',                                       [150,  10, 200])
    lut.add(429,   'IPL_PGa_l',                                       [175,  10, 176])
    lut.add(430,   'IPL_PGa_r',                                       [175,  10, 176])
    lut.add(431,   'IPL_PGp_l',                                       [ 10, 100, 255])
    lut.add(432,   'IPL_PGp_r',                                       [ 10, 100, 255])
    lut.add(433,   'Visual_V3d_l',                                    [150,  45,  70])
    lut.add(434,   'Visual_V3d_r',                                    [150,  45,  70])
    lut.add(435,   'Visual_V4_l',                                     [ 45, 200,  15])
    lut.add(436,   'Visual_V4_r',                                     [ 45, 200,  15])
    lut.add(437,   'Visual_V5_b',                                     [227,  45, 100])
    lut.add(438,   'Visual_VP_l',                                     [227,  45, 100])
    lut.add(439,   'Visual_VP_r',                                     [227,  45, 100])

    # wm lesions
    lut.add(498,   'wmsa',                                            [143, 188, 143])
    lut.add(499,   'other_wmsa',                                      [255, 248, 220])

    # hippocampus labeling
    lut.add(500,   'right_CA2_3',                                     [ 17,  85, 136])
    lut.add(501,   'right_alveus',                                    [119, 187, 102])
    lut.add(502,   'right_CA1',                                       [204,  68,  34])
    lut.add(503,   'right_fimbria',                                   [204,   0, 255])
    lut.add(504,   'right_presubiculum',                              [221, 187,  17])
    lut.add(505,   'right_hippocampal_fissure',                       [153, 221, 238])
    lut.add(506,   'right_CA4_DG',                                    [ 51,  17,  17])
    lut.add(507,   'right_subiculum',                                 [  0, 119,  85])
    lut.add(508,   'right_fornix',                                    [ 20, 100, 200])
    lut.add(550,   'left_CA2_3',                                      [ 17,  85, 137])
    lut.add(551,   'left_alveus',                                     [119, 187, 103])
    lut.add(552,   'left_CA1',                                        [204,  68,  35])
    lut.add(553,   'left_fimbria',                                    [204,   0, 254])
    lut.add(554,   'left_presubiculum',                               [221, 187,  16])
    lut.add(555,   'left_hippocampal_fissure',                        [153, 221, 239])
    lut.add(556,   'left_CA4_DG',                                     [ 51,  17,  18])
    lut.add(557,   'left_subiculum',                                  [  0, 119,  86])
    lut.add(558,   'left_fornix',                                     [ 20, 100, 201])

    lut.add(600,   'Tumor',                                           [254, 254, 254])

    # cerebellar parcellation labels
    lut.add(601,   'Cbm_Left_I_IV',                                   [ 70, 130, 180])
    lut.add(602,   'Cbm_Right_I_IV',                                  [245, 245, 245])
    lut.add(603,   'Cbm_Left_V',                                      [205,  62,  78])
    lut.add(604,   'Cbm_Right_V',                                     [120,  18, 134])
    lut.add(605,   'Cbm_Left_VI',                                     [196,  58, 250])
    lut.add(606,   'Cbm_Vermis_VI',                                   [  0, 148,   0])
    lut.add(607,   'Cbm_Right_VI',                                    [220, 248, 164])
    lut.add(608,   'Cbm_Left_CrusI',                                  [230, 148,  34])
    lut.add(609,   'Cbm_Vermis_CrusI',                                [  0, 118,  14])
    lut.add(610,   'Cbm_Right_CrusI',                                 [  0, 118,  14])
    lut.add(611,   'Cbm_Left_CrusII',                                 [122, 186, 220])
    lut.add(612,   'Cbm_Vermis_CrusII',                               [236,  13, 176])
    lut.add(613,   'Cbm_Right_CrusII',                                [ 12,  48, 255])
    lut.add(614,   'Cbm_Left_VIIb',                                   [204, 182, 142])
    lut.add(615,   'Cbm_Vermis_VIIb',                                 [ 42, 204, 164])
    lut.add(616,   'Cbm_Right_VIIb',                                  [119, 159, 176])
    lut.add(617,   'Cbm_Left_VIIIa',                                  [220, 216,  20])
    lut.add(618,   'Cbm_Vermis_VIIIa',                                [103, 255, 255])
    lut.add(619,   'Cbm_Right_VIIIa',                                 [ 80, 196,  98])
    lut.add(620,   'Cbm_Left_VIIIb',                                  [ 60,  58, 210])
    lut.add(621,   'Cbm_Vermis_VIIIb',                                [ 60,  58, 210])
    lut.add(622,   'Cbm_Right_VIIIb',                                 [ 60,  58, 210])
    lut.add(623,   'Cbm_Left_IX',                                     [ 60,  58, 210])
    lut.add(624,   'Cbm_Vermis_IX',                                   [ 60,  60,  60])
    lut.add(625,   'Cbm_Right_IX',                                    [255, 165,   0])
    lut.add(626,   'Cbm_Left_X',                                      [255, 165,   0])
    lut.add(627,   'Cbm_Vermis_X',                                    [  0, 255, 127])
    lut.add(628,   'Cbm_Right_X',                                     [165,  42,  42])
    
    # cerebellar lobule parcellations
    lut.add(640,   'Cbm_Right_I_V_med',                               [204,   0,   0])
    lut.add(641,   'Cbm_Right_I_V_mid',                               [255,   0,   0])
    lut.add(642,   'Cbm_Right_VI_med',                                [  0,   0, 255])
    lut.add(643,   'Cbm_Right_VI_mid',                                [ 30, 144, 255])
    lut.add(644,   'Cbm_Right_VI_lat',                                [100, 212, 237])
    lut.add(645,   'Cbm_Right_CrusI_med',                             [218, 165,  32])
    lut.add(646,   'Cbm_Right_CrusI_mid',                             [255, 215,   0])
    lut.add(647,   'Cbm_Right_CrusI_lat',                             [255, 255, 166])
    lut.add(648,   'Cbm_Right_CrusII_med',                            [153,   0, 204])
    lut.add(649,   'Cbm_Right_CrusII_mid',                            [153, 141, 209])
    lut.add(650,   'Cbm_Right_CrusII_lat',                            [204, 204, 255])
    lut.add(651,   'Cbm_Right_7med',                                  [ 31, 212, 194])
    lut.add(652,   'Cbm_Right_7mid',                                  [  3, 255, 237])
    lut.add(653,   'Cbm_Right_7lat',                                  [204, 255, 255])
    lut.add(654,   'Cbm_Right_8med',                                  [ 86,  74, 147])
    lut.add(655,   'Cbm_Right_8mid',                                  [114, 114, 190])
    lut.add(656,   'Cbm_Right_8lat',                                  [184, 178, 255])
    lut.add(657,   'Cbm_Right_PUNs',                                  [126, 138,  37])
    lut.add(658,   'Cbm_Right_TONs',                                  [189, 197, 117])
    lut.add(659,   'Cbm_Right_FLOs',                                  [240, 230, 140])
    lut.add(660,   'Cbm_Left_I_V_med',                                [204,   0,   0])
    lut.add(661,   'Cbm_Left_I_V_mid',                                [255,   0,   0])
    lut.add(662,   'Cbm_Left_VI_med',                                 [  0,   0, 255])
    lut.add(663,   'Cbm_Left_VI_mid',                                 [ 30, 144, 255])
    lut.add(664,   'Cbm_Left_VI_lat',                                 [100, 212, 237])
    lut.add(665,   'Cbm_Left_CrusI_med',                              [218, 165,  32])
    lut.add(666,   'Cbm_Left_CrusI_mid',                              [255, 215,   0])
    lut.add(667,   'Cbm_Left_CrusI_lat',                              [255, 255, 166])
    lut.add(668,   'Cbm_Left_CrusII_med',                             [153,   0, 204])
    lut.add(669,   'Cbm_Left_CrusII_mid',                             [153, 141, 209])
    lut.add(670,   'Cbm_Left_CrusII_lat',                             [204, 204, 255])
    lut.add(671,   'Cbm_Left_7med',                                   [ 31, 212, 194])
    lut.add(672,   'Cbm_Left_7mid',                                   [  3, 255, 237])
    lut.add(673,   'Cbm_Left_7lat',                                   [204, 255, 255])
    lut.add(674,   'Cbm_Left_8med',                                   [ 86,  74, 147])
    lut.add(675,   'Cbm_Left_8mid',                                   [114, 114, 190])
    lut.add(676,   'Cbm_Left_8lat',                                   [184, 178, 255])
    lut.add(677,   'Cbm_Left_PUNs',                                   [126, 138,  37])
    lut.add(678,   'Cbm_Left_TONs',                                   [189, 197, 117])
    lut.add(679,   'Cbm_Left_FLOs',                                   [240, 230, 140])
    lut.add(690,   'CbmWM_Gyri_Left',                                 [122, 135,  50])
    lut.add(691,   'CbmWM_Gyri_Right',                                [122, 135,  50])
    
    lut.add(701,   'CSF-FSL-FAST',                                    [120,  18, 134])
    lut.add(702,   'GrayMatter-FSL-FAST',                             [205,  62,  78])
    lut.add(703,   'WhiteMatter-FSL-FAST',                            [  0, 225,   0])
    
    # hypothalamus labels
    lut.add(801,   'L_hypothalamus_anterior_inferior',                [250, 255,  50])
    lut.add(802,   'L_hypothalamus_anterior_superior',                [ 80, 200, 255])
    lut.add(803,   'L_hypothalamus_posterior',                        [255, 160,   0])
    lut.add(804,   'L_hypothalamus_tubular_inferior',                 [255, 160, 200])
    lut.add(805,   'L_hypothalamus_tubular_superior',                 [ 20, 180, 130])
    lut.add(806,   'R_hypothalamus_anterior_inferior',                [250, 255,  50])
    lut.add(807,   'R_hypothalamus_anterior_superior',                [ 80, 200, 255])
    lut.add(808,   'R_hypothalamus_posterior',                        [255, 160,   0])
    lut.add(809,   'R_hypothalamus_tubular_inferior',                 [255, 160, 200])
    lut.add(810,   'R_hypothalamus_tubular_superior',                 [ 20, 180, 130])

    lut.add(999,   'SUSPICIOUS',                                      [255, 100, 100])
    
    # cortical labels
    lut.add(1000,  'ctx-lh-unknown',                                  [ 25,   5,  25])
    lut.add(1001,  'ctx-lh-bankssts',                                 [ 25, 100,  40])
    lut.add(1002,  'ctx-lh-caudalanteriorcingulate',                  [125, 100, 160])
    lut.add(1003,  'ctx-lh-caudalmiddlefrontal',                      [100,  25,   0])
    lut.add(1004,  'ctx-lh-corpuscallosum',                           [120,  70,  50])
    lut.add(1005,  'ctx-lh-cuneus',                                   [220,  20, 100])
    lut.add(1006,  'ctx-lh-entorhinal',                               [220,  20,  10])
    lut.add(1007,  'ctx-lh-fusiform',                                 [180, 220, 140])
    lut.add(1008,  'ctx-lh-inferiorparietal',                         [220,  60, 220])
    lut.add(1009,  'ctx-lh-inferiortemporal',                         [180,  40, 120])
    lut.add(1010,  'ctx-lh-isthmuscingulate',                         [140,  20, 140])
    lut.add(1011,  'ctx-lh-lateraloccipital',                         [ 20,  30, 140])
    lut.add(1012,  'ctx-lh-lateralorbitofrontal',                     [ 35,  75,  50])
    lut.add(1013,  'ctx-lh-lingual',                                  [225, 140, 140])
    lut.add(1014,  'ctx-lh-medialorbitofrontal',                      [200,  35,  75])
    lut.add(1015,  'ctx-lh-middletemporal',                           [160, 100,  50])
    lut.add(1016,  'ctx-lh-parahippocampal',                          [ 20, 220,  60])
    lut.add(1017,  'ctx-lh-paracentral',                              [ 60, 220,  60])
    lut.add(1018,  'ctx-lh-parsopercularis',                          [220, 180, 140])
    lut.add(1019,  'ctx-lh-parsorbitalis',                            [ 20, 100,  50])
    lut.add(1020,  'ctx-lh-parstriangularis',                         [220,  60,  20])
    lut.add(1021,  'ctx-lh-pericalcarine',                            [120, 100,  60])
    lut.add(1022,  'ctx-lh-postcentral',                              [220,  20,  20])
    lut.add(1023,  'ctx-lh-posteriorcingulate',                       [220, 180, 220])
    lut.add(1024,  'ctx-lh-precentral',                               [ 60,  20, 220])
    lut.add(1025,  'ctx-lh-precuneus',                                [160, 140, 180])
    lut.add(1026,  'ctx-lh-rostralanteriorcingulate',                 [ 80,  20, 140])
    lut.add(1027,  'ctx-lh-rostralmiddlefrontal',                     [ 75,  50, 125])
    lut.add(1028,  'ctx-lh-superiorfrontal',                          [ 20, 220, 160])
    lut.add(1029,  'ctx-lh-superiorparietal',                         [ 20, 180, 140])
    lut.add(1030,  'ctx-lh-superiortemporal',                         [140, 220, 220])
    lut.add(1031,  'ctx-lh-supramarginal',                            [ 80, 160,  20])
    lut.add(1032,  'ctx-lh-frontalpole',                              [100,   0, 100])
    lut.add(1033,  'ctx-lh-temporalpole',                             [ 70,  70,  70])
    lut.add(1034,  'ctx-lh-transversetemporal',                       [150, 150, 200])
    lut.add(1035,  'ctx-lh-insula',                                   [255, 192,  32])
    lut.add(1100,  'ctx-lh-Unknown',                                  [  0,   0,   0])
    lut.add(1101,  'ctx-lh-Corpus_callosum',                          [ 50,  50,  50])
    lut.add(1102,  'ctx-lh-G_and_S_Insula_ONLY_AVERAGE',              [180,  20,  30])
    lut.add(1103,  'ctx-lh-G_cingulate-Isthmus',                      [ 60,  25,  25])
    lut.add(1104,  'ctx-lh-G_cingulate-Main_part',                    [ 25,  60,  60])
    lut.add(1105,  'ctx-lh-G_cuneus',                                 [180,  20,  20])
    lut.add(1106,  'ctx-lh-G_frontal_inf-Opercular_part',             [220,  20, 100])
    lut.add(1107,  'ctx-lh-G_frontal_inf-Orbital_part',               [140,  60,  60])
    lut.add(1108,  'ctx-lh-G_frontal_inf-Triangular_part',            [180, 220, 140])
    lut.add(1109,  'ctx-lh-G_frontal_middle',                         [140, 100, 180])
    lut.add(1110,  'ctx-lh-G_frontal_superior',                       [180,  20, 140])
    lut.add(1111,  'ctx-lh-G_frontomarginal',                         [140,  20, 140])
    lut.add(1112,  'ctx-lh-G_insular_long',                           [ 21,  10,  10])
    lut.add(1113,  'ctx-lh-G_insular_short',                          [225, 140, 140])
    lut.add(1114,  'ctx-lh-G_and_S_occipital_inferior',               [ 23,  60, 180])
    lut.add(1115,  'ctx-lh-G_occipital_middle',                       [180,  60, 180])
    lut.add(1116,  'ctx-lh-G_occipital_superior',                     [ 20, 220,  60])
    lut.add(1117,  'ctx-lh-G_occipit-temp_lat-Or_fusiform',           [ 60,  20, 140])
    lut.add(1118,  'ctx-lh-G_occipit-temp_med-Lingual_part',          [220, 180, 140])
    lut.add(1119,  'ctx-lh-G_occipit-temp_med-Parahippocampal_part',  [ 65, 100,  20])
    lut.add(1120,  'ctx-lh-G_orbital',                                [220,  60,  20])
    lut.add(1121,  'ctx-lh-G_paracentral',                            [ 60, 100,  60])
    lut.add(1122,  'ctx-lh-G_parietal_inferior-Angular_part',         [ 20,  60, 220])
    lut.add(1123,  'ctx-lh-G_parietal_inferior-Supramarginal_part',   [100, 100,  60])
    lut.add(1124,  'ctx-lh-G_parietal_superior',                      [220, 180, 220])
    lut.add(1125,  'ctx-lh-G_postcentral',                            [ 20, 180, 140])
    lut.add(1126,  'ctx-lh-G_precentral',                             [ 60, 140, 180])
    lut.add(1127,  'ctx-lh-G_precuneus',                              [ 25,  20, 140])
    lut.add(1128,  'ctx-lh-G_rectus',                                 [ 20,  60, 100])
    lut.add(1129,  'ctx-lh-G_subcallosal',                            [ 60, 220,  20])
    lut.add(1130,  'ctx-lh-G_subcentral',                             [ 60,  20, 220])
    lut.add(1131,  'ctx-lh-G_temporal_inferior',                      [220, 220, 100])
    lut.add(1132,  'ctx-lh-G_temporal_middle',                        [180,  60,  60])
    lut.add(1133,  'ctx-lh-G_temp_sup-G_temp_transv_and_interm_S',    [ 60,  60, 220])
    lut.add(1134,  'ctx-lh-G_temp_sup-Lateral_aspect',                [220,  60, 220])
    lut.add(1135,  'ctx-lh-G_temp_sup-Planum_polare',                 [ 65, 220,  60])
    lut.add(1136,  'ctx-lh-G_temp_sup-Planum_tempolare',              [ 25, 140,  20])
    lut.add(1137,  'ctx-lh-G_and_S_transverse_frontopolar',           [ 13,   0, 250])
    lut.add(1138,  'ctx-lh-Lat_Fissure-ant_sgt-ramus_horizontal',     [ 61,  20, 220])
    lut.add(1139,  'ctx-lh-Lat_Fissure-ant_sgt-ramus_vertical',       [ 61,  20,  60])
    lut.add(1140,  'ctx-lh-Lat_Fissure-post_sgt',                     [ 61,  60, 100])
    lut.add(1141,  'ctx-lh-Medial_wall',                              [ 25,  25,  25])
    lut.add(1142,  'ctx-lh-Pole_occipital',                           [140,  20,  60])
    lut.add(1143,  'ctx-lh-Pole_temporal',                            [220, 180,  20])
    lut.add(1144,  'ctx-lh-S_calcarine',                              [ 63, 180, 180])
    lut.add(1145,  'ctx-lh-S_central',                                [221,  20,  10])
    lut.add(1146,  'ctx-lh-S_central_insula',                         [ 21, 220,  20])
    lut.add(1147,  'ctx-lh-S_cingulate-Main_part_and_Intracingulate', [183, 100,  20])
    lut.add(1148,  'ctx-lh-S_cingulate-Marginalis_part',              [221,  20, 100])
    lut.add(1149,  'ctx-lh-S_circular_insula_anterior',               [221,  60, 140])
    lut.add(1150,  'ctx-lh-S_circular_insula_inferior',               [221,  20, 220])
    lut.add(1151,  'ctx-lh-S_circular_insula_superior',               [ 61, 220, 220])
    lut.add(1152,  'ctx-lh-S_collateral_transverse_ant',              [100, 200, 200])
    lut.add(1153,  'ctx-lh-S_collateral_transverse_post',             [ 10, 200, 200])
    lut.add(1154,  'ctx-lh-S_frontal_inferior',                       [221, 220,  20])
    lut.add(1155,  'ctx-lh-S_frontal_middle',                         [141,  20, 100])
    lut.add(1156,  'ctx-lh-S_frontal_superior',                       [ 61, 220, 100])
    lut.add(1157,  'ctx-lh-S_frontomarginal',                         [ 21, 220,  60])
    lut.add(1158,  'ctx-lh-S_intermedius_primus-Jensen',              [141,  60,  20])
    lut.add(1159,  'ctx-lh-S_intraparietal-and_Parietal_transverse',  [143,  20, 220])
    lut.add(1160,  'ctx-lh-S_occipital_anterior',                     [ 61,  20, 180])
    lut.add(1161,  'ctx-lh-S_occipital_middle_and_Lunatus',           [101,  60, 220])
    lut.add(1162,  'ctx-lh-S_occipital_superior_and_transversalis',   [ 21,  20, 140])
    lut.add(1163,  'ctx-lh-S_occipito-temporal_lateral',              [221, 140,  20])
    lut.add(1164,  'ctx-lh-S_occipito-temporal_medial_and_S_Lingual', [141, 100, 220])
    lut.add(1165,  'ctx-lh-S_orbital-H_shapped',                      [101,  20,  20])
    lut.add(1166,  'ctx-lh-S_orbital_lateral',                        [221, 100,  20])
    lut.add(1167,  'ctx-lh-S_orbital_medial-Or_olfactory',            [181, 200,  20])
    lut.add(1168,  'ctx-lh-S_paracentral',                            [ 21, 180, 140])
    lut.add(1169,  'ctx-lh-S_parieto_occipital',                      [101, 100, 180])
    lut.add(1170,  'ctx-lh-S_pericallosal',                           [181, 220,  20])
    lut.add(1171,  'ctx-lh-S_postcentral',                            [ 21, 140, 200])
    lut.add(1172,  'ctx-lh-S_precentral-Inferior-part',               [ 21,  20, 240])
    lut.add(1173,  'ctx-lh-S_precentral-Superior-part',               [ 21,  20, 200])
    lut.add(1174,  'ctx-lh-S_subcentral_ant',                         [ 61, 180,  60])
    lut.add(1175,  'ctx-lh-S_subcentral_post',                        [ 61, 180, 250])
    lut.add(1176,  'ctx-lh-S_suborbital',                             [ 21,  20,  60])
    lut.add(1177,  'ctx-lh-S_subparietal',                            [101,  60,  60])
    lut.add(1178,  'ctx-lh-S_supracingulate',                         [ 21, 220, 220])
    lut.add(1179,  'ctx-lh-S_temporal_inferior',                      [ 21, 180, 180])
    lut.add(1180,  'ctx-lh-S_temporal_superior',                      [223, 220,  60])
    lut.add(1181,  'ctx-lh-S_temporal_transverse',                    [221,  60,  60])
    lut.add(1200,  'ctx-lh-G_cingulate-caudal_ACC',                   [ 25,  60,  61])
    lut.add(1201,  'ctx-lh-G_cingulate-rostral_ACC',                  [ 25,  90,  60])
    lut.add(1202,  'ctx-lh-G_cingulate-posterior',                    [ 25, 120,  60])
    lut.add(1205,  'ctx-lh-S_cingulate-caudal_ACC',                   [ 25, 150,  60])
    lut.add(1206,  'ctx-lh-S_cingulate-rostral_ACC',                  [ 25, 180,  60])
    lut.add(1207,  'ctx-lh-S_cingulate-posterior',                    [ 25, 210,  60])
    lut.add(1210,  'ctx-lh-S_pericallosal-caudal',                    [ 25, 150,  90])
    lut.add(1211,  'ctx-lh-S_pericallosal-rostral',                   [ 25, 180,  90])
    lut.add(1212,  'ctx-lh-S_pericallosal-posterior',                 [ 25, 210,  90])
    lut.add(1301,  'ctx-lh-frontal-lobe',                             [ 25, 100,  40])
    lut.add(1303,  'ctx-lh-cingulate-lobe',                           [100,  25,   0])
    lut.add(1304,  'ctx-lh-occipital-lobe',                           [120,  70,  50])
    lut.add(1305,  'ctx-lh-temporal-lobe',                            [220,  20, 100])
    lut.add(1306,  'ctx-lh-parietal-lobe',                            [220,  20,  10])
    lut.add(1307,  'ctx-lh-insula-lobe',                              [255, 192,  32])
    lut.add(2000,  'ctx-rh-unknown',                                  [ 25,   5,  25])
    lut.add(2001,  'ctx-rh-bankssts',                                 [ 25, 100,  40])
    lut.add(2002,  'ctx-rh-caudalanteriorcingulate',                  [125, 100, 160])
    lut.add(2003,  'ctx-rh-caudalmiddlefrontal',                      [100,  25,   0])
    lut.add(2004,  'ctx-rh-corpuscallosum',                           [120,  70,  50])
    lut.add(2005,  'ctx-rh-cuneus',                                   [220,  20, 100])
    lut.add(2006,  'ctx-rh-entorhinal',                               [220,  20,  10])
    lut.add(2007,  'ctx-rh-fusiform',                                 [180, 220, 140])
    lut.add(2008,  'ctx-rh-inferiorparietal',                         [220,  60, 220])
    lut.add(2009,  'ctx-rh-inferiortemporal',                         [180,  40, 120])
    lut.add(2010,  'ctx-rh-isthmuscingulate',                         [140,  20, 140])
    lut.add(2011,  'ctx-rh-lateraloccipital',                         [ 20,  30, 140])
    lut.add(2012,  'ctx-rh-lateralorbitofrontal',                     [ 35,  75,  50])
    lut.add(2013,  'ctx-rh-lingual',                                  [225, 140, 140])
    lut.add(2014,  'ctx-rh-medialorbitofrontal',                      [200,  35,  75])
    lut.add(2015,  'ctx-rh-middletemporal',                           [160, 100,  50])
    lut.add(2016,  'ctx-rh-parahippocampal',                          [ 20, 220,  60])
    lut.add(2017,  'ctx-rh-paracentral',                              [ 60, 220,  60])
    lut.add(2018,  'ctx-rh-parsopercularis',                          [220, 180, 140])
    lut.add(2019,  'ctx-rh-parsorbitalis',                            [ 20, 100,  50])
    lut.add(2020,  'ctx-rh-parstriangularis',                         [220,  60,  20])
    lut.add(2021,  'ctx-rh-pericalcarine',                            [120, 100,  60])
    lut.add(2022,  'ctx-rh-postcentral',                              [220,  20,  20])
    lut.add(2023,  'ctx-rh-posteriorcingulate',                       [220, 180, 220])
    lut.add(2024,  'ctx-rh-precentral',                               [ 60,  20, 220])
    lut.add(2025,  'ctx-rh-precuneus',                                [160, 140, 180])
    lut.add(2026,  'ctx-rh-rostralanteriorcingulate',                 [ 80,  20, 140])
    lut.add(2027,  'ctx-rh-rostralmiddlefrontal',                     [ 75,  50, 125])
    lut.add(2028,  'ctx-rh-superiorfrontal',                          [ 20, 220, 160])
    lut.add(2029,  'ctx-rh-superiorparietal',                         [ 20, 180, 140])
    lut.add(2030,  'ctx-rh-superiortemporal',                         [140, 220, 220])
    lut.add(2031,  'ctx-rh-supramarginal',                            [ 80, 160,  20])
    lut.add(2032,  'ctx-rh-frontalpole',                              [100,   0, 100])
    lut.add(2033,  'ctx-rh-temporalpole',                             [ 70,  70,  70])
    lut.add(2034,  'ctx-rh-transversetemporal',                       [150, 150, 200])
    lut.add(2035,  'ctx-rh-insula',                                   [255, 192,  32])
    lut.add(2100,  'ctx-rh-Unknown',                                  [  0,   0,   0])
    lut.add(2101,  'ctx-rh-Corpus_callosum',                          [ 50,  50,  50])
    lut.add(2102,  'ctx-rh-G_and_S_Insula_ONLY_AVERAGE',              [180,  20,  30])
    lut.add(2103,  'ctx-rh-G_cingulate-Isthmus',                      [ 60,  25,  25])
    lut.add(2104,  'ctx-rh-G_cingulate-Main_part',                    [ 25,  60,  60])
    lut.add(2105,  'ctx-rh-G_cuneus',                                 [180,  20,  20])
    lut.add(2106,  'ctx-rh-G_frontal_inf-Opercular_part',             [220,  20, 100])
    lut.add(2107,  'ctx-rh-G_frontal_inf-Orbital_part',               [140,  60,  60])
    lut.add(2108,  'ctx-rh-G_frontal_inf-Triangular_part',            [180, 220, 140])
    lut.add(2109,  'ctx-rh-G_frontal_middle',                         [140, 100, 180])
    lut.add(2110,  'ctx-rh-G_frontal_superior',                       [180,  20, 140])
    lut.add(2111,  'ctx-rh-G_frontomarginal',                         [140,  20, 140])
    lut.add(2112,  'ctx-rh-G_insular_long',                           [ 21,  10,  10])
    lut.add(2113,  'ctx-rh-G_insular_short',                          [225, 140, 140])
    lut.add(2114,  'ctx-rh-G_and_S_occipital_inferior',               [ 23,  60, 180])
    lut.add(2115,  'ctx-rh-G_occipital_middle',                       [180,  60, 180])
    lut.add(2116,  'ctx-rh-G_occipital_superior',                     [ 20, 220,  60])
    lut.add(2117,  'ctx-rh-G_occipit-temp_lat-Or_fusiform',           [ 60,  20, 140])
    lut.add(2118,  'ctx-rh-G_occipit-temp_med-Lingual_part',          [220, 180, 140])
    lut.add(2119,  'ctx-rh-G_occipit-temp_med-Parahippocampal_part',  [ 65, 100,  20])
    lut.add(2120,  'ctx-rh-G_orbital',                                [220,  60,  20])
    lut.add(2121,  'ctx-rh-G_paracentral',                            [ 60, 100,  60])
    lut.add(2122,  'ctx-rh-G_parietal_inferior-Angular_part',         [ 20,  60, 220])
    lut.add(2123,  'ctx-rh-G_parietal_inferior-Supramarginal_part',   [100, 100,  60])
    lut.add(2124,  'ctx-rh-G_parietal_superior',                      [220, 180, 220])
    lut.add(2125,  'ctx-rh-G_postcentral',                            [ 20, 180, 140])
    lut.add(2126,  'ctx-rh-G_precentral',                             [ 60, 140, 180])
    lut.add(2127,  'ctx-rh-G_precuneus',                              [ 25,  20, 140])
    lut.add(2128,  'ctx-rh-G_rectus',                                 [ 20,  60, 100])
    lut.add(2129,  'ctx-rh-G_subcallosal',                            [ 60, 220,  20])
    lut.add(2130,  'ctx-rh-G_subcentral',                             [ 60,  20, 220])
    lut.add(2131,  'ctx-rh-G_temporal_inferior',                      [220, 220, 100])
    lut.add(2132,  'ctx-rh-G_temporal_middle',                        [180,  60,  60])
    lut.add(2133,  'ctx-rh-G_temp_sup-G_temp_transv_and_interm_S',    [ 60,  60, 220])
    lut.add(2134,  'ctx-rh-G_temp_sup-Lateral_aspect',                [220,  60, 220])
    lut.add(2135,  'ctx-rh-G_temp_sup-Planum_polare',                 [ 65, 220,  60])
    lut.add(2136,  'ctx-rh-G_temp_sup-Planum_tempolare',              [ 25, 140,  20])
    lut.add(2137,  'ctx-rh-G_and_S_transverse_frontopolar',           [ 13,   0, 250])
    lut.add(2138,  'ctx-rh-Lat_Fissure-ant_sgt-ramus_horizontal',     [ 61,  20, 220])
    lut.add(2139,  'ctx-rh-Lat_Fissure-ant_sgt-ramus_vertical',       [ 61,  20,  60])
    lut.add(2140,  'ctx-rh-Lat_Fissure-post_sgt',                     [ 61,  60, 100])
    lut.add(2141,  'ctx-rh-Medial_wall',                              [ 25,  25,  25])
    lut.add(2142,  'ctx-rh-Pole_occipital',                           [140,  20,  60])
    lut.add(2143,  'ctx-rh-Pole_temporal',                            [220, 180,  20])
    lut.add(2144,  'ctx-rh-S_calcarine',                              [ 63, 180, 180])
    lut.add(2145,  'ctx-rh-S_central',                                [221,  20,  10])
    lut.add(2146,  'ctx-rh-S_central_insula',                         [ 21, 220,  20])
    lut.add(2147,  'ctx-rh-S_cingulate-Main_part_and_Intracingulate', [183, 100,  20])
    lut.add(2148,  'ctx-rh-S_cingulate-Marginalis_part',              [221,  20, 100])
    lut.add(2149,  'ctx-rh-S_circular_insula_anterior',               [221,  60, 140])
    lut.add(2150,  'ctx-rh-S_circular_insula_inferior',               [221,  20, 220])
    lut.add(2151,  'ctx-rh-S_circular_insula_superior',               [ 61, 220, 220])
    lut.add(2152,  'ctx-rh-S_collateral_transverse_ant',              [100, 200, 200])
    lut.add(2153,  'ctx-rh-S_collateral_transverse_post',             [ 10, 200, 200])
    lut.add(2154,  'ctx-rh-S_frontal_inferior',                       [221, 220,  20])
    lut.add(2155,  'ctx-rh-S_frontal_middle',                         [141,  20, 100])
    lut.add(2156,  'ctx-rh-S_frontal_superior',                       [ 61, 220, 100])
    lut.add(2157,  'ctx-rh-S_frontomarginal',                         [ 21, 220,  60])
    lut.add(2158,  'ctx-rh-S_intermedius_primus-Jensen',              [141,  60,  20])
    lut.add(2159,  'ctx-rh-S_intraparietal-and_Parietal_transverse',  [143,  20, 220])
    lut.add(2160,  'ctx-rh-S_occipital_anterior',                     [ 61,  20, 180])
    lut.add(2161,  'ctx-rh-S_occipital_middle_and_Lunatus',           [101,  60, 220])
    lut.add(2162,  'ctx-rh-S_occipital_superior_and_transversalis',   [ 21,  20, 140])
    lut.add(2163,  'ctx-rh-S_occipito-temporal_lateral',              [221, 140,  20])
    lut.add(2164,  'ctx-rh-S_occipito-temporal_medial_and_S_Lingual', [141, 100, 220])
    lut.add(2165,  'ctx-rh-S_orbital-H_shapped',                      [101,  20,  20])
    lut.add(2166,  'ctx-rh-S_orbital_lateral',                        [221, 100,  20])
    lut.add(2167,  'ctx-rh-S_orbital_medial-Or_olfactory',            [181, 200,  20])
    lut.add(2168,  'ctx-rh-S_paracentral',                            [ 21, 180, 140])
    lut.add(2169,  'ctx-rh-S_parieto_occipital',                      [101, 100, 180])
    lut.add(2170,  'ctx-rh-S_pericallosal',                           [181, 220,  20])
    lut.add(2171,  'ctx-rh-S_postcentral',                            [ 21, 140, 200])
    lut.add(2172,  'ctx-rh-S_precentral-Inferior-part',               [ 21,  20, 240])
    lut.add(2173,  'ctx-rh-S_precentral-Superior-part',               [ 21,  20, 200])
    lut.add(2174,  'ctx-rh-S_subcentral_ant',                         [ 61, 180,  60])
    lut.add(2175,  'ctx-rh-S_subcentral_post',                        [ 61, 180, 250])
    lut.add(2176,  'ctx-rh-S_suborbital',                             [ 21,  20,  60])
    lut.add(2177,  'ctx-rh-S_subparietal',                            [101,  60,  60])
    lut.add(2178,  'ctx-rh-S_supracingulate',                         [ 21, 220, 220])
    lut.add(2179,  'ctx-rh-S_temporal_inferior',                      [ 21, 180, 180])
    lut.add(2180,  'ctx-rh-S_temporal_superior',                      [223, 220,  60])
    lut.add(2181,  'ctx-rh-S_temporal_transverse',                    [221,  60,  60])
    lut.add(2200,  'ctx-rh-G_cingulate-caudal_ACC',                   [ 25,  60,  61])
    lut.add(2201,  'ctx-rh-G_cingulate-rostral_ACC',                  [ 25,  90,  60])
    lut.add(2202,  'ctx-rh-G_cingulate-posterior',                    [ 25, 120,  60])
    lut.add(2205,  'ctx-rh-S_cingulate-caudal_ACC',                   [ 25, 150,  60])
    lut.add(2206,  'ctx-rh-S_cingulate-rostral_ACC',                  [ 25, 180,  60])
    lut.add(2207,  'ctx-rh-S_cingulate-posterior',                    [ 25, 210,  60])
    lut.add(2210,  'ctx-rh-S_pericallosal-caudal',                    [ 25, 150,  90])
    lut.add(2211,  'ctx-rh-S_pericallosal-rostral',                   [ 25, 180,  90])
    lut.add(2212,  'ctx-rh-S_pericallosal-posterior',                 [ 25, 210,  90])
    lut.add(2301,  'ctx-rh-frontal-lobe',                             [ 25, 100,  40])
    lut.add(2303,  'ctx-rh-cingulate-lobe',                           [100,  25,   0])
    lut.add(2304,  'ctx-rh-occipital-lobe',                           [120,  70,  50])
    lut.add(2305,  'ctx-rh-temporal-lobe',                            [220,  20, 100])
    lut.add(2306,  'ctx-rh-parietal-lobe',                            [220,  20,  10])
    lut.add(2307,  'ctx-rh-insula-lobe',                              [255, 192,  32])

    # white matter labels
    lut.add(3000,  'wm-lh-unknown',                                   [230, 250, 230])
    lut.add(3001,  'wm-lh-bankssts',                                  [230, 155, 215])
    lut.add(3002,  'wm-lh-caudalanteriorcingulate',                   [130, 155,  95])
    lut.add(3003,  'wm-lh-caudalmiddlefrontal',                       [155, 230, 255])
    lut.add(3004,  'wm-lh-corpuscallosum',                            [135, 185, 205])
    lut.add(3005,  'wm-lh-cuneus',                                    [ 35, 235, 155])
    lut.add(3006,  'wm-lh-entorhinal',                                [ 35, 235, 245])
    lut.add(3007,  'wm-lh-fusiform',                                  [ 75,  35, 115])
    lut.add(3008,  'wm-lh-inferiorparietal',                          [ 35, 195,  35])
    lut.add(3009,  'wm-lh-inferiortemporal',                          [ 75, 215, 135])
    lut.add(3010,  'wm-lh-isthmuscingulate',                          [115, 235, 115])
    lut.add(3011,  'wm-lh-lateraloccipital',                          [235, 225, 115])
    lut.add(3012,  'wm-lh-lateralorbitofrontal',                      [220, 180, 205])
    lut.add(3013,  'wm-lh-lingual',                                   [ 30, 115, 115])
    lut.add(3014,  'wm-lh-medialorbitofrontal',                       [ 55, 220, 180])
    lut.add(3015,  'wm-lh-middletemporal',                            [ 95, 155, 205])
    lut.add(3016,  'wm-lh-parahippocampal',                           [235,  35, 195])
    lut.add(3017,  'wm-lh-paracentral',                               [195,  35, 195])
    lut.add(3018,  'wm-lh-parsopercularis',                           [ 35,  75, 115])
    lut.add(3019,  'wm-lh-parsorbitalis',                             [235, 155, 205])
    lut.add(3020,  'wm-lh-parstriangularis',                          [ 35, 195, 235])
    lut.add(3021,  'wm-lh-pericalcarine',                             [135, 155, 195])
    lut.add(3022,  'wm-lh-postcentral',                               [ 35, 235, 235])
    lut.add(3023,  'wm-lh-posteriorcingulate',                        [ 35,  75,  35])
    lut.add(3024,  'wm-lh-precentral',                                [195, 235,  35])
    lut.add(3025,  'wm-lh-precuneus',                                 [ 95, 115,  75])
    lut.add(3026,  'wm-lh-rostralanteriorcingulate',                  [175, 235, 115])
    lut.add(3027,  'wm-lh-rostralmiddlefrontal',                      [180, 205, 130])
    lut.add(3028,  'wm-lh-superiorfrontal',                           [235,  35,  95])
    lut.add(3029,  'wm-lh-superiorparietal',                          [235,  75, 115])
    lut.add(3030,  'wm-lh-superiortemporal',                          [115,  35,  35])
    lut.add(3031,  'wm-lh-supramarginal',                             [175,  95, 235])
    lut.add(3032,  'wm-lh-frontalpole',                               [155, 255, 155])
    lut.add(3033,  'wm-lh-temporalpole',                              [185, 185, 185])
    lut.add(3034,  'wm-lh-transversetemporal',                        [105, 105,  55])
    lut.add(3035,  'wm-lh-insula',                                    [ 20, 220, 160])
    lut.add(3100,  'wm-lh-Unknown',                                   [  0,   0,   0])
    lut.add(3101,  'wm-lh-Corpus_callosum',                           [ 50,  50,  50])
    lut.add(3102,  'wm-lh-G_and_S_Insula_ONLY_AVERAGE',               [180,  20,  30])
    lut.add(3103,  'wm-lh-G_cingulate-Isthmus',                       [ 60,  25,  25])
    lut.add(3104,  'wm-lh-G_cingulate-Main_part',                     [ 25,  60,  60])
    lut.add(3105,  'wm-lh-G_cuneus',                                  [180,  20,  20])
    lut.add(3106,  'wm-lh-G_frontal_inf-Opercular_part',              [220,  20, 100])
    lut.add(3107,  'wm-lh-G_frontal_inf-Orbital_part',                [140,  60,  60])
    lut.add(3108,  'wm-lh-G_frontal_inf-Triangular_part',             [180, 220, 140])
    lut.add(3109,  'wm-lh-G_frontal_middle',                          [140, 100, 180])
    lut.add(3110,  'wm-lh-G_frontal_superior',                        [180,  20, 140])
    lut.add(3111,  'wm-lh-G_frontomarginal',                          [140,  20, 140])
    lut.add(3112,  'wm-lh-G_insular_long',                            [ 21,  10,  10])
    lut.add(3113,  'wm-lh-G_insular_short',                           [225, 140, 140])
    lut.add(3114,  'wm-lh-G_and_S_occipital_inferior',                [ 23,  60, 180])
    lut.add(3115,  'wm-lh-G_occipital_middle',                        [180,  60, 180])
    lut.add(3116,  'wm-lh-G_occipital_superior',                      [ 20, 220,  60])
    lut.add(3117,  'wm-lh-G_occipit-temp_lat-Or_fusiform',            [ 60,  20, 140])
    lut.add(3118,  'wm-lh-G_occipit-temp_med-Lingual_part',           [220, 180, 140])
    lut.add(3119,  'wm-lh-G_occipit-temp_med-Parahippocampal_part',   [ 65, 100,  20])
    lut.add(3120,  'wm-lh-G_orbital',                                 [220,  60,  20])
    lut.add(3121,  'wm-lh-G_paracentral',                             [ 60, 100,  60])
    lut.add(3122,  'wm-lh-G_parietal_inferior-Angular_part',          [ 20,  60, 220])
    lut.add(3123,  'wm-lh-G_parietal_inferior-Supramarginal_part',    [100, 100,  60])
    lut.add(3124,  'wm-lh-G_parietal_superior',                       [220, 180, 220])
    lut.add(3125,  'wm-lh-G_postcentral',                             [ 20, 180, 140])
    lut.add(3126,  'wm-lh-G_precentral',                              [ 60, 140, 180])
    lut.add(3127,  'wm-lh-G_precuneus',                               [ 25,  20, 140])
    lut.add(3128,  'wm-lh-G_rectus',                                  [ 20,  60, 100])
    lut.add(3129,  'wm-lh-G_subcallosal',                             [ 60, 220,  20])
    lut.add(3130,  'wm-lh-G_subcentral',                              [ 60,  20, 220])
    lut.add(3131,  'wm-lh-G_temporal_inferior',                       [220, 220, 100])
    lut.add(3132,  'wm-lh-G_temporal_middle',                         [180,  60,  60])
    lut.add(3133,  'wm-lh-G_temp_sup-G_temp_transv_and_interm_S',     [ 60,  60, 220])
    lut.add(3134,  'wm-lh-G_temp_sup-Lateral_aspect',                 [220,  60, 220])
    lut.add(3135,  'wm-lh-G_temp_sup-Planum_polare',                  [ 65, 220,  60])
    lut.add(3136,  'wm-lh-G_temp_sup-Planum_tempolare',               [ 25, 140,  20])
    lut.add(3137,  'wm-lh-G_and_S_transverse_frontopolar',            [ 13,   0, 250])
    lut.add(3138,  'wm-lh-Lat_Fissure-ant_sgt-ramus_horizontal',      [ 61,  20, 220])
    lut.add(3139,  'wm-lh-Lat_Fissure-ant_sgt-ramus_vertical',        [ 61,  20,  60])
    lut.add(3140,  'wm-lh-Lat_Fissure-post_sgt',                      [ 61,  60, 100])
    lut.add(3141,  'wm-lh-Medial_wall',                               [ 25,  25,  25])
    lut.add(3142,  'wm-lh-Pole_occipital',                            [140,  20,  60])
    lut.add(3143,  'wm-lh-Pole_temporal',                             [220, 180,  20])
    lut.add(3144,  'wm-lh-S_calcarine',                               [ 63, 180, 180])
    lut.add(3145,  'wm-lh-S_central',                                 [221,  20,  10])
    lut.add(3146,  'wm-lh-S_central_insula',                          [ 21, 220,  20])
    lut.add(3147,  'wm-lh-S_cingulate-Main_part_and_Intracingulate',  [183, 100,  20])
    lut.add(3148,  'wm-lh-S_cingulate-Marginalis_part',               [221,  20, 100])
    lut.add(3149,  'wm-lh-S_circular_insula_anterior',                [221,  60, 140])
    lut.add(3150,  'wm-lh-S_circular_insula_inferior',                [221,  20, 220])
    lut.add(3151,  'wm-lh-S_circular_insula_superior',                [ 61, 220, 220])
    lut.add(3152,  'wm-lh-S_collateral_transverse_ant',               [100, 200, 200])
    lut.add(3153,  'wm-lh-S_collateral_transverse_post',              [ 10, 200, 200])
    lut.add(3154,  'wm-lh-S_frontal_inferior',                        [221, 220,  20])
    lut.add(3155,  'wm-lh-S_frontal_middle',                          [141,  20, 100])
    lut.add(3156,  'wm-lh-S_frontal_superior',                        [ 61, 220, 100])
    lut.add(3157,  'wm-lh-S_frontomarginal',                          [ 21, 220,  60])
    lut.add(3158,  'wm-lh-S_intermedius_primus-Jensen',               [141,  60,  20])
    lut.add(3159,  'wm-lh-S_intraparietal-and_Parietal_transverse',   [143,  20, 220])
    lut.add(3160,  'wm-lh-S_occipital_anterior',                      [ 61,  20, 180])
    lut.add(3161,  'wm-lh-S_occipital_middle_and_Lunatus',            [101,  60, 220])
    lut.add(3162,  'wm-lh-S_occipital_superior_and_transversalis',    [ 21,  20, 140])
    lut.add(3163,  'wm-lh-S_occipito-temporal_lateral',               [221, 140,  20])
    lut.add(3164,  'wm-lh-S_occipito-temporal_medial_and_S_Lingual',  [141, 100, 220])
    lut.add(3165,  'wm-lh-S_orbital-H_shapped',                       [101,  20,  20])
    lut.add(3166,  'wm-lh-S_orbital_lateral',                         [221, 100,  20])
    lut.add(3167,  'wm-lh-S_orbital_medial-Or_olfactory',             [181, 200,  20])
    lut.add(3168,  'wm-lh-S_paracentral',                             [ 21, 180, 140])
    lut.add(3169,  'wm-lh-S_parieto_occipital',                       [101, 100, 180])
    lut.add(3170,  'wm-lh-S_pericallosal',                            [181, 220,  20])
    lut.add(3171,  'wm-lh-S_postcentral',                             [ 21, 140, 200])
    lut.add(3172,  'wm-lh-S_precentral-Inferior-part',                [ 21,  20, 240])
    lut.add(3173,  'wm-lh-S_precentral-Superior-part',                [ 21,  20, 200])
    lut.add(3174,  'wm-lh-S_subcentral_ant',                          [ 61, 180,  60])
    lut.add(3175,  'wm-lh-S_subcentral_post',                         [ 61, 180, 250])
    lut.add(3176,  'wm-lh-S_suborbital',                              [ 21,  20,  60])
    lut.add(3177,  'wm-lh-S_subparietal',                             [101,  60,  60])
    lut.add(3178,  'wm-lh-S_supracingulate',                          [ 21, 220, 220])
    lut.add(3179,  'wm-lh-S_temporal_inferior',                       [ 21, 180, 180])
    lut.add(3180,  'wm-lh-S_temporal_superior',                       [223, 220,  60])
    lut.add(3181,  'wm-lh-S_temporal_transverse',                     [221,  60,  60])
    lut.add(3201,  'wm-lh-frontal-lobe',                              [235,  35,  95])
    lut.add(3203,  'wm-lh-cingulate-lobe',                            [ 35,  75,  35])
    lut.add(3204,  'wm-lh-occipital-lobe',                            [135, 155, 195])
    lut.add(3205,  'wm-lh-temporal-lobe',                             [115,  35,  35])
    lut.add(3206,  'wm-lh-parietal-lobe',                             [ 35, 195,  35])
    lut.add(3207,  'wm-lh-insula-lobe',                               [ 20, 220, 160])
    lut.add(4000,  'wm-rh-unknown',                                   [230, 250, 230])
    lut.add(4001,  'wm-rh-bankssts',                                  [230, 155, 215])
    lut.add(4002,  'wm-rh-caudalanteriorcingulate',                   [130, 155,  95])
    lut.add(4003,  'wm-rh-caudalmiddlefrontal',                       [155, 230, 255])
    lut.add(4004,  'wm-rh-corpuscallosum',                            [135, 185, 205])
    lut.add(4005,  'wm-rh-cuneus',                                    [ 35, 235, 155])
    lut.add(4006,  'wm-rh-entorhinal',                                [ 35, 235, 245])
    lut.add(4007,  'wm-rh-fusiform',                                  [ 75,  35, 115])
    lut.add(4008,  'wm-rh-inferiorparietal',                          [ 35, 195,  35])
    lut.add(4009,  'wm-rh-inferiortemporal',                          [ 75, 215, 135])
    lut.add(4010,  'wm-rh-isthmuscingulate',                          [115, 235, 115])
    lut.add(4011,  'wm-rh-lateraloccipital',                          [235, 225, 115])
    lut.add(4012,  'wm-rh-lateralorbitofrontal',                      [220, 180, 205])
    lut.add(4013,  'wm-rh-lingual',                                   [ 30, 115, 115])
    lut.add(4014,  'wm-rh-medialorbitofrontal',                       [ 55, 220, 180])
    lut.add(4015,  'wm-rh-middletemporal',                            [ 95, 155, 205])
    lut.add(4016,  'wm-rh-parahippocampal',                           [235,  35, 195])
    lut.add(4017,  'wm-rh-paracentral',                               [195,  35, 195])
    lut.add(4018,  'wm-rh-parsopercularis',                           [ 35,  75, 115])
    lut.add(4019,  'wm-rh-parsorbitalis',                             [235, 155, 205])
    lut.add(4020,  'wm-rh-parstriangularis',                          [ 35, 195, 235])
    lut.add(4021,  'wm-rh-pericalcarine',                             [135, 155, 195])
    lut.add(4022,  'wm-rh-postcentral',                               [ 35, 235, 235])
    lut.add(4023,  'wm-rh-posteriorcingulate',                        [ 35,  75,  35])
    lut.add(4024,  'wm-rh-precentral',                                [195, 235,  35])
    lut.add(4025,  'wm-rh-precuneus',                                 [ 95, 115,  75])
    lut.add(4026,  'wm-rh-rostralanteriorcingulate',                  [175, 235, 115])
    lut.add(4027,  'wm-rh-rostralmiddlefrontal',                      [180, 205, 130])
    lut.add(4028,  'wm-rh-superiorfrontal',                           [235,  35,  95])
    lut.add(4029,  'wm-rh-superiorparietal',                          [235,  75, 115])
    lut.add(4030,  'wm-rh-superiortemporal',                          [115,  35,  35])
    lut.add(4031,  'wm-rh-supramarginal',                             [175,  95, 235])
    lut.add(4032,  'wm-rh-frontalpole',                               [155, 255, 155])
    lut.add(4033,  'wm-rh-temporalpole',                              [185, 185, 185])
    lut.add(4034,  'wm-rh-transversetemporal',                        [105, 105,  55])
    lut.add(4035,  'wm-rh-insula',                                    [ 20, 220, 160])
    lut.add(4100,  'wm-rh-Unknown',                                   [  0,   0,   0])
    lut.add(4101,  'wm-rh-Corpus_callosum',                           [ 50,  50,  50])
    lut.add(4102,  'wm-rh-G_and_S_Insula_ONLY_AVERAGE',               [180,  20,  30])
    lut.add(4103,  'wm-rh-G_cingulate-Isthmus',                       [ 60,  25,  25])
    lut.add(4104,  'wm-rh-G_cingulate-Main_part',                     [ 25,  60,  60])
    lut.add(4105,  'wm-rh-G_cuneus',                                  [180,  20,  20])
    lut.add(4106,  'wm-rh-G_frontal_inf-Opercular_part',              [220,  20, 100])
    lut.add(4107,  'wm-rh-G_frontal_inf-Orbital_part',                [140,  60,  60])
    lut.add(4108,  'wm-rh-G_frontal_inf-Triangular_part',             [180, 220, 140])
    lut.add(4109,  'wm-rh-G_frontal_middle',                          [140, 100, 180])
    lut.add(4110,  'wm-rh-G_frontal_superior',                        [180,  20, 140])
    lut.add(4111,  'wm-rh-G_frontomarginal',                          [140,  20, 140])
    lut.add(4112,  'wm-rh-G_insular_long',                            [ 21,  10,  10])
    lut.add(4113,  'wm-rh-G_insular_short',                           [225, 140, 140])
    lut.add(4114,  'wm-rh-G_and_S_occipital_inferior',                [ 23,  60, 180])
    lut.add(4115,  'wm-rh-G_occipital_middle',                        [180,  60, 180])
    lut.add(4116,  'wm-rh-G_occipital_superior',                      [ 20, 220,  60])
    lut.add(4117,  'wm-rh-G_occipit-temp_lat-Or_fusiform',            [ 60,  20, 140])
    lut.add(4118,  'wm-rh-G_occipit-temp_med-Lingual_part',           [220, 180, 140])
    lut.add(4119,  'wm-rh-G_occipit-temp_med-Parahippocampal_part',   [ 65, 100,  20])
    lut.add(4120,  'wm-rh-G_orbital',                                 [220,  60,  20])
    lut.add(4121,  'wm-rh-G_paracentral',                             [ 60, 100,  60])
    lut.add(4122,  'wm-rh-G_parietal_inferior-Angular_part',          [ 20,  60, 220])
    lut.add(4123,  'wm-rh-G_parietal_inferior-Supramarginal_part',    [100, 100,  60])
    lut.add(4124,  'wm-rh-G_parietal_superior',                       [220, 180, 220])
    lut.add(4125,  'wm-rh-G_postcentral',                             [ 20, 180, 140])
    lut.add(4126,  'wm-rh-G_precentral',                              [ 60, 140, 180])
    lut.add(4127,  'wm-rh-G_precuneus',                               [ 25,  20, 140])
    lut.add(4128,  'wm-rh-G_rectus',                                  [ 20,  60, 100])
    lut.add(4129,  'wm-rh-G_subcallosal',                             [ 60, 220,  20])
    lut.add(4130,  'wm-rh-G_subcentral',                              [ 60,  20, 220])
    lut.add(4131,  'wm-rh-G_temporal_inferior',                       [220, 220, 100])
    lut.add(4132,  'wm-rh-G_temporal_middle',                         [180,  60,  60])
    lut.add(4133,  'wm-rh-G_temp_sup-G_temp_transv_and_interm_S',     [ 60,  60, 220])
    lut.add(4134,  'wm-rh-G_temp_sup-Lateral_aspect',                 [220,  60, 220])
    lut.add(4135,  'wm-rh-G_temp_sup-Planum_polare',                  [ 65, 220,  60])
    lut.add(4136,  'wm-rh-G_temp_sup-Planum_tempolare',               [ 25, 140,  20])
    lut.add(4137,  'wm-rh-G_and_S_transverse_frontopolar',            [ 13,   0, 250])
    lut.add(4138,  'wm-rh-Lat_Fissure-ant_sgt-ramus_horizontal',      [ 61,  20, 220])
    lut.add(4139,  'wm-rh-Lat_Fissure-ant_sgt-ramus_vertical',        [ 61,  20,  60])
    lut.add(4140,  'wm-rh-Lat_Fissure-post_sgt',                      [ 61,  60, 100])
    lut.add(4141,  'wm-rh-Medial_wall',                               [ 25,  25,  25])
    lut.add(4142,  'wm-rh-Pole_occipital',                            [140,  20,  60])
    lut.add(4143,  'wm-rh-Pole_temporal',                             [220, 180,  20])
    lut.add(4144,  'wm-rh-S_calcarine',                               [ 63, 180, 180])
    lut.add(4145,  'wm-rh-S_central',                                 [221,  20,  10])
    lut.add(4146,  'wm-rh-S_central_insula',                          [ 21, 220,  20])
    lut.add(4147,  'wm-rh-S_cingulate-Main_part_and_Intracingulate',  [183, 100,  20])
    lut.add(4148,  'wm-rh-S_cingulate-Marginalis_part',               [221,  20, 100])
    lut.add(4149,  'wm-rh-S_circular_insula_anterior',                [221,  60, 140])
    lut.add(4150,  'wm-rh-S_circular_insula_inferior',                [221,  20, 220])
    lut.add(4151,  'wm-rh-S_circular_insula_superior',                [ 61, 220, 220])
    lut.add(4152,  'wm-rh-S_collateral_transverse_ant',               [100, 200, 200])
    lut.add(4153,  'wm-rh-S_collateral_transverse_post',              [ 10, 200, 200])
    lut.add(4154,  'wm-rh-S_frontal_inferior',                        [221, 220,  20])
    lut.add(4155,  'wm-rh-S_frontal_middle',                          [141,  20, 100])
    lut.add(4156,  'wm-rh-S_frontal_superior',                        [ 61, 220, 100])
    lut.add(4157,  'wm-rh-S_frontomarginal',                          [ 21, 220,  60])
    lut.add(4158,  'wm-rh-S_intermedius_primus-Jensen',               [141,  60,  20])
    lut.add(4159,  'wm-rh-S_intraparietal-and_Parietal_transverse',   [143,  20, 220])
    lut.add(4160,  'wm-rh-S_occipital_anterior',                      [ 61,  20, 180])
    lut.add(4161,  'wm-rh-S_occipital_middle_and_Lunatus',            [101,  60, 220])
    lut.add(4162,  'wm-rh-S_occipital_superior_and_transversalis',    [ 21,  20, 140])
    lut.add(4163,  'wm-rh-S_occipito-temporal_lateral',               [221, 140,  20])
    lut.add(4164,  'wm-rh-S_occipito-temporal_medial_and_S_Lingual',  [141, 100, 220])
    lut.add(4165,  'wm-rh-S_orbital-H_shapped',                       [101,  20,  20])
    lut.add(4166,  'wm-rh-S_orbital_lateral',                         [221, 100,  20])
    lut.add(4167,  'wm-rh-S_orbital_medial-Or_olfactory',             [181, 200,  20])
    lut.add(4168,  'wm-rh-S_paracentral',                             [ 21, 180, 140])
    lut.add(4169,  'wm-rh-S_parieto_occipital',                       [101, 100, 180])
    lut.add(4170,  'wm-rh-S_pericallosal',                            [181, 220,  20])
    lut.add(4171,  'wm-rh-S_postcentral',                             [ 21, 140, 200])
    lut.add(4172,  'wm-rh-S_precentral-Inferior-part',                [ 21,  20, 240])
    lut.add(4173,  'wm-rh-S_precentral-Superior-part',                [ 21,  20, 200])
    lut.add(4174,  'wm-rh-S_subcentral_ant',                          [ 61, 180,  60])
    lut.add(4175,  'wm-rh-S_subcentral_post',                         [ 61, 180, 250])
    lut.add(4176,  'wm-rh-S_suborbital',                              [ 21,  20,  60])
    lut.add(4177,  'wm-rh-S_subparietal',                             [101,  60,  60])
    lut.add(4178,  'wm-rh-S_supracingulate',                          [ 21, 220, 220])
    lut.add(4179,  'wm-rh-S_temporal_inferior',                       [ 21, 180, 180])
    lut.add(4180,  'wm-rh-S_temporal_superior',                       [223, 220,  60])
    lut.add(4181,  'wm-rh-S_temporal_transverse',                     [221,  60,  60])
    lut.add(4201,  'wm-rh-frontal-lobe',                              [235,  35,  95])
    lut.add(4203,  'wm-rh-cingulate-lobe',                            [ 35,  75,  35])
    lut.add(4204,  'wm-rh-occipital-lobe',                            [135, 155, 195])
    lut.add(4205,  'wm-rh-temporal-lobe',                             [115,  35,  35])
    lut.add(4206,  'wm-rh-parietal-lobe',                             [ 35, 195,  35])
    lut.add(4207,  'wm-rh-insula-lobe',                               [ 20, 220, 160])

    lut.add(5001,  'Left-UnsegmentedWhiteMatter',                     [ 20,  30,  40])
    lut.add(5002,  'Right-UnsegmentedWhiteMatter',                    [ 20,  30,  40])

    lut.add(5100,  'fmajor',                                          [204, 102, 102])
    lut.add(5101,  'fminor',                                          [204, 102, 102])
    lut.add(5102,  'lh.atr',                                          [255, 255, 102])
    lut.add(5103,  'lh.cab',                                          [153, 204,   0])
    lut.add(5104,  'lh.ccg',                                          [  0, 153, 153])
    lut.add(5105,  'lh.cst',                                          [204, 153, 255])
    lut.add(5106,  'lh.ilf',                                          [255, 153,  51])
    lut.add(5107,  'lh.slfp',                                         [204, 204, 204])
    lut.add(5108,  'lh.slft',                                         [153, 255, 255])
    lut.add(5109,  'lh.unc',                                          [102, 153, 255])
    lut.add(5110,  'rh.atr',                                          [255, 255, 102])
    lut.add(5111,  'rh.cab',                                          [153, 204,   0])
    lut.add(5112,  'rh.ccg',                                          [  0, 153, 153])
    lut.add(5113,  'rh.cst',                                          [204, 153, 255])
    lut.add(5114,  'rh.ilf',                                          [255, 153,  51])
    lut.add(5115,  'rh.slfp',                                         [204, 204, 204])
    lut.add(5116,  'rh.slft',                                         [153, 255, 255])
    lut.add(5117,  'rh.unc',                                          [102, 153, 255])
    lut.add(5200,  'CC-ForcepsMajor',                                 [204, 102, 102])
    lut.add(5201,  'CC-ForcepsMinor',                                 [204, 102, 102])
    lut.add(5202,  'LAntThalRadiation',                               [255, 255, 102])
    lut.add(5203,  'LCingulumAngBundle',                              [153, 204,   0])
    lut.add(5204,  'LCingulumCingGyrus',                              [  0, 153, 153])
    lut.add(5205,  'LCorticospinalTract',                             [204, 153, 255])
    lut.add(5206,  'LInfLongFas',                                     [255, 153,  51])
    lut.add(5207,  'LSupLongFasParietal',                             [204, 204, 204])
    lut.add(5208,  'LSupLongFasTemporal',                             [153, 255, 255])
    lut.add(5209,  'LUncinateFas',                                    [102, 153, 255])
    lut.add(5210,  'RAntThalRadiation',                               [255, 255, 102])
    lut.add(5211,  'RCingulumAngBundle',                              [153, 204,   0])
    lut.add(5212,  'RCingulumCingGyrus',                              [  0, 153, 153])
    lut.add(5213,  'RCorticospinalTract',                             [204, 153, 255])
    lut.add(5214,  'RInfLongFas',                                     [255, 153,  51])
    lut.add(5215,  'RSupLongFasParietal',                             [204, 204, 204])
    lut.add(5216,  'RSupLongFasTemporal',                             [153, 255, 255])
    lut.add(5217,  'RUncinateFas',                                    [102, 153, 255])
    lut.add(6000,  'CST-orig',                                        [  0, 255,   0])
    lut.add(6001,  'CST-hammer',                                      [255, 255,   0])
    lut.add(6002,  'CST-CVS',                                         [  0, 255, 255])
    lut.add(6003,  'CST-flirt',                                       [  0,   0, 255])
    lut.add(6010,  'Left-SLF1',                                       [236,  16, 231])
    lut.add(6020,  'Right-SLF1',                                      [237,  18, 232])
    lut.add(6030,  'Left-SLF3',                                       [236,  13, 227])
    lut.add(6040,  'Right-SLF3',                                      [236,  17, 228])
    lut.add(6050,  'Left-CST',                                        [  1, 255,   1])
    lut.add(6060,  'Right-CST',                                       [  2, 255,   1])
    lut.add(6070,  'Left-SLF2',                                       [236,  14, 230])
    lut.add(6080,  'Right-SLF2',                                      [237,  14, 230])
    lut.add(7001,  'Lateral-nucleus',                                 [ 72, 132, 181])
    lut.add(7002,  'Basolateral-nucleus',                             [243, 243, 243])
    lut.add(7003,  'Basal-nucleus',                                   [207,  63,  79])
    lut.add(7004,  'Centromedial-nucleus',                            [121,  20, 135])
    lut.add(7005,  'Central-nucleus',                                 [197,  60, 248])
    lut.add(7006,  'Medial-nucleus',                                  [  2, 149,   2])
    lut.add(7007,  'Cortical-nucleus',                                [221, 249, 166])
    lut.add(7008,  'Accessory-Basal-nucleus',                         [232, 146,  35])
    lut.add(7009,  'Corticoamygdaloid-transitio',                     [ 20,  60, 120])
    lut.add(7010,  'Anterior-amygdaloid-area-AAA',                    [250, 250,   0])
    lut.add(7011,  'Fusion-amygdala-HP-FAH',                          [122, 187, 222])
    lut.add(7012,  'Hippocampal-amygdala-transition-HATA',            [237,  12, 177])
    lut.add(7013,  'Endopiriform-nucleus',                            [ 10,  49, 255])
    lut.add(7014,  'Lateral-nucleus-olfactory-tract',                 [205, 184, 144])
    lut.add(7015,  'Paralaminar-nucleus',                             [ 45, 205, 165])
    lut.add(7016,  'Intercalated-nucleus',                            [117, 160, 175])
    lut.add(7017,  'Prepiriform-cortex',                              [221, 217,  21])
    lut.add(7018,  'Periamygdaloid-cortex',                           [ 20,  60, 120])
    lut.add(7019,  'Envelope-Amygdala',                               [141,  21, 100])
    lut.add(7020,  'Extranuclear-Amydala',                            [225, 140, 141])
    lut.add(7100,  'Brainstem-inferior-colliculus',                   [ 42, 201, 168])
    lut.add(7101,  'Brainstem-cochlear-nucleus',                      [168, 104, 162])
    lut.add(7201,  'DR',                                              [121, 255, 250])
    lut.add(7202,  'MnR',                                             [  0, 255,   0])
    lut.add(7203,  'PAG',                                             [153, 153, 255])
    lut.add(7204,  'VTA',                                             [255,   0, 255])
    lut.add(7301,  'Left-LC',                                         [  0,   0, 255])
    lut.add(7302,  'Left-LDTg',                                       [255, 127,   0])
    lut.add(7303,  'Left-mRt',                                        [255,   0,   0])
    lut.add(7304,  'Left-PBC',                                        [255, 255,   0])
    lut.add(7305,  'Left-PnO',                                        [  0, 127, 255])
    lut.add(7306,  'Left-PTg',                                        [127,   0, 255])
    lut.add(7401,  'Right-LC',                                        [  0,   0, 255])
    lut.add(7402,  'Right-LDTg',                                      [255, 127,   0])
    lut.add(7403,  'Right-mRt',                                       [255,   0,   0])
    lut.add(7404,  'Right-PBC',                                       [255, 255,   0])
    lut.add(7405,  'Right-PnO',                                       [  0, 127, 255])
    lut.add(7406,  'Right-PTg',                                       [127,   0, 255])

    lut.add(8001,  'Thalamus-Anterior',                               [ 74, 130, 181])
    lut.add(8002,  'Thalamus-Ventral-anterior',                       [242, 241, 240])
    lut.add(8003,  'Thalamus-Lateral-dorsal',                         [206,  65,  78])
    lut.add(8004,  'Thalamus-Lateral-posterior',                      [120,  21, 133])
    lut.add(8005,  'Thalamus-Ventral-lateral',                        [195,  61, 246])
    lut.add(8006,  'Thalamus-Ventral-posterior-medial',               [  3, 147,   6])
    lut.add(8007,  'Thalamus-Ventral-posterior-lateral',              [220, 251, 163])
    lut.add(8008,  'Thalamus-intralaminar',                           [232, 146,  33])
    lut.add(8009,  'Thalamus-centromedian',                           [  4, 114,  14])
    lut.add(8010,  'Thalamus-mediodorsal',                            [121, 184, 220])
    lut.add(8011,  'Thalamus-medial',                                 [235,  11, 175])
    lut.add(8012,  'Thalamus-pulvinar',                               [ 12,  46, 250])
    lut.add(8013,  'Thalamus-lateral-geniculate',                     [203, 182, 143])
    lut.add(8014,  'Thalamus-medial-geniculate',                      [ 42, 204, 167])

    lut.add(8103,  'Left-AV',                                         [  0,  85,   0])
    lut.add(8104,  'Left-CeM',                                        [170,  85,   0])
    lut.add(8105,  'Left-CL',                                         [  0, 170,   0])
    lut.add(8106,  'Left-CM',                                         [170, 170,   0])
    lut.add(8108,  'Left-LD',                                         [170, 255,   0])
    lut.add(8109,  'Left-LGN',                                        [  0,   0, 127])
    lut.add(8110,  'Left-LP',                                         [  0,  85, 127])
    lut.add(8111,  'Left-L-Sg',                                       [170,  85, 127])
    lut.add(8112,  'Left-MDl',                                        [  0, 170, 127])
    lut.add(8113,  'Left-MDm',                                        [170, 170, 127])
    lut.add(8115,  'Left-MGN',                                        [170, 255, 127])
    lut.add(8116,  'Left-MV(Re)',                                     [  0,   0, 255])
    lut.add(8117,  'Left-Pc',                                         [170,   0, 255])
    lut.add(8118,  'Left-Pf',                                         [  0,  85, 255])
    lut.add(8119,  'Left-Pt',                                         [170,  85, 255])
    lut.add(8120,  'Left-PuA',                                        [  0, 170, 255])
    lut.add(8121,  'Left-PuI',                                        [170, 170, 255])
    lut.add(8122,  'Left-PuL',                                        [  0, 255, 255])
    lut.add(8123,  'Left-PuM',                                        [170, 255, 255])
    lut.add(8125,  'Left-R',                                          [255,   0,   0])
    lut.add(8126,  'Left-VA',                                         [ 85,  85,   0])
    lut.add(8127,  'Left-VAmc',                                       [255,  85,   0])
    lut.add(8128,  'Left-VLa',                                        [ 85, 170,   0])
    lut.add(8129,  'Left-VLp',                                        [255, 170,   0])
    lut.add(8130,  'Left-VM',                                         [ 85, 255,   0])
    lut.add(8133,  'Left-VPL',                                        [255,   0, 255])
    lut.add(8134,  'Left-PaV',                                        [120,  18, 134])
    lut.add(8203,  'Right-AV',                                        [  0,  85,   0])
    lut.add(8204,  'Right-CeM',                                       [170,  85,   0])
    lut.add(8205,  'Right-CL',                                        [  0, 170,   0])
    lut.add(8206,  'Right-CM',                                        [170, 170,   0])
    lut.add(8208,  'Right-LD',                                        [170, 255,   0])
    lut.add(8209,  'Right-LGN',                                       [  0,   0, 127])
    lut.add(8210,  'Right-LP',                                        [  0,  85, 127])
    lut.add(8211,  'Right-L-Sg',                                      [170,  85, 127])
    lut.add(8212,  'Right-MDl',                                       [  0, 170, 127])
    lut.add(8213,  'Right-MDm',                                       [170, 170, 127])
    lut.add(8215,  'Right-MGN',                                       [170, 255, 127])
    lut.add(8216,  'Right-MV(Re)',                                    [  0,   0, 255])
    lut.add(8217,  'Right-Pc',                                        [170,   0, 255])
    lut.add(8218,  'Right-Pf',                                        [  0,  85, 255])
    lut.add(8219,  'Right-Pt',                                        [170,  85, 255])
    lut.add(8220,  'Right-PuA',                                       [  0, 170, 255])
    lut.add(8221,  'Right-PuI',                                       [170, 170, 255])
    lut.add(8222,  'Right-PuL',                                       [  0, 255, 255])
    lut.add(8223,  'Right-PuM',                                       [170, 255, 255])
    lut.add(8225,  'Right-R',                                         [255,   0,   0])
    lut.add(8226,  'Right-VA',                                        [ 85,  85,   0])
    lut.add(8227,  'Right-VAmc',                                      [255,  85,   0])
    lut.add(8228,  'Right-VLa',                                       [ 85, 170,   0])
    lut.add(8229,  'Right-VLp',                                       [255, 170,   0])
    lut.add(8230,  'Right-VM',                                        [ 85, 255,   0])
    lut.add(8233,  'Right-VPL',                                       [255,   0, 255])
    lut.add(8234,  'Right-PaV',                                       [120,  18, 134])

    lut.add(9000,  'ctx-lh-prefrontal',                               [ 50, 100,  30])
    lut.add(9001,  'ctx-lh-primary-motor',                            [ 30, 100,  45])
    lut.add(9002,  'ctx-lh-premotor',                                 [130, 100, 165])
    lut.add(9003,  'ctx-lh-temporal',                                 [105,  25,   5])
    lut.add(9004,  'ctx-lh-posterior-parietal',                       [125,  70,  55])
    lut.add(9005,  'ctx-lh-prim-sec-somatosensory',                   [225,  20, 105])
    lut.add(9006,  'ctx-lh-occipital',                                [225,  20,  15])
    lut.add(9500,  'ctx-rh-prefrontal',                               [ 50, 200,  30])
    lut.add(9501,  'ctx-rh-primary-motor',                            [ 30, 150,  45])
    lut.add(9502,  'ctx-rh-premotor',                                 [130, 150, 165])
    lut.add(9503,  'ctx-rh-temporal',                                 [105,  75,   5])
    lut.add(9504,  'ctx-rh-posterior-parietal',                       [125, 120,  55])
    lut.add(9505,  'ctx-rh-prim-sec-somatosensory',                   [225,  70, 105])
    lut.add(9506,  'ctx-rh-occipital',                                [225,  70,  15])
    lut.add(11100, 'ctx_lh_Unknown',                                  [  0,   0,   0])
    lut.add(11101, 'ctx_lh_G_and_S_frontomargin',                     [ 23, 220,  60])
    lut.add(11102, 'ctx_lh_G_and_S_occipital_inf',                    [ 23,  60, 180])
    lut.add(11103, 'ctx_lh_G_and_S_paracentral',                      [ 63, 100,  60])
    lut.add(11104, 'ctx_lh_G_and_S_subcentral',                       [ 63,  20, 220])
    lut.add(11105, 'ctx_lh_G_and_S_transv_frontopol',                 [ 13,   0, 250])
    lut.add(11106, 'ctx_lh_G_and_S_cingul-Ant',                       [ 26,  60,   0])
    lut.add(11107, 'ctx_lh_G_and_S_cingul-Mid-Ant',                   [ 26,  60,  75])
    lut.add(11108, 'ctx_lh_G_and_S_cingul-Mid-Post',                  [ 26,  60, 150])
    lut.add(11109, 'ctx_lh_G_cingul-Post-dorsal',                     [ 25,  60, 250])
    lut.add(11110, 'ctx_lh_G_cingul-Post-ventral',                    [ 60,  25,  25])
    lut.add(11111, 'ctx_lh_G_cuneus',                                 [180,  20,  20])
    lut.add(11112, 'ctx_lh_G_front_inf-Opercular',                    [220,  20, 100])
    lut.add(11113, 'ctx_lh_G_front_inf-Orbital',                      [140,  60,  60])
    lut.add(11114, 'ctx_lh_G_front_inf-Triangul',                     [180, 220, 140])
    lut.add(11115, 'ctx_lh_G_front_middle',                           [140, 100, 180])
    lut.add(11116, 'ctx_lh_G_front_sup',                              [180,  20, 140])
    lut.add(11117, 'ctx_lh_G_Ins_lg_and_S_cent_ins',                  [ 23,  10,  10])
    lut.add(11118, 'ctx_lh_G_insular_short',                          [225, 140, 140])
    lut.add(11119, 'ctx_lh_G_occipital_middle',                       [180,  60, 180])
    lut.add(11120, 'ctx_lh_G_occipital_sup',                          [ 20, 220,  60])
    lut.add(11121, 'ctx_lh_G_oc-temp_lat-fusifor',                    [ 60,  20, 140])
    lut.add(11122, 'ctx_lh_G_oc-temp_med-Lingual',                    [220, 180, 140])
    lut.add(11123, 'ctx_lh_G_oc-temp_med-Parahip',                    [ 65, 100,  20])
    lut.add(11124, 'ctx_lh_G_orbital',                                [220,  60,  20])
    lut.add(11125, 'ctx_lh_G_pariet_inf-Angular',                     [ 20,  60, 220])
    lut.add(11126, 'ctx_lh_G_pariet_inf-Supramar',                    [100, 100,  60])
    lut.add(11127, 'ctx_lh_G_parietal_sup',                           [220, 180, 220])
    lut.add(11128, 'ctx_lh_G_postcentral',                            [ 20, 180, 140])
    lut.add(11129, 'ctx_lh_G_precentral',                             [ 60, 140, 180])
    lut.add(11130, 'ctx_lh_G_precuneus',                              [ 25,  20, 140])
    lut.add(11131, 'ctx_lh_G_rectus',                                 [ 20,  60, 100])
    lut.add(11132, 'ctx_lh_G_subcallosal',                            [ 60, 220,  20])
    lut.add(11133, 'ctx_lh_G_temp_sup-G_T_transv',                    [ 60,  60, 220])
    lut.add(11134, 'ctx_lh_G_temp_sup-Lateral',                       [220,  60, 220])
    lut.add(11135, 'ctx_lh_G_temp_sup-Plan_polar',                    [ 65, 220,  60])
    lut.add(11136, 'ctx_lh_G_temp_sup-Plan_tempo',                    [ 25, 140,  20])
    lut.add(11137, 'ctx_lh_G_temporal_inf',                           [220, 220, 100])
    lut.add(11138, 'ctx_lh_G_temporal_middle',                        [180,  60,  60])
    lut.add(11139, 'ctx_lh_Lat_Fis-ant-Horizont',                     [ 61,  20, 220])
    lut.add(11140, 'ctx_lh_Lat_Fis-ant-Vertical',                     [ 61,  20,  60])
    lut.add(11141, 'ctx_lh_Lat_Fis-post',                             [ 61,  60, 100])
    lut.add(11142, 'ctx_lh_Medial_wall',                              [ 25,  25,  25])
    lut.add(11143, 'ctx_lh_Pole_occipital',                           [140,  20,  60])
    lut.add(11144, 'ctx_lh_Pole_temporal',                            [220, 180,  20])
    lut.add(11145, 'ctx_lh_S_calcarine',                              [ 63, 180, 180])
    lut.add(11146, 'ctx_lh_S_central',                                [221,  20,  10])
    lut.add(11147, 'ctx_lh_S_cingul-Marginalis',                      [221,  20, 100])
    lut.add(11148, 'ctx_lh_S_circular_insula_ant',                    [221,  60, 140])
    lut.add(11149, 'ctx_lh_S_circular_insula_inf',                    [221,  20, 220])
    lut.add(11150, 'ctx_lh_S_circular_insula_sup',                    [ 61, 220, 220])
    lut.add(11151, 'ctx_lh_S_collat_transv_ant',                      [100, 200, 200])
    lut.add(11152, 'ctx_lh_S_collat_transv_post',                     [ 10, 200, 200])
    lut.add(11153, 'ctx_lh_S_front_inf',                              [221, 220,  20])
    lut.add(11154, 'ctx_lh_S_front_middle',                           [141,  20, 100])
    lut.add(11155, 'ctx_lh_S_front_sup',                              [ 61, 220, 100])
    lut.add(11156, 'ctx_lh_S_interm_prim-Jensen',                     [141,  60,  20])
    lut.add(11157, 'ctx_lh_S_intrapariet_and_P_trans',                [143,  20, 220])
    lut.add(11158, 'ctx_lh_S_oc_middle_and_Lunatus',                  [101,  60, 220])
    lut.add(11159, 'ctx_lh_S_oc_sup_and_transversal',                 [ 21,  20, 140])
    lut.add(11160, 'ctx_lh_S_occipital_ant',                          [ 61,  20, 180])
    lut.add(11161, 'ctx_lh_S_oc-temp_lat',                            [221, 140,  20])
    lut.add(11162, 'ctx_lh_S_oc-temp_med_and_Lingual',                [141, 100, 220])
    lut.add(11163, 'ctx_lh_S_orbital_lateral',                        [221, 100,  20])
    lut.add(11164, 'ctx_lh_S_orbital_med-olfact',                     [181, 200,  20])
    lut.add(11165, 'ctx_lh_S_orbital-H_Shaped',                       [101,  20,  20])
    lut.add(11166, 'ctx_lh_S_parieto_occipital',                      [101, 100, 180])
    lut.add(11167, 'ctx_lh_S_pericallosal',                           [181, 220,  20])
    lut.add(11168, 'ctx_lh_S_postcentral',                            [ 21, 140, 200])
    lut.add(11169, 'ctx_lh_S_precentral-inf-part',                    [ 21,  20, 240])
    lut.add(11170, 'ctx_lh_S_precentral-sup-part',                    [ 21,  20, 200])
    lut.add(11171, 'ctx_lh_S_suborbital',                             [ 21,  20,  60])
    lut.add(11172, 'ctx_lh_S_subparietal',                            [101,  60,  60])
    lut.add(11173, 'ctx_lh_S_temporal_inf',                           [ 21, 180, 180])
    lut.add(11174, 'ctx_lh_S_temporal_sup',                           [223, 220,  60])
    lut.add(11175, 'ctx_lh_S_temporal_transverse',                    [221,  60,  60])
    lut.add(12100, 'ctx_rh_Unknown',                                  [  0,   0,   0])
    lut.add(12101, 'ctx_rh_G_and_S_frontomargin',                     [ 23, 220,  60])
    lut.add(12102, 'ctx_rh_G_and_S_occipital_inf',                    [ 23,  60, 180])
    lut.add(12103, 'ctx_rh_G_and_S_paracentral',                      [ 63, 100,  60])
    lut.add(12104, 'ctx_rh_G_and_S_subcentral',                       [ 63,  20, 220])
    lut.add(12105, 'ctx_rh_G_and_S_transv_frontopol',                 [ 13,   0, 250])
    lut.add(12106, 'ctx_rh_G_and_S_cingul-Ant',                       [ 26,  60,   0])
    lut.add(12107, 'ctx_rh_G_and_S_cingul-Mid-Ant',                   [ 26,  60,  75])
    lut.add(12108, 'ctx_rh_G_and_S_cingul-Mid-Post',                  [ 26,  60, 150])
    lut.add(12109, 'ctx_rh_G_cingul-Post-dorsal',                     [ 25,  60, 250])
    lut.add(12110, 'ctx_rh_G_cingul-Post-ventral',                    [ 60,  25,  25])
    lut.add(12111, 'ctx_rh_G_cuneus',                                 [180,  20,  20])
    lut.add(12112, 'ctx_rh_G_front_inf-Opercular',                    [220,  20, 100])
    lut.add(12113, 'ctx_rh_G_front_inf-Orbital',                      [140,  60,  60])
    lut.add(12114, 'ctx_rh_G_front_inf-Triangul',                     [180, 220, 140])
    lut.add(12115, 'ctx_rh_G_front_middle',                           [140, 100, 180])
    lut.add(12116, 'ctx_rh_G_front_sup',                              [180,  20, 140])
    lut.add(12117, 'ctx_rh_G_Ins_lg_and_S_cent_ins',                  [ 23,  10,  10])
    lut.add(12118, 'ctx_rh_G_insular_short',                          [225, 140, 140])
    lut.add(12119, 'ctx_rh_G_occipital_middle',                       [180,  60, 180])
    lut.add(12120, 'ctx_rh_G_occipital_sup',                          [ 20, 220,  60])
    lut.add(12121, 'ctx_rh_G_oc-temp_lat-fusifor',                    [ 60,  20, 140])
    lut.add(12122, 'ctx_rh_G_oc-temp_med-Lingual',                    [220, 180, 140])
    lut.add(12123, 'ctx_rh_G_oc-temp_med-Parahip',                    [ 65, 100,  20])
    lut.add(12124, 'ctx_rh_G_orbital',                                [220,  60,  20])
    lut.add(12125, 'ctx_rh_G_pariet_inf-Angular',                     [ 20,  60, 220])
    lut.add(12126, 'ctx_rh_G_pariet_inf-Supramar',                    [100, 100,  60])
    lut.add(12127, 'ctx_rh_G_parietal_sup',                           [220, 180, 220])
    lut.add(12128, 'ctx_rh_G_postcentral',                            [ 20, 180, 140])
    lut.add(12129, 'ctx_rh_G_precentral',                             [ 60, 140, 180])
    lut.add(12130, 'ctx_rh_G_precuneus',                              [ 25,  20, 140])
    lut.add(12131, 'ctx_rh_G_rectus',                                 [ 20,  60, 100])
    lut.add(12132, 'ctx_rh_G_subcallosal',                            [ 60, 220,  20])
    lut.add(12133, 'ctx_rh_G_temp_sup-G_T_transv',                    [ 60,  60, 220])
    lut.add(12134, 'ctx_rh_G_temp_sup-Lateral',                       [220,  60, 220])
    lut.add(12135, 'ctx_rh_G_temp_sup-Plan_polar',                    [ 65, 220,  60])
    lut.add(12136, 'ctx_rh_G_temp_sup-Plan_tempo',                    [ 25, 140,  20])
    lut.add(12137, 'ctx_rh_G_temporal_inf',                           [220, 220, 100])
    lut.add(12138, 'ctx_rh_G_temporal_middle',                        [180,  60,  60])
    lut.add(12139, 'ctx_rh_Lat_Fis-ant-Horizont',                     [ 61,  20, 220])
    lut.add(12140, 'ctx_rh_Lat_Fis-ant-Vertical',                     [ 61,  20,  60])
    lut.add(12141, 'ctx_rh_Lat_Fis-post',                             [ 61,  60, 100])
    lut.add(12142, 'ctx_rh_Medial_wall',                              [ 25,  25,  25])
    lut.add(12143, 'ctx_rh_Pole_occipital',                           [140,  20,  60])
    lut.add(12144, 'ctx_rh_Pole_temporal',                            [220, 180,  20])
    lut.add(12145, 'ctx_rh_S_calcarine',                              [ 63, 180, 180])
    lut.add(12146, 'ctx_rh_S_central',                                [221,  20,  10])
    lut.add(12147, 'ctx_rh_S_cingul-Marginalis',                      [221,  20, 100])
    lut.add(12148, 'ctx_rh_S_circular_insula_ant',                    [221,  60, 140])
    lut.add(12149, 'ctx_rh_S_circular_insula_inf',                    [221,  20, 220])
    lut.add(12150, 'ctx_rh_S_circular_insula_sup',                    [ 61, 220, 220])
    lut.add(12151, 'ctx_rh_S_collat_transv_ant',                      [100, 200, 200])
    lut.add(12152, 'ctx_rh_S_collat_transv_post',                     [ 10, 200, 200])
    lut.add(12153, 'ctx_rh_S_front_inf',                              [221, 220,  20])
    lut.add(12154, 'ctx_rh_S_front_middle',                           [141,  20, 100])
    lut.add(12155, 'ctx_rh_S_front_sup',                              [ 61, 220, 100])
    lut.add(12156, 'ctx_rh_S_interm_prim-Jensen',                     [141,  60,  20])
    lut.add(12157, 'ctx_rh_S_intrapariet_and_P_trans',                [143,  20, 220])
    lut.add(12158, 'ctx_rh_S_oc_middle_and_Lunatus',                  [101,  60, 220])
    lut.add(12159, 'ctx_rh_S_oc_sup_and_transversal',                 [ 21,  20, 140])
    lut.add(12160, 'ctx_rh_S_occipital_ant',                          [ 61,  20, 180])
    lut.add(12161, 'ctx_rh_S_oc-temp_lat',                            [221, 140,  20])
    lut.add(12162, 'ctx_rh_S_oc-temp_med_and_Lingual',                [141, 100, 220])
    lut.add(12163, 'ctx_rh_S_orbital_lateral',                        [221, 100,  20])
    lut.add(12164, 'ctx_rh_S_orbital_med-olfact',                     [181, 200,  20])
    lut.add(12165, 'ctx_rh_S_orbital-H_Shaped',                       [101,  20,  20])
    lut.add(12166, 'ctx_rh_S_parieto_occipital',                      [101, 100, 180])
    lut.add(12167, 'ctx_rh_S_pericallosal',                           [181, 220,  20])
    lut.add(12168, 'ctx_rh_S_postcentral',                            [ 21, 140, 200])
    lut.add(12169, 'ctx_rh_S_precentral-inf-part',                    [ 21,  20, 240])
    lut.add(12170, 'ctx_rh_S_precentral-sup-part',                    [ 21,  20, 200])
    lut.add(12171, 'ctx_rh_S_suborbital',                             [ 21,  20,  60])
    lut.add(12172, 'ctx_rh_S_subparietal',                            [101,  60,  60])
    lut.add(12173, 'ctx_rh_S_temporal_inf',                           [ 21, 180, 180])
    lut.add(12174, 'ctx_rh_S_temporal_sup',                           [223, 220,  60])
    lut.add(12175, 'ctx_rh_S_temporal_transverse',                    [221,  60,  60])

    lut.add(13100, 'wm_lh_Unknown',                                   [  0,   0,   0])
    lut.add(13101, 'wm_lh_G_and_S_frontomargin',                      [ 23, 220,  60])
    lut.add(13102, 'wm_lh_G_and_S_occipital_inf',                     [ 23,  60, 180])
    lut.add(13103, 'wm_lh_G_and_S_paracentral',                       [ 63, 100,  60])
    lut.add(13104, 'wm_lh_G_and_S_subcentral',                        [ 63,  20, 220])
    lut.add(13105, 'wm_lh_G_and_S_transv_frontopol',                  [ 13,   0, 250])
    lut.add(13106, 'wm_lh_G_and_S_cingul-Ant',                        [ 26,  60,   0])
    lut.add(13107, 'wm_lh_G_and_S_cingul-Mid-Ant',                    [ 26,  60,  75])
    lut.add(13108, 'wm_lh_G_and_S_cingul-Mid-Post',                   [ 26,  60, 150])
    lut.add(13109, 'wm_lh_G_cingul-Post-dorsal',                      [ 25,  60, 250])
    lut.add(13110, 'wm_lh_G_cingul-Post-ventral',                     [ 60,  25,  25])
    lut.add(13111, 'wm_lh_G_cuneus',                                  [180,  20,  20])
    lut.add(13112, 'wm_lh_G_front_inf-Opercular',                     [220,  20, 100])
    lut.add(13113, 'wm_lh_G_front_inf-Orbital',                       [140,  60,  60])
    lut.add(13114, 'wm_lh_G_front_inf-Triangul',                      [180, 220, 140])
    lut.add(13115, 'wm_lh_G_front_middle',                            [140, 100, 180])
    lut.add(13116, 'wm_lh_G_front_sup',                               [180,  20, 140])
    lut.add(13117, 'wm_lh_G_Ins_lg_and_S_cent_ins',                   [ 23,  10,  10])
    lut.add(13118, 'wm_lh_G_insular_short',                           [225, 140, 140])
    lut.add(13119, 'wm_lh_G_occipital_middle',                        [180,  60, 180])
    lut.add(13120, 'wm_lh_G_occipital_sup',                           [ 20, 220,  60])
    lut.add(13121, 'wm_lh_G_oc-temp_lat-fusifor',                     [ 60,  20, 140])
    lut.add(13122, 'wm_lh_G_oc-temp_med-Lingual',                     [220, 180, 140])
    lut.add(13123, 'wm_lh_G_oc-temp_med-Parahip',                     [ 65, 100,  20])
    lut.add(13124, 'wm_lh_G_orbital',                                 [220,  60,  20])
    lut.add(13125, 'wm_lh_G_pariet_inf-Angular',                      [ 20,  60, 220])
    lut.add(13126, 'wm_lh_G_pariet_inf-Supramar',                     [100, 100,  60])
    lut.add(13127, 'wm_lh_G_parietal_sup',                            [220, 180, 220])
    lut.add(13128, 'wm_lh_G_postcentral',                             [ 20, 180, 140])
    lut.add(13129, 'wm_lh_G_precentral',                              [ 60, 140, 180])
    lut.add(13130, 'wm_lh_G_precuneus',                               [ 25,  20, 140])
    lut.add(13131, 'wm_lh_G_rectus',                                  [ 20,  60, 100])
    lut.add(13132, 'wm_lh_G_subcallosal',                             [ 60, 220,  20])
    lut.add(13133, 'wm_lh_G_temp_sup-G_T_transv',                     [ 60,  60, 220])
    lut.add(13134, 'wm_lh_G_temp_sup-Lateral',                        [220,  60, 220])
    lut.add(13135, 'wm_lh_G_temp_sup-Plan_polar',                     [ 65, 220,  60])
    lut.add(13136, 'wm_lh_G_temp_sup-Plan_tempo',                     [ 25, 140,  20])
    lut.add(13137, 'wm_lh_G_temporal_inf',                            [220, 220, 100])
    lut.add(13138, 'wm_lh_G_temporal_middle',                         [180,  60,  60])
    lut.add(13139, 'wm_lh_Lat_Fis-ant-Horizont',                      [ 61,  20, 220])
    lut.add(13140, 'wm_lh_Lat_Fis-ant-Vertical',                      [ 61,  20,  60])
    lut.add(13141, 'wm_lh_Lat_Fis-post',                              [ 61,  60, 100])
    lut.add(13142, 'wm_lh_Medial_wall',                               [ 25,  25,  25])
    lut.add(13143, 'wm_lh_Pole_occipital',                            [140,  20,  60])
    lut.add(13144, 'wm_lh_Pole_temporal',                             [220, 180,  20])
    lut.add(13145, 'wm_lh_S_calcarine',                               [ 63, 180, 180])
    lut.add(13146, 'wm_lh_S_central',                                 [221,  20,  10])
    lut.add(13147, 'wm_lh_S_cingul-Marginalis',                       [221,  20, 100])
    lut.add(13148, 'wm_lh_S_circular_insula_ant',                     [221,  60, 140])
    lut.add(13149, 'wm_lh_S_circular_insula_inf',                     [221,  20, 220])
    lut.add(13150, 'wm_lh_S_circular_insula_sup',                     [ 61, 220, 220])
    lut.add(13151, 'wm_lh_S_collat_transv_ant',                       [100, 200, 200])
    lut.add(13152, 'wm_lh_S_collat_transv_post',                      [ 10, 200, 200])
    lut.add(13153, 'wm_lh_S_front_inf',                               [221, 220,  20])
    lut.add(13154, 'wm_lh_S_front_middle',                            [141,  20, 100])
    lut.add(13155, 'wm_lh_S_front_sup',                               [ 61, 220, 100])
    lut.add(13156, 'wm_lh_S_interm_prim-Jensen',                      [141,  60,  20])
    lut.add(13157, 'wm_lh_S_intrapariet_and_P_trans',                 [143,  20, 220])
    lut.add(13158, 'wm_lh_S_oc_middle_and_Lunatus',                   [101,  60, 220])
    lut.add(13159, 'wm_lh_S_oc_sup_and_transversal',                  [ 21,  20, 140])
    lut.add(13160, 'wm_lh_S_occipital_ant',                           [ 61,  20, 180])
    lut.add(13161, 'wm_lh_S_oc-temp_lat',                             [221, 140,  20])
    lut.add(13162, 'wm_lh_S_oc-temp_med_and_Lingual',                 [141, 100, 220])
    lut.add(13163, 'wm_lh_S_orbital_lateral',                         [221, 100,  20])
    lut.add(13164, 'wm_lh_S_orbital_med-olfact',                      [181, 200,  20])
    lut.add(13165, 'wm_lh_S_orbital-H_Shaped',                        [101,  20,  20])
    lut.add(13166, 'wm_lh_S_parieto_occipital',                       [101, 100, 180])
    lut.add(13167, 'wm_lh_S_pericallosal',                            [181, 220,  20])
    lut.add(13168, 'wm_lh_S_postcentral',                             [ 21, 140, 200])
    lut.add(13169, 'wm_lh_S_precentral-inf-part',                     [ 21,  20, 240])
    lut.add(13170, 'wm_lh_S_precentral-sup-part',                     [ 21,  20, 200])
    lut.add(13171, 'wm_lh_S_suborbital',                              [ 21,  20,  60])
    lut.add(13172, 'wm_lh_S_subparietal',                             [101,  60,  60])
    lut.add(13173, 'wm_lh_S_temporal_inf',                            [ 21, 180, 180])
    lut.add(13174, 'wm_lh_S_temporal_sup',                            [223, 220,  60])
    lut.add(13175, 'wm_lh_S_temporal_transverse',                     [221,  60,  60])
    lut.add(14100, 'wm_rh_Unknown',                                   [  0,   0,   0])
    lut.add(14101, 'wm_rh_G_and_S_frontomargin',                      [ 23, 220,  60])
    lut.add(14102, 'wm_rh_G_and_S_occipital_inf',                     [ 23,  60, 180])
    lut.add(14103, 'wm_rh_G_and_S_paracentral',                       [ 63, 100,  60])
    lut.add(14104, 'wm_rh_G_and_S_subcentral',                        [ 63,  20, 220])
    lut.add(14105, 'wm_rh_G_and_S_transv_frontopol',                  [ 13,   0, 250])
    lut.add(14106, 'wm_rh_G_and_S_cingul-Ant',                        [ 26,  60,   0])
    lut.add(14107, 'wm_rh_G_and_S_cingul-Mid-Ant',                    [ 26,  60,  75])
    lut.add(14108, 'wm_rh_G_and_S_cingul-Mid-Post',                   [ 26,  60, 150])
    lut.add(14109, 'wm_rh_G_cingul-Post-dorsal',                      [ 25,  60, 250])
    lut.add(14110, 'wm_rh_G_cingul-Post-ventral',                     [ 60,  25,  25])
    lut.add(14111, 'wm_rh_G_cuneus',                                  [180,  20,  20])
    lut.add(14112, 'wm_rh_G_front_inf-Opercular',                     [220,  20, 100])
    lut.add(14113, 'wm_rh_G_front_inf-Orbital',                       [140,  60,  60])
    lut.add(14114, 'wm_rh_G_front_inf-Triangul',                      [180, 220, 140])
    lut.add(14115, 'wm_rh_G_front_middle',                            [140, 100, 180])
    lut.add(14116, 'wm_rh_G_front_sup',                               [180,  20, 140])
    lut.add(14117, 'wm_rh_G_Ins_lg_and_S_cent_ins',                   [ 23,  10,  10])
    lut.add(14118, 'wm_rh_G_insular_short',                           [225, 140, 140])
    lut.add(14119, 'wm_rh_G_occipital_middle',                        [180,  60, 180])
    lut.add(14120, 'wm_rh_G_occipital_sup',                           [ 20, 220,  60])
    lut.add(14121, 'wm_rh_G_oc-temp_lat-fusifor',                     [ 60,  20, 140])
    lut.add(14122, 'wm_rh_G_oc-temp_med-Lingual',                     [220, 180, 140])
    lut.add(14123, 'wm_rh_G_oc-temp_med-Parahip',                     [ 65, 100,  20])
    lut.add(14124, 'wm_rh_G_orbital',                                 [220,  60,  20])
    lut.add(14125, 'wm_rh_G_pariet_inf-Angular',                      [ 20,  60, 220])
    lut.add(14126, 'wm_rh_G_pariet_inf-Supramar',                     [100, 100,  60])
    lut.add(14127, 'wm_rh_G_parietal_sup',                            [220, 180, 220])
    lut.add(14128, 'wm_rh_G_postcentral',                             [ 20, 180, 140])
    lut.add(14129, 'wm_rh_G_precentral',                              [ 60, 140, 180])
    lut.add(14130, 'wm_rh_G_precuneus',                               [ 25,  20, 140])
    lut.add(14131, 'wm_rh_G_rectus',                                  [ 20,  60, 100])
    lut.add(14132, 'wm_rh_G_subcallosal',                             [ 60, 220,  20])
    lut.add(14133, 'wm_rh_G_temp_sup-G_T_transv',                     [ 60,  60, 220])
    lut.add(14134, 'wm_rh_G_temp_sup-Lateral',                        [220,  60, 220])
    lut.add(14135, 'wm_rh_G_temp_sup-Plan_polar',                     [ 65, 220,  60])
    lut.add(14136, 'wm_rh_G_temp_sup-Plan_tempo',                     [ 25, 140,  20])
    lut.add(14137, 'wm_rh_G_temporal_inf',                            [220, 220, 100])
    lut.add(14138, 'wm_rh_G_temporal_middle',                         [180,  60,  60])
    lut.add(14139, 'wm_rh_Lat_Fis-ant-Horizont',                      [ 61,  20, 220])
    lut.add(14140, 'wm_rh_Lat_Fis-ant-Vertical',                      [ 61,  20,  60])
    lut.add(14141, 'wm_rh_Lat_Fis-post',                              [ 61,  60, 100])
    lut.add(14142, 'wm_rh_Medial_wall',                               [ 25,  25,  25])
    lut.add(14143, 'wm_rh_Pole_occipital',                            [140,  20,  60])
    lut.add(14144, 'wm_rh_Pole_temporal',                             [220, 180,  20])
    lut.add(14145, 'wm_rh_S_calcarine',                               [ 63, 180, 180])
    lut.add(14146, 'wm_rh_S_central',                                 [221,  20,  10])
    lut.add(14147, 'wm_rh_S_cingul-Marginalis',                       [221,  20, 100])
    lut.add(14148, 'wm_rh_S_circular_insula_ant',                     [221,  60, 140])
    lut.add(14149, 'wm_rh_S_circular_insula_inf',                     [221,  20, 220])
    lut.add(14150, 'wm_rh_S_circular_insula_sup',                     [ 61, 220, 220])
    lut.add(14151, 'wm_rh_S_collat_transv_ant',                       [100, 200, 200])
    lut.add(14152, 'wm_rh_S_collat_transv_post',                      [ 10, 200, 200])
    lut.add(14153, 'wm_rh_S_front_inf',                               [221, 220,  20])
    lut.add(14154, 'wm_rh_S_front_middle',                            [141,  20, 100])
    lut.add(14155, 'wm_rh_S_front_sup',                               [ 61, 220, 100])
    lut.add(14156, 'wm_rh_S_interm_prim-Jensen',                      [141,  60,  20])
    lut.add(14157, 'wm_rh_S_intrapariet_and_P_trans',                 [143,  20, 220])
    lut.add(14158, 'wm_rh_S_oc_middle_and_Lunatus',                   [101,  60, 220])
    lut.add(14159, 'wm_rh_S_oc_sup_and_transversal',                  [ 21,  20, 140])
    lut.add(14160, 'wm_rh_S_occipital_ant',                           [ 61,  20, 180])
    lut.add(14161, 'wm_rh_S_oc-temp_lat',                             [221, 140,  20])
    lut.add(14162, 'wm_rh_S_oc-temp_med_and_Lingual',                 [141, 100, 220])
    lut.add(14163, 'wm_rh_S_orbital_lateral',                         [221, 100,  20])
    lut.add(14164, 'wm_rh_S_orbital_med-olfact',                      [181, 200,  20])
    lut.add(14165, 'wm_rh_S_orbital-H_Shaped',                        [101,  20,  20])
    lut.add(14166, 'wm_rh_S_parieto_occipital',                       [101, 100, 180])
    lut.add(14167, 'wm_rh_S_pericallosal',                            [181, 220,  20])
    lut.add(14168, 'wm_rh_S_postcentral',                             [ 21, 140, 200])
    lut.add(14169, 'wm_rh_S_precentral-inf-part',                     [ 21,  20, 240])
    lut.add(14170, 'wm_rh_S_precentral-sup-part',                     [ 21,  20, 200])
    lut.add(14171, 'wm_rh_S_suborbital',                              [ 21,  20,  60])
    lut.add(14172, 'wm_rh_S_subparietal',                             [101,  60,  60])
    lut.add(14173, 'wm_rh_S_temporal_inf',                            [ 21, 180, 180])
    lut.add(14174, 'wm_rh_S_temporal_sup',                            [223, 220,  60])
    lut.add(14175, 'wm_rh_S_temporal_transverse',                     [221,  60,  60])
    
    return lut


def destrieux():
    """
    """
    lut = LookupTable()
    lut.add(0,  'Unknown',                   [  0,   0,   0])
    lut.add(1,  'G_and_S_frontomargin',      [ 23, 220,  60])
    lut.add(2,  'G_and_S_occipital_inf',     [ 23,  60, 180])
    lut.add(3,  'G_and_S_paracentral',       [ 63, 100,  60])
    lut.add(4,  'G_and_S_subcentral',        [ 63,  20, 220])
    lut.add(5,  'G_and_S_transv_frontopol',  [ 13,   0, 250])
    lut.add(6,  'G_and_S_cingul-Ant',        [ 26,  60,   0])
    lut.add(7,  'G_and_S_cingul-Mid-Ant',    [ 26,  60,  75])
    lut.add(8,  'G_and_S_cingul-Mid-Post',   [ 26,  60, 150])
    lut.add(9,  'G_cingul-Post-dorsal',      [ 25,  60, 250])
    lut.add(10, 'G_cingul-Post-ventral',     [ 60,  25,  25])
    lut.add(11, 'G_cuneus',                  [180,  20,  20])
    lut.add(12, 'G_front_inf-Opercular',     [220,  20, 100])
    lut.add(13, 'G_front_inf-Orbital',       [140,  60,  60])
    lut.add(14, 'G_front_inf-Triangul',      [180, 220, 140])
    lut.add(15, 'G_front_middle',            [140, 100, 180])
    lut.add(16, 'G_front_sup',               [180,  20, 140])
    lut.add(17, 'G_Ins_lg_and_S_cent_ins',   [ 23,  10,  10])
    lut.add(18, 'G_insular_short',           [225, 140, 140])
    lut.add(19, 'G_occipital_middle',        [180,  60, 180])
    lut.add(20, 'G_occipital_sup',           [ 20, 220,  60])
    lut.add(21, 'G_oc-temp_lat-fusifor',     [ 60,  20, 140])
    lut.add(22, 'G_oc-temp_med-Lingual',     [220, 180, 140])
    lut.add(23, 'G_oc-temp_med-Parahip',     [ 65, 100,  20])
    lut.add(24, 'G_orbital',                 [220,  60,  20])
    lut.add(25, 'G_pariet_inf-Angular',      [ 20,  60, 220])
    lut.add(26, 'G_pariet_inf-Supramar',     [100, 100,  60])
    lut.add(27, 'G_parietal_sup',            [220, 180, 220])
    lut.add(28, 'G_postcentral',             [ 20, 180, 140])
    lut.add(29, 'G_precentral',              [ 60, 140, 180])
    lut.add(30, 'G_precuneus',               [ 25,  20, 140])
    lut.add(31, 'G_rectus',                  [ 20,  60, 100])
    lut.add(32, 'G_subcallosal',             [ 60, 220,  20])
    lut.add(33, 'G_temp_sup-G_T_transv',     [ 60,  60, 220])
    lut.add(34, 'G_temp_sup-Lateral',        [220,  60, 220])
    lut.add(35, 'G_temp_sup-Plan_polar',     [ 65, 220,  60])
    lut.add(36, 'G_temp_sup-Plan_tempo',     [ 25, 140,  20])
    lut.add(37, 'G_temporal_inf',            [220, 220, 100])
    lut.add(38, 'G_temporal_middle',         [180,  60,  60])
    lut.add(39, 'Lat_Fis-ant-Horizont',      [ 61,  20, 220])
    lut.add(40, 'Lat_Fis-ant-Vertical',      [ 61,  20,  60])
    lut.add(41, 'Lat_Fis-post',              [ 61,  60, 100])
    lut.add(42, 'Medial_wall',               [ 25,  25,  25])
    lut.add(43, 'Pole_occipital',            [140,  20,  60])
    lut.add(44, 'Pole_temporal',             [220, 180,  20])
    lut.add(45, 'S_calcarine',               [ 63, 180, 180])
    lut.add(46, 'S_central',                 [221,  20,  10])
    lut.add(47, 'S_cingul-Marginalis',       [221,  20, 100])
    lut.add(48, 'S_circular_insula_ant',     [221,  60, 140])
    lut.add(49, 'S_circular_insula_inf',     [221,  20, 220])
    lut.add(50, 'S_circular_insula_sup',     [ 61, 220, 220])
    lut.add(51, 'S_collat_transv_ant',       [100, 200, 200])
    lut.add(52, 'S_collat_transv_post',      [ 10, 200, 200])
    lut.add(53, 'S_front_inf',               [221, 220,  20])
    lut.add(54, 'S_front_middle',            [141,  20, 100])
    lut.add(55, 'S_front_sup',               [ 61, 220, 100])
    lut.add(56, 'S_interm_prim-Jensen',      [141,  60,  20])
    lut.add(57, 'S_intrapariet_and_P_trans', [143,  20, 220])
    lut.add(58, 'S_oc_middle_and_Lunatus',   [101,  60, 220])
    lut.add(59, 'S_oc_sup_and_transversal',  [ 21,  20, 140])
    lut.add(60, 'S_occipital_ant',           [ 61,  20, 180])
    lut.add(61, 'S_oc-temp_lat',             [221, 140,  20])
    lut.add(62, 'S_oc-temp_med_and_Lingual', [141, 100, 220])
    lut.add(63, 'S_orbital_lateral',         [221, 100,  20])
    lut.add(64, 'S_orbital_med-olfact',      [181, 200,  20])
    lut.add(65, 'S_orbital-H_Shaped',        [101,  20,  20])
    lut.add(66, 'S_parieto_occipital',       [101, 100, 180])
    lut.add(67, 'S_pericallosal',            [181, 220,  20])
    lut.add(68, 'S_postcentral',             [ 21, 140, 200])
    lut.add(69, 'S_precentral-inf-part',     [ 21,  20, 240])
    lut.add(70, 'S_precentral-sup-part',     [ 21,  20, 200])
    lut.add(71, 'S_suborbital',              [ 21,  20,  60])
    lut.add(72, 'S_subparietal',             [101,  60,  60])
    lut.add(73, 'S_temporal_inf',            [ 21, 180, 180])
    lut.add(74, 'S_temporal_sup',            [223, 220,  60])
    lut.add(75, 'S_temporal_transverse',     [221,  60,  60])
    return lut


def dkt():
    """
    """
    lut = LookupTable()
    lut.add(0,  'unknown',                  [ 25,   5,  25])
    lut.add(1,  'bankssts',                 [ 25, 100,  40])
    lut.add(2,  'caudalanteriorcingulate',  [125, 100, 160])
    lut.add(3,  'caudalmiddlefrontal',      [100,  25,   0])
    lut.add(4,  'corpuscallosum',           [120,  70,  50])
    lut.add(5,  'cuneus',                   [220,  20, 100])
    lut.add(6,  'entorhinal',               [220,  20,  10])
    lut.add(7,  'fusiform',                 [180, 220, 140])
    lut.add(8,  'inferiorparietal',         [220,  60, 220])
    lut.add(9,  'inferiortemporal',         [180,  40, 120])
    lut.add(10, 'isthmuscingulate',         [140,  20, 140])
    lut.add(11, 'lateraloccipital',         [ 20,  30, 140])
    lut.add(12, 'lateralorbitofrontal',     [ 35,  75,  50])
    lut.add(13, 'lingual',                  [225, 140, 140])
    lut.add(14, 'medialorbitofrontal',      [200,  35,  75])
    lut.add(15, 'middletemporal',           [160, 100,  50])
    lut.add(16, 'parahippocampal',          [ 20, 220,  60])
    lut.add(17, 'paracentral',              [ 60, 220,  60])
    lut.add(18, 'parsopercularis',          [220, 180, 140])
    lut.add(19, 'parsorbitalis',            [ 20, 100,  50])
    lut.add(20, 'parstriangularis',         [220,  60,  20])
    lut.add(21, 'pericalcarine',            [120, 100,  60])
    lut.add(22, 'postcentral',              [220,  20,  20])
    lut.add(23, 'posteriorcingulate',       [220, 180, 220])
    lut.add(24, 'precentral',               [ 60,  20, 220])
    lut.add(25, 'precuneus',                [160, 140, 180])
    lut.add(26, 'rostralanteriorcingulate', [ 80,  20, 140])
    lut.add(27, 'rostralmiddlefrontal',     [ 75,  50, 125])
    lut.add(28, 'superiorfrontal',          [ 20, 220, 160])
    lut.add(29, 'superiorparietal',         [ 20, 180, 140])
    lut.add(30, 'superiortemporal',         [140, 220, 220])
    lut.add(31, 'supramarginal',            [ 80, 160,  20])
    lut.add(32, 'frontalpole',              [100,   0, 100])
    lut.add(33, 'temporalpole',             [ 70,  20, 170])
    lut.add(34, 'transversetemporal',       [150, 150, 200])
    lut.add(35, 'insula',                   [255, 192,  32])
    return lut


def tissue_type():
    """
    """
    lut = LookupTable()
    lut.add(0, 'Unknown',                  [0,   0,   0])
    lut.add(1, 'Cortex',                   [205, 62,  78])
    lut.add(2, 'Subcortical-Gray-Matter',  [230, 148, 34])
    lut.add(3, 'White-Matter',             [245, 245, 245])
    lut.add(4, 'CSF',                      [120, 18,  134])
    lut.add(5, 'Head',                     [150, 150, 200])
    return lut


def tissue_type_no_skull():
    """
    """
    lut = LookupTable()
    lut.add(0, 'Unknown',                  [0,   0,   0])
    lut.add(1, 'Cortex',                   [205, 62,  78])
    lut.add(2, 'Subcortical-Gray-Matter',  [230, 148, 34])
    lut.add(3, 'White-Matter',             [245, 245, 245])
    lut.add(4, 'CSF',                      [120, 18,  134])
    return lut


def tissue_type_recoder():
    """
    Returns a recoding table that converts default brain labels to the
    corresponding tissue-type.

    Returns:
        RecodingLookupTable: .
    """
    rlut = RecodingLookupTable()
    rlut.target_lut = tissue_type()
    rlut.mapping = {
        0:    0,  # Unknown
        2:    3,  # Left-Cerebral-White-Matter
        3:    1,  # Left-Cerebral-Cortex
        4:    4,  # Left-Lateral-Ventricle
        5:    4,  # Left-Inf-Lat-Vent
        7:    3,  # Left-Cerebellum-White-Matter
        8:    2,  # Left-Cerebellum-Cortex
        10:   2,  # Left-Thalamus
        11:   2,  # Left-Caudate
        12:   2,  # Left-Putamen
        13:   2,  # Left-Pallidum
        14:   4,  # 3rd-Ventricle
        15:   4,  # 4th-Ventricle
        16:   3,  # Brain-Stem
        17:   2,  # Left-Hippocampus
        18:   2,  # Left-Amygdala
        24:   4,  # CSF
        26:   2,  # Left-Accumbens-Area
        28:   3,  # Left-VentralDC
        30:   4,  # Left-Vessel
        31:   4,  # Left-Choroid-Plexus
        41:   3,  # Right-Cerebral-White-Matter
        42:   1,  # Right-Cerebral-Cortex
        43:   4,  # Right-Lateral-Ventricle
        44:   4,  # Right-Inf-Lat-Vent
        46:   3,  # Right-Cerebellum-White-Matter
        47:   2,  # Right-Cerebellum-Cortex
        49:   2,  # Right-Thalamus
        50:   2,  # Right-Caudate
        51:   2,  # Right-Putamen
        52:   2,  # Right-Pallidum
        53:   2,  # Right-Hippocampus
        54:   2,  # Right-Amygdala
        58:   2,  # Right-Accumbens-Area
        60:   3,  # Right-VentralDC
        62:   4,  # Right-Vessel
        63:   4,  # Right-Choroid-Plexus
        77:   3,  # WM-Hypointensities
        78:   3,  # Left-WM-Hypointensities
        79:   3,  # Right-WM-Hypointensities
        80:   2,  # Non-WM-Hypointensities
        81:   2,  # Left-Non-WM-Hypointensities
        82:   2,  # Right-Non-WM-Hypointensities
        85:   3,  # Optic-Chiasm
        130:  5,  # Air
        165:  5,  # Skull
        172:  2,  # Vermis
        174:  3,  # Pons
        251:  3,  # CC_Posterior
        252:  3,  # CC_Mid_Posterior
        253:  3,  # CC_Central
        254:  3,  # CC_Mid_Anterior
        255:  3,  # CC_Anterior
        257:  4,  # CSF-ExtraCerebral
        258:  5,  # Head-ExtraCerebral
        1001: 1,  # ctx-lh-bankssts
        1002: 1,  # ctx-lh-caudalanteriorcingulate
        1003: 1,  # ctx-lh-caudalmiddlefrontal
        1005: 1,  # ctx-lh-cuneus
        1006: 1,  # ctx-lh-entorhinal
        1007: 1,  # ctx-lh-fusiform
        1008: 1,  # ctx-lh-inferiorparietal
        1009: 1,  # ctx-lh-inferiortemporal
        1010: 1,  # ctx-lh-isthmuscingulate
        1011: 1,  # ctx-lh-lateraloccipital
        1012: 1,  # ctx-lh-lateralorbitofrontal
        1013: 1,  # ctx-lh-lingual
        1014: 1,  # ctx-lh-medialorbitofrontal
        1015: 1,  # ctx-lh-middletemporal
        1016: 1,  # ctx-lh-parahippocampal
        1017: 1,  # ctx-lh-paracentral
        1018: 1,  # ctx-lh-parsopercularis
        1019: 1,  # ctx-lh-parsorbitalis
        1020: 1,  # ctx-lh-parstriangularis
        1021: 1,  # ctx-lh-pericalcarine
        1022: 1,  # ctx-lh-postcentral
        1023: 1,  # ctx-lh-posteriorcingulate
        1024: 1,  # ctx-lh-precentral
        1025: 1,  # ctx-lh-precuneus
        1026: 1,  # ctx-lh-rostralanteriorcingulate
        1027: 1,  # ctx-lh-rostralmiddlefrontal
        1028: 1,  # ctx-lh-superiorfrontal
        1029: 1,  # ctx-lh-superiorparietal
        1030: 1,  # ctx-lh-superiortemporal
        1031: 1,  # ctx-lh-supramarginal
        1032: 1,  # ctx-lh-frontalpole
        1033: 1,  # ctx-lh-temporalpole
        1034: 1,  # ctx-lh-transversetemporal
        1035: 1,  # ctx-lh-insula
        2001: 1,  # ctx-rh-bankssts
        2002: 1,  # ctx-rh-caudalanteriorcingulate
        2003: 1,  # ctx-rh-caudalmiddlefrontal
        2005: 1,  # ctx-rh-cuneus
        2006: 1,  # ctx-rh-entorhinal
        2007: 1,  # ctx-rh-fusiform
        2008: 1,  # ctx-rh-inferiorparietal
        2009: 1,  # ctx-rh-inferiortemporal
        2010: 1,  # ctx-rh-isthmuscingulate
        2011: 1,  # ctx-rh-lateraloccipital
        2012: 1,  # ctx-rh-lateralorbitofrontal
        2013: 1,  # ctx-rh-lingual
        2014: 1,  # ctx-rh-medialorbitofrontal
        2015: 1,  # ctx-rh-middletemporal
        2016: 1,  # ctx-rh-parahippocampal
        2017: 1,  # ctx-rh-paracentral
        2018: 1,  # ctx-rh-parsopercularis
        2019: 1,  # ctx-rh-parsorbitalis
        2020: 1,  # ctx-rh-parstriangularis
        2021: 1,  # ctx-rh-pericalcarine
        2022: 1,  # ctx-rh-postcentral
        2023: 1,  # ctx-rh-posteriorcingulate
        2024: 1,  # ctx-rh-precentral
        2025: 1,  # ctx-rh-precuneus
        2026: 1,  # ctx-rh-rostralanteriorcingulate
        2027: 1,  # ctx-rh-rostralmiddlefrontal
        2028: 1,  # ctx-rh-superiorfrontal
        2029: 1,  # ctx-rh-superiorparietal
        2030: 1,  # ctx-rh-superiortemporal
        2031: 1,  # ctx-rh-supramarginal
        2032: 1,  # ctx-rh-frontalpole
        2033: 1,  # ctx-rh-temporalpole
        2034: 1,  # ctx-rh-transversetemporal
        2035: 1,  # ctx-rh-insula
    }
    return rlut

def tissue_type_recoder_no_skull():
    """
    Returns a recoding table that converts default brain labels to the
    corresponding tissue-type.

    Returns:
        RecodingLookupTable: .
    """
    rlut = RecodingLookupTable()
    rlut.target_lut = tissue_type_no_skull()
    rlut.mapping = {
        0:    0,  # Unknown
        2:    3,  # Left-Cerebral-White-Matter
        3:    1,  # Left-Cerebral-Cortex
        4:    4,  # Left-Lateral-Ventricle
        5:    4,  # Left-Inf-Lat-Vent
        7:    3,  # Left-Cerebellum-White-Matter
        8:    2,  # Left-Cerebellum-Cortex
        10:   2,  # Left-Thalamus
        11:   2,  # Left-Caudate
        12:   2,  # Left-Putamen
        13:   2,  # Left-Pallidum
        14:   4,  # 3rd-Ventricle
        15:   4,  # 4th-Ventricle
        16:   3,  # Brain-Stem
        17:   2,  # Left-Hippocampus
        18:   2,  # Left-Amygdala
        24:   4,  # CSF
        26:   2,  # Left-Accumbens-Area
        28:   3,  # Left-VentralDC
        30:   4,  # Left-Vessel
        31:   4,  # Left-Choroid-Plexus
        41:   3,  # Right-Cerebral-White-Matter
        42:   1,  # Right-Cerebral-Cortex
        43:   4,  # Right-Lateral-Ventricle
        44:   4,  # Right-Inf-Lat-Vent
        46:   3,  # Right-Cerebellum-White-Matter
        47:   2,  # Right-Cerebellum-Cortex
        49:   2,  # Right-Thalamus
        50:   2,  # Right-Caudate
        51:   2,  # Right-Putamen
        52:   2,  # Right-Pallidum
        53:   2,  # Right-Hippocampus
        54:   2,  # Right-Amygdala
        58:   2,  # Right-Accumbens-Area
        60:   3,  # Right-VentralDC
        62:   4,  # Right-Vessel
        63:   4,  # Right-Choroid-Plexus
        77:   3,  # WM-Hypointensities
        78:   3,  # Left-WM-Hypointensities
        79:   3,  # Right-WM-Hypointensities
        80:   2,  # Non-WM-Hypointensities
        81:   2,  # Left-Non-WM-Hypointensities
        82:   2,  # Right-Non-WM-Hypointensities
        85:   3,  # Optic-Chiasm
        130:  0,  # Air
        165:  0,  # Skull
        172:  2,  # Vermis
        174:  3,  # Pons
        251:  3,  # CC_Posterior
        252:  3,  # CC_Mid_Posterior
        253:  3,  # CC_Central
        254:  3,  # CC_Mid_Anterior
        255:  3,  # CC_Anterior
        257:  4,  # CSF-ExtraCerebral
        258:  0,  # Head-ExtraCerebral
        1001: 1,  # ctx-lh-bankssts
        1002: 1,  # ctx-lh-caudalanteriorcingulate
        1003: 1,  # ctx-lh-caudalmiddlefrontal
        1005: 1,  # ctx-lh-cuneus
        1006: 1,  # ctx-lh-entorhinal
        1007: 1,  # ctx-lh-fusiform
        1008: 1,  # ctx-lh-inferiorparietal
        1009: 1,  # ctx-lh-inferiortemporal
        1010: 1,  # ctx-lh-isthmuscingulate
        1011: 1,  # ctx-lh-lateraloccipital
        1012: 1,  # ctx-lh-lateralorbitofrontal
        1013: 1,  # ctx-lh-lingual
        1014: 1,  # ctx-lh-medialorbitofrontal
        1015: 1,  # ctx-lh-middletemporal
        1016: 1,  # ctx-lh-parahippocampal
        1017: 1,  # ctx-lh-paracentral
        1018: 1,  # ctx-lh-parsopercularis
        1019: 1,  # ctx-lh-parsorbitalis
        1020: 1,  # ctx-lh-parstriangularis
        1021: 1,  # ctx-lh-pericalcarine
        1022: 1,  # ctx-lh-postcentral
        1023: 1,  # ctx-lh-posteriorcingulate
        1024: 1,  # ctx-lh-precentral
        1025: 1,  # ctx-lh-precuneus
        1026: 1,  # ctx-lh-rostralanteriorcingulate
        1027: 1,  # ctx-lh-rostralmiddlefrontal
        1028: 1,  # ctx-lh-superiorfrontal
        1029: 1,  # ctx-lh-superiorparietal
        1030: 1,  # ctx-lh-superiortemporal
        1031: 1,  # ctx-lh-supramarginal
        1032: 1,  # ctx-lh-frontalpole
        1033: 1,  # ctx-lh-temporalpole
        1034: 1,  # ctx-lh-transversetemporal
        1035: 1,  # ctx-lh-insula
        2001: 1,  # ctx-rh-bankssts
        2002: 1,  # ctx-rh-caudalanteriorcingulate
        2003: 1,  # ctx-rh-caudalmiddlefrontal
        2005: 1,  # ctx-rh-cuneus
        2006: 1,  # ctx-rh-entorhinal
        2007: 1,  # ctx-rh-fusiform
        2008: 1,  # ctx-rh-inferiorparietal
        2009: 1,  # ctx-rh-inferiortemporal
        2010: 1,  # ctx-rh-isthmuscingulate
        2011: 1,  # ctx-rh-lateraloccipital
        2012: 1,  # ctx-rh-lateralorbitofrontal
        2013: 1,  # ctx-rh-lingual
        2014: 1,  # ctx-rh-medialorbitofrontal
        2015: 1,  # ctx-rh-middletemporal
        2016: 1,  # ctx-rh-parahippocampal
        2017: 1,  # ctx-rh-paracentral
        2018: 1,  # ctx-rh-parsopercularis
        2019: 1,  # ctx-rh-parsorbitalis
        2020: 1,  # ctx-rh-parstriangularis
        2021: 1,  # ctx-rh-pericalcarine
        2022: 1,  # ctx-rh-postcentral
        2023: 1,  # ctx-rh-posteriorcingulate
        2024: 1,  # ctx-rh-precentral
        2025: 1,  # ctx-rh-precuneus
        2026: 1,  # ctx-rh-rostralanteriorcingulate
        2027: 1,  # ctx-rh-rostralmiddlefrontal
        2028: 1,  # ctx-rh-superiorfrontal
        2029: 1,  # ctx-rh-superiorparietal
        2030: 1,  # ctx-rh-superiortemporal
        2031: 1,  # ctx-rh-supramarginal
        2032: 1,  # ctx-rh-frontalpole
        2033: 1,  # ctx-rh-temporalpole
        2034: 1,  # ctx-rh-transversetemporal
        2035: 1,  # ctx-rh-insula
    }
    return rlut
