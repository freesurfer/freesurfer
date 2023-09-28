import yaml
import numpy as np

# Background
BG_LABELS = [0]
BG_LABELS_ASEG = [0]

# White matter labels
WM_LABELS = [2, 7, 68, 120, 130, 199, 210, 445, 513, 514, 602, 100, 161, 208, 209, 114]
WM_LABELS_ASEG = [2, 41]

# Gray matter
GM_LABELS = [
    3,
    103,
    106,
    108,
    111,
    117,
    127,
    128,
    129,
    146,
    158,
    99,
    113,
    125,
    155,
    156,
    157,
]
GM_LABELS_ASEG = [3, 42]

# Lateral and inferio lateral ventricle
VENT_LABELS = [50, 383, 620, 688, 328, 291]
VENT_LABELS_ASEG = [4, 5, 43, 44, 24, 14, 15]

# WM cerebellum
WM_CEREBELLUM = [611, 715, 721]
WM_CEREBELLUM_ASEG = [7, 46]

# GM cerebellum
GM_CEREBELLUM = [595, 597]
GM_CEREBELLUM_ASEG = [8, 47]

# Caudate and accumbens - GM nuclei
CAUDATE_AND_ACCUMBENS = [48, 118, 393, 101, 49, 184, 187, 69]
CAUDATE_AND_ACCUMBENS_ASEG = [11, 50, 26, 58]

# Putamen
PUTAMEN = [79, 349]
PUTAMEN_ASEG = [12, 51]

# Claustrum and temporal claustrum (the latter is more WM-ish)
CLAUSTRUM = [102, 174]
CLAUSTRUM_ASEG = []

# Pallidum
PALLIDUM = [119, 206]
PALLIDUM_ASEG = [13, 52]

# Hippocampus
HIPPOCAMPUS = [
    343,
    368,
    372,
    408,
    418,
    562,
    566,
    571,
    339,
    354,
    461,
    344,
    369,
    409,
    567,
    340,
    563,
    419,
    373,
    572,
    345,
    370,
    374,
    410,
    420,
    564,
    573,
    341,
    568,
    371,
    375,
    411,
    421,
    565,
    569,
    574,
    342,
    346,
    326,
    347,
    576,
    348,
    327,
    358,
    558,
    404,
    364,
    405,
    559,
    365,
    561,
    407,
    367,
    500,
    575,
    422,
    432,
    320,
]
HIPPOCAMPUS_ASEG = [17, 53]

# Thalamus
THALAMUS = [
    276,
    444,
    430,
    809,
    815,
    201,
    314,
    381,
    382,
    394,
    458,
    283,
    493,
    519,
    812,
    477,
    254,
    191,
    312,
    285,
    425,
    284,
    396,
    479,
    350,
    222,
    397,
    459,
    484,
    218,
    492,
    322,
    442,
    190,
    818,
    246,
    200,
    517,
    478,
    505,
    226,
    303,
    398,
    426,
    424,
    317,
    400,
    274,
    395,
    423,
    811,
    441,
    823,
    504,
    511,
    512,
    510,
    378,
    221,
    379,
    313,
    282,
    225,
    808,
    224,
    399,
    286,
    223,
    454,
    380,
    296,
    814,
    816,
    813,
    219,
    220,
    252,
    253,
    443,
    508,
    578,
]
THALAMUS_ASEG = [10, 49]

# Brain stem
BRAINSTEM = [
    1,
    385,
    472,
    525,
    451,
    526,
    384,
    401,
    431,
    467,
    468,
    469,
    470,
    496,
    498,
    403,
    352,
    449,
    507,
    386,
    310,
    465,
    412,
    414,
    415,
    435,
    436,
    447,
    464,
    435,
    471,
    499,
    506,
    515,
    520,
    521,
    527,
    528,
    531,
    535,
    541,
    542,
    543,
    546,
    580,
    581,
    586,
    589,
    591,
    592,
    598,
    599,
    600,
    610,
    617,
    624,
    630,
    645,
    647,
    649,
    654,
    655,
    659,
    660,
    662,
    666,
    669,
    671,
    677,
    682,
    683,
    685,
    687,
    689,
    690,
    693,
    694,
    695,
    697,
    703,
    717,
    724,
    725,
    726,
    731,
    751,
    752,
    755,
    756,
    765,
    771,
    779,
    827,
    828,
    829,
    842,
    843,
    844,
]
BRAINSTEM_ASEG = [16, 161, 162]

# Hypothalamus
HYPOTHALAMUS = [
    256,
    275,
    287,
    304,
    147,
    150,
    164,
    179,
    192,
    193,
    194,
    195,
    197,
    207,
    227,
    228,
    229,
    243,
    244,
    245,
    255,
    268,
    297,
    308,
    149,
    160,
    181,
    182,
    196,
    198,
    230,
    231,
    232,
    269,
    305,
    306,
    307,
    333,
    351,
    316,
    321,
    234,
    289,
    298,
    309,
    315,
    325,
    433,
]
HYPOTHALAMUS_ASEG = []

# Amygdala
AMYGDALA = [
    238,
    249,
    277,
    278,
    279,
    292,
    300,
    301,
    319,
    376,
    377,
    217,
    240,
    242,
    214,
    215,
    216,
    295,
    249,
    292,
    300,
    319,
    278,
    376,
]
AMYGDALA_ASEG = [18, 54]

# Cortex
CORTEX = list(range(1000, 3001))
CORTEX_ASEG = list(range(1000, 3001))

# Let's start with a handful of tissues, say WM, GM, CSF and background
if False:  # merge background and CSF (3D version)
    WM_TISSUES = WM_LABELS + WM_CEREBELLUM + BRAINSTEM + PALLIDUM
    GM_TISSUES = (
        GM_LABELS
        + GM_CEREBELLUM
        + CAUDATE_AND_ACCUMBENS
        + PUTAMEN
        + HIPPOCAMPUS
        + THALAMUS
        + HYPOTHALAMUS
        + AMYGDALA
        + CORTEX
        + CLAUSTRUM
    )
    BG_TISSUES = BG_LABELS + VENT_LABELS
    TISSUES = [BG_TISSUES, WM_TISSUES, GM_TISSUES]

    WM_TISSUES_ASEG = (
        WM_LABELS_ASEG + WM_CEREBELLUM_ASEG + BRAINSTEM_ASEG + PALLIDUM_ASEG
    )
    GM_TISSUES_ASEG = (
        GM_LABELS_ASEG
        + GM_CEREBELLUM_ASEG
        + CAUDATE_AND_ACCUMBENS_ASEG
        + PUTAMEN_ASEG
        + HIPPOCAMPUS_ASEG
        + THALAMUS_ASEG
        + HYPOTHALAMUS_ASEG
        + AMYGDALA_ASEG
        + CORTEX_ASEG
        + CLAUSTRUM_ASEG
    )
    BG_TISSUES_ASEG = BG_LABELS_ASEG + VENT_LABELS_ASEG
    TISSUES_ASEG = [BG_TISSUES_ASEG, WM_TISSUES_ASEG, GM_TISSUES_ASEG]

elif False:  # do not merge background and CSF (4d version)
    WM_TISSUES = WM_LABELS + WM_CEREBELLUM + BRAINSTEM + PALLIDUM
    GM_TISSUES = (
        GM_LABELS
        + GM_CEREBELLUM
        + CAUDATE_AND_ACCUMBENS
        + PUTAMEN
        + HIPPOCAMPUS
        + THALAMUS
        + HYPOTHALAMUS
        + AMYGDALA
        + CORTEX
        + CLAUSTRUM
    )
    CSF_TISSUES = VENT_LABELS
    BG_TISSUES = BG_LABELS
    TISSUES = [BG_TISSUES, WM_TISSUES, GM_TISSUES, CSF_TISSUES]

    WM_TISSUES_ASEG = (
        WM_LABELS_ASEG + WM_CEREBELLUM_ASEG + BRAINSTEM_ASEG + PALLIDUM_ASEG
    )
    GM_TISSUES_ASEG = (
        GM_LABELS_ASEG
        + GM_CEREBELLUM_ASEG
        + CAUDATE_AND_ACCUMBENS_ASEG
        + PUTAMEN_ASEG
        + HIPPOCAMPUS_ASEG
        + THALAMUS_ASEG
        + HYPOTHALAMUS_ASEG
        + AMYGDALA_ASEG
        + CORTEX_ASEG
        + CLAUSTRUM_ASEG
    )
    CSF_TISSUES_ASEG = VENT_LABELS_ASEG
    BG_TISSUES_ASEG = BG_LABELS_ASEG
    TISSUES_ASEG = [BG_TISSUES_ASEG, WM_TISSUES_ASEG, GM_TISSUES_ASEG, CSF_TISSUES_ASEG]
elif False:  # merge background and CSF, but split cerebellum classes( 5D atlas)
    WM_TISSUES = WM_LABELS + BRAINSTEM + PALLIDUM
    WM_CEREBELLUM_TISSUES = WM_CEREBELLUM
    GM_TISSUES = (
        GM_LABELS
        + CAUDATE_AND_ACCUMBENS
        + PUTAMEN
        + HIPPOCAMPUS
        + THALAMUS
        + HYPOTHALAMUS
        + AMYGDALA
        + CORTEX
        + CLAUSTRUM
    )
    GM_CEREBELLUM_TISSUES = GM_CEREBELLUM
    BG_TISSUES = BG_LABELS + VENT_LABELS
    TISSUES = [
        BG_TISSUES,
        WM_TISSUES,
        GM_TISSUES,
        WM_CEREBELLUM_TISSUES,
        GM_CEREBELLUM_TISSUES,
    ]

    WM_TISSUES_ASEG = WM_LABELS_ASEG + BRAINSTEM_ASEG + PALLIDUM_ASEG
    WM_CEREBELLUM_TISSUES_ASEG = WM_CEREBELLUM_ASEG
    GM_TISSUES_ASEG = (
        GM_LABELS_ASEG
        + CAUDATE_AND_ACCUMBENS_ASEG
        + PUTAMEN_ASEG
        + HIPPOCAMPUS_ASEG
        + THALAMUS_ASEG
        + HYPOTHALAMUS_ASEG
        + AMYGDALA_ASEG
        + CORTEX_ASEG
        + CLAUSTRUM_ASEG
    )
    GM_CEREBELLUM_TISSUES_ASEG = GM_CEREBELLUM_ASEG
    BG_TISSUES_ASEG = BG_LABELS_ASEG + VENT_LABELS_ASEG
    TISSUES_ASEG = [
        BG_TISSUES_ASEG,
        WM_TISSUES_ASEG,
        GM_TISSUES_ASEG,
        WM_CEREBELLUM_TISSUES_ASEG,
        GM_CEREBELLUM_TISSUES_ASEG,
    ]
else:
    WM_TISSUES = WM_LABELS + BRAINSTEM + PALLIDUM
    WM_CEREBELLUM_TISSUES = WM_CEREBELLUM
    GM_TISSUES = GM_LABELS + HIPPOCAMPUS + HYPOTHALAMUS + AMYGDALA + CORTEX
    CAU_ACC_TISSUES = CAUDATE_AND_ACCUMBENS
    PUT_TISSUES = PUTAMEN + CLAUSTRUM
    THAL_TISSUES = THALAMUS
    GM_CEREBELLUM_TISSUES = GM_CEREBELLUM
    BG_TISSUES = BG_LABELS + VENT_LABELS
    TISSUES = [
        BG_TISSUES,
        WM_TISSUES,
        GM_TISSUES,
        WM_CEREBELLUM_TISSUES,
        GM_CEREBELLUM_TISSUES,
        CAU_ACC_TISSUES,
        PUT_TISSUES,
        THAL_TISSUES,
    ]

    TISSUE_GROUPS = [
        {"Background": BG_TISSUES},
        {"White Matter": WM_TISSUES},
        {"Gray Matter": GM_TISSUES},
        {"Cerebellum White Matter": WM_CEREBELLUM_TISSUES},
        {"Cerebellum Gray Matter": GM_CEREBELLUM_TISSUES},
        {"Caudate and Accumbens": CAU_ACC_TISSUES},
        {"Putamen": PUT_TISSUES},
        {"Thalamus": THAL_TISSUES},
    ]

    WM_TISSUES_ASEG = WM_LABELS_ASEG + BRAINSTEM_ASEG + PALLIDUM_ASEG
    WM_CEREBELLUM_TISSUES_ASEG = WM_CEREBELLUM_ASEG
    GM_TISSUES_ASEG = (
        GM_LABELS_ASEG
        + HIPPOCAMPUS_ASEG
        + HYPOTHALAMUS_ASEG
        + AMYGDALA_ASEG
        + CORTEX_ASEG
    )
    CAU_ACC_TISSUES_ASEG = CAUDATE_AND_ACCUMBENS_ASEG
    PUT_TISSUES_ASEG = PUTAMEN_ASEG + CLAUSTRUM_ASEG
    THAL_TISSUES_ASEG = THALAMUS_ASEG
    GM_CEREBELLUM_TISSUES_ASEG = GM_CEREBELLUM_ASEG
    BG_TISSUES_ASEG = BG_LABELS_ASEG + VENT_LABELS_ASEG
    TISSUES_ASEG = [
        BG_TISSUES_ASEG,
        WM_TISSUES_ASEG,
        GM_TISSUES_ASEG,
        WM_CEREBELLUM_TISSUES_ASEG,
        GM_CEREBELLUM_TISSUES_ASEG,
        CAU_ACC_TISSUES_ASEG,
        PUT_TISSUES_ASEG,
        THAL_TISSUES_ASEG,
    ]

    TISSUE_GROUPS_ASEG = [
        {"Background": BG_TISSUES_ASEG},
        {"White Matter": WM_TISSUES_ASEG},
        {"Gray Matter": GM_TISSUES_ASEG},
        {"Cerebellum White Matter": WM_CEREBELLUM_TISSUES_ASEG},
        {"Cerebellum Gray Matter": GM_CEREBELLUM_TISSUES_ASEG},
        {"Caudate and Accumbens": CAU_ACC_TISSUES_ASEG},
        {"Putamen": PUT_TISSUES_ASEG},
        {"Thalamus": THAL_TISSUES_ASEG},
    ]


def assign_labels_in_atlas(atlas_labels):

    tissue_index = [-1] * len(atlas_labels)

    for it, T in enumerate(TISSUES):
        for ia, l in enumerate(atlas_labels):
            if l in T:
                tissue_index[ia] = it

    return tissue_index


def assign_labels_in_aseg(aseg_labels):

    tissue_index = [-1] * len(aseg_labels)

    for it, T in enumerate(TISSUES_ASEG):
        for ia, l in enumerate(aseg_labels):
            if l in T:
                tissue_index[ia] = it

    return tissue_index


def _get_allen_lut_names(allen_lut_file, labels):

    atlas_names_and_labels = []
    with open(allen_lut_file, "r") as f:
        line_tmp = f.readline()
        while line_tmp:
            split_line = line_tmp.split()
            label = int(split_line[0])
            name = split_line[1]
            if label in labels:
                atlas_names_and_labels.append({name: label})

            line_tmp = f.readline()

    return atlas_names_and_labels


def create_atlas_settings(allen_lut_file, label_array_file, tissue_groups, aseg_groups):

    labels = np.load(label_array_file).tolist()
    labels = [int(element) for label in labels for element in label]
    atlas_names_and_labels = _get_allen_lut_names(allen_lut_file, labels)

    dict_tmp = {"Atlas Names and Labels": atlas_names_and_labels}
    with open("atlas_names_and_labels.yaml", "w") as file:
        yaml.dump(dict_tmp, file)

    dict_tmp = {"Aseg Combined Structures": aseg_groups}
    with open("combined_aseg_labels.yaml", "w") as file:
        yaml.dump(dict_tmp, file)

    dict_tmp = {"Combined Structures": tissue_groups}
    with open("combined_atlas_labels.yaml", "w") as file:
        yaml.dump(dict_tmp, file)


def get_tissue_settings(
    atlas_label_file,
    atlas_combined_label_file,
    aseg_combined_label_file,
    gmm_component_file,
    aseg_label_list,
):

    with open(atlas_label_file) as f:
        atlas_labels = yaml.full_load(f)

    # Flatten the dict of all atlas labels
    all_labels = atlas_labels["Atlas Names and Labels"]
    all_labels_list = []
    for lab_dict in all_labels:
        for k in lab_dict.keys():
            all_labels_list.append(lab_dict[k])

    with open(atlas_combined_label_file) as f:
        atlas_combined_labels = yaml.full_load(f)

    atlas_groups = atlas_combined_labels["Combined Structures"]
    atlas_indices = [-1] * len(all_labels_list)

    for it, T in enumerate(atlas_groups):
        for ia, l in enumerate(all_labels_list):
            if l in list(T.values())[0]:
                atlas_indices[ia] = it

    with open(aseg_combined_label_file) as f:
        aseg_combined_labels = yaml.full_load(f)

    aseg_groups = aseg_combined_labels["Aseg Combined Structures"]
    # aseg_indices = [-1] * len(aseg_label_list)
    aseg_indices = []
    for it, T in enumerate(aseg_groups):
        list_of_labels = []
        for ia, l in enumerate(aseg_label_list):
            if l in list(T.values())[0]:
                list_of_labels.append(int(l))
        aseg_indices.append(list_of_labels)

    with open(gmm_component_file) as f:
        gmm_components = yaml.full_load(f)

    gmm_components = gmm_components["Number of components in a class"]
    gmm_components =[list(x.values())[0] for x in gmm_components]

    return atlas_indices, aseg_indices, np.array(all_labels_list), np.array(gmm_components)
