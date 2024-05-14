from typing import Tuple, Union

import numpy as np
from scipy.ndimage import distance_transform_edt, binary_dilation, binary_erosion
from skimage.morphology import ball, disk


def signed_edt(binary_input_array: np.ndarray, spacing: Union[float, Tuple[float, ...]] = 1.):
    mask = binary_input_array == 1

    if len(binary_input_array.shape) == 3:
        strel = ball(1)
    elif len(binary_input_array.shape) == 2:
        strel = disk(1)
    else:
        raise RuntimeError()

    contour = binary_dilation(np.copy(binary_input_array), strel) - binary_input_array
    distance_to_contour_dil = distance_transform_edt(1 - contour, spacing)
    distance_to_contour_dil[mask] = - distance_to_contour_dil[mask]

    contour = binary_input_array - binary_erosion(np.copy(binary_input_array), strel)
    distance_to_contour_ero = distance_transform_edt(1 - contour, spacing)
    distance_to_contour_ero[mask] = - distance_to_contour_ero[mask]

    return (distance_to_contour_dil + distance_to_contour_ero) / 2


def slice_interpolation_axial(binary_input_segmentation: np.ndarray):
    sums = np.sum(binary_input_segmentation, axis=(1, 2))
    slices_with_labels = np.where(sums)[0]
    interpolated_seg = np.zeros_like(binary_input_segmentation)
    lower = slices_with_labels[0]
    upper_idx = 1
    upper = slices_with_labels[upper_idx]
    edt_lower = signed_edt(binary_input_segmentation[lower])
    edt_upper = signed_edt(binary_input_segmentation[upper])
    for slice_id in range(min(slices_with_labels), max(slices_with_labels) + 1):
        if sums[slice_id] == 0:
            factor_lower = 1 - (slice_id - lower) / (upper - lower)
            factor_upper = 1 - (upper - slice_id) / (upper - lower)
            assert factor_upper + factor_lower == 1
            print(slice_id, lower, upper, factor_lower, factor_upper)
            interpolated_seg[slice_id] = (factor_lower * edt_lower + factor_upper * edt_upper) <= 0
        else:
            interpolated_seg[slice_id] = binary_input_segmentation[slice_id]
            if slice_id != lower:
                assert slice_id == upper
                lower = upper
                if lower == max(slices_with_labels):
                    break
                upper_idx += 1
                upper = slices_with_labels[upper_idx]
                edt_lower = edt_upper
                edt_upper = signed_edt(binary_input_segmentation[upper])
    return interpolated_seg


if __name__ == '__main__':
    import SimpleITK as sitk
    a = sitk.GetArrayFromImage(sitk.ReadImage('/media/fabian/data/nnUNet_raw/Dataset990_Daniel/imagesVal/20220330_154623_maize_bg_stack_3_80KeV_1mm_Al_height_138mm_tray_0_ix_2_iy_2_0000-labels.nrrd'))
    b = a.transpose((2, 0, 1))
    c = slice_interpolation_axial(b)
    d = c.transpose((1, 2, 0))
    sitk.WriteImage(sitk.GetImageFromArray(d), '/media/fabian/data/nnUNet_raw/Dataset990_Daniel/imagesVal/20220330_154623_maize_bg_stack_3_80KeV_1mm_Al_height_138mm_tray_0_ix_2_iy_2.nii.gz')