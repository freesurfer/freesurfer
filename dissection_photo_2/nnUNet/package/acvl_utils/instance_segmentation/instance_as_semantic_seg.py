from time import time
from multiprocessing import Pool
from typing import Union, Tuple, List

import SimpleITK
import numpy as np

from acvl_utils.cropping_and_padding.bounding_boxes import pad_bbox, bounding_box_to_slice
from acvl_utils.morphology.morphology_helper import generate_ball, label_with_component_sizes
from skimage.measure import label
from skimage.measure import regionprops
from skimage.morphology import ball, erosion, binary_dilation
from skimage.morphology import binary_erosion
from skimage.morphology import dilation
from acvl_utils.cropping_and_padding.bounding_boxes import regionprops_bbox_to_proper_bbox

BORDER_LABEL = 2
CENTER_LABEL = 1


def _internal_convert_semantic_to_instance_mp(cropped_core_instances, cropped_border, spacing):
    cropped_current = np.copy(cropped_core_instances)
    already_dilated_mm = np.array((0, 0, 0))
    cropped_final = np.copy(cropped_core_instances)
    background_mask = (cropped_border == 0) & (cropped_core_instances == 0)
    while np.sum(cropped_border) > 0:
        strel_size = [0, 0, 0]
        maximum_dilation = max(already_dilated_mm)
        for i in range(3):
            if spacing[i] == min(spacing):
                strel_size[i] = 1
                continue
            if already_dilated_mm[i] + spacing[i] / 2 < maximum_dilation:
                strel_size[i] = 1
        ball_here = ball(1)

        if strel_size[0] == 0: ball_here = ball_here[1:2]
        if strel_size[1] == 0: ball_here = ball_here[:, 1:2]
        if strel_size[2] == 0: ball_here = ball_here[:, :, 1:2]

        # print(1)
        dilated = dilation(cropped_current, ball_here)
        dilated[background_mask] = 0
        diff = (cropped_current == 0) & (dilated != cropped_current)
        cropped_final[diff & cropped_border] = dilated[diff & cropped_border]
        cropped_border[diff] = 0
        cropped_current = dilated
        already_dilated_mm = [already_dilated_mm[i] + spacing[i] if strel_size[i] == 1 else
                              already_dilated_mm[i] for i in range(3)]
    return cropped_final


def convert_semantic_to_instanceseg_mp(arr: np.ndarray,
                                       spacing: Union[Tuple[float, ...], List[float]] = (1, 1, 1),
                                       small_center_threshold: float = 30,
                                       isolated_border_as_separate_instance_threshold: float = 15,
                                       num_processes: int = 8) -> np.ndarray:
    """
    This function runs multiprocessing within one image. This is only useful for when you have one (or very few)
    very large images. For all other use cases you get much more throughput by parallelizing on a per image level
    (use multiprocessing pool to parallelize calls to convert_semantic_to_instanceseg)

    :param arr:
    :param spacing:
    :param small_center_threshold: volume, as dictated by spacing! If your spacing is (2, 2, 2) then a
    small_center_threshold of 16 would correspond to 2 pixels!
    :param isolated_border_as_separate_instance_threshold: volume, as dictated by spacing! If your spacing is (2, 2, 2) then a
    isolated_border_as_separate_instance_threshold of 16 would correspond to 2 pixels!
    :return:
    """
    assert np.issubdtype(arr.dtype, np.unsignedinteger), 'instance_segmentation must be an array of type unsigned ' \
                                                         'integer (can be uint8, uint16 etc)'
    pool = Pool(num_processes)

    spacing = np.array(spacing)
    small_center_threshold_in_pixels = round(small_center_threshold / np.prod(spacing))
    isolated_border_as_separate_instance_threshold_in_pixels = round(isolated_border_as_separate_instance_threshold / np.prod(
        spacing))

    # we first identify centers that are too small and set them to be border. This should remove false positive instances
    labeled_image, component_sizes = label_with_component_sizes(arr == CENTER_LABEL, connectivity=1)
    remove = np.array([i for i, j in component_sizes.items() if j <= small_center_threshold_in_pixels])
    remove = np.in1d(labeled_image.ravel(), remove).reshape(labeled_image.shape)
    arr[remove] = BORDER_LABEL

    # recompute core labels
    core_instances = label(arr == CENTER_LABEL, connectivity=1)

    # prepare empty array for results. uint32 for allowing HEAPS of instances. Otherwise we would cap at 65535
    final = np.zeros_like(core_instances, dtype=np.uint32)

    # besides the core instances we will need the borders
    border_mask = arr == BORDER_LABEL

    # we search for connected components and then convert each connected component into instance segmentation. This should
    # prevent bleeding into neighboring instances even if the instances don't touch
    connected_components, num_components = label(arr > 0, return_num=True, connectivity=1)
    max_component_idx = np.max(core_instances)
    rp = regionprops(connected_components)

    mp_results = []
    masks = []
    slicers = []
    for r in rp:
        bbox = regionprops_bbox_to_proper_bbox(r.bbox)
        slicer = bounding_box_to_slice(bbox)

        cropped_mask = connected_components[slicer] == r.label
        cropped_core_instances = np.copy(core_instances[slicer])
        cropped_border = np.copy(border_mask[slicer])

        # remove other objects from the current crop, only keep the current connected component
        cropped_core_instances[~cropped_mask] = 0
        cropped_border[~cropped_mask] = 0

        unique_core_idx = np.unique(cropped_core_instances)
        # do not use len(unique_core_idx) == 1 because there could be one code filling the entire thing
        if np.sum(cropped_core_instances) == 0:
            # special case no core
            if np.sum(cropped_border) > isolated_border_as_separate_instance_threshold_in_pixels:
                final[slicer][cropped_border] = max_component_idx + 1
                max_component_idx += 1
        elif len(unique_core_idx) == 2:
            # special case only one core = one object
            final[slicer][(cropped_core_instances > 0) | cropped_border] = unique_core_idx[1]
        else:
            mp_results.append(
                pool.starmap_async(
                    _internal_convert_semantic_to_instance_mp,
                    ((cropped_core_instances, cropped_border, spacing),)
                )
            )
            masks.append(cropped_mask)
            slicers.append(slicer)

    for mp_r, mask, slc in zip(mp_results, masks, slicers):
        cropped_final = mp_r.get()[0]
        # now place result back
        final[slc][mask] = cropped_final[mask]

    pool.close()
    pool.join()
    return final


def convert_semantic_to_instanceseg(arr: np.ndarray,
                                    spacing: Union[Tuple[float, ...], List[float]] = (1, 1, 1),
                                    small_center_threshold: float = 30,
                                    isolated_border_as_separate_instance_threshold: float = 15) -> np.ndarray:
    """
    :param arr:
    :param spacing:
    :param small_center_threshold: volume, as dictated by spacing! If your spacing is (2, 2, 2) then a
    small_center_threshold of 16 would correspond to 2 pixels!
    :param isolated_border_as_separate_instance_threshold: volume, as dictated by spacing! If your spacing is (2, 2, 2) then a
    isolated_border_as_separate_instance_threshold of 16 would correspond to 2 pixels!
    :return:
    """
    assert np.issubdtype(arr.dtype, np.unsignedinteger), 'instance_segmentation must be an array of type unsigned ' \
                                                         'integer (can be uint8, uint16 etc)'
    spacing = np.array(spacing)
    small_center_threshold_in_pixels = round(small_center_threshold / np.prod(spacing))
    isolated_border_as_separate_instance_threshold_in_pixels = round(isolated_border_as_separate_instance_threshold / np.prod(
        spacing))

    # we first identify centers that are too small and set them to be border. This should remove false positive instances
    labeled_image, component_sizes = label_with_component_sizes(arr == CENTER_LABEL, connectivity=1)
    remove = np.array([i for i, j in component_sizes.items() if j <= small_center_threshold_in_pixels])
    remove = np.in1d(labeled_image.ravel(), remove).reshape(labeled_image.shape)
    arr[remove] = BORDER_LABEL

    # recompute core labels
    core_instances = label(arr == CENTER_LABEL, connectivity=1)

    # prepare empty array for results
    final = np.zeros_like(core_instances, dtype=np.uint32)

    # besides the core instances we will need the borders
    border_mask = arr == BORDER_LABEL

    # we search for connected components and then convert each connected component into instance segmentation. This should
    # prevent bleeding into neighboring instances even if the instances don't touch
    connected_components, num_components = label(arr > 0, return_num=True, connectivity=1)
    max_component_idx = np.max(core_instances)
    rp = regionprops(connected_components)

    for r in rp:
        bbox = regionprops_bbox_to_proper_bbox(r.bbox)
        slicer = bounding_box_to_slice(bbox)

        cropped_mask = connected_components[slicer] == r.label
        cropped_core_instances = np.copy(core_instances[slicer])
        cropped_border = np.copy(border_mask[slicer])

        # remove other objects from the current crop, only keep the current connected component
        cropped_core_instances[~cropped_mask] = 0
        cropped_border[~cropped_mask] = 0
        cropped_current = np.copy(cropped_core_instances)

        unique_core_idx = np.unique(cropped_core_instances)
        # do not use len(unique_core_idx) == 1 because there could be one code filling the entire thing
        if np.sum(cropped_core_instances) == 0:
            # special case no core
            if np.sum(cropped_border) > isolated_border_as_separate_instance_threshold_in_pixels:
                final[slicer][cropped_border] = max_component_idx + 1
                max_component_idx += 1
        elif len(unique_core_idx) == 2:
            # special case only one core = one object
            final[slicer][(cropped_core_instances > 0) | cropped_border] = unique_core_idx[1]
        else:
            already_dilated_mm = np.array((0, 0, 0))
            cropped_final = np.copy(cropped_core_instances)
            while np.sum(cropped_border) > 0:
                strel_size = [0, 0, 0]
                maximum_dilation = max(already_dilated_mm)
                for i in range(3):
                    if spacing[i] == min(spacing):
                        strel_size[i] = 1
                        continue
                    if already_dilated_mm[i] + spacing[i] / 2 < maximum_dilation:
                        strel_size[i] = 1
                ball_here = ball(1)

                if strel_size[0] == 0: ball_here = ball_here[1:2]
                if strel_size[1] == 0: ball_here = ball_here[:, 1:2]
                if strel_size[2] == 0: ball_here = ball_here[:, :, 1:2]

                # print(1)
                dilated = dilation(cropped_current, ball_here)
                dilated[~cropped_mask] = 0
                diff = (cropped_current == 0) & (dilated != cropped_current)
                cropped_final[diff & cropped_border] = dilated[diff & cropped_border]
                cropped_border[diff] = 0
                cropped_current = dilated
                already_dilated_mm = [
                    already_dilated_mm[i] + spacing[i] if strel_size[i] == 1 else already_dilated_mm[i] for i in
                    range(3)]

            # now place result back
            final[slicer][cropped_mask] = cropped_final[cropped_mask]
    return final


def _internal_postprocess_instance_segmentation_mp(instance_mask, cropped_instance):
    # let's see if this instance is fragmented
    strel = ball(1)
    fragment_masks = []
    target_labels = []

    labeled_cropped_instance, fragment_sizes = label_with_component_sizes(instance_mask, connectivity=1)
    if len(fragment_sizes) > 1:
        max_fragment_size = max(fragment_sizes.values())
        # print(instance_id, fragment_sizes)
        for f in fragment_sizes.keys():
            if fragment_sizes[f] == max_fragment_size: continue
            fragment_mask = labeled_cropped_instance == f
            neighbor_mask = binary_dilation(fragment_mask, strel) != fragment_mask
            fragment_neighbors = cropped_instance[neighbor_mask]
            # remove background
            fragment_neighbors = fragment_neighbors[fragment_neighbors != 0]
            if len(fragment_neighbors) > 0:
                counts = np.bincount(fragment_neighbors)
                # replace fragment with most instance it shares the largest border with
                fragment_masks.append(fragment_mask)
                target_labels = np.argmax(counts)
    return fragment_masks, target_labels


def postprocess_instance_segmentation_mp(instance_segmentation: np.ndarray, num_processes: int = 8):
    """
    Sometimes the dilation used to convert sem seg back to inst seg can cause fragmented instances. This is more of an artifact
    rather than real. This function can fix this by merging all but the largest connected component of each fragment
    with the nearest neighboring instances.

    This function runs multiprocessing within one image. This is only useful for when you have one (or very few)
    very large images. For all other use cases you get much more throughput by parallelizing on a per image level
    (use multiprocessing pool to parallelize calls to postprocess_instance_segmentation)

    :param instance_segmentation:
    :return:
    """
    assert np.issubdtype(instance_segmentation.dtype,
                         np.unsignedinteger), 'instance_segmentation must be an array of type unsigned ' \
                                              'integer (can be uint8, uint16 etc)'
    pool = Pool(num_processes)

    results = []
    slicers = []

    rp = regionprops(instance_segmentation)
    for r in rp:
        bbox = regionprops_bbox_to_proper_bbox(r.bbox)
        slicer = bounding_box_to_slice(bbox)

        cropped_instance = instance_segmentation[slicer]
        instance_mask = instance_segmentation[slicer] == r.label

        results.append(
            pool.starmap_async(
                _internal_postprocess_instance_segmentation_mp,
                ((instance_mask, cropped_instance),)
            )
        )
        slicers.append(slicer)

    for slicer, result in zip(slicers, results):
        for fragment_mask, target_label in zip(*result.get()[0]):
            instance_segmentation[slicer][fragment_mask] = target_label

    pool.close()
    pool.join()
    return instance_segmentation


def postprocess_instance_segmentation(instance_segmentation: np.ndarray):
    """
    Sometimes the dilation used to convert sem seg back to inst seg can cause fragmented instances. This is more of an artifact
    rather than real. This function can fix this by merging all but the largest connected component of each fragment
    with the nearest neighboring instances.

    :param instance_segmentation:
    :return:
    """
    assert np.issubdtype(instance_segmentation.dtype,
                         np.unsignedinteger), 'instance_segmentation must be an array of type unsigned ' \
                                              'integer (can be uint8, uint16 etc)'
    rp = regionprops(instance_segmentation)
    strel = ball(1)
    for r in rp:
        bbox = regionprops_bbox_to_proper_bbox(r.bbox)
        slicer = bounding_box_to_slice(bbox)

        cropped_instance = instance_segmentation[slicer]
        instance_mask = instance_segmentation[slicer] == r.label

        # let's see if this instance is fragmented
        labeled_cropped_instance, fragment_sizes = label_with_component_sizes(instance_mask, connectivity=1)
        if len(fragment_sizes) > 1:
            max_fragment_size = max(fragment_sizes.values())
            # print(instance_id, fragment_sizes)
            for f in fragment_sizes.keys():
                if fragment_sizes[f] == max_fragment_size: continue
                fragment_mask = labeled_cropped_instance == f
                neighbor_mask = binary_dilation(fragment_mask, strel) != fragment_mask
                fragment_neighbors = cropped_instance[neighbor_mask]
                # remove background
                fragment_neighbors = fragment_neighbors[fragment_neighbors != 0]
                if len(fragment_neighbors) > 0:
                    counts = np.bincount(fragment_neighbors)
                    # replace fragment with most instance it shares the largest border with
                    instance_segmentation[slicer][fragment_mask] = np.argmax(counts)
    return instance_segmentation


def convert_instanceseg_to_semantic(instance_segmentation: np.ndarray, spacing: Union[Tuple, List] = (1, 1, 1),
                                    border_thickness: float = 2) -> np.ndarray:
    assert np.issubdtype(instance_segmentation.dtype,
                         np.unsignedinteger), 'instance_segmentation must be an array of type unsigned ' \
                                              'integer (can be uint8, uint16 etc)'
    border_semantic = np.zeros_like(instance_segmentation, dtype=np.uint8)
    selem = generate_ball([border_thickness] * 3, spacing)
    labels = np.unique(instance_segmentation)
    for label in labels:
        if label == 0:
            continue
        mask = instance_segmentation == label
        eroded = erosion(mask, selem)
        border_semantic[(~eroded) & mask] = BORDER_LABEL
        border_semantic[eroded & mask] = CENTER_LABEL
    return border_semantic


def convert_instanceseg_to_semantic_patched(instance_segmentation: np.ndarray, spacing: Union[Tuple, List] = (1, 1, 1),
                                            border_thickness: float = 2) -> np.ndarray:
    assert np.issubdtype(instance_segmentation.dtype,
                         np.unsignedinteger), 'instance_segmentation must be an array of type unsigned ' \
                                              'integer (can be uint8, uint16 etc)'
    border_semantic = np.zeros_like(instance_segmentation, dtype=np.uint8)
    selem = generate_ball([border_thickness] * 3, spacing)
    pad_amount = 1  # testing purposes only, erosion should not need padding
    instance_properties = regionprops(instance_segmentation)
    for ip in instance_properties:
        bbox = regionprops_bbox_to_proper_bbox(ip['bbox'])
        if pad_amount != 0:
            bbox = pad_bbox(bbox, pad_amount, instance_segmentation.shape)
        slicer = bounding_box_to_slice(bbox)
        instance_cropped = instance_segmentation[slicer]
        instance_mask = instance_cropped == ip["label"]
        instance_mask_eroded = binary_erosion(instance_mask, selem)
        border_semantic[slicer][(~instance_mask_eroded) & instance_mask] = BORDER_LABEL
        border_semantic[slicer][instance_mask_eroded & instance_mask] = CENTER_LABEL
    return border_semantic


def _internal_convert_instanceseg_to_semantic_patched_mp(instance_int: int, selem: np.ndarray, cropped_is: np.ndarray):
    instance_mask = cropped_is == instance_int
    result = np.zeros(instance_mask.shape, dtype=np.uint8)
    instance_mask_eroded = binary_erosion(instance_mask, selem)
    result[(~instance_mask_eroded) & instance_mask] = BORDER_LABEL
    result[instance_mask_eroded & instance_mask] = CENTER_LABEL
    return result, instance_mask


def convert_instanceseg_to_semantic_patched_mp(instance_segmentation: np.ndarray,
                                               spacing: Union[Tuple, List] = (1, 1, 1),
                                               border_thickness: float = 2,
                                               num_processes: int = 8) -> np.ndarray:
    """
    This function runs multiprocessing within one image. This is only useful for when you have one (or very few)
    very large images. For all other use cases you get much more throughput by parallelizing on a per image level
    (use multiprocessing pool to parallelize calls to convert_instanceseg_to_semantic_patched)

    :param instance_segmentation:
    :param spacing:
    :param border_thickness:
    :param num_processes:
    :return:
    """
    assert np.issubdtype(instance_segmentation.dtype,
                         np.unsignedinteger), 'instance_segmentation must be an array of type unsigned ' \
                                              'integer (can be uint8, uint16 etc)'
    pool = Pool(num_processes)
    border_semantic = np.zeros_like(instance_segmentation, dtype=np.uint8)
    selem = generate_ball([border_thickness] * 3, spacing)
    pad_amount = 0  # testing purposes only, erosion should not need padding
    instance_properties = regionprops(instance_segmentation)
    results = []
    slicers = []
    for ip in instance_properties:
        bbox = regionprops_bbox_to_proper_bbox(ip['bbox'])
        if pad_amount != 0:
            bbox = pad_bbox(bbox, pad_amount, instance_segmentation.shape)
        slicer = bounding_box_to_slice(bbox)
        instance_cropped = instance_segmentation[slicer]
        results.append(
            pool.starmap_async(
                _internal_convert_instanceseg_to_semantic_patched_mp,
                ((ip['label'], selem, instance_cropped),)
            )
        )
        slicers.append(slicer)

    for r, s in zip(results, slicers):
        semseg, instance_mask = r.get()[0]
        border_semantic[s][instance_mask] = semseg[instance_mask]

    pool.close()
    pool.join()
    return border_semantic


def main_instance_to_sem():
    source_file = '/media/isensee/data_stick/Task162_ShankTeeth/labelsTr/aarika-drunk-chimpanzee.nii.gz'
    target_normal = '/home/isensee/temp/aarika-drunk-chimpanzee_normal.nii.gz'
    target_patched = '/home/isensee/temp/aarika-drunk-chimpanzee_patched.nii.gz'
    target_patched_mp = '/home/isensee/temp/aarika-drunk-chimpanzee_patched_mp.nii.gz'
    import SimpleITK as sitk
    source_file_itk = sitk.ReadImage(source_file)
    spacing = list(source_file_itk.GetSpacing())[::-1]
    source_npy = sitk.GetArrayFromImage(source_file_itk).astype(np.uint8)
    source_npy[(source_npy <= 4) | (source_npy >= 37)] = 0  # only keep adult teeth

    start = time()
    semseg_patched = convert_instanceseg_to_semantic_patched(source_npy, spacing, border_thickness=1)
    print(f'convert_instanceseg_to_semantic_patched: {time() - start} s')
    semseg_itk = sitk.GetImageFromArray(semseg_patched)
    semseg_itk.CopyInformation(source_file_itk)
    sitk.WriteImage(semseg_itk, target_patched)

    start = time()
    semseg_patched_mp = convert_instanceseg_to_semantic_patched_mp(source_npy, spacing, border_thickness=1,
                                                                   num_processes=8)
    print(f'convert_instanceseg_to_semantic_patched_mp: {time() - start} s')
    semseg_itk = sitk.GetImageFromArray(semseg_patched_mp)
    semseg_itk.CopyInformation(source_file_itk)
    sitk.WriteImage(semseg_itk, target_patched_mp)

    start = time()
    semseg = convert_instanceseg_to_semantic(source_npy, spacing, border_thickness=1)
    print(f'convert_instanceseg_to_semantic: {time() - start} s')
    semseg_itk = sitk.GetImageFromArray(semseg)
    semseg_itk.CopyInformation(source_file_itk)
    sitk.WriteImage(semseg_itk, target_normal)


def main_sem_to_instance():
    source_file = '/home/isensee/temp/aarika-drunk-chimpanzee_patched.nii.gz'
    target_patched = '/home/isensee/temp/aarika-drunk-chimpanzee_patched_inst.nii.gz'
    target_patched_mp = '/home/isensee/temp/aarika-drunk-chimpanzee_patched_inst_mp.nii.gz'
    import SimpleITK as sitk
    source_file_itk = sitk.ReadImage(source_file)
    spacing = list(source_file_itk.GetSpacing())[::-1]
    source_npy = sitk.GetArrayFromImage(source_file_itk)

    start = time()
    instseg_patched = convert_semantic_to_instanceseg(source_npy, spacing, small_center_threshold=5,
                                                      isolated_border_as_separate_instance_threshold=5)
    print(f'convert_semantic_to_instanceseg: {time() - start} s')
    start = time()
    instseg_patched_postprocessed = postprocess_instance_segmentation(instseg_patched)
    print(f'postprocess_instance_segmentation: {time() - start} s')

    instseg_itk = sitk.GetImageFromArray(instseg_patched_postprocessed)
    instseg_itk.CopyInformation(source_file_itk)
    sitk.WriteImage(instseg_itk, target_patched)

    start = time()
    instseg_patched_mp = convert_semantic_to_instanceseg_mp(source_npy, spacing, small_center_threshold=5,
                                                            isolated_border_as_separate_instance_threshold=5,
                                                            num_processes=8)
    print(f'convert_semantic_to_instanceseg_mp: {time() - start} s')
    start = time()
    instseg_patched_postprocessed_mp = postprocess_instance_segmentation_mp(instseg_patched_mp, num_processes=8)
    print(f'postprocess_instance_segmentation_mp: {time() - start} s')
    instseg_itk = sitk.GetImageFromArray(instseg_patched_postprocessed_mp)
    instseg_itk.CopyInformation(source_file_itk)
    sitk.WriteImage(instseg_itk, target_patched_mp)


if __name__ == '__main__':
    main_sem_to_instance()

