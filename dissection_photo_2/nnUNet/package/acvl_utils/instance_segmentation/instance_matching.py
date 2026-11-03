from copy import deepcopy
from multiprocessing import Pool
from typing import List, Tuple, Callable

import numpy as np


def match_instances_dice(instances_gt: np.ndarray, instances_pred: np.ndarray, dice_cutoff: float = 0.1,
                         consume_instances: bool = True) -> List[Tuple]:
    """
    Takes a predicted and a ground truth instance maps and matches the instances between them via their Dice overlap.
    This is implemented by first computing the pair-wise Dice between all gt-pred instance combinations. Then, we
    iteratively take the pair with the highest Dice overlap and remove it from the set. If consume_instances=False,
    gt instances can be matched with multiple predicted instances and vice versa (Not recommended).

    This function returns a List of Tuples. Each entry in the list is a match, each match is represented as a tuple.
    This is what the matches look like:
    (instance_id_gt, instance_id_pred, dice, num_pixels_gt_instance, num_pixels_pred_instance)

    This is what an FP prediction looks like (instance in instances_pred that could not be matched):
    (None, instance_id_pred, 0, 0, num_pixels_pred_instance)

    This is what an FN prediction looks like (instance in instances_pred that could not be matched):
    (instance_id_gt, None, 0, num_pixels_gt_instance, 0)

    :param instances_gt: np.ndarray of some uint type (each integer value corresponds to one predicted instance)
    :param instances_pred: np.ndarray of some uint type (each integer value corresponds to one predicted instance)
    :param dice_cutoff: only instance pairs with an overlap > dice_cutoff will be considered for matching. All others
    are either FP or FN.
    :param consume_instances: If True, each gt/pred instance can only be paired with at most one pred/gt instance.
    :return:
    """
    assert np.issubdtype(instances_gt.dtype, np.unsignedinteger), "instances_gt must have some unsigned integer " \
                                                                  "dtype (uint8, utin16 etc)"
    assert np.issubdtype(instances_pred.dtype, np.unsignedinteger), "instances_pred must have some unsigned integer " \
                                                                    "dtype (uint8, utin16 etc)"
    instance_ids_gt = [i for i in np.unique(instances_gt) if i != 0]
    instance_ids_pred = [i for i in np.unique(instances_pred) if i != 0]

    pixel_mask = (instances_pred > 0) | (instances_gt > 0)
    pred_pixels = instances_pred[pixel_mask]
    gt_pixels = instances_gt[pixel_mask]

    # we match pred with gt instances based on Dice
    dice_matrix = np.zeros((len(instance_ids_gt), len(instance_ids_pred)))
    for i in range(len(instance_ids_gt)):
        instance_label = instance_ids_gt[i]
        mask_gt = gt_pixels == instance_label
        # we only need to check the dice for instance_ids_pred where there is any overlap with mask_gt
        instance_ids_pred_with_overlap = [i for i in np.unique(pred_pixels[mask_gt]) if i != 0]
        for j in range(len(instance_ids_pred)):
            instance_label_pred = instance_ids_pred[j]
            if instance_label_pred in instance_ids_pred_with_overlap:
                mask_pred = pred_pixels == instance_label_pred
                tp = np.sum(mask_gt & mask_pred)
                fp = np.sum((~mask_gt) & mask_pred)
                fn = np.sum(mask_gt & (~mask_pred))
                dc = 2 * tp / np.clip(2 * tp + fp + fn, a_min=1e-8, a_max=np.inf)
                dice_matrix[i, j] = dc

    remaining_gt = deepcopy(instance_ids_gt)
    remaining_pred = deepcopy(instance_ids_pred)
    matches = []
    while np.any(dice_matrix > dice_cutoff):
        gt_ind, pred_ind = np.unravel_index(np.argmax(dice_matrix), dice_matrix.shape)
        label_gt = instance_ids_gt[gt_ind]
        label_pred = instance_ids_pred[pred_ind]

        dice = dice_matrix[gt_ind, pred_ind]
        matches.append((label_gt, label_pred, dice, np.sum(gt_pixels == label_gt), np.sum(pred_pixels == label_pred)))

        if consume_instances:
            dice_matrix[gt_ind, :] = 0
            dice_matrix[:, pred_ind] = 0
        else:
            dice_matrix[gt_ind, pred_ind] = 0

        if label_gt in remaining_gt:
            remaining_gt.remove(label_gt)
        if label_pred in remaining_pred:
            remaining_pred.remove(label_pred)

    for lg in remaining_gt:
        matches.append((lg, None, 0, np.sum(gt_pixels == lg), 0))

    for lp in remaining_pred:
        matches.append((None, lp, 0, 0, np.sum(pred_pixels == lp)))

    return matches


def compute_all_matches(files_gt: List[str], files_pred: List[str], load_fn: Callable[[str], np.ndarray],
                        dice_cutoff: float = 0.1, consume_instances: bool = True,
                        num_processes: int = 6) -> List[Tuple]:
    """
    Wrapper for computing matches between instance segmentations using a multiprocessing.Pool. See match_instances_dice
    for documentation regarding dice_cutoff and consume_instances

    files_gt and files_pred msut have the same length and must be sorted properly (so that files_pred[i] and files_gt[i]
    are the correct pair)

    :param files_gt: list of files containing reference instance segmentations
    :param files_pred: list of files containing predicted instance segmentations
    :param load_fn: a callable that can load a file and returns a numpy ndarray of some uint type
    :param dice_cutoff: See match_instances_dice
    :param consume_instances: See match_instances_dice
    :param num_processes: number of CPU processes used for matching
    :return:
    """
    p = Pool(num_processes)
    results = []
    for file_gt, file_pred in zip(files_gt, files_pred):
        results.append(p.starmap_async(
            _load_compute_matches, ((file_gt, file_pred, load_fn, dice_cutoff, consume_instances), )
        ))

    results = [i.get()[0] for i in results]
    p.close()
    p.join()
    return results


def _load_compute_matches(image_gt: str, image_pred: str, load_fn, dice_cutoff, consume_instances):
    gt = load_fn(image_gt)
    pred = load_fn(image_pred)
    return match_instances_dice(gt, pred, dice_cutoff, consume_instances)
