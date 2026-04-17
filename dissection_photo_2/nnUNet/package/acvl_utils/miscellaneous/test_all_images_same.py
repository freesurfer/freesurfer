from typing import Callable

from batchgenerators.utilities.file_and_folder_operations import subfiles, join
from multiprocessing import Pool

import SimpleITK as sitk
import numpy as np


def test_same_generic(file_1: str, file_2: str, load_fn: Callable[[str], np.ndarray], verbose: bool = False):
    """
    reads two files with SimpleITK and returns whether theis pixelarrays have the same content
    :param file_1:
    :param file_2:
    :param load_fn: Callable, must take a string and return a numpy array
    :param verbose:
    :return:
    """
    file_1_np = load_fn(file_1)
    file_2_np = load_fn(file_2)

    if not len(file_1_np.shape) == len(file_2_np.shape):
        print(f'Dimension mismatch between {file_1} and {file_2}. Shapes {file_1_np.shape} and {file_2_np.shape}')
        return

    if not all([i == j for i, j in zip(file_1_np.shape, file_2_np.shape)]):
        print(f'Shape mismatch between {file_1} and {file_2}. Shapes {file_1_np.shape} and {file_2_np.shape}')
        return

    same = np.all(file_1_np == file_2_np)
    if not same:
        print(f'{file_1} and {file_2} do not have the same pixel array :-(')
    elif verbose:
        print(f'{file_1} and {file_2} OK!')


def test_all_images_in_folders_same_mp(folder_1: str, folder_2: str, load_fn: Callable[[str], np.ndarray],
                                       suffix: str = '.nii.gz', num_processes: int = 8,
                                       verbose: bool = False, strict: bool = True):
    """

    :param folder_1:
    :param folder_2:
    :param suffix:
    :param num_processes:
    :param verbose:
    :param strict: if True it will crash if not all files from folder 1 are in folder 2 and vice versa. If False it
    will just evaluate the common files
    :return:
    """
    pool = Pool(num_processes)
    files_1 = subfiles(folder_1, suffix=suffix, join=False)
    files_2 = subfiles(folder_2, suffix=suffix, join=False)
    if strict:
        assert all([i in files_2 for i in files_1]) and all([i in files_1 for i in files_2]), \
            "Not all files from folder_1 are in folder_2 or the other way around. Make sure this is fixed or " \
            "set strict=False"
        workon = files_1
    else:
        workon = [i for i in files_1 if i in files_2]
    if verbose:
        print(f'Found {len(workon)} files')
    res = pool.starmap_async(
        test_same_generic,
        zip(
            [join(folder_1, i) for i in workon],
            [join(folder_2, i) for i in workon],
            [load_fn] * len(workon),
            [verbose] * len(workon)
        )
    )
    _ = res.get()
    pool.close()
    pool.join()


def _load_sitk(file: str):
    return sitk.GetArrayFromImage(sitk.ReadImage(file))


def _load_npz_nnunet(file: str):
    return np.load(file)['data']


if __name__ == '__main__':
    base = './'
    test_all_images_in_folders_same_mp(join(base, 'nnUNetData_plans_v2.1_trgSp_05x05x05_stage0'), join(base, 'nnUNetData_plans_v2.1_trgSp_05x05x05_stage0_patched/'),
                                       _load_npz_nnunet, suffix='.npz',
                                       num_processes=8, verbose=False, strict=True)
