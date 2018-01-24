from as_python.samseg.dev_utils.compute_dice_measure import dice_compare_test_case
from as_python.samseg.dev_utils.debug_client import run_test_cases


def full_dice_compare(case_name, case_file_folder, savePath):
    dice_compare_test_case(case_name, case_file_folder, savePath, cross_compare=True)


if __name__ == '__main__':
    run_test_cases(action=full_dice_compare)
