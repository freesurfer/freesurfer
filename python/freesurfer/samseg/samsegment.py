import os

from .utilities import Specification
from .samsegment_part1 import samsegment_part1
from .samsegment_part2 import samsegment_part2
from .samsegment_part3 import samsegment_part3

# from .dev_utils.debug_client import run_test_cases, create_checkpoint_manager, load_starting_fixture


def samsegment(
    imageFileNames,
    transformedTemplateFileName,
    modelSpecifications,
    optimizationOptions,
    savePath,
    visualizer,
    checkpoint_manager=None
):

    # ------ Setup ------

    part0_results_dict = Specification({
        'imageFileNames': imageFileNames,
        'transformedTemplateFileName': transformedTemplateFileName,
        'modelSpecifications': modelSpecifications,
        'optimizationOptions': optimizationOptions,
        'savePath': savePath,
        'visualizer': visualizer,
    })

    if checkpoint_manager:
        checkpoint_manager.save_specification(part0_results_dict, 'part0', 1)

    # ------ Samsegment Part 1 ------

    print('calling part1...')
    part1_results_dict = samsegment_part1(
        imageFileNames,
        transformedTemplateFileName,
        modelSpecifications,
        optimizationOptions,
        savePath,
        visualizer,
        checkpoint_manager
    )

    if checkpoint_manager:
        checkpoint_manager.save(part1_results_dict, 'part1', 1)
    
    # ------ Samsegment Part 2 ------

    print('calling part2...')
    part2_results_dict = samsegment_part2(
        modelSpecifications,
        optimizationOptions,
        part1_results_dict,
        visualizer,
        checkpoint_manager
    )

    if checkpoint_manager:
        checkpoint_manager.save(part2_results_dict, 'part2', 1)
    
    # ------ Samsegment Part 3 ------

    print('calling part3...')
    part3_results_dict = samsegment_part3(
        modelSpecifications,
        optimizationOptions,
        part1_results_dict,
        part2_results_dict,
        imageFileNames,
        visualizer,
        checkpoint_manager
    )

    if checkpoint_manager:
        checkpoint_manager.save(part3_results_dict, 'part3', 1)

    # ------ Done ------

    print('...all parts completed')
    names = part1_results_dict['names']
    FreeSurferLabels = part3_results_dict['FreeSurferLabels']
    volumesInCubicMm = part3_results_dict['volumesInCubicMm']

    with open(os.path.join(savePath, 'samseg.stats'), 'w') as fid:
        for volume, name in zip(volumesInCubicMm, names):
            fid.write('# Measure {}, {:.6f}, mm^3\n'.format(name, volume))
    
    return [FreeSurferLabels, names, volumesInCubicMm]


def test_samsegment(case_name, case_file_folder, savePath):
    checkpoint_manager = create_checkpoint_manager(case_file_folder)
    fixture = load_starting_fixture()
    results = samsegment(
        fixture['imageFileNames'],
        fixture['transformedTemplateFileName'],
        fixture['modelSpecifications'],
        fixture['optimizationOptions'],
        savePath,
        fixture['visualizer'],
        checkpoint_manager
    )
    print(results)


if __name__ == '__main__':
    run_test_cases(action=test_samsegment)
