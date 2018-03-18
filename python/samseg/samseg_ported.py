import logging
import os.path

from samseg.dev_utils.debug_client import run_test_cases, create_checkpoint_manager, load_starting_fixture
from samseg.run_utilities import Specification
from samseg.samseg_ported_part1 import samsegment_part1
from samseg.samseg_ported_part2 import samsegment_part2
from samseg.samseg_ported_part3 import samsegment_part3

logger = logging.getLogger(__name__)

def samsegment(
        imageFileNames,
        transformedTemplateFileName,
        modelSpecifications,
        optimizationOptions,
        savePath,
        showFigures,
        checkpoint_manager=None
):
    part0_results_dict = Specification({
        'imageFileNames': imageFileNames,
        'transformedTemplateFileName': transformedTemplateFileName,
        'modelSpecifications': modelSpecifications,
        'optimizationOptions': optimizationOptions,
        'savePath': savePath,
        'showFigures': showFigures,
    })
    if checkpoint_manager:
        checkpoint_manager.save_specification(part0_results_dict, 'part0', 1)
    logger.info('calling part1...')
    part1_results_dict = samsegment_part1(
        imageFileNames,
        transformedTemplateFileName,
        modelSpecifications,
        optimizationOptions,
        savePath,
        showFigures,
        checkpoint_manager
    )
    if checkpoint_manager:
        checkpoint_manager.save(part1_results_dict, 'part1', 1)
    logger.info('calling part2...')

    part2_results_dict = samsegment_part2(
        modelSpecifications,
        optimizationOptions,
        part1_results_dict,
        checkpoint_manager
    )
    if checkpoint_manager:
        checkpoint_manager.save(part2_results_dict, 'part2', 1)

    logger.info('calling part3...')
    part3_results_dict = samsegment_part3(
        modelSpecifications,
        optimizationOptions,
        part1_results_dict,
        part2_results_dict,
        imageFileNames,
        checkpoint_manager
    )
    if checkpoint_manager:
        checkpoint_manager.save(part3_results_dict, 'part3', 1)
    logger.info('...all parts completed')
    names = part1_results_dict['names']
    FreeSurferLabels = part3_results_dict['FreeSurferLabels']
    volumesInCubicMm = part3_results_dict['volumesInCubicMm']

    volfp = open(os.path.join(savePath, 'samseg.stats'), 'w')
    for i, names in enumerate(names):
        vol =  "%.6f" % volumesInCubicMm[i]
        volfp.write('# Measure ' + names + ', ' + vol + ', mm^3\n')
    volfp.close()
    
    return [FreeSurferLabels, names, volumesInCubicMm]


def test_samseg_ported(case_name, case_file_folder, savePath):
    checkpoint_manager = create_checkpoint_manager(case_file_folder)
    fixture = load_starting_fixture()
    results = samsegment(
        fixture['imageFileNames'],
        fixture['transformedTemplateFileName'],
        fixture['modelSpecifications'],
        fixture['optimizationOptions'],
        savePath,
        fixture['showFigures'],
        checkpoint_manager
    )
    print(results)


if __name__ == '__main__':
    run_test_cases(action=test_samseg_ported)
