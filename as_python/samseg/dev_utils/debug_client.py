'''
USAGE: First start up the debug server in MATLAB. Please see debug_server.m for instructions
Then call `request_var` and pass in the name of the variable you want to auto-magically
transfer from MATLAB to Python. This is intended to be called with both Python and MATLAB
at a breakpoint to easily compare the values of variables for equality.
'''
import json
import logging
import os
import sys
import time
import traceback
from itertools import permutations

import numpy as np
import scipy.io

logger = logging.getLogger(__name__)

# MATLAB_DUMP_DIR = "/Users/ys/work/freesurfer/matlab_dumps"
# MATLAB_DUMP_DIR = '/home/willy/work/cm/my_tests/matlab_dump_dir'
MATLAB_DUMP_DIR = os.getenv('MATLAB_DUMP_DIR')


def request_var(varname):
    mat_filename = os.path.join(MATLAB_DUMP_DIR, varname + '.mat')
    if os.path.exists(mat_filename):
        os.remove(mat_filename)
    request_file = os.path.join(MATLAB_DUMP_DIR, varname + '.request')
    open(request_file, 'w')
    for i in range(200):
        if not os.path.exists(request_file):
            var = scipy.io.loadmat(mat_filename, struct_as_record=False, squeeze_me=True)[varname]
            os.remove(mat_filename)
            return var
        else:
            time.sleep(0.1)
    else:
        raise Exception("Timed out waiting for matlab. Are you at a breakpoint and have the debug server running?")


def compare_vars(varname):
    '''
    Assumes Python is at a breakpoint with a REPL shell. Also assumes MATLAB is at a breakpoint with debug_server
    running. This function returns values for both variable names so they can be compared.
    '''
    import inspect
    frame = inspect.currentframe()
    try:
        return {'python': frame.f_back.f_locals[varname],
                'matlab': request_var(varname)}
    finally:
        del frame


def compare_lists(expected, actual, name=''):
    if all(expected == actual):
        print('{0}: matches at {1}'.format(name, expected))
    else:
        print('{0}: {1} != {2}'.format(name, expected, actual))


def compare_scalars(expected, actual, name=''):
    if expected == actual:
        print('{0}: matches at {1}'.format(name, expected))
    else:
        absolute_difference = abs(expected - actual)
        expected_absolute_value = abs(expected)
        relative_difference = absolute_difference / expected_absolute_value if expected_absolute_value > 0 else absolute_difference
        print('{0}: {1} != {2} absolute difference={3} relative_difference={4}'.format(
            name, expected, actual, absolute_difference, relative_difference))


def compare_ndarray_closeness(expected, actual, name=''):
    def show(message, prefix=''):
        print("{2}{0}: {1}".format(name, message, prefix))

    if not hasattr(expected, 'shape'):
        show("expected has no shape")
        return
    expected_shape = expected.shape
    if not hasattr(actual, 'shape'):
        show('expected.shape={0} but actual has no shape'.format(expected_shape))
        return
    actual_shape = actual.shape
    if actual.shape != expected.shape:
        show('shapes differ ref={0} actual={1}'.format(expected_shape, actual_shape))
        if (len(actual_shape) == len(expected_shape)):
            expected, actual = match_on_overlap(expected, actual, expected_shape, actual_shape)
        else:
            return
    else:
        show('shapes identical {0}'.format(expected_shape))
    try:
        total_size = np.prod(expected_shape)
        expected_max = np.max(expected)
        expected_min = np.min(expected)
        actual_max = np.max(actual)
        actual_min = np.min(actual)
        difference = np.double(actual) - np.double(expected)
        max_difference = np.max(difference)
        max_location = np.argmax(difference)
        min_difference = np.min(difference)
        min_location = np.argmin(difference)
        total_difference = np.sum(np.abs(difference))
        if total_difference == 0 and expected_shape == actual_shape:
            show('are identical', '    ')
            return
        average_difference = total_difference / total_size
        expected_total = np.sum(np.abs(expected))
        relative_difference = total_difference / expected_total if expected_total != 0 else total_difference
        if relative_difference < 0.001:
            show('small differences...', '    ')
        elif relative_difference > 0.1:
            show('LARGE DIFFERENCES...', '    ')
        else:
            show('differences...', '    ')

        show('expected_max={0} expected_min={1} actual_max={2} actual_min={3}'.format(expected_max, expected_min, actual_max, actual_min), '    ')
        show('max_diff={0} min_diff={1}'.format(max_difference, min_difference), '    ')
        show('max_location={0} min_location={1}'.format(max_location, min_location), '    ')
        show('average difference={0} relative_difference={1}'.format(average_difference, relative_difference), '    ')
    except Exception as flaw:
        show('flaw = {0}'.format(str(flaw)))
        traceback.print_exc()


def match_on_overlap(expected, actual, expected_shape, actual_shape):
    overlap_slices = [slice(0, min(expected_dim, actual_dim)) for expected_dim, actual_dim in zip(expected_shape, actual_shape)]
    expected = expected[overlap_slices]
    actual = actual[overlap_slices]
    return expected, actual


def fixup_part1(python_dict, matlab_dict):
    python_dict['numberOfGaussiansPerClass'] = [
        int(number_of_gaussians)
        for number_of_gaussians in python_dict['numberOfGaussiansPerClass']]
    python_dict['numberOfBasisFunctions'] = [
        int(number_of_basis_functions)
        for number_of_basis_functions in python_dict['numberOfBasisFunctions']]
    return python_dict

def fixup_part0(python_dict, matlab_dict):
    imageFileNames = python_dict['imageFileNames']
    if imageFileNames[0] == '/':
        imageFileNames = [imageFileNames]
    python_dict['imageFileNames'] = imageFileNames
    return python_dict


def fixup_part2(python_dict, matlab_dict):
    python_dict['nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel'] = \
        matlab_dict['historyWithinEachMultiResolutionLevel'][-1].finalNodePositionsInTemplateSpace - \
        matlab_dict['historyWithinEachMultiResolutionLevel'][-1].initialNodePositionsInTemplateSpace
    return python_dict


def pass_unchanged(python_dict, matlab_dict):
    return python_dict


FIXUPS = {
    'part0': fixup_part0,
    'part1': fixup_part1,
    'part2': fixup_part2,
}


class CheckpointManager:
    def __init__(self, matlab_dump_dir=None, python_dump_dir=None, detailed=False):
        if matlab_dump_dir is None:
            matlab_dump_dir = MATLAB_DUMP_DIR
        if python_dump_dir is None:
            python_dump_dir = matlab_dump_dir
        self.matlab_dump_dir = matlab_dump_dir
        self.python_dump_dir = python_dump_dir
        self.counts = {}
        self.detailed = detailed

    def increment(self, checkpoint_name):
        self.counts.setdefault(checkpoint_name, 0)
        self.counts[checkpoint_name] += 1
        return self.counts[checkpoint_name]

    def file_name_for_checkpoint(self, dump_dir, checkpoint_name, checkpoint_number=None, suffix='.mat'):
        if checkpoint_number is None and checkpoint_name not in self.counts:
            raise Exception('You must either specify the checkpoint number or call '
                            'increment in the same logical location as in the MATLAB code.')
        checkpoint_number = checkpoint_number or self.counts[checkpoint_name]
        return os.path.join(dump_dir, '{}_{}{}'.format(checkpoint_name, checkpoint_number, suffix))

    def load(self, checkpoint_name, checkpoint_number=None):
        matlab_load_path = self.file_name_for_checkpoint(self.matlab_dump_dir, checkpoint_name, checkpoint_number)
        logger.info('loading from %s', matlab_load_path)
        return scipy.io.loadmat(matlab_load_path, struct_as_record=False, squeeze_me=True)


    def substitute(self, checkpoint_name, checkpoint_number=None, python_dict=None):
        if python_dict is None:
            python_dict = self.load_python(checkpoint_name, checkpoint_number)
        matlab_dict = self.load(checkpoint_name, checkpoint_number)
        result_dict = {key: matlab_dict.get(key, value) for key, value in python_dict.items()}
        fixup = FIXUPS.get(checkpoint_name, pass_unchanged)
        result_dict = fixup(result_dict, matlab_dict)
        return result_dict, python_dict, matlab_dict

    def load_python(self, checkpoint_name, checkpoint_number=None):
        python_load_path = self.file_name_for_checkpoint(self.python_dump_dir, checkpoint_name, checkpoint_number)
        logger.info('loading from %s', python_load_path)
        return scipy.io.loadmat(python_load_path, struct_as_record=False, squeeze_me=True)

    def save(self, value_dict, checkpoint_name, checkpoint_number=None):
        mat_path = self.file_name_for_checkpoint(self.python_dump_dir, checkpoint_name, checkpoint_number)
        logger.info('saving at %s', mat_path)
        scipy.io.savemat(mat_path, value_dict, long_field_names=True)

    def increment_and_save(self, value_dict, checkpoint_name):
        self.increment(checkpoint_name)
        self.save(value_dict, checkpoint_name)

    def save_specification(self, specification, checkpoint_name, checkpoint_number=None):
        json_path = self.file_name_for_checkpoint(
            self.python_dump_dir,
            checkpoint_name,
            checkpoint_number,
            suffix='.json'
        )
        logger.info('saving specification at %s', json_path)
        with open(json_path, 'w') as outfile:
            outfile.write(specification.toJSON())

    def load_specification(self, checkpoint_name, checkpoint_number=None):
        json_path = self.file_name_for_checkpoint(
            self.python_dump_dir,
            checkpoint_name,
            checkpoint_number,
            suffix='.json'
        )
        logger.info('loading specification at %s', json_path)
        with open(json_path, 'r') as infile:
            specification = json.load(infile)
        return specification


class Inspector:
    def __init__(self, target_names):
        self.target_names = target_names

    def inspect(self, prefix, expected_dictionary, actual_dictionary):
        for name in self.target_names:
            message = ':'.join((prefix, name))
            expected_value = expected_dictionary.get(name)
            if expected_value is None:
                print('could not find {0}[{1}] in matlab'.format(prefix, name))
                continue
            actual_value = actual_dictionary.get(name)
            if actual_value is None:
                print('could not find {0}[{1}] in python'.format(prefix, name))
                continue
            self.compare(expected_value, actual_value, message)


class InspectionTeam:
    def __init__(self, checkpoint_name, inspectors):
        self.checkpoint_name = checkpoint_name
        self.inspectors = inspectors

    def inspect(self, checkpoint_manager, actual_dictionary=None):
        if checkpoint_manager is None:
            return
        prefix = ":".join((self.checkpoint_name, str(checkpoint_manager.increment(self.checkpoint_name))))
        try:
            expected_dictionary = checkpoint_manager.load(self.checkpoint_name)
        except Exception as flaw:
            print('no matlab checkpoint {0} for {1}'.format(prefix, checkpoint_manager.matlab_dump_dir))
            return False
        if actual_dictionary is None:
            try:
                actual_dictionary = checkpoint_manager.load_python(self.checkpoint_name)
            except Exception as flaw:
                print('no python checkpoint {0} for {1}'.format(prefix, checkpoint_manager.python_dump_dir))
                return False
        for inspector in self.inspectors:
            inspector.inspect(prefix, expected_dictionary, actual_dictionary)
        return True

    def inspect_all(self, checkpoint_manager):
        while (self.inspect(checkpoint_manager)):
            print('')
        print('')
        print('')


class NdArrayInspector(Inspector):
    def compare(self, expected_value, actual_value, message):
        compare_ndarray_closeness(expected_value, actual_value, message)


def compare_ndarray_dice(expected_value, actual_value, name):
    def show(message, prefix=''):
        print("{2}{0}: {1}".format(name, message, prefix))

    if not hasattr(expected_value, 'shape'):
        show("expected has no shape")
        return
    ref_shape = expected_value.shape
    if not hasattr(actual_value, 'shape'):
        show('expected.shape={0} but actual has no shape'.format(ref_shape))
        return
    item_shape = actual_value.shape
    if actual_value.shape != expected_value.shape:
        show('shapes differ expected={0} actual={1}'.format(ref_shape, item_shape))
        if (len(item_shape) == len(ref_shape)):
            expected_value, actual_value = match_on_overlap(expected_value, actual_value, ref_shape, item_shape)
        else:
            return
    else:
        show('shapes identical {0}'.format(ref_shape))
    try:
        total_size = np.prod(ref_shape)
        matches = actual_value == expected_value
        total_matches = np.sum(matches)
        if total_matches == total_size and ref_shape == item_shape:
            show('are identical', '    ')
        matching_rate = total_matches / total_size
        show('matching rate={0} = {1}/{2}'.format(matching_rate, total_matches, total_size), '    ')
        expected_interior = expected_value != 0
        expected_interior_count = np.sum(expected_interior)
        actual_interior_count = np.sum(actual_value != 0)
        interior_match_count = np.sum(np.logical_and(expected_interior, matches))
        dice = 2 * interior_match_count / (actual_interior_count + expected_interior_count)
        show('dice={0} = 2*{1}/({2} + {3})'.format(
            dice, interior_match_count, expected_interior_count, actual_interior_count), '    ')
    except Exception as flaw:
        show('flaw = {0}'.format(str(flaw)))
        traceback.print_exc()


class DiceInspector(Inspector):
    def compare(self, expected_value, actual_value, message):
        compare_ndarray_dice(expected_value, actual_value, message)


class ScalarInspector(Inspector):
    def compare(self, expected_value, actual_value, message):
        compare_scalars(expected_value, actual_value, message)


class ListInspector(Inspector):
    def compare(self, expected_value, actual_value, message):
        compare_lists(expected_value, actual_value, message)


class ScrambleInspector(Inspector):
    def compare(self, expected_value, actual_value, message):
        expected_sorted = sorted_index_array(expected_value)
        actual_sorted = sorted_index_array(actual_value)
        how_many = min(25, expected_value.shape[0])
        width, height = expected_value.shape
        mapping = [[[0, 0] for i in range(width)] for j in range(height)]
        unmapping = [[[0, 0] for i in range(width)] for j in range(height)]
        for index in range(how_many):
            print('    {0} {1} {2}'.format(message, expected_sorted[index], actual_sorted[index]))
        for index in range(how_many):
            print('    {0} {1} {2}'.format(message, expected_sorted[-1 - index], actual_sorted[-1 - index]))
        for index, ref_val in enumerate(expected_sorted):
            actual_sorted_value = actual_sorted[index]
            r_i = ref_val[1]
            r_j = ref_val[2]
            t_i = actual_sorted_value[1]
            t_j = actual_sorted_value[2]
            mapping[r_i][r_j] = [t_i, t_j]
            unmapping[t_i][t_j] = [r_i, r_j]
        for k in range(width):
            t_i, t_j = mapping[k][k]
            r_i, r_j = unmapping[k][k]
            print('    {0} [{1},{1}] => [{2},{3}]={7}  [{4},{5}] => [{1},{1}]{6}<->{7} {8}<->{9}'.format(
                message, k,
                t_i, t_j,
                r_i, r_j,
                expected_value[k][k],
                actual_value[t_i][t_j],
                expected_value[r_i][r_j],
                actual_value[k][k],
            ))
        for k in range(width):
            t_i, t_j = mapping[0][k]
            r_i, r_j = unmapping[0][k]
            print(
                '    {0} [0,{1}] => [{2},{3}]  [{4},{5}] => [0,{1}] {11}:{12} {10} {13:14} {6}<->{7} {8}<->{9}'.format(
                    message, k,
                    t_i, t_j,
                    r_i, r_j,
                    expected_value[0][k],
                    actual_value[t_i][t_j],
                    expected_value[r_i][r_j],
                    actual_value[0][k],
                    basis(k),
                    basis(r_i),
                    basis(r_j),
                    basis(t_i),
                    basis(t_j),
                ))
        for k in range(width):
            t_i, t_j = mapping[k][0]
            r_i, r_j = unmapping[k][0]
            print(
                '    {0} [{1},0] => [{2},{3}]  [{4},{5}] => [{1},0] {11}:{12} {10} {13}:{14} {6}<->{7} {8}<->{9}'.format(
                    message, k,
                    t_i, t_j,
                    r_i, r_j,
                    expected_value[k][0],
                    actual_value[t_i][t_j],
                    expected_value[r_i][r_j],
                    actual_value[k][0],
                    basis(k),
                    basis(r_i),
                    basis(r_j),
                    basis(t_i),
                    basis(t_j),
                ))
        pass


def basis(i):
    a = i % 5
    i = i // 5
    b = i % 5
    c = i // 5
    return str([a, b, c])


def indexed_array(data):
    width, height = data.shape
    return [(data[i, j], i, j, width * j + 1, height * i + j) for i in range(width) for j in range(height)]


def take_first(a):
    return a[0]


def sorted_index_array(data):
    return sorted(indexed_array(data), key=take_first)


def create_part1_inspection_team():
    return InspectionTeam('part1', [
        NdArrayInspector([
            'biasFieldCoefficients',
            'colors',
            'croppingOffset',
            'FreeSurferLabels',
            'imageBuffers',
            'imageSize',
            'imageToWorldTransformMatrix',
            # 'kroneckerProductBasisFunctions',
            'mask',
            # 'names',
            'nonCroppedImageSize',
            'reducingLookupTable',
            # 'savePath',
            'transformMatrix',
            'voxelSpacing',
        ]),
        ScalarInspector([
            'numberOfClasses',
            'numberOfContrasts',
            'numberOfGaussians',
        ]),
        ListInspector([
            'numberOfGaussiansPerClass',
        ]),
        DiceInspector(['mask']),
    ])


def create_part2_inspection_team():
    return InspectionTeam('part2', [
        NdArrayInspector([
            'biasFieldCoefficients',
            'imageBuffers',
            'means',
            'mixtureWeights',
            'transformMatrix',
            'variances',
        ]),
    ])


def create_part3_inspection_team():
    return InspectionTeam('part3', [
        NdArrayInspector([
            'FreeSurferLabels',
            'freeSurferSegmentation',
            'uncroppedFreeSurferSegmentation',
            'volumesInCubicMm',
        ]),
        DiceInspector([
            'freeSurferSegmentation',
            'uncroppedFreeSurferSegmentation',
        ]),
    ])


def create_multiresWarp_inspection_team():
    return InspectionTeam('multiresWarp', [
        NdArrayInspector([
            'desiredNodePositions',
            'desiredNodePositionsInTemplateSpace',
            'initialNodeDeformationInTemplateSpace',
            'nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel',
            'tmp',
        ]),
    ])


def create_reduced_alphas_inspection_team():
    return InspectionTeam('reducedAlphas', [NdArrayInspector(['reducedAlphas'])])


def create_optimizer_inspection_team():
    return InspectionTeam('optimizer', [
        NdArrayInspector(['nodePositionsAfterDeformation']),
        ScalarInspector(['maximalDeformationApplied']),
    ])


def create_optimizer_exit_inspection_team():
    return InspectionTeam('optimizerPerVoxelExit', [
        ScalarInspector([
            'activeVoxelCount',
            'priorCost',
            'currentCost',
            'costChange',
            'perVoxelDecrease',
            'perVoxelDecreaseThreshold',
            'minLogLikelihoodTimesDeformationPrior',
            'intensityModelParameterCost',
        ]),
    ])


def create_optimizer_em_exit_inspection_team():
    return InspectionTeam('optimizerEmExit', [
        ScalarInspector([
            'activeVoxelCount',
            'priorEMCost',
            'currentEMCost',
            'costChangeEM',
            'changeCostEMPerVoxel',
            'changeCostEMPerVoxelThreshold',
            'minLogLikelihood',
            'intensityModelParameterCost',
        ]),
    ])


def create_bias_correction_inspection_team():
    return InspectionTeam('estimateBiasField', [
        NdArrayInspector([
            'biasFieldCoefficients',
            'lhs',
            'rhs',
            'downSampledBiasField',
            'downSampledBiasCorrectedImageBuffers',
            'biasCorrectedData',
        ]),
        # ScrambleInspector(['lhs']),
        # UnScrambleInspector(['lhs']),
    ])


def distance(ref, item):
    difference = np.double(item) - np.double(ref)
    return np.sum(np.abs(difference))


class UnScrambleInspector(Inspector):
    def compare(self, expected_value, actual_value, message):
        lowest_distance = 150 * 150 * 10000000
        for new_shape in possible_shapes():
            for permutation in permutations(range(6)):
                target_array = np.reshape(actual_value, new_shape)
                permuted_target = np.transpose(target_array, permutation)
                reshaped_target = permuted_target.reshape((150, 150))
                permuted_distance = distance(expected_value, reshaped_target) / (150 * 150)
                if permuted_distance <= lowest_distance:
                    best_permutation = permutation
                    lowest_distance = permuted_distance
                    best_shape = new_shape
                if permuted_distance < 0.01:
                    print('perm={1} shape={0} has distance {2}'.format(
                        str(new_shape), str(permutation), permuted_distance))
        print('*** perm={1} shape={0} has distance {2}'.format(
            str(best_shape), str(best_permutation), lowest_distance))


def possible_shapes():
    for recipe in permutations(range(6), 2):
        i, j = recipe
        if i < j:
            permutation = [5 for k in range(6)]
            permutation[5 - i] = 6
            permutation[5 - j] = 6
            yield permutation


def run_test_cases(case_names=None, testing_dir=None, action=None):
    testing_dir = find_testing_dir(testing_dir)
    if action is None:
        raise ValueError("Nothing to do!!")
    if not case_names:
        case_names = sys.argv[1:]
    if not case_names:
        case_names = ['004']
    for case_file_folder, savePath in valid_case_folders_and_save_paths(case_names, testing_dir):
        action(case_file_folder, savePath)


def make_checkpoint_dir(case_file_folder, subdir_name):
    checkpoint_dir = os.path.join(case_file_folder, subdir_name)
    os.makedirs(checkpoint_dir, exist_ok=True)
    return checkpoint_dir


def find_testing_dir(testing_dir=None):
    if testing_dir is None:
        testing_dir = os.getenv('TESTING_DIR')
    if not testing_dir or not os.path.isdir(testing_dir):
        raise ValueError("testing_dir={0}: must be a valid directory".format(testing_dir))
    return testing_dir


def valid_case_folders_and_save_paths(case_names=None, testing_dir=None):
    testing_dir = find_testing_dir(testing_dir)
    for case_name in valid_case_names(case_names, testing_dir):
        case_file_folder, savePath = \
            find_case_file_folder_and_save_path(case_name, testing_dir)
        if os.path.isdir(case_file_folder):
            yield (case_file_folder, savePath)
        else:
            print("{0} is not a folder".format(case_file_folder))


def valid_case_names(case_names=None, testing_dir=None):
    if not case_names:
        testing_dir = find_testing_dir(testing_dir)
        case_names = sorted([f for f in os.listdir(testing_dir)])
    for case_name in case_names:
        if not 'temp_data' in case_name:
            yield case_name


def find_case_file_folder_and_save_path(case_name, testing_dir=None):
    testing_dir = find_testing_dir(testing_dir)
    case_file_folder = os.path.join(testing_dir, case_name)
    savePath = os.path.join(testing_dir, 'python_temp_data', case_name)
    return case_file_folder, savePath


def create_checkpoint_manager(case_file_folder):
    matlab_dump_dir = make_checkpoint_dir(case_file_folder, 'matlab_checkpoints')
    python_dump_dir = make_checkpoint_dir(case_file_folder, 'python_checkpoints')
    return CheckpointManager(matlab_dump_dir=matlab_dump_dir, python_dump_dir=python_dump_dir)


if __name__ == '__main__':
    x = [5, 5, 6]
    y = [1, 1]
    z = list(np.kron(x, y))
    print(z)
    z.reverse()
    print(z)
