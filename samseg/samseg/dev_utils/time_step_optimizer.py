import numpy as np
from gems2python import GEMS2Python

from samseg.dev_utils.rasterize_mesh import show_result, create_mesh
from samseg.process_timer import ProcessTimer


def require_np_array(np_array):
    return np.require(np_array, requirements=['F_CONTIGUOUS', 'ALIGNED'])


IMAGE_FILE_NAME = '/home/willy/work/cm_m/innolitics_testing/buckner40/004/orig.mgz'


def full_optimizer_test():
    mesh = create_mesh()
    results = [measure_optimizer_step_time(mesh, thread_count) for thread_count in range(1, 7)]
    for index, elapsed_time in enumerate(results):
        thread_count = 1 + index
        show_result(thread_count, elapsed_time)


def measure_optimizer_step_time(mesh=None, thread_count=1, repeat_count=3):
    GEMS2Python.setGlobalDefaultNumberOfThreads(thread_count)
    test_parts = make_test_parts(mesh)
    optimizer = test_parts['optimizer']
    process_timer = ProcessTimer()
    for _ in range(repeat_count):
        minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer()
    elapsed_time = process_timer.elapsed_time / repeat_count
    print(str(elapsed_time))
    return elapsed_time



def make_test_parts(mesh=None):
    if mesh is None:
        mesh = create_mesh()
    image = create_image()
    calculator = create_calculator(image)
    optimization_parameters = create_optimization_parameters()
    optimizer = GEMS2Python.KvlOptimizer(
        'L-BFGS',
        mesh,
        calculator,
        optimization_parameters
    )
    return {
        'image': image,
        'calculator': calculator,
        'mesh': mesh,
        'optimizer': optimizer,
    }


def create_calculator(image):
    return GEMS2Python.KvlCostAndGradientCalculator('MutualInformation', [image], 'Affine')


def create_image():
    image_file_name = IMAGE_FILE_NAME
    return GEMS2Python.KvlImage(image_file_name)


def create_optimization_parameters():
    return {
        'Verbose': 1.0,
        'MaximalDeformationStopCriterion': 0.005,
        'LineSearchMaximalDeformationIntervalStopCriterion':  0.005,
        'BFGS-MaximumMemoryLength': 12.0  # Affine registration only has 12 DOF
    }


def profile_optimizer_test(thread_count):
    show_result(thread_count, measure_optimizer_step_time(thread_count=thread_count))


if __name__ == '__main__':
    full_optimizer_test()
    # profile_optimizer_test(1)
