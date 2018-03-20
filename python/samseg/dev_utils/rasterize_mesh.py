import numpy as np
from gems2python import GEMS2Python

from samseg.process_timer import ProcessTimer

MESH_COLLECTION_NAME = '/home/willy/work/cm_p/innolitics_testing/atlas/20Subjects_smoothing2_down2_smoothingForAffine2/atlasForAffineRegistration.txt.gz'


def require_np_array(np_array):
    return np.require(np_array, requirements=['F_CONTIGUOUS', 'ALIGNED'])


def create_mesh(mesh_collection_file_name=MESH_COLLECTION_NAME):
    scaling = 0.9
    scalingMatrix = np.diag([scaling, scaling, scaling, 1.0])
    initialWorldToWorldTransformMatrix = scalingMatrix
    mesh_collection = GEMS2Python.KvlMeshCollection()
    mesh_collection.read(mesh_collection_file_name)
    K = 1e-7
    K = K / (scaling * scaling * scaling)
    mesh_collection.k = K
    mesh_collection.transform(GEMS2Python.KvlTransform(require_np_array(initialWorldToWorldTransformMatrix)))
    mesh = mesh_collection.reference_mesh
    return mesh

def rasterize_a_mesh(mesh_collection_file_name, thread_count=1, repeat_count=10):
    mesh = create_mesh(mesh_collection_file_name)
    GEMS2Python.setGlobalDefaultNumberOfThreads(thread_count)
    imageSize = [200, 200, 200]
    process_timer = ProcessTimer()
    for _ in range(repeat_count):
        priors = mesh.rasterize(imageSize)
    elapsed_time = process_timer.elapsed_time / repeat_count
    print(str(elapsed_time))
    return elapsed_time


def full_mesh_test():
    mesh_collection_file_name = MESH_COLLECTION_NAME
    results = [rasterize_a_mesh(mesh_collection_file_name, thread_count) for thread_count in range(1, 7)]
    for index, elapsed_time in enumerate(results):
        thread_count = 1 + index
        show_result(thread_count, elapsed_time)


def show_result(thread_count, elapsed_time):
    total_seconds = elapsed_time.total_seconds()
    work = total_seconds * thread_count
    productivity = 1.0 / total_seconds
    thread_productivity = productivity / thread_count
    print('threads={0}   time={1}   threads*time={2}   throughput={3}   throughput per thread={4}'.format(thread_count,
                                                                                                          total_seconds,
                                                                                                          work,
                                                                                                          productivity,
                                                                                                          thread_productivity))


def profile_test(thread_count):
    mesh_collection_file_name = MESH_COLLECTION_NAME
    elapsed_time = rasterize_a_mesh(mesh_collection_file_name, thread_count)
    show_result(thread_count, elapsed_time)


if __name__ == '__main__':
    # full_mesh_test()
    profile_test(6)
