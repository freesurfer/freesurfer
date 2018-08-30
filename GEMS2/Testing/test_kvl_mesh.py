import pytest
import numpy as np
import freesurfer.gems as gems

MESH_COLLECTION_TEST_FILE_NAME = 'Testing/test.txt'
TEST_FILE_POINT_COUNT = 58307
TEST_FILE_MESH_COUNT = 20
TEST_FILE_LABEL_COUNT = 17
CONSTRUCTED_MESH_SIZE = (3, 5, 7)
CONSTRUCTED_DOMAIN_SIZE = (10, 11, 13)
CONSTRUCTED_STIFFNESS = 0.25
CONSTRUCTED_NUMBER_OF_CLASSES = 6
CONSTRUCTED_POINT_COUNT = 105


@pytest.fixture(scope='session')
def writeable_mesh_file_name(tmpdir_factory):
    file = tmpdir_factory.mktemp('data').join('test_mesh.txt')
    return file.strpath


def make_a_test_mesh_collection(number_of_meshes):
    mesh_collection = gems.KvlMeshCollection()
    mesh_collection.construct(
        CONSTRUCTED_MESH_SIZE,
        CONSTRUCTED_DOMAIN_SIZE,
        CONSTRUCTED_STIFFNESS,
        CONSTRUCTED_NUMBER_OF_CLASSES,
        number_of_meshes)
    return mesh_collection


@pytest.fixture
def single_mesh_collection():
    return make_a_test_mesh_collection(1)


@pytest.fixture
def simple_mesh():
    return make_a_test_mesh_collection(1).reference_mesh


@pytest.fixture
def triple_mesh_collection():
    return make_a_test_mesh_collection(3)


class TestMeshCollection:
    _collection = None
    # pytest gathers attributes and we do not want those calls to collection to trigger a read
    _okay_to_read = False

    def setup(self):
        TestMeshCollection._okay_to_read = True

    def teardown(self):
        TestMeshCollection._okay_to_read = False

    @property
    def collection(self):
        # Read is very slow, so many tests operate on the same mesh collection.
        # TODO: separate these tests by reading once and then operating on copies of original
        if TestMeshCollection._collection is None and TestMeshCollection._okay_to_read:
            TestMeshCollection._collection = gems.KvlMeshCollection()
            print('reading mesh collection...')
            TestMeshCollection._collection.read(MESH_COLLECTION_TEST_FILE_NAME)
            print('..done reading mesh collection')
        return TestMeshCollection._collection

    def test_construction(self):
        number_of_meshes = 3
        mesh_collection = make_a_test_mesh_collection(number_of_meshes)
        mesh_count = mesh_collection.mesh_count
        assert mesh_count == number_of_meshes
        point_count = mesh_collection.reference_mesh.point_count
        assert point_count == 105

    def test_bad_construction_mesh_size(self):
        mesh_collection = gems.KvlMeshCollection()
        mesh_size = (3, 5)  # not 3d
        domain_size = (10, 11, 13)
        stiffness = 0.25
        number_of_classes = 6
        number_of_meshes = 3
        with pytest.raises(Exception):
            mesh_collection.construct(mesh_size, domain_size, stiffness, number_of_classes, number_of_meshes)

    def test_bad_construction_domain_size(self):
        mesh_collection = gems.KvlMeshCollection()
        mesh_size = (3, 5, 7)
        domain_size = (10, 11, 13, 14)  # not 3d
        stiffness = 0.25
        number_of_classes = 6
        number_of_meshes = 3
        with pytest.raises(Exception):
            mesh_collection.construct(mesh_size, domain_size, stiffness, number_of_classes, number_of_meshes)

    def test_access_k(self, single_mesh_collection):
        single_mesh_collection.k = 123
        actual_k = single_mesh_collection.k
        assert 123 == actual_k

    def test_get_mesh_count(self, triple_mesh_collection):
        actual_count = triple_mesh_collection.mesh_count
        assert actual_count == 3

    def test_reference_mesh(self, triple_mesh_collection):
        reference_mesh = triple_mesh_collection.reference_mesh
        actual_count = reference_mesh.point_count
        assert actual_count == CONSTRUCTED_POINT_COUNT

    def test_get_first_mesh(self, triple_mesh_collection):
        first_mesh = triple_mesh_collection.get_mesh(0)
        actual_count = first_mesh.point_count
        assert actual_count == CONSTRUCTED_POINT_COUNT

    def test_get_reference_mesh_using_magic_flag(self, triple_mesh_collection):
        reference_mesh = triple_mesh_collection.get_mesh(-1)
        actual_count = reference_mesh.point_count
        assert actual_count == CONSTRUCTED_POINT_COUNT

    def test_get_mesh_with_index_way_too_small(self, triple_mesh_collection):
        with pytest.raises(Exception):
            mesh = triple_mesh_collection.get_mesh(-99)

    def test_get_mesh_with_index_too_small_by_one(self, triple_mesh_collection):
        with pytest.raises(Exception):
            mesh = triple_mesh_collection.get_mesh(-2)

    def test_get_mesh_with_huge_index(self, triple_mesh_collection):
        with pytest.raises(Exception):
            mesh = triple_mesh_collection.get_mesh(99999)

    def test_get_mesh_with_index_too_big_by_one(self, triple_mesh_collection):
        with pytest.raises(Exception):
            mesh = triple_mesh_collection.get_mesh(3)

    @pytest.mark.slowtest
    def test_read(self):
        collection = self.collection
        mesh_count = collection.mesh_count
        assert mesh_count == TEST_FILE_MESH_COUNT
        mesh = collection.reference_mesh
        points = mesh.points
        [point_count, point_dimensions] = points.shape
        assert point_dimensions == 3
        assert point_count == TEST_FILE_POINT_COUNT
        alphas = mesh.alphas
        [point_count, label_count] = alphas.shape
        assert label_count == TEST_FILE_LABEL_COUNT
        assert point_count == TEST_FILE_POINT_COUNT
        first_mesh = self.collection.get_mesh(0)
        assert first_mesh.point_count == TEST_FILE_POINT_COUNT

    def test_bad_read(self):
        with pytest.raises(Exception):
            bad_collection = gems.KvlMeshCollection()
            bad_collection.read('nada_nada_nada')

    @pytest.mark.slowtest
    def test_write(self, writeable_mesh_file_name):
        print('writing mesh collection...')
        self.collection.write(writeable_mesh_file_name)
        print('...done writing mesh collection')

    def test_empty_write(self, writeable_mesh_file_name):
        empty_collection = gems.KvlMeshCollection()
        with pytest.raises(Exception):
            empty_collection.write(writeable_mesh_file_name)

    def test_transform_mesh_collection(self, single_mesh_collection):
        initial_transform = np.array([
            [2, 0, 0, 20],
            [0, 3, 0, 30],
            [0, 0, 5, 40],
            [0, 0, 0, 1],
        ], dtype=np.double, order='F')
        kvl_transform = gems.KvlTransform(initial_transform)

        current_points = single_mesh_collection.reference_mesh.points;
        [x, y, z] = current_points[13]
        single_mesh_collection.transform(kvl_transform)

        transformed_points = single_mesh_collection.reference_mesh.points;
        [sx, sy, sz] = transformed_points[13]
        assert 20 + x * 2 == sx
        assert 30 + y * 3 == sy
        assert 40 + z * 5 == sz

class TestMesh:

    def test_get_points(self, simple_mesh):
        points = simple_mesh.points
        [point_count, point_dimensions] = points.shape
        assert point_dimensions == 3
        assert point_count == CONSTRUCTED_POINT_COUNT

    def test_get_alphas(self, simple_mesh):
        alphas = simple_mesh.alphas
        [point_count, label_count] = alphas.shape
        assert label_count == CONSTRUCTED_NUMBER_OF_CLASSES
        assert point_count == CONSTRUCTED_POINT_COUNT

    def test_mesh_set_points(self, simple_mesh):
        points = simple_mesh.points
        [x, y, z] = points[77]
        assert x != 99
        assert y != 100
        assert z != 101
        points[77] = [99, 100, 101]  # Change local points
        refetched_points = simple_mesh.points

        # refetched points are same as before
        [xx, yy, zz] = refetched_points[77]
        assert x == xx
        assert y == yy
        assert z == zz

        # save back, however, will change the points
        simple_mesh.points = points
        again_points = simple_mesh.points
        [xxx, yyy, zzz] = again_points[77]
        assert xxx == 99
        assert yyy == 100
        assert zzz == 101

    def test_set_alphas(self, simple_mesh):
        alphas = simple_mesh.alphas
        [a, b, c, d, e, f] = alphas[33]
        assert a != 99
        assert b != 100
        assert c != 101
        assert d != 102
        assert e != 104
        assert f != 105
        alphas[33] = [99, 100, 101, 102, 103, 104]  # Change local alphas
        refetched_alphas = simple_mesh.alphas

        # refetched alphas are same as before
        [aa, bb, cc, dd, ee, ff] = refetched_alphas[33]
        assert aa == a
        assert bb == b
        assert cc == c
        assert dd == d
        assert ee == e
        assert ff == f

        # save back, however, will change the alphas
        simple_mesh.alphas = alphas
        again_alphas = simple_mesh.alphas
        [aaa, bbb, ccc, ddd, eee, fff] = again_alphas[33]
        assert aaa == 99
        assert bbb == 100
        assert ccc == 101
        assert ddd == 102
        assert eee == 103
        assert fff == 104

    def test_scale_mesh(self, simple_mesh):
        current_points = simple_mesh.points;
        [x, y, z] = current_points[13]
        simple_mesh.scale([2, 3, 5])
        scaled_points = simple_mesh.points;
        [sx, sy, sz] = scaled_points[13]
        assert x * 2 == sx
        assert y * 3 == sy
        assert z * 5 == sz
        # TODO: test cell data content (inverted matrix) has scaled

