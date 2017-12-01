import GEMS2Python
import pytest

MESH_COLLECTION_TEST_FILE = 'Testing/test.txt'
MESH_COLLECTION_TEST_POINT_COUNT = 58307
MESH_COLLECTION_TEST_MESH_COUNT = 20
MESH_COLLECTION_TEST_LABEL_COUNT = 17


@pytest.fixture(scope='session')
def mesh_file_name(tmpdir_factory):
    file = tmpdir_factory.mktemp('data').join('test_mesh.txt')
    return file.strpath


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
            TestMeshCollection._collection = GEMS2Python.KvlMeshCollection()
            print('reading mesh collection...')
            TestMeshCollection._collection.read(MESH_COLLECTION_TEST_FILE)
            print('..done reading mesh collection')
        return TestMeshCollection._collection

    def test_construction(self):
        mesh_collection = GEMS2Python.KvlMeshCollection()
        mesh_size = (3, 5, 7)
        domain_size = (10, 11, 13)
        stiffness = 0.25
        number_of_classes = 6
        number_of_meshes = 3
        mesh_collection.construct(mesh_size, domain_size, stiffness, number_of_classes, number_of_meshes)
        mesh_count = mesh_collection.mesh_count
        assert mesh_count == number_of_meshes

    def test_bad_construction_mesh_size(self):
        mesh_collection = GEMS2Python.KvlMeshCollection()
        mesh_size = (3, 5)  # not 3d
        domain_size = (10, 11, 13)
        stiffness = 0.25
        number_of_classes = 6
        number_of_meshes = 3
        with pytest.raises(Exception):
            mesh_collection.construct(mesh_size, domain_size, stiffness, number_of_classes, number_of_meshes)

    def test_bad_construction_domain_size(self):
        mesh_collection = GEMS2Python.KvlMeshCollection()
        mesh_size = (3, 5, 7)
        domain_size = (10, 11, 13, 14)  # not 3d
        stiffness = 0.25
        number_of_classes = 6
        number_of_meshes = 3
        with pytest.raises(Exception):
            mesh_collection.construct(mesh_size, domain_size, stiffness, number_of_classes, number_of_meshes)

    def test_access_k(self):
        empty_collection = GEMS2Python.KvlMeshCollection()
        empty_collection.k = 123
        actual_k = empty_collection.k
        assert 123 == actual_k

    def test_bad_read(self):
        with pytest.raises(Exception):
            bad_collection = GEMS2Python.KvlMeshCollection()
            bad_collection.read('nada_nada_nada')

    @pytest.mark.slowtest
    def test_get_mesh_count(self):
        reference_mesh = self.collection.get_reference_mesh()
        actual_count = self.collection.mesh_count
        assert actual_count == MESH_COLLECTION_TEST_MESH_COUNT

    @pytest.mark.slowtest
    def test_get_reference_mesh(self):
        reference_mesh = self.collection.get_reference_mesh()
        actual_count = reference_mesh.point_count
        assert actual_count == MESH_COLLECTION_TEST_POINT_COUNT

    @pytest.mark.slowtest
    def test_get_first_mesh(self):
        first_mesh = self.collection.get_mesh(0)
        actual_count = first_mesh.point_count
        assert actual_count == MESH_COLLECTION_TEST_POINT_COUNT

    @pytest.mark.slowtest
    def test_get_reference_mesh_using_magic_flag(self):
        reference_mesh = self.collection.get_mesh(-1)
        actual_count = reference_mesh.point_count
        assert actual_count == MESH_COLLECTION_TEST_POINT_COUNT

    @pytest.mark.slowtest
    def test_get_mesh_with_index_way_too_small(self):
        with pytest.raises(Exception):
            mesh = self.collection.get_mesh(-99)

    @pytest.mark.slowtest
    def test_get_mesh_with_index_too_small_by_one(self):
        with pytest.raises(Exception):
            mesh = self.collection.get_mesh(-2)

    @pytest.mark.slowtest
    def test_get_mesh_with_huge_index(self):
        with pytest.raises(Exception):
            mesh = self.collection.get_mesh(99999)

    @pytest.mark.slowtest
    def test_get_mesh_with_index_too_big_by_one(self):
        with pytest.raises(Exception):
            mesh = self.collection.get_mesh(MESH_COLLECTION_TEST_MESH_COUNT)

    @pytest.mark.slowtest
    def test_get_points(self):
        mesh = self.collection.get_reference_mesh()
        points = mesh.points
        [point_count, point_dimensions] = points.shape
        assert point_dimensions == 3
        assert point_count == MESH_COLLECTION_TEST_POINT_COUNT

    @pytest.mark.slowtest
    def test_get_alphas(self):
        mesh = self.collection.get_reference_mesh()
        alphas = mesh.alphas
        [point_count, label_count] = alphas.shape
        assert label_count == MESH_COLLECTION_TEST_LABEL_COUNT
        assert point_count == MESH_COLLECTION_TEST_POINT_COUNT

    @pytest.mark.slowtest
    def test_write(self, mesh_file_name):
        print('writing mesh collection...')
        self.collection.write(mesh_file_name)
        print('...done writing mesh collection')

    def test_empty_write(self, mesh_file_name):
        empty_collection = GEMS2Python.KvlMeshCollection()
        with pytest.raises(Exception):
            empty_collection.write(mesh_file_name)
