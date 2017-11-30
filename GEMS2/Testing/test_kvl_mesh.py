import GEMS2Python
from py._path.local import LocalPath
import pytest

MESH_COLLECTION_TEST_FILE = 'Testing/test.txt'


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
        self.empty_collection = GEMS2Python.KvlMeshCollection()

    def teardown(self):
        self.empty_collection = None
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

    def test_access_k(self):
        self.empty_collection.set_k(123)
        actual_k = self.empty_collection.get_k()
        assert 123 == actual_k

    def test_bad_read(self):
        with pytest.raises(Exception):
            bad_collection = GEMS2Python.KvlMeshCollection()
            bad_collection.read('nada_nada_nada')

    @pytest.mark.slowtest
    def test_write(self, mesh_file_name):
        print('writing mesh collection...')
        self.collection.write(mesh_file_name)
        print('...done writing mesh collection')

    def test_empty_write(self, mesh_file_name):
        with pytest.raises(Exception):
            self.empty_collection.write(mesh_file_name)

    @pytest.mark.slowtest
    def test_get_mesh(self):
        reference_mesh = self.collection.get_reference_mesh()
        first_mesh = self.collection.get_mesh(0)
