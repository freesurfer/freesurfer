import GEMS2Python
from py._path.local import LocalPath
import pytest

@pytest.fixture(scope='session')
def mesh_file_name(tmpdir_factory):
    file = tmpdir_factory.mktemp('data').join('test_mesh.txt')
    return file.strpath

class TestMeshCollection:
    def setup(self):
        self.collection = GEMS2Python.KvlMeshCollection()

    def test_access_k(self):
        self.collection.set_k(123)
        actual_k = self.collection.get_k()
        assert 123 == actual_k

    def test_bad_read(self):
        with pytest.raises(Exception):
            self.collection.read('nada_nada_nada')

    @pytest.mark.slowtest
    def test_good_read_write(self, mesh_file_name):
        self.collection.read('Testing/test.txt')
        self.collection.write(mesh_file_name)

    def test_empty_write(self, mesh_file_name):
        with pytest.raises(Exception):
            self.collection.write(mesh_file_name)

    @pytest.mark.slowtest
    def test_get_mesh(self):
        self.collection.read('Testing/test.txt')
        reference_mesh = self.collection.get_reference_mesh()
        first_mesh = self.collection.get_mesh(0)