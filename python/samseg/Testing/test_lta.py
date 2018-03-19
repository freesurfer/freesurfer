import os
import numpy as np
import pytest

from samseg.lta import LTA, filter_until_type, parse_assignment, parse_expected

TEST_LEAF_NAME = 'samseg.talairach.lta'


class TestIta:
    def setup(self):
        self.test_fixture_file_name = os.path.join(os.path.dirname(__file__), TEST_LEAF_NAME)

    def test_the_test(self):
        assert "abc" == "abc"

    def test_filter_until_type(self):
        input_lines = ['nada', '', 'nope', ' type with indent', 'type is here', 'more', 'and more']
        filtered_lines = filter_until_type(input_lines)
        assert filtered_lines is not None

    def test_parse_assignment(self):
        assert parse_assignment('') is None
        assert parse_assignment('fruit is banana') is None
        assert parse_assignment(' fruit is apple banana') is None
        assert parse_assignment(' fruit = apple') == ['fruit', ['apple']]
        assert parse_assignment(' fruit = apple ') == ['fruit', ['apple']]
        assert parse_assignment(' fruit = apple banana') == ['fruit', ['apple', 'banana']]
        assert parse_assignment(' fruit =     apple banana') == ['fruit', ['apple', 'banana']]
        assert parse_assignment(' fruit = ') == ['fruit', []]

    def test_parse_expected(self):
        assert parse_expected('fruit', 'fruit = apple') == ['apple']
        with pytest.raises(Exception):
            parse_expected('fruit', 'stuff = apple')

    def test_read(self):
        lta = LTA()
        lta.read(self.test_fixture_file_name)
        self.check_is_valid(lta)

    def check_is_valid(self, lta):
        assert lta is not None
        assert lta.type == 0
        assert lta.nxforms == 1
        assert lta.mean == [0.0, 0.0, 0.0]
        assert lta.sigma == 0.0
        assert lta.dims == [1, 4, 4]
        assert lta.xform.shape == (4, 4)
        assert lta.xform[1, 2] == 0.154134709291283
        assert lta.srcmri.valid == 1
        assert lta.srcfile == '/autofs/cluster/fsm/users/samseg/subjects/Buckner40/004/mri/in_vol.mgz'
        assert lta.srcmri.volume == [256, 256, 256]
        assert lta.srcmri.voxelsize == [1.0, 1.0, 1.0]
        assert lta.srcmri.xras == [-1.0, 0.0, 0.0]
        assert lta.srcmri.yras == [0.0, 0.0, -1.0]
        assert lta.srcmri.zras == [0.0, 1.0, 0.0]
        assert lta.srcmri.cras == [0.0, 0.0, 0.0]
        assert lta.dstmri.valid == 1
        assert lta.dstfile == '/autofs/cluster/fsm/users/samseg/freesurfer/dist/subjects/fsaverage/mri/orig.mgz'
        assert lta.dstmri.volume == [256, 256, 256]
        assert lta.dstmri.voxelsize == [1.0, 1.0, 1.0]
        assert lta.dstmri.xras == [-1.0, 0.0, 0.0]
        assert lta.dstmri.yras == [0.0, 0.0, -1.0]
        assert lta.dstmri.zras == [0.0, 1.0, 0.0]
        assert lta.dstmri.cras == [0.0, 0.0, 0.0]
        assert lta.subject == 'fsaverage'

    def test_write(self, tmpdir):
        lta = LTA()
        lta.read(self.test_fixture_file_name)
        outfile_name = str(tmpdir.mkdir('sub').join('/result.lta'))
        lta.write(outfile_name)
        result_lta = LTA()
        result_lta.read(outfile_name)
        self.check_is_valid(result_lta)

