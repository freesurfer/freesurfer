import os

import numpy as np
import pytest

from samseg.lta import LTA, filter_until_type, parse_assignment, parse_expected, MRI
from samseg.mri_util import load_mgh_header
from samseg.register_atlas_ported import compute_talairach

TEST_LEAF_NAME = 'samseg.talairach.lta'


class TestIta:
    def setup(self):
        self.test_fixture_file_name = os.path.join(os.path.dirname(__file__), TEST_LEAF_NAME)
        self.fs_home = os.path.abspath(os.path.join(__file__, '../../../../'))
        self.expected_dstfile = os.path.join(self.fs_home, 'subjects/fsaverage/mri/orig.mgz')
        self.test_image_folder = os.path.join(os.path.dirname(self.fs_home), 'innolitics_testing')

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

    def test_lta_defaults(self):
        lta = LTA()
        assert lta.nxforms == 1
        assert lta.srcmri is not None
        assert lta.dstmri is not None
        assert lta.sigma == 0.0
        assert lta.mean == [0, 0, 0]
        assert lta.dims == [1, 4, 4]

    def test_mri_defaults(self):
        mri = MRI()
        assert mri.analyzehdr is None
        assert mri.bhdr is None
        assert not mri.valid

    def test_read(self):
        lta = LTA().read(self.test_fixture_file_name)
        self.check_is_valid(lta)

    def check_is_valid(self, lta):
        assert lta is not None
        assert lta.type == 0
        assert lta.nxforms == 1
        assert list(lta.mean) == [0.0, 0.0, 0.0]
        assert lta.sigma == 0.0
        assert list(lta.dims) == [1, 4, 4]
        assert lta.xform.shape == (4, 4)
        assert lta.xform[1, 2] == 0.154134709291283
        assert lta.srcmri.valid == 1
        assert lta.srcfile == '/autofs/cluster/fsm/users/samseg/subjects/Buckner40/004/mri/in_vol.mgz'
        assert list(lta.srcmri.volsize) == [256, 256, 256]
        assert list(lta.srcmri.volres) == [1.0, 1.0, 1.0]
        assert list(lta.srcmri.xras) == [-1.0, 0.0, 0.0]
        assert list(lta.srcmri.yras) == [0.0, 0.0, -1.0]
        assert list(lta.srcmri.zras) == [0.0, 1.0, 0.0]
        assert list(lta.srcmri.cras) == [0.0, 0.0, 0.0]
        assert lta.dstmri.valid == 1
        assert lta.dstfile == '/autofs/cluster/fsm/users/samseg/freesurfer/dist/subjects/fsaverage/mri/orig.mgz'
        assert list(lta.dstmri.volsize) == [256, 256, 256]
        assert list(lta.dstmri.volres) == [1.0, 1.0, 1.0]
        assert list(lta.dstmri.xras) == [-1.0, 0.0, 0.0]
        assert list(lta.dstmri.yras) == [0.0, 0.0, -1.0]
        assert list(lta.dstmri.zras) == [0.0, 1.0, 0.0]
        assert list(lta.dstmri.cras) == [0.0, 0.0, 0.0]
        assert lta.subject == 'fsaverage'
        assert lta.srcmri.valid
        assert lta.dstmri.valid

    def test_write(self, tmpdir):
        lta = LTA()
        lta.read(self.test_fixture_file_name)
        outfile_name = str(tmpdir.mkdir('sub').join('/result.lta'))
        lta.write(outfile_name)
        result_lta = LTA()
        result_lta.read(outfile_name)
        self.check_is_valid(result_lta)

    def test_calculate(self):
        imageFileName = self.test_image_folder + '/buckner40/004/orig.mgz'
        imageToImageTransformMatrix = np.diag([4.0, -3.0, 2.0, 1.0])
        templateImageToWorldTransformMatrix = np.diag([-1.0, 2.0, 3.0, 1.0])
        lta = compute_talairach(
            imageFileName,
            imageToImageTransformMatrix,
            templateImageToWorldTransformMatrix,
            fshome=self.fs_home)
        assert lta is not None
        assert lta.type == 0
        assert lta.subject == 'fsaverage'
        assert lta.srcmri.vol == []
        assert lta.dstmri.vol == []
        assert lta.srcmri.valid
        assert lta.dstmri.valid
        assert lta.srcfile == imageFileName
        assert lta.dstfile == self.expected_dstfile

    def test_load_mgh_header(self):
        [vol, M, mr_parms, volsz] = load_mgh_header(self.expected_dstfile)
        assert list(volsz) == [256, 256, 256, 1]
        assert vol is None
        assert M is not None
        assert mr_parms is not None

    def test_MRI_read_header(self):
        mri = MRI()
        mri.read_header(self.expected_dstfile)
        assert mri is not None
        assert mri.fspec == self.expected_dstfile
        assert mri.vol is None  # reading only header
        assert mri.pwd == os.getcwd()
        assert mri.flip_angle == 0
        assert mri.tr == 0
        assert mri.te == 0
        assert mri.ti == 0
        assert mri.vox2ras is not None
        assert mri.vox2ras0 is not None
        assert list(mri.volsize) == [256, 256, 256]
        assert mri.height == 256
        assert mri.width == 256
        assert mri.depth == 256
        assert mri.nframes == 1
        assert mri.nvoxels == 256 ** 3
        assert mri.xsize == 1.0
        assert mri.ysize == 1.0
        assert mri.zsize == 1.0
        assert mri.volres == [1, 1, 1]
        assert mri.x_r == -1.0
        assert mri.x_a == 0.0
        assert mri.x_s == 0.0
        assert mri.y_r == 0.0
        assert mri.y_a == 0.0
        assert mri.y_s == -1.0
        assert mri.z_r == 0.0
        assert mri.z_a == 1.0
        assert mri.z_s == 0.0
        assert mri.c_r == 0.0
        assert mri.c_a == 0.0
        assert mri.c_s == 0.0
        assert mri.xras == [-1.0, 0.0, 0.0]
        assert mri.yras == [0.0, 0.0, -1.0]
        assert mri.zras == [0.0, 1.0, 0.0]
        assert mri.cras == [0, 0, 0]
        assert mri.vox2ras0[0] == pytest.approx([-1, 0, 0, 128])
        assert mri.vox2ras0[1] == pytest.approx([0, 0, 1, -128])
        assert mri.vox2ras0[2] == pytest.approx([0, -1, 0, 128])
        assert mri.vox2ras0[3] == pytest.approx([0, 0, 0, 1])
