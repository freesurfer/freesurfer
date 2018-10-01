import os
import math
import gzip
import numpy as np
import nibabel as nib


def filter_until_type(lines):
    yielding = False
    for line in lines:
        if not yielding:
            split_line = line.split(' ')
            yielding = len(split_line) and split_line[0] == 'type'
        if yielding:
            yield line


def parse_assignment(line):
    parts = line.split('=')
    if len(parts):
        var_name = parts[0].strip()
        if var_name:
            if len(parts) == 2:
                var_value = [item for item in parts[1].strip().split(' ') if item != '']
                return [var_name, var_value]
    return None


def parse_expected(expected_var, line):
    var_name, var_value = parse_assignment(line)
    if var_name == expected_var:
        return var_value
    else:
        message = 'Expected {0} but found {1}'.format(expected_var, var_name)
        raise Exception(message)


def parse_expected_string(expected_var, gen):
    value_list = parse_expected(expected_var, next(gen))
    if len(value_list) == 1:
        return value_list[0]


def parse_expected_string_no_equals(expected_var, gen):
    line_parts = next(gen).strip().split(' ')
    if len(line_parts) == 2:
        var_name = line_parts[0].strip()
        if var_name == expected_var:
            return line_parts[1].strip()
        else:
            message = 'Expected {0} but found {1}'.format(expected_var, var_name)
            raise Exception(message)
    else:
        message = 'Expected {0} but found {1}'.format(expected_var, str(line_parts))
        raise Exception(message)


def parse_expected_int(expected_var, gen):
    value_list = parse_expected(expected_var, next(gen))
    if len(value_list) == 1:
        return int(value_list[0])


def parse_expected_float(expected_var, gen):
    value_list = parse_expected(expected_var, next(gen))
    if len(value_list) == 1:
        return float(value_list[0])


def parse_expected_float_array(expected_var, gen):
    return np.array([float(value) for value in parse_expected(expected_var, next(gen))])


def parse_expected_int_array(expected_var, gen):
    return np.array([int(value) for value in parse_expected(expected_var, next(gen))])


def parse_int_array(gen):
    line = next(gen).strip()
    return np.array([int(value) for value in line.split(' ') if value != ''])


def parse_float_array(gen):
    line = next(gen).strip()
    return np.array([float(value) for value in line.split(' ') if value != ''])


def parse_2d_float_array(gen, rows=4):
    return np.array([parse_float_array(gen) for row in range(rows)])


def nice_array_format(data):
    return ' '.join(['{: 18.15f}'.format(val) for val in data])


class LTA:
    def __init__(self):
        self.srcmri = MRI()
        self.dstmri = MRI()
        self.sigma = 0.0
        self.mean = [0.0, 0.0, 0.0]
        self.dims = [1, 4, 4]
        self.nxforms = 1

    def read(self, filename):
        with open(filename, 'r') as file:
            content = iter(filter_until_type(file.readlines()))
            self.type = parse_expected_int('type', content)
            self.nxforms = parse_expected_int('nxforms', content)
            self.mean = parse_expected_float_array('mean', content)
            self.sigma = parse_expected_float('sigma', content)
            self.dims = parse_int_array(content)
            self.xform = np.array(parse_2d_float_array(content))
            self.srcfile = self.srcmri.read(content)
            self.dstfile = self.dstmri.read(content)
            self.subject = parse_expected_string_no_equals('subject', content)
        return self

    def write(self, filename):
        header_lines = [
            '# transform file {0}'.format(filename),
            '# created by LTA.write',
        ]
        description_lines = [
            'type     = %d' % self.type,
            'nxforms  = %d' % self.nxforms,
            'mean     = {}'.format(' '.join([str(val) for val in self.mean])),
            'sigma    = %d' % self.sigma,
            ' '.join(str(val) for val in self.dims),
        ]
        transform_lines = self.transformLines()
        source_lines = self.srcmri.formatted_lines('src', self.srcfile)
        destination_lines = self.dstmri.formatted_lines('dst', self.dstfile)
        subject_line = [' '.join(['subject', self.subject or 'unknown'])]
        output_lines = header_lines + description_lines + transform_lines + \
                       source_lines + destination_lines + subject_line
        with open(filename, 'w') as file:
            file.writelines('\n'.join(output_lines))

    def transformLines(self):
        return [nice_array_format(row) for row in self.xform]

    def calculate(self):
        pass


def load_mgh_header(filename):
    image = nib.load(filename)
    header = image.header
    volsz = header['dims']
    if header['goodRASFlag']:
        M = header.get_affine()
    else:
        M = None
    mrparms = [header['tr'], header['flip_angle'], header['te'], header['ti']]
    return [None, M, mrparms, np.array(volsz)]


def construct_affine(mat, offset):
    M = np.zeros([4, 4])
    for row in range(3):
        for column in range(3):
            M[row, column] = mat[row, column]
        M[row, 3] = offset[row]
    M[3, 3] = 1.0
    return M


class MRI:
    def __init__(self):
        self.analyzehdr = None
        self.bhdr = None
        self.valid = False

    def read_header(self, filename):
        self.fspec = filename
        self.pwd = os.getcwd()
        [self.vol, M, mr_parms, volsz] = load_mgh_header(filename)
        [self.tr, self.flip_angle, self.te, self.ti] = mr_parms
        volsz[0], volsz[1] = volsz[1], volsz[0]  # first two dimensions get swapped
        self.vox2ras0 = M
        self.volsize = volsz[0:3]

        self.nframes = volsz[3]

        # Everything below is redundant in that they can be derivied from
        # stuff above, but they are in the MRI struct defined in mri.h, so I
        # thought I would add them here for completeness.  Note: MRIwrite()
        # derives all geometry information (ie, direction cosines, voxel
        # resolution, and P0 from vox2ras0. If you change other geometry
        # elements below, it will not be reflected in the output volume.

        M2 = M * M

        self.volres = [
            math.sqrt(sum(M2[:, 0])),  # Column
            math.sqrt(sum(M2[:, 1])),  # Row
            math.sqrt(sum(M2[:, 2])),  # Slice
        ]
        self.xras = [M[row, 0] / self.xsize for row in range(3)]
        self.yras = [M[row, 1] / self.xsize for row in range(3)]
        self.zras = [M[row, 2] / self.xsize for row in range(3)]
        ic = np.array([self.width / 2, self.height / 2, self.depth / 2, 1]).T
        c = M @ ic
        self.cras = list(c[0:3])
        self.valid = True
        self.upate_vox2ras()
        return self

    def read(self, content):
        next(content)  # "xxx volume info"
        self.valid = parse_expected_int('valid', content)
        filename = parse_expected_string('filename', content)
        self.volsize = parse_expected_int_array('volume', content)
        self.volres = parse_expected_float_array('voxelsize', content)
        self.xras = parse_expected_float_array('xras', content)
        self.yras = parse_expected_float_array('yras', content)
        self.zras = parse_expected_float_array('zras', content)
        self.cras = parse_expected_float_array('cras', content)
        self.upate_vox2ras()
        return filename

    def upate_vox2ras(self):
        # Compute the vox2ras matrix
        Mdc = np.array([self.xras, self.yras, self.zras]).T
        Nvox2 = self.volsize / 2
        D = np.diag(self.volres)
        MdcD = Mdc @ D
        P0 = self.cras - MdcD @ Nvox2
        self.vox2ras0 = construct_affine(MdcD, P0)
        pass

    def formatted_lines(self, prefix, filename):
        return [
            '{} volume info'.format(prefix),
            'valid = %d' % 1 if self.valid else 0,
            'filename = %s' % filename,
            'volume = ' + ' '.join([str(val) for val in self.volsize]),
            'voxelsize = ' + nice_array_format(self.volres),
            'xras = ' + nice_array_format(self.xras),
            'yras = ' + nice_array_format(self.yras),
            'zras = ' + nice_array_format(self.zras),
            'cras = ' + nice_array_format(self.cras),
        ]

    @property
    def c_r(self):
        return self.cras[0]

    @property
    def c_a(self):
        return self.cras[1]

    @property
    def c_s(self):
        return self.cras[2]

    @property
    def x_r(self):
        return self.xras[0]

    @property
    def x_a(self):
        return self.xras[1]

    @property
    def x_s(self):
        return self.xras[2]

    @property
    def y_r(self):
        return self.yras[0]

    @property
    def y_a(self):
        return self.yras[1]

    @property
    def y_s(self):
        return self.yras[2]

    @property
    def z_r(self):
        return self.zras[0]

    @property
    def z_a(self):
        return self.zras[1]

    @property
    def z_s(self):
        return self.zras[2]

    @property
    def xsize(self):
        return self.volres[1]

    @property
    def ysize(self):
        return self.volres[0]

    @property
    def zsize(self):
        return self.volres[2]

    @property
    def height(self):
        return self.volsize[0]

    @property
    def width(self):
        return self.volsize[1]

    @property
    def depth(self):
        return self.volsize[2]

    @property
    def vox2ras(self):
        return self.vox2ras0

    @property
    def nvoxels(self):
        return self.height * self.width * self.depth
