import numpy as np


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
    return [float(value) for value in parse_expected(expected_var, next(gen))]


def parse_expected_int_array(expected_var, gen):
    return [int(value) for value in parse_expected(expected_var, next(gen))]


def parse_int_array(gen):
    line = next(gen).strip()
    return [int(value) for value in line.split(' ') if value != '']


def parse_float_array(gen):
    line = next(gen).strip()
    return [float(value) for value in line.split(' ') if value != '']


def parse_2d_float_array(gen, rows=4):
    return [parse_float_array(gen) for row in range(rows)]


class LTA:
    def __init__(self):
        self.srcmri = LTA_MRI()
        self.dstmri = LTA_MRI()
        self.sigma = 0.0
        self.mean = [0.0, 0.0, 0.0]
        self.dims = [1, 4, 4]

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


    def write(self, filename):
        with open(filename, 'w') as file:
            file.writelines('\n'.join([
                '# transform file {0}'.format(filename),
                '# created by LTA.write',
                'type     = %d' % self.type,
                'nxforms  = %d' % self.nxforms,
                'mean     = {}'.format(' '.join([str(val) for val in self.mean])),
                'sigma     = %d' % self.sigma,
                ' '.join(str(val) for val in self.dims),
                ] + self.transformLines() + [
                ] + self.srcmri.formatted_lines('src', self.srcfile) + [
                ] + self.dstmri.formatted_lines('dst', self.dstfile) + [
            ' '.join(['subject', self.subject or 'unknown']),
            ]))
    def transformLines(self):
        return [' '.join([str(col) for col in row]) for row in self.xform]


class LTA_MRI:
    def read(self, content):
        next(content)  # "xxx volume info"
        self.valid = parse_expected_int('valid', content)
        filename = parse_expected_string('filename', content)
        self.volume = parse_expected_int_array('volume', content)
        self.voxelsize = parse_expected_float_array('voxelsize', content)
        self.xras = parse_expected_float_array('xras', content)
        self.yras = parse_expected_float_array('yras', content)
        self.zras = parse_expected_float_array('zras', content)
        self.cras = parse_expected_float_array('cras', content)
        return filename

    def formatted_lines(self, prefix, filename):
        return [
            '{} volume info'.format(prefix),
            'valid = %d' % self.valid,
            'filename = %s' % filename,
            'volume = ' + ' '.join([str(val) for val in self.volume]),
            'voxelsize = ' + ' '.join([str(val) for val in self.voxelsize]),
            'xras = ' + ' '.join([str(val) for val in self.xras]),
            'yras = ' + ' '.join([str(val) for val in self.yras]),
            'zras = ' + ' '.join([str(val) for val in self.zras]),
            'cras = ' + ' '.join([str(val) for val in self.cras]),
        ]
