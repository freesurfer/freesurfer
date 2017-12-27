'''
USAGE: First start up the debug server in MATLAB. Please see debug_server.m for instructions
Then call `request_var` and pass in the name of the variable you want to auto-magically
transfer from MATLAB to Python. This is intended to be called with both Python and MATLAB
at a breakpoint to easily compare the values of variables for equality.
'''
import scipy.io
import os
import time

MATLAB_DUMP_DIR = "/Users/ys/work/freesurfer/matlab_dumps" #os.getenv('MATLAB_DUMP_DIR')


def request_var(varname):
    mat_filename = os.path.join(MATLAB_DUMP_DIR, varname + '.mat')
    if os.path.exists(mat_filename):
        os.remove(mat_filename)
    request_file = os.path.join(MATLAB_DUMP_DIR, varname + '.request')
    open(request_file, 'w')
    for i in range(200):
        if not os.path.exists(request_file):
            var = scipy.io.loadmat(mat_filename, struct_as_record=False, squeeze_me=True)[varname]
            os.remove(mat_filename)
            return var
        else:
            time.sleep(0.1)
    else:
        raise Exception("Timed out waiting for matlab. Are you at a breakpoint and have the debug server running?")


def compare_vars(varname):
    '''
    Assumes Python is at a breakpoint with a REPL shell. Also assumes MATLAB is at a breakpoint with debug_server
    running. This function returns values for both variable names so they can be compared.
    '''
    import inspect
    frame = inspect.currentframe()
    try:
        return {'python': frame.f_back.f_locals[varname],
                'matlab': request_var(varname)}
    finally:
        del frame


class CheckpointManager:
    def __init__(self):
        self.counts = {}

    def increment(self, checkpoint_name):
        self.counts.setdefault(checkpoint_name, 0)
        self.counts[checkpoint_name] += 1

    def load(self, checkpoint_name, checkpoint_number=None):
        if checkpoint_name not in self.counts:
            raise Exception('You must either specify the checkpoint number or call '
                            'increment in the same logical location as in the MATLAB code.')
        checkpoint_number = checkpoint_number or self.counts[checkpoint_name]
        mat_path = os.path.join(MATLAB_DUMP_DIR, '{}_{}.mat'.format(checkpoint_name, checkpoint_number))
        return scipy.io.loadmat(mat_path, struct_as_record=False, squeeze_me=True)
