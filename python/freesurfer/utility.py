import os
import sys
import re
import platform
import datetime as dt
import subprocess as sp


def run(command, silent=False, background=False, executable='/bin/bash'):
    '''Runs a shell command and returns the exit code.

    Note:
        A command run in the background will always return `None`.
    Args:
        command (str): Command to run.
        silent (bool): Send output to devnull. Defaults to False.
        background (bool): Run command as a background process. Defaults to False.
        executable (str): Shell executable. Defaults to `/bin/bash`.
    Returns:
        Command exit code.
    '''
    # redirect the standard output appropriately
    if silent:
        std = {'stdout': sp.DEVNULL, 'stderr': sp.DEVNULL}
    elif not background:
        std = {'stdout': sp.PIPE, 'stderr': sp.STDOUT}
    else:
        std = {}  # do not redirect
    # run the command
    process = sp.Popen(command, **std, shell=True, executable=executable)
    if not background:
        # write the standard output stream
        if process.stdout:
            for line in process.stdout:
                sys.stdout.write(line.decode('utf-8'))
        # wait for process to finish
        process.wait()
    return process.returncode


def collectOutput(command, executable='/bin/bash'):
    '''Collects the output of a shell command.

    Args:
        command (str): Command to run.
        executable (str): Shell executable. Defaults to `/bin/bash`.
    Returns:
        Tuple containing the command output and its exit code.
    '''
    result = sp.run(command, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True, executable=executable)
    return result.stdout.decode('utf-8'), result.returncode


def source(filename):
    '''Sources a bash file and updates the current environment.

    Args:
        filename: Bash script to source.
    Raises:
        IOError: If file does not exist.
        RuntimeError: If sourcing fails.
    '''
    if not os.path.isfile(filename):
        raise IOError('%s is not a valid file' % filename)
    # source the file and return the resulting environment
    output, retcode = collectOutput('. %s && env' % filename)
    if retcode != 0:
        raise RuntimeError('could not source %s' % filename)
    # parse the env output into a list of pairs
    env = {}
    for line in output.splitlines():
        pair = line.split('=', 1)
        if len(pair) == 1:
            # this line is just a continuation of the previous variable
            env[previous] += '\n' + line
        else:
            key, value = pair
            env[key] = value
            previous = key
    # update the environment
    os.environ.update(env)


def hostname(short=True):
    '''Gets the system hostname.

    Args:
        short (bool): Provide the short hostname. Defaults to True.
    Returns:
        Hostname string.
    '''
    node = platform.node()
    if short:
        return node.split('.')[0]
    return node


def fshome():
    '''Returns the freesurfer home directory.'''
    return os.environ.get('FREESURFER_HOME')


def vmpeak():
    '''Returns the peak memory usage of the process in kilobytes. This only works
    on linux machines because it requires `/proc/self/status`.'''
    procstatus = '/proc/self/status'
    if os.path.exists(procstatus):
        with open(procstatus, 'r') as file:
            for line in file:
                if 'VmPeak' in line:
                    return int(line.split()[1])
    return None


def printPeakMemory(prefix=''):
    '''Prints the peak memory usage of the process if available.'''
    peak = vmpeak()
    if peak is not None:
        if prefix:
            prefix += ' '
        print('%sVmPeak: %d kB' % (prefix, peak))


class Timer:
    '''A simple timer class to track process speed.'''
    def __init__(self, message=None):
        if message: print(message)
        self.start_time = dt.datetime.now()

    @property
    def elapsed(self):
        return dt.datetime.now() - self.start_time

    def mark(self, message):
        print('%s: %s' % (message, str(self.elapsed)))
