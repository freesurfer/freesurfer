import os
import sys
import platform
import subprocess as sp


def run(command, silent=False, background=False, executable='/bin/bash', logfile=None):
    '''Runs a shell command and returns the exit code.

    Note:
        A command run in the background will always return `None`.
    Args:
        command: Command to run.
        silent: Send output to devnull. Defaults to False.
        background: Run command as a background process. Defaults to False.
        executable: Shell executable. Defaults to `/bin/bash`.
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
                decoded = line.decode('utf-8')
                if logfile is not None:
                    with open(logfile, 'a') as file:
                        file.write(decoded)
                sys.stdout.write(decoded)
        # wait for process to finish
        process.wait()

    return process.returncode


def collect_output(command, executable='/bin/bash'):
    '''Collects the output of a shell command.

    Args:
        command: Command to run.
        executable: Shell executable. Defaults to `/bin/bash`.
    Returns:
        Tuple containing the command output and its exit code.
    '''
    result = sp.run(command, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True, executable=executable)
    return (result.stdout.decode('utf-8'), result.returncode)


def hostname(short=True):
    '''Gets the system hostname.

    Args:
        short: Provide the short hostname. Defaults to True.
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

    # TODO: switch to this (portable across platforms)
    # return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    procstatus = '/proc/self/status'
    if os.path.exists(procstatus):
        with open(procstatus, 'r') as file:
            for line in file:
                if 'VmPeak' in line:
                    return int(line.split()[1])
    return None
