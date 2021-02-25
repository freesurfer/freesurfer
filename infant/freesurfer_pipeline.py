import os
import sys
import datetime
import argparse
import platform
import getpass
import freesurfer as fs


class CommandPipeline:
    """
    Framework for making a dependency-checking command-line pipline.
    """

    def __init__(self, name, log=None, keep_going=False):
        
        # pipeline name
        self.name = name

        # set default for now
        self._scanning_for_outdated = False

        # configure log file
        self.log = os.path.abspath(log)
        if self.log is not None:
            os.makedirs(os.path.dirname(self.log), exist_ok=True)
            # add a few empty lines if the log already exists
            if os.path.isfile(self.log):
                with open(self.log, 'a') as file:
                    file.write('\n\n\n\n')

        # cache start time
        self.start_time = datetime.datetime.now()

        # print some useful information
        self.info('%s pipeline' % self.name)
        if keep_going:
            self.info(f'Picking-up processing on {datetime.datetime.now()}')
        else:
            self.info(f'New invocation on {datetime.datetime.now()}')
        cmdline = ' '.join(sys.argv[1:])
        self.info(f'Command options: {cmdline}')
        self.info(f'User: {getpass.getuser()}')
        self.info(f'Host: {platform.node()}')
        self.info(f'System: {platform.system()} {platform.release()}')

        # self.info('Total Virtual Memory: %.2f GB' % (psutil.virtual_memory().total * 1e-9))
        # self.info('Avail Virtual Memory: %.2f GB' % (psutil.virtual_memory().available * 1e-9))

        # make sure to define at the end
        self._scanning_for_outdated = keep_going
        if self._scanning_for_outdated:
            self.info('Scanning for changes or missing files...')

    def total_time(self):
        """
        Returns the total time since start.
        """
        return datetime.datetime.now() - self.start_time

    def total_time_str(self):
        """
        Returns the total time since start in string form.
        """
        return str(self.total_time()).split('.')[0]

    def _print_message(self, tag, message, color=None):
        """
        Prints a message to the console and log with timestamp and tag information.
        """

        # convert to time-stamped message
        dt = self.total_time_str()
        tag = tag.rjust(5)
        message = f'{dt} {tag} | {message}'

        # write to log file
        if self.log:
            with open(self.log, 'a') as file:
                file.write(message + '\n')

        # colorize when writing to console
        if color is not None:
            color_func = getattr(fs.utils.term, color)
            message = color_func(message)
        print(message)

    def info(self, message):
        """
        Print a message.
        """
        if self._scanning_for_outdated:
            return
        self._print_message('INFO', message, 'yellow')

    def print(self, message):
        """
        Print a message. Alias for `info()`.
        """
        self.info(message)

    def done(self):
        """
        Exit the pipeline successfully.
        """
        dt = self.total_time_str()
        if self._scanning_for_outdated:
            self._print_message('EXIT', 'Everything up-to-date', 'green')
        else:
            self._print_message('EXIT', '%s pipeline finished successfully' % self.name, 'green')
        sys.exit(0)

    def fatal(self, message=None, code=1):
        """
        Through an error and exit the pipeline.
        """
        dt = self.total_time_str()
        if message is None:
            self._print_message('ERROR', 'Fatal error', 'red')
        else:
            self._print_message('ERROR', 'Fatal: ' + message, 'red')
        sys.exit(code)

    def _is_outdated(self, inputs, outputs):
        """
        Checks whether the output files are outdated given the input files.
        """

        inputs = inputs if isinstance(inputs, list) else [inputs]
        outputs = outputs if isinstance(outputs, list) else [outputs]

        if not inputs and not outputs:
            return True

        # first check if the output files even exist
        if any([not os.path.exists(f) for f in outputs]):
            return True

        if inputs:
            # now compare timestamps across inputs and outputs
            input_mtime = max([os.path.getmtime(f) for f in inputs])
            output_mtime = min([os.path.getmtime(f) for f in outputs])
            if output_mtime < input_mtime:
                return True

        # outputs are up to date
        return False

    def run(self, commands, inputs=[], outputs=[]):
        """
        Run a command or set of commands. Input and output files are used for dependency-checking.
        """

        # check for outdates file (if searching)
        if self._scanning_for_outdated:
            if self._is_outdated(inputs, outputs):
                self._scanning_for_outdated = False
                self.info('Found outdated files, picking-up processing')
            else:
                return

        # force command as a list
        if not isinstance(commands, (list, tuple)):
            commands = [commands]

        # run commands and check for errors
        for cmd in commands:
            self._print_message('CMD', cmd, 'green')
            retcode = fs.run(cmd, logfile=self.log)
            if retcode != 0:
                self.fatal('Command "%s" failed with exit code %d' % (cmd, retcode), retcode)

    def copy(self, source, target):
        """
        Alias to copy a file.
        """
        self.run(f'cp {source} {target}', inputs=source, outputs=target)

    def mkdir(self, directory):
        """
        Alias to make a directory if it does not exist.
        """
        self.run(f'mkdir -p {directory}', outputs=directory)
