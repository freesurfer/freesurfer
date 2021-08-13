import sys
import re
import numpy as np
import datetime as dt

from . import term


def log(message):
    '''Prints a timestamped message.'''
    print(term.cyan(dt.datetime.now().strftime('%H:%M:%S') + ' | ') + message)


def warning(message):
    '''Prints a warning message.'''
    print(term.yellow('warning: ') + message)


def error(message):
    '''Prints an error message.'''
    print(term.red('error: ') + message)


def fatal(message, retcode=1):
    '''Prints an error message and exits (or raises an exception if in interactive mode).

    Args:
        message (str): Message to print.
        retcode (int): Exit code. Defaults to 1.
    '''
    import __main__ as main
    if hasattr(main, '__file__'):
        error(message)
        sys.exit(retcode)
    else:
        raise Exception(message)


def assertion(condition, message, retcode=1):
    '''Exits with error if condition is not met.'''
    if not condition:
        fatal(message, retcode)


class Logger:
    '''A base logging class that writes to a file and stdout (if tee'ing is turned on).

    Note:
        To automatically direct all standard output to a `Logger`,
        use the `RedirectingLogger` subclass.

    Args:
        filename (str): Path to log file.
        tee (bool): Write to stdout as well. Defaults to True.
        mode (str): File opening mode. Overwrites by default.
    '''
    def __init__(self, filename, tee=True, mode='w'):
        self.tee = tee
        self.filename = filename
        self.file = open(self.filename, mode)
        self.stdout = sys.stdout

    def __del__(self):
        self.close()  # be sure to close the file upon deletion

    def write(self, string):
        '''Writes a string to the log file and, if tee'd, to the system stdout.'''
        if term.istty:
            if self.tee:
                # write to stdout
                self.stdout.write(string)
            # write to file
            string = term.removecolor(string)
        self.file.write(string)

    def flush(self):
        '''Flushes the file and, if tee'd, stdout.'''
        if self.tee:
            self.stdout.flush()
        self.file.flush()

    def close(self):
        '''Closes the log file.'''
        self.file.close()


class RedirectingLogger(Logger):
    '''
    A `Logger` that will redirect stdout and stderr to itself after
    initialization. Make sure to call `close` when you would no longer like to
    redirect the output.
    '''
    def __init__(self, *args, **kwargs):
        super(RedirectingLogger, self).__init__(*args, **kwargs)
        sys.stdout = sys.stderr = self

    def close(self):
        '''Closes the log file and resets stdout and stderr to the system defaults.'''
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        super(RedirectingLogger, self).close()


class Table:
    '''Facilitates simple table formatting.

    Args:
        columns: List of column header names.
    '''
    def __init__(self, columns=[], precision=2):
        self.rows = []
        self.columns = columns
        self.precision = precision

    def add_row(self, cells):
        self.addRow(cells)

    def addRow(self, cells):
        '''Adds a row of cells to the table. The number of cells provided must match the
        table's number of columns.'''
        if len(cells) != len(self.columns):
            raise ValueError('number of input cells does not match number of columns')
        self.rows.append([self._format(cell) for cell in cells])

    def get_column(self, idx):
        return [row[idx] for row in self.rows]

    def getColumn(self, idx):
        '''Returns a list of cells for a given column index.'''
        return [self.columns[idx]] + [row[idx] for row in self.rows]

    def _format(self, element):
        if isinstance(element, np.floating):
            return '%.*f' % (self.precision, element)
        return str(element)

    def _pad(self, string, width, rjust=False):
        # width should never be smaller than the length of the string
        string = term.removecolor(string)
        return string.rjust(width) if rjust else string.ljust(width)

    def __str__(self):
        # compute column widths
        widths = []
        for col in range(len(self.columns)):
            colorless = [term.removecolor(cell) for cell in self.getColumn(col)]
            widths.append(len(max(colorless, key=len)))
        strings = []
        # format header
        header = '  '.join([self._pad(name, widths[n]) for n, name in enumerate(self.columns)])
        strings.append(term.bold(header))
        # format rows
        for row in self.rows:
            strings.append('  '.join([self._pad(cell, widths[n]) for n, cell in enumerate(row)]))
        return '\n'.join(strings)

    def write(self, filename):
        with open(filename, 'w') as file:
            if isinstance(self.header, str):
                self.header = [self.header]
            for line in self.header:
                file.write('# %s\n' % line)
            file.write('# NRows %d\n' % len(self.rows))
            file.write('# NTableCols %d\n' % (len(self.columns) + 1))
            file.write('# ColHeaders Index %s\n' % ' '.join(self.columns))

            widths = []
            rjusts = []
            for col in range(len(self.columns)):
                column = self.get_column(col)
                widths.append(len(max(column, key=len)))
                rjusts.append(all([s.replace('-','').replace('.','').isnumeric() for s in column]))

            idxwidth = len(str(len(self.rows)))
            for idx, row in enumerate(self.rows):
                idx = str(idx + 1).rjust(idxwidth)
                row_string = '  '.join([self._pad(cell, widths[n], rjusts[n]) for n, cell in enumerate(row)])
                file.write('%s  %s\n' % (idx, row_string))
