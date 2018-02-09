#!/usr/bin/env python3
# vim:fileencoding=utf-8:ft=python
#
# Author: R.F. Smith <rsmith@xs4all.nl>
# Last modified: 2015-09-23 21:19:02 +0200
#
# To the extent possible under law, Roland Smith has waived all copyright and
# related or neighboring rights to update-modified-keywords.py. This work is
# published from the Netherlands.
# See http://creativecommons.org/publicdomain/zero/1.0/

"""Remove and check out those files that that contain keywords and have
changed since in the last commit in the current working directory."""

from base64 import b64decode
import mmap
import os
import subprocess
import sys


def main(args):
    """Main program.

    Arguments:
        args: command line arguments
    """
    # Check if git is available.
    checkfor(['git', '--version'])
    # Check if .git exists
    if not os.access('.git', os.F_OK):
        print('No .git directory found!')
        sys.exit(1)
    print('{}: Updating modified files.'.format(args[0]))
    # Get modified files
    files = modifiedfiles()
    if not files:
        print('{}: No modified files.'.format(args[0]))
        sys.exit(0)
    files.sort()
    # Find files that have keywords in them
    kwfn = keywordfiles(files)
    if not kwfn:
        print('{}: No keyword files modified.'.format(args[0]))
        sys.exit(0)
    for fn in kwfn:
        os.remove(fn)
    sargs = ['git', 'checkout', '-f'] + kwfn
    subprocess.call(sargs)


def checkfor(args):
    """Make sure that a program necessary for using this script is
    available.

    Arguments:
    args -- string or list of strings of commands. A single string may
            not contain spaces.
    """
    if isinstance(args, str):
        if ' ' in args:
            raise ValueError('No spaces in single command allowed.')
        args = [args]
    try:
        subprocess.check_call(args, stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print("Required program '{}' not found! exiting.".format(args[0]))
        sys.exit(1)


def modifiedfiles():
    """Find files that have been modified in the last commit.

    Returns:
        A list of filenames.
    """
    fnl = []
    try:
        args = ['git', 'diff-tree', 'HEAD~1', 'HEAD', '--name-only', '-r',
                '--diff-filter=ACMRT']
        fnl = subprocess.check_output(args, stderr=subprocess.DEVNULL)
        fnl = fnl.decode('utf8').splitlines()
        # Deal with unmodified repositories
        if len(fnl) == 1 and fnl[0] is 'clean':
            return []
    except subprocess.CalledProcessError as e:
        if e.returncode == 128:  # new repository
            args = ['git', 'ls-files']
            fnl = subprocess.check_output(args, stderr=subprocess.DEVNULL)
            fnl = fnl.decode('utf8').splitlines()
    # Only return regular files.
    fnl = [i for i in fnl if os.path.isfile(i)]
    return fnl


def keywordfiles(fns):
    """Filter those files that have keywords in them

    Arguments:
        fns: A list of filenames.

    Returns:
        A list for filenames for files that contain keywords.
    """
    # These lines are encoded otherwise they would be mangled if this file
    # is checked in my git repo!
    datekw = b64decode('JERhdGU=')
    revkw = b64decode('JFJldmlzaW9u')
    rv = []
    for fn in fns:
        with open(fn, 'rb') as f:
            try:
                mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                if mm.find(datekw) > -1 or mm.find(revkw) > -1:
                    rv.append(fn)
                mm.close()
            except ValueError:
                pass
    return rv


if __name__ == '__main__':
    main(sys.argv)
