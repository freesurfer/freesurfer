#!/usr/bin/env python3
# vim:fileencoding=utf-8:ft=python
#
# Author: R.F. Smith <rsmith@xs4all.nl>
# Last modified: 2015-09-23 21:17:05 +0200
#
# To the extent possible under law, Roland Smith has waived all copyright and
# related or neighboring rights to update-all-keywords.py. This work is
# published from the Netherlands.
# See http://creativecommons.org/publicdomain/zero/1.0/

"""Remove and check out all files under git's control that contain keywords in
the current working directory."""

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
    # Get all files that are controlled by git.
    files = git_ls_files()
    # Remove those that aren't checked in
    mod = git_not_checkedin()
    if mod:
        files = [f for f in files if f not in mod]
    if not files:
        print('{}: Only uncommitted changes, nothing to do.'.format(args[0]))
        sys.exit(0)
    files.sort()
    # Find files that have keywords in them
    kwfn = keywordfiles(files)
    if kwfn:
        print('{}: Updating all files.'.format(args[0]))
        for fn in kwfn:
            os.remove(fn)
        sargs = ['git', 'checkout', '-f'] + kwfn
        subprocess.call(sargs)
    else:
        print('{}: Nothing to update.'.format(args[0]))


def checkfor(args):
    """Make sure that a program necessary for using this script is
    available.

    Arguments:
        args: String or list of strings of commands. A single string may
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


def git_ls_files():
    """Find ordinary files that are controlled by git.

    Returns:
        A list of files
    """
    args = ['git', 'ls-files']
    flist = subprocess.check_output(args).decode('utf8').splitlines()
    return flist


def git_not_checkedin():
    """Find files that are modified but are not checked in.

    Returns:
        A list of modified files that are not checked in.
    """
    lns = subprocess.check_output(['git', 'status', '-s'])
    lns.decode('utf8').splitlines()
    lns = [l.split()[-1] for l in lns]
    return lns


def keywordfiles(fns):
    """Filter those files that have keywords in them

    Arguments:
        fns: A list of filenames.

    Returns:
        A list for filenames for files that contain keywords.
    """
    # These lines are encoded otherwise they would be mangled if this file
    # is checked in!
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
