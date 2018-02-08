#!/usr/bin/env python3
# vim:fileencoding=utf-8:ft=python
#
# Author: R.F. Smith <rsmith@xs4all.nl>
# Last modified: 2015-05-03 22:06:55 +0200
#
# To the extent possible under law, Roland Smith has waived all copyright and
# related or neighboring rights to kwclean.py. This work is published from the
# Netherlands. See http://creativecommons.org/publicdomain/zero/1.0/

"""Remove the Date and Revision keyword contents from the standard input."""

import io
import re
import sys

if __name__ == '__main__':
    dre = re.compile(''.join([r'\$', r'Date.*\$']))
    drep = ''.join(['$', 'Date', '$'])
    rre = re.compile(''.join([r'\$', r'Revision.*\$']))
    rrep = ''.join(['$', 'Revision', '$'])
    input_stream = io.TextIOWrapper(sys.stdin.buffer, encoding='utf-8')
    for line in input_stream:
        line = dre.sub(drep, line)
        print(rre.sub(rrep, line), end="")
