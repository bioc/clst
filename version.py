#!/usr/bin/env python

import sys
import commands
import re
import time
import os
import fileinput

R_VERSION = commands.getoutput("R --version").split()[2]
VERSION = re.findall(r'Revision: (\d+)',commands.getoutput('svn info'))[0]
DATE = time.strftime('%Y-%m-%d')

print 'R_VERSION',R_VERSION
print 'VERSION', VERSION
print 'DATE', DATE

write = lambda k,v: sys.stdout.write('%s: %s\n' % (k, v.strip()))

for line in fileinput.input('DESCRIPTION',inplace=True):
    key, val = line.split(':')

    if key == 'Version':
        val = '.'.join(val.split('.')[:2] + [VERSION])
    elif key == 'Date':
        val = DATE
    write(key, val)

for line in fileinput.input('man/clst-package.Rd',inplace=True):
    if line.startswith('Version:'):
        line = re.sub(r'\d+(?=\\cr)', VERSION, line)
    elif line.startswith('Date:'):
        line = 'Date: \\tab %s\\cr\n' % DATE
    elif line.startswith('Built:'):
        line = 'Built: \\tab R %s; ; %s; unix\\cr\n' % \
            (R_VERSION, DATE)

    sys.stdout.write(line)


