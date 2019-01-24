#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

# Strip out system, host, machine, user, platform information whjich will vary
#
# Get only first 3 decimal places for line with "Average face Area.." as results may 
# vary on consecutive runs starting at 4th decimal place (at least on the Mac)
rt.run('./mris_decimate -a 0.5 lh.orig.nofix.predec lh.orig.nofix | grep -vi "*Id:*" | grep -vi sysname | grep -vi hostname | grep -vi machine | grep -vi user | sed \'s;^Average Face Area of output.*;&END;\' | sed \'s;[0-9][0-9][0-9]END;;\' > out.mgh')

rt.diff('out.mgh', 'out.ref.mgh')

rt.cleanup()
