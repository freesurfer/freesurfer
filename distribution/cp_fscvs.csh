#! /bin/tcsh -ef
set echo=1
cd /space/repo/1/dev
tar -X /autofs/space/minerva_001/users/nicks/CVSROOT/fscvs_exclude \
    -cvzf /space/outgoing/fsdev/nicks/tmp/fscvs.tgz dev CVSROOT
