
##
## unpack_ima.tcl
##
## Original Author: Tony Harris
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2011/03/02 00:04:36 $
##    $Revision: 1.12 $
##
## Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer Software License Agreement' contained
## in the file 'LICENSE' found in the FreeSurfer distribution, and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##

#-------------------- NOTICE ------------------------------------#
# This program is under revision control. Do not edit it without
# going through the proper checkin/checkout steps!
#----------------------------------------------------------------#

# this script looks at the headers of ima files in a predetermined 
# archive directory. The default archive directory is over-ridden 
# by the environment variable ARCHIVE_DIR. the user selects a session 
# and the path to that session is proveded to the script that copies 
# the relevant files via nfs to the local machine and unpacks them 
# locally into b shorts.

% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
