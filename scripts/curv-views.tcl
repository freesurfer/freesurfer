#! /usr/bin/tclsh

##
## curv-views.tcl
## surfer script: curv-views  [display curvature on requested views]
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2007/01/05 00:21:02 $
##    $Revision: 1.4 $
##
## Copyright (C) 2002-2007,
## The General Hospital Corporation (Boston, MA). 
## All rights reserved.
##
## Distribution, usage and copying of this software is covered under the
## terms found in the License Agreement file named 'COPYING' found in the
## FreeSurfer source code root directory, and duplicated here:
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
##
## General inquiries: freesurfer@nmr.mgh.harvard.edu
## Bug reports: analysis-bugs@nmr.mgh.harvard.edu
##
#############################################################################
##############################################################################

#### file defaults: can reset in csh script with setenv
set rgbname curv         ;# curv,sulc -- rgbfile name + sets display type

#### parm defaults: can reset in csh script with setenv
puts "tksurfer: [file tail $script]: set flags"
set overlayflag 0       ;# overlay data on gray brain
set surfcolor 1         ;# draw the curvature under data
set avgflag 1           ;# make half convex/concave
set cslope 5.0          ;# curv sigmoid steepness
set offset 0.20         ;# default lighting offset

#### read non-cap setenv vars (or ext w/correct rgbname) to override defaults
source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl

#### read curvature (or sulc)
puts "tksurfer: [file tail $script]: read curvature"
if { [string match *curv* $rgbname] } {
  read_binary_curv
} elseif { [string match *sulc* $rgbname] } {
  read_binary_sulc
} else {
  puts \
  "### tksurfer: [file tail $script]: $rgbname: bad type: curv/sulc not in name"
  exit
}

#### scale and position brain
puts "tksurfer: [file tail $script]: scale, position brain"
open_window
make_lateral_view
do_lighting_model -1 -1 -1 -1 $offset   ;# -1:nochange offset:diffusecurv

#### save requested rgbs
puts "tksurfer: [file tail $script]: save rgb's"
source $env(FREESURFER_HOME)/lib/tcl/saveviews.tcl

if ![info exists noexit] { exit }

