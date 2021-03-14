#! /usr/bin/tclsh

##
## curv-flat.tcl
## surfer script: curv-flat  [display curvature on flat]
##
##
## Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer Software License Agreement' contained
## in the file 'LICENSE' found in the FreeSurfer distribution, and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##
#############################################################################
##############################################################################

#### file defaults: can resent in csh script with setenv
set rgbname curv         ;# curv,sulc -- rgbfile name + sets display type
set patchname patch

#### parm defaults: can reset in csh script with setenv
puts "tksurfer: [file tail $script]: set flags"
set overlayflag 0       ;# overlay data on gray brain
set surfcolor 1         ;# draw the curvature under data
set avgflag 1           ;# make half convex/concave
set cslope 5.0         ;# curv sigmoid steepness
set flatzrot 0
set flatscale 1.0

#### read default patch view if there
source $env(FREESURFER_HOME)/lib/tcl/setdefpatchview.tcl

#### read non-cap setenv vars (or ext w/correct rgbname) to override defaults
source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl

#### read curvature (or sulc)
puts "tksurfer: [file tail $script]: read curvature"
if { $rgbname == "curv" } {
  read_binary_curv
} elseif { $rgbname == "sulc"} {
  read_binary_sulc
} else {
  puts "### tksurfer: [file tail $script]: $rgbname: bad coloring type"
  exit
}

#### read 2D patch; calc,write fieldsign and mask
puts "tksurfer: [file tail $script]: read patch"
setfile patch ~/surf/$hemi.$patchname
read_binary_patch       ;# overwrites initial surface read in

#### initial scale and position patch
puts "tksurfer: [file tail $script]: scale, position brain"
open_window
restore_zero_position   ;# undo initial centering
rotate_brain_x -90

#### save requested rgbs (transforms done here)
puts "tksurfer: [file tail $script]: save rgb's"
source $env(FREESURFER_HOME)/lib/tcl/saveflat.tcl

if ![info exists noexit] { exit }

