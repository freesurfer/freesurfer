#! /usr/bin/tclsh
##
## setdefpatchview.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2007/01/05 00:21:24 $
##    $Revision: 1.3 $
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
# tksurfer script: setdefpatchview
#############################################################################

# before readenv (just check for patchname; don't blow away aliased patchname)
if { ![info exists aliasedscript] } { set patchname patch } ;# default
foreach var [array names env] {
  if { $var == "patchname" } { set patchname $env($var) }
  if { $var == "nomidzrot" } { set nomidzrot $env($var) } 
  if { $var == "nomidscale" } { set nomidscale $env($var) } 
}

# source position.tcl from *subjectsdir* if there 
if { [file exists $env(SUBJECTS_DIR)/$subject/scripts/position.tcl] } {
  puts "setdefpatchview: ==> read default patchview from position.tcl"
  source $env(SUBJECTS_DIR)/$subject/scripts/position.tcl
} else {
  puts "setdefpatchview: ### position.tcl for $subject does not exist"
  puts "setdefpatchview: ### setting flatzrot to 90 (default)"
  set flatzrot 90
  puts "setdefpatchview: ### setting flatscale to 1.0 (default)"
  set flatscale 1.0
  puts "setdefpatchview: ### setting flatxtrans to 0.0 (default)"
  set flatxtrans 1.0
  puts "setdefpatchview: ### setting flatytrans to 0.0 (default)"
  set flatytrans 1.0
  return
} 

# set rot/trans/scale vars for this hemi/patch combo
if { [info exists $hemi.$patchname.flatzrot] } {
  set flatzrot [set $hemi.$patchname.flatzrot]
  puts "setdefpatchview: set flatzrot $flatzrot"
} else {
  puts "setdefpatchview: ### $hemi.$patchname.flatzrot not set in position.tcl"
  puts "setdefpatchview: ### setting flatzrot to 0 (default)"
  set flatzrot 0
}

if { [info exists $hemi.$patchname.flatscale] } {
  set flatscale [set $hemi.$patchname.flatscale]
  puts "setdefpatchview: set flatscale $flatscale"
} else {
  puts "setdefpatchview: ### $hemi.$patchname.flatscale not set in position.tcl"
  puts "setdefpatchview: ### setting flatscale to 1.0 (default)"
  set flatscale 1.0
}

if { [info exists $hemi.$patchname.flatxtrans] } {
  set flatxtrans [set $hemi.$patchname.flatxtrans]
  puts "setdefpatchview: set flatxtrans $flatxtrans"
} else {
  puts \
      "setdefpatchview: ### $hemi.$patchname.flatxtrans not set in position.tcl"
  puts "setdefpatchview: ### setting flatxtrans to 0.0 (default)"
  set flatxtrans 0.0
}

if { [info exists $hemi.$patchname.flatytrans] } {
  set flatytrans [set $hemi.$patchname.flatytrans]
  puts "setdefpatchview: set flatytrans $flatytrans"
} else {
  puts \
      "setdefpatchview: ### $hemi.$patchname.flatytrans not set in position.tcl"
  puts "setdefpatchview: ### setting flatytrans to 0.0 (default)"
  set flatytrans 0.0
}

