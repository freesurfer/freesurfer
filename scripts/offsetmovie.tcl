#! /usr/bin/tclsh

##
## offsetmovie.tcl
## surfer script: offsetmovie     [save offset movie; interactive or batch]
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2007/01/05 00:21:14 $
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
##############################################################################
##############################################################################

#### defaults
if ![info exists rgbname] { set rgbname current }     ;# else use current
set offsetsteps 25

### maybe just print defaults
if [info exists justvars] {
  puts "offsetmovie.tcl: ==> print defaults"
  source $env(FREESURFER_HOME)/lib/tcl/printdef.tcl
  return
}

#### two ways to override defaults
if [winfo viewable .] {    ;# make popup
  tmpcontrols "OFFSET MOVIE" { rgbname offsetsteps }
  if {!$userok} { return }
} else {                   ;# batch scripts; re-read env to override defaults
  source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl
}

puts "offsetmovie.tcl: making movie of angle_offsets"
puts "     files: $rgbname-offmov-$hemi-$ext-??.rgb"
prompt

#### save offsetmovie
set angle_offset 0.0
set i 0
while {$i < $offsetsteps} {
  redraw
  save_rgb_named $rgbname-offmov-$hemi-$ext-[format "%02d" $i].rgb
  set angle_offset [expr $angle_offset + [expr 1.0/$offsetsteps]]
  incr i
}

