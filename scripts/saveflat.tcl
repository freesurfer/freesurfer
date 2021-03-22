#! /usr/bin/tclsh

##
## saveflat.tcl
## called by flat scripts: saveflat  [save flat image, construct filename]
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

#### make infix for patch rgb name
puts "tksurfer: [file tail $script]: save rgb"
if {$patchname == "patch"} {
  set flatname flat       ;# default rgb infix
} elseif { [string match *nomid* $patchname] } {
  set flatname nomid      ;# default rgb infix of allflat/no midline
} else {
  set flatname $patchname ;# insert other patch names verbatim
}

#### transforms (order: faceview, position.tcl, setenv)
rotate_brain_z $flatzrot
translate_brain_x $flatxtrans
translate_brain_y $flatytrans
scale_brain $flatscale
redraw

#### be nice
if { [info exists nosave] } {
  puts "tksurfer: [file tail $script]: nosave setenv'd => nothing will be saved"
  dontsave_rgb_named $rgbname-$hemi-$flatname.rgb
} elseif { [info exists phasemovie] } {   ;# still get transforms
  puts "tksurfer: [file tail $script]: phasemovie setenv'd => no single rgbs"
  dontsave_rgb_named $rgbname-$hemi-$flatname.rgb
} else {
  save_rgb_named $rgbname-$hemi-$flatname.rgb
}
