#! /usr/bin/tclsh

##
## saveviews.tcl
## called by views scripts: saveviews  [save requested subset of views]
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

#### setenv some subset of: {lat,med,ven,pos,dor,all or custom<transform>}

#### maybe just print defaults
if [info exists justvars] {
  puts "saveviews.tcl: ==> print defaults"
  puts "  pre-set views: setenv {lat,med,ven,pos,dor,all}"
  puts "  customview (one additional view): setenv custom{transform(s)}"
  puts \
 "  order: \[lat\] => xrot => yrot => zrot => xtrans => ytrans => scalepercent"
  source $env(FREESURFER_HOME)/lib/tcl/printdef.tcl
  return
}

#### show lat even when nothing set
if { ![info exists lat] && ![info exists med] && ![info exists ven] && \
     ![info exists dor] && ![info exists pos] && ![info exists all] && \
     [info vars custom*] == "" } {
  set lat 1
}

#### be nice
if { [info exists nosave] } {
  puts "saveviews.tcl: nosave setenv'd => nothing will be saved"
}
if { [info exists phasemovie] } {
  puts "saveviews.tcl: phasemovie setenv'd => no single rgbs"
}

#### save lateral view
if { [info exists lat] || [info exists all] } {
  make_lateral_view
  if { [info exists viewsscale] } { scale_brain $viewsscale }
  redraw
  if { ![info exists nosave] && ![info exists phasemovie] } {
    save_rgb_named $rgbname-$hemi-$ext-lat.rgb
  } else {
    dontsave_rgb_named $rgbname-$hemi-$ext-lat.rgb
  }
}

#### save medial view
if { [info exists med] || [info exists all] } {
  make_lateral_view
  rotate_brain_y 180
  if [info exists viewsscale] { scale_brain $viewsscale }
  redraw
  if { ![info exists nosave] && ![info exists phasemovie] } {
    save_rgb_named $rgbname-$hemi-$ext-med.rgb
  } else {
    dontsave_rgb_named $rgbname-$hemi-$ext-med.rgb
  }
}

#### save ventral view
if { [info exists ven] || [info exists all] } {
  make_lateral_view
  rotate_brain_x 90
  if [info exists viewsscale] { scale_brain $viewsscale }
  redraw
  if { ![info exists nosave] && ![info exists phasemovie] } {
    save_rgb_named $rgbname-$hemi-$ext-ven.rgb
  } else {
    dontsave_rgb_named $rgbname-$hemi-$ext-ven.rgb
  }
}

#### save dorsal view
if { [info exists dor] || [info exists all] } {
  make_lateral_view
  rotate_brain_x -90
  if [info exists viewsscale] { scale_brain $viewsscale }
  redraw
  if { ![info exists nosave] && ![info exists phasemovie] } {
    save_rgb_named $rgbname-$hemi-$ext-dor.rgb
  } else {
    dontsave_rgb_named $rgbname-$hemi-$ext-dor.rgb
  }
}

#### save posterior view
if { [info exists pos] || [info exists all] } {
  make_lateral_view
  if {$hemi == "lh"} {
    rotate_brain_y 90
  } elseif {$hemi == "rh"} {
    rotate_brain_y -90
  }
  redraw
  if [info exists viewsscale] { scale_brain $viewsscale }
  if { ![info exists nosave] && ![info exists phasemovie] } {
    save_rgb_named $rgbname-$hemi-$ext-pos.rgb
  } else {
    dontsave_rgb_named $rgbname-$hemi-$ext-pos.rgb
  }
}

#### save custom view if any custom transforms
if { [info exists customxrot] || \
     [info exists customyrot] || \
     [info exists customzrot] || \
     [info exists customxtrans] || \
     [info exists customytrans] || \
     [info exists customscalepercent] } {
  make_lateral_view
  if { [info exists customxrot] } { rotate_brain_x $customxrot }
  if { [info exists customyrot] } { rotate_brain_y $customyrot }
  if { [info exists customzrot] } { rotate_brain_z $customzrot }
  if { [info exists customxtrans] } { translate_brain_x $customxtrans }
  if { [info exists customytrans] } { translate_brain_y $customytrans }
  if { [info exists customscalepercent] } {
    scale_brain [expr $customscalepercent/100.0]
  }
  redraw
  if { ![info exists nosave] && ![info exists phasemovie] } {
    save_rgb_named $rgbname-$hemi-$ext-cus.rgb
  } else {
    dontsave_rgb_named $rgbname-$hemi-$ext-cus.rgb
  }
}

