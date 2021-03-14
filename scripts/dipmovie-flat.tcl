#! /usr/bin/tclsh

##
## dipmovie-flat.tcl
## tksurfer script: dipmovie-flat [render dipole estimates on flat surface]
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
#############################################################################
#
#### file defaults: can reset in csh script with setenv
set floatstem sig                   ;# float file stem

#### parm defaults: can reset in csh script with setenv
puts "tksurfer: [file tail $script]: set flags"
set overlayflag 1       ;# overlay data on gray brain
set surfcolor 1         ;# draw the curvature under data
set avgflag 1           ;# make half convex/concave
set colscale 0          ;# 0=wheel,1=heat,2=BR,3=BGR,4=twocondGR,5=gray
set fthresh 0.0         ;# val/curv sigmoid zero (neg=>0)
set fslope 3.0          ;# contast (was fsquash 2.5)
set fmid   0.3          ;# set linear region
set flatzrot 0
set flatscale 1.0
set smoothsteps 15
set pridir prior
set cond_num 0
set lat0 0
set lat1 400
set dlat 50
set zeropad 10
set rgbname dipmovie
set solstem undefined
set prifile undefined
set dipspacing undefined
set offset 0.20

#### read non-cap setenv vars (or ext w/correct rgbname) to override defaults
source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl

#### read curvature (or sulc)
puts "tksurfer: [file tail $script]: read curvature"
read_binary_curv

#### read patch (overwrites initial surface read in)
puts "tksurfer: [file tail $script]: read patch"
read_binary_patch

#### read sol file 
puts "tksurfer: [file tail $script]: read sol file"
read_sols $solstem-$hemi-$dipspacing

#### read pri file 
puts "tksurfer: [file tail $script]: read pri file"
setfile val */$pridir/$prifile-$hemi.pri
read_binary_values

#### read dec file
puts "tksurfer: [file tail $script]: read dec file"
setfile dec $hemi-$dipspacing.dec
read_binary_decimation

#### read dip file
puts "tksurfer: [file tail $script]: read dip file"
read_binary_dipoles

#### normalize time courses
puts "tksurfer: [file tail $script]: normalize time courses"
normalize_time_courses

#### scale and position brain
puts "tksurfer: [file tail $script]: scale, position brain"
open_window
restore_zero_position   ;# undo initial centering
translate_brain_x 0.0
translate_brain_y 0.0
translate_brain_z 0.0
rotate_brain_x -90
rotate_brain_z $flatzrot
scale_brain $flatscale
do_lighting_model -1 -1 -1 -1 $offset  ;# -1(def); offset=curvdiffuse (def=0.15)

#### backward compatibility
if { [info exists nomidzrot] && [info exists nomidscale] } {
  restore_zero_position
  rotate_brain_x -90
  rotate_brain_z $nomidzrot
  scale_brain $nomidscale
}

#### draw and save images
open_rgb_cmp_named $rgbname-$hemi-flat.cmp
for {set lat $lat0} {$lat <= $lat1} {set lat [expr $lat + $dlat]} {
  load_vals_from_sol $lat $dlat $cond_num
  smooth_val_sparse $smoothsteps
  redraw
  if {$lat == $lat0} {
    for {set i 0} {$i < $zeropad} {incr i} {save_rgb_cmp_frame_named $lat}
  }
  save_rgb_cmp_frame_named $lat
}
close_rgb_cmp_named
exit
