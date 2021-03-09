#! /usr/bin/tclsh

##
## solmovie-flat.tcl
## tksurfer script: solmovie-flat [render dipole estimates on flat surface]
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
# NOTE: just setenv inpatch for nomid
#setenv flat              ;# savergb flag, like med,ven,.. flags in saveviews
#setenv flatzrot 105      ;# (nomidzrot still recognized)
#setenv flatscale 1.2
#setenv patchname patch-nomid.1000     ;# explicit here (else: patch)

#### file defaults: can reset in csh script with setenv
set patchname patch

#### parm defaults: can reset in csh script with setenv
puts "tksurfer: [file tail $script]: set flags"
set offset 0.20         ;# default lighting offset
set overlayflag 1       ;# overlay data on gray brain
set surfcolor 1         ;# draw the curvature under data
set avgflag 1           ;# make half convex/concave
set colscale 1          ;# 0=wheel,1=heat,2=BR,3=BGR,4=twocondGR,5=gray
set fthresh 1.0
set fslope 0.2
set fmid   3.0
set flatzrot 0
set flatscale 1.0
set smoothsteps 15
set lat0 0
set lat1 400
set dlat 5
set normtype 1

#### read non-cap setenv vars (or ext w/correct rgbname) to override defaults
source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl

set sol_lat0 $lat0
set sol_lat1 $lat1

#### read curvature (or sulc)
puts "tksurfer: [file tail $script]: read curvature"
read_binary_curv

#### read patch (overwrites initial surface read in)
puts "tksurfer: [file tail $script]: read patch"
setfile patch ~/surf/$hemi.$patchname
read_binary_patch

#### read dec file
puts "tksurfer: [file tail $script]: read dec file"
setfile dec $hemi-$dipspacing.dec
read_binary_decimation

#### read dip file
puts "tksurfer: [file tail $script]: read dip file"
setfile dip $hemi.dip
read_binary_dipoles

#### set hemi_num to offset in iop file for hemi
set hemi_num 1         ;# default rh
if {$hemi == "lh"} {
  set hemi_num 2       ;# lh
}

#### read iop file 
puts "tksurfer: [file tail $script]: read iop file"
read_iop ../sol/${iopstem}.iop $hemi_num

#### read rec file 
puts "tksurfer: [file tail $script]: read rec file"
read_rec ../data/${recstem}.rec

#### filter rec data 
puts "tksurfer: [file tail $script]: filter rec data"
filter_recs

#### normalize inverse 
puts "tksurfer: [file tail $script]: normalize inverse"
normalize_inverse

#### read noise covariance matrix
if {$normtype == "3"} {
  read_ncov ../data/${ncovstem}.ncov
}

#### compute timecourses 
puts "tksurfer: [file tail $script]: compute timecourses"
compute_timecourses

#### normalize timecourses
puts "tksurfer: [file tail $script]: normalize timecourses"
normalize_time_courses $normtype

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

#### draw and save collapsed
load_var_from_sol 0
smooth_val_sparse $smoothsteps
if {$sol_plot_type == "2"} {
  smooth_val $smoothsteps
}
redraw
save_rgb_named $rgbstem-$hemi-nomid.rgb

#### draw and save compressed movie
open_rgb_cmp_named $rgbstem-$hemi-nomid.cmp
for {set lat $lat0} {$lat <= $lat1} {set lat [expr $lat + $dlat]} {
  load_vals_from_sol $lat $dlat 0
  smooth_val_sparse $smoothsteps
  redraw
  save_rgb_cmp_frame_named $lat
}

### draw first frame again
load_vals_from_sol $lat0 $dlat 0
smooth_val_sparse $smoothsteps
if {$sol_plot_type == "2"} {
  smooth_val $smoothsteps
}
redraw
save_rgb_cmp_frame_named $lat0

close_rgb_cmp_named

if ![info exists noexit] { exit }
