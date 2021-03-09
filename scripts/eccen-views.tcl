#! /usr/bin/tclsh

##
## eccen-views.tcl
## surfer script: eccen-views  [disp phase-encoded data on 3D folded/unfolded]
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

#### session dir autoset to $cwd/.. when cwd=scripts
#setsession ~/fmri/DALE0308/08798

#### files defaults: can reset in csh script with setenv
set dir 004/image7                  ;# usually set in calling script setenv
set floatstem sig                   ;# float file stem
set realname 2                      ;# analyse infix
set complexname 3                   ;# analyse infix
set rgbname eccen                   ;# name of rgbfiles

#### parm defaults: can reset in csh script with setenv
puts "tksurfer: [file tail $script]: set flags"
set overlayflag 1       ;# overlay data on gray brain
set surfcolor 1         ;# draw the curvature under data
set avgflag 1           ;# make half convex/concave
set complexvalflag 1    ;# two-component data
set colscale 0          ;# 0=wheel,1=heat,2=BR,3=BGR,4=twocondGR,5=gray
set angle_offset 0.82   ;# phase offset
set angle_cycles 1.0    ;# adjust range
set fthresh 0.3         ;# val/curv sigmoid zero (neg=>0)
set fslope 1.5          ;# contast (was fsquash 2.5)
set fmid   0.8          ;# set linear region
set smoothsteps 20
set offset 0.20    ;# default lighting offset 

#### read non-cap setenv vars (or ext w/correct rgbname) to override defaults
source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl

#### print defaults
if [info exists justvars] {
  source $env(FREESURFER_HOME)/lib/tcl/printdef.tcl
  return
}

#### read curvature (or sulc)
puts "tksurfer: [file tail $script]: read curvature"
read_binary_curv

#### setenv eccendir overrides setenv dir           (TODO: omit--force use $dir)
if [info exists eccendir] { set dir $eccendir }  ;# (breaks 2nd after fs sets)

#### read and smooth complex component MRI Fourier transform of data
puts "tksurfer: [file tail $script]: read and smooth complex Fourier comp"
setfile val */$dir/${floatstem}${complexname}-$hemi.w     ;# eccentricity
read_binary_values
smooth_val $smoothsteps
shift_values     ;# shift complex component out of way

#### read and smooth real component MRI Fourier transform of data
puts "tksurfer: [file tail $script]: read and smooth real Fourier comp"
setfile val */$dir/${floatstem}${realname}-$hemi.w     ;# eccentricity
read_binary_values
smooth_val $smoothsteps

#### scale and position brain
puts "tksurfer: [file tail $script]: scale, position brain"
open_window
make_lateral_view       ;# rotate either hemisphere
do_lighting_model -1 -1 -1 -1 $offset ;# -1 => nochange; diffuse curv (def=0.15)

#### save requested rgbs (transforms done here)
puts "tksurfer: [file tail $script]: save rgb's"
source $env(FREESURFER_HOME)/lib/tcl/saveviews.tcl

#### save phasemovie
if [info exists phasemovie] {
  puts "tksurfer: [file tail $script]: save phasemovie"
  source $env(FREESURFER_HOME)/lib/tcl/phasemovie.tcl
}

if ![info exists noexit] { exit }
