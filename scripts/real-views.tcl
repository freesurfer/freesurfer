#! /usr/bin/tclsh

##
## real-views.tcl
## tksurfer script: real-views [disp two-condition data on 3D folded/unfolded]
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
##############################################################################
##############################################################################

#### session dir autoset to $cwd/.. when cwd=scripts
#setsession ~/fmri/DALE0308/08798

#### file defaults: can reset in csh script with setenv
set floatstem sig                   ;# float file stem
set statname ""                     ;# TODO
set rgbname real                    ;# name of rgbfiles

#### parm defaults: can reset in csh script with setenv
puts "tksurfer: [file tail $script]: set flags"
set overlayflag 1       ;# overlay data on gray brain
set surfcolor 1         ;# draw the curvature under data
set avgflag 1           ;# make half convex/concave
set complexvalflag 0    ;# two-component data
set colscale 1          ;# 0=wheel,1=heat,2=BR,3=BGR,4=twocondGR,5=gray
set angle_offset 0.0    ;# phase offset
set angle_cycles 1.0    ;# adjust range
set fthresh 0.3         ;# val/curv sigmoid zero (neg=>0)
set fslope 1.5          ;# contast (was fsquash 2.5)
set fmid   1.5          ;# set linear region
set flatzrot 0
set flatscale 1.0
set smoothsteps 0
set offset 0.2

#### read non-cap setenv vars (or ext w/correct rgbname) to override defaults
source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl

#### read curvature (or sulc)
puts "tksurfer: [file tail $script]: read curvature"
read_binary_curv

#### read and smooth real data
puts "tksurfer: [file tail $script]: read and smooth real data"
setfile val */$dir/${floatstem}${statname}-$hemi.w
read_binary_values
smooth_val $smoothsteps

#### scale and position brain
puts "tksurfer: [file tail $script]: scale, position brain"
open_window
make_lateral_view       ;# rotate either hemisphere
do_lighting_model -1 -1 -1 -1 $offset ;# -1(no change),offset=curvdiffuse (0.15)

#### save requested rgbs
puts "tksurfer: [file tail $script]: save rgb's"
source $env(FREESURFER_HOME)/lib/tcl/saveviews.tcl

if ![info exists noexit] { exit }

