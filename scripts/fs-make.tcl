#! /usr/bin/tclsh

##
## fs-make.tcl
## surfer script: fs-make    [calc,write fieldsign using patch--optional disp]
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

#### set default surface-related files (here for reference)
#setfile curv   ~/surf/$hemi.curv
#setfile patch  ~/surf/$hemi.patch
#setfile fs     */fs/$hemi.fs
#setfile fm     */fs/$hemi.fm

#### file defaults: can reset in csh script with setenv
set eccendir eccen   ;# override w/calling script setenv
set polardir polar   ;# override w/calling script setenv

#### parm defaults: can reset in csh script with setenv
set floatstem sig               ;# float file stem
set realname 2                  ;# analyse infix
set complexname 3               ;# analyse infix
set smoothsteps 50

#### read non-cap setenv vars (or ext w/correct rgbname) to override defaults
source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl

#### read curvature (or sulc)
puts "tksurfer: [file tail $script]: read curvature"
read_binary_curv

#### deal gracefully with old style fs,fm location
if [file exists [file dirname $fs]] {
  # OK: fs dir exists
} else {
  puts "tksurfer: [file tail $script]: new style fs dir not found  ...making it"
  puts "tksurfer:     [file dirname $fs]"
  exec mkdir [file dirname $fs]
}

#### N.B.: overwrites existing (e.g., first) lh.fs, rh.fs

#### fieldsign functions compute junk without patch
if ![file exists $patch] {
  puts "tksurfer: [file tail $script]: ### can't compute fieldsign without patch"
  puts "tksurfer: [file tail $script]: ######## compute fieldsign failed ########"
  exit
}

#### ECCENTRICITY ####
#### read and smooth complex component MRI Fourier transform of *eccen* data
puts "tksurfer: [file tail $script]: read, smooth complex Fourier comp: eccen"
setfile val */$eccendir/${floatstem}${complexname}-$hemi.w
read_binary_values
smooth_val $smoothsteps
shift_values            ;# shift complex component out of way

#### read and smooth real component MRI Fourier transform of *eccen* data
puts "tksurfer: [file tail $script]: read, smooth real Fourier comp: eccen"
setfile val */$eccendir/${floatstem}${realname}-$hemi.w
read_binary_values
smooth_val $smoothsteps
swap_values             ;# swap both components eccentricity out of way

#### POLAR ANGLE ####
#### read and smooth complex component MRI Fourier transform of *theta* data
puts "tksurfer: [file tail $script]: read, smooth complex Fourier comp: polar"
setfile val */$polardir/${floatstem}${complexname}-$hemi.w
read_binary_values
smooth_val $smoothsteps
shift_values            ;# shift complex component out of way

#### read and smooth real component MRI Fourier transform of *theta* data
puts "tksurfer: [file tail $script]: read, smooth real Fourier comp: polar"
setfile val */$polardir/${floatstem}${realname}-$hemi.w
read_binary_values
smooth_val $smoothsteps
swap_values             ;# swap again (r,th)

#### read 2D patch; calc,write fieldsign and mask
setfile patch  ~/surf/${hemi}.${patchname}
puts "tksurfer: [file tail $script]: read patch"
read_binary_patch
puts "tksurfer: [file tail $script]: x-y to polar"
compute_angles          ;# real/complex => ampl/phase
puts "tksurfer: [file tail $script]: compute fieldsign"
compute_fieldsign       ;# gradients, cross prod, scale by geom mean r,th pow
puts "tksurfer: [file tail $script]: write fieldsign"
write_fieldsign         ;# -1,0,1
puts "tksurfer: [file tail $script]: write fieldsign mask"
write_fsmask            ;# thresh, now depends on eccen *and* polar angle

if ![info exists noexit] { exit }

