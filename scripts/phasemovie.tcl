#! /usr/bin/tclsh

##
## phasemovie.tcl
## surfer script: phasemovie.tcl  [save phase movie; interactive or batch]
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

#### defaults
if ![info exists rgbname] { set rgbname current }
set phasemin -0.5
set phasemax 0.5
set phasewidth 0.03
if {$rgbname == "eccen"} {
  set phasemin -0.3
  set phasemax 0.5
  set phasewidth 0.05
}
if {$rgbname == "polar"} {
  set phasemin -0.5
  set phasemax 0.2
  set phasewidth 0.03
}
if {$rgbname == "downsweep"} {
  set phasemin 0.2
  set phasemax 0.52
  set phasewidth 0.03
}
if {$rgbname == "upsweep"} {
  set phasemin 0.3
  set phasemax 0.7
  set phasewidth 0.03
}
set phasesteps 50
set phasecontour_bright 255

### maybe just print defaults
if [info exists justvars] {
  puts "phasemovie.tcl: ==> print defaults"
  puts "  setenv phasemovie: save default phasemovie using most recent view"
  puts "  (sequence of) single views still shown but not saved"
  source $env(FREESURFER_HOME)/lib/tcl/printdef.tcl
  return
}

#### two ways to override defaults
if [winfo viewable .] {    ;# make popup
  tmpcontrols "PHASE MOVIE" { rgbname phasemin phasemax phasewidth phasesteps }
  if {!$userok} { return }
} else {                   ;# batch scripts; re-read env to override defaults
  source $env(FREESURFER_HOME)/lib/tcl/readenv.tcl
}

if { [info exists nosave] } {
  puts "phasemovie.tcl: nosave is set (csh:unsetenv or tcl:unset to save)"
  return
} else {
  puts "phasemovie.tcl: making phasemovie"
  puts "     files: $rgbname-phmov-$hemi-$ext-??.rgb"
  prompt
}

#### save phasemovie (phase in {-0.5,0.5}; before rev, inv)
set flatext flat
set incr [expr ($phasemax-$phasemin)/$phasesteps]
set phasecontourflag TRUE
set i 0
if {$incr > 0.0} { 
  set phasetest {$phase < $phasemax}
} else {
  set phasetest {$phase > $phasemax}
}
for {set phase $phasemin} $phasetest {set phase [expr $phase+$incr]} {
  # surfer divides phase across max/min boundary--just wrap into [-0.5,0.5]
  set phasecontour_min [expr $phase - $phasewidth/2.0]
  if { $phasecontour_min > 0.5 } {
    set phasecontour_min [expr $phasecontour_min - 1.0]
  }
  if { $phasecontour_min < -0.5 } {
    set phasecontour_min [expr $phasecontour_min + 1.0]
  }
  set phasecontour_max [expr $phase + $phasewidth/2.0]
  if { $phasecontour_max > 0.5 } {
    set phasecontour_max [expr $phasecontour_max - 1.0]
  }
  if { $phasecontour_max < -0.5 } {
    set phasecontour_max [expr $phasecontour_max + 1.0]
  }
  redraw
  if {$flag2d} {
    if { [string match *nomid* $patchname] } { set flatext nomid }
    save_rgb_named $rgbname-phmov-$hemi-$flatext-[format "%02d" $i].rgb
  } else {
    save_rgb_named $rgbname-phmov-$hemi-$ext-[format "%02d" $i].rgb
  }
  incr i
}

set phasecontourflag FALSE

