#! /usr/bin/wish

##
## test_tkUtils.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2011/05/02 21:02:44 $
##    $Revision: 1.12 $
##
## Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer Software License Agreement' contained
## in the file 'LICENSE' found in the FreeSurfer distribution, and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##
## General inquiries: freesurfer@nmr.mgh.harvard.edu
## Bug reports: analysis-bugs@nmr.mgh.harvard.edu
##

source ../scripts/tkUtils.tcl
package require Tix

if { [catch {
    set font [tkuLabelFont]
} sResult] } {
    puts "tkuLabelFont failed"
    puts "$sResult"
}

if { [catch {
    set font [tkuNormalFont]
} sResult] } {
    puts "tkuNormalFont failed"
    puts "$sResult"
}

set entryVar "Callback is puts"
tkuMakeEntry .ew \
    -variable entryVar -command {puts $entryVar} -width 20
pack .ew

set entryVar2 "No callback, notify"
tkuMakeEntry .ew2 \
    -variable entryVar2 -width 20 -notify 1
pack .ew2

tkuMakeCheckboxes .cb \
    -orientation h \
    -checkboxes {
	{ -type text -label "label1" -variable bValue1 -command "puts 1" }
	{ -type text -label "label2" -variable bValue2 -command "puts 2" }
    }
pack .cb

tkuMakeCheckboxes .cb2 \
    -orientation v \
    -checkboxes {
	{ -type text -label "label1" -variable bValue1 -command "puts 1" }
	{ -type text -label "label2" -variable bValue2 -command "puts 2" }
    }
pack .cb2

tkuMakeSliders .sw \
    -sliders {
	{ -label "10-20" -variable slider1 -min 10 -max 20
	    -command { puts $slider1 } }
	{ -label "00-20(.5)" -variable slider2 -min 0 -max 20 -resolution 0.5
	    -command { puts $slider2 } }
	{ -label "-10-10e" -variable slider3 -min -10 -max 10
	    -command { puts $slider3 } -entry 1 }
	{ -label "0-10 nolimit" -variable slider4 -min 0 -max 10 
	    -limitentry 0 
	    -command { puts $slider4 } -entry 1 }
    }
pack .sw

set red1 0
set green1 0
set blue1 0
set red2 0
set green2 0
set blue2 0
tkuMakeColorPickers .cp \
    -pickers {
	{ -label "color1" -command {puts "color1: $red1 $green1 $blue1"}
	    -redVariable red1 -blueVariable blue1 -greenVariable green1 }
	{ -label "color2" -command {puts "color2: $red2 $green2 $blue2"}
	    -redVariable red2 -blueVariable blue2 -greenVariable green2 }
    }
pack .cp

tkuMakeToolbar .tb \
    -allowzero false \
    -radio true \
    -variable toolbar \
    -command {tbWrapper} \
    -buttons {
	{ -type text -name tb1 -label "tb1" }
	{ -type text -name tb2 -label "tb2" -balloon "hi there" }
    }
set toolbar tb2

pack .tb

tkuMakeFileSelector .fsw \
    -text "Choose a file" \
    -variable fileName \
    -command {puts "got $fileName"}
set fileName /tmp/blah

pack .fsw


for {set n 0} {$n < 20} {incr n} {
    lappend low1 $n
}

for {set n 0} {$n < 100} {incr n} {
    lappend low2 $n
}

tkuMakeOptionMenu .ow1 \
    -entries $low1 \
    -label "No submenus" \
    -command puts

tkuMakeOptionMenu .ow2 \
    -entries $low2 \
    -labelwidth 10 \
    -label "Submenus" \
    -command puts

pack .ow1 .ow2 -anchor w

proc tbWrapper { isName iValue } {
    puts "tbWrapper: $isName = $iValue"
}

button .bwFile -text "Test File Dlog" -command "TestFileDlog"
button .bwError -text "Test Error Dlog" -command "TestErrorDlog"
button .bwFormattedError -text "Test Formatted Error Dlog" -command "TestFormattedErrorDlog"

pack .bwFile .bwError .bwFormattedError -side bottom

proc TestFileDlog {} {
    
    tkuDoFileDlog -title "Window title" \
	-type1 file \
	-prompt1 "file1 " \
	-note1 "note1, default dir should be /usr/bin" \
	-defaultdir1 "/usr/bin" \
	-shortcutdirs1 "/shortcut/dir/1" \
	-type2 dir \
	-prompt2 "dir2 " \
	-note2 "note2, default dir should be /usr/local/bin" \
	-defaultvalue2 "default value 2" \
	-defaultdir2 "/usr/local/bin" \
	-shortcutdirs2 "/shortcut/dir/2" \
	-type3 checkbox \
	-prompt3 "cb3 " \
	-defaultvalue3 1 \
	-type4 menu \
	-prompt4 "Choose an item:" \
	-menu4 {{0 "Zero"} {1 "One"} {2 "Two"}} \
	-type5 menu \
	-prompt5 "Choose an item (default should be item 3 \"Two\"):" \
	-menu5 {{0 "Zero"} {1 "One"} {2 "Two"}} \
	-defaultitem5 3 \
	-okCmd "puts %s1; puts %s2; puts %s3; puts %s4; puts %s5"
}
 
proc TestFormattedErrorDlog {} {
    tkuFormattedErrorDlog "Error" \
	"Error has occured" \
	"This is a long description of the error"
}

proc TestErrorDlog {} {
    tkuErrorDlog "hi this is an error message whee"
}


tkuFinish
