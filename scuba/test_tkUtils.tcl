#! /usr/bin/wish
#

source $env(DEV)/scripts/tkUtils.tcl
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

tkuMakeToolbar .tb \
    -allowzero false \
    -radio true \
    -variable toolbar \
    -command {tbWrapper} \
    -buttons {
	{ -type text -name tb1 -label "tb1" }
	{ -type text -name tb2 -label "tb2" }
    }
set toolbar tb2

pack .tb

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
	-okCmd "puts %s1; puts %2; puts %3"
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
