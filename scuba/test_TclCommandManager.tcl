

set value 0
set err [catch { set value [ReturnSingleInt] } sResult]
if { $err != 0 } {
    puts "Error on ReturnSingleInt: $sResult"
}
if { $value != 5 } {
    puts "Error on ReturnSingleInt:"
    puts "Return value was incorrect: $value"
}

set value 999
set err [catch { set value [ReturnSingleFloat] } sResult]
if { $err != 0 } {
    puts "Error on ReturnSingleFloat: $sResult"
}
if { $value != 5.5 } {
    puts "Error on ReturnSingleFloat:"
    puts "Return value was incorrect: $value"
}

set sValue "nope"
set err [catch { set sValue [ReturnSingleString] } sResult]
if { $err != 0 } {
    puts "Error on ReturnSingleString: $sResult"
}
if { $sValue != "hello" } {
    puts "Error on ReturnSingleString:"
    puts "Return value was incorrect: $sValue"
}

set sValue "nope"
set err [catch { set sValue [ReturnSingleTwoWordString] } sResult]
if { $err != 0 } {
    puts "Error on ReturnSingleTwoWordString: $sResult"
}
if { $sValue != "hello world" } {
    puts "Error on ReturnSingleTwoWordString::"
    puts "Return value was incorrect: $sValue"
}

set sValue "nope"
set err [catch { set sValue [ReturnSingleLongString] } sResult]
if { $err != 0 } {
    puts "Error on ReturnSingleLongString: $sResult"
}
if { $sValue != "to neither love nor reverence will thou be tied" } {
    puts "Error on ReturnSingleString:"
    puts "Return value was incorrect: $sValue"
}

set lValue {}
set err [catch { set lValue [ReturnSingleList] } sResult]
if { $err != 0 } {
    puts "Error on ReturnSingleList: $sResult"
}
if { [lindex $lValue 0] != 5 || 
     [lindex $lValue 1] != 5.5 || 
     [lindex $lValue 2] != "hello" } {
    puts "Error on ReturnSingleList:"
    puts "Return value was incorrect: $lValue"
}

set lValue {}
set err [catch { set lValue [ReturnNestedList] } sResult]
if { $err != 0 } {
    puts "Error on ReturnNestedList:: $sResult"
}
set lNested [lindex $lValue 3]
if { [lindex $lValue 0] != 5 || 
     [lindex $lValue 1] != 5.5 || 
     [lindex $lValue 2] != "hello" ||
     [lindex $lNested 0] != 6 || 
     [lindex $lNested 1] != 6.6 || 
     [lindex $lNested 2] != "hello world" } {
    puts "Error on ReturnNestedList:"
    puts "Return value was incorrect: $lValue"
}
