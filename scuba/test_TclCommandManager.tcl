

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

set lValue {}
set err [catch { set lValue [ReturnNestedList2] } sResult]
if { $err != 0 } {
    puts "Error on ReturnNestedList2:: $sResult"
}
set lNested1 [lindex $lValue 0]
set lNested2 [lindex $lValue 1]
if { [lindex $lNested1 0] != "Label 1" || 
     [lindex $lNested1 1] != "Value 1" || 
     [lindex $lNested2 0] != "Label 2" || 
     [lindex $lNested2 1] != "Value 2" } {
    puts "Error on ReturnNestedList2:"
    puts "Return value was incorrect: $lValue"
}

set err [catch { ReturnError } sResult]
if { $err == 0 } {
    puts "ReturnError did not return an error: $sResult"
}
if { $sResult != "This is an error string." } {
    puts "Error on ReturnError:"
    puts "Result value was incorrect: $sResult"
}
