
# $Id: tkm_common.tcl,v 1.3 2003/09/03 18:09:36 kteich Exp $

proc dputs { isString } {

    if { [catch {DebugPrint $isString}] != 0 } {
	puts $isString
    }
}

