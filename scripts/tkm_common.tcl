
# $Id: tkm_common.tcl,v 1.2 2003/05/07 23:00:35 kteich Exp $

proc dputs { isString } {

    if { [catch {DebugPrint $isString}] != 0 } {
	puts $isString
    }
}

proc tkm_SourceTclFile { isLibraryPath isFileName } {

    set lPath [list "" "$isLibraryPath/lib/tcl"]

    foreach sPath $lPath {
	
	set sFileName [ file join $sPath $isFileName ]
	dputs "trying $sFileName"
	set nErr [catch { source $sFileName } sResult]
	if { $nErr == 0 } {
	    dputs "read $sFileName"
	    return;
	} else {
	    dputs "error loading $sFileName: $sResult"
	}
    }

    dputs "Couldn't load $isFileName: Not found in $lPath"
}