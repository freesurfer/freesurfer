proc Dialog_Create { iwwTop isTitle iArgs } {

    global gDialog

    if [winfo exists $iwwTop] {
	switch -- [wm state $iwwTop] {
	    normal {
		raise $iwwTop 
	    }
	    withdrawn - iconic {
		wm deiconify $iwwTop
		catch {
		    wm geometry $iwwTop $gDialog(geo,$iwwTop)
		}
	    }
	}
	return 0
    } else {
	eval {toplevel $iwwTop} $iArgs
	wm title $iwwTop $isTitle
	return 1
    }
}

proc Dialog_Wait { iwwTop iVarName {iFocus {}} } {

    upvar $iVarName var

    bind $iwwTop <Destroy> {list set $iVarName cancel}

    if {[string length $iFocus] == 0} {
	set focus $iwwTop
    }
l    set $saveFocus [focus -displayof $iwwTop]
    focus $iFocus
    catch {
	tkwait visibilty $iwwTop
    }
    catch {
	grab $iwwTop
    }

    tkwait variable $iVarName
    catch {
	grab release $iwwTop
    }
    focus $saveFocus
}

proc Dialog_Close { iwwTop } {

    global gDialog
    
    catch {
	set gDialog(geo,$iwwTop) [wm geometry $iwwTop]
	wm withdraw $iwwTop
    }
}
