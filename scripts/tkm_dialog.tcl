##
## tkm_dialog.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2007/01/11 20:15:15 $
##    $Revision: 1.5 $
##
## Copyright (C) 2002-2007, CorTechs Labs, Inc. (La Jolla, CA) and
## The General Hospital Corporation (Boston, MA). 
## All rights reserved.
##
## Distribution, usage and copying of this software is covered under the
## terms found in the License Agreement file named 'COPYING' found in the
## FreeSurfer source code root directory, and duplicated here:
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
##
## General inquiries: freesurfer@nmr.mgh.harvard.edu
## Bug reports: analysis-bugs@nmr.mgh.harvard.edu
##

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
