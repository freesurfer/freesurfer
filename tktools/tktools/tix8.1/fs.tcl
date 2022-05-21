# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: fs.tcl,v 1.4.2.1 2001/11/03 07:26:10 idiscovery Exp $
#
# tixAssert --
#
#	Debugging routine. Evaluates the test script in the context of the
#	caller. The test script is responsible for generating the error.
#	
proc tixAssert {script} {
    uplevel $script
}

proc tixAssertNorm {path} {
    if {![tixFSIsNorm $path]} {
	error "\"$path\" is not a normalized path"
    }
}

proc tixAssertVPath {vpath} {
    if {![tixFSIsVPath $vpath]} {
	error "\"$vpath\" is not a VPATH"
    }
}

# tixFSAbsPath --
#
#	Converts $path into an normalized absolute path
#
proc tixFSAbsPath {path} {
    return [lindex [tixFSNorm [tixFSVPWD] $path] 0]
}

# tixFSVPWD --
#
#	Returns the VPATH of the current directory.
#
proc tixFSVPWD {} {
    return [tixFSVPath [tixFSPWD]]
}

if {![info exists tcl_platform] || $tcl_platform(platform) == "unix"} {

# tixFSPWD --
#
#	Return the current directory
#
proc tixFSPWD {} {
    return [pwd]
}

# tixFSDisplayName --
#
#	Returns the name of a normalized path which is usually displayed by
#	the OS
#
proc tixFSDisplayName {normpath} {
    tixAssert {
	tixAssertNorm $normpath
    }
    return $normpath
}

proc tixFSIsAbsPath {path} {
    return [tixStrEq [string index $path 0] /]
}

# tixFSIsNorm_os --
#
#	Returns true iff this pathname is normalized, in the OS native name
#	format
#
proc tixFSIsNorm_os {path} {
    return [tixFSIsNorm $path]
}

proc tixFSIsNorm {path} {
    if {[tixStrEq $path /]} {
	return 1
    }

    # relative path
    #
    if {![regexp -- {^/} $path]} {
	return 0
    }

    if {[regexp -- {/[.]$} $path]} {
	return 0
    }
    if {[regexp -- {/[.][.]$} $path]} {
	return 0
    }
    if {[regexp -- {/[.]/} $path]} {
	return 0
    }
    if {[regexp -- {/[.][.]/} $path]} {
	return 0
    }
    if {[tixStrEq $path .]} {
	return 0
    }
    if {[tixStrEq $path ..]} {
	return 0
    }

    # Tilde
    #
    if {[regexp -- {^~} $path]} {
	return 0
    }

    # Double slashes
    #
    if {[regexp -- {//} $path]} {
	return 0
    }

    # Trailing slashes
    #
    if {[regexp -- {/$} $path]} {
	return 0
    }

    return 1
}

# tixFSIsValid --
#
#	Checks whether a native pathname contains invalid characters.
#
proc tixFSIsValid {path} {
    return 1
}

proc tixFSIsVPath {vpath} {
    return [tixFSIsNorm $vpath]
}

# tixFSVPath --
#
#	Converts a native pathname to its VPATH
#
proc tixFSVPath {path} {
    tixAssert {
	tixAssertNorm $path
    }
    return $path
}

# tixFSPath --
#
#	Converts a vpath to a native pathname
proc tixFSPath {vpath} {
    tixAssert {
	tixAssertVPath $vpath
    }
    return $vpath
}

# tixFSTildeSubst -- [Unix only]
#
#	Substitutes any leading tilde characters if possible. No error is
#	generated if the user doesn't exist.
#
proc tixFSTildeSubst {text} {
    if {[tixStrEq [string index $text 0] ~]} {
	# The following will report if the user doesn't exist
	if {[catch {
	    file isdir $text
	}]} {
	    return ./$text
	}
	return [tixFile tilde $text]
    } else {
	return $text
    }
}

# tixFSNorm --
#
#	Interprets the user's input and return file information about this
#	input.
#
# Arguments:
#	See documentation (docs/Files.txt)
#
proc tixFSNorm {context text {defFile ""} {flagsVar ""} {errorMsgVar ""}} {
    tixAssert {
	tixAssertVPath $context
    }

    if {![tixStrEq $errorMsgVar ""]} {
	upvar $errorMsgVar errorMsg
    }
    if {![tixStrEq $flagsVar ""]} {
	upvar $flagsVar flags
    }

    set hasDirSuffix [regexp -- {/$} $text]
    set text [tixFSTildeSubst $text]
    set text [_tixJoin $context $text]

    if {$hasDirSuffix || [file isdir $text]} {
	set dir $text
	set tail $defFile
    } else {
	set dir [file dirname $text]
	set tail [file tail $text]
    }

    set norm $dir/$tail
    regsub -all -- /+ $norm / norm
    if {![tixStrEq $norm /]} {
	regsub -- {/$} $norm "" norm
    }

    if {![info exists flag(noPattern)]} {
	set isPat 0
	foreach char [split $tail ""] {
	    if {$char == "*" || $char == "?"} {
		set isPat 1
		break
	    }
	}
	if {$isPat} {
	    return [list $norm $dir "" $tail]
	}
    }

    return [list $norm $dir $tail ""]
}

# _tixJoin -- [Internal]
# 
#	Joins two native pathnames.
#
proc _tixJoin {p1 p2} {
    if {[tixStrEq [string index $p2 0] /]} {
	return [_tixNormalize $p2]
    } else {
	return [_tixNormalize $p1/$p2]
    }
}

# tixFSNormDir --
#
#	Normalizes an absolute path.
#
proc tixFSNormDir {dir} {
    set dir [tixFile tilde $dir]
    if {![tixStrEq [string index $dir 0] /]} {
	error "\"$dir\" must be an absolute pathname"
    }
    if {![file isdir $dir]} {
	error "\"$dir\" is not a directory"
    }
    return [_tixNormalize $dir]
}

# _tixNormalize --
#
#	Normalizes an absolute pathname.
#
# 	$dir must be an absolute pathname
#
proc _tixNormalize {path} {
    tixAssert {
	if {![tixStrEq [string index $path 0] /]} {
	    error "\"$path\" must be an absolute pathname"
	}
    }

    # Don't be fooled: $path doesn't need to be a directory. The following
    # code just makes it easy to get rid of trailing . and ..
    #
    set path $path/
    regsub -all -- /+ $path / path
    while {1} {
	if {![regsub -- {/\./} $path "/" path]} {break}
    }
    while {1} {
	if {![regsub -- {/\.$} $path "" path]} {break}
    }

    while {1} {
	if {![regsub -- {/[^/]+/\.\./} $path "/" path]} {break}
	while {1} {
	    if {![regsub -- {^/\.\./} $path "/" path]} {break}
	}
    }
    while {1} {
	if {![regsub -- {^/\.\./} $path "/" path]} {break}
    }

    regsub -- {/$} $path "" path
    if {[tixStrEq $path ""]} {
	return /
    } else {
	return $path
    }
}

# tixFSCreateDirs
#
#
# 
proc tixFSCreateDirs {path} {
    tixAssert {
	error "Procedure tixFSCreateDirs not implemented on all platforms"
    }
    if {[tixStrEq $path /]} {
	return 1
    }
    if {[file exists $path]} {
	return 1
    }
    if {![tixFSCreateDirs [file dirname $path]]} {
	return 0
    }
    if {[catch {exec mkdir $path}]} {
	return 0
    }
    return 1
}

} else {

#-Win--------------------------------------------------------------------

# (Win) tixFSPWD --
#
#	Return the current directory
#
proc tixFSPWD {} {
    set p [pwd]
    regsub -all -- / $p \\ p
    return $p
}
# Win
#
proc tixFSIsNorm {path} {

    # Drive root directory
    # CYGNUS: drive can be immediately followed by directory separator.
    #
    if {[regexp -- {^[A-z]:\\?$} $path]} {
	return 1
    }

    # If it is not a drive root directory, it must
    # have a leading [drive letter:]\\[non empty string]
    # CYGNUS: A UNC path (\\host\dir) is also OK.
    if {![regexp -- {^[A-z]:\\.} $path]} {
	if {![regexp -- {^\\\\.*\\.} $path]} {
	return 0
    }
    }

    # relative path
    #
    if {[regexp -- {\\[.]$} $path]} {
	return 0
    }
    if {[regexp -- {\\[.][.]$} $path]} {
	return 0
    }
    if {[regexp -- {\\[.]\\} $path]} {
	return 0
    }
    if {[regexp -- {\\[.][.]\\} $path]} {
	return 0
    }
    if {[tixStrEq $path .]} {
	return 0
    }
    if {[tixStrEq $path ..]} {
	return 0
    }

    # Double slashes
    # CYGNUS: Double slashes at the front are OK.
    #
    if {[regexp -- {.\\\\} $path]} {
	return 0
    }

    # Trailing slashes
    #
    if {[regexp -- {[\\]$} $path]} {
	return 0
    }

    return 1
}

# (Win) tixFSNorm --
#
#	Interprets the user's input and return file information about this
#	input.
#
# Arguments:
#	See documentation (docs/Files.txt)
#
proc tixFSNorm {context text {defFile ""} {flagsVar ""} {errorMsgVar ""}} {
    tixAssert {
	tixAssertVPath $context
    }

    if {![tixStrEq $errorMsgVar ""]} {
	upvar $errorMsgVar errorMsg
    }
    if {![tixStrEq $flagsVar ""]} {
	upvar $flagsVar flags
    }

    set isDir [regexp -- {[\\]$} $text]
    set text [_tixJoin $context $text]
    set path [tixFSPath $text]

    if {$isDir || [file isdir $path]} {
	set vpath $text
	set tail $defFile
    } else {
	set list [split $text \\]
	set tail [lindex $list end]
	set len [string length $tail]
	set vpath [string range $text 0 [expr [string len $text]-$len-1]]
	regsub -- {[\\]$} $vpath "" vpath
    }

    set path [tixFSPath $vpath]

    if {![info exists flag(noPattern)]} {
	set isPat 0
	foreach char [split $tail ""] {
	    if {$char == "*" || $char == "?"} {
		set isPat 1
		break
	    }
	}
	if {$isPat} {
	    return [list $path $vpath "" $tail]
	}
    }

    return [list $path $vpath $tail ""]
}

# Win
#
# _tixJoin -- [internal]
#
#	Joins a pathname to a VPATH
#
proc _tixJoin {vp1 p2} {
    if {[tixFSIsAbsPath $p2]} {
	return [tixFSVPath [_tixNormalize $p2]]
    } else {
	return [tixFSVPath [_tixNormalize [tixFSPath $vp1]\\$p2]]
    }
}

# (Win) tixFSIsAbsPath
#
#	The Tcl "file pathtype" is buggy. E.g. C:\.\..\. is absolute, but
#	"file pathtype" thinks that it isn't
#

proc tixFSIsAbsPath {path} {
    # CYGNUS: Handle a UNC path (\\host\dir)
    if {[regexp -- {^\\\\.*\\.} $path]} {
	return 1
    }
    return [regexp -- {^[A-z]:\\} $path]
}

# (Win) tixFSIsNorm_os
#
#	Returns true iff this pathname is normalized, in the OS native name
#	format
#
proc tixFSIsNorm_os {path} {
    if {[regexp -- {^[A-z]:[\\]$} $path]} {
	return 1
    }
    if {[regexp -- {^[A-z]:$} $path]} {
	return 0
    }

    return [tixFSIsNorm $path]

}

# Win
#
# _tixNormalize --
#
#	Normalizes an absolute pathname.
#
# 	$dir must be an absolute native pathname
#
proc _tixNormalize {abpath} {
    tixAssert {
	if {![tixFSIsAbsPath $abpath]} {
	    error "\"$abpath\" must be an absolute pathname"
	}
    }

    if {![regexp -- {^[A-z]:} $abpath drive]} {
	tixPanic "\"$abpath\" does not contain a drive letter"
    }
    set drive [string toupper $drive]

    # CYGNUS: Handle UNC paths (\\host\dir)
    if {[regexp -- {^\\\\.*\\.} $abpath]} {
	set drive "\\"
	regsub -- {^\\} $abpath "" path
    } else {
	if {![regexp -- {^[A-z]:} $abpath drive]} {
	    tixPanic "\"$abpath\" does not contain a drive letter"
	}
	set drive [string toupper $drive]

	regsub -- {^[A-z]:} $abpath "" path
    }

    # Don't be fooled: $path doesn't need to be a directory. The following
    # code "set path $path\\" just makes it easy to get rid of trailing
    # . and ..
    #
    set path $path\\
    regsub -all -- {[\\]+} $path \\ path
    while {1} {
	if {![regsub -- {\\[.]\\} $path "\\" path]} {break}
    }
    while {1} {
	if {![regsub -- {\\[.]$} $path "" path]} {break}
    }

    while {1} {
	if {![regsub -- {\\[^\\]+\\[.][.]\\} $path "\\" path]} {break}
	while {1} {
	    if {![regsub -- {^\\[.][.]\\} $path "\\" path]} {break}
	}
    }
    while {1} {
	if {![regsub -- {^\\[.][.]\\} $path "\\" path]} {break}
    }

    regsub -- {[\\]+$} $path "" path
    return $drive$path
}

# Win
#
# tixFSNormDir --
#
#	Normalizes a directory
#
proc tixFSNormDir {dir} {
    if {![tixFSIsAbsPath $dir]} {
	error "\"$dir\" must be an absolute pathname"
    }
    if {![file isdir $dir]} {
	error "\"$dir\" is not a directory"
    }
    return [_tixNormalize $dir]
}


proc tixPanic {message} {
    error $message
}

# tixFSIsValid --
#
#	Checks whether a native pathname contains invalid characters.
#
proc tixFSIsValid {path} {
    return 1
}

# Win
#
#
proc tixFSIsVPath {vpath} {
    global tixPriv
    if {$tixPriv(isWin95)} {
	# CYGNUS: Accept UNC path (\\host\dir)
	if {[string match {xx\\xx\\\\\\*\\*} $vpath]} {
	    return 1
	}
	return [string match {xx\\xx\\[A-z]:*} $vpath]
    } else {
	return [string match {xx\\[A-z]:*} $vpath]
    }
}

# Win
#
# tixFSVPath --
#
#	Converts a normalized native pathname to its VPATH
#
proc tixFSVPath {path} {
    global tixPriv

    tixAssert {
	tixAssertNorm $path
    }
    return $tixPriv(WinPrefix)\\$path
}

# tixFSPath --
#
#	Converts a vpath to a native pathname
proc tixFSPath {vpath} {
    global tixPriv
    tixAssert {
	tixAssertVPath $vpath
    }
    if {$tixPriv(isWin95)} {
	set path [string range $vpath 6 end]
    } else {
	set path [string range $vpath 3 end]
    }
    regsub -- {:$} $path :\\ path

    return $path
}

# tixFSDisplayName --
#
#	Returns the name of a normalized path which is usually displayed by
#	the OS
#
proc tixFSDisplayName {normpath} {
    tixAssert {
	tixAssertNorm $normpath
    }

    if {[regexp -- {^[A-z]:$} $normpath]} {
	return $normpath\\
    } else {
	return $normpath
    }
}


tixInitFileCmpt:Win 

}
