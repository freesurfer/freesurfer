# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: WinFile.tcl,v 1.4.2.1 2001/11/03 07:26:10 idiscovery Exp $
#
# WinFile.tcl --
#
#	MS Window file access portibility routines.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tixInitFileCmpt:Win {} {
    global tixPriv tcl_platform

    if {$tcl_platform(osVersion) >= 4.0} {
	set tixPriv(isWin95) 1
    } else {
	set tixPriv(isWin95) 0
    }
    if {$tixPriv(isWin95)} {
	set tixPriv(WinPrefix) xx\\xx
    } else {
	set tixPriv(WinPrefix) xx
    }

#----------------------------------------------------------------------
#
#		MS Windows
#
#----------------------------------------------------------------------

# splits a Windows directory into its hierarchical components
#
proc tixFSSplit {vpath} {
    global tixPriv

    set path ""
    if $tixPriv(isWin95) {
	if {![string compare $vpath xx]} {
	    lappend path [list xx     "Desktop"     "C:\\Windows\\Desktop" ]
	    return $path
	}
	if {![string compare $vpath xx\\xx]} {
	    lappend path [list xx     "Desktop"     "C:\\Windows\\Desktop" ]
	    lappend path [list xx\\xx "My Computer" "C:\\"]
	    return $path
	}

	set prefix "xx\\xx"
	if {![regsub -- {^xx\\xx\\} $vpath "" dir]} {
	    if {[regsub -- {^xx\\} $vpath "" dir]} {
		lappend path [list xx     "Desktop"     "C:\\Windows\\Desktop" ]
		set v "xx"
		set p "C:\\Windows\\Desktop"
		foreach d [split $dir \\] {
		    append v \\$d
		    append p \\$d
		    lappend path [list $v $d $p]
		}
		return $path
	    }
	}
	regsub -- {:$} $dir :/ dir
	lappend path [list xx     "Desktop"     "C:\\Windows\\Desktop" ]
	lappend path [list xx\\xx "My Computer" "C:\\"]
    } else {
	if {![string compare $vpath xx]} {
	    lappend path [list xx     "My Computer" "C:\\"]
	    return $path
	}
	lappend path [list xx     "My Computer" "C:\\"]

	set prefix xx
	regsub -- {^xx\\} $vpath "" dir
	regsub -- {:$} $dir :/ dir
    }

    if {![string compare $dir ""]} {
	return $path
    }
    if {[string compare [file pathtype $dir] "absolute"]} {
	error "$dir must be an absolute path"
    }

    set dirs [file split $dir]
    set p ""
    foreach d $dirs {
	set p [file join $p $d]
	regsub -all / $p \\ p
	set vpath $prefix\\$p
	regsub -- {[\\]$} $vpath "" vpath
	regsub -- {:/$} $d ":" d
	lappend path [list $vpath $d $p]
    }

    return $path
}

# returns true if $dir is an valid path (not equal to "")
#
proc tixFSValid {dir} {
    return [expr ![string compare $dir ""]]
}

# tixFSIntName
#
#	Returns the "virtual path" of a filename
#
proc tixFSIntName {dir} {
    global tixPriv

    if {![string compare $dir ""]} {
	if $tixPriv(isWin95) {
	    return "xx\\xx"
	} else {
	    return xx
	}
    }
        
    if {[string compare [file pathtype $dir] "absolute"]} {
	error "$dir must be an absolute path"
    }

    if $tixPriv(isWin95) {
        set vpath "xx\\xx\\$dir"
    } else {
        set vpath "xx\\$dir"
    }
    regsub -- {:/$} $vpath ":" vpath
    regsub -- {[\\]$} $vpath "" vpath
    return $vpath
}

proc tixFSIntJoin {dir sub} {
    set vpath $dir\\$sub
    regsub -all {\\\\} $vpath \\ vpath
    regsub -- {:/$} $vpath : vpath
    regsub -- {[\\]$} $vpath "" vpath
    return $vpath
}

proc tixFSJoin {dir sub} {
    set p [file join $dir $sub]
    regsub -all / $p \\ p
    return $p
}

proc tixFSResolveName {p} {
    regsub -all / $p \\ p
    if {[regexp -- {:([^\\]|$)} $p]} {
	regsub : $p :\\ p
    }
    return $p
}

# dir:		Make a listing of this directory
# showSubDir:	Want to list the subdirectories?
# showFile:	Want to list the non-directory files in this directory?
# showPrevDir:	Want to list ".." as well?
# showHidden:	Want to list the hidden files? (%% is ignored)
#
# return value:	a list of files and/or subdirectories
#
proc tixFSListDir {vpath showSubDir showFile showPrevDir showHidden {pattern ""}} {
    global tixPriv
    set appPWD [pwd]
    set list ""

    if $tixPriv(isWin95) {
	if {![string compare $vpath xx]} {
	    set dir C:\\Windows\\Desktop
	    if {$showSubDir} {
		lappend list xx:
	    }
	} elseif {![string compare $vpath xx\\xx]} {
	    if {$showSubDir} {
		return [tixFSGetDrives]
	    } else {
		return ""
	    }
	} else {
	    if {![regsub -- {^xx\\xx\\} $vpath "" dir]} {
		regsub -- {^xx\\} $vpath C:\\Windows\\Desktop\\ dir
	    }
	    regsub -- {:$} $dir :\\ dir
	}
    } else {
	if {![string compare $vpath xx]} {
	    if {$showSubDir} {
		return [tixFSGetDrives]
	    } else {
		return ""
	    }
	}

	regsub -- {^xx\\} $vpath "" dir
	regsub -- {:$} $dir :\\ dir
    }

    if {[catch {cd $dir} err]} {
	# The user has entered an invalid directory
	# %% todo: prompt error, go back to last succeed directory
	cd $appPWD
	return ""
    }

    if {$pattern == ""} {
	set pattern "*"
    }

    if {[catch {set names [lsort [eval glob -nocomplain $pattern]]} err]} {
	# Cannot read directory
	# %% todo: show directory permission denied
	cd $appPWD
	return ""
    }

    catch {
	# We are catch'ing, just in case the "file" command returns unexpected
	# errors
	#
	foreach fname $names {
	    if {![string compare . $fname]} {
		continue
	    }
	    if {![string compare ".." $fname]} {
		continue
	    }
	    if {[file isdirectory $fname]} {
		if $showSubDir {
		    lappend list [file tail $fname]
		}
	    } else {
		if $showFile {
		    lappend list [file tail $fname]
		}
	    }
	}
    }
    cd $appPWD

    if {$showSubDir && $showPrevDir && $dir != "/"} {
	return [tixFSMakeList $vpath $dir [lsort [concat .. $list]]]
    } else {
	return [tixFSMakeList $vpath $dir $list]
    }
}

proc tixFSMakeList {vpath dir list} {
    global tixPriv

    if $tixPriv(isWin95) {
	set prefix xx\\xx
    } else {
	set prefix xx
    }
    set l ""
    foreach file $list {
	if {![string compare $file xx:]} {
	     lappend l [list xx\\xx "My Computer" "C:\\"]
	} else {
	    set path [tixFSJoin $dir $file]
	    lappend l [list $vpath\\$file $file $path]
	}
    }

    return $l
}

proc tixFSSep {} {
    return "\\"
}

proc tixFSGetDrives {} {
    global tixPriv

    if {[info exists tixPriv(drives)]} {
	return $tixPriv(drives)
    } else {
	set drives [list A: B:]
	foreach d {c d e f g h i j k l m n o p q r s t u v w x y z} {
	    if {[file exists $d:\\]} {
		lappend drives [string toupper $d:]
	    }
	}

	set tixPriv(drives) ""
	foreach d $drives {
	     lappend tixPriv(drives) [list $tixPriv(WinPrefix)\\$d $d $d\\]
	}
    }
    return $tixPriv(drives)
}

#----------------------------------------------------------------------
#
#		OBSOLETE
#
#----------------------------------------------------------------------



# Directory separator
#
proc tixDirSep {} {
    return "\\"
}

# returns the "root directory" of this operating system
#
# out:	intName
proc tixRootDir {} {
    return "/"
}

# is an absoulte path only if it starts with a baclskash
# or starts with "<drive letter>:"
#
# in: nativeName
#
proc tixIsAbsPath {nativeName} {
    set c [string index $nativeName 0]
    if {$c == "\\"} {
	return 1
    }

    if {[string compare [string toupper $c] A] < 0} {
	return 0
    }
    if {[string compare [string toupper $c] Z] > 0} {
	return 0
    }
    if {[string index $nativeName 1] != ":"} {
	return 0
    }
    return 1
}

# returns <drive>:
#
proc tixWinGetFileDrive {nativeName} {
    set c [string index $nativeName 0]
    if {$c == "\\"} {
	return [string toupper [string range [pwd] 0 1]]
    }

    if {[string compare [string toupper $c] A] < 0} {
	return [string toupper [string range [pwd] 0 1]]
    }
    if {[string compare [string toupper $c] Z] > 0} {
	return [string toupper [string range [pwd] 0 1]]
    }
    if {[string index $nativeName 1] != ":"} {
	return [string toupper [string range [pwd] 0 1]]
    }
    return [string toupper [string range $nativeName 0 1]]
}

# returns the absolute pathname of the file 
# (not including the drive letter or the first backslash)
#
# [tixWinGetFileDrive]\\[tixWinGetFilePath] gives the complete
# drive and pathname
#
proc tixWinGetFilePath {nativeName} {
    set c [string index $nativeName 0]
    if {$c == "\\"} {
	return ""
    }

    if {[string compare [string toupper $c] A] < 0} {
	return [tixWinGetPathFromDrive $nativeName]
    }
    if {[string compare [string toupper $c] Z] > 0} {
	return [tixWinGetPathFromDrive $nativeName]
    }
    if {[string index $nativeName 1] != ":"} {
	return [tixWinGetPathFromDrive $nativeName]
    }
    if {[string index $nativeName 2] != "\\"} {
        regexp -- {[A-z]:} $nativeName drive
	regsub -- {[A-z]:} $nativeName "" path
	return [tixWinGetPathFromDrive $path $drive]
    }

    regsub -- {[A-z]:[\\]} $nativeName "" path
    return $path
}

proc tixWinCurrentDrive {} {
    return [string range [pwd] 0 1]
}

proc tixWinGetPathFromDrive {path {drive ""}} {
    if {$drive == ""} {
        set drive [tixWinCurrentDrive]
    }

    #
    # %% currently TCL (7.5b3) does not tell what the current path
    #    on a particular drive is

    return $path
}

#
#
# nativeName:	native filename used in this OS, comes from the user or
#		application programmer
# defParent:	(intName) if the filename is not an absolute path,
#		treat it as a subfolder of $defParent
#		(must be an intName, must be absolute)
proc tixFileIntName {nativeName {defParent ""}} {
    if {![tixIsAbsPath $nativeName]} {
        if {$defParent != ""} {
	    if {[string index $defParent 0] != "/"} {
	        error "Tix toolkit error: \"$defParent\" is not an absolute internal file name"
	    }
	    set path [tixSubFolder $defParent $nativeName]
	} else {
	    set path $nativeName
	}
    } else {
	set path /[tixWinGetFileDrive $nativeName]\\[tixWinGetFilePath $nativeName]
    }

    set intName ""
    foreach name [tixFileSplit $path] {
	set intName [tixSubFolder $intName $name]
    }

    return $intName
}

# in:	internal name
# out:	native name
proc tixNativeName {intName {mustBeAbs 1}} {
    if {[string index $intName 0] != "/"} {
        if {$mustBeAbs} {
            error "Tix internal error: \"$intName\" is not an intName"
	} else {
	    return $intName
	}
    }
    if {$intName == "/"} {
        return C:\\
    }
    regsub -- {/[\\]} $intName "" nativeName
    if {[string length $nativeName] == 2} {
        return $nativeName\\
    } else {
        return $nativeName
    }
}

# how a filename should be displayed
# 
# e.g. /\C: becomes C:\\
#      /\   becomes "My Computer"
#      /\C:\\Windows is Windows
proc tixFileDisplayName {intName} {
    if {[string index $intName 0] != "/"} {
        error "Tix internal error: \"$intName\" is not an intName"
    }

    if {$intName == "/"} {
        return "My Computer"
    }

    regsub -- {/[\\]} $intName "" nativeName

    if {[string length $nativeName] == 2} {
        return [string toupper $nativeName\\]
    } else {
        return [file tail $nativeName]
    }
}

# in:	internal name
# out:	a list of paths
proc tixFileSplit {intName} {

    set l ""
    foreach n [split $intName /\\] {
	if {$n == ""} {
	    continue
	}
	if {$n == "."} {
	    continue
	}

	lappend l $n
    }
    

    while {1} {
	set idx [lsearch $l ".."]
	if {$idx == -1} {
	    break;
	}
	set l [lreplace $l [expr $idx -1] $idx]
    }


    if {[string index $intName 0] == "/"} {
	return [concat "/" $l]
    } else {
	return $l
    }
}

# parent, sub:	intName
#
proc tixSubFolder {parent sub} {
    if {$parent == ""} {
	return $sub
    }
    return $parent\\$sub
}

proc tixWinGetDrives {} {
    global tixPriv

    if {[info exists tixPriv(drives)]} {
	return $tixPriv(drives)
    } else {
	set tixPriv(drives) {A: B:}
        foreach d {c e d f g h i j k l m n o p q r s t u v w x y z} {
	    if {[file exists $d:]} {
		lappend tixPriv(drives) [string toupper $d:]
	    }
        }
    }
    return $tixPriv(drives)
}

# dir:		Make a listing of this directory
# showSubDir:	Want to list the subdirectories?
# showFile:	Want to list the non-directory files in this directory?
# showPrevDir:	Want to list ".." as well?
# showHidden:	Want to list the hidden files? (%% is ignored)
#
# return value:	a list of files and/or subdirectories
#
proc tixListDir {dir showSubDir showFile showPrevDir showHidden {pattern ""}} { 
    set appPWD [pwd]

    if {$dir == "/"} {
	if {$showSubDir} {
	    return [tixWinGetDrives]
        } else {
	    return ""
	}
    }

    if {[catch {cd [tixNativeName $dir]} err]} {
	# The user has entered an invalid directory
	# %% todo: prompt error, go back to last succeed directory
	cd $appPWD
	return ""
    }

    if {$pattern == ""} {
	set pattern "*"
    }

    if {[catch {set names [lsort [eval glob -nocomplain $pattern]]} err]} {
	# Cannot read directory
	# %% todo: show directory permission denied
	cd $appPWD
	return ""
    }

    set list ""
    catch {
	# We are catch'ing, just in case the "file" command returns unexpected
	# errors
	#
 	foreach fname $names {
	    if {![string compare . $fname]} {
		continue
	    }
 	    if {![string compare ".." $fname]} {
	        continue
	    }
	    if {[file isdirectory $fname]} {
		if $showSubDir {
		    lappend list [file tail $fname]
		}
	    } else {
		if $showFile {
		    lappend list [file tail $fname]
		}
	    }
	}
    }
    cd $appPWD

    if {$showSubDir && $showPrevDir && $dir != "/"} {
	return [lsort [concat .. $list]]
    } else {
        return $list
    }
}

proc tixVerifyFile {file} {
    return [tixFileIntName $file]
}

proc tixFilePattern {args} {
    if {[lsearch $args allFiles] != -1} {
	return *
    }
    return *
}

}

# tixWinFileEmu --
#
#	Emulates a MS Windows file system environemnt inside Unix
#
proc tixWinFileEmu {} {
    cd /mnt/c
    rename pwd __pwd
    rename cd  __cd
    proc EmuConvert {path} {
	if {[regsub ^/mnt/c/ $path c:/ path]} {
	    return $path
	}
	if {[regsub ^/mnt/d/ $path d:/ path]} {
	    return $path
	}
	if {[regsub ^/mnt/c\$ $path c:/ path]} {
	    return $path
	}
	if {[regsub ^/mnt/d\$ $path d:/ path]} {
	    return $path
	}
	return c:/windows
    }

    proc pwd {} {
	return [EmuConvert [__pwd]]
    }
    proc glob {args} {

    }
}
