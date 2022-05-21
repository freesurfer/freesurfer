# -*- mode: TCL; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: UnixFile.tcl,v 1.4.2.1 2001/11/03 07:25:12 idiscovery Exp $
#
# UnixFile.tcl --
#
#	Unix file access portibility routines.
#
# Copyright (c) 1993-1999 Ioi Kim Lam.
# Copyright (c) 2000-2001 Tix Project Group.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc tixInitFileCmpt:Unix {} {

# tixFSSplit --
# 
# Splits a directory into its hierarchical components
#
# "hlist-type hierachical path"		<- "vpath"
# "name"
# "directory name"			<- "path"
#
proc tixFSSplit {dir} {
    if {[string compare [tixFSPathType $dir] "absolute"]} {
	error "$dir must be an absolute path"
    }

    set path ""
    set p ""
    foreach d [tixFileSplit $dir] {
	set p [tixFSJoin $p $d]
	lappend path [list $p $d $p]
    }
    return $path
}

# returns true if $dir is an valid path (always true in Unix)
#
proc tixFSValid {dir} {
    return 1
}

# Directory separator
#
proc tixFSSep {} {
    return "/"
}

# tixFSIntName
#
#	Returns the "virtual path" of a filename
#
proc tixFSIntName {dir} {
    if {[string compare [tixFSPathType $dir] "absolute"]} {
	error "$dir must be an absolute path"
    }

    return $dir
}

proc tixFSResolveName {p} {
    return $p
}


# These subcommands of "file" only exist in Tcl 7.5+. We define the following
# wrappers so that the code also works under Tcl 7.4
#
global tcl_version
if {![string compare $tcl_version 7.4]} {

    proc tixFSPathType {dir} {
	if {![string compare [string index $dir 0] /]} {
	    return "absolute"
	}
	if {![string compare [string index $dir 0] ~]} {
	    return "absolute"
	}
	return "relative"
    }

    proc tixFSJoin {dir sub} {
	set joined $dir/$sub

	regsub -all {[/]+} $joined / joined
	return $joined
    }

} else {
    proc tixFSPathType {dir} {
	return [file pathtype $dir]
    }

    proc tixFSJoin {dir sub} {
	return [file join $dir $sub]
    }
}

# dir:		Make a listing of this directory
# showSubDir:	Want to list the subdirectories?
# showFile:	Want to list the non-directory files in this directory?
# showPrevDir:	Want to list ".." as well?
# showHidden:	Want to list the hidden files?
#
# return value:	a list of files and/or subdirectories
#
proc tixFSListDir {dir showSubDir showFile showPrevDir showHidden {pattern ""}} {
    set appPWD [pwd]

    if {[catch {cd $dir} err]} {
	# The user has entered an invalid directory
	# %% todo: prompt error, go back to last succeed directory
	cd $appPWD
	return ""
    }

    if {$pattern == ""} {
	if $showHidden {
	    set pattern "* .*"
	} else {
	    set pattern *
	}
    } elseif {$pattern == "*"} {
	if $showHidden {
	    set pattern "* .*"
	}
    }

    set list ""
    foreach pat $pattern {
	if {[catch {set names [lsort [glob -nocomplain $pat]]} err]} {
	    # Cannot read directory
	    # %% todo: show directory permission denied
	    continue
	}

	catch {
	    # We are catch'ing, just in case the "file" command
	    # returns unexpected errors
	    #
	    foreach fname $names {
		if {![string compare . $fname]} {
		    continue
		}
		if {[file isdirectory $fname]} {
		    if {![string compare ".." $fname] && !$showPrevDir} {
			continue
		    }
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
    }

    cd $appPWD

    if {[llength $pattern] > 1} {
	# get rid of duplicated names
	#
	set list1 ""
	set oldfile ""
	foreach name [lsort $list] {
	    if {$name == $oldfile} {
		continue
	    }
	    lappend list1 $name
	    set oldfile $name
	}
	return [_tixFSMakeList $dir $list1]
    } else {
	return [_tixFSMakeList $dir $list]
    }
}

# _tixFSMakeList -
#
#	Internal procedure. Used only by tixFSListDir
proc _tixFSMakeList {dir list} {
    set l ""
    foreach file $list {
	set path [tixFSJoin $dir $file]
	lappend l [list $path $file $path]
    }

    return $l
}

# Directory separator
#
proc tixDirSep {} {
    return "/"
}


# tixFSInfo --
#
#	Returns information about the file system of this OS
#
# hasdrives: Boolean
#	Does this file system support seperate disk drives?
#
proc tixFSInfo {args} {
    case [lindex $args 0] {
	hasdrives {
	    return 0
	}
    }
}

#----------------------------------------------------------------------
# Obsolete
#----------------------------------------------------------------------

# nativeName:	native filename used in this OS, comes from the user or
#		application programmer
# defParent:	if the filename is not an absolute path, treat it as a
#		subfolder of $defParent
proc tixFileIntName {nativeName {defParent ""}} {
    if {![tixIsAbsPath $nativeName]} {
	if {$defParent != ""} {
	    set path [tixSubFolder $defParent $nativeName]
	} else {
	    set path $nativeName
	}
    } else {
	set path $nativeName
    }

    set intName ""
    set path [tixFile trimslash [tixFile tildesubst $path]]
    foreach name [tixFileSplit $path] {
	set intName [tixSubFolder $intName $name]
    }
    return $intName
}

proc tixNativeName {name {mustBeAbs ""}} {
    return $name
}

proc tixFileDisplayName {intName} {
    if {$intName == "/"} {
	return "/"
    } else {
	return [file tail $intName]
    }
}


proc tixFileSplit {intName} {

    set l ""
    foreach n [split $intName /] {
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

proc tixSubFolder {parent sub} {
    if {$parent == ""} {
	return $sub
    }
    if {$parent == "/"} {
	return /$sub
    } else {
	return $parent/$sub
    }
}

# dir:		Make a listing of this directory
# showSubDir:	Want to list the subdirectories?
# showFile:	Want to list the non-directory files in this directory?
# showPrevDir:	Want to list ".." as well?
# showHidden:	Want to list the hidden files?
#
# return value:	a list of files and/or subdirectories
#
proc tixListDir {dir showSubDir showFile showPrevDir showHidden {pattern ""}} { 

    set appPWD [pwd]

    if {[catch {cd $dir} err]} {
	# The user has entered an invalid directory
	# %% todo: prompt error, go back to last succeed directory
	cd $appPWD
	return ""
    }

    if {$pattern == ""} {
	if $showHidden {
	    set pattern "* .*"
	} else {
	    set pattern *
	}
    } elseif {$pattern == "*"} {
	if $showHidden {
	    set pattern "* .*"
	}
    }

    set list ""
    foreach pat $pattern {
	if {[catch {set names [lsort [glob -nocomplain $pat]]} err]} {
	    # Cannot read directory
	    # %% todo: show directory permission denied
	    continue
	}

	catch {
	    # We are catch'ing, just in case the "file" command
	    # returns unexpected errors
	    #
	    foreach fname $names {
		if {![string compare . $fname]} {
		    continue
		}
		if {[file isdirectory $fname]} {
		    if {![string compare ".." $fname] && !$showPrevDir} {
			continue
		    }
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
    }

    cd $appPWD

    if {[llength $pattern] > 1} {
	set list1 ""
	set oldfile ""
	foreach name [lsort $list] {
	    if {$name == $oldfile} {
		continue
	    }
	    lappend list1 $name
	    set oldfile $name
	}
	return $list1
    } else {
	return $list
    }
}

# returns the "root directory" of this operating system
#
proc tixRootDir {} {
    return "/"
}

proc tixIsAbsPath {nativeName} {
    set c [string index $nativeName 0]
    if {$c == "~" || $c == "/"} {
	return 1
    } else {
	return 0
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





