# Tcl package index file, version 1.0

proc LoadBLT { version dir } {

    set prefix "lib"
    set suffix [info sharedlibextension]
    regsub {\.} $version {} version_no_dots

    # Determine whether to load the full BLT library or
    # the "lite" tcl-only version.
    
    if { [info commands tk] == "tk" } {
        set name libBLT.2.dylib
    } else {
        set name libBLTlite.2.dylib
    }
    
    global tcl_platform
    if { $tcl_platform(platform) == "unix" } {
	set library [file join $dir $name]
	if { ![file exists $library] } {
	    # Try the parent directory.
	    set library [file join [file dirname $dir] $name]
	}
	if { ![file exists $library] } {
	    # Default to the path generated at compilation.
	    set library [file join "/usr/pubsw/packages/tcltktixblt/current/lib" $name]
	}
    } else {
	set library $name
    }
    load $library BLT
}

set version "2.4"

package ifneeded BLT $version [list LoadBLT $version $dir]

# End of package index file
