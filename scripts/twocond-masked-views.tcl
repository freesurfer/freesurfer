if [info exists noexit] { unsetenv condition0 }
set mask ${hemi}-cortex
source $env(MRI_DIR)/lib/tcl/twocond-views.tcl
#mask_label ${hemi}-cortex
