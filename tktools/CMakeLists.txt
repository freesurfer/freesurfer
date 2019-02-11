project(tktools)

# The tk GUIs - tkmedit, tkregister2, tksurfer, and qdec - are no
# no longer maintained and cannot be built due to incompatibles across
# systems and with modern vtk (which no longer offers tcl wrapping).
# Eventually these programs can be completely deprecated, but for now
# we have to option to install precompiled binaries in linux distributions.

if(NOT APPLE)

  install_symlinks(
    qdec
    tkmedit
    tkregister2
    tksurfer
    TYPE programs
    DESTINATION bin
  )

  install(PROGRAMS
    scripts/tkanalyse.tcl
    scripts/tkm_common.tcl
    scripts/tkm_dialog.tcl
    scripts/tkmedit.tcl
    scripts/tkm_functional.tcl
    scripts/tkm_graph.tcl
    scripts/tkm_wrappers.tcl
    scripts/tkregister2.tcl
    scripts/tksurfer.tcl
    scripts/tkUtils.tcl
    DESTINATION tktools
  )

  install_tarball(resources.tar.gz DESTINATION lib)

endif()
