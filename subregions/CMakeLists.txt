project(subregions)

if(NOT INTEGRATE_SAMSEG)
   message(STATUS "No install of segment_subregions INTEGRATE_SAMSEG=OFF")
   # if (FSPYTHON_INSTALL_TREE)
   #    install_pyscript_fspython_tree(segment_subregions)
   # else()
   #    install_pyscript(segment_subregions)
   # endif()
endif()  

install_symlinks(ConcatenateSubregionsResults.sh TYPE programs DESTINATION bin)
