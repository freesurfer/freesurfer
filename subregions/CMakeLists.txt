project(subregions)

if(NOT INTEGRATE_SAMSEG)
  install_pyscript(segment_subregions)
endif()  

install_symlinks(ConcatenateSubregionsResults.sh TYPE programs DESTINATION bin)
