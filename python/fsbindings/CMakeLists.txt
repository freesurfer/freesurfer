include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${FS_INCLUDE_DIRS})

# temporary fix so that -g doesn't produce linker errors when binding cxx/python code
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")

pybind11_add_module(fsbindings
  module.cpp
  bindings_numpy.cpp
  bindings_surface.cpp
  bindings_mri.cpp
  bindings_transform.cpp 
)

# link utilities
target_link_libraries(fsbindings PRIVATE utils nifti ${ITK_LIBRARIES})

if(LINK_SHARED_AS_NEEDED)
   target_link_options(fsbindings PRIVATE "LINKER:--as-needed")
endif()

# make sure the bindings library gets built into the repository even in out-of-source builds
set_target_properties(fsbindings PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/fsbindings)

## FIX ME - this workaround copies in fsbindings libraries that work with an fspython 3.8 from an ubuntu 22 system
## on an ubuntu 24 system where the installed python is an newer/different version, i.e., 3.12.

if(HOST_OS MATCHES "Ubuntu24")
   add_custom_command(TARGET fsbindings POST_BUILD COMMAND
     # cp -p -f /space/freesurfer/ubuntu22_x86_64/dev/python/packages/fsbindings/*.cpython-38-x86_64-linux-gnu.so ${CMAKE_CURRENT_SOURCE_DIR}/. && ls -lt *.cpython-38-x86_64-linux-gnu.so
     cp -p -f /space/freesurfer/ubuntu22_x86_64/8.2.0/python/packages/fsbindings/*.cpython-38-x86_64-linux-gnu.so ${CMAKE_CURRENT_SOURCE_DIR}/. && ls -lt *.cpython-38-x86_64-linux-gnu.so
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
     COMMENT "Copy in ubuntu22 fsbindings libraries for ubuntu24"
   )
   # remove when "make clean" is run - wildcard expressions do not work
   set_property(TARGET fsbindings APPEND PROPERTY ADDITIONAL_CLEAN_FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/fsbindings.cpython-38-x86_64-linux-gnu.so"
      "${CMAKE_CURRENT_SOURCE_DIR}/labelfusion.cpython-38-x86_64-linux-gnu.so"
   )
endif()

