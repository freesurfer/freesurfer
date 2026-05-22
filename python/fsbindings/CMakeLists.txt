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

## While python 3.8 can be installed on Ubuntu24, switching the
## system to py38 causes other cmake conflicts where 3.12 (the
## default native version) is expected for the tool chain.

# if(HOST_OS MATCHES "Ubuntu24")
#   set(PYTHON_EXECUTABLE "/usr/bin/python3.8" CACHE FILEPATH "" FORCE)
#   set(PYBIND11_PYTHON_VERSION "3.8" CACHE STRING "" FORCE)
#   find_package(pybind11 REQUIRED)
# endif()

# link utilities
target_link_libraries(fsbindings PRIVATE utils nifti ${ITK_LIBRARIES})

if(LINK_SHARED_AS_NEEDED)
   target_link_options(fsbindings PRIVATE "LINKER:--as-needed")
endif()

# make sure the bindings library gets built into the repository even in out-of-source builds
set_target_properties(fsbindings PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/fsbindings)

## the workaround for ubuntu 24 is to copy the python3.8 fsbindings libraries as built on and saved from an ubuntu 22 system

if(HOST_OS MATCHES "Ubuntu24")
   set(VAST_BASE "/space/freesurfer/ubuntu22_x86_64/dev/python/packages/fsbindings")
   if("${FSPYTHON_VERSION}" MATCHES "3.8")
      set(FSBIND_LIB "fsbindings.cpython-38-x86_64-linux-gnu.so")
   elseif("${FSPYTHON_VERSION}" MATCHES "3.9")
      set(FSBIND_LIB "fsbindings.cpython-39-x86_64-linux-gnu.so")
   endif()
   add_custom_command(TARGET fsbindings POST_BUILD COMMAND
     cp -p -f ${VAST_BASE}/${FSBIND_LIB} ${CMAKE_CURRENT_SOURCE_DIR}/. && ls -lt ${CMAKE_CURRENT_SOURCE_DIR}/${FSBIND_LIB}
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
     COMMENT "Copy in ubuntu22 fsbindings libraries for ubuntu24"
   )
   install(FILES ${FSBIND_LIB} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION ${FSPYTHON_INSTALL_PREFIX}/python/packages/fsbindings)
   # remove when "make clean" is run - wildcard expressions do not work
   set_property(TARGET fsbindings APPEND PROPERTY ADDITIONAL_CLEAN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${FSBIND_LIB}")
endif()

# PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE

