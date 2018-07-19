# PETSC Find Module

if(NOT PETSC_DIR)
  set(PETSC_DIR ${FS_PACKAGES_DIR}/petsc/2.3.3-p13)
endif()

find_path(PETSC_INCLUDE_DIR PATHS ${PETSC_DIR} NAMES petsc.h PATH_SUFFIXES include)

# the petsc libs we want
set(LIBS petscts petscsnes petscksp petscdm petscmat petscvec petsc petsccontrib fmpich pmpich mpich)

# find the (many) libraries
foreach(LIB ${LIBS})
  find_library(tmp PATHS ${PETSC_DIR} NAMES lib${LIB}.a PATH_SUFFIXES lib)
  set(PETSC_LIBRARIES ${PETSC_LIBRARIES} ${tmp})
  unset(tmp CACHE)  # this is necessary for find_library to work (plus it clears it from the cache)
endforeach()

find_package_handle_standard_args(PETSC DEFAULT_MSG PETSC_INCLUDE_DIR PETSC_LIBRARIES)

# make sure the petsc path gets linked from the lib directory during install
symlink(${PETSC_DIR} ${CMAKE_INSTALL_PREFIX}/lib/petsc)
