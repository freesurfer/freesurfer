include(FindPackageHandleStandardArgs)

set(MARTINOS_PETSC /usr/pubsw/packages/petsc/2.3.3-p13)

# find the include dir
find_path(PETSC_INCLUDE_DIR PATHS ${PETSC_DIR} ${MARTINOS_PETSC} NAMES petsc.h PATH_SUFFIXES include)
# find the libraries
foreach(lib petscts petscsnes petscksp petscdm petscmat petscvec petsc petsccontrib fmpich pmpich mpich)
  find_library(${lib}_LIBRARY PATHS ${PETSC_DIR} ${MARTINOS_PETSC} NAMES lib${lib}.a PATH_SUFFIXES lib)
  set(PETSC_LIBRARIES ${PETSC_LIBRARIES} ${${lib}_LIBRARY})
endforeach()
find_package_handle_standard_args(PETSC DEFAULT_MSG PETSC_LIBRARIES PETSC_INCLUDE_DIR)
