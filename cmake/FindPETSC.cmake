# PETSC Find Module

if(NOT PETSC_DIR)
  set(PETSC_DIR ${FS_PACKAGES_DIR}/petsc/2.3.3)
endif()

find_path(PETSC_INCLUDE_DIR HINTS ${PETSC_DIR} NAMES petsc.h PATH_SUFFIXES include)

find_package_handle_standard_args(PETSC DEFAULT_MSG PETSC_INCLUDE_DIR)

library_paths(
  NAME PETSC_LIBRARIES
  LIBDIR ${PETSC_DIR}/lib
  LIBRARIES
  petscts
  petscsnes
  petscksp
  petscdm
  petscmat
  petscvec
  petsc
  petsccontrib
  mpich
  pmpich
)
