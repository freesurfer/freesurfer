#
# Find the Petsc includes and dependent libraries
#
# This works assuming gnu is used.
#
# PETSC_DIR	
# PETSC_LINK_DIR
# PETSC_LIBS
# PETSC_FOUND

IF(NOT PETSC_DIR)
	SET(PETSC_DIR $ENV{PETSC_DIR} CACHE PATH "path to Petsc files")
	#SET(MPI_DIR $ENV{MPI_DIR} CACHE PATH "path to MPI")
	
	IF(EXISTS ${PETSC_DIR})
		SET(PETSC_DIR_FOUND yes)
	ENDIF(EXISTS ${PETSC_DIR})

	IF(PETSC_DIR_FOUND)
		SET(PETSC_FOUND yes)
	ENDIF(PETSC_DIR_FOUND)
ELSE (NOT PETSC_DIR)	
	SET(PETSC_FOUND yes)
ENDIF(NOT PETSC_DIR)

IF(PETSC_FOUND)
	SET(PETSC_INCLUDE_DIR ${PETSC_DIR}/include ${PETSC_DIR}/include/mpiuni
	)
	SET(PETSC_LINK_DIR ${PETSC_DIR}/lib)
	SET(PETSC_LIBS  petscts petscsnes petscksp petscdm petscmat petscvec petsc
			mpiuni nsl aio rt blas lapack)

ENDIF(PETSC_FOUND)
