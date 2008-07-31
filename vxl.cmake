IF (WIN32)
  FIND_PATH(vxl_ROOT_DIR 
    NAMES lib/vcl.lib
    PATHS 
      $ENV{vxl_ROOT_DIR}
      $ENV{ProgramFiles}/gnuwin32/
      C:/
      D:/
      $ENV{ProgramFiles}
    PATH_SUFFIXES 
      vxl-1.10.0
    vxl-1.10
      vxl
    DOC "vxl base/installation directory?"
    )
ENDIF (WIN32)


IF (WIN32)
  FIND_PATH(vxl_ROOT_INCLUDE_DIR 
    NAMES core/vxl_version.h
    PATHS
      $ENV{vxl_ROOT_INCLUDE_DIR}
    ${vxl_ROOT_DIR}/include
    ${vxl_ROOT_DIR}
    PATH_SUFFIXES 
      vxl
    DOC "vxl include directory?"
    )
ENDIF (WIN32)


FIND_LIBRARY (VNL_LIBRARIES 
  NAMES vnl vnl.lib
  PATHS
    ${vxl_ROOT_DIR}/lib
    $ENV{DEFAULT_LIB}
    $ENV{LIBRARY_PATH}
) 
MARK_AS_ADVANCED(VNL_LIBRARY)

FIND_LIBRARY (VNL_ALGO_LIBRARIES 
  NAMES vnl_algo 
  PATHS
    ${vxl_ROOT_DIR}/lib
    $ENV{DEFAULT_LIB}
    $ENV{LIBRARY_PATH}
)
MARK_AS_ADVANCED(VNL_ALGO_LIBRARY)

FIND_LIBRARY (VCL_LIBRARIES
  NAMES vcl
  PATHS
    ${vxl_ROOT_DIR}/lib
    $ENV{DEFAULT_LIB}
    $ENV{LIBRARY_PATH}
)
MARK_AS_ADVANCED(VCL_LIBRARY)

FIND_LIBRARY (VPL_LIBRARIES 
  NAMES vpl
  PATHS
    ${vxl_ROOT_DIR}/lib
    $ENV{DEFAULT_LIB}
    $ENV{LIBRARY_PATH}
)
MARK_AS_ADVANCED(VPL_LIBRARY)

FIND_LIBRARY (V3P_LIBRARIES 
  NAMES v3p_netlib
  PATHS
    ${vxl_ROOT_DIR}/lib
    $ENV{DEFAULT_LIB}
    $ENV{LIBRARY_PATH}
)
MARK_AS_ADVANCED(V3P_LIBRARY)


SET (VXL_INCLUDES
  ${vxl_ROOT_INCLUDE_DIR}
  ${vxl_ROOT_INCLUDE_DIR}/v3p
  ${vxl_ROOT_INCLUDE_DIR}/v3p/netlib
  ${vxl_ROOT_INCLUDE_DIR}/vcl
  ${vxl_ROOT_INCLUDE_DIR}/core
  ${vxl_ROOT_INCLUDE_DIR}/core/vpl
  ${vxl_ROOT_INCLUDE_DIR}/core/vnl
  ${vxl_ROOT_INCLUDE_DIR}/core/vnl/algo
)
MARK_AS_ADVANCED(VXL_INCLUDES)

SET (VXL_LIBS 
  ${VNL_LIBRARIES} 
  ${VNL_ALGO_LIBRARIES} 
  ${V3P_LIBRARIES} 
  ${VCL_LIBRARIES} 
  ${VPL_LIBRARIES}
)
