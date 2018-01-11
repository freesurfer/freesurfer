# note: right now this ONLY compiles the libutils.a library on
# a martinos centOS7 machine

project(freesurfer)
cmake_minimum_required(VERSION 3.5)

include_directories(${CMAKE_SOURCE_DIR}/include
                    ${CMAKE_SOURCE_DIR}/include/dicom
                    ${CMAKE_SOURCE_DIR}/include/NrrdIO
                    ${CMAKE_SOURCE_DIR}/xml2
                    ${CMAKE_SOURCE_DIR}/jpeg
                    ${CMAKE_SOURCE_DIR}/tiff
                    ${CMAKE_SOURCE_DIR}/expat)

include_directories(/usr/pubsw/packages/mni/1.5/include
                    /usr/pubsw/packages/vxl/current/include/vxl/core
                    /usr/pubsw/packages/vxl/current/include/vxl/vcl
                    /usr/pubsw/packages/vxl/current/include/vxl/v3p/netlib
                    /usr/pubsw/packages/vxl/current/include/vxl/v3p/netlib/opt)

include_directories(/usr/pubsw/packages/itk/current/include/InsightToolkit
                    /usr/pubsw/packages/itk/current/include/InsightToolkit/Algorithms
                    /usr/pubsw/packages/itk/current/include/InsightToolkit/BasicFilters
                    /usr/pubsw/packages/itk/current/include/InsightToolkit/Common
                    /usr/pubsw/packages/itk/current/include/InsightToolkit/IO
                    /usr/pubsw/packages/itk/current/include/InsightToolkit/Numerics
                    /usr/pubsw/packages/itk/current/include/InsightToolkit/Numerics/Statistics
                    /usr/pubsw/packages/itk/current/include/InsightToolkit/Review
                    /usr/pubsw/packages/itk/current/include/InsightToolkit/Review/Statistics
                    /usr/pubsw/packages/itk/current/include/InsightToolkit/SpatialObject
                    /usr/pubsw/packages/itk/current/include/InsightToolkit/Utilities)


# --------------------------------------------------
#                  operating system
# --------------------------------------------------

add_definitions(-D${CMAKE_SYSTEM_NAME})

# NrrdIO-required definitions
add_definitions(-DTEEM_DIO=0 -DTEEM_QNANHIBIT=1)

# check for big endian
include(TestBigEndian)
test_big_endian(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
  add_definitions(-DTEEM_ENDIAN=4321)
else()
  add_definitions(-DTEEM_ENDIAN=1234)
endif()

# non 32-bit
add_definitions(-DTEEM_32BIT=0)


# --------------------------------------------------
#                      packages
# --------------------------------------------------

# designate module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# ZLIB
find_package(ZLIB REQUIRED)
add_definitions(-DHAVE_ZLIB)

# mni
find_package(MNI REQUIRED)

# itk
# find_package(ITK REQUIRED)
add_definitions(-DHAVE_ITK_LIBS)

# OpenGL
# find_package(OpenGL)

# GLUT
find_package(GLUT)
if(GLUT_FOUND)
  add_definitions(-DHAVE_GLUT_LIBS)
  include_directories(${CMAKE_SOURCE_DIR}/glut)
endif()

# vxl
# find_package(VXL REQUIRED)

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  add_definitions(-DHAVE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# SSE matrix and math functions (affine.h, sse_mathfun.h)
add_definitions(-DUSE_SSE_MATHFUN)


# --------------------------------------------------
#                      options
# --------------------------------------------------

# flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3 -msse2 -mfpmath=sse -Wall -Wno-unused-but-set-variable")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -msse2 -mfpmath=sse -Wall -Wno-unused-but-set-variable")

# for stripping code
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdata-sections -ffunction-sections -Wl,--gc-sections")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdata-sections -ffunction-sections -Wl,--gc-sections")

add_definitions(-DANSI)

# create utils c library
add_subdirectory(utils)
