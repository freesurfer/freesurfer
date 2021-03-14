if(FS_BUILD_MATRIX)
  if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
  endif()
  unset(FS_BUILD_MATRIX)
  project(freesurfer_super)
  include(CTest)
  set(CONF_MATRIX
      "coverage;profile;Debug;Release;RelWithDebInfo;MinSizeRel;asan;lsan;ubsan;tsan;cfisan;msan;sssan;thinlto;fulllto"
      )
  set(COMPILER_MATRIX "clang09;clang10;clang11;clang12;clang13;gcc09;gcc10")
  if(APPLE)
    list(APPEND COMPILER_MATRIX "appleclang")
    set(FS_LLVM_PREFIX "/Users/aboualiaa/Desktop/.temp/install")
    set(FS_MATRIX_CACHE "${CMAKE_SOURCE_DIR}/cmake/caches/clion_cache.cmake")
    set(CMAKE_YARD_PREFIX "/Users/aboualiaa/Desktop/.temp/scratch/cmake-build-yard")
    execute_process(COMMAND "sw_vers" "-productVersion" OUTPUT_VARIABLE CURRENT_SYSTEM_VERSION)
    set(CURRENT_SYSTEM_NAME "macOS")
    string(STRIP ${CURRENT_SYSTEM_VERSION} CURRENT_SYSTEM_VERSION)
  else()
    set(FS_LLVM_PREFIX "/media/psf/Home/Desktop/.temp/install/ubuntu")
    set(FS_MATRIX_CACHE "${CMAKE_SOURCE_DIR}/cmake/caches/focal_cache.cmake")
    set(CMAKE_YARD_PREFIX "/media/psf/Home/Desktop/.temp/scratch/cmake-build-yard")
    execute_process(COMMAND "lsb_release" "-is" OUTPUT_VARIABLE CURRENT_SYSTEM_NAME)
    execute_process(COMMAND "lsb_release" "-rs" OUTPUT_VARIABLE CURRENT_SYSTEM_VERSION)
    string(STRIP ${CURRENT_SYSTEM_NAME} CURRENT_SYSTEM_NAME)
    string(STRIP ${CURRENT_SYSTEM_VERSION} CURRENT_SYSTEM_VERSION)
  endif()

  foreach(Z_COMPILER IN LISTS COMPILER_MATRIX)
    foreach(TYPE IN LISTS CONF_MATRIX)

      if(Z_COMPILER STREQUAL "clang09" AND APPLE)
        set(COMPILER_PROXY_OPTION "-DFS_USE_LLVM=${FS_LLVM_PREFIX}/llvm-09/bin")
      elseif(Z_COMPILER STREQUAL "clang09")
        set(COMPILER_PROXY_OPTION "-DCMAKE_CXX_COMPILER=/usr/bin/clang++-9")
        set(COMPILER_C_PROXY_OPTION "-DCMAKE_C_COMPILER=/usr/bin/clang-9")

      elseif(Z_COMPILER STREQUAL "clang10" AND APPLE)
        set(COMPILER_PROXY_OPTION "-DFS_USE_LLVM=${FS_LLVM_PREFIX}/llvm-10/bin")
      elseif(Z_COMPILER STREQUAL "clang10")
        set(COMPILER_PROXY_OPTION "-DCMAKE_CXX_COMPILER=/usr/bin/clang++-10")
        set(COMPILER_C_PROXY_OPTION "-DCMAKE_C_COMPILER=/usr/bin/clang-10")

      elseif(Z_COMPILER STREQUAL "clang11" AND APPLE)
        set(COMPILER_PROXY_OPTION "-DFS_USE_LLVM=${FS_LLVM_PREFIX}/llvm-11/bin")
      elseif(Z_COMPILER STREQUAL "clang11")
        set(COMPILER_PROXY_OPTION "-DCMAKE_CXX_COMPILER=/usr/bin/clang++-11")
        set(COMPILER_C_PROXY_OPTION "-DCMAKE_C_COMPILER=/usr/bin/clang-11")

      elseif(Z_COMPILER STREQUAL "clang12" AND APPLE)
        set(COMPILER_PROXY_OPTION "-DFS_USE_LLVM=${FS_LLVM_PREFIX}/llvm-12/bin")
      elseif(Z_COMPILER STREQUAL "clang12" AND APPLE)
        set(COMPILER_PROXY_OPTION "-DCMAKE_CXX_COMPILER=/usr/bin/clang++-12")
        set(COMPILER_C_PROXY_OPTION "-DCMAKE_C_COMPILER=/usr/bin/clang-12")

      elseif(Z_COMPILER STREQUAL "clang13" AND APPLE)
        set(COMPILER_PROXY_OPTION "-DFS_USE_LLVM=${FS_LLVM_PREFIX}/llvm-13/bin")
      elseif(Z_COMPILER STREQUAL "clang13")
        set(COMPILER_PROXY_OPTION "-DCMAKE_CXX_COMPILER=/usr/bin/clang++-13")
        set(COMPILER_C_PROXY_OPTION "-DCMAKE_C_COMPILER=/usr/bin/clang-13")

      elseif(Z_COMPILER STREQUAL "gcc10" AND APPLE)
        set(COMPILER_PROXY_OPTION "-DFS_USE_GCC=/usr/local/Cellar/gcc/10.2.0_4")
      elseif(Z_COMPILER STREQUAL "gcc10")
        set(COMPILER_PROXY_OPTION "-DCMAKE_CXX_COMPILER=/usr/bin/g++-10")
        set(COMPILER_C_PROXY_OPTION "-DCMAKE_C_COMPILER=/usr/bin/gcc-10")

      elseif(Z_COMPILER STREQUAL "gcc09" AND APPLE)
        set(COMPILER_PROXY_OPTION "-DFS_USE_GCC=/usr/local/Cellar/gcc@9/9.3.0_1")
      elseif(Z_COMPILER STREQUAL "gcc09")
        set(COMPILER_PROXY_OPTION "-DCMAKE_CXX_COMPILER=/usr/bin/g++-9")
        set(COMPILER_C_PROXY_OPTION "-DCMAKE_C_COMPILER=/usr/bin/gcc-9")

      elseif(Z_COMPILER STREQUAL "appleclang")
        set(COMPILER_PROXY_OPTION "")
      endif()

      string(TOLOWER ${TYPE} FOLDER_NAME)

      if(Z_COMPILER MATCHES "clang" AND TYPE STREQUAL "profile")
        message(NOTICE "skipping profile build with llvm clang")
        continue()
      endif()

      if(Z_COMPILER STREQUAL "appleclang" AND TYPE STREQUAL "profile")
        message(NOTICE "skipping profile build with llvm clang")
        continue()
      endif()

      if(Z_COMPILER STREQUAL "gcc" AND TYPE STREQUAL "cfisan")
        message(NOTICE "skipping cfisan build with gcc")
        continue()
      endif()

      if(Z_COMPILER STREQUAL "appleclang" AND TYPE STREQUAL "cfisan")
        message(NOTICE "skipping cfisan build with apple clang")
        continue()
      endif()

      if(NOT LINUX AND TYPE STREQUAL "msan")
        message(NOTICE "skipping msan build with on non linux os")
        continue()
      endif()

      if(NOT LINUX AND TYPE STREQUAL "sssan")
        message(NOTICE "skipping sssan build with on non linux os")
        continue()
      endif()

      if(Z_COMPILER STREQUAL "appleclang" AND TYPE STREQUAL "lsan")
        message(NOTICE "skipping lsan build with apple clang")
        continue()
      endif()

      message(NOTICE "configuring ${TYPE} build with ${Z_COMPILER}")

      add_custom_target(mkdir-${Z_COMPILER}-${TYPE}
                        COMMAND bash -c \" mkdir -p
                                ${CMAKE_YARD_PREFIX}/freesurfer-${CURRENT_SYSTEM_NAME}-${CURRENT_SYSTEM_VERSION}-${FOLDER_NAME}-${Z_COMPILER}
                                \"
                        )

      add_custom_target(cmake-${Z_COMPILER}-${TYPE} ALL
                        COMMAND bash -c \" ${CMAKE_COMMAND} -Wno-dev ${COMPILER_PROXY_OPTION}
                                ${COMPILER_C_PROXY_OPTION} -DCMAKE_BUILD_TYPE=${TYPE} -C ${FS_MATRIX_CACHE} -G Ninja -S
                                ${CMAKE_SOURCE_DIR} -B
                                ${CMAKE_YARD_PREFIX}/freesurfer-${CURRENT_SYSTEM_NAME}-${CURRENT_SYSTEM_VERSION}-${FOLDER_NAME}-${Z_COMPILER}
                                \"
                        )

      add_custom_target(ninja-${Z_COMPILER}-${TYPE}
                        COMMAND bash -c \" ninja -k 0 -l 10 -C
                                ${CMAKE_YARD_PREFIX}/freesurfer-${CURRENT_SYSTEM_NAME}-${CURRENT_SYSTEM_VERSION}-${FOLDER_NAME}-${Z_COMPILER}/
                                \"
                        )

      add_dependencies(cmake-${Z_COMPILER}-${TYPE} mkdir-${Z_COMPILER}-${TYPE})

      if(NOT TARGET compile-${Z_COMPILER})
        add_custom_target(compile-${Z_COMPILER})
      endif()

      if(NOT TARGET compile-${Z_COMPILER}-${TYPE})
        add_custom_target(compile-${Z_COMPILER}-${TYPE})
      endif()

      add_dependencies(compile-${Z_COMPILER} ninja-${Z_COMPILER}-${TYPE})
      add_dependencies(compile-${Z_COMPILER}-${TYPE} ninja-${Z_COMPILER}-${TYPE})

      if(NOT TARGET compile-${TYPE})
        add_custom_target(compile-${TYPE})
      endif()

      add_dependencies(compile-${TYPE} ninja-${Z_COMPILER}-${TYPE})

      if(NOT TARGET compile-all)
        add_custom_target(compile-all)
      endif()

      add_dependencies(compile-all compile-${TYPE})

      add_custom_target(ctest-${Z_COMPILER}-${TYPE}
                        COMMAND bash -c \" ${CMAKE_CTEST_COMMAND} -T Test -L Unit --progress --test-load 10 \"
                        WORKING_DIRECTORY "${CMAKE_YARD_PREFIX}/freesurfer-${CURRENT_SYSTEM_NAME}-${CURRENT_SYSTEM_VERSION}-${FOLDER_NAME}-${Z_COMPILER}"
                        )

      add_test(ctest-${Z_COMPILER}-${TYPE}
               COMMAND ${CMAKE_CTEST_COMMAND} -T Test -L Unit --progress --test-load 10
               WORKING_DIRECTORY "${CMAKE_YARD_PREFIX}/freesurfer-${CURRENT_SYSTEM_NAME}-${CURRENT_SYSTEM_VERSION}-${FOLDER_NAME}-${Z_COMPILER}"
               )

      if(NOT TARGET check-${Z_COMPILER})
        add_custom_target(check-${Z_COMPILER})
      endif()

      add_dependencies(ctest-${Z_COMPILER}-${TYPE} ninja-${Z_COMPILER}-${TYPE})
      add_dependencies(check-${Z_COMPILER} ctest-${Z_COMPILER}-${TYPE})

      if(NOT TARGET check-${TYPE})
        add_custom_target(check-${TYPE})
      endif()

      add_dependencies(check-${TYPE} ctest-${Z_COMPILER}-${TYPE})

      if(NOT TARGET check-all)
        add_custom_target(check-all)
      endif()

      add_dependencies(check-all check-${TYPE})

      if(NOT TARGET freshstart)
        add_custom_target(freshstart)
      endif()

      if(NOT TARGET freshstart-${Z_COMPILER})
        add_custom_target(freshstart-${Z_COMPILER})
      endif()

      if(NOT TARGET freshstart-${TYPE})
        add_custom_target(freshstart-${TYPE})
      endif()

      add_custom_target(freshstart-${Z_COMPILER}-${TYPE}
                        COMMAND bash -c \" rm -rf
                                ${CMAKE_YARD_PREFIX}/freesurfer-${CURRENT_SYSTEM_NAME}-${CURRENT_SYSTEM_VERSION}-${TYPE}-${Z_COMPILER}
                                \"
                        )

      add_dependencies(freshstart-${Z_COMPILER} freshstart-${Z_COMPILER}-${TYPE})
      add_dependencies(freshstart-${TYPE} freshstart-${Z_COMPILER}-${TYPE})
      add_dependencies(freshstart freshstart-${Z_COMPILER} freshstart-${TYPE})

      if(NOT TARGET clean-${Z_COMPILER})
        add_custom_target(clean-${Z_COMPILER})
      endif()

      if(NOT TARGET clean-${TYPE})
        add_custom_target(clean-${TYPE})
      endif()

      if(NOT TARGET clean-all)
        add_custom_target(clean-all)
      endif()

      add_custom_target(clean-${Z_COMPILER}-${TYPE}
                        COMMAND bash -c \" ninja -C
                                ${CMAKE_YARD_PREFIX}/freesurfer-${CURRENT_SYSTEM_NAME}-${CURRENT_SYSTEM_VERSION}-${FOLDER_NAME}-${Z_COMPILER}
                                clean || true \"
                        )

      add_dependencies(clean-${Z_COMPILER} clean-${Z_COMPILER}-${TYPE})
      add_dependencies(clean-${TYPE} clean-${Z_COMPILER}-${TYPE})
      add_dependencies(clean-all clean-${Z_COMPILER}-${TYPE})
    endforeach()
  endforeach()
  return()
endif()
if(UNIX AND NOT APPLE)
  set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran")
elseif(APPLE)
  set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran")
endif()
if(DEFINED FS_USE_GCC AND NOT FS_USE_GCC STREQUAL "")
  if(FS_USE_GCC MATCHES "10")
    set(CMAKE_CXX_COMPILER "${FS_USE_GCC}/bin/g++-10")
    set(CMAKE_C_COMPILER "${FS_USE_GCC}/bin/gcc-10")
    set(CMAKE_AR "${FS_USE_GCC}/bin/gcc-ar-10")
    set(CMAKE_NM "${FS_USE_GCC}/bin/gcc-nm-10")
    set(CMAKE_RANLIB "${FS_USE_GCC}/bin/gcc-ranlib-10")
  else()
    set(CMAKE_CXX_COMPILER "${FS_USE_GCC}/bin/g++-9")
    set(CMAKE_C_COMPILER "${FS_USE_GCC}/bin/gcc-9")
    set(CMAKE_AR "${FS_USE_GCC}/bin/gcc-ar-9")
    set(CMAKE_NM "${FS_USE_GCC}/bin/gcc-nm-9")
    set(CMAKE_RANLIB "${FS_USE_GCC}/bin/gcc-ranlib-9")
  endif()
  #    add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
  #    set(ITK_DIR "/Users/aboualiaa/Downloads/packages/itk/4.13.0")
  #    set(VTK_DIR "/Users/aboualiaa/Downloads/packages/vtk/5.10.1")
  #    add_definitions(-DFS_ITK_LEGACY_REMOVE)
endif()
if(DEFINED FS_USE_LLVM
   AND NOT FS_USE_LLVM STREQUAL ""
   AND APPLE
   )
  unset(CMAKE_C_COMPILER CACHE)
  unset(CMAKE_CXX_COMPILER CACHE)
  unset(CMAKE_AR CACHE)
  unset(CMAKE_NM CACHE)
  unset(CMAKE_RANLIB CACHE)
  set(CMAKE_C_COMPILER "${FS_USE_LLVM}/clang")
  set(CMAKE_CXX_COMPILER "${FS_USE_LLVM}/clang++")
  set(CMAKE_AR "${FS_USE_LLVM}/llvm-ar")
  set(CMAKE_NM "${FS_USE_LLVM}/llvm-nm")
  set(CMAKE_RANLIB "${FS_USE_LLVM}/llvm-ranlib")
endif()

include(cmake/caches/prepare_build_types.cmake)

if(POLICY CMP0025) # AppleClang and Clang are different compilers.
  cmake_policy(SET CMP0025 NEW)
endif()

if(POLICY CMP0048) # The project() command manages VERSION variables.
  cmake_policy(SET CMP0048 NEW)
endif()

if(NOT DEFINED FS_VERSION_MAJOR)
  set(FS_VERSION_MAJOR 9)
endif()
if(NOT DEFINED FS_VERSION_MINOR)
  set(FS_VERSION_MINOR 0)
endif()
if(NOT DEFINED FS_VERSION_PATCH)
  set(FS_VERSION_PATCH 0)
endif()
if(NOT DEFINED FS_VERSION_SUFFIX)
  set(FS_VERSION_SUFFIX "")
endif()

if(NOT PACKAGE_VERSION)
  set(PACKAGE_VERSION "${FS_VERSION_MAJOR}.${FS_VERSION_MINOR}.${FS_VERSION_PATCH}${FS_VERSION_SUFFIX}")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "[Cc]lang"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "10"
   AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
   )
  set(CMAKE_RANLIB "${FS_USE_LLVM}/llvm-ranlib")
endif()

if(APPLE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)[Cc]lang")
    set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
    set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
    set(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
    set(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
  endif()
  mark_as_advanced(CMAKE_OSX_SYSROOT CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_ARCHITECTURES CMAKE_EXECUTABLE_FORMAT)
endif()

find_program(CCACHE_PROGRAM ccache)
mark_as_advanced(CCACHE_PROGRAM)
if(CCACHE_PROGRAM AND FS_USE_CCACHE #   AND NOT CMAKE_BUILD_TYPE STREQUAL "coverage"
   )
  message(STATUS "Using ccache at ${CCACHE_PROGRAM}")
  set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
endif()

function(prevent_in_source_builds_and_installs)
  get_filename_component(SRCDIR "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(BINDIR "${CMAKE_BINARY_DIR}" REALPATH)

  # disallow in-source builds
  if("${SRCDIR}" STREQUAL "${BINDIR}")
    message(FATAL_ERROR "In-source builds/installs are disallowed. source dir: ${SRCDIR}; build dir: ${BINDIR}")
  endif()
endfunction()

# project info
cmake_minimum_required(VERSION 3.9)

# requiring at least c++11 standard
if(CMAKE_CXX_STANDARD EQUAL "98")
  message(FATAL_ERROR "CMAKE_CXX_STANDARD:STRING=98 is not supported in.")
endif()

# set the default language dialect settings (This must be done BEFORE the project call)
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
endif()
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(NOT CMAKE_CXX_EXTENSIONS)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

project(freesurfer
        VERSION ${PACKAGE_VERSION}
        DESCRIPTION "An open source software suite for processing and analyzing (human) brain MRI images."
        LANGUAGES C CXX Fortran
        )

if(NOT APPLE)
  cmake_minimum_required(VERSION 3.9)
  set(HOST_OS $ENV{HOSTING_OS})
  # TODO: upgrade to itk 5 then remove
  if(HOST_OS STREQUAL "trusty")
    add_definitions(-DFS_ITK_LEGACY_REMOVE)
  endif()
  add_definitions(-DVCL_GCC_80) # trick itk 4 into believing the compiler is gcc 8
  add_definitions(-DVXL_FULLCXX11SUPPORT)
else()
  cmake_minimum_required(VERSION 3.10.2) # preparing for itk5
endif()

if(FS_SUPRESS_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wno-everything)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-w)
  endif()
  add_link_options(-w)
endif()

if(CMAKE_CONFIGURATION_TYPES AND FS_EXP_BUILD_CONFIGURATIONS)
  if(NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    list(REMOVE_ITEM FS_EXP_BUILD_CONFIGURATIONS "profile")
  endif()
  if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    list(REMOVE_ITEM FS_EXP_BUILD_CONFIGURATIONS "cfisan")
  endif()
  if(NOT LINUX)
    list(REMOVE_ITEM FS_EXP_BUILD_CONFIGURATIONS "msan")
  endif()
  if("${CMAKE_C_COMPILER_ID}" MATCHES "(Apple)[Cc]lang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "(Apple)[Cc]lang")
    list(REMOVE_ITEM FS_EXP_BUILD_CONFIGURATIONS "lsan")
  endif()
  set(CMAKE_CONFIGURATION_TYPES ${FS_EXP_BUILD_CONFIGURATIONS})
endif()

if(FS_BUILD_TESTING)
  include(CTest)
endif()

include(cmake/functions.cmake)
include(cmake/HandleFSOptions.cmake)
include(TestBigEndian)
include(CheckFunctionExists)

# prevent third-party packages from importing as a system
set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)

# if an install prefix is not provided, check the FS_INSTALL_DIR env var, and
# if that is not defined, set the default path to /usr/local/freesurfer
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(NOT " $ENV{FS_INSTALL_DIR}" STREQUAL "")
    set(CMAKE_INSTALL_PREFIX
        "$ENV{FS_INSTALL_DIR}"
        CACHE PATH "Copied from FS_INSTALL_DIR env variable" FORCE
        )
  else()
    set(CMAKE_INSTALL_PREFIX
        "/usr/local/freesurfer"
        CACHE PATH "Default install path" FORCE
        )
  endif()
endif()

# version stamp
set(FS_VERSION
    "$ENV{USER}-local"
    CACHE STRING "Distribution version"
    )

# build stamp
string(TIMESTAMP TODAY "%Y%m%d")
set(BUILD_STAMP
    "freesurfer-local-build-${TODAY}"
    CACHE STRING "Distribution build stamp"
    )
install(CODE "file (WRITE ${CMAKE_INSTALL_PREFIX}/build-stamp.txt ${BUILD_STAMP}\\n) ")

# set the default build type to 'Release' for optimization purposes
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS " No build type selected - defaulting to Release")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Default build type" FORCE
      )
endif()

# create a custom target "make install-nmr" to only install nmr-component files
add_custom_target(install-nmr COMMAND ${CMAKE_COMMAND} -D COMPONENT=nmr -P cmake_install.cmake)

# --------------------------------------------------
#               external dependencies
# --------------------------------------------------

set(ABSL_RUN_TESTS OFF)
if(BUILD_TESTING)
  set(BUILD_TESTING OFF)
endif()
add_subdirectory(third_party/abseil-cpp)
set(BUILD_TESTING ${_BUILD_TESTING})

# xxd is used to generate the helptext headers
find_program(XXD xxd)
if(NOT XXD)
  message(WARNING "MISSING xxd - required to build executable help texts")
endif()
mark_as_advanced(XXD)

# Most of the packages required by freesurfer are located by custom find-modules stored in the
# 'cmake' subdir. The find-modules expect each package to be installed under a common
# path defined by FS_PACKAGES_DIR. On Martinos machines, this variable automatically defaults
# to /usr/pubsw/packages, but external developers must provide this path manually. External developers
# can run the packages/build_packages.py script to compile the dependencies locally. If a package
# is not found under FS_PACKAGES_DIR, cmake will continue to look through the default search paths.
# Additionally, alternative paths to package installs can be specified with the <PACKAGE>_DIR variables
# Note: FS_PACKAGES_DIR can be defined in an environment variable

if(NOT FS_PACKAGES_DIR)
  if(NOT "$ENV{FS_PACKAGES_DIR}" STREQUAL "")
    # check if FS_PACKAGES_DIR has been set
    set(FS_PACKAGES_DIR
        "$ENV{FS_PACKAGES_DIR}"
        CACHE INTERNAL "Copied from FS_PACKAGES_DIR environment variable"
        )
  elseif(EXISTS /usr/pubsw/packages)
    # if the user is at Martinos, default to /usr/pubsw/packages
    set(FS_PACKAGES_DIR
        /usr/pubsw/packages
        CACHE INTERNAL "Default Martinos packages dir"
        )
  else()
    # setting FS_PACKAGES_DIR is required if the user is outside Martinos
    # for more information visit https://surfer.nmr.mgh.harvard.edu/fswiki/BuildRequirements
    message(WARNING "FS_PACKAGES_DIR is NOT defined - build will be limited")
  endif()
endif()

# all custom find-modules are stored in the cmake subdir
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(FS_USE_ARRAYFIRE)
  find_package(ArrayFire)
  mark_as_advanced(ArrayFire_DIR)
endif()
find_package(MPI REQUIRED)

# -------- zlib --------
find_package(ZLIB REQUIRED)
add_definitions(-DHAVE_ZLIB)

# -------- gfortran/blas/lapack --------
if(NOT APPLE)
  find_library(GFORTRAN_LIBRARIES
               HINTS /usr/lib64 /usr/lib /usr/lib/gcc/*/*/
               NAMES libgfortran.a gfortran
               )
  find_library(BLAS_LIBRARIES NAMES libblas.a libopenblas.a)
  find_library(LAPACK_LIBRARIES NAMES liblapack.a)
  find_library(QUADMATH_LIBRARIES
               HINTS /usr/lib64 /usr/lib /usr/lib/gcc/*/*/
               NAMES libquadmath.a libquadmath.so libquadmath.so.0
               )
  if(NOT QUADMATH_LIBRARIES)
    set(QUADMATH_LIBRARIES "")
  endif()
endif()

# -------- armadillo --------
find_package(ARMADILLO)
mark_as_advanced(ARMADILLO_LIBRARIES ARMADILLO_LIBRARIES ARMADILLO_INCLUDE_DIRS)

# -------- petsc --------
find_package(PETSC REQUIRED)

# --- boost ---
set(Boost_NO_BOOST_CMAKE ON)
set(Boost_USE_MULTITHREADED TRUE)
find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework filesystem log atomic)

# -------- itk --------
find_package(ITK HINTS ${ITK_DIR} REQUIRED)
mark_as_advanced(ITK_DIR)

if(APPLE)
  list(REMOVE_ITEM ITK_INCLUDE_DIRS "/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/usr/include")
endif()

# -------- ann --------
# ANN is only used by mris_resample
set(ANN_DIR "${FS_PACKAGES_DIR}/ann/1.1.2/")
find_package(ANN)

# -------- gts --------
find_package(GTS)
mark_as_advanced(ICONV_LIBRARY)

#endif()

# -------- matlab --------
if(FS_BUILD_MATLAB)
  if(UNIX AND NOT APPLE)
    set(Matlab_ROOT_DIR /usr/local/MATLAB/MATLAB_Compiler_Runtime/v84)
    find_package(Matlab)
  else()
    set(Matlab_ROOT_DIR /Applications/Polyspace/MATLAB.app)
    find_package(Matlab)
  endif()
  mark_as_advanced(Matlab_ROOT_DIR)
endif()
# -------- hpx --------
if(FS_USE_HPX)
  set(HPX_IGNORE_COMPILER_COMPATIBILITY ON)
  set(HPX_IGNORE_CMAKE_BUILD_TYPE_COMPATIBILITY OFF FORCE)
  find_package(HPX HINTS "/Users/aboualiaa/Desktop/.temp/install/hpx/")
  mark_as_advanced(HPX_GIT_COMMIT
                   HPX_Fortran_COMPILER_VERSION
                   HPX_Fortran_COMPILER_ID
                   HPX_Fortran_COMPILER
                   HPX_DIR
                   HPX_C_COMPILER_VERSION
                   HPX_C_COMPILER_ID
                   HPX_C_COMPILER
                   HPX_CXX_COMPILER_VERSION
                   HPX_CXX_COMPILER_ID
                   HPX_CXX_COMPILER
                   LLVM_FILECHECK_EXE
                   )
endif()

# -------- openCV --------
if(NOT APPLE)
  find_package(OpenCV)
endif()

# -------- X11 ---------
find_package(X11)
mark_as_advanced(XT_LIBRARY XMU_LIBRARY X11_LIBRARY X11_INCLUDE_DIR GL_LIBRARY GLU_LIBRARY)

# search for GUI libraries
if(FS_BUILD_GUIS)

  # -------- openGL --------
  if(EXISTS /usr/lib64/nvidia/libGL.so)
    # on some martinos centos7 machines, linking to the default libGL.so fails due to missing symbols
    # in the nvidia drivers, so as a temporary fix, link directly to the nvidia GL library
    list(APPEND CMAKE_PREFIX_PATH "/usr/lib64/nvidia")
  endif()
  find_package(OpenGL)
  if(OPENGL_FOUND)
    add_definitions(-DHAVE_OPENGL)
  endif()

  # -------- qt --------
  find_package(QT)

  mark_as_advanced(Qt5Widgets_DIR Qt5Sql_DIR Qt5Gui_DIR Qt5Core_DIR)
endif()

# -------- vtk --------
#set(VTK_DIR ${FS_PACKAGES_DIR}/vtk/7.1 CACHE PATH "VTK install directory")
find_package(VTK HINTS ${VTK_DIR})
if(VTK_FOUND AND VTK_BUILD_SHARED_LIBS)
  # locate all shared vtk libraries
  if(${VTK_VERSION} LESS 7)
    file(GLOB _vtk_library_list "${VTK_LIBRARY_DIRS}/lib*.so*")
  else()
    foreach(_vtk_module ${VTK_MODULES_ENABLED})
      foreach(_vtk_library ${${_vtk_module}_LIBRARIES})
        file(GLOB _found_library "${${_vtk_module}_RUNTIME_LIBRARY_DIRS}/lib${_vtk_library}*.so*")
        list(APPEND _vtk_library_list ${_found_library})
      endforeach()
      foreach(_vtk_library ${${_vtk_module}_DEPENDS})
        file(GLOB _found_library "${${_vtk_module}_RUNTIME_LIBRARY_DIRS}/lib${_vtk_library}*.so*")
        list(APPEND _vtk_library_list ${_found_library})
      endforeach()
    endforeach()
  endif()
  # install found libraries to the freesurfer lib/vtk dir and add this directory to rpath
  if(_vtk_library_list)
    list(REMOVE_DUPLICATES _vtk_library_list)
    install(PROGRAMS ${_vtk_library_list} DESTINATION lib/vtk)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/vtk:${CMAKE_INSTALL_RPATH}")
  endif()
endif()
list(REMOVE_ITEM VTK_LIBRARIES "/usr/lib/libxml2.dylib")

# --------------------------------------------------
#              global system information
# --------------------------------------------------

add_definitions(-D${CMAKE_SYSTEM_NAME})

# general check for big endian
test_big_endian(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
  set(BYTEORDER 4321)
else()
  set(BYTEORDER 1234)
endif()

# --------------------------------------------------
#                  third-party code
# --------------------------------------------------
set(BENCHMARK_ENABLE_TESTING
    OFF
    CACHE BOOL "Suppressing benchmark's tests" FORCE
    )

if(APPLE)
  find_package(EXPAT REQUIRED)
  find_package(LibXml2 REQUIRED)
endif()
#find_package(tiff REQUIRED)
if(UNIX AND NOT APPLE)
  set(EXPAT_LIBRARIES "/usr/lib/x86_64-linux-gnu/libexpat.so")
  set(LibXml2_DIR "/usr/lib/x86_64-linux-gnu")
  find_package(LibXml2 REQUIRED)
endif()
find_package(fmt REQUIRED HINTS "/usr/local/opt/fmt")
find_package(spdlog REQUIRED HINTS "/usr/local/opt/spdlog")
#find_package(Microsoft.GSL REQUIRED HINTS "/usr/local/opt/cpp-gsl")
#find_package(gtest REQUIRED)
#find_package(benchmark REQUIRED)
add_subdirectory(packages)
add_subdirectory(third_party)

#if(APPLE)
#  # this is an ugly hack to ensure that cmake looks for the
#  # static version of libomp on mac, since we don't want
#  # to distribute the shared library
#  set(ORIG_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
#  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#  find_package(OpenMP)
#  set(CMAKE_FIND_LIBRARY_SUFFIXES ${ORIG_LIBRARY_SUFFIXES})
#else()
#  find_package(OpenMP)
#endif()
#
#if(OPENMP_FOUND)
#  add_definitions(-DHAVE_OPENMP)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#  set(OMP_CXX_LIBRARIES OpenMP::OpenMP_CXX)
#endif()

# -------- openMP --------
if(APPLE)
  if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(OpenMP_C "${CMAKE_C_COMPILER}")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -lomp -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include")
    set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
    set(OpenMP_CXX_FLAGS
        "-Xpreprocessor -fopenmp -lomp -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include -L/usr/local/opt/hdf5/lib -L/usr/local/opt/szip/lib"
        )
    set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  endif()
  add_link_options(-Xpreprocessor
                   -fopenmp
                   -lomp
                   -L/usr/local/opt/libomp/lib
                   -I/usr/local/opt/libomp/include
                   -L/usr/local/opt/hdf5/lib
                   -L/usr/local/opt/szip/lib
                   )
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
  add_definitions(-DHAVE_OPENMP -D_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(OMP_CXX_LIBRARIES OpenMP::OpenMP_CXX)
endif()

# --------------------------------------------------
#                    setup python
# --------------------------------------------------

if(FS_DISTRIBUTE_FSPYTHON)
  set(FSPYTHON_DIR ${FS_PACKAGES_DIR}/fspython/3.6)
  if(EXISTS ${FSPYTHON_DIR}/bin/python3)
    # specify the python to use
    set(PYTHON_EXECUTABLE ${FSPYTHON_DIR}/bin/python3)
    set(PYTHON_LIBRARY ${FSPYTHON_DIR}/lib/libpython3.6m)
    # install the barebones python distribution
    install(DIRECTORY ${FSPYTHON_DIR}/bin ${FSPYTHON_DIR}/lib ${FSPYTHON_DIR}/include
            DESTINATION ${CMAKE_INSTALL_PREFIX}/python
            USE_SOURCE_PERMISSIONS
            )
  else()
    message(FATAL_ERROR "Cannot find fspython distribution in FS_PACKAGES_DIR")
  endif()
endif()

# initialize pybind for python wrapping
set(PYBIND11_PYTHON_VERSION 3.5)
add_subdirectory(packages/pybind11)
mark_as_advanced(PYBIND11_TEST PYBIND11_INSTALL USE_PYTHON_INCLUDE_DIR)

if(NOT FS_DISTRIBUTE_FSPYTHON)
  # link to an external python binary if fspython won't be distributed
  symlink(${PYTHON_EXECUTABLE} ${CMAKE_INSTALL_PREFIX}/python/bin/python3)
endif()

# --------------------------------------------------
#             freesurfer build settings
# --------------------------------------------------

# linker options
if(APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -dead_strip")
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
  endif()
else()
  set(STRIP_FLAGS "-fdata-sections -ffunction-sections -Wl,--gc-sections")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STRIP_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${STRIP_FLAGS} -Wl,-Map,ld_map.txt -Wl,--no-demangle")
endif()

# --------------------------------------------------
#                 build freesurfer
# --------------------------------------------------

# prepare the freesurfer distribution
add_subdirectory(distribution)

# the top-level include dir contain the most commonly included freesurfer header files
set(FS_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/include
    ${CMAKE_SOURCE_DIR}/packages/minc
    ${CMAKE_SOURCE_DIR}/packages/netcdf
    ${CMAKE_SOURCE_DIR}/packages/nifti
    ${CMAKE_SOURCE_DIR}/include/mri_tools
    ${CMAKE_SOURCE_DIR}/packages/benchmark/include/
    ${CMAKE_SOURCE_DIR}/packages/fmt/include/
    ${CMAKE_SOURCE_DIR}/packages/GSL/include
    ${CMAKE_SOURCE_DIR}/packages/spdlog/include
    ${CMAKE_SOURCE_DIR}/packages/tiff/
    ${CMAKE_SOURCE_DIR}/third_party/abseil-cpp/
    )

if(NOT APPLE)
  set(FS_INCLUDE_DIRS ${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/xml2/)
endif()
# build the freesurfer static libraries
add_subdirectory(utils)

# the following utility libraries are required for freeview
add_subdirectories(vtkutils lineprof)

# TODO: solve build issues on linux/gcc: looser throw spec...
# the fem elastic subdir also builds a library, so it should be added first
if("${CMAKE_C_COMPILER_ID}" MATCHES "(Apple)?[Cc]lang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "(Apple)?[Cc]lang")
  if(APPLE)
    add_subdirectory(fem_elastic)
  endif()
endif()

# build the gems library
add_subdirectory(gems)

# --------------------------------------------------
#                      programs
# --------------------------------------------------

# the following subdirectories contain programs required for recon-all
# and should be included during a "minimal " install. Extra programs that aren't
# used in the standard stream should be added further down
add_subdirectories(AntsDenoiseImageFs
                   AntsN4BiasFieldCorrectionFs
                   lta_convert
                   mri_and
                   mri_add_xform_to_header
                   mri_annotation2label
                   mri_aparc2aseg
                   mri_aparc2wmseg
                   mri_binarize
                   mri_brainvol_stats
                   mri_ca_label
                   mri_ca_normalize
                   mri_ca_register
                   mri_cc
                   mri_compute_overlap
                   mri_compute_seg_overlap
                   mri_concat
                   mri_concatenate_lta
                   mri_convert
                   mri_coreg
                   mri_deface
                   mri_diff
                   mri_edit_segmentation
                   mri_edit_segmentation_with_surfaces
                   mri_edit_wm_with_aseg
                   mri_em_register
                   mri_extract_largest_CC
                   mri_fill
                   mri_fuse_segmentations
                   mri_fwhm
                   mri_gcut
                   mri_info
                   mri_label2label
                   mri_label2vol
                   mri_log_likelihood
                   mri_make_bem_surfaces
                   mri_mask
                   mri_matrix_multiply
                   mri_mc
                   mri_normalize
                   mri_normalize_tp2
                   mri_probedicom
                   mri_refine_seg
                   mri_relabel_hypointensities
                   mri_relabel_nonwm_hypos
                   mri_remove_neck
                   mri_robust_register
                   mri_seg_overlap
                   mri_seg_diff
                   mri_segment
                   mri_segreg
                   mri_segstats
                   mri_stats2seg
                   mri_surf2surf
                   mri_surf2vol
                   mri_surfcluster
                   mri_tessellate
                   mri_vol2surf
                   mri_vol2vol
                   mri_voldiff
                   mri_watershed
                   mris_anatomical_stats
                   mris_annot_diff
                   mris_ca_label
                   mris_calc
                   mris_convert
                   mris_curvature
                   mris_curvature_stats
                   mris_defects_pointset
                   mris_diff
                   mris_divide_parcellation
                   mris_euler_number
                   mris_fix_topology
                   mris_inflate
                   mris_info
                   mris_jacobian
                   mris_label2annot
                   mris_left_right_register
                   mris_make_surfaces
                   mris_register
                   mris_smooth
                   mris_sphere
                   mris_surface_stats
                   mris_thickness
                   mris_thickness_diff
                   mris_topo_fixer
                   mris_volmask
                   mrisp_paint
                   python
                   samseg
                   scripts
                   talairach_afd
                   talairach_avi
                   tkregister2
                   )

# the following program subdirectories aren't required in the standard recon-all stream.
# they will be built by default, but not if a minimal build is configured
if(NOT FS_MINIMAL_BUILD)
  add_subdirectories(anatomicuts
                     attic
                     BrainstemSS
                     diffusion_tool
                     dummy
                     freeview
                     fsfast
                     fslutils
                     hiam_make_surfaces
                     hiam_make_template
                     hiam_register
                     HippoSF
                     histo_compute_joint_density
                     histo_register_block
                     histo_synthesize
                     label2flat
                     label2patch
                     mkxsubjreg
                     matlab
                     mri_average
                     mri_bias
                     mri_ca_tissue_parms
                     mri_ca_train
                     mri_cal_renormalize_gca
                     mri_cnr
                     mri_compile_edits
                     mri_compute_change_map
                     mri_compute_volume_fractions
                     mri_concatenate_gcam
                     mri_copy_values
                     mri_cor2label
                     mri_correct_segmentations
                     mri_cvs_register
                     mri_dct_align
                     mri_distance_transform
                     mri_dist_surf_label
                     mri_evaluate_morph
                     mri_extract
                     mri_extract_fcd_features
                     mri_extract_label
                     mri_fdr
                     mri_fieldsign
                     mri_fit_bias
                     mri_fslmat_to_lta
                     mri_fuse_intensity_images
                     mri_gca_ambiguous
                     mri_glmfit
                     mri_gtmpvc
                     mri_gtmseg
                     mri_hausdorff_dist
                     mri_head
                     mri_hires_register
                     mri_histo_eq
                     mri_jacobian
                     mri_joint_density
                     mri_label_histo
                     mri_label_vals
                     mri_label_volume
                     mri_linear_register
                     mri_map_cpdat
                     mri_mark_temporal_lobe
                     mri_mcsim
                     mri_mi
                     mri_modify
                     mri_morphology
                     mri_ms_fitparms
                     mri_nlfilter
                     mri_otl
                     mri_paint
                     mri_parse_sdcmdir
                     mri_path2label
                     mri_polv
                     mri_probe_ima
                     mri_reduce
                     mri_rf_label
                     mri_rf_long_train
                     mri_rf_train
                     mri_ribbon
                     mri_rigid_register
                     mri_sbbr
                     mri_segcentroids
                     mri_seghead
                     mri_segment_hypothalamic_subunits
                     mri_strip_nonwhite
                     mri_strip_subject_info
                     mri_surfacemask
                     mri_synthesize
                     mri_threshold
                     mri_topologycorrection
                     mri_train
                     mri_transform
                     mri_twoclass
                     mri_volcluster
                     mri_volsynth
                     mri_warp_convert
                     mri_wbc
                     mri_xvolavg
                     mri_z2p
                     mris2rgb
                     mris_annot_to_segmentation
                     mris_apply_reg
                     mris_average_curvature
                     mris_BA_segment
                     mris_ca_train
                     mris_compute_acorr
                     mris_compute_lgi
                     mris_compute_overlap
                     mris_compute_parc_overlap
                     mris_compute_volume_fractions
                     mris_congeal
                     mris_copy_header
                     mris_deform
                     mris_distance_map
                     mris_distance_to_label
                     mris_distance_transform
                     mris_entropy
                     mris_errors
                     mris_expand
                     mris_extract_patches
                     mris_fill
                     mris_find_flat_regions
                     mris_flatten
                     mris_fwhm
                     mris_hausdorff_dist
                     mris_init_global_tractography
                     mris_interpolate_warp
                     mris_label_area
                     mris_label_calc
                     mris_label_mode
                     mris_make_average_surface
                     mris_make_face_parcellation
                     mris_make_template
                     mris_map_cuts
                     mris_merge_parcellations
                     mris_mesh_subdivide
                     mris_morph_stats
                     mris_ms_refine
                     mris_multiscale_stats
                     mris_niters2fwhm
                     mris_parcellate_connectivity
                     mris_pmake
                     mris_register_label_map
                     mris_register_to_volume
                     mris_remove_variance
                     mris_remesh
                     mris_resample
                     mris_rescale
                     mris_reverse
                     mris_rf_label
                     mris_rf_train
                     mris_rotate
                     mris_sample_label
                     mris_sample_parc
                     mris_seg2annot
                     mris_segment
                     mris_segment_vals
                     mris_segmentation_stats
                     mris_shrinkwrap
                     mris_simulate_atrophy
                     mris_smooth_intracortical
                     mris_spherical_average
                     mris_surface_to_vol_distances
                     mris_talairach
                     mris_thickness_comparison
                     mris_transform
                     mris_translate_annotation
                     mris_transmantle_dysplasia_paths
                     mris_volume
                     mris_warp
                     mris_watershed
                     mrisp_write
                     nmovie_qt
                     oct_register_mosaic
                     optseq2
                     qdecproject
                     qdec_glmfit
                     resurf
                     spline3
                     stat_normalize
                     stim_polar
                     swi_processing
                     template
                     test_makevol
                     ThalamicNuclei
                     trc
                     tridec
                     )
endif()

host_os()

message(STATUS "Building on host OS ${HOST_OS}")

if(FS_INFANT_MODULE)
  add_subdirectory(infant)
  add_subdirectory(sscnn_skullstripping)
endif()

if(FS_TKTOOLS_MODULE)
  add_subdirectory(tktools)
endif()

if(FS_QATOOLS_MODULE)
  add_subdirectory(qatools)
endif()

if(FS_BUILD_DNG)
  add_subdirectory(dngtester)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "coverage")

  add_custom_target(check-all
                    COMMAND bash -c \" ctest -T Test -j 12 --schedule-random -L Unit || true \"
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
                    )

  add_custom_target(generate-coverage-reports)
  get_directory_property(CURR_ADDIDITIONAL_CLEAN_FILES DIRECTORY ${CMAKE_BINARY_DIR} ADDITIONAL_CLEAN_FILES)
  set(GCOV_ARTIFACTS "coverage.info;gcov-report")
  set(LLVM_COV_ARTIFACTS "profiles.manifest;profiles;merged.profdata;llvm-cov-report")

  if(FS_COVERAGE_STYLE STREQUAL "gcov" OR FS_COVERAGE_STYLE STREQUAL "all")
    add_custom_target(generate-gcov-coverage-report
                      COMMAND bash -c \" lcov --capture --directory ${CMAKE_BINARY_DIR} --output-file
                              ${CMAKE_BINARY_DIR}/coverage.info && lcov --remove ${CMAKE_BINARY_DIR}/coverage.info -o
                              ${CMAKE_BINARY_DIR}/coverage.info '${CMAKE_SOURCE_DIR}/packages/*' '/usr/local/*'
                              '${CMAKE_SOURCE_DIR}/third_party/*' '/Applications/Xcode.app/*'
                              '/Library/Developer/CommandLineTools/*' && genhtml ${CMAKE_BINARY_DIR}/coverage.info
                              --output-directory ${CMAKE_BINARY_DIR}/gcov-report \"
                      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
                      )

    add_custom_target(upload-gcov-coverage-report
                      COMMAND bash -c \" ${CMAKE_SOURCE_DIR}/repo_tools/codecov.sh -f ${CMAKE_BINARY_DIR}/coverage.info
                              -t ${FS_CODECOV_TOKEN} \"
                      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                      )

    add_dependencies(generate-gcov-coverage-report check-all)
    add_dependencies(upload-gcov-coverage-report generate-gcov-coverage-report)
    add_dependencies(generate-coverage-reports generate-gcov-coverage-report)
    list(APPEND CURR_ADDIDITIONAL_CLEAN_FILES "GCOV_ARTIFACTS")

  endif(coverage_style_gcov)

  if(FS_COVERAGE_STYLE STREQUAL "llvm" OR FS_COVERAGE_STYLE STREQUAL "all")
    get_filename_component(CMAKE_CXX_COMPILER_PATH ${CMAKE_CXX_COMPILER} DIRECTORY)
    set(FS_LLVM_COV_PATH "${CMAKE_CXX_COMPILER_PATH}/llvm-cov")
    if(EXISTS "${FS_LLVM_COV_PATH}")
      add_custom_target(generate-llvm-coverage-report
                        COMMAND bash -c \" /bin/ls -1 -d ${CMAKE_BINARY_DIR}/profiles/* >
                                ${CMAKE_BINARY_DIR}/profiles.manifest && ${CMAKE_CXX_COMPILER_PATH}/llvm-profdata merge
                                -sparse -f ${CMAKE_BINARY_DIR}/profiles.manifest -o ${CMAKE_BINARY_DIR}/merged.profdata
                                && ${CMAKE_CXX_COMPILER_PATH}/llvm-cov show ${CMAKE_BINARY_DIR}/utils/libutils.dylib
                                -format html -instr-profile ${CMAKE_BINARY_DIR}/merged.profdata -o
                                ${CMAKE_BINARY_DIR}/llvm-cov-report -show-line-counts-or-regions -Xdemangler c++filt
                                -Xdemangler -n -show-instantiation-summary -show-instantiations
                                --ignore-filename-regex="${CMAKE_SOURCE_DIR}/packages/*"
                                --ignore-filename-regex="${CMAKE_SOURCE_DIR}/third_party/*"
                                --ignore-filename-regex="/usr/local/*"
                                --ignore-filename-regex="/Applications/Xcode.app/*"
                                --ignore-filename-regex="/Library/Developer/CommandLineTools/*" \"
                        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
                        )

      add_dependencies(generate-llvm-coverage-report check-all)
      add_dependencies(generate-coverage-reports generate-llvm-coverage-report)

      list(APPEND CURR_ADDIDITIONAL_CLEAN_FILES ${LLVM_COV_ARTIFACTS})

    else()
      message(WARNING "compiler ${CMAKE_CXX_COMPILER} doesn't support clang style source based coverage")
    endif(llvm_cov_exists)
  endif(coverage_style_llvm)
  set_property(DIRECTORY ${CMAKE_BINARY_DIR} PROPERTY ADDITIONAL_CLEAN_FILES ${CURR_ADDIDITIONAL_CLEAN_FILES})
endif()

set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
configure_file("${CMAKE_SOURCE_DIR}/goback.sh" "${CMAKE_BINARY_DIR}/goback.sh" @ONLY)
