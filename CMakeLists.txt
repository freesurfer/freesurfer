# project info
cmake_minimum_required(VERSION 3.9)

# requiring at least c++11 standard
if(CMAKE_CXX_STANDARD EQUAL "98" )
   message(FATAL_ERROR "CMAKE_CXX_STANDARD:STRING=98 is not supported in.")
endif()

# set the default language dialect settings (This must be done BEFORE the project call)
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(NOT CMAKE_CXX_EXTENSIONS)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

project(freesurfer 
  DESCRIPTION "An open source software suite for processing and analyzing (human) brain MRI images."
  LANGUAGES C CXX Fortran
)

# a few build options
option(MINIMAL "Only build core components" OFF)
option(BUILD_GUIS "Build GUI tools" ON)
option(BUILD_ATTIC "Build deprecated tools from the attic" OFF)
option(INFANT_MODULE "Include infant recon-all" OFF)
option(QATOOLS_MODULE "Install QA tools" ON)
option(BUILD_DNG "Build Doug's testing tools" OFF)
option(FREEVIEW_LINEPROF "Build FreeView with lineprof enabled" OFF)
option(PROFILING "Complile binaries for profiling with gprof" OFF)
option(INSTALL_PYTHON_DEPENDENCIES "Install python package dependencies" ON)
option(BUILD_FORTRAN "Build subdirs with source using gfortran" ON)
option(PATCH_FSPYTHON "Build subdirs with source using gfortran" OFF)

#  Modern compilers on new platforms often provide many more warnings than older compilers
#  for example the gcc 7.5.0 compiler warnings about falling through switch statements
#  without explicit breaks on every case.  This desired fallthrough behavior causes false
#  compiler errors with "-Werror" turned on.  The WARNING_AS_ERROR option allows temporarily
#  disabling warnings as errors to allow a more gradual move to supporting new compilers.
option(WARNING_AS_ERROR "Convert build warnings to errors" ON)
option(MARTINOS_BUILD "Does the build have access to the martinos network" ON)

if(NOT APPLE)
  # linux-only build options
  option(TKTOOLS_MODULE "Install old Linux TK GUIs" OFF)
endif()

# enable ctest
enable_testing()

# include our custom cmake functions (and others)
include(cmake/functions.cmake)
include(TestBigEndian)
include(CheckFunctionExists)

# prevent third-party packages from importing as a system
set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)

# if an install prefix is not provided, check the FS_INSTALL_DIR env var, and
# if that is not defined, set the default path to /usr/local/freesurfer
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(NOT "$ENV{FS_INSTALL_DIR}" STREQUAL "")
    set(CMAKE_INSTALL_PREFIX "$ENV{FS_INSTALL_DIR}" CACHE PATH "Copied from FS_INSTALL_DIR env variable" FORCE)
  else()
    set(CMAKE_INSTALL_PREFIX "/usr/local/freesurfer" CACHE PATH "Default install path" FORCE)
  endif()
endif()

# version stamp
set(FS_VERSION "$ENV{USER}-local" CACHE STRING "Distribution version")

# build stamp
string(TIMESTAMP TODAY "%Y%m%d")
set(BUILD_STAMP "freesurfer-local-build-${TODAY}" CACHE STRING "Distribution build stamp")
install(CODE "file(WRITE ${CMAKE_INSTALL_PREFIX}/build-stamp.txt ${BUILD_STAMP}\\n)")

# set the default build type to 'Release' for optimization purposes
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type selected - defaulting to Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Default build type" FORCE)
endif()

# create a custom target "make install-nmr" to only install nmr-component files 
add_custom_target(install-nmr COMMAND ${CMAKE_COMMAND} -D COMPONENT=nmr -P cmake_install.cmake)

# --------------------------------------------------
#               external dependencies
# --------------------------------------------------

# xxd is used to generate the helptext headers
find_program(XXD xxd)
if(NOT XXD)
  message(WARNING "MISSING xxd - required to build executable help texts")
endif()

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
    set(FS_PACKAGES_DIR "$ENV{FS_PACKAGES_DIR}" CACHE INTERNAL "Copied from FS_PACKAGES_DIR environment variable")
  elseif(EXISTS /usr/pubsw/packages)
    # if the user is at Martinos, default to /usr/pubsw/packages
    set(FS_PACKAGES_DIR /usr/pubsw/packages CACHE INTERNAL "Default Martinos packages dir")
  else()
    # setting FS_PACKAGES_DIR is required if the user is outside Martinos
    # for more information visit https://surfer.nmr.mgh.harvard.edu/fswiki/BuildRequirements
    message(WARNING "FS_PACKAGES_DIR is NOT defined - build will be limited")
  endif()
endif()

# all custom find-modules are stored in the cmake subdir
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Use host_os() to get HOST_OS to reference the correct 3rd party package build
host_os()

# -------- zlib --------
find_package(ZLIB REQUIRED)
add_definitions(-DHAVE_ZLIB)

# -------- gfortran/blas/lapack --------
if(NOT APPLE)
  if(HOST_OS MATCHES "Ubuntu20" OR HOST_OS MATCHES "Ubuntu22")
     # assume gcc8 tools distribution on Ubuntu 20 which contains its own libquadmath and libgfortran
     find_library(GFORTRAN_LIBRARIES HINTS /usr/lib/gcc/x86_64-linux-gnu/8 NAMES libgfortran.a gfortran)
     find_library(BLAS_LIBRARIES NAMES libopenblas.a libblas.a)
     find_library(LAPACK_LIBRARIES NAMES liblapack.a)
     find_library(QUADMATH_LIBRARIES HINTS /usr/lib/gcc/x86_64-linux-gnu/8 NAMES libquadmath.a libquadmath.so libquadmath.so.0)
  else()
     find_library(GFORTRAN_LIBRARIES HINTS /usr/lib64 /usr/lib /usr/lib/gcc/*/*/ NAMES libgfortran.a gfortran)
     find_library(BLAS_LIBRARIES NAMES libopenblas.a libblas.a)
     find_library(LAPACK_LIBRARIES NAMES liblapack.a)
     find_library(QUADMATH_LIBRARIES HINTS /usr/lib64 /usr/lib /usr/lib/gcc/*/*/ NAMES libquadmath.a libquadmath.so libquadmath.so.0)
  endif()
  if(NOT QUADMATH_LIBRARIES)
    set(QUADMATH_LIBRARIES "")
  endif()
endif()

# commented out on 10/2022
# -------- armadillo --------
# find_package(ARMADILLO)

# -------- petsc --------
find_package(PETSC)

# -------- itk --------
if(NOT ITK_DIR)
   if(HOST_OS MATCHES "CentOS8" OR HOST_OS MATCHES "CentOS8stream")
      set(ITK_DIR ${FS_PACKAGES_DIR}/itk/4.13.2_centos8 CACHE PATH "ITK install directory")
   elseif(HOST_OS MATCHES "Ubuntu20" OR HOST_OS MATCHES "Ubuntu22")
      set(ITK_DIR ${FS_PACKAGES_DIR}/itk/4.13.2_ubuntu20 CACHE PATH "ITK install directory")
   elseif(HOST_OS MATCHES "CentOS7")
      set(ITK_DIR ${FS_PACKAGES_DIR}/itk/4.13.2_centos7 CACHE PATH "ITK install directory")
   elseif(HOST_OS MATCHES "Ubuntu18")
      set(ITK_DIR ${FS_PACKAGES_DIR}/itk/4.13.2_ubuntu18 CACHE PATH "ITK install directory")
   elseif(APPLE)
      set(ITK_DIR ${FS_PACKAGES_DIR}/itk/4.13.2_darwin CACHE PATH "ITK install directory")
   else()
      # default
      set(ITK_DIR ${FS_PACKAGES_DIR}/itk/4.13.2 CACHE PATH "ITK install directory")
   endif()
endif()

find_package(ITK HINTS ${ITK_DIR} REQUIRED)
message(STATUS "For HOST_OS=${HOST_OS} ITK_DIR=${ITK_DIR}")

# -------- ann --------
# ANN is only used by mris_resample
find_package(ANN)

# -------- openMP --------
if(APPLE)
  # this is an ugly hack to ensure that cmake looks for the
  # static version of libomp on mac, since we don't want
  # to distribute the shared library
  set(ORIG_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  find_package(OpenMP)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${ORIG_LIBRARY_SUFFIXES})

  # With added security in MacOS 12.X (Monterey)
  # and 13.X (Ventura), permissions may be broken under
  # /usr/local for accessing gfortran/libs via Homebrew. So for
  # now, turn off build/install of subdirs using gfortran.
  if(HOST_OS MATCHES "MacOS-12" OR HOST_OS MATCHES "MacOS-13")
     set(BUILD_FORTRAN OFF)
  endif()
else()
  find_package(OpenMP)
endif()

if(OPENMP_FOUND)
  add_definitions(-DHAVE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(OMP_CXX_LIBRARIES OpenMP::OpenMP_CXX)
endif()

# Avoid error from __powf_finite undefined.
# The 2.31 version of the GNU C Library removed symbols from LLVM libc++ math.h including,
# __powf_finite __asinf_finite U_atan2f_finite __exp2f_finite __exp2_finite __expf_finite __log2f_finite
if(HOST_OS MATCHES "Ubuntu20" OR HOST_OS MATCHES "Ubuntu22")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-builtin")
endif()

# -------- openCV --------
if(NOT APPLE)
  find_package(OpenCV)
endif()

# -------- X11 ---------
find_package(X11)

# search for GUI libraries
if(BUILD_GUIS)

  # -------- openGL --------
  if(EXISTS /usr/lib64/nvidia/libGL.so)
    # on some martinos centos7 machines, linking to the default libGL.so fails due to missing symbols
    # in the nvidia drivers, so as a temporary fix, link directly to the nvidia GL library
    list(APPEND CMAKE_PREFIX_PATH "/usr/lib64/nvidia")
  endif()

  # Output from ...
  # $ cmake --help-policy CMP0072
  # 
  # CMake 3.11 and above prefer to choose GLVND libraries.  This policy provides
  # compatibility with projects that expect the legacy GL library to be used.
  # The ``OLD`` behavior for this policy is to set ``OpenGL_GL_PREFERENCE`` to
  # ``LEGACY``.  The ``NEW`` behavior for this policy is to set
  # ``OpenGL_GL_PREFERENCE`` to ``GLVND``.

  set(OpenGL_GL_PREFERENCE "LEGACY")
  # set(OpenGL_GL_PREFERENCE "GLVND")

  find_package(OpenGL)
  if(OPENGL_FOUND)
    add_definitions(-DHAVE_OPENGL)
  endif()

  # -------- qt --------

if(NOT Qt6_DIR)
  if(NOT Qt5_DIR)
     if(HOST_OS MATCHES "CentOS7" OR HOST_OS MATCHES "CentOS8" OR HOST_OS MATCHES "CentOS8stream")
        set(Qt5_DIR ${FS_PACKAGES_DIR}/qt/5.12.11_LTS_cos7/lib/cmake/Qt5 CACHE PATH "Qt5 install directory")
     elseif(HOST_OS MATCHES "Ubuntu18" OR HOST_OS MATCHES "Ubuntu20" OR HOST_OS MATCHES "Ubuntu22")
        set(Qt5_DIR ${FS_PACKAGES_DIR}/qt/5.12.11_LTS_ub18/lib/cmake/Qt5 CACHE PATH "Qt5 install directory")
     else()
        # default
        set(Qt5_DIR ${FS_PACKAGES_DIR}/qt/5.12.11_LTS_cos7/lib/cmake/Qt5 CACHE PATH "Qt5 install directory")
     endif()
  endif()
endif()

  find_package(QT)

  if(Qt6_DIR)
     message(STATUS "For HOST_OS=${HOST_OS} Qt6_DIR=${Qt6_DIR}")
  elseif(Qt5_DIR)
     message(STATUS "For HOST_OS=${HOST_OS} Qt5_DIR=${Qt5_DIR}")
  endif()

endif()

# -------- vtk --------

if(NOT VTK_DIR)
   set(VTK_DIR ${FS_PACKAGES_DIR}/vtk/7.1 CACHE PATH "VTK install directory")
endif()

find_package(VTK HINTS ${VTK_DIR})
message(STATUS "For HOST_OS=${HOST_OS} VTK_DIR=${VTK_DIR}")

if(VTK_FOUND AND VTK_BUILD_SHARED_LIBS)
  # locate all shared vtk libraries
  if(${VTK_VERSION} LESS 7)
    file(GLOB _vtk_library_list "${VTK_LIBRARY_DIRS}/lib*.so*")
  else()
    # special libvtkproj
    file(GLOB _found_library "${VTK_LIBRARY_DIRS}/libvtkproj*.so*")
    list(APPEND _vtk_library_list ${_found_library})
    # special libvtkNetCDF
    file(GLOB _found_library "${VTK_LIBRARY_DIRS}/libvtkNetCDF*.so*")
    list(APPEND _vtk_library_list ${_found_library})
    # the rest of the VTK libraries
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

add_subdirectory(packages)

# --------------------------------------------------
#                    setup python
# --------------------------------------------------

# Unfortunately, the python version used to run pybind c-libraries must be equivalent to
# the version used to build the libraries. The easiest and least intrusive way of making freesurfer
# python scripts run out-of-the-box (and to help guarantee reproducibility) requires
# distributing a minimal, custom python installation called fspython. This fspython package
# will get installed to "freesurfer/python", but an external python can be used instead when the
# DISTRIBUTE_FSPYTHON option is turned off. Turning this off will create a freesurfer/bin/fspythonlink
# symlink during install that points to the external python executable located by pybind
option(DISTRIBUTE_FSPYTHON "Include the fspython distribution in the installation" OFF)

if(DISTRIBUTE_FSPYTHON)
  set(FSPYTHON_DIR ${FS_PACKAGES_DIR}/fspython/3.8)
  if(EXISTS ${FSPYTHON_DIR}/bin/python3)
    # specify the python to use
    set(PYTHON_EXECUTABLE ${FSPYTHON_DIR}/bin/python3)
    set(PYTHON_LIBRARY ${FSPYTHON_DIR}/lib/libpython3.8m)
    # install the barebones python distribution
    INSTALL(DIRECTORY ${FSPYTHON_DIR}/bin ${FSPYTHON_DIR}/lib ${FSPYTHON_DIR}/include
      DESTINATION ${CMAKE_INSTALL_PREFIX}/python USE_SOURCE_PERMISSIONS
    )
  else()
    message(FATAL_ERROR "Cannot find fspython distribution in FS_PACKAGES_DIR")
  endif()
endif()

# initialize pybind for python wrapping
set(PYBIND11_PYTHON_VERSION 3)
add_subdirectory(packages/pybind11)

if(NOT DISTRIBUTE_FSPYTHON)
  # link to an external python binary if fspython won't be distributed
  symlink(${PYTHON_EXECUTABLE} ${CMAKE_INSTALL_PREFIX}/python/bin/python3)
endif()

add_custom_target(install-fspy COMMAND ${CMAKE_COMMAND} -D COMPONENT=fspy -P cmake_install.cmake)

# cmake apparently generates incorrect dependencies/paths
if(APPLE)
   if(HOST_OS MATCHES "MacOS-12" OR HOST_OS MATCHES "MacOS-13")
      add_custom_target(check_macos_depends ALL COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/macos_depends.sh ${CMAKE_SOURCE_DIR})
   endif()
endif()

# --------------------------------------------------
#             freesurfer build settings
# --------------------------------------------------

# warnings
if(NOT NO_CXX_WARN)
  if(WARNING_AS_ERROR)
    set(WARN_AS_ERROR_FLAG "-Werror")
    if(HOST_OS MATCHES "CentOS8" OR HOST_OS MATCHES "CentOS8stream")
       # -Wno-cpp needed with gcc8 for /usr/pubse/packages/vtk 7.1 compiled with older compilers
       # -Wno-restrict needed with gcc8 for ./freesurfer/mris_compute_acorr/mris_compute_acorr.cpp
       # -Wno-format-overflow needed with gcc8 for ./fresurfer/mris2rgb/mris2rgb.cpp
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra ${WARN_AS_ERROR_FLAG} -Wno-cpp -Wno-restrict -Wno-format-overflow -Wno-sign-compare -Wno-unused-result -Wno-unused-parameter")
    elseif(HOST_OS MATCHES "Ubuntu20" OR HOST_OS MATCHES "Ubuntu22")
       # uses gcc8 as well
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin -Wall -Wextra ${WARN_AS_ERROR_FLAG} -Wno-cpp -Wno-restrict -Wno-format-overflow -Wno-sign-compare -Wno-unused-result -Wno-unused-parameter")
    elseif(HOST_OS MATCHES "MacOS-12" OR HOST_OS MATCHES "MacOS-13")
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-non-c-typedef-for-linkage")
    else ()
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra ${WARN_AS_ERROR_FLAG} -Wno-sign-compare -Wno-unused-result -Wno-unused-parameter")
    endif()
  endif()
else ()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error -Wno-deprecated-declarations -fpermissive") 
endif()

# clang complains about -Wno-unused-but-set-variable and says to use -Wno-unused-const-variable
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-const-variable -Wno-inconsistent-missing-override -Wno-self-assign-field -Wno-tautological-overlap-compare -Wno-tautological-compare -Wno-unused-value -Wno-range-loop-analysis -Wno-return-stack-address -Wno-dangling-gsl")
else()
  if(NOT NO_CXX_WARN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-but-set-variable")
  endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# linker options
if(APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -dead_strip")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
  endif()
else()
  set(STRIP_FLAGS "-fdata-sections -ffunction-sections -Wl,--gc-sections")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STRIP_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${STRIP_FLAGS} -Wl,-Map,ld_map.txt -Wl,--no-demangle")
endif()

if(PROFILING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

message(STATUS "For HOST_OS=${HOST_OS} CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
message(STATUS "For HOST_OS=${HOST_OS} CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")

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
)

# build the freesurfer static libraries
add_subdirectory(utils)

# the following utility libraries are required for freeview
add_subdirectories(vtkutils lineprof)

# the fem elastic subdir also builds a library, so it should be added first
add_subdirectory(fem_elastic)

# build the gems library
add_subdirectory(gems)

# --------------------------------------------------
#                      programs
# --------------------------------------------------

# the following subdirectories contain programs required for recon-all
# and should be included during a "minimal" install. Extra programs that aren't
# used in the standard stream should be added further down
add_subdirectories(
  AntsDenoiseImageFs
  AntsN4BiasFieldCorrectionFs
  lta_convert
  io
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
  mris_skeletonize
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
  samsegmesh2surf
  talairach_afd
  scripts
  tkregister2
)
if(BUILD_FORTRAN)
  add_subdirectories(
    talairach_avi
  )
else()
  message(STATUS "For HOST_OS=${HOST_OS} NOT BUILDING talairach subdirs requiring gfortran.")
endif()

# the following program subdirectories aren't required in the standard recon-all stream.
# they will be built by default, but not if a minimal build is configured
if(NOT MINIMAL)
  add_subdirectories(
    anatomicuts
    attic
    BrainstemSS
    diffusion_tool
    dissection_photo
    dummy
    exvivo
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
    mri_easyreg
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
    mri_gradunwarp
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
    mri_sclimbic_seg
    mri_segcentroids
    mri_seghead
    mri_segment_hypothalamic_subunits
    mri_strip_nonwhite
    mri_strip_subject_info
    mri_surfacemask
    mri_synthesize
    mri_synthmorph
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
    mris_estimate_wm
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
    mris_reposition_surface
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
    subregions
    stat_normalize
    stim_polar
    mri_synthstrip
    mri_synthseg
    mri_synthsr
    mri_3d_photo_recon
    swi_processing
    template
    test_makevol
    ThalamicNuclei
    trc
    tridec
    ventfix
    pointset2label
  )
endif()

# # Only build some old subdirectories if they exist
# # AND the build host is a certain version of Cent/Redhat OS
# if(EXISTS "/etc/redhat-release")
#    execute_process(COMMAND cat \/etc\/redhat\-release OUTPUT_VARIABLE REDHAT_VERSION)
#    string(STRIP ${REDHAT_VERSION} REDHAT_VERSION)
#    # message(STATUS "Host is ${REDHAT_VERSION}")
#    if(EXISTS "./mri_make_bem_surfaces")
#       if(REDHAT_VERSION MATCHES "release 7")
#          add_subdirectories(mri_make_bem_surfaces)
#       elseif(REDHAT_VERSION MATCHES "release 8")
#          add_subdirectories(mri_make_bem_surfaces)
#       endif()
#    endif()
# endif()

# Only build some old subdirectories if they exist
# AND the build host is a certain version of Cent/Redhat OS
if(EXISTS "./mri_make_bem_surfaces")
   if(HOST_OS MATCHES "CentOS7")
      add_subdirectories(mri_make_bem_surfaces)
      message(STATUS "Adding mri_make_bem_surfaces to build on ${HOST_OS}")
   elseif(HOST_OS MATCHES "CentOS8")
      add_subdirectories(mri_make_bem_surfaces)
      message(STATUS "Adding mri_make_bem_surfaces to build on ${HOST_OS}")
   elseif(HOST_OS MATCHES "Ubuntu20" OR HOST_OS MATCHES "Ubuntu22")
      add_subdirectories(mri_make_bem_surfaces)
      message(STATUS "Adding mri_make_bem_surfaces to build on ${HOST_OS}")
   elseif(APPLE)
      add_subdirectories(mri_make_bem_surfaces)
      message(STATUS "Adding mri_make_bem_surfaces to build on ${HOST_OS}")
   endif()
endif()

if(INFANT_MODULE)
  add_subdirectory(infant)
  add_subdirectory(sscnn_skullstripping)
endif()

if(INFANT_MODULE AND MINIMAL)
  add_subdirectory(mri_distance_transform)
  add_subdirectory(mri_correct_segmentations)
  add_subdirectory(mri_morphology)
endif()

if(TKTOOLS_MODULE)
  add_subdirectory(tktools)
endif()

if(QATOOLS_MODULE)
  add_subdirectory(qatools)
endif()

if(BUILD_DNG)
  add_subdirectory(dngtester)
endif()

