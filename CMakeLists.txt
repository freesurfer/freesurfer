# before starting the project, set gcc/g++ as default
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

# project info
project(freesurfer)
cmake_minimum_required(VERSION 3.9)
enable_language(C CXX Fortran)

# a few build options
option(MINIMAL "Only build core components" OFF)
option(BUILD_QDEC "Build QDEC" OFF)
option(BUILD_GUIS "Build GUIs" ON)
option(USE_OPENCL "Build OpenCL algorithms" OFF)
option(INFANT_MODULE "Include infant recon-all" OFF)
option(USE_QT4 "Use Qt4 instead of Qt5 for GUIs" OFF)
if(APPLE)
  option(BUILD_OSX_TK_GUIS "Compile OSX Tk GUIs" OFF)
endif()

# enable ctest
enable_testing()

# include our custom cmake functions (and others)
include(cmake/functions.cmake)
include(TestBigEndian)
include(CheckFunctionExists)

# prevents itk from importing as a system
set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)

# xxd is used to generate the helptext headers - make sure it exists
find_program(XXD xxd)
if(NOT XXD)
  message(FATAL_ERROR "The xxd program is required to build freesurfer")
endif()

# allow reading installation dir from environment var, but if not defined,
# then override default install dir /usr/local (too messy to install there)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(NOT "$ENV{FS_INSTALL_DIR}" STREQUAL "")
    set(CMAKE_INSTALL_PREFIX "$ENV{FS_INSTALL_DIR}" CACHE PATH "Copied from FS_INSTALL_DIR env variable" FORCE)
  else()
    set(CMAKE_INSTALL_PREFIX "/usr/local/freesurfer" CACHE PATH "Default install path" FORCE)
  endif()
endif()

# set the default build type to 'Release' for optimization purposes
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type selected - defaulting to Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Default build type" FORCE)
endif()

# create a custom target "make install-nmr" to only install nmr-component files 
add_custom_target(install-nmr COMMAND ${CMAKE_COMMAND} -D COMPONENT=nmr -P cmake_install.cmake)

# create build stamp during install
if(NOT FS_BUILD_STAMP)
  set(FS_BUILD_STAMP "freesurfer")
endif()
install(CODE "file(WRITE ${CMAKE_INSTALL_PREFIX}/build-stamp.txt ${FS_BUILD_STAMP}\\n)")

# --------------------------------------------------
#   library dependencies and third-party packages
# --------------------------------------------------

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
    if(APPLE)
      set(FS_PACKAGES_URL "http://surfer.nmr.mgh.harvard.edu/pub/data/fspackages/prebuilt/osx10.11-packages.tar.gz")
    else()
      set(FS_PACKAGES_URL "http://surfer.nmr.mgh.harvard.edu/pub/data/fspackages/prebuilt/centos7-packages.tar.gz")
    endif()
    message(FATAL_ERROR "External FreeSurfer developers must provide a valid FS_PACKAGES_DIR path! "
                        "Required FS dependencies can be built locally by running:\n"
                        "    packages/build_packages.sh \"/desired/packages/install/dir\"\n"
                        "Or they can be downloaded from here:\n"
                        "    ${FS_PACKAGES_URL}\n"
                        "After the dependencies have been setup, cmake should be configured with:\n"
                        "    -DFS_PACKAGES_DIR=\"/path/to/fs-packages\"")
  endif()
endif()

# all custom find-modules are stored in the cmake subdir
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# -------- zlib --------
find_package(ZLIB REQUIRED)
add_definitions(-DHAVE_ZLIB)

if(NOT APPLE)

  # -------- gfortran --------
  find_library(GFORTRAN_LIBRARIES HINTS /usr/lib/gcc/x86_64-linux-gnu/4.8 NAMES libgfortran.a gfortran)

  # -------- blas and lapack --------
  find_library(BLAS_LIBRARIES NAMES libblas.a)
  find_library(LAPACK_LIBRARIES NAMES liblapack.a)

endif()

# -------- jpeg --------
find_package(JPEG REQUIRED)

# -------- tiff --------
find_package(TIFF REQUIRED)

# -------- expat --------
find_package(Expat REQUIRED)

# -------- xml2 --------
find_package(XML2 REQUIRED)

# -------- tetgen --------
find_package(Tetgen REQUIRED)

# -------- petsc --------
find_package(PETSC)

# -------- itk --------
if(NOT ITK_DIR)
  set(ITK_DIR ${FS_PACKAGES_DIR}/itk/4.13.0)
endif()
find_package(ITK HINTS ${ITK_DIR} REQUIRED)
# todo: since ITK is now required, we can get rid of this macro
add_definitions(-DHAVE_ITK_LIBS)

# -------- boost --------
set(BOOST_ROOT ${FS_PACKAGES_DIR}/boost/1.53)
find_package(Boost COMPONENTS program_options)

# -------- ann --------
# ANN is only used by mris_resample
find_package(ANN)

# -------- gts --------
# only used by mris_decimate and mris_decimate_gui
find_package(GTS)

# -------- openMP --------
find_package(OpenMP)
if(OPENMP_FOUND)
  add_definitions(-DHAVE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

# -------- openCV --------
# version issue for ./freesurfer/histo_register
if(NOT APPLE)
  find_package(OpenCV)
endif()

# -------- openCL --------
if(USE_OPENCL)
  find_package(OpenCL)
endif()

# -------- X11 ---------
find_package(X11 REQUIRED)

# -------- vtk --------
find_package(VTK)

# search for GUI libraries
if(BUILD_GUIS)

  # -------- openGL --------
  find_package(OpenGL)
  if(OPENGL_FOUND)
    add_definitions(-DHAVE_OPENGL)
    # since we have OpenGL, make sure we have glut too
    find_package(GLUT REQUIRED)
    add_definitions(-DHAVE_GLUT_LIBS)
  endif()

  # -------- qt --------
  find_package(QT)

  # -------- tcltktixblt --------
  find_package(TCLTKTIXBLT)

  # -------- kwwidgets --------
  # once qdec is deprecated, KWWidgets can be removed
  if(BUILD_QDEC)
    find_package(KWWidgets)
  endif()

  # -------- wxwidgets --------
  # only used by mris_decimate_gui
  find_package(WXWidgets)

endif()

# --------------------------------------------------
#             global build configuration
# --------------------------------------------------

add_definitions(-D${CMAKE_SYSTEM_NAME})

# general check for big endian
test_big_endian(IS_BIG_ENDIAN)

# compiler warnings
set(C_CXX_FLAGS "${C_CXX_FLAGS} -Wall -Wno-unused-result -Wno-unused-local-typedefs")
# clang complains about -Wno-unused-but-set-variable and says to use -Wno-unused-const-variable
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(C_CXX_FLAGS "${C_CXX_FLAGS} -Wno-unused-const-variable")
else()
  set(C_CXX_FLAGS "${C_CXX_FLAGS} -Wno-unused-but-set-variable")
endif()

# ANSI
add_definitions(-DANSI)

# SSE matrix and math functions (affine.h and sse_mathfun.h)
# todo: write some logic to set this correctly
add_definitions(-DUSE_SSE_MATHFUN)
set(C_CXX_FLAGS "${C_CXX_FLAGS} -msse2 -mfpmath=sse")

# large file support: support for files > 2GB
if(APPLE)
  add_definitions(-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE)
endif()

# todo: implement processor-based logic to determine this
set(C_CXX_FLAGS "${C_CXX_FLAGS} -m64")

# for stripping unused code
if(APPLE)
  add_subdirectory(histo_register)
  set(STRIP_FLAGS "-dead_strip")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${STRIP_FLAGS} -static-libgcc")
else()
  set(STRIP_FLAGS "-fdata-sections -ffunction-sections -Wl,--gc-sections")
  set(C_CXX_FLAGS "${C_CXX_FLAGS} ${STRIP_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${STRIP_FLAGS}")
endif()

if(NOT APPLE)
  # link map (only for linux)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-Map,ld_map.txt -Wl,--no-demangle")
  # CentOS/Redhat
  if(EXISTS "/etc/redhat-release")
     add_subdirectory(histo_register)
     # Try to use native libGLdispatch.so as versions from /usr/lib64/nvidia may be missing symbols
     # Alternately add /lib64 to global CMAKE_EXE_LINKER_FLAGS
     if(EXISTS "/lib64/libGLdispatch.so")
       set(CMAKE_EXE_LINKER_FLAGS "-lGLdispatch ${CMAKE_EXE_LINKER_FLAGS}")
       message("CentOS/Redhat build will use the native libGLdispatch.so if needed")
     endif()
  # Ubuntu
  elseif(EXISTS "/etc/os-release")
     # add_subdirectory(histo_register EXCLUDE_FROM_ALL)
     # list(REMOVE_ITEM histo_register ????)
     message("Ubuntu build will exclude subdirectory ./freesurfer/histo_register")
  else()
     add_subdirectory(histo_register)
  endif()
endif()

# check for clock_gettime on mac
if(APPLE)
  check_function_exists(clock_gettime HAVE_CLOCK_GETTIME)
  if(HAVE_CLOCK_GETTIME)
    add_definitions(-DHAVE_CLOCK_GETTIME)
  endif()
endif()

# apply C and CXX flags
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${C_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${C_CXX_FLAGS}")

# enable std c++11
set(CMAKE_CXX_STANDARD 11)

# set executable rpath
if(NOT APPLE)
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:${CMAKE_INSTALL_RPATH}")
endif()

# --------------------------------------------------
#                   setup python
# --------------------------------------------------

# Unfortunately, the python version used to run wrapped c-libraries must be equal to
# the version used to build the libraries. The easiest and least intrusive way of making freesurfer
# python scripts run out-of-the-box (and to help guarantee reproducibility) requires
# distributing a minimal, custom python installation called fspython. This fspython package
# will get installed by default to "freesurfer/lib/python", but an external python can be
# used instead by turning off the DISTRIBUTE_FSPYTHON option. Turning this off will create a
# freesurfer/bin/fspythonlink symlink that points to the external python executable located by pybind
option(DISTRIBUTE_FSPYTHON "Include the fspython distribution in the installation" ON)

if(DISTRIBUTE_FSPYTHON)
  set(FSPYTHON_DIR ${FS_PACKAGES_DIR}/fspython/3.6/base)
  if(EXISTS ${FSPYTHON_DIR}/bin/python)
    # link the fspython distribution
    symlink(${FSPYTHON_DIR} ${CMAKE_INSTALL_PREFIX}/lib/python)
    # specify the python to use
    set(PYTHON_EXECUTABLE ${FSPYTHON_DIR}/bin/python)
    set(PYTHON_LIBRARY ${FSPYTHON_DIR}/lib/libpython3.6m.so)
  else()
    message(WARNING "Cannot find fspython distribution in FS_PACKAGES_DIR, "
                    "so must use external python executable")
    set(DISTRIBUTE_FSPYTHON OFF)
  endif()
endif()

# initialize pybind for python wrapping
set(PYBIND11_PYTHON_VERSION 3.5)  # minimum version required
add_subdirectory(packages/pybind11)

if(NOT DISTRIBUTE_FSPYTHON)
  # link to an external python binary if fspython won't be distributed
  symlink(${PYTHON_EXECUTABLE} ${CMAKE_INSTALL_PREFIX}/bin/fspythonlink)
endif()

# --------------------------------------------------
#                 build freesurfer
# --------------------------------------------------

# first, prepare the freesurfer distribution
add_subdirectory(distribution)
if(NOT MINIMAL)
  add_subdirectory(images)
endif()

# the top-level include dir contain the most commonly
# included freesurfer header files
set(FS_INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/minc_1_5_1
  ${CMAKE_SOURCE_DIR}/netcdf_3_6_0_p1
)

# make sure to use the local minc source
add_definitions(-DUSE_LOCAL_MINC)

# now build the freesurfer static libraries (order here matters!)
add_subdirectories(
  minc_1_5_1
  netcdf_3_6_0_p1
  unix
  log
  rgb
  hipsstubs
  dicom
  utils
  utilscpp
  fsgdf
  svm-lib
  opencl_algorithms
)

# the following utility libraries are required for freeview and qdec
add_subdirectories(vtkfsio vtkutils lineprof)

# the fem elastic subdir also builds a library, so it should be added first
add_subdirectory(fem_elastic)

# build the gems python library
add_subdirectories(GEMS2)

# --------------------------------------------------
#                      programs
# --------------------------------------------------

# the following subdirectories contain programs required for recon-all
# and should be included during a "minimal" install. Extra programs that aren't
# used in the standard stream should be added further down
add_subdirectories(
  lta_convert
  mri_and
  mri_add_xform_to_header
  mri_annotation2label
  mri_aparc2aseg
  mri_aparc2wmseg
  mri_binarize
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
  mris_ca_label
  mris_calc
  mris_convert
  mris_curvature
  mris_curvature_stats
  mris_decimate
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

# the following program subdirectories aren't required in the standard
# recon-all stream - they will be built by default, but not if a minimal build
# is configured
if(NOT MINIMAL)
  add_subdirectories(
    anatomicuts
    BrainstemSS
    check_siemens_dir
    connectgraph
    diffusion_tool
    dmri_tensoreig
    dngtester
    freeview
    fsfast
    fslutils
    hiam_make_surfaces
    hiam_make_template
    hiam_register
    HippoSF
    histo_compute_joint_density
    histo_fix_topology
    histo_register_block
    histo_segment
    histo_synthesize
    label2flat
    label2patch
    label_area
    label_border
    mkxsubjreg
    matlab
    mri_apply_inu_correction
    mri_aseg_edit_reclassify
    mri_aseg_edit_train
    mri_auto_fill
    mri_average
    mri_bc_sc_bias_correct
    mri_bias
    mri_build_priors
    mri_ca_tissue_parms
    mri_ca_train
    mri_cal_renormalize_gca
    mri_cht2p
    mri_classify
    mri_cnr
    mri_compile_edits
    mri_compute_change_map
    mri_compute_structure_transforms
    mri_compute_volume_fractions
    mri_concatenate_gcam
    mri_convert_mdh
    mri_copy_values
    mri_cor2label
    mri_correct_segmentations
    mri_cvs_register
    mri_dct_align
    mri_distance_transform
    mri_divide_segmentation
    mri_elastic_energy
    mri_estimate_tissue_parms
    mri_evaluate_morph
    mri_extract
    mri_extract_conditions
    mri_extract_fcd_features
    mri_extract_label
    mri_fcili
    mri_fdr
    mri_fieldsign
    mri_fit_bias
    mri_fslmat_to_lta
    mri_fuse_intensity_images
    mri_gca_ambiguous
    mri_gcab_train
    mri_gdfglm
    mri_glmfit
    mri_gtmpvc
    mri_gtmseg
    mri_hausdorff_dist
    mri_head
    mri_hires_register
    mri_histo_eq
    mri_histo_normalize
    mri_ibmc
    mri_interpolate
    mri_jacobian
    mri_joint_density
    mri_label_accuracy
    mri_label_histo
    mri_label_vals
    mri_label_volume
    mri_linear_register
    mri_make_bem_surfaces
    mri_make_density_map
    mri_make_labels
    mri_make_register
    mri_make_template
    mri_map_cpdat
    mri_mark_temporal_lobe
    mri_mcsim
    mri_mi
    mri_modify
    mri_morphology
    mri_mosaic
    mri_ms_EM
    mri_ms_EM_with_atlas
    mri_ms_fitparms
    mri_ms_LDA
    mri_multiscale_segment
    mri_multispectral_segment
    mri_nlfilter
    mri_otl
    mri_paint
    mri_parse_sdcmdir
    mri_parselabel
    mri_partial_ribbon
    mri_path2label
    mri_polv
    mri_probe_ima
    mri_reduce
    mri_register
    mri_rf_label
    mri_rf_long_label
    mri_rf_long_train
    mri_rf_train
    mri_ribbon
    mri_rigid_register
    mri_sbbr
    mri_seg_diff
    mri_segcentroids
    mri_seghead
    mri_segment_tumor
    mri_segment_wm_damage
    mri_simulate_atrophy
    mri_strip_nonwhite
    mri_strip_subject_info
    mri_surfacemask
    mri_synthesize
    mri_threshold
    mri_topologycorrection
    mri_train
    mri_train_autoencoder
    mri_transform
    mri_transform_to_COR
    mri_twoclass
    mri_update_gca
    mri_vol2roi
    mri_volcluster
    mri_volsynth
    mri_warp_convert
    mri_wbc
    mri_wmfilter
    mri_xcorr
    mri_xvolavg
    mri_z2p
    mris2rgb
    mris_add_template
    mris_annot_to_segmentation
    mris_apply_reg
    mris_aseg_distance
    mris_average_curvature
    mris_average_parcellation
    mris_BA_segment
    mris_ca_train
    mris_classify_thickness
    mris_compute_acorr
    mris_compute_lgi
    mris_compute_optimal_kernel
    mris_compute_overlap
    mris_compute_parc_overlap
    mris_compute_volume_fractions
    mris_congeal
    mris_copy_header
    mris_decimate_gui
    mris_defects_pointset
    mris_deform
    mris_density
    mris_distance_map
    mris_distance_to_label
    mris_distance_transform
    mris_entropy
    mris_errors
    mris_expand
    mris_extract_patches
    mris_fbirn_annot
    mris_fill
    mris_find_flat_regions
    mris_flatten
    mris_fwhm
    mris_glm
    mris_hausdorff_dist
    mris_init_global_tractography
    mris_interpolate_warp
    mris_label_area
    mris_label_calc
    mris_label_mode
    mris_longitudinal_surfaces
    mris_make_average_surface
    mris_make_face_parcellation
    mris_make_map_surfaces
    mris_make_template
    mris_map_cuts
    mris_merge_parcellations
    mris_mesh_subdivide
    mris_morph_stats
    mris_ms_refine
    mris_ms_surface_CNR
    mris_multiscale_stats
    mris_niters2fwhm
    mris_parcellate_connectivity
    mris_pmake
    mris_register_label_map
    mris_register_to_volume
    mris_remove_variance
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
    mris_spherical_average
    mris_surface_change
    mris_surface_to_vol_distances
    mris_svm_classify
    mris_svm_train
    mris_talairach
    mris_thickness_comparison
    mris_transform
    mris_translate_annotation
    mris_transmantle_dysplasia_paths
    mris_twoclass
    mris_volume
    mris_w_to_curv
    mris_warp
    mris_watershed
    mrisp_write
    nmovie_qt
    oct_register_mosaic
    oct_train
    optseq2
    qdecproject
    qdec
    qdec_glmfit
    rbftest
    regdat2xfm
    repair_siemens_file
    resurf
    spherical_stats
    spline3
    stat_normalize
    stim_polar
    swi_processing
    template
    test_makevol
    ThalamicNuclei
    tkmedit
    # tksurfer
    trc
    tridec
  )
endif()

# infant freesurfer
if(INFANT_MODULE)
  add_subdirectories(infant)
endif()
