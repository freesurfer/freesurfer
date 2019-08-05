# before starting the project, set gcc/g++ as the default compiler
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

# project info
project(freesurfer)
cmake_minimum_required(VERSION 3.9)
enable_language(C CXX Fortran)

# a few build options
option(MINIMAL "Only build core components" OFF)
option(BUILD_GUIS "Build GUI tools" ON)
option(BUILD_ATTIC "Build subdirs moved to ./freesurfer/attic" OFF)
option(INFANT_MODULE "Include infant recon-all" OFF)
option(SUPPRESS_WARNINGS "Suppress some selected warnings" ON)
option(PROFILING "Complile binaries for profiling with gprof" OFF)

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

# create build stamp during install
string(TIMESTAMP TODAY "%Y%m%d")
set(BUILD_STAMP "freesurfer-local-build-dev-${TODAY}" CACHE STRING "Distribution build stamp")
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
  message(FATAL_ERROR "The xxd program is required to build freesurfer")
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

# -------- gfortran/blas/lapack --------
if(NOT APPLE)
  find_library(GFORTRAN_LIBRARIES
    HINTS /usr/lib/gcc/x86_64-linux-gnu/4.8 /usr/lib/gcc/x86_64-redhat-linux/4.8.2/
    NAMES libgfortran.a gfortran)
  find_library(BLAS_LIBRARIES NAMES libblas.a)
  find_library(LAPACK_LIBRARIES NAMES liblapack.a)
endif()

# -------- armadillo --------
find_package(Armadillo)

# -------- petsc --------
find_package(PETSC)

# -------- itk --------
set(ITK_DIR ${FS_PACKAGES_DIR}/itk/4.13.0 CACHE PATH "ITK install directory")
find_package(ITK HINTS ${ITK_DIR} REQUIRED)

# -------- ann --------
# ANN is only used by mris_resample
find_package(ANN)

# -------- gts --------
# only used by mris_decimate
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
if(NOT APPLE)
  find_package(OpenCV)
endif()

# -------- X11 ---------
find_package(X11 REQUIRED)

# search for GUI libraries
if(BUILD_GUIS)

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

endif()

# -------- vtk --------
set(VTK_DIR ${FS_PACKAGES_DIR}/vtk/5.10.1 CACHE PATH "VTK install directory")
find_package(VTK HINTS ${VTK_DIR})
if(VTK_FOUND AND NOT APPLE)
  # install the shared libraries to the freesurfer lib directory
  file(GLOB VTK_LIBRARIES_TO_INSTALL "${VTK_LIBRARY_DIRS}/lib*.so*")
  if(VTK_LIBRARIES_TO_INSTALL)
    install(PROGRAMS ${VTK_LIBRARIES_TO_INSTALL} DESTINATION lib/vtk)
    # add vtk library directory to rpath
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

# enable std c++11
set(CMAKE_CXX_STANDARD 11)

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
  set(FSPYTHON_DIR ${FS_PACKAGES_DIR}/fspython/3.6)
  if(EXISTS ${FSPYTHON_DIR}/bin/python3)
    # specify the python to use
    set(PYTHON_EXECUTABLE ${FSPYTHON_DIR}/bin/python3)
    set(PYTHON_LIBRARY ${FSPYTHON_DIR}/lib/libpython3.6m)
    # install the barebones python distribution
    INSTALL(DIRECTORY ${FSPYTHON_DIR}/bin ${FSPYTHON_DIR}/lib ${FSPYTHON_DIR}/include
      DESTINATION ${CMAKE_INSTALL_PREFIX}/python USE_SOURCE_PERMISSIONS
    )
  else()
    message(FATAL_ERROR "Cannot find fspython distribution in FS_PACKAGES_DIR")
  endif()
endif()

# initialize pybind for python wrapping
set(PYBIND11_PYTHON_VERSION 3.5)
add_subdirectory(packages/pybind11)

if(NOT DISTRIBUTE_FSPYTHON)
  # link to an external python binary if fspython won't be distributed
  symlink(${PYTHON_EXECUTABLE} ${CMAKE_INSTALL_PREFIX}/python/bin/python3)
endif()

# --------------------------------------------------
#             freesurfer build settings
# --------------------------------------------------

# warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wno-absolute-value -Wno-sign-compare -Wno-write-strings -Wno-unused-result -Wno-unused-parameter")

# clang complains about -Wno-unused-but-set-variable and says to use -Wno-unused-const-variable
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-const-variable")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-but-set-variable")
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

# check for clock_gettime on mac
if(APPLE)
  check_function_exists(clock_gettime HAVE_CLOCK_GETTIME)
  if(HAVE_CLOCK_GETTIME)
    add_definitions(-DHAVE_CLOCK_GETTIME)
  endif()
endif()

# --------------------------------------------------
#                 build freesurfer
# --------------------------------------------------

# prepare the freesurfer distribution
add_subdirectory(distribution)

# the top-level include dir contain the most commonly included freesurfer header files
set(FS_INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/packages/minc
  ${CMAKE_SOURCE_DIR}/packages/netcdf
  ${CMAKE_SOURCE_DIR}/packages/nifti
)

# build the freesurfer static libraries
add_subdirectories(utils)

# the following utility libraries are required for freeview
add_subdirectories(vtkutils lineprof)

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
  mris_annot_diff
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

# the following program subdirectories aren't required in the standard recon-all stream.
# they will be built by default, but not if a minimal build is configured
if(NOT MINIMAL)
  add_subdirectories(
    anatomicuts
    BrainstemSS
    diffusion_tool
    dummy
    dngtester
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
    mri_copy_values
    mri_cor2label
    mri_correct_segmentations
    mri_cvs_register
    mri_dct_align
    mri_distance_transform
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
    mri_seg_diff
    mri_segcentroids
    mri_seghead
    mri_strip_nonwhite
    mri_strip_subject_info
    mri_surfacemask
    mri_synthesize
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
    mris_defects_pointset
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

# No CMakeLists.txt in subdirs
# mris_spherical_wavelets

# infant freesurfer
if(INFANT_MODULE)
  add_subdirectories(infant)
endif()

# tk tools
if(TKTOOLS_MODULE)
  add_subdirectory(tktools)
endif()

# 5/2019 - The attic subdirectories were previously listed above to be built by default under
# ./freesurfer. But because the code is old/seldom used - or never built - these subdirs have been
# moved to ./freesurfer/attic (attic is an old CVS term).  They are no longer built and installed by 
# default, but can be built by adding -DBUILD_ATTIC:BOOL=ON to the cmake command line.
#
# To see the complete log history of any file in under ./attic prior to its move there, add the
# --follow option to the git log command, e.g., cd attic/mri_mosaic and do "git log --follow mri_mosaic.cpp"

if(BUILD_ATTIC)
  add_subdirectories(
    attic/check_siemens_dir
    attic/connectgraph
    attic/dmri_tensoreig
    attic/histo_fix_topology
    attic/histo_segment
    attic/label_area
    attic/label_border
    attic/mri_apply_inu_correction
    attic/mri_aseg_edit_reclassify
    attic/mri_aseg_edit_train
    attic/mri_auto_fill
    attic/mri_bc_sc_bias_correct
    attic/mri_build_priors
    attic/mri_cht2p
    attic/mri_classify
    attic/mri_compute_structure_transforms
    attic/mri_concatenate_gcam
    attic/mri_convert_mdh
    attic/mri_divide_segmentation
    attic/mri_elastic_energy
    attic/mri_estimate_tissue_parms
    attic/mri_extract_conditions
    attic/mri_fcili
    attic/mri_gcab_train
    attic/mri_gdfglm
    attic/mri_histo_normalize
    attic/mri_ibmc
    attic/mri_interpolate
    attic/mri_label_accuracy
    attic/mri_make_bem_surfaces
    attic/mri_make_density_map
    attic/mri_make_labels
    attic/mri_make_register
    attic/mri_make_template
    attic/mri_mosaic
    attic/mri_ms_EM
    attic/mri_ms_EM_with_atlas
    attic/mri_ms_LDA
    attic/mri_multiscale_segment
    attic/mri_multispectral_segment
    attic/mri_parselabel
    attic/mri_partial_ribbon
    attic/mri_register
    attic/mri_rf_long_label
    attic/mri_segment_tumor
    attic/mri_segment_wm_damage
    attic/mri_simulate_atrophy
    attic/mri_threshold
    attic/mri_train_autoencoder
    attic/mri_transform_to_COR
    attic/mri_update_gca
    attic/mri_vol2roi
    attic/mri_wmfilter
    attic/mri_xcorr
    attic/mris_add_template
    attic/mris_aseg_distance
    attic/mris_average_parcellation
    attic/mris_classify_thickness
    attic/mris_compute_optimal_kernel
    attic/mris_density
    attic/mris_fbirn_annot
    attic/mris_glm
    attic/mris_longitudinal_surfaces
    attic/mris_make_map_surfaces
    attic/mris_ms_surface_CNR
    attic/mris_surface_change
    attic/mris_svm_classify
    attic/mris_svm_train
    attic/mris_twoclass
    attic/mris_w_to_curv
    attic/oct_train
    attic/rbftest
    attic/regdat2xfm
    attic/repair_siemens_file
    attic/spherical_stats
  )
endif()

# These never built, i.e. have no CMakeLists.txt, so they are not listed above in attic list
#
# attic/mri_cc_ma_fill
# attic/mri_cc_medial_axis
# attic/mri_extract_ma_intensity
# attic/mri_map_atrophy
# attic/mris_curvature_vectors
# attic/mris_pval_fill

# This subdirectory does the install step to add a script named after
# binaries in the attic when BUILD_ATTIC is false (the default). The
# script is a placeholder to notify the user the command is no longer
# built (and echoes that). If this is the last subdir added, then this
# should be the last step in the install pass.
add_subdirectory(${CMAKE_SOURCE_DIR}/attic/ATTIC_POSTINSTALL)

