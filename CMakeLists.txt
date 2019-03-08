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
option(INFANT_MODULE "Include infant recon-all" OFF)
option(SUPPRESS_WARNINGS "Suppress some selected warnings" ON)

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

# create build stamp during install
set(BUILD_STAMP "freesurfer" CACHE STRING "Distribution build-stamp")
install(CODE "file(WRITE ${CMAKE_INSTALL_PREFIX}/build-stamp.txt ${BUILD_STAMP}\\n)")

# if an install prefix is not provided, check the FS_INSTALL_DIR env var, and
# if that is not defined, set the default path to /usr/local/freesurfer
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
    set(PYTHON_LIBRARY ${FSPYTHON_DIR}/lib/libpython3.6m.a)
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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wno-sign-compare -Wno-write-strings -Wno-unused-result")

# clang complains about -Wno-unused-but-set-variable and says to use -Wno-unused-const-variable
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-const-variable")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-but-set-variable")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# linker options
if(APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -dead_strip -static-libgcc")
else()
  set(STRIP_FLAGS "-fdata-sections -ffunction-sections -Wl,--gc-sections")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STRIP_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${STRIP_FLAGS} -Wl,-Map,ld_map.txt -Wl,--no-demangle")
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
    check_siemens_dir
    connectgraph
    diffusion_tool
    dmri_tensoreig
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
    trc
    tridec
  )
endif()

# infant freesurfer
if(INFANT_MODULE)
  add_subdirectories(infant)
endif()

# tk tools
if(TKTOOLS_MODULE)
  add_subdirectory(tktools)
endif()
