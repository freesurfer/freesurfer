# WARNING: this is in development and should not yet be used to replace automake

# before starting the project, set gcc as the default compiler
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project(freesurfer)
cmake_minimum_required(VERSION 3.5)

set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)


# --------------------------------------------------
#                 operating system
# --------------------------------------------------

add_definitions(-D${CMAKE_SYSTEM_NAME})

# not sure what this ANSI is used for
# todo: look into removing it
add_definitions(-DANSI)

# todo: logic for SSE matrix and math functions (affine.h and sse_mathfun.h)
add_definitions(-DUSE_SSE_MATHFUN)
set(COMPILE_FLAGS "${COMPILE_FLAGS} -msse2 -mfpmath=sse")

# check for big endian
include(TestBigEndian)
test_big_endian(IS_BIG_ENDIAN)

# todo: -DHAVE_CONFIG_H


enable_language(Fortran)

# --------------------------------------------------
#                  TMP functions
# --------------------------------------------------

# maybe should move this somewhere else
include(CheckFunctionExists)
include(CheckIncludeFiles)

function(add_help target helptext)
  add_custom_command(COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} &&
                             xxd -i ${helptext} ${CMAKE_CURRENT_BINARY_DIR}/${helptext}.h
                     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${helptext}.h
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${helptext})
  include_directories(${CMAKE_CURRENT_BINARY_DIR})
  target_sources(${target} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/${helptext}.h)
endfunction(add_help)


# todo: is target included in ARGN???
function(add_external_dependencies target)
  foreach(i ${ARGN})
    if(TARGET ${i})
      add_dependencies(${target} ${i})
    endif()
  endforeach()
endfunction()

include(cmake/functions.cmake)

# --------------------------------------------------
#                      options
# --------------------------------------------------

# set the compiler warnings
set(COMPILE_FLAGS "${COMPILE_FLAGS} -Wall -Wno-unused-but-set-variable -Wno-unused-result -Wno-unused-local-typedefs")

# -g compiles debugging info for use by debuggers like gdb. Normally we want
# that in the build, but here we add the DEBUGGING option, which greatly
# reduces the size of compiled binaries when turned OFF
option(DEBUGGING "Compile with the -g flag" ON)
if(DEBUGGING)
  set(COMPILE_FLAGS "${COMPILE_FLAGS} -g")
endif()

# strip unused code
if(APPLE)
  set(COMPILE_FLAGS "${COMPILE_FLAGS} -dead_strip")
else()
  set(COMPILE_FLAGS "${COMPILE_FLAGS} -fdata-sections -ffunction-sections -Wl,--gc-sections")
endif()

# optimization
# todo: there needs to be processor-based logic for this
set(COMPILE_FLAGS "${COMPILE_FLAGS} -O3 -m64")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILE_FLAGS}")

# link map
if(NOT APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-Map,ld_map.txt -Wl,--no-demangle")
endif()

# --------------------------------------------------
#                     libraries
# --------------------------------------------------

# -------- zlib --------
find_package(ZLIB REQUIRED)
add_definitions(-DHAVE_ZLIB)

# -------- gfortran --------
find_library(GFORTRAN_LIBRARIES NAMES libgfortran.a)
find_library(BLAS_LIBRARIES NAMES libblas.a)
find_library(LAPACK_LIBRARIES NAMES liblapack.a)
message(STATUS "GFORTRAN_LIBRARIES: ${GFORTRAN_LIBRARIES}")
message(STATUS "BLAS_LIBRARIES: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")


# --------------------------------------------------
#                third-party packages
# --------------------------------------------------

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# -------- expat --------
find_package(Expat REQUIRED)

# -------- jpeg --------
find_package(JPEG REQUIRED)

# -------- tiff --------
find_package(TIFF REQUIRED)

# -------- xml2 --------
find_package(XML2 REQUIRED)

# -------- netcdf --------
find_package(NetCDF REQUIRED)

# -------- minc --------
find_package(MINC REQUIRED)
add_definitions(-DUSE_LOCAL_MINC)

# -------- itk --------
find_package(ITK REQUIRED)
add_definitions(-DHAVE_ITK_LIBS)

# -------- petsc --------
find_package(PETSC REQUIRED)

# -------- ann --------
find_package(ANN REQUIRED)

# -------- tetgen --------
find_package(Tetgen REQUIRED)

# -------- vtk --------
find_package(VTK REQUIRED)

# -------- openMP --------
find_package(OpenMP)
if(OPENMP_FOUND)
  add_definitions(-DHAVE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# -------- openGL --------
find_package(OpenGL)
if(OPENGL_FOUND)
  add_definitions(-DHAVE_OPENGL)
endif()

# -------- glut --------
# todo fix the include dirs in utils/fix the source so that include "GL/glut.h"
find_package(GLUT)
if(GLUT_FOUND)
  add_definitions(-DHAVE_GLUT_LIBS)
endif()

# -------- gts --------
find_package(GLIB2)

# -------- gts --------
find_package(GTS)

# -------- boost program options --------
find_package(Boost COMPONENTS program_options)

# --------------------------------------------------
#                  freesurfer libs
# --------------------------------------------------

# unix
add_subdirectory(unix)

# log
add_subdirectory(log)

# rgb
add_subdirectory(rgb)

# hipsstubs
add_subdirectory(hipsstubs)

# dicom
add_subdirectory(dicom)

# c utils
add_subdirectory(utils)

# cpp utils
# add_subdirectory(utilscpp)

# fsgdf
add_subdirectory(fsgdf)


# --------------------------------------------------
#                      programs
# --------------------------------------------------

set(FS_SUBDIRS
  # mris_decimate  # uses GTS and crazy glib stuff
  # mris_topo_fixer  # uses the CPP utils
  # qdec
  # qdec_glmfit  # qdec has its own thing goin on
  # tkmedit
  # tkregister2
  # tksurfer
  BrainstemSS
  check_siemens_dir
  connectgraph
  dmri_tensoreig
  dummy
  fem_elastic
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
  lta_convert
  mkxsubjreg
  matlab
  mri_add_xform_to_header
  mri_and
  mri_annotation2label
  mri_aparc2aseg
  mri_aparc2wmseg
  mri_apply_inu_correction
  mri_aseg_edit_reclassify
  mri_aseg_edit_train
  mri_auto_fill
  mri_average
  mri_bc_sc_bias_correct
  mri_bias
  mri_binarize
  mri_build_priors
  mri_ca_label
  mri_ca_normalize
  mri_ca_register
  mri_ca_tissue_parms
  mri_ca_train
  mri_cal_renormalize_gca
  mri_cc
  mri_cht2p
  mri_classify
  mri_cnr
  mri_compile_edits
  mri_compute_change_map
  mri_compute_overlap
  mri_compute_seg_overlap
  mri_compute_structure_transforms
  mri_compute_volume_fractions
  mri_concat
  mri_concatenate_gcam
  mri_concatenate_lta
  mri_convert
  mri_convert_mdh
  mri_copy_values
  mri_cor2label
  mri_coreg
  mri_correct_segmentations
  mri_cvs_register
  mri_dct_align
  mri_deface
  mri_diff
  mri_distance_transform
  mri_divide_segmentation
  mri_edit_segmentation
  mri_edit_segmentation_with_surfaces
  mri_edit_wm_with_aseg
  mri_elastic_energy
  mri_em_register
  mri_estimate_tissue_parms
  mri_evaluate_morph
  mri_extract
  mri_extract_conditions
  mri_extract_fcd_features
  mri_extract_label
  mri_extract_largest_CC
  mri_fcili
  mri_fdr
  mri_fieldsign
  mri_fill
  mri_fit_bias
  mri_fslmat_to_lta
  mri_fuse_intensity_images
  mri_fuse_segmentations
  mri_fwhm
  mri_gca_ambiguous
  mri_gcab_train
  mri_gcut
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
  mri_info
  mri_interpolate
  mri_jacobian
  mri_joint_density
  mri_label2label
  mri_label2vol
  mri_label_accuracy
  mri_label_histo
  mri_label_vals
  mri_label_volume
  mri_linear_register
  mri_log_likelihood
  mri_make_bem_surfaces
  mri_make_density_map
  mri_make_labels
  mri_make_register
  mri_make_template
  mri_map_cpdat
  mri_mark_temporal_lobe
  mri_mask
  mri_matrix_multiply
  mri_mc
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
  mri_normalize
  mri_normalize_tp2
  mri_otl
  mri_paint
  mri_parse_sdcmdir
  mri_parselabel
  mri_partial_ribbon
  mri_path2label
  mri_polv
  mri_probe_ima
  mri_probedicom
  mri_reduce
  mri_register
  mri_relabel_hypointensities
  mri_relabel_nonwm_hypos
  mri_remove_neck
  mri_rf_label
  mri_rf_long_label
  mri_rf_long_train
  mri_rf_train
  mri_ribbon
  mri_rigid_register
  mri_robust_register
  mri_sbbr
  mri_seg_diff
  mri_segcentroids
  mri_seghead
  mri_segment
  mri_segment_tumor
  mri_segment_wm_damage
  mri_segreg
  mri_segstats
  mri_simulate_atrophy
  mri_stats2seg
  mri_strip_nonwhite
  mri_strip_subject_info
  mri_surf2surf
  mri_surf2vol
  mri_surfacemask
  mri_surfcluster
  mri_synthesize
  mri_tessellate
  mri_threshold
  mri_topologycorrection
  mri_train
  mri_train_autoencoder
  mri_transform
  mri_transform_to_COR
  mri_twoclass
  mri_update_gca
  mri_vol2roi
  mri_vol2surf
  mri_vol2vol
  mri_volcluster
  mri_voldiff
  mri_volsynth
  mri_warp_concat
  mri_warp_convert
  mri_watershed
  mri_wbc
  mri_wmfilter
  mri_xcorr
  mri_xvolavg
  mri_z2p
  mris2rgb
  mris_add_template
  mris_anatomical_stats
  mris_annot_to_segmentation
  mris_apply_reg
  mris_aseg_distance
  mris_average_curvature
  mris_average_parcellation
  mris_BA_segment
  mris_ca_label
  mris_ca_train
  mris_calc
  mris_classify_thickness
  mris_compute_acorr
  mris_compute_lgi
  mris_compute_optimal_kernel
  mris_compute_overlap
  mris_compute_parc_overlap
  mris_compute_volume_fractions
  mris_congeal
  mris_convert
  mris_copy_header
  mris_curvature
  mris_curvature_stats
  mris_deform
  mris_density
  mris_diff
  mris_distance_map
  mris_distance_to_label
  mris_distance_transform
  mris_divide_parcellation
  mris_entropy
  mris_errors
  mris_euler_number
  mris_expand
  mris_extract_patches
  mris_fbirn_annot
  mris_fill
  mris_find_flat_regions
  mris_fix_topology
  mris_flatten
  mris_fwhm
  mris_glm
  mris_hausdorff_dist
  mris_inflate
  mris_info
  mris_init_global_tractography
  mris_interpolate_warp
  mris_jacobian
  mris_label2annot
  mris_label_area
  mris_label_calc
  mris_label_mode
  mris_left_right_register
  mris_longitudinal_surfaces
  mris_make_average_surface
  mris_make_face_parcellation
  mris_make_map_surfaces
  mris_make_surfaces
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
  # mris_pmake
  mris_register
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
  mris_smooth
  mris_sphere
  mris_spherical_average
  mris_surface_change
  mris_surface_stats
  mris_surface_to_vol_distances
  mris_svm_classify
  mris_svm_train
  mris_talairach
  mris_thickness
  mris_thickness_comparison
  mris_thickness_diff
  mris_transform
  mris_translate_annotation
  mris_transmantle_dysplasia_paths
  mris_twoclass
  mris_volmask
  mris_volume
  mris_w_to_curv
  mris_warp
  mris_watershed
  mrisp_paint
  mrisp_write
  oct_register_mosaic
  oct_train
  optseq2
  rbftest
  regdat2xfm
  repair_siemens_file
  scripts
  spherical_stats
  spline3
  stat_normalize
  stim_polar
  swi_processing
  talairach_afd
  talairach_avi
  template
  test_makevol
  ThalamicNuclei
  trc
  tridec
)

foreach(SUBDIR ${FS_SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()


# --------------------------------------------------
#                      TMP to be named
# --------------------------------------------------

add_subdirectory(distribution)
