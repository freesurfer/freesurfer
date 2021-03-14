/*
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef UTILS_AFNI_HPP
#define UTILS_AFNI_HPP

#include "mri.h"

namespace fs::utils {

struct AFNI_HEADER {
  // MARK: Mandatory Mttributes

  int dataset_rank[2];
  // Two values that determine the dimensionality of the dataset:
  //     [0] = Number of spatial dimensions (must be 3).
  //     [1] = Number of sub-bricks in the dataset (in most programs, this is
  // called "nvals").
  // At one time I thought I might extend AFNI to support n-dimensional
  // datasets, but as time went one, I decided to support the fourth dimension
  // not by increasing the "rank" of a dataset, but by adding the time axis
  // instead. Thus, the dataset rank is always set to 3.

  int dataset_dimensions[3];
  // Three values that determine the size of each spatial axis of the dataset:
  //                      [0] = number of voxels along the x-axis (nx).
  //                      [1] = number of voxels along the y-axis (ny).
  //                      [2] = number of voxels along the z-axis (nz).
  // The voxel with 3-index (i,j,k) in a sub-brick is located at position
  // (i+j*nx+k*nx*ny), for i=0..nx-1, j=0..ny-1, k=0..nz-1.  Each axis must have
  // at least 2 points!

  char typestring[16];
  // One of "3DIM_HEAD_ANAT" or "3DIM_HEAD_FUNC" or "3DIM_GEN_ANAT"  or
  // "3DIM_GEN_FUNC". Determines if the dataset is of Anat or Func type
  // (grayscale underlay or color overlay).  If Anat type, and if it is a _HEAD_
  // dataset in the +orig view, then Talairach markers might be attached to it
  // (if it was created by to3d).

  int scene_data[3];
  // Three integer codes describing the dataset type
  // [0] = view type: 0=+orig, 1=+acpc, 2=+tlrc
  // [1] = func type:
  // If dataset is Anat type, then this is one of the following codes:
  //                      #define ANAT_SPGR_TYPE   0
  //                      #define ANAT_FSE_TYPE    1
  //                      #define ANAT_EPI_TYPE    2
  //                      #define ANAT_MRAN_TYPE   3
  //                      #define ANAT_CT_TYPE     4
  //                      #define ANAT_SPECT_TYPE  5
  //                      #define ANAT_PET_TYPE    6
  //                      #define ANAT_MRA_TYPE    7
  //                      #define ANAT_BMAP_TYPE   8
  //                      #define ANAT_DIFF_TYPE   9
  //                      #define ANAT_OMRI_TYPE   10
  //                      #define ANAT_BUCK_TYPE   11
  //
  // At this time, Anat codes 0..10 are treated identically by all AFNI
  // programs.  Code 11 marks the dataset as a "bucket" type, which is treated
  // differently in the display; the "Define Overlay" control panel will have a
  // chooser that allows you to specify which sub-brick from the bucket should
  // be used to make the underlay image.
  //
  // If dataset is Func type, then this is one of the following codes (Please
  // modify @statauxcode if you make additions or changes here):
  //
  //                      #define FUNC_FIM_TYPE   0  /* 1 value           */
  //                      #define FUNC_THR_TYPE   1  /* obsolete          */
  //                      #define FUNC_COR_TYPE   2  /* fico: correlation */
  //                      #define FUNC_TT_TYPE    3  /* fitt: t-statistic */
  //                      #define FUNC_FT_TYPE    4  /* fift: F-statistic */
  //                      #define FUNC_ZT_TYPE    5  /* fizt: z-score     */
  //                      #define FUNC_CT_TYPE    6  /* fict: Chi squared */
  //                      #define FUNC_BT_TYPE    7  /* fibt: Beta stat   */
  //                      #define FUNC_BN_TYPE    8  /* fibn: Binomial    */
  //                      #define FUNC_GT_TYPE    9  /* figt: Gamma       */
  //                      #define FUNC_PT_TYPE    10 /* fipt: Poisson     */
  //                      #define FUNC_BUCK_TYPE  11 /* fbuc: bucket      */
  //
  // These types are defined more fully in README.func_types.
  //
  // Unfortunately, the func type codes overlap for Func and Anat datasets. This
  // means that one cannot tell the contents of a dataset from a single
  // attribute. However, this bad design choice (from 1994) is now enshrined in
  // the .HEAD files of thousands of datasets, so it will be hard to change.
  //
  // [2] = 0 or 1 or 2 or 3, corresponding to the TYPESTRING values given above.
  // If this value does not match the typestring value, then the dataset is
  // malformed and AFNI will reject it!

  int orient_specific[3];
  // Three integer codes describing the spatial orientation of the dataset axes;
  // [0] for the x-axis, [1] for the y-axis, and [2] for the z-axis.  The
  // possible codes are:
  //
  //                    #define ORI_R2L_TYPE  0  /* Right to Left         */
  //                    #define ORI_L2R_TYPE  1  /* Left to Right         */
  //                    #define ORI_P2A_TYPE  2  /* Posterior to Anterior */
  //                    #define ORI_A2P_TYPE  3  /* Anterior to Posterior */
  //                    #define ORI_I2S_TYPE  4  /* Inferior to Superior  */
  //                    #define ORI_S2I_TYPE  5  /* Superior to Inferior  */
  //
  // Note that these codes must make sense (e.g., they can't all be 4).  Only
  // program to3d enforces this restriction, but if you create a nonsensical
  // dataset, then bad things will happen at some point.
  //
  // Spatial xyz-coordinates in AFNI are sometimes used in dataset order, which
  // refers to the order given here. They are also sometimes used in Dicom
  // order, in which x=R-L, y=A-P, and z=I-S (R,A,I are < 0; L,P,S are > 0).
  // There are utility functions for converting dataset ordered 3-vectors to and
  // from Dicom ordered 3-vectors -- see the functions in file thd_coords.c.
  // Distances in AFNI are always encoded in millimeters.

  float origin[3];
  // Three numbers giving the xyz-coordinates of the center of the (0,0,0) voxel
  // in the dataset.  The order of these numbers is the same as the order of the
  // xyz-axes (cf. ORIENT_SPECIFIC). However, the AFNI convention is that R-L,
  // A-P, and I-S are negative-to-positive.  Thus, if the y-axis is P-A (say),
  // then the y-origin is likely to be positive (and the y-delta, below, would
  // be negative).  These numbers are usually computed from the centering
  // controls in to3d.

  float delta[3];
  //  Three numbers giving the (x,y,z) voxel sizes, in the same order as
  //  ORIENT_SPECIFIC.  That is, [0] = x-delta, [1] = y-delta, and [2] =
  //  z-delta.  These values may be negative; in the example above, where the
  //  y-axis is P-A, then y-delta would be negative. The center of the (i,j,k)
  //  voxel is located at xyz-coordinates ORIGIN[0]+i*DELTA[0],
  //  ORIGIN[1]+j*DELTA[1], ORIGIN[2]+k*DELTA[2]

  // MARK: Time-Dependent Dataset Attributes
  // These attributes are mandatory if the .HEAD file describes a 3D+time
  // dataset.

  int taxis_nums[3];
  // [0] = Number of points in time (at present, must be equal to
  // nvals=DATASET_RANK[1], or AFNI programs will not be happy; that is, each
  // time point can only have a single numerical value per voxel).
  //
  // [1] = Number of slices with time offsets.  If zero, then no slice-dependent
  // time offsets are present (all slices are presumed to be acquired at the
  // same time).  If positive, specifies the number of values to read from
  // TAXIS_OFFSETS.  Normally, this would either be 0 or be equal to
  // DATASET_DIMENSIONS[2].
  //
  // [2] = Units codes for TAXIS_FLOATS[1]; one of the following
  //                     #define UNITS_MSEC_TYPE  77001  /* don't ask me */
  //                     #define UNITS_SEC_TYPE   77002  /* where these */
  //                     #define UNITS_HZ_TYPE    77003  /* came from! */

  float taxis_floats[5];
  // [0] = Time origin (in units given by TAXIS_NUMS[2]). This is 0 in datasets
  // created by to3d (at present).
  //
  // [1] = Time step (TR).
  //
  // [2] = Duration of acquisition.  This is 0 in datasets created by to3d (at
  // present)
  //
  // [3] = If TAXIS_NUMS[1] > 0, then this is the z-axis offset for the
  // slice-dependent time offsets.  This will be equal to ORIGIN[2] in datasets
  // created by to3d.c.
  //
  // [4] = If TAXIS_NUMS[1] > 0, then this is the z-axis step for the
  // slice-dependent time offsets.  This will be equal to DELTA[2] in datasets
  // created by to3d.c.

  float *taxis_offsets;
  // If TAXIS_NUMS[1] > 0, then this array gives the time offsets of the slices
  // defined by TAXIS_FLOATS[3..4].
  //
  // The time offset at z = TAXIS_FLOATS[3] + k*TAXIS_FLOATS[4] is
  // TAXIS_OFFSETS[k], for k=0..TAXIS_NUMS[1]-1.
  //
  // If TAXIS_NUMS[1] == 0, the this attribute is not used.

  // MARK: Almost Mandatory Attributes
  // The following useful attributes are present in most AFNI datasets created
  // by AFNI package programs.  However, if they are not present, then the
  // function that assembles a dataset struct will get by.
  char * idcode_string;
  char * idcode_date;
  char   byteorder_string[10];
  float *brick_stats;
  int *  brick_types;
  float *brick_float_facs;
  char * brick_labs;
  float *brick_statusux;
  float *statu_aux;

  int numchars;
  int numstats;
  int numtypes;
  int numfacs;

  // MARK: Notes Attributes
  char * history_note;
  int    notes_count;
  char **notes; // NOTE_NUMBER_001 to NOTE_NUMBER_999

  // MARK: Registration Attributes
  float * tagalign_matved;
  float **volreg_matved; // for sub-brick 001 to 999
  char ** volreg_rotcom; // for sub-brick 001 to 999
  float * volreg_center_old;
  float * volreg_center_base;
  char *  volreg_rotparent_idcode;
  char *  volreg_rotparent_name;
  char *  volreg_gridparend_idcode;
  char *  volreg_gridparent_name;
  char *  volreg_input_idcode;
  char *  volreg_input_name;
  char *  volreg_base_idcode;
  char *  volreg_base_name;
  int     volreg_rotcom_num;

  // MARK: Miscellaneous Attributes

  char *idcode_anat_parent;
  int * to3d_zpad;

  // MARK: Warping
  char * idcode_warp_parent;
  int    warp_type[2];
  float *warp_data;

  // MARK: Talairach Markers Attributes
  float *marks_xyz;
  char * marks_lab;
  char * marks_help;
  int    marks_flags[2];

  // MARK: Attributes for User-Defined Tags
  int    tagset_num[2];
  float *tagset_floats;
  char * tagset_labels;

  // MARK: Nearly Useless Attributes

  char *label_1;
  char *label_2;
  char *dataset_name;
  char *dataset_keywords;
  char *brick_keywords;
};

using AF = struct AFNI_HEADER;

auto afniRead(const char *fname, int read_volume) -> MRI *;
auto afniWrite(MRI *mri, const char *fname) -> int;
auto readAFNIHeader(FILE *fp, AF *paf) -> int;
void AFinit(AF &pAF);
void AFclean(AF *pAF);
void printAFNIHeader(AF &pAF);

} // namespace fs::utils

#endif
