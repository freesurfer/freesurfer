/**
 * @brief utilities for reading/writing MRI data structure
 *
 * Reading and writing most of the major MRI volume types, to and from
 * the Freesurfer .mgz format, is provided here.
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/
#define _MRIIO_SRC

#include <sstream>
#include <iomanip>
#include <vector>

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "machine.h"
#include "mri.h"
#include "AFNI.h"
#include "Bruker.h"
#include "DICOMRead.h"
#include "analyze.h"
#include "autoencoder.h"
#include "bfileio.h"
#include "chklc.h"
#include "cma.h"
#include "diag.h"
#include "dti.h"
#include "error.h"
#include "fio.h"
#include "gcamorph.h"
#include "gifti.h"
#include "imautils.h"
#include "macros.h"
#include "matfile.h"
#include "math.h"
#include "matrix.h"
#include "mghendian.h"
#include "mri2.h"
#include "mri_circulars.h"
#include "mri_identify.h"
#include "nifti1.h"
#include "nifti1_io.h"
#include "proto.h"
#include "region.h"
#include "signa.h"
#include "tags.h"
#include "utils.h"
#include "znzlib.h"
#include "romp_support.h"
#include "NrrdIO.h"

#include "nii_dicom.h"
#include "nii_dicom_batch.h"

static int niiPrintHdr(FILE *fp, struct nifti_1_header *hdr);

// unix director separator
#define DIR_SEPARATOR '/'
#define CURDIR "./"

#define MM_PER_METER 1000.0f
#define INFO_FNAME "COR-.info"


MRI *mri_read(const char *fname, int type, int volume_flag, int start_frame, int end_frame);
static MRI *corRead(const char *fname, int read_volume);
static int corWrite(MRI *mri, const char *fname);
static MRI *siemensRead(const char *fname, int read_volume);
static MRI *readGCA(const char *fname, int start_frame, int end_frame);

static MRI *mincRead(const char *fname, int read_volume);
static int mincWrite(MRI *mri, const char *fname);

static int bvolumeWrite(MRI *vol, const char *fname_passed, int type);
// static int bshortWrite(MRI *mri, const char *fname_passed);
// static int bfloatWrite(MRI *mri, const char *stem);
static int write_bhdr(MRI *mri, FILE *fp);
static int read_bhdr(MRI *mri, FILE *fp);

static MRI *bvolumeRead(const char *fname_passed, int read_volume, int type);
// static MRI *bshortRead(const char *fname_passed, int read_volume);
// static MRI *bfloatRead(const char *fname_passed, int read_volume);
static MRI *genesisRead(const char *stem, int read_volume);
static MRI *gelxRead(const char *stem, int read_volume);

static int CountAnalyzeFiles(const char *analyzefname, int nzpad, char **ppstem);
static MRI *analyzeRead(const char *fname, int read_volume);
static dsr *ReadAnalyzeHeader(const char *hdrfile, int *swap, int *mritype, int *bytes_per_voxel);
static int DumpAnalyzeHeader(FILE *fp, dsr *hdr);

static int analyzeWrite(MRI *mri, const char *fname);
static int analyzeWriteFrame(MRI *mri, const char *fname, int frame);
static int analyzeWriteSeries(MRI *mri, const char *fname);
static int analyzeWrite4D(MRI *mri, const char *fname);

static void swap_analyze_header(dsr *hdr);

static int nan_inf_check(MRI *mri);
#ifdef VC_TO_CV
static int voxel_center_to_center_voxel(MRI *mri, float *x, float *y, float *z);
#endif
static MRI *gdfRead(const char *fname, int read_volume);
static int gdfWrite(MRI *mri, const char *fname);

static MRI *ximgRead(const char *fname, int read_volume);

static MRI *nifti1Read(const char *fname, int read_volume);
static int nifti1Write(MRI *mri, const char *fname);
static MRI *niiRead(const char *fname, int read_volume);
static MRI *niiReadFromMriFsStruct(MRIFSSTRUCT *mrifsStruct);
static int niiWrite(MRI *mri, const char *fname);
static int itkMorphWrite(MRI *mri, const char *fname);
static int niftiQformToMri(MRI *mri, struct nifti_1_header *hdr);
static int mriToNiftiQform(MRI *mri, struct nifti_1_header *hdr);
static int mriToNiftiSform(MRI *mri, struct nifti_1_header *hdr);
static int niftiSformToMri(MRI *mri, struct nifti_1_header *hdr);
static void swap_nifti_1_header(struct nifti_1_header *hdr);
static MRI *MRISreadCurvAsMRI(const char *curvfile, int read_volume);

static MRI *mriNrrdRead(const char *fname, int read_volume);
static int mriNrrdWrite(MRI *mri, const char *fname);

/********************************************/

static void short_local_buffer_to_image(short *buf, MRI *mri, int slice, int frame);

static void int_local_buffer_to_image(int *buf, MRI *mri, int slice, int frame);

static void long32_local_buffer_to_image(long32 *buf, MRI *mri, int slice, int frame);

static void float_local_buffer_to_image(float *buf, MRI *mri, int slice, int frame);

static void local_buffer_to_image(BUFTYPE *buf, MRI *mri, int slice, int frame);

static MRI *sdtRead(const char *fname, int read_volume);
static MRI *mghRead(const char *fname, int read_volume, int frame);
static int mghWrite(MRI *mri, const char *fname, int frame);
static int mghAppend(MRI *mri, const char *fname, int frame);

/********************************************/

extern int errno;

extern const char *Progname;

static char *command_line;
static char *subject_name;
static int gdf_crop_flag = FALSE;

#define MAX_UNKNOWN_LABELS 100

static short cma_field[512][512];
static char unknown_labels[MAX_UNKNOWN_LABELS][STRLEN];
static int n_unknown_labels;

//////////////////////////////////////////////////////////////////////
// this is a one way of setting direction cosine
// when the direction cosine is not provided in the volume.
// may not agree with the volume. what can we do?  Let them set by themselves.
int setDirectionCosine(MRI *mri, int orientation)
{
  switch (orientation) {
    case MRI_CORONAL:  // x is from right to left.
      // y is from top to neck, z is from back to front
      mri->x_r = -1;
      mri->y_r = 0;
      mri->z_r = 0;
      mri->c_r = 0;
      mri->x_a = 0;
      mri->y_a = 0;
      mri->z_a = 1;
      mri->c_a = 0;
      mri->x_s = 0;
      mri->y_s = -1;
      mri->z_s = 0;
      mri->c_s = 0;
      break;
    case MRI_SAGITTAL:  // x is frp, back to front,
      // y is from top to neck, z is from left to right
      mri->x_r = 0;
      mri->y_r = 0;
      mri->z_r = 1;
      mri->c_r = 0;
      mri->x_a = 1;
      mri->y_a = 0;
      mri->z_a = 0;
      mri->c_a = 0;
      mri->x_s = 0;
      mri->y_s = -1;
      mri->z_s = 0;
      mri->c_s = 0;
      break;
    case MRI_HORIZONTAL:  // x is from right to left,
      // y is from front to back, z is from neck to top
      mri->x_r = -1;
      mri->y_r = 0;
      mri->z_r = 0;
      mri->c_r = 0;
      mri->x_a = 0;
      mri->y_a = -1;
      mri->z_a = 0;
      mri->c_a = 0;
      mri->x_s = 0;
      mri->y_s = 0;
      mri->z_s = 1;
      mri->c_s = 0;
      break;
    default:
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "setDirectionCosine():unknown slice direction"));
      break;  // should not reach here (handled at the conversion)
  }
  return (NO_ERROR);
}

#define isOne(a) (fabs(fabs(a) - 1) < 0.00001)
#define isCloseToOne(a) (fabs(fabs(a) - 1) < 0.1)

// here I take the narrow view of slice_direction
int getSliceDirection(MRI *mri)
{
  int direction = MRI_UNDEFINED;

  if (!strcmp(MRIsliceDirectionName(mri), "coronal"))
    direction = MRI_CORONAL;
  else if (!strcmp(MRIsliceDirectionName(mri), "sagittal"))
    direction = MRI_SAGITTAL;
  else if (!strcmp(MRIsliceDirectionName(mri), "axial"))
    direction = MRI_AXIAL;

  /*  if (isCloseToOne(mri->x_r) &&
      isCloseToOne(mri->y_s) &&
      isCloseToOne(mri->z_a))
    direction = MRI_CORONAL;
  else if (isCloseToOne(mri->x_a) &&
           isCloseToOne(mri->y_s) &&
           isCloseToOne(mri->z_r))
    direction = MRI_SAGITTAL;
  else if (isCloseToOne(mri->x_r) &&
           isCloseToOne(mri->y_a) &&
           isCloseToOne( mri->z_s))
           direction = MRI_HORIZONTAL;*/
  return direction;
}

// For surface, we currently cannot handle volumes with general slice direction
// nor we cannot handle non-conformed volumes
int mriOKforSurface(MRI *mri)
{
  // first check slice direction
  if (getSliceDirection(mri) != MRI_CORONAL) return 0;
  // remove slice size limitation
  // else if (mri->width != 256 || mri->height != 256 || mri->depth != 256)
  //  return 0;
  // remove check for 1 mm size
  //  else if (mri->xsize != 1 || mri->ysize != 1 || mri->zsize != 1)
  //  return 0;
  else
    return 1;
}

int mriConformed(MRI *mri)
{
  // first check slice direction
  if (getSliceDirection(mri) != MRI_CORONAL)
    return 0;
  else if (mri->width != 256 || mri->height != 256 || mri->depth != 256)
    return 0;
  // check to four decimals of precision only
  else if (nint(mri->xsize * 10000) != (int)10000 || nint(mri->ysize * 10000) != (int)10000 ||
           nint(mri->zsize * 10000) != (int)10000)
    return 0;
  else
    return 1;
}

float MRIfindMinSize(MRI *mri, int *conform_width)
{
  double xsize, ysize, zsize, minsize;
  double fwidth, fheight, fdepth, fmax;
  xsize = mri->xsize;
  ysize = mri->ysize;
  zsize = mri->zsize;
  // there are 3! = 6 ways of ordering
  //             xy  yz  zx
  // x > y > z    z min
  // x > z > y    y min
  // z > x > y    y min
  //////////////////////////
  // y > x > z    z min
  // y > z > x    x min
  // z > y > x    x min
  if (xsize > ysize) {
    minsize = (ysize > zsize) ? zsize : ysize;
  }
  else {
    minsize = (zsize > xsize) ? xsize : zsize;
  }

  // now decide the conformed_width
  // algorighm ----------------------------------------------
  // calculate the size in mm for all three directions
  fwidth = mri->xsize * mri->width;
  fheight = mri->ysize * mri->height;
  fdepth = mri->zsize * mri->depth;
  // pick the largest
  if (fwidth > fheight) {
    fmax = (fwidth > fdepth) ? fwidth : fdepth;
  }
  else {
    fmax = (fdepth > fheight) ? fdepth : fheight;
  }

  *conform_width = (int)ceil(nint(fmax / minsize * 10000) / 10000);
  // used to be: *conform_width = (int) ceil(fmax/minsize);
  // just to make sure that if smaller than 256, use 256 anyway
  if (*conform_width < 256) {
    *conform_width = 256;
  }

  return (float)minsize;
}

// this function is called when conform is done
int MRIfindRightSize(MRI *mri, float conform_size)
{
  // user gave the conform_size
  double xsize, ysize, zsize;
  double fwidth, fheight, fdepth, fmax;
  int conform_width;

  xsize = mri->xsize;
  ysize = mri->ysize;
  zsize = mri->zsize;

  // now decide the conformed_width
  // calculate the size in mm for all three directions
  fwidth = mri->xsize * mri->width;
  fheight = mri->ysize * mri->height;
  fdepth = mri->zsize * mri->depth;
  // pick the largest
  if (fwidth > fheight) {
    fmax = (fwidth > fdepth) ? fwidth : fdepth;
  }
  else {
    fmax = (fdepth > fheight) ? fdepth : fheight;
  }
  // get the width with conform_size
  conform_width = (int)ceil(fmax / conform_size);

  // just to make sure that if smaller than 256, use 256 anyway
  if (conform_width < 256) {
    conform_width = 256;
  }
  // conform_width >= 256.   allow 10% leeway
  else if ((conform_width - 256.) / 256. < 0.1) {
    conform_width = 256;
  }

  // if more than 256, warn users
  if (conform_width > 256) {
    fprintf(stderr,
            "WARNING =================="
            "++++++++++++++++++++++++"
            "=======================================\n");
    fprintf(stderr,
            "The physical sizes are "
            "(%.2f mm, %.2f mm, %.2f mm), "
            "which cannot fit in 256^3 mm^3 volume.\n",
            fwidth,
            fheight,
            fdepth);
    fprintf(stderr, "The resulting volume will have %d slices.\n", conform_width);
    fprintf(stderr,
            "If you find problems, please let us know "
            "(freesurfer@nmr.mgh.harvard.edu).\n");
    fprintf(stderr,
            "=================================================="
            "++++++++++++++++++++++++"
            "===============\n\n");
  }
  return conform_width;
}

void setMRIforSurface(MRI *mri)
{
  if (!mriOKforSurface(mri))
    ErrorExit(ERROR_BADPARM,
              "%s: the volume is not conformed, that is, "
              "the volume must be in CORONAL direction.\n",
              Progname);
}

int mriio_command_line(int argc, char *argv[])
{
  int i;
  int length;
  char *c;

  length = 0;
  for (i = 0; i < argc; i++) length += strlen(argv[i]);

  /* --- space for spaces and \0 --- */
  length += argc;

  command_line = (char *)malloc(length);

  c = command_line;
  for (i = 0; i < argc; i++) {
    strcpy(c, argv[i]);
    c += strlen(argv[i]);
    *c = (i == argc - 1 ? '\0' : ' ');
    c++;
  }

  return (NO_ERROR);

} /* end mriio_command_line() */

void mriio_set_gdf_crop_flag(int new_gdf_crop_flag)
{
  gdf_crop_flag = new_gdf_crop_flag;

  return;

} /* end mriio_set_gdf_crop_flag() */

int MRIgetVolumeName(const char *string, char *name_only)
{
  char *at, *pound;

  strcpy(name_only, string);

  if ((at = strrchr(name_only, '@')) != NULL) *at = '\0';

  if (MRIIO_Strip_Pound) {
    if ((pound = strrchr(name_only, '#')) != NULL) *pound = '\0';
  }

  return (NO_ERROR);

} /* end MRIgetVolumeName() */

MRI *mri_read(const char *fname, int type, int volume_flag, int start_frame, int end_frame)
{
  MRI *mri, *mri2;
  IMAGE *I;
  char fname_copy[STRLEN];
  char *ptmpstr;
  char *at, *pound, *colon;
  char *ep;
  int i, j, k, t;
  int volume_frames;

  // sanity-checks
  if (fname == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "mri_read(): null fname!\n"));
  }
  if (fname[0] == 0) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "mri_read(): empty fname!\n"));
  }

  // if filename does not contain any directory separator, then add cwd
  if (!strchr(fname, DIR_SEPARATOR)) {
    char *cwd = getcwd(NULL, 0);  // posix 1 extension
    // (allocate as much space needed)
    if (cwd) {
      strcpy(fname_copy, cwd);
      strcat(fname_copy, "/");
      strcat(fname_copy, fname);
      free(cwd);
    }
    else  // why fail?
      strcpy(fname_copy, fname);
  }
  else
    strcpy(fname_copy, fname);

  at = strrchr(fname_copy, '@');
  if (MRIIO_Strip_Pound)
    pound = strrchr(fname_copy, '#');
  else
    pound = NULL;

  if (at != NULL) {
    *at = '\0';
    at++;
  }

  if (pound != NULL) {
    *pound = '\0';
    pound++;
  }

  if (at != NULL) {
    type = string_to_type(at);
    if (type == MRI_VOLUME_TYPE_UNKNOWN) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADPARM, "mri_read(): unknown type '%s'\n", at));
    }
  }
  else if (type == MRI_VOLUME_TYPE_UNKNOWN) {
    type = mri_identify(fname_copy);
    if (type == MRI_VOLUME_TYPE_UNKNOWN) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADPARM, "mri_read(): couldn't determine type of file %s", fname_copy));
    }
  }

  if (pound != NULL) {
    colon = strchr(pound, ':');
    if (colon != NULL) {
      *colon = '\0';
      colon++;
      if (*colon == '\0') {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "mri_read(): bad frame specification ('%s:')\n", pound));
      }

      start_frame = strtol(pound, &ep, 10);
      if (*ep != '\0') {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "mri_read(): bad start frame ('%s')\n", pound));
      }

      end_frame = strtol(colon, &ep, 10);
      if (*ep != '\0') {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "mri_read(): bad end frame ('%s')\n", colon));
      }
    }
    else {
      start_frame = end_frame = strtol(pound, &ep, 10);
      if (*ep != '\0') {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "mri_read(): bad frame ('%s')\n", pound));
      }
    }

    if (start_frame < 0) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADPARM, "mri_read(): frame (%d) is less than zero\n", start_frame));
    }

    if (end_frame < 0) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADPARM, "mri_read(): frame (%d) is less than zero\n", end_frame));
    }

    if (start_frame > end_frame) {
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_BADPARM,
                   "mri_read(): the start frame (%d) is "
                   "greater than the end frame (%d)\n",
                   start_frame,
                   end_frame));
    }
  }

  if (type == MRI_CORONAL_SLICE_DIRECTORY) {
    mri = corRead(fname_copy, volume_flag);
  }
  else if (type == SIEMENS_FILE) {
    mri = siemensRead(fname_copy, volume_flag);
  }
  else if (type == MRI_GCA_FILE) {
    mri = readGCA(fname_copy, start_frame, end_frame);
    start_frame = -1;
  }
  else if (type == BHDR) {
    ptmpstr = bhdr_firstslicefname(fname_copy);
    t = bhdr_precision(fname_copy);
    mri = bvolumeRead(ptmpstr, volume_flag, t);
    free(ptmpstr);
  }
  else if (type == BSHORT_FILE) {
    // mri = bshortRead(fname_copy, volume_flag);
    mri = bvolumeRead(fname_copy, volume_flag, MRI_SHORT);
  }
  else if (type == BFLOAT_FILE) {
    // mri = bfloatRead(fname_copy, volume_flag);
    mri = bvolumeRead(fname_copy, volume_flag, MRI_FLOAT);
  }
  else if (type == GENESIS_FILE) {
    mri = genesisRead(fname_copy, volume_flag);
  }
  else if (type == SIGNA_FILE) {
    mri = signaRead(fname_copy, volume_flag);
  }
  else if (type == GE_LX_FILE) {
    mri = gelxRead(fname_copy, volume_flag);
  }
  else if (type == MRI_ANALYZE_FILE || type == MRI_ANALYZE4D_FILE) {
    mri = analyzeRead(fname_copy, volume_flag);
  }
  else if (type == BRIK_FILE) {
    mri = afniRead(fname_copy, volume_flag);
  }
  else if (type == MRI_MINC_FILE) {
    // mri = mincRead2(fname_copy, volume_flag);
    mri = mincRead(fname_copy, volume_flag);
  }
  else if (type == SDT_FILE) {
    mri = sdtRead(fname_copy, volume_flag);
  }
  else if (type == MRI_MGH_FILE) {
    mri = mghRead(fname_copy, volume_flag, -1);
  }
  else if (type == MGH_MORPH) {
    int which = start_frame ;
    GCA_MORPH *gcam;
    gcam = GCAMread(fname_copy);
    if (gcam == NULL) ErrorReturn(NULL, (ERROR_BADPARM, "MRIread(%s): could not read .m3z\n", fname_copy));
    if (gcam->type == GCAM_RAS) GCAMrasToVox(gcam, NULL);
    if (start_frame < 0)
      mri = GCAMwriteWarpToMRI(gcam, NULL);
    else {
      printf("reading 'frame' # %d from gcam (see gcamorph.h for definitions)\n", start_frame);
      if (which == GCAM_NODEX || which == GCAM_NODEY || which == GCAM_NODEZ)
      {
	MRI *mri_template ;

	mri_template = MRIalloc(gcam->image.width, gcam->image.height, gcam->image.depth, MRI_FLOAT);
	MRIcopyVolGeomToMRI(mri_template, &gcam->image);
	MRIsetResolution(mri_template, gcam->image.xsize, gcam->image.ysize, gcam->image.zsize);
	MRIreInitCache(mri_template);
	if (!mri_template)
	  ErrorExit(ERROR_NOMEMORY, "gcamWrite: could not allocate %dx%dx%d MRI\n", gcam->width, gcam->height, gcam->depth);
	
	mri = GCAMmorphFieldFromAtlas(gcam, mri_template, which, 0, 0) ;
	MRIfree(&mri_template) ;
      }
      else
	mri = GCAMwriteMRI(gcam, NULL, start_frame);
      start_frame = end_frame = 0;
    }

    GCAMfree(&gcam);
  }
  else if (type == MGH_AUTOENCODER) {
    SAE *sae;
    sae = SAEread(fname_copy);
    if (sae == NULL) ErrorReturn(NULL, (ERROR_BADPARM, "MRIread(%s): could not read autoencoder\n", fname_copy));

    mri = SAElayerWeightsToMRI(sae, start_frame);
    start_frame = 0;
    end_frame = mri->nframes - 1;
    SAEfree(&sae);
  }
  else if (type == GDF_FILE) {
    mri = gdfRead(fname_copy, volume_flag);
  }
  else if (type == DICOM_FILE) {  
    if (UseDICOMRead3) {
      MRIFSSTRUCT *mrifsStruct = DICOMRead3(fname_copy, volume_flag);
      if (mrifsStruct == NULL) 
	return NULL;

      mri = niiReadFromMriFsStruct(mrifsStruct);
      if (mri != NULL && mri->ti < 0)
	mri->ti = 0;

      free(mrifsStruct->imgM);
      free(mrifsStruct->tdti);
    }
    else if (!UseDICOMRead2){
      DICOMRead(fname_copy, &mri, volume_flag);
      printf("mriio.cpp: starting DICOMRead()\n");
    }
    else{
      mri = DICOMRead2(fname_copy, volume_flag);
      printf("mriio.cpp: starting DICOMRead2()\n");
    }
  }
  else if (type == SIEMENS_DICOM_FILE) {
    if (UseDICOMRead3) {
      printf("mriio.cpp: starting DICOMRead3()\n");
      MRIFSSTRUCT *mrifsStruct = DICOMRead3(fname_copy, volume_flag);
      if (mrifsStruct == NULL) 
	return NULL;
      mri = niiReadFromMriFsStruct(mrifsStruct);
      free(mrifsStruct->imgM);
      free(mrifsStruct->tdti);
    } 
    else {
      printf("mriio.cpp: starting sdcmLoadVolume()\n");
      // mri_convert -nth option sets start_frame = nth.  otherwise -1
      mri = sdcmLoadVolume(fname_copy, volume_flag, start_frame);
      start_frame = -1;
      // in order to avoid the later processing on start_frame and end_frame
      // read the comment later on
      end_frame = 0;
    }
  }
  else if (type == BRUKER_FILE) {
    mri = brukerRead(fname_copy, volume_flag);
  }
  else if (type == XIMG_FILE) {
    mri = ximgRead(fname_copy, volume_flag);
  }
  else if (type == NIFTI1_FILE) {
    mri = nifti1Read(fname_copy, volume_flag);
  }
  else if (type == NII_FILE) {
    mri = niiRead(fname_copy, volume_flag);
  }
  else if (type == NRRD_FILE) {
    mri = mriNrrdRead(fname_copy, volume_flag);
  }
  else if (type == MRI_CURV_FILE)
    mri = MRISreadCurvAsMRI(fname_copy, volume_flag);
  else if (type == GIFTI_FILE)
    mri = MRISreadGiftiAsMRI(fname_copy, volume_flag);
  else if (type == IMAGE_FILE) {
    I = ImageRead(fname_copy);
    mri = ImageToMRI(I);
    ImageFree(&I);
  }
  else {
    fprintf(stderr, "mri_read(): type = %d\n", type);
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "mri_read(): code inconsistency "
                 "(file type recognized but not caught)"));
  }

  if (mri == NULL) return (NULL);
  strcpy(mri->fname, fname);  // added by dng 11/16/2010

  // update/cache the transform
  MATRIX *tmp;
  if (mri->i_to_r__) {
    AffineMatrixFree(&(mri->i_to_r__));
  }
  AffineMatrixAlloc(&(mri->i_to_r__));
  tmp = extract_i_to_r(mri);
  SetAffineMatrix(mri->i_to_r__, tmp);
  MatrixFree(&tmp);

  if (mri->r_to_i__) MatrixFree(&(mri->r_to_i__));
  mri->r_to_i__ = extract_r_to_i(mri);

  /* Compute the FOV from the vox2ras matrix (don't rely on what
     may or may not be in the file).*/

  if (start_frame == -1) return (mri);

  /* --- select frames --- */

  if (start_frame >= mri->nframes) {
    volume_frames = mri->nframes;
    MRIfree(&mri);
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "mri_read(): start frame (%d) is "
                 "out of range (%d frames in volume)",
                 start_frame,
                 volume_frames));
  }

  if (end_frame >= mri->nframes) {
    volume_frames = mri->nframes;
    MRIfree(&mri);
    errno = 0;
    ErrorReturn(
        NULL,
        (ERROR_BADPARM, "mri_read(): end frame (%d) is out of range (%d frames in volume)", end_frame, volume_frames));
  }

  if (!volume_flag) {
    if (nan_inf_check(mri) != NO_ERROR) {
      MRIfree(&mri);
      return (NULL);
    }
    mri2 = MRIcopy(mri, NULL);
    MRIfree(&mri);
    mri2->nframes = (end_frame - start_frame + 1);
    return (mri2);
  }

  if (start_frame == 0 && end_frame == mri->nframes - 1) {
    if (nan_inf_check(mri) != NO_ERROR) {
      MRIfree(&mri);
      return (NULL);
    }
    return (mri);
  }

  /* reading the whole volume and then copying only
     some frames is a bit inelegant,
     but it's easier for now
     to do this than to insert the appropriate
     code into the read function for each format... */
  /////////////////////////////////////////////////////////////////
  // This method does not work, because
  // 1. each frame has different acq parameters
  // 2. when making frames, it takes the info only from the first frame
  // Thus taking a frame out must be done at each frame extraction
  ////////////////////////////////////////////////////////////////////
  mri2 = MRIallocSequence(mri->width, mri->height, mri->depth, mri->type, (end_frame - start_frame + 1));
  MRIcopyHeader(mri, mri2);
  mri2->nframes = (end_frame - start_frame + 1);
  mri2->imnr0 = 1;
  mri2->imnr0 = mri2->nframes;

  if (0)   // disable it for now, don't know a case to execute this code path // (getenv("FS_MGZIO_USEVOXELBUF"))
  {
    printf("mri_read(): copy voxel by row ...\n");
    for (t = 0; t < mri2->nframes; t++)
      for (k = 0; k < mri2->depth; k++) {
        //for (j = 0; j < mri2->height; j++) {
          int bytes_to_copy = MRIsizeof(mri2->type) * mri2->width * mri2->height;
          void *pbufsrc = &MRIseq_vox(mri, 0, 0, k, t + start_frame);
          void *pbufdst = &MRIseq_vox(mri2, 0, 0, k, t);
	  memcpy(pbufdst, pbufsrc, bytes_to_copy);
	  //}
      }
  }
  else
  {
    if (mri2->type == MRI_UCHAR)
      for (t = 0; t < mri2->nframes; t++)
        for (i = 0; i < mri2->width; i++)
          for (j = 0; j < mri2->height; j++)
            for (k = 0; k < mri2->depth; k++) MRIseq_vox(mri2, i, j, k, t) = MRIseq_vox(mri, i, j, k, t + start_frame);

    if (mri2->type == MRI_SHORT)
      for (t = 0; t < mri2->nframes; t++)
        for (i = 0; i < mri2->width; i++)
          for (j = 0; j < mri2->height; j++)
            for (k = 0; k < mri2->depth; k++) MRISseq_vox(mri2, i, j, k, t) = MRISseq_vox(mri, i, j, k, t + start_frame);

    if (mri2->type == MRI_USHRT)
      for (t = 0; t < mri2->nframes; t++)
        for (i = 0; i < mri2->width; i++)
          for (j = 0; j < mri2->height; j++)
            for (k = 0; k < mri2->depth; k++) MRIUSseq_vox(mri2, i, j, k, t) = MRIUSseq_vox(mri, i, j, k, t + start_frame);

    if (mri2->type == MRI_INT)
      for (t = 0; t < mri2->nframes; t++)
        for (i = 0; i < mri2->width; i++)
          for (j = 0; j < mri2->height; j++)
            for (k = 0; k < mri2->depth; k++) MRIIseq_vox(mri2, i, j, k, t) = MRIIseq_vox(mri, i, j, k, t + start_frame);

    if (mri2->type == MRI_LONG)
      for (t = 0; t < mri2->nframes; t++)
        for (i = 0; i < mri2->width; i++)
          for (j = 0; j < mri2->height; j++)
            for (k = 0; k < mri2->depth; k++) MRILseq_vox(mri2, i, j, k, t) = MRILseq_vox(mri, i, j, k, t + start_frame);

    if (mri2->type == MRI_FLOAT)
      for (t = 0; t < mri2->nframes; t++)
        for (i = 0; i < mri2->width; i++)
          for (j = 0; j < mri2->height; j++)
            for (k = 0; k < mri2->depth; k++) MRIFseq_vox(mri2, i, j, k, t) = MRIFseq_vox(mri, i, j, k, t + start_frame);
  }

  if (nan_inf_check(mri) != NO_ERROR) {
    MRIfree(&mri);
    return (NULL);
  }

  MRIfree(&mri);

  return (mri2);

} /* end mri_read() */

static int nan_inf_check(MRI *mri)
{
  int i, j, k, t;

  if (mri->type != MRI_FLOAT) return (NO_ERROR);

  for (i = 0; i < mri->width; i++)
    for (j = 0; j < mri->height; j++)
      for (k = 0; k < mri->depth; k++)
        for (t = 0; t < mri->nframes; t++)
          if (!devFinite((MRIFseq_vox(mri, i, j, k, t)))) {
            if (devIsinf((MRIFseq_vox(mri, i, j, k, t))) != 0) {
              errno = 0;
              ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "nan_inf_check(): Inf at voxel %d, %d, %d, %d", i, j, k, t));
            }
            else if (devIsnan((MRIFseq_vox(mri, i, j, k, t)))) {
              errno = 0;
              ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "nan_inf_check(): NaN at voxel %d, %d, %d, %d", i, j, k, t));
            }
            else {
              errno = 0;
              ErrorReturn(ERROR_BADPARM,
                          (ERROR_BADPARM,
                           "nan_inf_check(): bizarre value (not Inf, "
                           "not NaN, but not finite) at %d, %d, %d, %d",
                           i,
                           j,
                           k,
                           t));
            }
          }

  return (NO_ERROR);

} /* end nan_inf_check() */

MRI *MRIreadType(const char *fname, int type)
{
  MRI *mri;

  chklc();

  mri = mri_read(fname, type, TRUE, -1, -1);

  return (mri);

} /* end MRIreadType() */

MRI *MRIread(const char *fname)
{
  char buf[STRLEN];
  MRI *mri = NULL;

  chklc();

  FileNameFromWildcard(fname, buf);
  fname = buf;
  int nstart = global_progress_range[0];
  int nend = global_progress_range[1];
  global_progress_range[1] = nstart + (nend - nstart) * 2 / 3;
  mri = mri_read(fname, MRI_VOLUME_TYPE_UNKNOWN, TRUE, -1, -1);

  /* some volume format needs to read many
     different files for slices (GE DICOM or COR).
     we make sure that mri_read() read the slices, not just one   */
  if (mri == NULL) return NULL;

  /* MRIremoveNaNs also takes long time so we split existing progress range
     between it and mri_read */
  global_progress_range[0] = global_progress_range[1];
  global_progress_range[1] = nend;
  MRIremoveNaNs(mri, mri);
  return (mri);

} /* end MRIread() */

// allow picking one frame out of many frame
// currently implemented only for Siemens dicom file
MRI *MRIreadEx(const char *fname, int nthframe)
{
  char buf[STRLEN];
  MRI *mri = NULL;

  chklc();

  FileNameFromWildcard(fname, buf);
  fname = buf;
  mri = mri_read(fname, MRI_VOLUME_TYPE_UNKNOWN, TRUE, nthframe, nthframe);

  /* some volume format needs to read many
     different files for slices (GE DICOM or COR).
     we make sure that mri_read() read the slices, not just one   */
  if (mri == NULL) return NULL;

  MRIremoveNaNs(mri, mri);
  return (mri);

} /* end MRIread() */

MRI *MRIreadInfo(const char *fname)
{
  MRI *mri = NULL;

  mri = mri_read(fname, MRI_VOLUME_TYPE_UNKNOWN, FALSE, -1, -1);

  return (mri);

} /* end MRIreadInfo() */

/*---------------------------------------------------------------
  MRIreadHeader() - reads the MRI header of the given file name.
  If type is MRI_VOLUME_TYPE_UNKNOWN, then the type will be
  inferred from the file name.
  ---------------------------------------------------------------*/
MRI *MRIreadHeader(const char *fname, int type)
{
  int usetype;
  MRI *mri = NULL;
  char modFname[STRLEN];
  struct stat stat_buf;

  usetype = type;

  if (!strchr(fname, DIR_SEPARATOR)) {
    char *cwd = getcwd(NULL, 0);  // posix 1 extension
    // (allocate as much space needed)
    if (cwd) {
      strcpy(modFname, cwd);
      strcat(modFname, "/");
      strcat(modFname, fname);
      free(cwd);
    }
    else  // why fail?
      strcpy(modFname, fname);
  }
  else
    strcpy(modFname, fname);

  if (usetype == MRI_VOLUME_TYPE_UNKNOWN) {
    usetype = mri_identify(modFname);
    if (usetype == MRI_VOLUME_TYPE_UNKNOWN) {
      // just check again
      if (stat(fname, &stat_buf) < 0)
        printf("ERROR: could not find volume %s.  Does it exist?\n", fname);
      else
        printf("ERROR: could not determine type of %s\n", fname);
      return (NULL);
    }
  }
  mri = mri_read(modFname, usetype, FALSE, -1, -1);

  return (mri);

} /* end MRIreadInfo() */

int MRIwriteType(MRI *mri, const char *fname, int type)
{
  struct stat stat_buf;
  int error = 0;
  char *fstem;
  char tmpstr[STRLEN];

  if (type == MRI_CORONAL_SLICE_DIRECTORY) {
    error = corWrite(mri, fname);
  }
  else {
    /* ----- all remaining types should write to
       a filename, not to within a directory,
       so check that it isn't an existing directory name we've been passed.
       failure to check will result in a segfault
       when attempting to write file data
       to a 'file' that is actually an existing directory ----- */
    if (stat(fname, &stat_buf) == 0) {  // if can stat, then fname exists
      if (S_ISDIR(stat_buf.st_mode)) {  // if is directory...
        errno = 0;
        ErrorReturn(ERROR_BADFILE,
                    (ERROR_BADFILE,
                     "MRIwriteType(): %s is an existing directory. (type=%d)\n"
                     "                A filename should be specified instead.",
                     fname,
                     type));
      }
    }
  }

  /* ----- continue file writing... ----- */

  if (type == MRI_MINC_FILE) {
    error = mincWrite(mri, fname);
  }
  else 
  if (type == IMAGE_FILE) {
    IMAGE *image;
    if (mri->depth != 1) ErrorExit(ERROR_BADPARM, "MRIwriteType(%s): image files cannnot have depth > 1\n", fname);
    image = MRItoImage(mri, NULL, 0);
    ImageWrite(image, fname);
    ImageFree(&image);
  }
  else if (type == BHDR) {
    fstem = bhdr_stem(fname);
    if (mri->type == MRI_SHORT) {
      sprintf(tmpstr, "%s_000.bshort", fstem);
      error = bvolumeWrite(mri, tmpstr, MRI_SHORT);
    }
    else {
      sprintf(tmpstr, "%s_000.bfloat", fstem);
      error = bvolumeWrite(mri, tmpstr, MRI_FLOAT);
    }
    free(fstem);
  }
  else if (type == BSHORT_FILE) {
    error = bvolumeWrite(mri, fname, MRI_SHORT);
  }
  else if (type == BFLOAT_FILE) {
    error = bvolumeWrite(mri, fname, MRI_FLOAT);
  }
  else if (type == MRI_ANALYZE_FILE) {
    error = analyzeWrite(mri, fname);
  }
  else if (type == MRI_ANALYZE4D_FILE) {
    error = analyzeWrite4D(mri, fname);
  }
  else if (type == BRIK_FILE) {
    error = afniWrite(mri, fname);
  }
  else if (type == MGH_MORPH) {
    GCA_MORPH *gcam;

    if (mri->nframes != 3)
      ErrorReturn(-1, (ERROR_UNSUPPORTED, "MRIwriteType: gcam - input mri must have 3 frames, not %d", mri->nframes));
    gcam = GCAMalloc(mri->width, mri->height, mri->depth);
    GCAMinitVolGeom(gcam, mri, mri);
    GCAMinit(gcam, mri, NULL, NULL, 0);
    GCAMreadWarpFromMRI(gcam, mri);
    GCAMwrite(gcam, fname);
    GCAMfree(&gcam);
  }
  else if (type == ITK_MORPH) {
    if (mri->nframes != 3)
      ErrorReturn(-1,
                  (ERROR_UNSUPPORTED, "MRIwriteType: itk warp - input mri must have 3 frames, not %d", mri->nframes));
    error = itkMorphWrite(mri, fname);
  }
  else if (type == MRI_MGH_FILE) {
    error = mghWrite(mri, fname, -1);
  }
  else if (type == GDF_FILE) {
    error = gdfWrite(mri, fname);
  }
  else if (type == NIFTI1_FILE) {
    error = nifti1Write(mri, fname);
  }
  else if (type == NII_FILE) {
    // printf("Before writing nii file \n");
    error = niiWrite(mri, fname);
    // printf("The error code is: %d\n", error);
  }
  else if (type == NRRD_FILE) {
    error = mriNrrdWrite(mri, fname);
  }
  else if (type == GIFTI_FILE) {
    error = mriWriteGifti(mri, fname);
  }
  else if (type == GENESIS_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of GENESIS file type not supported"));
  }
  else if (type == GE_LX_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of GE LX file type not supported"));
  }
  else if (type == SIEMENS_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of SIEMENS file type not supported"));
  }
  else if (type == SDT_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of SDT file type not supported"));
  }
  else if (type == OTL_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of OTL file type not supported"));
  }
  else if (type == RAW_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of RAW file type not supported"));
  }
  else if (type == SIGNA_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of SIGNA file type not supported"));
  }
  else if (type == DICOM_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of DICOM file type not supported"));
  }
  else if (type == SIEMENS_DICOM_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of SIEMENS DICOM file type not supported"));
  }
  else if (type == BRUKER_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of BRUKER file type not supported"));
  }
  else if (type == XIMG_FILE) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIwriteType(): writing of XIMG file type not supported"));
  }
  else if (type == MRI_CORONAL_SLICE_DIRECTORY) {
    // already processed above
  }
  else {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "MRIwriteType(): code inconsistency "
                 "(file type recognized but not caught)"));
  }
  if (error || mri->bvals == NULL) return (error);

  fstem = IDstemFromName(fname);
  if (fstem == NULL) return (error);

  printf("Saving bvals and bvecs\n");
  sprintf(tmpstr, "%s.bvals", fstem);
  DTIwriteBValues(mri->bvals, tmpstr);
  switch (mri->bvec_space) {
    case BVEC_SPACE_SCANNER:
      sprintf(tmpstr, "%s.scanner_space.bvecs", fstem);
      break;
    case BVEC_SPACE_VOXEL:
      sprintf(tmpstr, "%s.voxel_space.bvecs", fstem);
      break;
    default:
      sprintf(tmpstr, "%s.unknown_space.bvecs", fstem);
      break;
  }
  DTIwriteBVectors(mri->bvecs, tmpstr);
  free(fstem);

  return (error);

} /* end MRIwriteType() */

int MRIwriteFrame(MRI *mri, const char *fname, int frame)
{
  MRI *mri_tmp;

  if (frame >= mri->nframes)
    ErrorExit(ERROR_BADPARM, "MRIwriteFrame(%d) - frame out of bounds (%d)", frame, mri->nframes);
  mri_tmp = MRIcopyFrame(mri, NULL, frame, 0);
  MRIwrite(mri_tmp, fname);
  MRIfree(&mri_tmp);
  return (NO_ERROR);
}

int MRIwrite(MRI *mri, const char *fname)
{
  int int_type = -1;
  int error;

  chklc();
  if ((int_type = mri_identify(fname)) < 0) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "unknown file type for file (%s)", fname));
  }

  error = MRIwriteType(mri, fname, int_type);
  return (error);

} /* end MRIwrite() */

/* ----- required header fields ----- */

#define COR_ALL_REQUIRED 0x00001fff

#define IMNR0_FLAG 0x00000001
#define IMNR1_FLAG 0x00000002
#define PTYPE_FLAG 0x00000004
#define X_FLAG 0x00000008
#define Y_FLAG 0x00000010
#define THICK_FLAG 0x00000020
#define PSIZ_FLAG 0x00000040
#define STRTX_FLAG 0x00000080
#define ENDX_FLAG 0x00000100
#define STRTY_FLAG 0x00000200
#define ENDY_FLAG 0x00000400
#define STRTZ_FLAG 0x00000800
#define ENDZ_FLAG 0x00001000

/* trivially time course clean */
static MRI *corRead(const char *fname, int read_volume)
{
  MRI *mri;
  struct stat stat_buf;
  char fname_use[STRLEN];
  char *fbase;
  FILE *fp;
  int i, j;
  char line[STRLEN];
  int imnr0, imnr1, x, y, ptype;
  double fov, thick, psiz, locatn; /* using floats to read creates problems
                                              when checking values
                                              (e.g. thick = 0.00100000005) */
  float strtx, endx, strty, endy, strtz, endz;
  float tr, te, ti, flip_angle;
  int ras_good_flag;
  float x_r, x_a, x_s;
  float y_r, y_a, y_s;
  float z_r, z_a, z_s;
  float c_r, c_a, c_s;
  char xform[STRLEN];
  long gotten;
  char *cur_char;
  char *last_slash;

  /* ----- check that it is a directory we've been passed ----- */
  strcpy(fname_use, fname);
  if (stat(fname_use, &stat_buf) < 0) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "corRead(): can't stat %s", fname_use));
  }

  if (!S_ISDIR(stat_buf.st_mode)) {
    /* remove the last bit and try again */
    cur_char = fname_use;
    last_slash = cur_char;
    while (*(cur_char + 1) != '\0') {
      if (*cur_char == '/') last_slash = cur_char;
      cur_char++;
    }
    *last_slash = '\0';
    if (stat(fname_use, &stat_buf) < 0) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADFILE, "corRead(): can't stat %s", fname_use));
    }
    if (!S_ISDIR(stat_buf.st_mode)) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADFILE, "corRead(): %s isn't a directory", fname_use));
    }
  }

  /* ----- copy the directory name and remove any trailing '/' ----- */
  fbase = &fname_use[strlen(fname_use)];
  if (*(fbase - 1) != '/') {
    *fbase = '/';
    fbase++;
  }

  /* ----- read the header file ----- */
  sprintf(fbase, "COR-.info");
  if ((fp = fopen(fname_use, "r")) == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "corRead(): can't open file %s", fname_use));
  }

  /* ----- defaults (a good idea for non-required values...) ----- */
  xform[0] = '\0';
  ras_good_flag = 0;
  x_r = x_a = x_s = 0.0;
  y_r = y_a = y_s = 0.0;
  z_r = z_a = z_s = 0.0;
  c_r = c_a = c_s = 0.0;
  flip_angle = tr = te = ti = 0.0;
  fov = 0.0;
  locatn = 0.0;

  gotten = 0x00;

  while (fgets(line, STRLEN, fp) != NULL) {
    if (strncmp(line, "imnr0 ", 6) == 0) {
      sscanf(line, "%*s %d", &imnr0);
      gotten = gotten | IMNR0_FLAG;
    }
    else if (strncmp(line, "imnr1 ", 6) == 0) {
      sscanf(line, "%*s %d", &imnr1);
      gotten = gotten | IMNR1_FLAG;
    }
    else if (strncmp(line, "ptype ", 6) == 0) {
      sscanf(line, "%*s %d", &ptype);
      gotten = gotten | PTYPE_FLAG;
    }
    else if (strncmp(line, "x ", 2) == 0) {
      sscanf(line, "%*s %d", &x);
      gotten = gotten | X_FLAG;
    }
    else if (strncmp(line, "y ", 2) == 0) {
      sscanf(line, "%*s %d", &y);
      gotten = gotten | Y_FLAG;
    }
    else if (strncmp(line, "fov ", 4) == 0) {
      sscanf(line, "%*s %lf", &fov);
    }
    else if (strncmp(line, "thick ", 6) == 0) {
      sscanf(line, "%*s %lf", &thick);
      gotten = gotten | THICK_FLAG;
    }
    else if (strncmp(line, "flip ", 5) == 0) {
      sscanf(line + 11, "%f", &flip_angle);
      flip_angle = RADIANS(flip_angle);
    }
    else if (strncmp(line, "psiz ", 5) == 0) {
      sscanf(line, "%*s %lf", &psiz);
      gotten = gotten | PSIZ_FLAG;
    }
    else if (strncmp(line, "locatn ", 7) == 0) {
      sscanf(line, "%*s %lf", &locatn);
    }
    else if (strncmp(line, "strtx ", 6) == 0) {
      sscanf(line, "%*s %f", &strtx);
      gotten = gotten | STRTX_FLAG;
    }
    else if (strncmp(line, "endx ", 5) == 0) {
      sscanf(line, "%*s %f", &endx);
      gotten = gotten | ENDX_FLAG;
    }
    else if (strncmp(line, "strty ", 6) == 0) {
      sscanf(line, "%*s %f", &strty);
      gotten = gotten | STRTY_FLAG;
    }
    else if (strncmp(line, "endy ", 5) == 0) {
      sscanf(line, "%*s %f", &endy);
      gotten = gotten | ENDY_FLAG;
    }
    else if (strncmp(line, "strtz ", 6) == 0) {
      sscanf(line, "%*s %f", &strtz);
      gotten = gotten | STRTZ_FLAG;
    }
    else if (strncmp(line, "endz ", 5) == 0) {
      sscanf(line, "%*s %f", &endz);
      gotten = gotten | ENDZ_FLAG;
    }
    else if (strncmp(line, "tr ", 3) == 0) {
      sscanf(line, "%*s %f", &tr);
    }
    else if (strncmp(line, "te ", 3) == 0) {
      sscanf(line, "%*s %f", &te);
    }
    else if (strncmp(line, "ti ", 3) == 0) {
      sscanf(line, "%*s %f", &ti);
    }
    else if (strncmp(line, "ras_good_flag ", 14) == 0) {
      sscanf(line, "%*s %d", &ras_good_flag);
    }
    else if (strncmp(line, "x_ras ", 6) == 0) {
      sscanf(line, "%*s %f %f %f", &x_r, &x_a, &x_s);
    }
    else if (strncmp(line, "y_ras ", 6) == 0) {
      sscanf(line, "%*s %f %f %f", &y_r, &y_a, &y_s);
    }
    else if (strncmp(line, "z_ras ", 6) == 0) {
      sscanf(line, "%*s %f %f %f", &z_r, &z_a, &z_s);
    }
    else if (strncmp(line, "c_ras ", 6) == 0) {
      sscanf(line, "%*s %f %f %f", &c_r, &c_a, &c_s);
    }
    else if (strncmp(line, "xform", 5) == 0 || strncmp(line, "transform", 9) == 0) {
      sscanf(line, "%*s %s", xform);
    }
  }

  fclose(fp);

  /* ----- check for required fields ----- */
  if ((gotten & COR_ALL_REQUIRED) != COR_ALL_REQUIRED) {
    errno = 0;
    ErrorPrintf(ERROR_BADFILE, "missing fields in file %s:", fname_use);
    if (!(gotten & IMNR0_FLAG)) ErrorPrintf(ERROR_BADFILE, "  imnr0 field missing");
    if (!(gotten & IMNR1_FLAG)) ErrorPrintf(ERROR_BADFILE, "  imnr1 field missing");
    if (!(gotten & PTYPE_FLAG)) ErrorPrintf(ERROR_BADFILE, "  ptype field missing");
    if (!(gotten & X_FLAG)) ErrorPrintf(ERROR_BADFILE, "  x field missing");
    if (!(gotten & Y_FLAG)) ErrorPrintf(ERROR_BADFILE, "  y field missing");
    if (!(gotten & THICK_FLAG)) ErrorPrintf(ERROR_BADFILE, "  thick field missing");
    if (!(gotten & PSIZ_FLAG)) ErrorPrintf(ERROR_BADFILE, "  psiz field missing");
    if (!(gotten & STRTX_FLAG)) ErrorPrintf(ERROR_BADFILE, "  strtx field missing");
    if (!(gotten & ENDX_FLAG)) ErrorPrintf(ERROR_BADFILE, "  endx field missing");
    if (!(gotten & STRTY_FLAG)) ErrorPrintf(ERROR_BADFILE, "  strty field missing");
    if (!(gotten & ENDY_FLAG)) ErrorPrintf(ERROR_BADFILE, "  endy field missing");
    if (!(gotten & STRTZ_FLAG)) ErrorPrintf(ERROR_BADFILE, "  strtz field missing");
    if (!(gotten & ENDZ_FLAG)) ErrorPrintf(ERROR_BADFILE, "  endz field missing");
    return (NULL);
  }

  /* ----- check for required but forced (constant) values ----- */

  if (imnr0 != 1) {
    printf(
        "warning: non-standard value for imnr0 "
        "(%d, usually 1) in file %s\n",
        imnr0,
        fname_use);
  }

  if (imnr1 != 256) {
    printf(
        "warning: non-standard value for imnr1 "
        "(%d, usually 256) in file %s\n",
        imnr1,
        fname_use);
  }

  if (ptype != 2) {
    printf(
        "warning: non-standard value for ptype "
        "(%d, usually 2) in file %s\n",
        ptype,
        fname_use);
  }

  if (x != 256) {
    printf(
        "warning: non-standard value for x "
        "(%d, usually 256) in file %s\n",
        x,
        fname_use);
  }

  if (y != 256) {
    printf(
        "warning: non-standard value for y "
        "(%d, usually 256) in file %s\n",
        y,
        fname_use);
  }

  if (thick != 0.001) {
    printf(
        "warning: non-standard value for thick "
        "(%g, usually 0.001) in file %s\n",
        thick,
        fname_use);
  }

  if (psiz != 0.001) {
    printf(
        "warning: non-standard value for psiz "
        "(%g, usually 0.001) in file %s\n",
        psiz,
        fname_use);
  }

  /* ----- copy header information to an mri structure ----- */

  if (read_volume)
    mri = MRIalloc(x, y, imnr1 - imnr0 + 1, MRI_UCHAR);
  else
    mri = MRIallocHeader(x, y, imnr1 - imnr0 + 1, MRI_UCHAR, 1);

  /* hack */
  /*
    printf("%g, %g, %g\n", x_r, x_a, x_s);
    printf("%g, %g, %g\n", y_r, y_a, y_s);
    printf("%g, %g, %g\n", z_r, z_a, z_s);
  */
  if (x_r == 0.0 && x_a == 0.0 && x_s == 0.0 && y_r == 0.0 && y_a == 0.0 && y_s == 0.0 && z_r == 0.0 && z_a == 0.0 &&
      z_s == 0.0 && Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stderr,
            "----------------------------------------------------\n"
            "Could not find the direction cosine information.\n"
            "Will use the CORONAL orientation.\n"
            "If not suitable, please provide the information "
            "in COR-.info file\n"
            "----------------------------------------------------\n");
    x_r = -1.0;
    y_s = -1.0;
    z_a = 1.0;
    ras_good_flag = 0;
  }

  mri->imnr0 = imnr0;
  mri->imnr1 = imnr1;
  mri->fov = (float)(fov * 1000);
  mri->thick = (float)(thick * 1000);
  mri->ps = (float)(psiz * 1000);
  mri->xsize = mri->ps;
  mri->ysize = mri->ps;
  mri->zsize = (float)(mri->thick);
  mri->xstart = strtx * 1000;
  mri->xend = endx * 1000;
  mri->ystart = strty * 1000;
  mri->yend = endy * 1000;
  mri->zstart = strtz * 1000;
  mri->zend = endz * 1000;
  strcpy(mri->fname, fname);
  mri->tr = tr;
  mri->te = te;
  mri->ti = ti;
  mri->flip_angle = flip_angle;
  mri->ras_good_flag = ras_good_flag;
  mri->x_r = x_r;
  mri->x_a = x_a;
  mri->x_s = x_s;
  mri->y_r = y_r;
  mri->y_a = y_a;
  mri->y_s = y_s;
  mri->z_r = z_r;
  mri->z_a = z_a;
  mri->z_s = z_s;
  mri->c_r = c_r;
  mri->c_a = c_a;
  mri->c_s = c_s;

  if (strlen(xform) > 0) {
    char xform_use[STRLEN];
    if (xform[0] == '/') {
      strcpy(mri->transform_fname, xform);
    } else {
      int req = snprintf(mri->transform_fname, STRLEN, "%s/%s", fname, xform);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }

    strcpy(xform_use, mri->transform_fname);

    if (!FileExists(xform_use)) {
      sprintf(xform_use, "%s/../transforms/talairach.xfm", fname);

      if (!FileExists(xform_use)) printf("can't find talairach file '%s'\n", xform_use);
    }

    if (FileExists(xform_use)) {
      if (input_transform_file(xform_use, &(mri->transform)) == NO_ERROR) {
        mri->linear_transform = get_linear_transform_ptr(&mri->transform);
        mri->inverse_linear_transform = get_inverse_linear_transform_ptr(&mri->transform);
        mri->free_transform = 1;
        strcpy(mri->transform_fname, xform_use);
        if (DIAG_VERBOSE_ON) fprintf(stderr, "INFO: loaded talairach xform : %s\n", mri->transform_fname);
      }
      else {
        errno = 0;
        ErrorPrintf(ERROR_BAD_FILE, "error loading transform from %s", xform_use);
        mri->linear_transform = NULL;
        mri->inverse_linear_transform = NULL;
        mri->free_transform = 1;
        (mri->transform_fname)[0] = '\0';
      }
    }
  }

  if (!read_volume) return (mri);

  /* ----- read the data files ----- */
  for (i = mri->imnr0; i <= imnr1; i++) {
    sprintf(fbase, "COR-%03d", i);
    if ((fp = fopen(fname_use, "r")) == NULL) {
      MRIfree(&mri);
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADFILE, "corRead(): can't open file %s", fname_use));
    }
    for (j = 0; j < mri->height; j++) {
      if ((int)fread(mri->slices[i - mri->imnr0][j], 1, mri->width, fp) < mri->width) {
        MRIfree(&mri);
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADFILE, "corRead(): error reading from file %s", fname_use));
      }
    }
    fclose(fp);
    exec_progress_callback(i - mri->imnr0, imnr1 - imnr0 + 1, 0, 1);
  }

  return (mri);

} /* end corRead() */

static int corWrite(MRI *mri, const char *fname)
{
  struct stat stat_buf;
  char fname_use[STRLEN];
  char *fbase;
  FILE *fp;
  int i, j;
  int rv;

  /* ----- check the mri structure for COR file compliance ----- */

  if (mri->slices == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "corWrite(): mri structure to be "
                 "written contains no voxel data"));
  }

  if (mri->imnr0 != 1) {
    printf("non-standard value for imnr0 (%d, usually 1) in volume structure\n", mri->imnr0);
  }

  if (mri->imnr1 != 256) {
    printf(
        "non-standard value for imnr1 (%d, "
        "usually 256) in volume structure\n",
        mri->imnr1);
  }

  if (mri->type != MRI_UCHAR) {
    printf(
        "non-standard value for type "
        "(%d, usually %d) in volume structure\n",
        mri->type,
        MRI_UCHAR);
  }

  if (mri->width != 256) {
    printf(
        "non-standard value for width "
        "(%d, usually 256) in volume structure\n",
        mri->width);
  }

  if (mri->height != 256) {
    printf(
        "non-standard value for height "
        "(%d, usually 256) in volume structure\n",
        mri->height);
  }

  if (mri->thick != 1) {
    printf(
        "non-standard value for thick "
        "(%g, usually 1) in volume structure\n",
        mri->thick);
  }

  if (mri->ps != 1) {
    printf(
        "non-standard value for ps "
        "(%g, usually 1) in volume structure\n",
        mri->ps);
  }

  /* ----- copy the directory name and remove any trailing '/' ----- */
  strcpy(fname_use, fname);
  fbase = &fname_use[strlen(fname_use)];
  if (*(fbase - 1) != '/') {
    *fbase = '/';
    fbase++;
    *fbase = '\0';
  }

  /* Create the directory */
  rv = mkdir(fname_use, 0777);
  if (rv != 0 && errno != EEXIST) {
    printf("ERROR: creating directory %s\n", fname_use);
    perror(NULL);
    return (1);
  }

  /* ----- check that it is a directory we've been passed ----- */
  if (stat(fname, &stat_buf) < 0) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "corWrite(): can't stat %s", fname));
  }

  if (!S_ISDIR(stat_buf.st_mode)) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "corWrite(): %s isn't a directory", fname));
  }

  sprintf(fbase, "COR-.info");
  if ((fp = fopen(fname_use, "w")) == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "corWrite(): can't open file %s for writing", fname_use));
  }

  fprintf(fp, "imnr0 %d\n", mri->imnr0);
  fprintf(fp, "imnr1 %d\n", mri->imnr1);
  fprintf(fp, "ptype %d\n", 2);
  fprintf(fp, "x %d\n", mri->width);
  fprintf(fp, "y %d\n", mri->height);
  fprintf(fp, "fov %g\n", mri->fov / 1000.0);
  fprintf(fp, "thick %g\n", mri->zsize / 1000.0); /* was xsize 11/1/01 */
  fprintf(fp, "psiz %g\n", mri->xsize / 1000.0);  /* was zsize 11/1/01 */
  fprintf(fp, "locatn %g\n", 0.0);
  fprintf(fp, "strtx %g\n", mri->xstart / 1000.0);
  fprintf(fp, "endx %g\n", mri->xend / 1000.0);
  fprintf(fp, "strty %g\n", mri->ystart / 1000.0);
  fprintf(fp, "endy %g\n", mri->yend / 1000.0);
  fprintf(fp, "strtz %g\n", mri->zstart / 1000.0);
  fprintf(fp, "endz %g\n", mri->zend / 1000.0);
  fprintf(fp, "tr %f\n", mri->tr);
  fprintf(fp, "te %f\n", mri->te);
  fprintf(fp, "ti %f\n", mri->ti);
  fprintf(fp, "flip angle %f\n", DEGREES(mri->flip_angle));
  fprintf(fp, "ras_good_flag %d\n", mri->ras_good_flag);
  fprintf(fp, "x_ras %f %f %f\n", mri->x_r, mri->x_a, mri->x_s);
  fprintf(fp, "y_ras %f %f %f\n", mri->y_r, mri->y_a, mri->y_s);
  fprintf(fp, "z_ras %f %f %f\n", mri->z_r, mri->z_a, mri->z_s);
  fprintf(fp, "c_ras %f %f %f\n", mri->c_r, mri->c_a, mri->c_s);
  if (strlen(mri->transform_fname) > 0) fprintf(fp, "xform %s\n", mri->transform_fname);

  fclose(fp);

  for (i = mri->imnr0; i <= mri->imnr1; i++) {
    sprintf(fbase, "COR-%03d", i);
    if ((fp = fopen(fname_use, "w")) == NULL) {
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "corWrite(): can't open file %s for writing", fname_use));
    }
    for (j = 0; j < mri->height; j++) {
      if ((int)fwrite(mri->slices[i - mri->imnr0][j], 1, mri->width, fp) < mri->width) {
        fclose(fp);
        errno = 0;
        ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "corWrite(): error writing to file %s ", fname_use));
      }
    }
    fclose(fp);
    exec_progress_callback(i - mri->imnr0, mri->imnr1 - mri->imnr0 + 1, 0, 1);
  }

  return (NO_ERROR);

} /* end corWrite() */

static MRI *siemensRead(const char *fname, int read_volume_flag)
{
  int file_n, n_low, n_high;
  char fname_use[STRLEN];
  MRI *mri;
  FILE *fp;
  char *c, *c2;
  short rows, cols;
  int n_slices;
  double d, d2;
  double im_c_r, im_c_a, im_c_s;
  int i, j;
  int n_files, base_raw_matrix_size, number_of_averages;
  int mosaic_size;
  int mosaic;
  int mos_r, mos_c;
  char pulse_sequence_name[STRLEN], ps2[STRLEN];
  int n_t;
  int n_dangling_images, n_full_mosaics, mosaics_per_volume;
  int br, bc;
  MRI *mri_raw;
  int t, s;
  int slice_in_mosaic;
  int file;
  char ima[4];
  IMAFILEINFO *ifi;

  /* ----- stop compiler complaints ----- */
  mri = NULL;
  mosaic_size = 0;

  printf("Starting siemensRead()\n");

  /* Check whether it is really a dicom file with ima extension*/
  if (IsDICOM(fname)) {
    printf("INFO: this looks like a dicom, not an old Siemens .ima file, ");
    printf("so I'm going to unpack it as a dicom.\n");
    printf("If this fails, run mri_convert with -it dicom\n");
    // mri = DICOMRead2(fname, read_volume_flag); // generic dicom
    mri = sdcmLoadVolume(fname, read_volume_flag, -1);  // siemens dicom
    return (mri);
  }

  ifi = imaLoadFileInfo(fname);
  if (ifi == NULL) {
    printf("ERROR: siemensRead(): %s\n", fname);
    return (NULL);
  }

  strcpy(fname_use, fname);

  /* ----- point to ".ima" ----- */
  c = fname_use + (strlen(fname_use) - 4);

  if (c < fname_use) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "siemensRead(): bad file name %s (must end in '.ima' or '.IMA')", fname_use));
  }
  if (strcmp(".ima", c) == 0)
    sprintf(ima, "ima");
  else if (strcmp(".IMA", c) == 0)
    sprintf(ima, "IMA");
  else {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "siemensRead(): bad file name %s (must end in '.ima' or '.IMA')", fname_use));
  }

  c2 = c;
  for (c--; isdigit((int)(*c)); c--)
    ;
  c++;

  /* ----- c now points to the first digit in the last number set
     (e.g. to the "5" in 123-4-567.ima) ----- */

  /* ----- count down and up from this file --
     max and min image number within the series ----- */
  *c2 = '\0';
  file_n = atol(c);
  *c2 = '.';

  if (!FileExists(fname_use)) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "siemensRead(): file %s doesn't exist", fname_use));
  }

  /* --- get the low image number --- */
  n_low = file_n - 1;
  sprintf(c, "%d.%s", n_low, ima);
  while (FileExists(fname_use)) {
    n_low--;
    sprintf(c, "%d.%s", n_low, ima);
  }
  n_low++;

  /* --- get the high image number --- */
  n_high = file_n + 1;
  sprintf(c, "%d.%s", n_high, ima);
  while (FileExists(fname_use)) {
    n_high++;
    sprintf(c, "%d.%s", n_high, ima);
  }
  n_high--;

  n_files = n_high - n_low + 1;

  sprintf(c, "%d.%s", n_low, ima);
  if ((fp = fopen(fname_use, "r")) == NULL) {
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "siemensRead(): can't open file %s (low = %d, this = %d, high = %d)",
                 fname_use,
                 n_low,
                 file_n,
                 n_high));
  }

  fseek(fp, 4994, SEEK_SET);
  if(fread(&rows, 2, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  rows = orderShortBytes(rows);
  printf("rows = %d\n", rows);
  fseek(fp, 4996, SEEK_SET);
  if (fread(&cols, 2, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  cols = orderShortBytes(cols);
  printf("cols = %d\n", cols);
  fseek(fp, 4004, SEEK_SET);
  if (fread(&n_slices, 4, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  
  n_slices = orderIntBytes(n_slices);
  if (n_slices == 0) {
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "\n\nPlease try with the option '-it siemens_dicom'.\n "
                 "The handling failed assuming the old siemens format.\n"))
  }
  fseek(fp, 2864, SEEK_SET);
  if (fread(&base_raw_matrix_size, 4, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  base_raw_matrix_size = orderIntBytes(base_raw_matrix_size);
  fseek(fp, 1584, SEEK_SET);
  if (fread(&number_of_averages, 4, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  number_of_averages = orderIntBytes(number_of_averages);
  memset(pulse_sequence_name, 0x00, STRLEN);
  fseek(fp, 3009, SEEK_SET);
  if (fread(&pulse_sequence_name, 1, 65, fp) != 65){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }

  /* --- scout --- */
  strcpy(ps2, pulse_sequence_name);
  StrLower(ps2);
  if (strstr(ps2, "scout") != NULL) {
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "siemensRead(): series appears to be a scout "
                 "(sequence file name is %s)",
                 pulse_sequence_name));
  }

  /* --- structural --- */
  if (n_slices == 1 && !ifi->IsMosaic) {
    n_slices = n_files;
    n_t = 1;
    if (base_raw_matrix_size != rows || base_raw_matrix_size != cols) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADPARM, "siemensRead(): bad file/base matrix sizes"));
    }
    mos_r = mos_c = 1;
    mosaic_size = 1;
  }
  else {
    if (rows % base_raw_matrix_size != 0) {
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_BADPARM,
                   "siemensRead(): file rows (%hd) not"
                   " divisible by image rows (%d)",
                   rows,
                   base_raw_matrix_size));
    }
    if (cols % base_raw_matrix_size != 0) {
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_BADPARM,
                   "siemensRead(): file cols (%hd)"
                   " not divisible by image cols (%d)",
                   cols,
                   base_raw_matrix_size));
    }

    mos_r = rows / base_raw_matrix_size;
    mos_c = cols / base_raw_matrix_size;
    mosaic_size = mos_r * mos_c;
    if (mosaic_size == 0) {
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_BADPARM,
                   "siemensRead(): mosaic_size==0, "
                   "Try '-it siemens_dicom' flag with mri_convert"));
    }

    n_dangling_images = n_slices % mosaic_size;
    n_full_mosaics = (n_slices - n_dangling_images) / mosaic_size;

    mosaics_per_volume = n_full_mosaics + (n_dangling_images == 0 ? 0 : 1);

    if (n_files % mosaics_per_volume != 0) {
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_BADPARM,
                   "siemensRead(): files in volume "
                   "(%d) not divisible by mosaics per volume (%d)",
                   n_files,
                   mosaics_per_volume));
    }

    n_t = n_files / mosaics_per_volume;
  }

  if (read_volume_flag)
    mri = MRIallocSequence(base_raw_matrix_size, base_raw_matrix_size, n_slices, MRI_SHORT, n_t);
  else {
    mri = MRIallocHeader(base_raw_matrix_size, base_raw_matrix_size, n_slices, MRI_SHORT, n_t);
    mri->nframes = n_t;
  }

  /* --- pixel sizes --- */
  /* --- mos_r and mos_c factors are strange, but they're there... --- */
  fseek(fp, 5000, SEEK_SET);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->xsize = mos_r * orderDoubleBytes(d);
  fseek(fp, 5008, SEEK_SET);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->ysize = mos_c * orderDoubleBytes(d);

  /* --- slice distance factor --- */
  fseek(fp, 4136, SEEK_SET);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  d = orderDoubleBytes(d);
  if (d == -19222) /* undefined (Siemens code) -- I assume this to mean 0 */
    d = 0.0;
  /* --- slice thickness --- */
  fseek(fp, 1544, SEEK_SET);
  if (fread(&d2, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  d2 = orderDoubleBytes(d2);
  /* --- distance between slices --- */
  mri->zsize = (1.0 + d) * d2;

  /* --- field of view (larger of height, width fov) --- */
  fseek(fp, 3744, SEEK_SET);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  d = orderDoubleBytes(d);
  fseek(fp, 3752, SEEK_SET);
  if (fread(&d2, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  d2 = orderDoubleBytes(d);
  mri->fov = (d > d2 ? d : d2);

  mri->thick = mri->zsize;
  mri->ps = mri->xsize;

  strcpy(mri->fname, fname);

  mri->location = 0.0;

  fseek(fp, 2112, SEEK_SET);
  mri->flip_angle = RADIANS(freadDouble(fp)); /* was in degrees */

  fseek(fp, 1560, SEEK_SET);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->tr = orderDoubleBytes(d);
  fseek(fp, 1568, SEEK_SET);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->te = orderDoubleBytes(d);
  fseek(fp, 1576, SEEK_SET);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->ti = orderDoubleBytes(d);

  fseek(fp, 3792, SEEK_SET);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->z_r = -orderDoubleBytes(d);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->z_a = orderDoubleBytes(d);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->z_s = -orderDoubleBytes(d);

  fseek(fp, 3832, SEEK_SET);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->x_r = -orderDoubleBytes(d);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->x_a = orderDoubleBytes(d);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->x_s = -orderDoubleBytes(d);

  fseek(fp, 3856, SEEK_SET);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->y_r = -orderDoubleBytes(d);
 if (fread(&d, 8, 1, fp) != 1){
   ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
 }
  mri->y_a = orderDoubleBytes(d);
  if (fread(&d, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  mri->y_s = -orderDoubleBytes(d);

  fseek(fp, 3768, SEEK_SET);
  if (fread(&im_c_r, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  im_c_r = -orderDoubleBytes(im_c_r);
  if (fread(&im_c_a, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  im_c_a = orderDoubleBytes(im_c_a);
  if (fread(&im_c_s, 8, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  im_c_s = -orderDoubleBytes(im_c_s);

  mri->c_r = im_c_r - (mosaic_size - 1) * mri->z_r * mri->zsize + ((mri->depth) / 2.0) * mri->z_r * mri->zsize;
  mri->c_a = im_c_a - (mosaic_size - 1) * mri->z_a * mri->zsize + ((mri->depth) / 2.0) * mri->z_a * mri->zsize;
  mri->c_s = im_c_s - (mosaic_size - 1) * mri->z_s * mri->zsize + ((mri->depth) / 2.0) * mri->z_s * mri->zsize;

  mri->ras_good_flag = 1;

  fseek(fp, 3760, SEEK_SET);
  if (fread(&i, 4, 1, fp) != 1){
    ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
  }
  i = orderIntBytes(i);

  if (i < 1 || i > 6) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "siemensRead(): bad slice direction (%d) in file %s", i, fname_use));
  }

  mri->xstart = -mri->width * mri->xsize / 2.0;
  mri->xend = -mri->xstart;
  mri->ystart = -mri->height * mri->ysize / 2.0;
  mri->yend = -mri->ystart;
  mri->zstart = -mri->depth * mri->zsize / 2.0;
  mri->zend = -mri->zstart;

  /*
    printf("%d, %d; %d, %hd, %hd; %d\n", n_files, number_of_averages,
    base_raw_matrix_size, rows, cols,
    slices);
  */
  /*
    rows, cols, brms, mosaic i, j, mosaic size, n slices, n files, n_t
  */
  fclose(fp);

  mri->imnr0 = 1;
  mri->imnr1 = mri->depth;
  /*
    printf("%d, %d, %d, %d\n",
    mri->width, mri->height, mri->depth, mri->nframes);
  */
  if (read_volume_flag) {
    mri_raw = MRIalloc(rows, cols, n_files, MRI_SHORT);

    for (file_n = n_low; file_n <= n_high; file_n++) {
      sprintf(c, "%d.%s", file_n, ima);
      if ((fp = fopen(fname_use, "r")) == NULL) {
        MRIfree(&mri);
        errno = 0;
        ErrorReturn(NULL,
                    (ERROR_BADFILE,
                     "siemensRead(): can't open file %s "
                     "(low = %d, this = %d, high = %d)",
                     fname_use,
                     n_low,
                     file_n,
                     n_high));
      }

      fseek(fp, 6144, SEEK_SET);

      for (i = 0; i < rows; i++) {
        if (fread(&MRISvox(mri_raw, 0, i, file_n - n_low), sizeof(short), cols, fp) != (unsigned)cols){
          ErrorPrintf(ERROR_BADFILE, "siemensRead(): could not read file");
        }
#if (BYTE_ORDER == LITTLE_ENDIAN)
#if defined(SunOS)
        swab((const char *)&MRISvox(mri_raw, 0, i, file_n - n_low),
             (char *)&MRISvox(mri_raw, 0, i, file_n - n_low),
             sizeof(short) * cols);
#else
	{
	  std::vector<short> temp(cols);
	  // Note:
	  // void swab(const void *from, void *to, ssize_t n);
	  // void *memcpy(void *dest, const void *src, size_t n);
	  // Because consistency is the hobgoblin of small minds...
	  swab(&MRISvox(mri_raw, 0, i, file_n - n_low), temp.data(), sizeof(short)*cols );
	  memcpy(&MRISvox(mri_raw, 0, i, file_n - n_low), temp.data(), sizeof(short)*cols );
	}
#endif
#endif
      }

      fclose(fp);
    }

    for (t = 0; t < mri->nframes; t++) {
      for (s = 0; s < mri->depth; s++) {
        slice_in_mosaic = s % mosaic_size;
        mosaic = (s - slice_in_mosaic) / mosaic_size;
        file = mri->nframes * mosaic + t;
        /*
          printf("s, t = %d, %d; f, sm = %d, %d\n",
          s, t, file, slice_in_mosaic);
        */
        bc = slice_in_mosaic % mos_r;
        br = (slice_in_mosaic - bc) / mos_r;

        for (i = 0; i < mri->width; i++) {
          for (j = 0; j < mri->height; j++) {
            MRISseq_vox(mri, i, j, s, t) = MRISvox(mri_raw, mri->width * bc + i, mri->height * br + j, file);
          }
        }
        exec_progress_callback(s, mri->depth, t, mri->nframes);
      }
    }

    MRIfree(&mri_raw);
  }

  return (mri);

} /* end siemensRead() */
/*-----------------------------------------------------------*/
static MRI *mincRead(const char *fname, int read_volume)
{
  // double wx, wy, wz;
  MRI *mri;
  Volume vol;
  VIO_Status status;
  const char *dim_names[4];
  int dim_sizes[4];
  int ndims;
  int dtype;
  volume_input_struct input_info;
  double separations[4];
  double voxel[4];
  double worldr, worlda, worlds;
  double val;
  int i, j, k, t;
  float xfov, yfov, zfov;
  double f;
  BOOLEAN sflag = TRUE;
  Transform *pVox2WorldLin;
  General_transform *pVox2WorldGen;

  /* ----- read the volume ----- */

  /*   we no longer much around axes and thus */
  /*   it is safe to read in the standard way    */
  dim_names[0] = MIxspace;
  dim_names[1] = MIyspace;
  dim_names[2] = MIzspace;
  dim_names[3] = MItime;

  if (!FileExists(fname)) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "mincRead(): can't find file %s", fname));
  }

  char *tmp = strcpyalloc(fname);
  status = start_volume_input(tmp, 0, const_cast<char**>(dim_names), NC_UNSPECIFIED, 0, 0, 0, TRUE, &vol, NULL, &input_info);
  free(tmp);

  if (Gdiag & DIAG_VERBOSE_ON & DIAG_SHOW) {
    printf("status = %d\n", status);
    printf("n_dimensions = %d\n", get_volume_n_dimensions(vol));
    printf("nc_data_type = %d\n", vol->nc_data_type);
  }

  if (status != OK) {
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "mincRead(): error reading "
                 "volume from file %s",
                 fname));
  }

  /* ----- check the number of dimensions ----- */
  ndims = get_volume_n_dimensions(vol);
  if (ndims != 3 && ndims != 4) {
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "mincRead(): %d dimensions "
                 "in file; expecting 3 or 4",
                 ndims));
  }

  /* ----- get the dimension sizes ----- */
  get_volume_sizes(vol, dim_sizes);

  /* --- one time point if there are only three dimensions in the file --- */
  if (ndims == 3) dim_sizes[3] = 1;

  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
    printf("DimSizes: %d, %d, %d, %d, %d\n", ndims, dim_sizes[0], dim_sizes[1], dim_sizes[2], dim_sizes[3]);
  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON) printf("DataType: %d\n", vol->nc_data_type);

  dtype = get_volume_nc_data_type(vol, &sflag);
  dtype = orderIntBytes(vol->nc_data_type);

  /* ----- get the data type ----- */
  if (vol->nc_data_type == NC_BYTE || vol->nc_data_type == NC_CHAR)
    dtype = MRI_UCHAR;
  else if (vol->nc_data_type == NC_SHORT)
    dtype = MRI_SHORT;
  else if (vol->nc_data_type == NC_LONG)
    dtype = MRI_LONG;
  else if (vol->nc_data_type == NC_FLOAT || vol->nc_data_type == NC_DOUBLE)
    dtype = MRI_FLOAT;
  else {
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "mincRead(): bad data type "
                 "(%d) in input file %s",
                 vol->nc_data_type,
                 fname));
  }

  /* ----- allocate the mri structure ----- */
  if (read_volume)
    mri = MRIallocSequence(dim_sizes[0], dim_sizes[1], dim_sizes[2], dtype, dim_sizes[3]);
  else {
    mri = MRIallocHeader(dim_sizes[0], dim_sizes[1], dim_sizes[2], dtype, dim_sizes[3]);
    mri->nframes = dim_sizes[3];
  }

  /* ----- set up the mri structure ----- */
  get_volume_separations(vol, separations);
  mri->xsize = fabs(separations[0]);
  mri->ysize = fabs(separations[1]);
  mri->zsize = fabs(separations[2]);
  mri->ps = mri->xsize;
  mri->thick = mri->zsize;

  mri->x_r = vol->direction_cosines[0][0];
  mri->x_a = vol->direction_cosines[0][1];
  mri->x_s = vol->direction_cosines[0][2];

  mri->y_r = vol->direction_cosines[1][0];
  mri->y_a = vol->direction_cosines[1][1];
  mri->y_s = vol->direction_cosines[1][2];

  mri->z_r = vol->direction_cosines[2][0];
  mri->z_a = vol->direction_cosines[2][1];
  mri->z_s = vol->direction_cosines[2][2];

  if (separations[0] < 0) {
    mri->x_r = -mri->x_r;
    mri->x_a = -mri->x_a;
    mri->x_s = -mri->x_s;
  }
  if (separations[1] < 0) {
    mri->y_r = -mri->y_r;
    mri->y_a = -mri->y_a;
    mri->y_s = -mri->y_s;
  }
  if (separations[2] < 0) {
    mri->z_r = -mri->z_r;
    mri->z_a = -mri->z_a;
    mri->z_s = -mri->z_s;
  }

  voxel[0] = (mri->width) / 2.0;
  voxel[1] = (mri->height) / 2.0;
  voxel[2] = (mri->depth) / 2.0;
  voxel[3] = 0.0;
  convert_voxel_to_world(vol, voxel, &worldr, &worlda, &worlds);
  mri->c_r = worldr;
  mri->c_a = worlda;
  mri->c_s = worlds;

  mri->ras_good_flag = 1;

  mri->xend = mri->xsize * mri->width / 2.0;
  mri->xstart = -mri->xend;
  mri->yend = mri->ysize * mri->height / 2.0;
  mri->ystart = -mri->yend;
  mri->zend = mri->zsize * mri->depth / 2.0;
  mri->zstart = -mri->zend;

  xfov = mri->xend - mri->xstart;
  yfov = mri->yend - mri->ystart;
  zfov = mri->zend - mri->zstart;

  mri->fov = (xfov > yfov ? (xfov > zfov ? xfov : zfov) : (yfov > zfov ? yfov : zfov));

  strcpy(mri->fname, fname);

  //////////////////////////////////////////////////////////////////////////
  // test transform their way and our way:
  // MRIvoxelToWorld(mri, voxel[0], voxel[1], voxel[2], &wx, &wy, &wz);
  // printf("MNI calculated %.2f, %.2f, %.2f
  //     vs. MRIvoxelToWorld: %.2f, %.2f, %.2f\n",
  //     worldr, worlda, worlds, wx, wy, wz);

  /* ----- copy the data from the file to the mri structure ----- */
  if (read_volume) {
    while (input_more_of_volume(vol, &input_info, &f))
      ;

    for (i = 0; i < mri->width; i++) {
      for (j = 0; j < mri->height; j++) {
        for (k = 0; k < mri->depth; k++) {
          for (t = 0; t < mri->nframes; t++) {
            val = get_volume_voxel_value(vol, i, j, k, t, 0);
            if (mri->type == MRI_UCHAR) MRIseq_vox(mri, i, j, k, t) = (unsigned char)val;
            if (mri->type == MRI_SHORT) MRISseq_vox(mri, i, j, k, t) = (short)val;
            if (mri->type == MRI_LONG) MRILseq_vox(mri, i, j, k, t) = (long)val;
            if (mri->type == MRI_FLOAT) MRIFseq_vox(mri, i, j, k, t) = (float)val;
          }
        }
      }
      exec_progress_callback(i, mri->width, 0, 1);
    }
  }

  pVox2WorldGen = get_voxel_to_world_transform(vol);
  pVox2WorldLin = get_linear_transform_ptr(pVox2WorldGen);
  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON) {
    printf("MINC Linear Transform\n");
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) printf("%7.4f ", pVox2WorldLin->m[j][i]);
      printf("\n");
    }
  }

  delete_volume_input(&input_info);
  delete_volume(vol);

  return (mri);

} /* end mincRead() */

/*----------------------------------------------------------*/
/* time course clean */
static int mincWrite(MRI *mri, const char *fname)
{
  Volume minc_volume;
  const char* dimension_names[4] = {"xspace", "yspace", "zspace", "time"};
  nc_type nc_data_type;
  double min, max;
  float fmin, fmax;
  int dimension_sizes[4];
  double separations[4];
  double dir_cos[4];
  int return_value;
  double voxel[4], world[4];
  int signed_flag;
  int di_x, di_y, di_z;
  int vi[4];
  /*   int r, a, s; */
  /*   float r_max; */
  VIO_Status status;

  /* di gives the bogus minc index
     di[0] is width, 1 is height, 2 is depth, 3 is time if
     there is a time dimension of length > 1
     minc wants these to be ras */

  /* here: minc volume index 0 is r, 1 is a, 2 is s */
  /* mri is lia */ /* (not true for some volumes *)*/

  /* ----- get the orientation of the volume ----- */
  if (mri->ras_good_flag == 0) {
    setDirectionCosine(mri, MRI_CORONAL);
  }

  /* we don't muck around axes anymore */
  di_x = 0;
  di_y = 1;
  di_z = 2;


  /* ----- set the data type ----- */
  if (mri->type == MRI_UCHAR) {
    nc_data_type = NC_BYTE;
    signed_flag = 0;
  }
  else if (mri->type == MRI_SHORT) {
    nc_data_type = NC_SHORT;
    signed_flag = 1;
  }
  else if (mri->type == MRI_INT) {
    nc_data_type = NC_LONG;
    signed_flag = 1;
  }
  else if (mri->type == MRI_LONG) {
    nc_data_type = NC_LONG;
    signed_flag = 1;
  }
  else if (mri->type == MRI_FLOAT) {
    nc_data_type = NC_FLOAT;
    signed_flag = 1;
  }
  else {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "mincWrite(): bad data type (%d) in mri structure", mri->type));
  }

  if ((return_value = MRIlimits(mri, &fmin, &fmax)) != NO_ERROR) return (return_value);

  min = (double)fmin;
  max = (double)fmax;

  if (mri->nframes == 1)
    minc_volume = create_volume(3, const_cast<char**>(dimension_names), nc_data_type, signed_flag, min, max);
  else
    minc_volume = create_volume(4, const_cast<char**>(dimension_names), nc_data_type, signed_flag, min, max);

  /* di_(x,y,z) is the map from minc to orig */
  /* minc dimension size is that of di_x, etc. */
  dimension_sizes[di_x] = mri->width;
  dimension_sizes[di_y] = mri->height;
  dimension_sizes[di_z] = mri->depth;
  dimension_sizes[3] = mri->nframes;

  set_volume_sizes(minc_volume, dimension_sizes);

  alloc_volume_data(minc_volume);

  separations[di_x] = (double)(mri->xsize);
  separations[di_y] = (double)(mri->ysize);
  separations[di_z] = (double)(mri->zsize);
  separations[3] = 1.0;  // appears to do nothing
  set_volume_separations(minc_volume, separations);
  /* has side effect to change transform and thus must be set first */

  dir_cos[0] = (double)mri->x_r;
  dir_cos[1] = (double)mri->x_a;
  dir_cos[2] = (double)mri->x_s;
  set_volume_direction_cosine(minc_volume, di_x, dir_cos);

  dir_cos[0] = (double)mri->y_r;
  dir_cos[1] = (double)mri->y_a;
  dir_cos[2] = (double)mri->y_s;
  set_volume_direction_cosine(minc_volume, di_y, dir_cos);

  dir_cos[0] = (double)mri->z_r;
  dir_cos[1] = (double)mri->z_a;
  dir_cos[2] = (double)mri->z_s;
  set_volume_direction_cosine(minc_volume, di_z, dir_cos);

  voxel[di_x] = mri->width / 2.0;  // promoted to double
  voxel[di_y] = mri->height / 2.0;
  voxel[di_z] = mri->depth / 2.0;
  voxel[3] = 0.0;
  world[0] = (double)(mri->c_r);
  world[1] = (double)(mri->c_a);
  world[2] = (double)(mri->c_s);
  world[3] = 0.0;
  set_volume_translation(minc_volume, voxel, world);

  /* get the position from (vi[di_x], vi[di_y], vi[di_z]) orig position     */
  /*      put the value to (vi[0], vi[1], vi[2]) minc volume                */
  /* this will keep the orientation of axes the same as the original        */

  /* vi[n] gives the index of the variable along minc axis x */
  /* vi[di_x] gives the index of the variable
     along minc axis di_x, or along mri axis x */
  for (vi[3] = 0; vi[3] < mri->nframes; vi[3]++) {              /* frames */
    for (vi[di_x] = 0; vi[di_x] < mri->width; vi[di_x]++) {     /* columns */
      for (vi[di_y] = 0; vi[di_y] < mri->height; vi[di_y]++) {  /* rows */
        for (vi[di_z] = 0; vi[di_z] < mri->depth; vi[di_z]++) { /* slices */

	  double voxel;
	  switch (mri->type) {
	  case MRI_UCHAR: voxel = (double)MRIseq_vox (mri, vi[di_x], vi[di_y], vi[di_z], vi[3]); break;
	  case MRI_SHORT: voxel = (double)MRISseq_vox(mri, vi[di_x], vi[di_y], vi[di_z], vi[3]); break;
          case MRI_INT:   voxel = (double)MRIIseq_vox(mri, vi[di_x], vi[di_y], vi[di_z], vi[3]); break;
          case MRI_LONG:  voxel = (double)MRILseq_vox(mri, vi[di_x], vi[di_y], vi[di_z], vi[3]); break;
          case MRI_FLOAT: voxel = (double)MRIFseq_vox(mri, vi[di_x], vi[di_y], vi[di_z], vi[3]); break;
	  default: fprintf(stderr, "%s:%d bad type", __FILE__, __LINE__); exit(1);
          }

	  
	  set_volume_voxel_value(minc_volume,
                                   vi[0],
                                   vi[1],
                                   vi[2],
                                   vi[3],
                                   0,
                                   voxel);
	}
      }
      exec_progress_callback(vi[di_x], mri->width, vi[3], mri->nframes);
    }
  }


  status = output_volume((char*)fname, nc_data_type, signed_flag, min, max, minc_volume, "", NULL);
  delete_volume(minc_volume);

  if (status) {
    printf("ERROR: mincWrite: output_volume exited with %d\n", status);
    return (1);
  }

  return (NO_ERROR);

} /* end mincWrite() */

/*----------------------------------------------------------------
  bvolumeWrite() - replaces bshortWrite and bfloatWrite.
  Bug: fname_passed must the the stem, not the full file name.
  -----------------------------------------------------------------*/
static int bvolumeWrite(MRI *vol, const char *fname_passed, int type)
{
  int i, j, t;
  std::string fname;;
  short *bufshort;
  float *buffloat;
  FILE *fp;
  int result;
  MRI *subject_info = NULL;
  char subject_volume_dir[STRLEN];
  char *subjects_dir;
  char *sn;
  char analyse_fname[STRLEN], register_fname[STRLEN];
  char output_dir[STRLEN];
  char *c;
  int od_length;
  char t1_path[STRLEN];
  MATRIX *cf, *bf, *ibf, *af, *iaf, *as, *bs, *cs, *ics, *r;
  MATRIX *r1, *r2, *r3, *r4;
  float det;
  int bad_flag;
  int l;
  char stem[STRLEN];
  const char *c1, *c2, *c3;
  struct stat stat_buf;
  char subject_dir[STRLEN];
  int dealloc, nslices, nframes;
  MRI *mri;
  float min, max;
  int swap_bytes_flag, size, bufsize, endian = 0;
  const char *ext;
  void *buf;

  /* check the type and set the extension and size*/
  switch (type) {
    case MRI_SHORT:
      ext = "bshort";
      size = sizeof(short);
      break;
    case MRI_FLOAT:
      ext = "bfloat";
      size = sizeof(float);
      break;
    default:
      fprintf(stderr, "ERROR: bvolumeWrite: type (%d) is not short or float\n", type);
      return (1);
  }

  if (vol->type != type) {
    if (DIAG_VERBOSE_ON) printf("INFO: bvolumeWrite: changing type\n");
    nslices = vol->depth;
    nframes = vol->nframes;
    vol->depth = nslices * nframes;
    vol->nframes = 1;
    MRIlimits(vol, &min, &max);
    if (DIAG_VERBOSE_ON) printf("INFO: bvolumeWrite: range %g %g\n", min, max);
    mri = MRIchangeType(vol, type, min, max, 1);
    if (mri == NULL) {
      fprintf(stderr, "ERROR: bvolumeWrite: MRIchangeType\n");
      return (1);
    }
    vol->depth = nslices;
    vol->nframes = nframes;
    mri->depth = nslices;
    mri->nframes = nframes;
    dealloc = 1;
  }
  else {
    mri = vol;
    dealloc = 0;
  }

  /* ----- get the stem from the passed file name ----- */
  /*
    four options:
    1. stem_xxx.bshort
    2. stem.bshort
    3. stem_xxx
    4. stem
    other possibles:
    stem_.bshort
  */
  l = strlen(fname_passed);
  c1 = fname_passed + l - 11;
  c2 = fname_passed + l - 7;
  c3 = fname_passed + l - 4;
  strcpy(stem, fname_passed);
  if (c1 > fname_passed) {
    if ((*c1 == '_' && strcmp(c1 + 4, ".bshort") == 0) || (*c1 == '_' && strcmp(c1 + 4, ".bfloat") == 0))
      stem[(int)(c1 - fname_passed)] = '\0';
  }
  if (c2 > fname_passed) {
    if ((*c2 == '_' && strcmp(c2 + 4, ".bshort") == 0) || (*c2 == '_' && strcmp(c2 + 4, ".bfloat") == 0))
      stem[(int)(c2 - fname_passed)] = '\0';
  }
  if (c3 > fname_passed) {
    if (*c3 == '_') stem[(int)(c3 - fname_passed)] = '\0';
  }
  printf("INFO: bvolumeWrite: stem = %s\n", stem);

  c = strrchr(stem, '/');
  if (c == NULL)
    output_dir[0] = '\0';
  else {
    od_length = (int)(c - stem);
    strncpy(output_dir, stem, od_length);
    /* -- leaving the trailing '/' on a directory is not my
       usual convention, but here it's a load easier if there's
       no directory in stem... -ch -- */
    output_dir[od_length] = '/';
    output_dir[od_length + 1] = '\0';
  }

  int needed = snprintf(analyse_fname, STRLEN, "%s%s", output_dir, "analyse.dat");
  if( needed >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": analyse_fname truncated" << std::endl;
  }
  needed = snprintf(register_fname, STRLEN, "%s%s", output_dir, "register.dat");
  if( needed >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": register_fname truncated" << std::endl;
  }

  bufsize = mri->width * size;
  bufshort = (short *)malloc(bufsize);
  buffloat = (float *)malloc(bufsize);

  buf = NULL; /* shuts up compiler */
  if (type == MRI_SHORT)
    buf = bufshort;
  else
    buf = buffloat;

  swap_bytes_flag = 0;
#if (BYTE_ORDER == LITTLE_ENDIAN)
  swap_bytes_flag = 1;  // make it big endian
#endif

  if (getenv("BFILE_LITTLE_ENDIAN") != NULL) {
    // This is a mechanism to force bfloat/bshort files to be written as little endian
    printf("INFO: BFILE_LITTLE_ENDIAN is set, so writing as little endian\n");
    endian = 1;
    if (BYTE_ORDER == LITTLE_ENDIAN) swap_bytes_flag = 0;
    printf("swap = %d, order=%d, le=%d be=%d, endian=%d\n",
           swap_bytes_flag,
           BYTE_ORDER,
           LITTLE_ENDIAN,
           BIG_ENDIAN,
           endian);
  }

  // printf("--------------------------------\n");
  // MRIdump(mri,stdout);
  // MRIdumpBuffer(mri,stdout);
  // printf("--------------------------------\n");

  for (i = 0; i < mri->depth; i++) {
    /* ----- write the header file ----- */
    std::stringstream tmp;
    tmp << stem << '_' << std::setw(3) << std::setfill('0') << i << ".hdr";
    fname = tmp.str();
    if ((fp = fopen(fname.c_str(), "w")) == NULL) {
      if (dealloc) MRIfree(&mri);
      free(bufshort);
      free(buffloat);
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "bvolumeWrite(): can't open file %s", fname.c_str()));
    }
    fprintf(fp, "%d %d %d %d\n", mri->height, mri->width, mri->nframes, endian);
    fclose(fp);

    /* ----- write the data file ----- */
    tmp.str("");
    tmp << stem << '_' << std::setw(3) << std::setfill('0') << i << '.' << ext;
    fname = tmp.str();
    if ((fp = fopen(fname.c_str(), "w")) == NULL) {
      if (dealloc) MRIfree(&mri);
      free(bufshort);
      free(buffloat);
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "bvolumeWrite(): can't open file %s", fname.c_str()));
    }

    for (t = 0; t < mri->nframes; t++) {
      for (j = 0; j < mri->height; j++) {
        memmove(buf, mri->slices[t * mri->depth + i][j], bufsize);

        if (swap_bytes_flag) {
          if (type == MRI_SHORT)
            byteswapbufshort((void *)buf, bufsize);
          else
            byteswapbuffloat((void *)buf, bufsize);
        }
        fwrite(buf, size, mri->width, fp);
      }
    }

    fclose(fp);
    exec_progress_callback(i, mri->depth, 0, 1);
  }

  free(bufshort);
  free(buffloat);

  sn = subject_name;
  if (mri->subject_name[0] != '\0') sn = mri->subject_name;

  if (sn != NULL) {
    if ((subjects_dir = getenv("SUBJECTS_DIR")) == NULL) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "bvolumeWrite(): environment variable SUBJECTS_DIR unset");
      if (dealloc) MRIfree(&mri);
    }
    else {
      int req = snprintf(subject_dir, STRLEN, "%s/%s", subjects_dir, sn);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      if (stat(subject_dir, &stat_buf) < 0) {
        fprintf(stderr, "can't stat %s; writing to bhdr instead\n", subject_dir);
      }
      else {
        if (!S_ISDIR(stat_buf.st_mode)) {
          fprintf(stderr, "%s is not a directory; writing to bhdr instead\n", subject_dir);
        }
        else {
          int needed = snprintf(subject_volume_dir, STRLEN, "%s/mri/T1", subject_dir);
	  if( needed >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": subject_volume_dir truncated" << std::endl;
	  }
          subject_info = MRIreadInfo(subject_volume_dir);
          if (subject_info == NULL) {
            int needed = snprintf(subject_volume_dir, STRLEN, "%s/mri/orig", subject_dir);
	    if( needed >= STRLEN ) {
	      std::cerr << __FUNCTION__ << ": subject_volume_dir truncated" << std::endl;
	    }
            subject_info = MRIreadInfo(subject_volume_dir);
            if (subject_info == NULL) {
              fprintf(stderr,
                      "can't read the subject's orig or T1 volumes; "
                      "writing to bhdr instead\n");
	    }
          }
        }
      }
    }
  }

  if (subject_info != NULL) {
    if (subject_info->ras_good_flag == 0) {
      setDirectionCosine(subject_info, MRI_CORONAL);
    }
  }

  cf = bf = ibf = af = iaf = as = bs = cs = ics = r = NULL;
  r1 = r2 = r3 = r4 = NULL;

  /* ----- write the register.dat and analyse.dat  or bhdr files ----- */
  if (subject_info != NULL) {
    bad_flag = FALSE;

    if ((as = MatrixAlloc(4, 4, MATRIX_REAL)) == NULL) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "bvolumeWrite(): error creating matrix");
      bad_flag = TRUE;
    }
    stuff_four_by_four(as,
                       subject_info->x_r,
                       subject_info->y_r,
                       subject_info->z_r,
                       subject_info->c_r,
                       subject_info->y_r,
                       subject_info->y_r,
                       subject_info->y_r,
                       subject_info->c_r,
                       subject_info->z_r,
                       subject_info->z_r,
                       subject_info->z_r,
                       subject_info->c_r,
                       0.0,
                       0.0,
                       0.0,
                       1.0);

    if ((af = MatrixAlloc(4, 4, MATRIX_REAL)) == NULL) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "bvolumeWrite(): error creating matrix");
      bad_flag = TRUE;
    }
    stuff_four_by_four(af,
                       mri->x_r,
                       mri->y_r,
                       mri->z_r,
                       mri->c_r,
                       mri->y_r,
                       mri->y_r,
                       mri->y_r,
                       mri->c_r,
                       mri->z_r,
                       mri->z_r,
                       mri->z_r,
                       mri->c_r,
                       0.0,
                       0.0,
                       0.0,
                       1.0);

    if ((bs = MatrixAlloc(4, 4, MATRIX_REAL)) == NULL) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "bvolumeWrite(): error creating matrix");
      bad_flag = TRUE;
    }
    stuff_four_by_four(bs,
                       1,
                       0,
                       0,
                       (subject_info->width - 1) / 2.0,
                       0,
                       1,
                       0,
                       (subject_info->height - 1) / 2.0,
                       0,
                       0,
                       1,
                       (subject_info->depth - 1) / 2.0,
                       0,
                       0,
                       0,
                       1.0);

    if ((bf = MatrixAlloc(4, 4, MATRIX_REAL)) == NULL) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "bvolumeWrite(): error creating matrix");
      bad_flag = TRUE;
    }
    stuff_four_by_four(bf,
                       1,
                       0,
                       0,
                       (mri->width - 1) / 2.0,
                       0,
                       1,
                       0,
                       (mri->height - 1) / 2.0,
                       0,
                       0,
                       1,
                       (mri->depth - 1) / 2.0,
                       0,
                       0,
                       0,
                       1.0);

    if ((cs = MatrixAlloc(4, 4, MATRIX_REAL)) == NULL) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "bvolumeWrite(): error creating matrix");
      bad_flag = TRUE;
    }
    stuff_four_by_four(cs,
                       -subject_info->xsize,
                       0,
                       0,
                       (subject_info->width * mri->xsize) / 2.0,
                       0,
                       0,
                       subject_info->zsize,
                       -(subject_info->depth * mri->zsize) / 2.0,
                       0,
                       -subject_info->ysize,
                       0,
                       (subject_info->height * mri->ysize) / 2.0,
                       0,
                       0,
                       0,
                       1);

    if ((cf = MatrixAlloc(4, 4, MATRIX_REAL)) == NULL) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "bvolumeWrite(): error creating matrix");
      bad_flag = TRUE;
    }
    stuff_four_by_four(cf,
                       -mri->xsize,
                       0,
                       0,
                       (mri->width * mri->xsize) / 2.0,
                       0,
                       0,
                       mri->zsize,
                       -(mri->depth * mri->zsize) / 2.0,
                       0,
                       -mri->ysize,
                       0,
                       (mri->height * mri->ysize) / 2.0,
                       0,
                       0,
                       0,
                       1);

    if (bad_flag) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "bvolumeWrite(): error creating one "
                  "or more matrices; aborting register.dat "
                  "write and writing bhdr instead");
      MRIfree(&subject_info);
    }
  }

  if (subject_info != NULL) {
    bad_flag = FALSE;

    if ((det = MatrixDeterminant(as)) == 0.0) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "bvolumeWrite(): bad determinant in "
                  "matrix (check structural volume)");
      bad_flag = TRUE;
    }
    if ((det = MatrixDeterminant(bs)) == 0.0) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "bvolumeWrite(): bad determinant in "
                  "matrix (check structural volume)");
      bad_flag = TRUE;
    }
    if ((det = MatrixDeterminant(cs)) == 0.0) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "bvolumeWrite(): bad determinant in "
                  "matrix (check structural volume)");
      bad_flag = TRUE;
    }

    if ((det = MatrixDeterminant(af)) == 0.0) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "bvolumeWrite(): bad determinant in "
                  "matrix (check functional volume)");
      bad_flag = TRUE;
    }
    if ((det = MatrixDeterminant(bf)) == 0.0) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "bvolumeWrite(): bad determinant in "
                  "matrix (check functional volume)");
      bad_flag = TRUE;
    }
    if ((det = MatrixDeterminant(cf)) == 0.0) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "bvolumeWrite(): bad determinant in "
                  "matrix (check functional volume)");
      bad_flag = TRUE;
    }

    if (bad_flag) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "bvolumeWrite(): one or more zero determinants;"
                  " aborting register.dat write and writing bhdr instead");
      MRIfree(&subject_info);
    }
  }

  if (subject_info != NULL) {
    bad_flag = FALSE;

    if ((iaf = MatrixInverse(af, NULL)) == NULL) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "bvolumeWrite(): error inverting matrix");
      bad_flag = TRUE;
    }
    if ((ibf = MatrixInverse(bf, NULL)) == NULL) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "bvolumeWrite(): error inverting matrix");
      bad_flag = TRUE;
    }
    if ((ics = MatrixInverse(cs, NULL)) == NULL) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "bvolumeWrite(): error inverting matrix");
      bad_flag = TRUE;
    }

    if (bad_flag) {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM,
                  "bvolumeWrite(): one or more zero determinants; "
                  "aborting register.dat write and writing bhdr instead");
      MRIfree(&subject_info);
    }
  }

  bad_flag = FALSE;

  if (subject_info != NULL) {
    if ((r1 = MatrixMultiply(bs, ics, NULL)) == NULL) {
      bad_flag = TRUE;
      MRIfree(&subject_info);
    }
  }

  if (subject_info != NULL) {
    if ((r2 = MatrixMultiply(as, r1, NULL)) == NULL) {
      bad_flag = TRUE;
      MRIfree(&subject_info);
    }
  }

  if (subject_info != NULL) {
    if ((r3 = MatrixMultiply(iaf, r2, NULL)) == NULL) {
      bad_flag = TRUE;
      MRIfree(&subject_info);
    }
  }

  if (subject_info != NULL) {
    if ((r4 = MatrixMultiply(ibf, r3, NULL)) == NULL) {
      bad_flag = TRUE;
      MRIfree(&subject_info);
    }
  }

  if (subject_info != NULL) {
    if ((r = MatrixMultiply(cf, r4, NULL)) == NULL) {
      bad_flag = TRUE;
      MRIfree(&subject_info);
    }
  }

  if (bad_flag) {
    errno = 0;
    ErrorPrintf(ERROR_BADPARM,
                "bvolumeWrite(): error during matrix multiplications; "
                "aborting register.dat write and writing bhdr instead");
  }

  if (as != NULL) MatrixFree(&as);
  if (bs != NULL) MatrixFree(&bs);
  if (cs != NULL) MatrixFree(&cs);
  if (af != NULL) MatrixFree(&af);
  if (bf != NULL) MatrixFree(&bf);
  if (cf != NULL) MatrixFree(&cf);
  if (iaf != NULL) MatrixFree(&iaf);
  if (ibf != NULL) MatrixFree(&ibf);
  if (ics != NULL) MatrixFree(&ics);
  if (r1 != NULL) MatrixFree(&r1);
  if (r2 != NULL) MatrixFree(&r2);
  if (r3 != NULL) MatrixFree(&r3);
  if (r4 != NULL) MatrixFree(&r4);

  if (subject_info != NULL) {
    if (mri->path_to_t1[0] == '\0')
      sprintf(t1_path, ".");
    else
      strcpy(t1_path, mri->path_to_t1);

    if (FileExists(analyse_fname)) fprintf(stderr, "warning: overwriting file %s\n", analyse_fname);

    if ((fp = fopen(analyse_fname, "w")) == NULL) {
      MRIfree(&subject_info);
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "bvolumeWrite(): couldn't open file %s for writing", analyse_fname));
    }

    fprintf(fp, "%s\n", t1_path);
    fprintf(fp, "%s_%%03d.bshort\n", stem);
    fprintf(fp, "%d %d\n", mri->depth, mri->nframes);
    fprintf(fp, "%d %d\n", mri->width, mri->height);

    fclose(fp);

    if (FileExists(analyse_fname)) fprintf(stderr, "warning: overwriting file %s\n", register_fname);

    if ((fp = fopen(register_fname, "w")) == NULL) {
      MRIfree(&subject_info);
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "bvolumeWrite(): couldn't open file %s for writing", register_fname));
    }

    fprintf(fp, "%s\n", sn);
    fprintf(fp, "%g\n", mri->xsize);
    fprintf(fp, "%g\n", mri->zsize);
    fprintf(fp, "%g\n", 1.0);
    fprintf(fp,
            "%g %g %g %g\n",
            *MATRIX_RELT(r, 1, 1),
            *MATRIX_RELT(r, 1, 2),
            *MATRIX_RELT(r, 1, 3),
            *MATRIX_RELT(r, 1, 4));
    fprintf(fp,
            "%g %g %g %g\n",
            *MATRIX_RELT(r, 2, 1),
            *MATRIX_RELT(r, 2, 2),
            *MATRIX_RELT(r, 2, 3),
            *MATRIX_RELT(r, 2, 4));
    fprintf(fp,
            "%g %g %g %g\n",
            *MATRIX_RELT(r, 3, 1),
            *MATRIX_RELT(r, 3, 2),
            *MATRIX_RELT(r, 3, 3),
            *MATRIX_RELT(r, 3, 4));
    fprintf(fp,
            "%g %g %g %g\n",
            *MATRIX_RELT(r, 4, 1),
            *MATRIX_RELT(r, 4, 2),
            *MATRIX_RELT(r, 4, 3),
            *MATRIX_RELT(r, 4, 4));

    fclose(fp);

    MatrixFree(&r);
  }

  if (subject_info == NULL) {
    fname = std::string(stem) + ".hdr";
    if ((fp = fopen(fname.c_str(), "w")) == NULL) {
      if (dealloc) MRIfree(&mri);
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "bvolumeWrite(): can't open file %s", fname.c_str()));
    }

    result = write_bhdr(mri, fp);

    fclose(fp);

    if (result != NO_ERROR) return (result);
  }
  else
    MRIfree(&subject_info);

  if (dealloc) MRIfree(&mri);

  return (NO_ERROR);

} /* end bvolumeWrite() */

static MRI *get_b_info(const char *fname_passed, int read_volume, char *directory, char *stem, int type)
{
  MRI *mri, *mri2;
  FILE *fp;
  int nslices = 0, nt;
  int nx, ny, i;
  ;
  int result;
  char fname[STRLEN];
  char extension[STRLEN];
  char bhdr_name[STRLEN];

  if (type == MRI_SHORT)
    sprintf(extension, "bshort");
  else if (type == MRI_FLOAT)
    sprintf(extension, "bfloat");
  else {
    errno = 0;
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "internal error: get_b_info() passed type %d", type));
  }

  result = decompose_b_fname(fname_passed, directory, stem);
  if (result != NO_ERROR) return (NULL);

  if (directory[0] == '\0') sprintf(directory, ".");


  mri = MRIallocHeader(1, 1, 1, type, 1);

  /* ----- try to read the stem.bhdr ----- */
  int required = snprintf(bhdr_name, STRLEN, "%s/%s.bhdr", directory, stem);
  if( required >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": bhdr_name truncated" << std::endl;
  }
  if ((fp = fopen(bhdr_name, "r")) != NULL) {
    read_bhdr(mri, fp);
    int required = snprintf(fname, STRLEN, "%s/%s_000.hdr", directory, stem);
    if( required >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": fname truncated" << std::endl;
    }
    if ((fp = fopen(fname, "r")) == NULL) {
      MRIfree(&mri);
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADFILE, "cannot open %s", fname));
    }
    if (fscanf(fp, "%d %d %d %*d", &ny, &nx, &nt) != 3){
      ErrorPrintf(ERROR_BADFILE, "bvolumeRead(): could not read file"); 
    }
    mri->nframes = nt;
    free(mri->frames);
    mri->frames = (MRI_FRAME *)calloc(nt, sizeof(MRI_FRAME));
    if (mri->frames == NULL) ErrorExit(ERROR_NOMEMORY, "get_b_info: could not allocate %d frames", nt);
    fclose(fp);
    {
      int i;
      for (i = 0; i < nt; i++) mri->frames->m_ras2vox = MatrixIdentity(4, NULL);
    }

    strcpy(mri->fname, fname_passed);
  }
  else {
    /* ----- get defaults ----- */
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr,
              "-----------------------------------------------------------------\n"
              "Could not find the direction cosine information.\n"
              "Will use the CORONAL orientation.\n"
              "If not suitable, please provide the information in %s file\n"
              "-----------------------------------------------------------------\n",
              bhdr_name);
    int required = snprintf(fname, STRLEN, "%s/%s_000.hdr", directory, stem);
    if( required >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if ((fp = fopen(fname, "r")) == NULL) {
      MRIfree(&mri);
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "can't find file %s (last resort);"
                   "bailing out on read",
                   fname));
    }
    if (fscanf(fp, "%d %d %d %*d", &ny, &nx, &nt) != 3){
      ErrorPrintf(ERROR_BADFILE, "get_b_info(): could not read file");
    }
    fclose(fp);

    /* --- get the number of slices --- */
    int req = snprintf(fname, STRLEN, "%s/%s_000.%s", directory, stem, extension);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (!FileExists(fname)) {
      MRIfree(&mri);
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "can't find file %s; "
                   "bailing out on read",
                   fname));
    }
    for (nslices = 0; FileExists(fname); nslices++) {
      int req = snprintf(fname, STRLEN, "%s/%s_%03d.%s", directory, stem, nslices, extension);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    nslices--;

    mri->width = nx;
    mri->height = ny;
    mri->depth = nslices;
    mri->nframes = nt;

    mri->imnr0 = 1;
    mri->imnr1 = nslices;

    mri->thick = mri->ps = 1.0;
    mri->xsize = mri->ysize = mri->zsize = 1.0;

    setDirectionCosine(mri, MRI_CORONAL);

    mri->ras_good_flag = 0;

    strcpy(mri->fname, fname_passed);
  }

  mri->imnr0 = 1;
  mri->imnr1 = nslices;

  mri->xend = mri->width * mri->xsize / 2.0;
  mri->xstart = -mri->xend;
  mri->yend = mri->height * mri->ysize / 2.0;
  mri->ystart = -mri->yend;
  mri->zend = mri->depth * mri->zsize / 2.0;
  mri->zstart = -mri->zend;

  mri->fov =
      ((mri->xend - mri->xstart) > (mri->yend - mri->ystart) ? (mri->xend - mri->xstart) : (mri->yend - mri->ystart));

  mri->frames = (MRI_FRAME *)calloc(mri->nframes, sizeof(MRI_FRAME));
  for (i = 0; i < mri->nframes; i++) mri->frames[i].m_ras2vox = NULL;

  if (read_volume) {
    mri2 = MRIallocSequence(mri->width, mri->height, mri->depth, mri->type, mri->nframes);
    MRIcopyHeader(mri, mri2);
    MRIfree(&mri);
    mri = mri2;
  }

  return (mri);

} /* end get_b_info() */

/*-------------------------------------------------------------------
  bvolumeRead() - this replaces bshortRead and bfloatRead.
  -------------------------------------------------------------------*/
static MRI *bvolumeRead(const char *fname_passed, int read_volume, int type)
{
  MRI *mri;
  FILE *fp;
  std::string fname;
  char directory[STRLEN];
  char stem[STRLEN];
  int swap_bytes_flag;
  int slice, frame, row, k;
  int nread;
  const char *ext;
  int size;
  float min, max;

  /* check the type and set the extension and size*/
  switch (type) {
    case MRI_SHORT:
      ext = "bshort";
      size = sizeof(short);
      break;
    case MRI_FLOAT:
      ext = "bfloat";
      size = sizeof(float);
      break;
    default:
      fprintf(stderr,
              "ERROR: bvolumeRead: type (%d) is not "
              "short or float\n",
              type);
      return (NULL);
  }

  /* Get the header info (also allocs if needed) */
  mri = get_b_info(fname_passed, read_volume, directory, stem, type);
  if (mri == NULL) return (NULL);

  /* If not reading the volume, return now */
  if (!read_volume) return (mri);

  /* Read in the header of the first slice to get the endianness */
  std::stringstream tmp;
  tmp << directory << '/' << stem << '_'
      << std::setw(3) << std::setfill('0') << 0
      << ".hdr";
  fname = tmp.str();
  if ((fp = fopen(fname.c_str(), "r")) == NULL) {
    fprintf(stderr, "ERROR: can't open file %s; assuming big-endian bvolume\n", fname.c_str());
    swap_bytes_flag = 0;
  }
  else {
    if (fscanf(fp, "%*d %*d %*d %d", &swap_bytes_flag) != 1){
      ErrorPrintf(ERROR_BADFILE, "bvolumeRead(): could not read file");
    }
#if (BYTE_ORDER == LITTLE_ENDIAN)
    swap_bytes_flag = !swap_bytes_flag;
#endif
    fclose(fp);
  }

  /* Go through each slice */
  for (slice = 0; slice < mri->depth; slice++) {
    /* Open the file for this slice */ 
    std::stringstream tmp;
    tmp << directory << '/' << stem << '_'
	<< std::setw(3) << std::setfill('0') << slice
	<< '.' << ext;
    fname = tmp.str();
    if ((fp = fopen(fname.c_str(), "r")) == NULL) {
      MRIfree(&mri);
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADFILE, "bvolumeRead(): error opening file %s", fname.c_str()));
    }
    // fprintf(stderr, "Reading %s ... \n", fname);
    /* Loop through the frames */
    for (frame = 0; frame < mri->nframes; frame++) {
      k = slice + mri->depth * frame;

      /* Loop through the rows */
      for (row = 0; row < mri->height; row++) {
        /* read in all the columns for a row */
        nread = fread(mri->slices[k][row], size, mri->width, fp);
        if (nread != mri->width) {
          fclose(fp);
          MRIfree(&mri);
          errno = 0;
          ErrorReturn(NULL,
                      (ERROR_BADFILE,
                       "bvolumeRead(): "
                       "error reading from file %s",
                       fname.c_str()));
        }

        if (swap_bytes_flag) {
          if (type == MRI_SHORT) {
	    std::vector<short> temp(mri->width);
	    // Note:
	    // void swab(const void *from, void *to, ssize_t n);
	    // void *memcpy(void *dest, const void *src, size_t n);
	    // Because consistency is the hobgoblin of small minds...
	    swab(mri->slices[k][row], temp.data(), (size_t)(mri->width * size));
	    memcpy(mri->slices[k][row], temp.data(), (size_t)(mri->width * size));
          } else {
            byteswapbuffloat((void *)mri->slices[k][row], size * mri->width);
	  }
        }
      } /* row loop */
    }   /* frame loop */

    fclose(fp);
    exec_progress_callback(slice, mri->depth, 0, 1);
  } /* slice loop */

  MRIlimits(mri, &min, &max);
  printf("INFO: bvolumeRead: min = %g, max = %g\n", min, max);

  mri->imnr0 = 1;
  mri->imnr1 = mri->depth;
  mri->thick = mri->zsize;

  return (mri);

} /* end bvolumeRead() */


int decompose_b_fname(const char *fname, char *dir, char *stem)
{
  char *slash, *dot, *stem_start, *underscore;
  int fname_length;
  int und_pos;
  char fname_copy[STRLEN];

  /*

  options for file names:

  stem_xxx.bshort
  stem_xxx
  stem_.bshort
  stem_
  stem.bshort
  stem

  with or without a preceding directory

  */

  strcpy(fname_copy, fname);

  slash = strrchr(fname_copy, '/');
  dot = strrchr(fname_copy, '.');
  underscore = strrchr(fname_copy, '_');

  if (slash == NULL) {
    stem_start = fname_copy;
    sprintf(dir, ".");
  }
  else {
    *slash = '\0';
    strcpy(dir, fname_copy);
    stem_start = slash + 1;
  }

  if (*stem_start == '\0') {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "decompose_b_fname(): bad bshort/bfloat file specifier \"%s\"", fname));
  }

  if (dot != NULL)
    if (strcmp(dot, ".bshort") == 0 || strcmp(dot, ".bfloat") == 0) *dot = '\0';

  fname_length = strlen(stem_start);

  if (underscore == NULL)
    strcpy(stem, stem_start);
  else {
    und_pos = (underscore - stem_start);
    if (und_pos == fname_length - 1 || und_pos == fname_length - 4) *underscore = '\0';
    strcpy(stem, stem_start);
  }

  return (NO_ERROR);

} /* end decompose_b_fname() */

/*-------------------------------------------------------------*/
static int write_bhdr(MRI *mri, FILE *fp)
{
  float vl;            /* vector length */
  float tlr, tla, tls; /* top left coordinates */
  float trr, tra, trs; /* top right coordinates */
  float brr, bra, brs; /* bottom right coordinates */
  float nr, na, ns;    /* normal coordinates */
  MATRIX *T, *crs1, *ras1;

  crs1 = MatrixAlloc(4, 1, MATRIX_REAL);
  crs1->rptr[1][4] = 1;
  ras1 = MatrixAlloc(4, 1, MATRIX_REAL);

  /* Construct the matrix to convert CRS to XYZ, assuming
     that CRS is 1-based */
  T = MRIxfmCRS2XYZ(mri, 0);
  // printf("------- write_bhdr: T ---------\n");
  // MatrixPrint(stdout,T);
  // printf("------------------------------\n");

  /* The "top left" is the center of the first voxel */
  tlr = T->rptr[1][4];
  tla = T->rptr[2][4];
  tls = T->rptr[3][4];

  /* The "top right" is the center of the last voxel + 1 in the
     column direction - this unfortunate situation is historical.
     It makes the difference between the TL and TR equal to the
     edge-to-edge FOV in the column direction. */
  crs1->rptr[1][1] = mri->width;
  crs1->rptr[1][2] = 0;
  crs1->rptr[1][3] = 0;
  MatrixMultiply(T, crs1, ras1);
  trr = ras1->rptr[1][1];
  tra = ras1->rptr[1][2];
  trs = ras1->rptr[1][3];

  /* The "bottom right" is the center of the last voxel + 1 in the
     column and row directions - this unfortunate situation is
     historical. It makes the difference between the TR and BR equal
     to the edge-to-edge FOV in the row direction. */
  crs1->rptr[1][1] = mri->width;
  crs1->rptr[1][2] = mri->height;
  crs1->rptr[1][3] = 0;
  MatrixMultiply(T, crs1, ras1);
  brr = ras1->rptr[1][1];
  bra = ras1->rptr[1][2];
  brs = ras1->rptr[1][3];

  MatrixFree(&T);
  MatrixFree(&crs1);
  MatrixFree(&ras1);

  /* ----- normalize this just in case ----- */
  vl = sqrt(mri->z_r * mri->z_r + mri->z_a * mri->z_a + mri->z_s * mri->z_s);
  nr = mri->z_r / vl;
  na = mri->z_a / vl;
  ns = mri->z_s / vl;

  fprintf(fp, "          cols: %d\n", mri->width);
  fprintf(fp, "          rows: %d\n", mri->height);
  fprintf(fp, "       nslices: %d\n", mri->depth);
  fprintf(fp, " n_time_points: %d\n", mri->nframes);
  fprintf(fp, "   slice_thick: %f\n", mri->zsize);
  fprintf(fp, "    top_left_r: %f\n", tlr);
  fprintf(fp, "    top_left_a: %f\n", tla);
  fprintf(fp, "    top_left_s: %f\n", tls);
  fprintf(fp, "   top_right_r: %f\n", trr);
  fprintf(fp, "   top_right_a: %f\n", tra);
  fprintf(fp, "   top_right_s: %f\n", trs);
  fprintf(fp, "bottom_right_r: %f\n", brr);
  fprintf(fp, "bottom_right_a: %f\n", bra);
  fprintf(fp, "bottom_right_s: %f\n", brs);
  fprintf(fp, "      normal_r: %f\n", nr);
  fprintf(fp, "      normal_a: %f\n", na);
  fprintf(fp, "      normal_s: %f\n", ns);
  fprintf(fp, "      image_te: %f\n", mri->te);
  fprintf(fp, "      image_tr: %f\n", mri->tr / 1000.0);  // convert to sec
  fprintf(fp, "      image_ti: %f\n", mri->ti);
  fprintf(fp, "    flip_angle: %f\n", mri->flip_angle);

  return (NO_ERROR);

} /* end write_bhdr() */

/*------------------------------------------------------*/
int read_bhdr(MRI *mri, FILE *fp)
{
  char line[STRLEN];
  char *l;
  float tlr = 0.;
  float tla = 0.;
  float tls = 0.; /* top left coordinates */
  float trr = 0.;
  float tra = 0.;
  float trs = 0.; /* top right coordinates */
  float brr = 0.;
  float bra = 0.;
  float brs = 0.; /* bottom right coordinates */
  float xr = 0.;
  float xa = 0.;
  float xs = 0.;
  float yr = 0.;
  float ya = 0.;
  float ys = 0.;
  MATRIX *T, *CRSCenter, *RASCenter;

  while (1) {  // don't use   "while (!feof(fp))"

    /* --- read the line --- */
    if (!fgets(line, STRLEN, fp) && ferror(fp)){
      ErrorPrintf(ERROR_BADFILE, "read_bhdr(): could not read file");
    }

    if (feof(fp))  // wow, it was beyound the end of the file. get out.
      break;

    /* --- remove the newline --- */
    if (line[strlen(line) - 1] == '\n') line[strlen(line) - 1] = '\0';

    /* --- skip the initial spaces --- */
    for (l = line; isspace((int)(*l)); l++)
      ;

    /* --- get the varible name and value(s) --- */
    if (strlen(l) > 0) {
      if (strncmp(l, "cols: ", 6) == 0)
        sscanf(l, "%*s %d", &mri->width);
      else if (strncmp(l, "rows: ", 6) == 0)
        sscanf(l, "%*s %d", &mri->height);
      else if (strncmp(l, "nslices: ", 9) == 0)
        sscanf(l, "%*s %d", &mri->depth);
      else if (strncmp(l, "n_time_points: ", 15) == 0)
        sscanf(l, "%*s %d", &mri->nframes);
      else if (strncmp(l, "slice_thick: ", 13) == 0)
        sscanf(l, "%*s %f", &mri->zsize);
      else if (strncmp(l, "image_te: ", 10) == 0)
        sscanf(l, "%*s %f", &mri->te);
      else if (strncmp(l, "image_tr: ", 10) == 0) {
        sscanf(l, "%*s %f", &mri->tr);
        mri->tr = 1000.0 * mri->tr;  // convert from sec to msec
      }
      else if (strncmp(l, "image_ti: ", 10) == 0)
        sscanf(l, "%*s %f", &mri->ti);
      else if (strncmp(l, "flip_angle: ", 10) == 0)
        sscanf(l, "%*s %lf", &mri->flip_angle);
      else if (strncmp(l, "top_left_r: ", 12) == 0)
        sscanf(l, "%*s %g", &tlr);
      else if (strncmp(l, "top_left_a: ", 12) == 0)
        sscanf(l, "%*s %g", &tla);
      else if (strncmp(l, "top_left_s: ", 12) == 0)
        sscanf(l, "%*s %g", &tls);
      else if (strncmp(l, "top_right_r: ", 13) == 0)
        sscanf(l, "%*s %g", &trr);
      else if (strncmp(l, "top_right_a: ", 13) == 0)
        sscanf(l, "%*s %g", &tra);
      else if (strncmp(l, "top_right_s: ", 13) == 0)
        sscanf(l, "%*s %g", &trs);
      else if (strncmp(l, "bottom_right_r: ", 16) == 0)
        sscanf(l, "%*s %g", &brr);
      else if (strncmp(l, "bottom_right_a: ", 16) == 0)
        sscanf(l, "%*s %g", &bra);
      else if (strncmp(l, "bottom_right_s: ", 16) == 0)
        sscanf(l, "%*s %g", &brs);
      else if (strncmp(l, "normal_r: ", 10) == 0)
        sscanf(l, "%*s %g", &mri->z_r);
      else if (strncmp(l, "normal_a: ", 10) == 0)
        sscanf(l, "%*s %g", &mri->z_a);
      else if (strncmp(l, "normal_s: ", 10) == 0)
        sscanf(l, "%*s %g", &mri->z_s);
      else { /* --- ignore it --- */
      }
    }
  } /* end while(!feof()) */
  // forget to close file
  fclose(fp);

  //  vox to ras matrix is
  //
  //  Xr*Sx  Yr*Sy  Zr*Sz  Tr
  //  Xa*Sx  Ya*Sy  Za*Sz  Ta
  //  Xs*Sx  Ys*Sy  Zs*Sz  Ts
  //    0      0      0    1
  //
  //  Therefore
  //
  //  trr = Xr*Sx*W + Tr, tlr = Tr
  //  tra = Xa*Sx*W + Ta, tla = Ta
  //  trs = Xs*Sx*W + Ts, tls = Ts
  //
  //  Thus
  //
  //  Sx = sqrt ( ((trr-tlr)/W)^2 + ((tra-tla)/W)^2 + ((trs-tls)/W)^2)
  //     since Xr^2 + Xa^2 + Xs^2 = 1
  //  Xr = (trr-tlr)/(W*Sx)
  //  Xa = (tra-tla)/(W*Sx)
  //  Xs = (trs-tls)/(W*Sx)
  //
  //  Similar things for others
  //
  xr = (trr - tlr) / mri->width;
  xa = (tra - tla) / mri->width;
  xs = (trs - tls) / mri->width;
  mri->xsize = sqrt(xr * xr + xa * xa + xs * xs);
  if (mri->xsize)  // avoid nan
  {
    mri->x_r = xr / mri->xsize;
    mri->x_a = xa / mri->xsize;
    mri->x_s = xs / mri->xsize;
  }
  else  // fake values
  {
    mri->xsize = 1;
    mri->x_r = -1;
    mri->x_a = 0;
    mri->x_s = 0;
  }
  yr = (brr - trr) / mri->height;
  ya = (bra - tra) / mri->height;
  ys = (brs - trs) / mri->height;
  mri->ysize = sqrt(yr * yr + ya * ya + ys * ys);
  if (mri->ysize)  // avoid nan
  {
    mri->y_r = yr / mri->ysize;
    mri->y_a = ya / mri->ysize;
    mri->y_s = ys / mri->ysize;
  }
  else  // fake values
  {
    mri->ysize = 1;
    mri->y_r = 0;
    mri->y_a = 0;
    mri->y_s = -1;
  }
  T = MRIxfmCRS2XYZ(mri, 0);

  T->rptr[1][4] = tlr;
  T->rptr[2][4] = tla;
  T->rptr[3][4] = tls;

  // printf("------- read_bhdr: T ---------\n");
  // MatrixPrint(stdout,T);
  // printf("------------------------------\n");

  CRSCenter = MatrixAlloc(4, 1, MATRIX_REAL);
  CRSCenter->rptr[1][1] = (mri->width) / 2.0;
  CRSCenter->rptr[2][1] = (mri->height) / 2.0;
  CRSCenter->rptr[3][1] = (mri->depth) / 2.0;
  CRSCenter->rptr[4][1] = 1;

  RASCenter = MatrixMultiply(T, CRSCenter, NULL);
  mri->c_r = RASCenter->rptr[1][1];
  mri->c_a = RASCenter->rptr[2][1];
  mri->c_s = RASCenter->rptr[3][1];

  MatrixFree(&T);
  MatrixFree(&CRSCenter);
  MatrixFree(&RASCenter);


  mri->ras_good_flag = 1;

  mri->thick = mri->zsize;
  mri->ps = mri->xsize;

  return (NO_ERROR);

} /* end read_bhdr() */

static MRI *genesisRead(const char *fname, int read_volume)
{
  char fname_format[STRLEN];
  char fname_format2[STRLEN];
  char fname_dir[STRLEN];
  char fname_base[STRLEN];
  MRI *mri = NULL;
  int im_init;
  int im_low, im_high;
  int im_low2, im_high2;
  char fname_use[STRLEN];
  char temp_string[STRLEN];
  FILE *fp;
  int width, height;
  int pixel_data_offset;
  int image_header_offset;
  float tl_r, tl_a, tl_s;
  float tr_r, tr_a, tr_s;
  float br_r, br_a, br_s;
  float c_r, c_a, c_s;
  float n_r, n_a, n_s;
  float xlength, ylength, zlength;
  int i, y;
  MRI *header = nullptr;
  float xfov, yfov, zfov;
  float nlength;
  int twoformats = 0, odd_only, even_only;

  odd_only = even_only = 0;
  if (getenv("GE_ODD")) {
    odd_only = 1;
    printf("only using odd # GE files\n");
  }
  else if (getenv("GE_EVEN")) {
    even_only = 1;
    printf("only using even # GE files\n");
  }

  /* ----- check the first (passed) file ----- */
  if (!FileExists(fname)) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error opening file %s", fname));
  }

  /* ----- split the file name into name and directory ----- */
  const char *cc = strrchr(fname, '/');
  if (cc == NULL) {
    fname_dir[0] = '\0';
    strcpy(fname_base, fname);
  }
  else {
    strncpy(fname_dir, fname, (cc - fname + 1));
    fname_dir[cc - fname + 1] = '\0';
    strcpy(fname_base, cc + 1);
  }

  /* ----- derive the file name format (for sprintf) ----- */
  // this one fix fname_format only
  if (strncmp(fname_base, "I.", 2) == 0) {
    twoformats = 0;
    im_init = atoi(&fname_base[2]);
    sprintf(fname_format, "I.%%03d");
  }
  // this one fix both fname_format and fname_format2
  else if (strlen(fname_base) >= 3) /* avoid core dumps below... */
  {
    twoformats = 1;
    char *c = &fname_base[strlen(fname_base) - 3];
    if (strcmp(c, ".MR") == 0) {
      *c = '\0';
      for (c--; isdigit(*c) && c >= fname_base; c--)
        ;
      c++;
      im_init = atoi(c);
      *c = '\0';
      // this is too quick to assume of this type
      // another type %s%%03d.MR" must be examined
      int req = snprintf(fname_format, STRLEN, "%s%%d.MR", fname_base);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      req = snprintf(fname_format2, STRLEN, "%s%%03d.MR", fname_base);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    else {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADPARM, "genesisRead(): can't determine file name format for %s", fname));
    }
  }
  else {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "genesisRead(): can't determine file name format for %s", fname));
  }

  if (strlen(fname_format) != 0) {
    strcpy(temp_string, fname_format);
    int req = snprintf(fname_format, STRLEN, "%s%s", fname_dir, temp_string);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("fname_format  : %s\n", fname_format);
  }
  if (strlen(fname_format2) != 0) {
    strcpy(temp_string, fname_format2);
    int req = snprintf(fname_format2, STRLEN, "%s%s", fname_dir, temp_string);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("fname_format2 : %s\n", fname_format2);
  }
  /* ----- find the low and high files ----- */
  if (odd_only || even_only) {
    if ((odd_only && ISEVEN(im_init)) || (even_only && ISODD(im_init))) im_init++;
    im_low = im_init;
    do {
      im_low -= 2;
      sprintf(fname_use, fname_format, im_low);
    } while (FileExists(fname_use));
    im_low += 2;

    im_high = im_init;
    do {
      im_high += 2;
      sprintf(fname_use, fname_format, im_high);
    } while (FileExists(fname_use));
    im_high -= 2;
  }
  else {
    im_low = im_init;
    do {
      im_low--;
      sprintf(fname_use, fname_format, im_low);
    } while (FileExists(fname_use));
    im_low++;

    im_high = im_init;
    do {
      im_high++;
      sprintf(fname_use, fname_format, im_high);
    } while (FileExists(fname_use));
    im_high--;
  }

  if (twoformats) {
    // now test fname_format2
    im_low2 = im_init;
    do {
      im_low2--;
      sprintf(fname_use, fname_format2, im_low2);
    } while (FileExists(fname_use));
    im_low2++;

    im_high2 = im_init;
    do {
      im_high2++;
      sprintf(fname_use, fname_format2, im_high2);
    } while (FileExists(fname_use));
    im_high2--;
  }
  else {
    im_high2 = im_low2 = 0;
  }
  // now decide which one to pick
  if ((im_high2 - im_low2) > (im_high - im_low)) {
    // we have to use fname_format2
    strcpy(fname_format, fname_format2);
    im_high = im_high2;
    im_low = im_low2;
  }
  // otherwise the same

  /* ----- allocate the mri structure ----- */
  header = MRIallocHeader(1, 1, 1, MRI_SHORT, header->nframes);

  if (odd_only || even_only)
    header->depth = (im_high - im_low) / 2 + 1;
  else
    header->depth = im_high - im_low + 1;
  header->imnr0 = 1;
  header->imnr1 = header->depth;

  /* ----- get the header information from the first file ----- */
  sprintf(fname_use, fname_format, im_low);
  if ((fp = fopen(fname_use, "r")) == NULL) {
    MRIfree(&header);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error opening file %s\n", fname_use));
  }

  fseek(fp, 8, SEEK_SET);
  if (fread(&width, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  width = orderIntBytes(width);
  if (fread(&height, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  height = orderIntBytes(height);
  fseek(fp, 148, SEEK_SET);
  if (fread(&image_header_offset, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  image_header_offset = orderIntBytes(image_header_offset);

  header->width = width;
  header->height = height;

  strcpy(header->fname, fname);

  fseek(fp, image_header_offset + 26, SEEK_SET);
  if (fread(&(header->thick), 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  header->thick = orderFloatBytes(header->thick);
  header->zsize = header->thick;

  fseek(fp, image_header_offset + 50, SEEK_SET);
  if (fread(&(header->xsize), 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  header->xsize = orderFloatBytes(header->xsize);
  if (fread(&(header->ysize), 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  header->ysize = orderFloatBytes(header->ysize);
  header->ps = header->xsize;

/* all in micro-seconds */
#define MICROSECONDS_PER_MILLISECOND 1e3
  fseek(fp, image_header_offset + 194, SEEK_SET);
  header->tr = freadInt(fp) / MICROSECONDS_PER_MILLISECOND;
  fseek(fp, image_header_offset + 198, SEEK_SET);
  header->ti = freadInt(fp) / MICROSECONDS_PER_MILLISECOND;
  fseek(fp, image_header_offset + 202, SEEK_SET);
  header->te = freadInt(fp) / MICROSECONDS_PER_MILLISECOND;
  fseek(fp, image_header_offset + 254, SEEK_SET);
  header->flip_angle = RADIANS(freadShort(fp));    /* was in degrees */
  fseek(fp, image_header_offset + 210, SEEK_SET);  // # of echoes
  header->nframes = freadShort(fp);
  if (header->nframes > 1) {
    printf("multi-echo genesis file detected (%d echoes)...\n", header->nframes);
  }
  else if (header->nframes == 0) {
    printf("zero frames specified in file - setting to 1\n");
    header->nframes = 1;
  }

  fseek(fp, image_header_offset + 130, SEEK_SET);
  if (fread(&c_r, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  c_r = orderFloatBytes(c_r);
  if (fread(&c_a, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  c_a = orderFloatBytes(c_a);
  if (fread(&c_s, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  c_s = orderFloatBytes(c_s);
  if (fread(&n_r, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  n_r = orderFloatBytes(n_r);
  if (fread(&n_a, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  n_a = orderFloatBytes(n_a);
  if (fread(&n_s, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  n_s = orderFloatBytes(n_s);
  if (fread(&tl_r, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tl_r = orderFloatBytes(tl_r);
  if (fread(&tl_a, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tl_a = orderFloatBytes(tl_a);
  if (fread(&tl_s, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tl_s = orderFloatBytes(tl_s);
  if (fread(&tr_r, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tr_r = orderFloatBytes(tr_r);
  if (fread(&tr_a, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tr_a = orderFloatBytes(tr_a);
  if (fread(&tr_s, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tr_s = orderFloatBytes(tr_s);
  if (fread(&br_r, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  br_r = orderFloatBytes(br_r);
  if (fread(&br_a, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  br_a = orderFloatBytes(br_a);
  if (fread(&br_s, 4, 1, fp) != 1){
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  br_s = orderFloatBytes(br_s);

  nlength = sqrt(n_r * n_r + n_a * n_a + n_s * n_s);
  n_r = n_r / nlength;
  n_a = n_a / nlength;
  n_s = n_s / nlength;

  if (getenv("KILLIANY_SWAP") != NULL) {
    printf("WARNING - swapping normal direction!\n");
    n_a *= -1;
  }

  header->x_r = (tr_r - tl_r);
  header->x_a = (tr_a - tl_a);
  header->x_s = (tr_s - tl_s);
  header->y_r = (br_r - tr_r);
  header->y_a = (br_a - tr_a);
  header->y_s = (br_s - tr_s);

  /* --- normalize -- the normal vector from the file
     should have length 1, but just in case... --- */
  xlength = sqrt(header->x_r * header->x_r + header->x_a * header->x_a + header->x_s * header->x_s);
  ylength = sqrt(header->y_r * header->y_r + header->y_a * header->y_a + header->y_s * header->y_s);
  zlength = sqrt(n_r * n_r + n_a * n_a + n_s * n_s);

  header->x_r = header->x_r / xlength;
  header->x_a = header->x_a / xlength;
  header->x_s = header->x_s / xlength;
  header->y_r = header->y_r / ylength;
  header->y_a = header->y_a / ylength;
  header->y_s = header->y_s / ylength;
  header->z_r = n_r / zlength;
  header->z_a = n_a / zlength;
  header->z_s = n_s / zlength;

  header->c_r = (tl_r + br_r) / 2.0 + n_r * header->zsize * (header->depth - 1.0) / 2.0;
  header->c_a = (tl_a + br_a) / 2.0 + n_a * header->zsize * (header->depth - 1.0) / 2.0;
  header->c_s = (tl_s + br_s) / 2.0 + n_s * header->zsize * (header->depth - 1.0) / 2.0;

  header->ras_good_flag = 1;

  header->xend = header->xsize * (double)header->width / 2.0;
  header->xstart = -header->xend;
  header->yend = header->ysize * (double)header->height / 2.0;
  header->ystart = -header->yend;
  header->zend = header->zsize * (double)header->depth / 2.0;
  header->zstart = -header->zend;

  xfov = header->xend - header->xstart;
  yfov = header->yend - header->ystart;
  zfov = header->zend - header->zstart;

  header->fov = (xfov > yfov ? (xfov > zfov ? xfov : zfov) : (yfov > zfov ? yfov : zfov));

  fclose(fp);

  if (read_volume)
    mri = MRIallocSequence(header->width, header->height, header->depth, header->type, header->nframes);
  else
    mri = MRIallocHeader(header->width, header->height, header->depth, header->type, header->nframes);

  MRIcopyHeader(header, mri);
  MRIfree(&header);

  /* ----- read the volume if required ----- */
  if (read_volume) {
    int slice, frame;
    for (slice = 0, i = im_low; i <= im_high; (odd_only || even_only) ? i += 2 : i++, slice++) {
      frame = (i - im_low) % mri->nframes;
      sprintf(fname_use, fname_format, i);
      if ((fp = fopen(fname_use, "r")) == NULL) {
        MRIfree(&mri);
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error opening file %s", fname_use));
      }

      fseek(fp, 4, SEEK_SET);
      if (fread(&pixel_data_offset, 4, 1, fp) != 1) {
         ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
      }
      pixel_data_offset = orderIntBytes(pixel_data_offset);
      fseek(fp, pixel_data_offset, SEEK_SET);

      for (y = 0; y < mri->height; y++) {
        if ((int)fread(&MRISseq_vox(mri, 0, y, slice, frame), sizeof(short), mri->width, fp) != mri->width) {
          fclose(fp);
          MRIfree(&mri);
          errno = 0;
          ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error reading from file file %s", fname_use));
        }
#if (BYTE_ORDER == LITTLE_ENDIAN)
// swab(mri->slices[slice][y],
// mri->slices[slice][y], 2 * mri->width);
#if defined(SunOS)
        swab((const char *)&MRISseq_vox(mri, 0, y, slice, frame),
             (char *)&MRISseq_vox(mri, 0, y, slice, frame),
             sizeof(short) * mri->width);
#else
        {
	  std::vector<short> temp(mri->width);
	  // Note:
	  // void swab(const void *from, void *to, ssize_t n);
	  // void *memcpy(void *dest, const void *src, size_t n);
	  // Because consistency is the hobgoblin of small minds...
	  swab(&MRISseq_vox(mri, 0, y, slice, frame), temp.data(), sizeof(short) * mri->width);
	  memcpy(&MRISseq_vox(mri, 0, y, slice, frame), temp.data(), sizeof(short) * mri->width);
	}
#endif
#endif
      }

      fclose(fp);

      if (frame != (mri->nframes - 1)) slice--;

      exec_progress_callback(i - im_low, im_high - im_low + 1, 0, 1);
    }
  }

  return (mri);

} /* end genesisRead() */

static MRI *gelxRead(const char *fname, int read_volume)
{
  char fname_format[STRLEN];
  char fname_dir[STRLEN];
  char fname_base[STRLEN];
  char *c;
  MRI *mri = NULL;
  int im_init;
  int im_low, im_high;
  char fname_use[STRLEN];
  char temp_string[STRLEN];
  FILE *fp;
  int width, height;
  float tl_r, tl_a, tl_s;
  float tr_r, tr_a, tr_s;
  float br_r, br_a, br_s;
  float c_r, c_a, c_s;
  float n_r, n_a, n_s;
  float xlength, ylength, zlength;
  int i, y;
  int ecount, scount, icount;
  int good_flag;
  MRI *header = nullptr;
  float xfov, yfov, zfov;

  /* ----- check the first (passed) file ----- */
  if (!FileExists(fname)) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error opening file %s", fname));
  }

  /* ----- split the file name into name and directory ----- */
  const char *cc = strrchr(fname, '/');
  if (cc == NULL) {
    fname_dir[0] = '\0';
    strcpy(fname_base, fname);
  }
  else {
    strncpy(fname_dir, fname, (cc - fname + 1));
    fname_dir[cc - fname + 1] = '\0';
    strcpy(fname_base, cc + 1);
  }

  ecount = scount = icount = 0;
  good_flag = TRUE;
  for (c = fname_base; *c != '\0'; c++) {
    if (*c == 'e')
      ecount++;
    else if (*c == 's')
      scount++;
    else if (*c == 'i')
      icount++;
    else if (!isdigit(*c))
      good_flag = FALSE;
  }
  if (good_flag && ecount == 1 && scount == 1 && icount == 1) {
    c = strrchr(fname_base, 'i');
    im_init = atoi(c + 1);
    *c = '\0';
    int req = snprintf(fname_format, STRLEN, "%si%%d", fname_base);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  }
  else {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "genesisRead(): can't determine file name format for %s", fname));
  }

  int req = snprintf(temp_string, STRLEN, "%s", fname_format);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  req = snprintf(fname_format, STRLEN, "%s%s", fname_dir, temp_string);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  /* ----- find the low and high files ----- */
  im_low = im_init;
  do {
    im_low--;
    sprintf(fname_use, fname_format, im_low);
  } while (FileExists(fname_use));
  im_low++;

  im_high = im_init;
  do {
    im_high++;
    sprintf(fname_use, fname_format, im_high);
  } while (FileExists(fname_use));
  im_high--;

  /* ----- allocate the mri structure ----- */
  header = MRIallocHeader(1, 1, 1, MRI_SHORT, header->nframes);

  header->depth = im_high - im_low + 1;
  header->imnr0 = 1;
  header->imnr1 = header->depth;

  /* ----- get the header information from the first file ----- */
  sprintf(fname_use, fname_format, im_low);
  if ((fp = fopen(fname_use, "r")) == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error opening file %s\n", fname_use));
  }

  fseek(fp, 3236, SEEK_SET);
  if (fread(&width, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  width = orderIntBytes(width);
  if (fread(&height, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  height = orderIntBytes(height);
  header->width = width;
  header->height = height;

  strcpy(header->fname, fname);

  fseek(fp, 2184 + 28, SEEK_SET);
  if (fread(&(header->thick), 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  header->thick = orderFloatBytes(header->thick);
  header->zsize = header->thick;

  fseek(fp, 2184 + 52, SEEK_SET);
  if (fread(&(header->xsize), 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  header->xsize = orderFloatBytes(header->xsize);
  if (fread(&(header->ysize), 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  header->ysize = orderFloatBytes(header->ysize);
  header->ps = header->xsize;

  fseek(fp, 2184 + 136, SEEK_SET);
  if (fread(&c_r, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  c_r = orderFloatBytes(c_r);
  if (fread(&c_a, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  c_a = orderFloatBytes(c_a);
  if (fread(&c_s, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  c_s = orderFloatBytes(c_s);
  if (fread(&n_r, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  n_r = orderFloatBytes(n_r);
  if (fread(&n_a, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  n_a = orderFloatBytes(n_a);
  if (fread(&n_s, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  n_s = orderFloatBytes(n_s);
  if (fread(&tl_r, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tl_r = orderFloatBytes(tl_r);
  if (fread(&tl_a, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tl_a = orderFloatBytes(tl_a);
  if (fread(&tl_s, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tl_s = orderFloatBytes(tl_s);
  if (fread(&tr_r, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tr_r = orderFloatBytes(tr_r);
  if (fread(&tr_a, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tr_a = orderFloatBytes(tr_a);
  if (fread(&tr_s, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tr_s = orderFloatBytes(tr_s);
  if (fread(&br_r, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  br_r = orderFloatBytes(br_r);
  if (fread(&br_a, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  br_a = orderFloatBytes(br_a);
  if (fread(&br_s, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  br_s = orderFloatBytes(br_s);

  header->x_r = (tr_r - tl_r);
  header->x_a = (tr_a - tl_a);
  header->x_s = (tr_s - tl_s);
  header->y_r = (br_r - tr_r);
  header->y_a = (br_a - tr_a);
  header->y_s = (br_s - tr_s);

  /* --- normalize -- the normal vector from the file
     should have length 1, but just in case... --- */
  xlength = sqrt(header->x_r * header->x_r + header->x_a * header->x_a + header->x_s * header->x_s);
  ylength = sqrt(header->y_r * header->y_r + header->y_a * header->y_a + header->y_s * header->y_s);
  zlength = sqrt(n_r * n_r + n_a * n_a + n_s * n_s);

  header->x_r = header->x_r / xlength;
  header->x_a = header->x_a / xlength;
  header->x_s = header->x_s / xlength;
  header->y_r = header->y_r / ylength;
  header->y_a = header->y_a / ylength;
  header->y_s = header->y_s / ylength;
  header->z_r = n_r / zlength;
  header->z_a = n_a / zlength;
  header->z_s = n_s / zlength;

  header->c_r = (tl_r + br_r) / 2.0 + n_r * header->zsize * (header->depth - 1.0) / 2.0;
  header->c_a = (tl_a + br_a) / 2.0 + n_a * header->zsize * (header->depth - 1.0) / 2.0;
  header->c_s = (tl_s + br_s) / 2.0 + n_s * header->zsize * (header->depth - 1.0) / 2.0;

  header->ras_good_flag = 1;

  header->xend = header->xsize * (double)header->width / 2.0;
  header->xstart = -header->xend;
  header->yend = header->ysize * (double)header->height / 2.0;
  header->ystart = -header->yend;
  header->zend = header->zsize * (double)header->depth / 2.0;
  header->zstart = -header->zend;

  xfov = header->xend - header->xstart;
  yfov = header->yend - header->ystart;
  zfov = header->zend - header->zstart;

  header->fov = (xfov > yfov ? (xfov > zfov ? xfov : zfov) : (yfov > zfov ? yfov : zfov));

  fclose(fp);

  if (read_volume)
    mri = MRIalloc(header->width, header->height, header->depth, MRI_SHORT);
  else
    mri = MRIallocHeader(header->width, header->height, header->depth, MRI_SHORT, header->nframes);

  MRIcopyHeader(header, mri);
  MRIfree(&header);

  /* ----- read the volume if required ----- */
  if (read_volume) {
    for (i = im_low; i <= im_high; i++) {
      sprintf(fname_use, fname_format, i);
      if ((fp = fopen(fname_use, "r")) == NULL) {
        MRIfree(&mri);
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error opening file %s", fname_use));
      }

      fseek(fp, 8432, SEEK_SET);

      for (y = 0; y < mri->height; y++) {
        if ((int)fread(mri->slices[i - im_low][y], 2, mri->width, fp) != mri->width) {
          fclose(fp);
          MRIfree(&mri);
          errno = 0;
          ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error reading from file file %s", fname_use));
        }
#if (BYTE_ORDER == LITTLE_ENDIAN)
	{
	  std::vector<short> temp(2*mri->width);
	  // Note:
	  // void swab(const void *from, void *to, ssize_t n);
	  // void *memcpy(void *dest, const void *src, size_t n);
	  // Because consistency is the hobgoblin of small minds...
	  swab(mri->slices[i - im_low][y], temp.data(), (size_t)(2 * mri->width));
	  memcpy(mri->slices[i - im_low][y], temp.data(), (size_t)(2 * mri->width));
	}
#endif
      }

      fclose(fp);

      exec_progress_callback(i - im_low, im_high - im_low + 1, 0, 1);
    }
  }

  return (mri);

} /* end gelxRead() */

/*----------------------------------------------------------------------
  GetSPMStartFrame() - gets an environment variable called
  SPM_START_FRAME and uses its value as the number of the first frame
  for an SPM series.  If this variable does not exist, then uses  1.
  ----------------------------------------------------------------------*/
int GetSPMStartFrame(void)
{
  char *s;
  int startframe;
  s = getenv("SPM_START_FRAME");
  if (s == NULL) return (1);
  sscanf(s, "%d", &startframe);
  printf("Using env var SPM_START_FRAME = %d\n", startframe);

  return (startframe);
}

/*-------------------------------------------------------------------------
  CountAnalyzeFiles() - counts the number analyze files associated with
  the given analyze file name.  The name can take several forms:
  stem.img - a single file
  stem     - stem. If nzpad < 0, then an extension of .img is implied.
  Otherwise, it looks for a series of files named stemXXX.img
  where XXX is a zero padded integer. The width of the padding
  is controlled by nzpad. The stem is returned in ppstem (unless
  ppstem == NULL).
  Note: the files must actually exist.
  -------------------------------------------------------------------------*/
int CountAnalyzeFiles(const char *analyzefname, int nzpad, char **ppstem)
{
  int len, ncopy;
  char *stem, fmt[1000], fname[1000];
  int nfiles, keepcounting, startframe;
  FILE *fp;
  nfiles = 0;

  len = strlen(analyzefname);

  /* Determine whether the file name has a .img extension */
  if (len > 4 && strcmp(&(analyzefname[len - 4]), ".img") == 0) {
    ncopy = len - 4;
    nfiles = 1;
    if (nzpad >= 0) {
      printf(
          "ERROR: CountAnalyzeFiles: file with .img extension specified "
          "       with zero pad variable.\n");
      return (-1);
    }
  }
  else {
    ncopy = len;
    if (nzpad < 0) nfiles = 1;
  }

  /* Get the stem (ie, everything without the .img */
  stem = (char *)calloc(len + 1, sizeof(char));
  memmove(stem, analyzefname, ncopy);

  if (ppstem != NULL) *ppstem = stem;

  /* If there's only one file, check that it's there */
  if (nfiles == 1) {
    sprintf(fname, "%s.img", stem);
    if (ppstem == NULL) free(stem);
    fp = fopen(fname, "r");
    if (fp == NULL) return (0);
    fclose(fp);
    return (1);
  }

  /* If there are multiple files, count them, starting at 1 */
  sprintf(fmt, "%s%%0%dd.img", stem, nzpad);
  keepcounting = 1;
  startframe = GetSPMStartFrame();
  nfiles = 0;
  while (keepcounting) {
    sprintf(fname, fmt, nfiles + startframe);
    fp = fopen(fname, "r");
    if (fp == NULL)
      keepcounting = 0;
    else {
      fclose(fp);
      nfiles++;
    }
  }

  if (ppstem == NULL) free(stem);
  return (nfiles);
}
/*-------------------------------------------------------------------------*/
static int DumpAnalyzeHeader(FILE *fp, dsr *hdr)
{
  fprintf(fp, "Header Key\n");
  fprintf(fp, "  sizeof_hdr    %d\n", hdr->hk.sizeof_hdr);
  fprintf(fp, "  data_type     %s\n", hdr->hk.data_type);
  fprintf(fp, "  db_name       %s\n", hdr->hk.db_name);
  fprintf(fp, "  extents       %d\n", hdr->hk.extents);
  fprintf(fp, "  session_error %d\n", hdr->hk.session_error);
  fprintf(fp, "  regular       %c\n", hdr->hk.regular);
  fprintf(fp, "  hkey_un0      %c\n", hdr->hk.hkey_un0);
  fprintf(fp, "Image Dimension \n");
  fprintf(fp,
          "  dim %d %d %d %d %d %d %d %d \n",
          hdr->dime.dim[0],
          hdr->dime.dim[1],
          hdr->dime.dim[2],
          hdr->dime.dim[3],
          hdr->dime.dim[4],
          hdr->dime.dim[5],
          hdr->dime.dim[6],
          hdr->dime.dim[7]);
  fprintf(fp,
          "  pixdim %f %f %f %f %f %f %f %f \n",
          hdr->dime.pixdim[0],
          hdr->dime.pixdim[1],
          hdr->dime.pixdim[2],
          hdr->dime.pixdim[3],
          hdr->dime.pixdim[4],
          hdr->dime.pixdim[5],
          hdr->dime.pixdim[6],
          hdr->dime.pixdim[7]);
  fprintf(fp, "  vox_units     %s\n", hdr->dime.vox_units);
  fprintf(fp, "  cal_units     %s\n", hdr->dime.cal_units);
  fprintf(fp, "  datatype      %d\n", hdr->dime.datatype);
  fprintf(fp, "  bitpix        %d\n", hdr->dime.bitpix);
  fprintf(fp, "  glmax         %g\n", (float)hdr->dime.glmax);
  fprintf(fp, "  glmin         %g\n", (float)hdr->dime.glmin);
  fprintf(fp, "  orient        %d\n", hdr->hist.orient);

  return (0);
}
/*-------------------------------------------------------------------------*/
static dsr *ReadAnalyzeHeader(const char *hdrfile, int *swap, int *mritype, int *bytes_per_voxel)
{
  FILE *fp;
  dsr *hdr;

  /* Open and read the header */
  fp = fopen(hdrfile, "r");
  if (fp == NULL) {
    printf("ERROR: ReadAnalyzeHeader(): cannot open %s\n", hdrfile);
    return (NULL);
  }

  /* Read the header file */
  hdr = (dsr *)calloc(1, sizeof(dsr));
  if (fread(hdr, sizeof(dsr), 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "ReadAnalyzeHeader(): could not read file");
  }
  fclose(fp);

  *swap = 0;
  if (hdr->hk.sizeof_hdr != sizeof(dsr)) {
    *swap = 1;
    swap_analyze_header(hdr);
  }

  if (hdr->dime.datatype == DT_UNSIGNED_CHAR) {
    *mritype = MRI_UCHAR;
    *bytes_per_voxel = 1;
  }
  else if (hdr->dime.datatype == DT_SIGNED_SHORT) {
    *mritype = MRI_SHORT;
    *bytes_per_voxel = 2;
  }
  else if (hdr->dime.datatype == DT_UINT16) {
    // Can happen if this is nifti
    printf("INFO: This is an unsigned short.\n");
    //    "Unsigned short not supported, but trying to read it \n"
    //    "in as a signed short. Will be ok if no vals >= 32k.\n");
    *mritype = MRI_USHRT;
    *bytes_per_voxel = 2;
  }
  else if (hdr->dime.datatype == DT_SIGNED_INT) {
    *mritype = MRI_INT;
    *bytes_per_voxel = 4;
  }
  else if (hdr->dime.datatype == DT_FLOAT) {
    *mritype = MRI_FLOAT;
    *bytes_per_voxel = 4;
  }
  else if (hdr->dime.datatype == DT_DOUBLE) {
    *mritype = MRI_FLOAT;
    *bytes_per_voxel = 8;
  }
  else {
    free(hdr);
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,
                 "ReadAnalyzeHeader: "
                 "unsupported data type %d",
                 hdr->dime.datatype));
  }

  return (hdr);
}

/*-------------------------------------------------------------------------
  analyzeRead() - see the end of file for the old (pre-10/11/01) analyzeRead().
  The fname can take one of several forms.
  -------------------------------------------------------------------------*/
static MRI *analyzeRead(const char *fname, int read_volume)
{
  extern int N_Zero_Pad_Input;
  int nfiles, k, nread;
  char *stem;
  char imgfile[1000], hdrfile[1000], matfile[1000], fmt[1000];
  char *buf;
  FILE *fp;
  dsr *hdr;
  int swap = 0, mritype = 0, bytes_per_voxel = 0, cantreadmatfile = 0;
  int ncols, nrows, nslcs, nframes, row, slice, frame, startframe;
  MATRIX *T = NULL, *PcrsCenter, *PxyzCenter, *T1 = NULL, *Q = NULL;
  MRI *mri, *mritmp;
  float min, max;
  struct stat StatBuf;
  int nv, nreal;
  char direction[64];
  int signX, signY, signZ;
  int thiserrno, nifticode;

  fp = NULL;
  startframe = GetSPMStartFrame();

  /* Count the number of files associated with this file name,
     and get the stem. */
  nfiles = CountAnalyzeFiles(fname, N_Zero_Pad_Input, &stem);

  /* If there are no files, return NULL */
  if (nfiles < 0) return (NULL);
  if (nfiles == 0) {
    printf("ERROR: analyzeRead(): cannot find any files for %s\n", fname);
    return (NULL);
  }
  // printf("INFO: analyzeRead(): found %d files for %s\n",nfiles,fname);

  /* Create file names of header and mat files */
  if (N_Zero_Pad_Input > -1) {
    sprintf(fmt, "%s%%0%dd.%%s", stem, N_Zero_Pad_Input);
    sprintf(hdrfile, fmt, startframe, "hdr");
    sprintf(matfile, fmt, startframe, "mat");
    sprintf(imgfile, fmt, startframe, "img");
  }
  else {
    sprintf(hdrfile, "%s.hdr", stem);
    sprintf(matfile, "%s.mat", stem);
    sprintf(imgfile, "%s.img", stem);
  }

  // here, nifticode can only be 0 (ANALYZE) or 2 (Two-File NIFTI)
  nifticode = is_nifti_file(hdrfile);
  if (Gdiag > 0) printf("nifticode = %d\n", nifticode);
  if (nifticode == 2) {
    if (nfiles == 1) {
      // don't print this info, as it messes-up the output
      // when --tr, --te etc are used with mri_info
      // printf("INFO: reading as a two-file NIFTI\n");
      mri = nifti1Read(hdrfile, read_volume);
      return (mri);
    }
  }
  // It can get here if it is a multi-frame two-file nifti in
  // which each frame is stored as a separate file. In this
  // case, it reads the header with the analyze reader, but
  // then reads in the vox2ras matrix with the nifti reader.
  // It does not look like it makes a difference.
  hdr = ReadAnalyzeHeader(hdrfile, &swap, &mritype, &bytes_per_voxel);
  if (mritype == MRI_FLOAT && getenv("ATROPHY_SIMULATOR") != NULL) {
    printf("adjusting analyze type to uchar to correct simulator bug\n");
    mritype = MRI_UCHAR;
    bytes_per_voxel = 1;
  }
  if (hdr == NULL) return (NULL);
  if (Gdiag_no > 0) DumpAnalyzeHeader(stdout, hdr);

  /* Get the number of frames as either the fourth dimension or
     the number of files. */
  if (nfiles == 1) {
    nframes = hdr->dime.dim[4];
    if (nframes == 0) nframes = 1;
  }
  else
    nframes = nfiles;

  ncols = hdr->dime.dim[1];
  nrows = hdr->dime.dim[2];
  nslcs = hdr->dime.dim[3];
  nv = ncols * nrows * nslcs * nframes;

  if (ncols == 1 || nrows == 1) {
    lstat(imgfile, &StatBuf);
    if (StatBuf.st_size != nv * bytes_per_voxel) {
      nreal = StatBuf.st_size / bytes_per_voxel;
      if (ncols == 1)
        nrows = nreal;
      else
        ncols = nreal;
      printf("INFO: %s really has %d rows/cols\n", imgfile, nreal);
    }
  }

  /* Alloc the Header and/or Volume */
  if (read_volume)
    mri = MRIallocSequence(ncols, nrows, nslcs, mritype, nframes);
  else {
    mri = MRIallocHeader(ncols, nrows, nslcs, mritype, nframes);
    mri->nframes = nframes;
  }

  /* Load Variables into header */
  mri->xsize = fabs(hdr->dime.pixdim[1]); /* col res */
  mri->ysize = fabs(hdr->dime.pixdim[2]); /* row res */
  mri->zsize = fabs(hdr->dime.pixdim[3]); /* slice res */
  mri->tr = 1000 * hdr->dime.pixdim[4];   /* time  res */

  signX = (hdr->dime.pixdim[1] > 0) ? 1 : -1;
  signY = (hdr->dime.pixdim[2] > 0) ? 1 : -1;
  signZ = (hdr->dime.pixdim[3] > 0) ? 1 : -1;

  // Handel the vox2ras matrix
  if (nifticode != 2) {  // Not a nifti file
    /* Read the matfile, if there */
    if (FileExists(matfile)) {
      T1 = MatlabRead(matfile);  // orientation info
      if (T1 == NULL) {
        printf("WARNING: analyzeRead(): matfile %s exists but could not read ... \n", matfile);
        printf("  may not be matlab4 mat file ... proceeding without it.\n");
        fflush(stdout);
        cantreadmatfile = 1;
      }
      else {
        /* Convert from 1-based to 0-based */
        Q = MtxCRS1toCRS0(Q);
        T = MatrixMultiply(T1, Q, T);
        // printf("------- Analyze Input Matrix (zero-based) --------\n");
        // MatrixPrint(stdout,T);
        // printf("-------------------------------------\n");
        mri->ras_good_flag = 1;
        MatrixFree(&Q);
        MatrixFree(&T1);
      }
    }
    if (!FileExists(matfile) || cantreadmatfile) {
      /* when not found, it is a fun exercise.                      */
      /* see http://wideman-one.com/gw/brain/analyze/formatdoc.htm  */
      /* for amgibuities.                                           */
      /* I will follow his advise                                   */
      /* hist.orient  Mayo name   Voxel[Index0, Index1, Index2] */
      /*                          Index0  Index1  Index2              */
      /* 0 transverse unflipped   R-L     P-A     I-S     LAS */
      /* 3 transverse flipped     R-L     A-P     I-S     LPS */
      /* 1 coronal unflipped      R-L     I-S     P-A     LSA */
      /* 4 coronal flipped        R-L     S-I     P-A     LIA */
      /* 2 sagittal unflipped     P-A     I-S     R-L     ASL */
      /* 5 sagittal flipped       P-A     S-I     R-L     AIL */
      //   P->A I->S L->R

      /* FLIRT distributes analyze format image which has a marked LR */
      /* in fls/etc/standard/avg152T1_LR-marked.img.  The convention    */
      /* is tested with this image for hdr->hist.orient== 0.              */
      /* note that flirt image has negative pixdim[1], but their software */
      /* ignores the sign anyway.                                         */
      if (hdr->hist.orient == 0) /* x = - r, y = a, z = s */
      {
        strcpy(direction, "transverse unflipped (default)");
        T = MatrixAlloc(4, 4, MATRIX_REAL);
        T->rptr[1][1] = -mri->xsize;
        T->rptr[2][2] = mri->ysize;
        T->rptr[3][3] = mri->zsize;
        T->rptr[1][4] = mri->xsize * (mri->width / 2.0);
        T->rptr[2][4] = -mri->ysize * (mri->height / 2.0);
        T->rptr[3][4] = -mri->zsize * (mri->depth / 2.0);
        T->rptr[4][4] = 1.;
      }
      else if (hdr->hist.orient == 1) /* x = -r, y = s, z = a */
      {
        strcpy(direction, "coronal unflipped");
        T = MatrixAlloc(4, 4, MATRIX_REAL);
        T->rptr[1][1] = -mri->xsize;
        T->rptr[2][3] = mri->zsize;
        T->rptr[3][2] = mri->ysize;
        T->rptr[1][4] = mri->xsize * (mri->width / 2.0);
        T->rptr[2][4] = -mri->zsize * (mri->depth / 2.0);
        T->rptr[3][4] = -mri->ysize * (mri->height / 2.0);
        T->rptr[4][4] = 1.;
      }
      else if (hdr->hist.orient == 2) /* x = a, y = s, z = -r */
      {
        strcpy(direction, "sagittal unflipped");
        T = MatrixAlloc(4, 4, MATRIX_REAL);
        T->rptr[1][3] = -mri->zsize;
        T->rptr[2][1] = mri->xsize;
        T->rptr[3][2] = mri->ysize;
        T->rptr[1][4] = mri->zsize * (mri->depth / 2.0);
        T->rptr[2][4] = -mri->xsize * (mri->width / 2.0);
        T->rptr[3][4] = -mri->ysize * (mri->height / 2.0);
        T->rptr[4][4] = 1.;
      }
      else if (hdr->hist.orient == 3) /* x = -r, y = -a, z = s */
      {
        strcpy(direction, "transverse flipped");
        T = MatrixAlloc(4, 4, MATRIX_REAL);
        T->rptr[1][1] = -mri->xsize;
        T->rptr[2][2] = -mri->ysize;
        T->rptr[3][3] = mri->zsize;
        T->rptr[1][4] = mri->xsize * (mri->width / 2.0);
        T->rptr[2][4] = mri->ysize * (mri->height / 2.0);
        T->rptr[3][4] = -mri->zsize * (mri->depth / 2.0);
        T->rptr[4][4] = 1.;
      }
      else if (hdr->hist.orient == 4) /* x = -r, y = -s, z = a */
      {
        strcpy(direction, "coronal flipped");
        T = MatrixAlloc(4, 4, MATRIX_REAL);
        T->rptr[1][1] = -mri->xsize;
        T->rptr[2][3] = mri->zsize;
        T->rptr[3][2] = -mri->ysize;
        T->rptr[1][4] = mri->xsize * (mri->width / 2.0);
        T->rptr[2][4] = -mri->zsize * (mri->depth / 2.0);
        T->rptr[3][4] = mri->ysize * (mri->height / 2.0);
        T->rptr[4][4] = 1.;
      }
      else if (hdr->hist.orient == 5) /* x = a, y = -s, z = -r */
      {
        strcpy(direction, "sagittal flipped");
        T = MatrixAlloc(4, 4, MATRIX_REAL);
        T->rptr[1][3] = -mri->zsize;
        T->rptr[2][1] = mri->xsize;
        T->rptr[3][2] = -mri->ysize;
        T->rptr[1][4] = mri->zsize * (mri->depth / 2.0);
        T->rptr[2][4] = -mri->xsize * (mri->width / 2.0);
        T->rptr[3][4] = mri->ysize * (mri->height / 2.0);
        T->rptr[4][4] = 1.;
      }
      if (hdr->hist.orient == (char)-1) {
        /* Unknown, so assume: x = -r, y = -a, z = s */
        // This is incompatible with mghRead() when rasgood=0.
        // mghRead() uses coronal dircos, not transverse.
        // I'm not changing it because don't know what will happen.
        strcpy(direction, "transverse flipped");
        T = MatrixAlloc(4, 4, MATRIX_REAL);
        T->rptr[1][1] = -mri->xsize;
        T->rptr[2][2] = -mri->ysize;
        T->rptr[3][3] = mri->zsize;
        T->rptr[1][4] = mri->xsize * (mri->width / 2.0);
        T->rptr[2][4] = mri->ysize * (mri->height / 2.0);
        T->rptr[3][4] = -mri->zsize * (mri->depth / 2.0);
        T->rptr[4][4] = 1.;
        mri->ras_good_flag = 0;
        fprintf(stderr,
                "WARNING: could not find %s file for direction cosine info.\n"
                "WARNING: Analyze 7.5 hdr->hist.orient value = -1, not valid.\n"
                "WARNING: assuming %s\n",
                matfile,
                direction);
      }
      else {
        // It's probably not a good idea to set this to 1, but setting it
        // to 0 created all kinds of problems with mghRead() which will
        // force dir cos to be Coronal if rasgood<0.
        mri->ras_good_flag = 1;
        fprintf(stderr,
                "-----------------------------------------------------------------\n"
                "INFO: could not find %s file for direction cosine info.\n"
                "INFO: use Analyze 7.5 hdr->hist.orient value: %d, %s.\n"
                "INFO: if not valid, please provide the information in %s file\n"
                "-----------------------------------------------------------------\n",
                matfile,
                hdr->hist.orient,
                direction,
                matfile);
      }
    }
  }
  else {
    // Just read in this one file as nifti to get vox2ras matrix
    // What a hack.
    mritmp = MRIreadHeader(imgfile, NIFTI1_FILE);
    T = MRIxfmCRS2XYZ(mritmp, 0);
    MRIfree(&mritmp);
  }

  /* ---- Assign the Geometric Paramaters -----*/
  mri->x_r = T->rptr[1][1] / mri->xsize;
  mri->x_a = T->rptr[2][1] / mri->xsize;
  mri->x_s = T->rptr[3][1] / mri->xsize;

  /* Row Direction Cosines */
  mri->y_r = T->rptr[1][2] / mri->ysize;
  mri->y_a = T->rptr[2][2] / mri->ysize;
  mri->y_s = T->rptr[3][2] / mri->ysize;

  /* Slice Direction Cosines */
  mri->z_r = T->rptr[1][3] / mri->zsize;
  mri->z_a = T->rptr[2][3] / mri->zsize;
  mri->z_s = T->rptr[3][3] / mri->zsize;

  /* Center of the FOV in voxels */
  PcrsCenter = MatrixAlloc(4, 1, MATRIX_REAL);
  PcrsCenter->rptr[1][1] = mri->width / 2.0;
  PcrsCenter->rptr[2][1] = mri->height / 2.0;
  PcrsCenter->rptr[3][1] = mri->depth / 2.0;
  PcrsCenter->rptr[4][1] = 1.0;

  /* Center of the FOV in XYZ */
  PxyzCenter = MatrixMultiply(T, PcrsCenter, NULL);
  mri->c_r = PxyzCenter->rptr[1][1];
  mri->c_a = PxyzCenter->rptr[2][1];
  mri->c_s = PxyzCenter->rptr[3][1];

  MatrixFree(&PcrsCenter);
  MatrixFree(&PxyzCenter);
  MatrixFree(&T);

  if (!read_volume) return (mri);

  /* Alloc the maximum amount of memory that a row could need */
  buf = (char *)malloc(mri->width * 8);
  if (NULL == buf) {
    printf("ERROR: analyzeRead(): malloc failure\n");
    MRIfree(&mri);
    return (NULL);
  }

  /* Open the one file, if there is one file */
  if (nfiles == 1) {
    fp = fopen(imgfile, "r");
    if (fp == NULL) {
      printf("ERROR: analyzeRead(): could not open %s\n", imgfile);
      MRIfree(&mri);
      return (NULL);
    }
    fseek(fp, (int)(hdr->dime.vox_offset), SEEK_SET);
  }

  /*--------------------- Frame Loop --------------------------*/
  for (frame = 0; frame < nframes; frame++) {
    /* Open the frame file if there is more than one file */
    if (N_Zero_Pad_Input > -1) {
      sprintf(imgfile, fmt, frame + startframe, "img");
      fp = fopen(imgfile, "r");
      if (fp == NULL) {
        printf("ERROR: analyzeRead(): could not open %s\n", imgfile);
        MRIfree(&mri);
        return (NULL);
      }
      fseek(fp, (int)(hdr->dime.vox_offset), SEEK_SET);
    }

    /* ----------- Slice Loop ----------------*/
    for (slice = 0; slice < mri->depth; slice++) {
      k = slice + mri->depth * frame;

      /* --------- Row Loop ------------------*/
      for (row = 0; row < mri->height; row++) {
        nread = fread(buf, bytes_per_voxel, mri->width, fp);
        thiserrno = errno;
        if (nread != mri->width) {
          if (feof(fp)) printf("ERROR: premature end of file\n");
          printf("ERROR: %s (%d)\n", strerror(thiserrno), thiserrno);
          errno = 0;
          printf("frame = %d, slice = %d, k=%d, row = %d\n", frame, slice, k, row);
          printf("nread = %d, nexpected = %d\n", nread, mri->width);
          MRIdump(mri, stdout);
          fflush(stdout);
          MRIfree(&mri);
          free(buf);
          fclose(fp);
          ErrorReturn(NULL,
                      (ERROR_BADFILE,
                       "analyzeRead2(): error reading "
                       "from file %s\n",
                       imgfile));
        }

        if (swap) {
          if (bytes_per_voxel == 2) byteswapbufshort((void *)buf, bytes_per_voxel * mri->width);
          if (bytes_per_voxel == 4) byteswapbuffloat((void *)buf, bytes_per_voxel * mri->width);
          if (bytes_per_voxel == 8) byteswapbufdouble((void *)buf, bytes_per_voxel * mri->width);
          // nflip(buf, bytes_per_voxel, mri->width); /* byte swap */
        }

        // check if this is 64bit double data, which we have to convert
        // to float, since the mri struct only supports float
        if ((bytes_per_voxel == 8) && (hdr->dime.datatype == DT_DOUBLE) && (mritype == MRI_FLOAT)) {
          // convert double to float
          int width;
          double *dp = (double *)buf;
          for (width = 0; width < mri->width; width++) {
            MRIFvox(mri, width, row, slice) = (float)*dp++;
          }
        }
        else {
          /*copy*/
          memmove(mri->slices[k][row], buf, bytes_per_voxel * mri->width);
        }

      } /* End Row Loop */

      exec_progress_callback(slice, mri->depth, frame, nframes);
    } /* End Slice Loop */

    /* Close file if each frame is in a separate file */
    if (N_Zero_Pad_Input > -1) fclose(fp);

  } /* End Frame Loop */

  if (N_Zero_Pad_Input < 0) fclose(fp);

  printf("  analyzeRead() roi_scale %13.9f\n", hdr->dime.roi_scale);
  if (getenv("FS_ANALYZE_NO_RESCALE") != NULL) printf("FS_ANALYZE_NO_RESCALE set, so not rescaling\n");
  if (fabs(hdr->dime.roi_scale - 1) > FLT_EPSILON && fabs(hdr->dime.roi_scale) > FLT_EPSILON &&
      getenv("FS_ANALYZE_NO_RESCALE") == NULL) {
    // Rescale if it is neither 1 nor 0
    if (mri->type != MRI_FLOAT) {
      MRI *mritmp;
      printf("  analyzeRead() changing type to float\n");
      mritmp = MRISeqchangeType(mri, MRI_FLOAT, 0, 0, 0);
      MRIfree(&mri);
      mri = mritmp;
    }
    printf("  analyzeRead() scaling by %13.9f\n", hdr->dime.roi_scale);
    mri = MRImultiplyConst(mri, hdr->dime.roi_scale, mri);
  }
  free(buf);
  free(hdr);

  if (Gdiag_no > 0) {
    MRIlimits(mri, &min, &max);
    printf("INFO: analyzeRead(): min = %g, max = %g\n", min, max);
  }
  return (mri);
}

/*---------------------------------------------------------------
  analyzeWrite() - writes out a volume in SPM analyze format. If the
  file name has a .img extension, then the first frame is stored into
  the given file name. If it does not include the extension, then the
  volume is saved as a series of frame files using fname as a base
  followed by the zero-padded frame number (see analyzeWriteSeries).  A
  series can be saved from mri_convert by specifying the basename and
  adding "--out_type spm". See also analyzeWrite4D(). DNG
  ---------------------------------------------------------------*/
static int analyzeWrite(MRI *mri, const char *fname)
{
  int len;
  int error_value;

  /* Check for .img extension */
  len = strlen(fname);
  if (len > 4) {
    if (fname[len - 4] == '.' && fname[len - 3] == 'i' && fname[len - 2] == 'm' && fname[len - 1] == 'g') {
      /* There is a trailing .img - save frame 0 into fname */
      error_value = analyzeWriteFrame(mri, fname, 0);
      return (error_value);
    }
  }
  /* There is no trailing .img - save as a series */
  error_value = analyzeWriteSeries(mri, fname);
  return (error_value);
}
/*---------------------------------------------------------------
  analyzeWriteFrame() - this used to be analyzeWrite() but was modified
  by DNG to be able to save a particular frame so that it could be
  used to write out an entire series.
  ---------------------------------------------------------------*/
static void printDirCos(MRI *mri)
{
  fprintf(stderr, "Direction cosines for %s are:\n", mri->fname);
  fprintf(stderr, "  x_r = %8.4f, y_r = %8.4f, z_r = %8.4f, c_r = %10.4f\n", mri->x_r, mri->y_r, mri->z_r, mri->c_r);
  fprintf(stderr, "  x_a = %8.4f, y_a = %8.4f, z_a = %8.4f, c_a = %10.4f\n", mri->x_a, mri->y_a, mri->z_a, mri->c_a);
  fprintf(stderr, "  x_s = %8.4f, y_s = %8.4f, z_s = %8.4f, c_s = %10.4f\n", mri->x_s, mri->y_s, mri->z_s, mri->c_s);
}

static int analyzeWriteFrame(MRI *mri, const char *fname, int frame)
{
  dsr hdr;
  float max, min, det;
  MATRIX *T, *invT;
  char hdr_fname[STRLEN];
  char mat_fname[STRLEN];
  const char *c;
  FILE *fp;
  int error_value;
  int i, j, k;
  int bytes_per_voxel;
  short i1, i2, i3;
  int shortmax;
  const char *orientname[7] = {"transverse unflipped",
                         "coronal unflipped",
                         "sagittal unflipped",
                         "transverse flipped",
                         "coronal flipped",
                         "sagittal flipped",
                         "unknown"};

  if (frame >= mri->nframes) {
    fprintf(stderr,
            "ERROR: analyzeWriteFrame(): frame number (%d) exceeds "
            "number of frames (%d)\n",
            frame,
            mri->nframes);
    return (1);
  }

  shortmax = (int)(pow(2.0, 15.0));
  if (mri->width > shortmax) {
    printf("ANALYZE FORMAT ERROR: ncols %d in volume exceeds %d\n", mri->width, shortmax);
    exit(1);
  }
  if (mri->height > shortmax) {
    printf("ANALYZE FORMAT ERROR: nrows %d in volume exceeds %d\n", mri->height, shortmax);
    exit(1);
  }
  if (mri->depth > shortmax) {
    printf("ANALYZE FORMAT ERROR: nslices %d in volume exceeds %d\n", mri->depth, shortmax);
    exit(1);
  }
  if (mri->nframes > shortmax) {
    printf("ANALYZE FORMAT ERROR:  nframes %d in volume exceeds %d\n", mri->nframes, shortmax);
    exit(1);
  }

  c = strrchr(fname, '.');
  if (c == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "analyzeWriteFrame(): "
                 "bad file name %s",
                 fname));
  }
  if (strcmp(c, ".img") != 0) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "analyzeWriteFrame(): "
                 "bad file name %s",
                 fname));
  }

  strcpy(hdr_fname, fname);
  sprintf(hdr_fname + (c - fname), ".hdr");

  strcpy(mat_fname, fname);
  sprintf(mat_fname + (c - fname), ".mat");

  memset(&hdr, 0x00, sizeof(hdr));

  hdr.hk.sizeof_hdr = sizeof(hdr);

  hdr.dime.vox_offset = 0.0;

  if (mri->type == MRI_UCHAR) {
    hdr.dime.datatype = DT_UNSIGNED_CHAR;
    bytes_per_voxel = 1;
  }
  else if (mri->type == MRI_SHORT) {
    hdr.dime.datatype = DT_SIGNED_SHORT;
    bytes_per_voxel = 2;
  }
  else if (mri->type == MRI_USHRT) {
    hdr.dime.datatype = DT_UINT16;
    bytes_per_voxel = 2;
  }
  /* --- assuming long and int are identical --- */
  else if (mri->type == MRI_INT || mri->type == MRI_LONG) {
    hdr.dime.datatype = DT_SIGNED_INT;
    bytes_per_voxel = 4;
  }
  else if (mri->type == MRI_FLOAT) {
    hdr.dime.datatype = DT_FLOAT;
    bytes_per_voxel = 4;
  }
  else {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "analyzeWriteFrame(): bad data type %d", mri->type));
  }

  /* Added by DNG 10/2/01 */
  hdr.dime.dim[0] = 4;                    /* number of dimensions */
  hdr.dime.dim[4] = 1;                    /* time */
  hdr.dime.pixdim[4] = mri->tr / 1000.0;  // convert to sec
  hdr.dime.bitpix = 8 * bytes_per_voxel;
  memmove(hdr.dime.vox_units, "mm\0", 3);
  /*----------------------------*/

  hdr.dime.dim[1] = mri->width;  /* ncols */
  hdr.dime.dim[2] = mri->height; /* nrows */
  hdr.dime.dim[3] = mri->depth;  /* nslices */

  hdr.dime.pixdim[1] = mri->xsize; /* col res */
  hdr.dime.pixdim[2] = mri->ysize; /* row res */
  hdr.dime.pixdim[3] = mri->zsize; /* slice res */

  MRIlimits(mri, &min, &max);
  hdr.dime.glmin = (int)min;
  hdr.dime.glmax = (int)max;

  /* Construct the matrix to convert CRS to XYZ, assuming
     that CRS is 1-based */
  T = MRIxfmCRS2XYZ(mri, 1);

  det = MatrixDeterminant(T);
  if (det == 0) {
    printf(
        "WARNING: cannot determine volume orientation, "
        "assuming identity. It's ok if the output is a surface.\n");
    T = MatrixIdentity(4, T);
    if (mri->xsize > 0) T->rptr[1][1] = mri->xsize;
    if (mri->ysize > 0) T->rptr[2][2] = mri->ysize;
    if (mri->zsize > 0) T->rptr[3][3] = mri->zsize;
  }
  if (frame == 0) {
    printf("Analyze Output Matrix\n");
    MatrixPrint(stdout, T);
    printf("--------------------\n");
  }

  /* ----- write the matrix to the .mat file ----- */
  error_value = MatlabWrite(T, mat_fname, "M");
  if (error_value != NO_ERROR) return (error_value);

  /* Compute inverse ot T (converts from XYZ to CRS) */
  invT = MatrixInverse(T, NULL);
  /* If the inverse cannot be computed, set to the identity */
  if (invT == NULL) invT = MatrixIdentity(4, NULL);

  /* Load the CRS into the originator field as 3 shorts for SPM */
  /* These come from the last column of invT */
  i1 = (short)(rint(*MATRIX_RELT(invT, 1, 4)));
  memmove((&hdr.hist.originator[0]), &i1, sizeof(short));
  i2 = (short)(rint(*MATRIX_RELT(invT, 2, 4)));
  memmove((&hdr.hist.originator[0] + sizeof(short)), &i2, sizeof(short));
  i3 = (short)(rint(*MATRIX_RELT(invT, 3, 4)));
  memmove((&hdr.hist.originator[0] + 2 * sizeof(short)), &i3, sizeof(short));

  MatrixFree(&T);
  MatrixFree(&invT);

  /* Set the hist.orient field -- this is not always correct */
  /* see http://wideman-one.com/gw/brain/analyze/formatdoc.htm  */
  if (fabs(mri->z_s) > fabs(mri->z_r) && fabs(mri->z_s) > fabs(mri->z_a)) {
    // Transverse: Superior/Inferior > both Right and Anterior
    if (mri->x_r < 0 && mri->y_a > 0 && mri->z_s > 0)
      hdr.hist.orient = 0;  // transverse unflipped  LAS
    else if (mri->x_r < 0 && mri->y_a < 0 && mri->z_s > 0)
      hdr.hist.orient = 3;  // transverse flipped    LPS
    else {
      fprintf(stderr,
              "No such orientation specified in Analyze7.5. "
              "Set orient to 0.\n");
      fprintf(stderr, "Not a problem as long as you use the .mat file.\n");
      printDirCos(mri);
      hdr.hist.orient = 0;
    }
  }
  if (fabs(mri->z_a) > fabs(mri->z_r) && fabs(mri->z_a) > fabs(mri->z_s)) {
    // Cor: Anterior/Post > both Right and Superior
    if (mri->x_r < 0 && mri->y_s > 0 && mri->z_a > 0)
      hdr.hist.orient = 1;  // cor unflipped   LSA
    else if (mri->x_r < 0 && mri->y_s < 0 && mri->z_a > 0)
      hdr.hist.orient = 4;  // cor flipped     LIA
    else {
      fprintf(stderr,
              "No such orientation specified in Analyze7.5. "
              "Set orient to 0\n");
      printDirCos(mri);
      hdr.hist.orient = 0;
    }
  }
  if (fabs(mri->z_r) > fabs(mri->z_a) && fabs(mri->z_r) > fabs(mri->z_s)) {
    // Sag: Righ/Left > both Anterior and Superior
    if (mri->x_a > 0 && mri->y_s > 0 && mri->z_r < 0)
      hdr.hist.orient = 2;  // sag unflipped   ASL
    else if (mri->x_a > 0 && mri->y_s < 0 && mri->z_r < 0)
      hdr.hist.orient = 5;  // sag flipped     AIL
    else {
      fprintf(stderr, "No such orientation specified in Analyze7.5. Set orient to 0\n");
      printDirCos(mri);
      hdr.hist.orient = 0;
    }
  }
  printf("INFO: set hdr.hist.orient to '%s'\n", orientname[(int)hdr.hist.orient]);

  /* ----- open the header file ----- */
  if ((fp = fopen(hdr_fname, "w")) == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "analyzeWriteFrame(): error opening file %s for writing", hdr_fname));
  }
  /* ----- write the header ----- */
  /* should write element-by-element */
  if (fwrite(&hdr, sizeof(hdr), 1, fp) != 1) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "analyzeWriteFrame(): error writing to file %s", hdr_fname));
  }
  fclose(fp);
  if (Gdiag_no > 0) DumpAnalyzeHeader(stdout, &hdr);

  /* ----- open the data file ----- */
  if ((fp = fopen(fname, "w")) == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "analyzeWriteFrame(): error opening file %s for writing", fname));
  }

  /* ----- write the data ----- */
  for (i = 0; i < mri->depth; i++) {
    /* this is the change to make it save a given frame */
    k = i + mri->depth * frame;
    for (j = 0; j < mri->height; j++) {
      if ((int)fwrite(mri->slices[k][j], bytes_per_voxel, mri->width, fp) != mri->width) {
        errno = 0;
        ErrorReturn(ERROR_BADFILE,
                    (ERROR_BADFILE,
                     "analyzeWriteFrame(): "
                     " error writing to file %s",
                     fname));
      } /* end if */
    }   /* end row loop */
    if (mri->nframes > 1) exec_progress_callback(i, mri->depth, 0, 1);
  } /* end col loop */

  fclose(fp);

  return (0);

} /* end analyzeWriteFrame() */

/*-------------------------------------------------------------*/
static int analyzeWriteSeries(MRI *mri, const char *fname)
{
  extern int N_Zero_Pad_Output;
  int frame;
  int err;
  char framename[STRLEN];
  char spmnamefmt[STRLEN];

  /* NOTE: This function assumes that fname does not have a .img extension. */

  /* N_Zero_Pad_Output is a global variable used to control the name of      */
  /* files in successive frames within the series. It can be set from     */
  /* mri_convert using "--spmnzeropad N" where N is the width of the pad. */
  if (N_Zero_Pad_Output < 0) N_Zero_Pad_Output = 3;

  /* Create the format string used to create the filename for each frame.   */
  /* The frame file name format will be fname%0Nd.img where N is the amount */
  /* of zero padding. The (padded) frame numbers will go from 1 to nframes. */
  sprintf(spmnamefmt, "%%s%%0%dd.img", N_Zero_Pad_Output);

  /* loop over slices */
  for (frame = 0; frame < mri->nframes; frame++) {
    sprintf(framename, spmnamefmt, fname, frame + 1);
    // printf("%3d %s\n",frame,framename);
    err = analyzeWriteFrame(mri, framename, frame);
    if (err) return (err);

    exec_progress_callback(frame, mri->nframes, 0, 1);
  }

  return (0);
}

/*---------------------------------------------------------------
  analyzeWrite4D() - saves data in analyze 4D format.
  ---------------------------------------------------------------*/
static int analyzeWrite4D(MRI *mri, const char *fname)
{
  dsr hdr;
  float max, min, det;
  MATRIX *T, *invT;
  char hdr_fname[STRLEN];
  char mat_fname[STRLEN];
  const char *c;
  FILE *fp;
  int error_value;
  int i, j, k, frame;
  int bytes_per_voxel;
  short i1, i2, i3;
  int shortmax;

  shortmax = (int)(pow(2.0, 15.0));
  if (mri->width > shortmax) {
    printf("ANALYZE FORMAT ERROR: ncols %d in volume exceeds %d\n", mri->width, shortmax);
    exit(1);
  }
  if (mri->height > shortmax) {
    printf("ANALYZE FORMAT ERROR: nrows %d in volume exceeds %d\n", mri->height, shortmax);
    exit(1);
  }
  if (mri->depth > shortmax) {
    printf("ANALYZE FORMAT ERROR: nslices %d in volume exceeds %d\n", mri->depth, shortmax);
    exit(1);
  }
  if (mri->nframes > shortmax) {
    printf("ANALYZE FORMAT ERROR:  nframes %d in volume exceeds %d\n", mri->nframes, shortmax);
    exit(1);
  }

  c = strrchr(fname, '.');
  if (c == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "analyzeWrite4D(): "
                 "bad file name %s",
                 fname));
  }
  if (strcmp(c, ".img") != 0) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "analyzeWrite4D(): "
                 "bad file name %s",
                 fname));
  }

  /* create the file name for the header */
  strcpy(hdr_fname, fname);
  sprintf(hdr_fname + (c - fname), ".hdr");

  /* create the file name for the mat file */
  strcpy(mat_fname, fname);
  sprintf(mat_fname + (c - fname), ".mat");

  memset(&hdr, 0x00, sizeof(hdr));
  hdr.hk.sizeof_hdr = sizeof(hdr);

  hdr.dime.vox_offset = 0.0;

  if (mri->type == MRI_UCHAR) {
    hdr.dime.datatype = DT_UNSIGNED_CHAR;
    bytes_per_voxel = 1;
  }
  else if (mri->type == MRI_SHORT) {
    hdr.dime.datatype = DT_SIGNED_SHORT;
    bytes_per_voxel = 2;
  }
  else if (mri->type == MRI_USHRT) {
    hdr.dime.datatype = DT_UINT16;
    bytes_per_voxel = 2;
  }
  /* --- assuming long and int are identical --- */
  else if (mri->type == MRI_INT || mri->type == MRI_LONG) {
    hdr.dime.datatype = DT_SIGNED_INT;
    bytes_per_voxel = 4;
  }
  else if (mri->type == MRI_FLOAT) {
    hdr.dime.datatype = DT_FLOAT;
    bytes_per_voxel = 4;
  }
  else {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "analyzeWrite4D(): bad data type %d", mri->type));
  }

  hdr.dime.bitpix = 8 * bytes_per_voxel;
  memmove(hdr.dime.vox_units, "mm\0", 3);

  hdr.dime.dim[1] = mri->width;  /* ncols */
  hdr.dime.dim[2] = mri->height; /* nrows */
  hdr.dime.dim[3] = mri->depth;  /* nslices */
  hdr.dime.dim[4] = mri->nframes;
  hdr.dime.dim[0] = 4; /* flirt expects to be always 4 */

  hdr.dime.pixdim[1] = mri->xsize;       /* col res */
  hdr.dime.pixdim[2] = mri->ysize;       /* row res */
  hdr.dime.pixdim[3] = mri->zsize;       /* slice res */
  hdr.dime.pixdim[4] = mri->tr / 1000.0; /* time res in sec*/

  MRIlimits(mri, &min, &max);
  hdr.dime.glmin = (int)min;
  hdr.dime.glmax = (int)max;

  /* Construct the matrix to convert CRS to XYZ, assuming
     that CRS is 1-based */
  T = MRIxfmCRS2XYZ(mri, 1);

  det = MatrixDeterminant(T);
  if (det == 0) {
    printf(
        "WARNING: cannot determine volume orientation, "
        "assuming identity. It's ok if the output is a surface.\n");
    T = MatrixIdentity(4, T);
    if (mri->xsize > 0) T->rptr[1][1] = mri->xsize;
    if (mri->ysize > 0) T->rptr[2][2] = mri->ysize;
    if (mri->zsize > 0) T->rptr[3][3] = mri->zsize;
  }
  printf("Analyze Output Matrix\n");
  MatrixPrint(stdout, T);
  printf("--------------------\n");

  /* Set the hist.orient field -- this is not always correct */
  /* see http://wideman-one.com/gw/brain/analyze/formatdoc.htm  */
  if (fabs(mri->z_s) > fabs(mri->z_r) && fabs(mri->z_s) > fabs(mri->z_a)) {
    // Transverse: Superior/Inferior > both Right and Anterior
    if (mri->x_r > 0)
      hdr.hist.orient = 0;  // transverse unflipped
    else
      hdr.hist.orient = 3;  // transverse flipped
  }
  if (fabs(mri->z_a) > fabs(mri->z_r) && fabs(mri->z_a) > fabs(mri->z_s)) {
    // Cor: Anterior/Post > both Right and Superior
    if (mri->x_r > 0)
      hdr.hist.orient = 1;  // cor unflipped
    else
      hdr.hist.orient = 4;  // cor flipped
  }
  if (fabs(mri->z_r) > fabs(mri->z_a) && fabs(mri->z_r) > fabs(mri->z_s)) {
    // Sag: Righ/Left > both Anterior and Superior
    if (mri->x_r > 0)
      hdr.hist.orient = 2;  // sag unflipped
    else
      hdr.hist.orient = 5;  // sag flipped
  }
  hdr.hist.orient = -1;
  printf("INFO: set hdr.hist.orient to %d\n", hdr.hist.orient);

  /* ----- write T to the  .mat file ----- */
  error_value = MatlabWrite(T, mat_fname, "M");
  if (error_value != NO_ERROR) return (error_value);

  /* This matrix converts from XYZ to CRS */
  invT = MatrixInverse(T, NULL);
  /* If the inverse cannot be computed, set to the identity */
  if (invT == NULL) invT = MatrixIdentity(4, NULL);

  /* Load the CRS into the originator field as 3 shorts for SPM */
  /* These come from the last column of invT */
  i1 = (short)(rint(*MATRIX_RELT(invT, 1, 4)));
  memmove((&hdr.hist.originator[0]), &i1, sizeof(short));
  i2 = (short)(rint(*MATRIX_RELT(invT, 2, 4)));
  memmove((&hdr.hist.originator[0] + sizeof(short)), &i2, sizeof(short));
  i3 = (short)(rint(*MATRIX_RELT(invT, 3, 4)));
  memmove((&hdr.hist.originator[0] + 2 * sizeof(short)), &i3, sizeof(short));

  MatrixFree(&T);
  MatrixFree(&invT);

  /* ----- write the header ----- */
  if ((fp = fopen(hdr_fname, "w")) == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "analyzeWrite4D(): error opening file %s for writing", hdr_fname));
  }
  if (fwrite(&hdr, sizeof(hdr), 1, fp) != 1) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "analyzeWrite4D(): error writing to file %s", hdr_fname));
  }
  fclose(fp);
  if (Gdiag_no > 0) DumpAnalyzeHeader(stdout, &hdr);

  /* ----- write the data ----- */
  if ((fp = fopen(fname, "w")) == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "analyzeWrite4D(): error opening file %s for writing", fname));
  }

  for (frame = 0; frame < mri->nframes; frame++) {
    for (i = 0; i < mri->depth; i++) {
      k = i + mri->depth * frame;
      for (j = 0; j < mri->height; j++) {
        if ((int)fwrite(mri->slices[k][j], bytes_per_voxel, mri->width, fp) != mri->width) {
          errno = 0;
          ErrorReturn(ERROR_BADFILE,
                      (ERROR_BADFILE,
                       "analyzeWrite4D(): "
                       "error writing to file %s",
                       fname));
        }
      }
      exec_progress_callback(i, mri->depth, frame, mri->nframes);
    }
  }

  fclose(fp);

  return (0);

} /* end analyzeWrite4D() */

/*-------------------------------------------------------------*/
static void swap_analyze_header(dsr *hdr)
{
  int i;
  char c;

  hdr->hk.sizeof_hdr = swapInt(hdr->hk.sizeof_hdr);
  hdr->hk.extents = swapShort(hdr->hk.extents);
  hdr->hk.session_error = swapShort(hdr->hk.session_error);

  for (i = 0; i < 5; i++) hdr->dime.dim[i] = swapShort(hdr->dime.dim[i]);
  hdr->dime.unused1 = swapShort(hdr->dime.unused1);
  hdr->dime.datatype = swapShort(hdr->dime.datatype);
  hdr->dime.bitpix = swapShort(hdr->dime.bitpix);
  hdr->dime.dim_un0 = swapShort(hdr->dime.dim_un0);
  hdr->dime.vox_offset = swapFloat(hdr->dime.vox_offset);
  hdr->dime.roi_scale = swapFloat(hdr->dime.roi_scale);
  hdr->dime.funused1 = swapFloat(hdr->dime.funused1);
  hdr->dime.funused2 = swapFloat(hdr->dime.funused2);
  hdr->dime.cal_max = swapFloat(hdr->dime.cal_max);
  hdr->dime.cal_min = swapFloat(hdr->dime.cal_min);
  hdr->dime.compressed = swapInt(hdr->dime.compressed);
  hdr->dime.verified = swapInt(hdr->dime.verified);
  hdr->dime.glmin = swapInt(hdr->dime.glmin);
  hdr->dime.glmax = swapInt(hdr->dime.glmax);
  for (i = 0; i < 8; i++) hdr->dime.pixdim[i] = swapFloat(hdr->dime.pixdim[i]);

  hdr->hist.views = swapInt(hdr->hist.views);
  hdr->hist.vols_added = swapInt(hdr->hist.vols_added);
  hdr->hist.start_field = swapInt(hdr->hist.start_field);
  hdr->hist.field_skip = swapInt(hdr->hist.field_skip);
  hdr->hist.omax = swapInt(hdr->hist.omax);
  hdr->hist.omin = swapInt(hdr->hist.omin);
  hdr->hist.smax = swapInt(hdr->hist.smax);
  hdr->hist.smin = swapInt(hdr->hist.smin);

  /* spm uses the originator char[10] as shorts */
  for (i = 0; i < 5; i++) {
    c = hdr->hist.originator[2 * i + 1];
    hdr->hist.originator[2 * i + 1] = hdr->hist.originator[2 * i];
    hdr->hist.originator[2 * i] = c;
  }

} /* end swap_analyze_header */
/*------------------------------------------------------*/



static MRI *gdfRead(const char *fname, int read_volume)
{
  MRI *mri;
  FILE *fp;
  char line[STRLEN];
  char *c;
  char file_path[STRLEN];
  float ipr[2];
  float st;
  char units_string[STRLEN], orientation_string[STRLEN], data_type_string[STRLEN];
  int size[2];
  int path_d, ipr_d, st_d, u_d, dt_d, o_d, s_d, x_ras_d, y_ras_d, z_ras_d, c_ras_d;
  int data_type;
  int orientation = MRI_UNDEFINED;
  char *or_ptr;
  char os_orig[STRLEN];
  float units_factor;
  char file_path_1[STRLEN], file_path_2[STRLEN];
  int i, j, k;
  short *sbuf = NULL;
  float *fbuf = NULL;
  unsigned char *ucbuf = NULL;
  int n_files;
  char fname_use[STRLEN];
  int pad_zeros_flag;
  int file_offset = 0;
  float x_r, x_a, x_s;
  float y_r, y_a, y_s;
  float z_r, z_a, z_s;
  float c_r, c_a, c_s;
  int have_min_crop = FALSE;
  int have_max_crop = FALSE;
  float min_crop[3], max_crop[3];

  if ((fp = fopen(fname, "r")) == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "gdfRead(): error opening file %s", fname));
  }

  /* --- defined flags --- */
  path_d = ipr_d = st_d = u_d = dt_d = o_d = s_d = x_ras_d = y_ras_d = z_ras_d = c_ras_d = FALSE;

  while (fgets(line, STRLEN, fp) != NULL) {
    /* --- strip the newline --- */
    if ((c = strrchr(line, '\n')) != NULL) *c = '\0';

    if (strncmp(line, "IMAGE_FILE_PATH", 15) == 0) {
      sscanf(line, "%*s %s", file_path);
      path_d = TRUE;
    }
    else if (strncmp(line, "IP_RES", 6) == 0) {
      sscanf(line, "%*s %f %f", &ipr[0], &ipr[1]);
      ipr_d = TRUE;
    }
    else if (strncmp(line, "SL_THICK", 8) == 0) {
      sscanf(line, "%*s %f", &st);
      st_d = TRUE;
    }
    else if (strncmp(line, "UNITS", 5) == 0) {
      sscanf(line, "%*s %s", units_string);
      u_d = TRUE;
    }
    else if (strncmp(line, "SIZE", 4) == 0) {
      sscanf(line, "%*s %d %d", &size[0], &size[1]);
      s_d = TRUE;
    }
    else if (strncmp(line, "FS_X_RAS", 8) == 0) {
      sscanf(line, "%*s %f %f %f", &x_r, &x_a, &x_s);
      x_ras_d = TRUE;
    }
    else if (strncmp(line, "FS_Y_RAS", 8) == 0) {
      sscanf(line, "%*s %f %f %f", &y_r, &y_a, &y_s);
      y_ras_d = TRUE;
    }
    else if (strncmp(line, "FS_Z_RAS", 8) == 0) {
      sscanf(line, "%*s %f %f %f", &z_r, &z_a, &z_s);
      z_ras_d = TRUE;
    }
    else if (strncmp(line, "FS_C_RAS", 8) == 0) {
      sscanf(line, "%*s %f %f %f", &c_r, &c_a, &c_s);
      c_ras_d = TRUE;
    }
    else if (strncmp(line, "DATA_TYPE", 9) == 0) {
      strcpy(data_type_string, &line[10]);
      dt_d = TRUE;
    }
    else if (strncmp(line, "ORIENTATION", 11) == 0) {
      sscanf(line, "%*s %s", orientation_string);
      strcpy(os_orig, orientation_string);
      o_d = TRUE;
    }
    else if (strncmp(line, "FILE_OFFSET", 11) == 0) {
      sscanf(line, "%*s %d", &file_offset);
    }
    else if (strncmp(line, "MIN_CROP", 8) == 0) {
      sscanf(line, "%*s %f %f %f", &min_crop[0], &min_crop[1], &min_crop[2]);
      have_min_crop = TRUE;
    }
    else if (strncmp(line, "MAX_CROP", 8) == 0) {
      sscanf(line, "%*s %f %f %f", &max_crop[0], &max_crop[1], &max_crop[2]);
      have_max_crop = TRUE;
    }
    else {
    }
  }

  fclose(fp);

  if (!(path_d && ipr_d && st_d && s_d)) {
    errno = 0;
    ErrorPrintf(ERROR_BADPARM, "gdfRead(): missing field(s) from %s:", fname);
    if (!path_d) ErrorPrintf(ERROR_BADPARM, "  IMAGE_FILE_PATH");
    if (!ipr_d) ErrorPrintf(ERROR_BADPARM, "  IP_RES");
    if (!st_d) ErrorPrintf(ERROR_BADPARM, "  SL_THICK");
    if (!s_d) ErrorPrintf(ERROR_BADPARM, "  SIZE");
    return (NULL);
  }

  if (!(o_d)) {
    if (x_ras_d && y_ras_d && z_ras_d) {
      printf(
          "missing field ORIENTATION in file %s, "
          "but you've got {xyz}_{ras}, so never mind\n",
          fname);
    }
    else {
      printf("missing field ORIENTATION in file %s; assuming 'coronal'\n", fname);
      sprintf(orientation_string, "coronal");
    }
  }

  if (!(dt_d)) {
    printf("missing field DATA_TYPE in file %s; assuming 'short'\n", fname);
    sprintf(data_type_string, "short");
  }

  if (!(u_d)) {
    printf("missing field UNITS in file %s; assuming 'mm'\n", fname);
    sprintf(units_string, "mm");
  }

  StrLower(data_type_string);
  StrLower(orientation_string);
  StrLower(units_string);

  /* --- data type --- */
  if (strcmp(data_type_string, "short") == 0)
    data_type = MRI_SHORT;
  else if (strcmp(data_type_string, "float") == 0)
    data_type = MRI_FLOAT;
  else if (strcmp(data_type_string, "unsigned char") == 0)
    data_type = MRI_UCHAR;
  else {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "gdfRead(): unknown data type '%s'", data_type_string));
  }

  /* --- orientation --- */
  or_ptr = strrchr(orientation_string, ' ');
  or_ptr = (or_ptr == NULL ? orientation_string : or_ptr + 1);
  if (strncmp(or_ptr, "cor", 3) == 0)
    orientation = MRI_CORONAL;
  else if (strncmp(or_ptr, "sag", 3) == 0)
    orientation = MRI_SAGITTAL;
  else if (strncmp(or_ptr, "ax", 2) == 0 || strncmp(or_ptr, "hor", 3) == 0)
    orientation = MRI_HORIZONTAL;
  else if (!(x_ras_d && y_ras_d && z_ras_d)) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "gdfRead(): can't determine orientation from string '%s'", os_orig));
  }

  if (strcmp(units_string, "mm") == 0)
    units_factor = 1.0;
  else if (strcmp(units_string, "cm") == 0)
    units_factor = 1.0;
  else if (strcmp(units_string, "m") == 0)
    units_factor = 100.0;
  else {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "gdfRead(): unknown units '%s'", units_string));
  }

  ipr[0] /= units_factor;
  ipr[1] /= units_factor;
  st /= units_factor;

  if (gdf_crop_flag && !(have_min_crop && have_max_crop)) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "gdfRead(): cropping desired but missing MIN_CROP or MAX_CROP\n"));
  }

  strcpy(file_path_1, file_path);
  c = strrchr(file_path_1, '*');
  if (c == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "gdfRead(): file path %s does not contain '*'\n", file_path));
  }

  *c = '\0';
  c++;
  strcpy(file_path_2, c);


  pad_zeros_flag = FALSE;

  n_files = 0;
  do {
    n_files++;
    int req = snprintf(fname_use, STRLEN, "%s%d%s", file_path_1, n_files, file_path_2);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  } while (FileExists(fname_use));

  /* ----- try padding the zeros if no files are found ----- */
  if (n_files == 1) {
    pad_zeros_flag = TRUE;

    n_files = 0;
    do {
      n_files++;
      int req = snprintf(fname_use, STRLEN, "%s%03d%s", file_path_1, n_files, file_path_2);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } while (FileExists(fname_use));

    /* ----- still a problem? ----- */
    if (n_files == 1) {
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_BADFILE, "gdfRead(): can't find file %s%d%s or %s\n", file_path_1, 1, file_path_2, fname_use));
    }
  }

  n_files--;

  if (read_volume)
    mri = MRIallocSequence(size[0], size[1], n_files, data_type, 1);
  else
    mri = MRIallocHeader(size[0], size[1], n_files, data_type, 1);

  mri->xsize = ipr[0];
  mri->ysize = ipr[1];
  mri->zsize = st;

  mri->xend = mri->width * mri->xsize / 2.0;
  mri->xstart = -mri->xend;
  mri->yend = mri->height * mri->ysize / 2.0;
  mri->ystart = -mri->yend;
  mri->zend = mri->depth * mri->zsize / 2.0;
  mri->zstart = -mri->zend;

  strcpy(mri->fname, fname);

  /* --- set volume orientation --- */
  if (x_ras_d && y_ras_d && z_ras_d) {
    mri->x_r = x_r;
    mri->x_a = x_a;
    mri->x_s = x_s;
    mri->y_r = y_r;
    mri->y_a = y_a;
    mri->y_s = y_s;
    mri->z_r = z_r;
    mri->z_a = z_a;
    mri->z_s = z_s;
    mri->ras_good_flag = TRUE;
  }
  else {
    /*
       direction cosine is not set.  we pick a particular kind
       of direction cosine.  If the volume is different you have
       to modify (how?)
    */
    if (setDirectionCosine(mri, orientation) != NO_ERROR) {
      MRIfree(&mri);
      return NULL;
    }
    printf("warning: gdf volume may be incorrectly oriented\n");
  }

  /* --- set volume center --- */
  if (c_ras_d) {
    mri->c_r = c_r;
    mri->c_a = c_a;
    mri->c_s = c_s;
  }
  else {
    mri->c_r = mri->c_a = mri->c_s = 0.0;
    printf("warning: gdf volume may be incorrectly centered\n");
  }

  if (!read_volume) return (mri);

  if (mri->type == MRI_UCHAR) {
    ucbuf = (unsigned char *)malloc(mri->width);
    if (ucbuf == NULL) {
      MRIfree(&mri);
      errno = 0;
      ErrorReturn(NULL, (ERROR_NOMEMORY, "gdfRead(): error allocating %d bytes for read buffer", mri->width));
    }
  }
  else if (mri->type == MRI_SHORT) {
    sbuf = (short *)malloc(mri->width * sizeof(short));
    if (sbuf == NULL) {
      MRIfree(&mri);
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_NOMEMORY, "gdfRead(): error allocating %d bytes for read buffer", mri->width * sizeof(short)));
    }
  }
  else if (mri->type == MRI_FLOAT) {
    fbuf = (float *)malloc(mri->width * sizeof(float));
    if (fbuf == NULL) {
      MRIfree(&mri);
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_NOMEMORY, "gdfRead(): error allocating %d bytes for read buffer", mri->width * sizeof(float)));
    }
  }
  else {
    MRIfree(&mri);
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "gdfRead(): internal error; data type %d "
                 "accepted but not supported in read",
                 mri->type));
  }

  for (i = 1; i <= n_files; i++) {
    if (pad_zeros_flag) {
      int req = snprintf(fname_use, STRLEN, "%s%03d%s", file_path_1, i, file_path_2);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname_use, STRLEN, "%s%d%s", file_path_1, i, file_path_2);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }

    fp = fopen(fname_use, "r");
    if (fp == NULL) {
      if (mri->type == MRI_UCHAR) free(ucbuf);
      if (mri->type == MRI_SHORT) free(sbuf);
      if (mri->type == MRI_FLOAT) free(fbuf);
      MRIfree(&mri);
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADFILE, "gdfRead(): error opening file %s", fname_use));
    }

    fseek(fp, file_offset, SEEK_SET);

    if (mri->type == MRI_UCHAR) {
      for (j = 0; j < mri->height; j++) {
        if ((int)fread(ucbuf, 1, mri->width, fp) != mri->width) {
          free(ucbuf);
          MRIfree(&mri);
          errno = 0;
          ErrorReturn(NULL, (ERROR_BADFILE, "gdfRead(): error reading from file %s", fname_use));
        }
        for (k = 0; k < mri->width; k++) MRIvox(mri, k, j, i - 1) = ucbuf[k];
      }
    }

    if (mri->type == MRI_SHORT) {
      for (j = 0; j < mri->height; j++) {
        if ((int)fread(sbuf, sizeof(short), mri->width, fp) != mri->width) {
          free(sbuf);
          MRIfree(&mri);
          errno = 0;
          ErrorReturn(NULL, (ERROR_BADFILE, "gdfRead(): error reading from file %s", fname_use));
        }
#if (BYTE_ORDER == LITTLE_ENDIAN)
        for (k = 0; k < mri->width; k++) MRISvox(mri, k, j, i - 1) = orderShortBytes(sbuf[k]);
#else
        for (k = 0; k < mri->width; k++) MRISvox(mri, k, j, i - 1) = sbuf[k];
#endif
      }
    }

    if (mri->type == MRI_FLOAT) {
      for (j = 0; j < mri->height; j++) {
        if ((int)fread(fbuf, sizeof(float), mri->width, fp) != mri->width) {
          free(fbuf);
          MRIfree(&mri);
          errno = 0;
          ErrorReturn(NULL, (ERROR_BADFILE, "gdfRead(): error reading from file %s", fname_use));
        }
#if (BYTE_ORDER == LITTLE_ENDIAN)
        for (k = 0; k < mri->width; k++) MRIFvox(mri, k, j, i - 1) = orderFloatBytes(fbuf[k]);
#else
        for (k = 0; k < mri->width; k++) MRIFvox(mri, k, j, i - 1) = fbuf[k];
#endif
      }
    }

    fclose(fp);

    exec_progress_callback(i - 1, n_files, 0, 1);
  }

  if (mri->type == MRI_UCHAR) free(ucbuf);
  if (mri->type == MRI_SHORT) free(sbuf);
  if (mri->type == MRI_FLOAT) free(fbuf);

  /* ----- crop if desired (and possible!) ----- */

  if (gdf_crop_flag) {
    /* we've checked for MIN_CROP and MAX_CROP above, */
    /* before reading the data... */

    printf("NOTICE: cropping GDF data\n");

    for (k = 0; k < mri->depth; k++) {
      for (j = 0; j < mri->height; j++) {
        for (i = 0; i < mri->width; i++) {
          if (i < (int)min_crop[0] || i > (int)max_crop[0] || j < (int)min_crop[2] || j > (int)max_crop[2] ||
              k < (int)min_crop[1] - 1 || k > (int)max_crop[1] - 1) {
            if (mri->type == MRI_UCHAR)
              MRIvox(mri, i, j, k) = 0;
            else if (mri->type == MRI_SHORT)
              MRISvox(mri, i, j, k) = 0;
            else if (mri->type == MRI_FLOAT)
              MRIFvox(mri, i, j, k) = 0.0;
            else {
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "gdfRead(): internal error (data type while cropping"));
            }
          }
        }
      }
    }
  }
  else {
    printf("NOTICE: not cropping GDF data\n");
  }

  return (mri);

} /* end gdfRead() */

static int gdfWrite(MRI *mri, const char *fname)
{
  FILE *fp;
  int i, j;
  std::string im_fname;
  unsigned char *buf;
  int buf_size = 0;

  if (strlen(mri->gdf_image_stem) == 0) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "GDF write attempted without GDF image stem"));
  }

  if (mri->nframes != 1) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "GDF write attempted with %d frames "
                 "(supported only for 1 frame)",
                 mri->nframes));
  }

  /* ----- type checks first ----- */

  if (mri->type != MRI_UCHAR && mri->type != MRI_FLOAT && mri->type != MRI_SHORT) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "gdfWrite(): bad data type (%d) for GDF write "
                 "(only uchar, float, short supported)",
                 mri->type));
  }

  /* ----- then write the image files ----- */

  printf("writing GDF image files...\n");

  if (mri->type == MRI_UCHAR) buf_size = mri->width * sizeof(unsigned char);
  if (mri->type == MRI_FLOAT) buf_size = mri->width * sizeof(float);
  if (mri->type == MRI_SHORT) buf_size = mri->width * sizeof(short);

  buf = (unsigned char *)malloc(buf_size);
  if (buf == NULL) {
    errno = 0;
    ErrorReturn(ERROR_NO_MEMORY, (ERROR_NO_MEMORY, "gdfWrite(): no memory for voxel write buffer"));
  }

  for (i = 0; i < mri->depth; i++) {
    im_fname = std::string(mri->gdf_image_stem) + '_' + std::to_string(i+1) + ".img";
    fp = fopen(im_fname.c_str(), "w");
    if (fp == NULL) {
      free(buf);
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "gdfWrite(): error opening file %s", im_fname.c_str()));
    }

    for (j = 0; j < mri->height; j++) {
      memmove(buf, mri->slices[i][j], buf_size);
#if (BYTE_ORDER == LITTLE_ENDIAN)
      if (mri->type == MRI_FLOAT) byteswapbuffloat(buf, buf_size);
      if (mri->type == MRI_SHORT) byteswapbufshort(buf, buf_size);
#endif
      fwrite(buf, 1, buf_size, fp);
    }

    fclose(fp);

    exec_progress_callback(i, mri->depth, 0, 1);
  }

  free(buf);

  /* ----- and finally the info file ----- */

  printf("writing GDF info file...\n");

  fp = fopen(fname, "w");
  if (fp == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "gdfWrite(): error opening file %s", fname));
  }

  fprintf(fp, "GDF FILE VERSION3\n");
  fprintf(fp, "START MAIN HEADER\n");
  fprintf(fp, "IMAGE_FILE_PATH %s_*.img\n", mri->gdf_image_stem);
  fprintf(fp, "SIZE %d %d\n", mri->width, mri->height);
  fprintf(fp, "IP_RES %g %g\n", mri->xsize, mri->ysize);
  fprintf(fp, "SL_THICK %g\n", mri->zsize);
  if (mri->ras_good_flag) {
    fprintf(fp, "FS_X_RAS %f %f %f\n", mri->x_r, mri->x_a, mri->x_s);
    fprintf(fp, "FS_Y_RAS %f %f %f\n", mri->y_r, mri->y_a, mri->y_s);
    fprintf(fp, "FS_Z_RAS %f %f %f\n", mri->z_r, mri->z_a, mri->z_s);
    fprintf(fp, "FS_C_RAS %f %f %f\n", mri->c_r, mri->c_a, mri->c_s);
  }
  fprintf(fp, "FILE_OFFSET 0\n");
  if (mri->type == MRI_UCHAR) fprintf(fp, "DATA_TYPE unsigned char\n");
  if (mri->type == MRI_FLOAT) fprintf(fp, "DATA_TYPE float\n");
  if (mri->type == MRI_SHORT) fprintf(fp, "DATA_TYPE short\n");
  fprintf(fp, "END MAIN HEADER\n");

  fclose(fp);

  return (NO_ERROR);
}


static int register_unknown_label(const char *label)
{
  int i;

  if (n_unknown_labels == MAX_UNKNOWN_LABELS) return (NO_ERROR);

  for (i = 0; i < n_unknown_labels; i++) {
    if (strcmp(unknown_labels[i], label) == 0) return (NO_ERROR);
  }

  strcpy(unknown_labels[n_unknown_labels], label);
  n_unknown_labels++;

  return (NO_ERROR);

} /* end register_unknown_label() */

static int clear_unknown_labels(void)
{
  n_unknown_labels = 0;

  return (NO_ERROR);

} /* end clear_unknown_labels() */

static int print_unknown_labels(const char *prefix)
{
  int i;

  for (i = 0; i < n_unknown_labels; i++) printf("%s%s\n", prefix, unknown_labels[i]);

  return (NO_ERROR);

} /* end print_unknown_labels() */

static int read_otl_file(
    FILE *fp, MRI *mri, int slice, COLOR_TABLE *ctab, int fill_flag, int translate_label_flag, int zero_outlines_flag)
{
  int n_outlines = -1;
  int n_rows, n_cols;
  char label[STRLEN], label_to_compare[STRLEN];
  int seed_x, seed_y;
  char line[STRLEN];
  int main_header_flag;
  int i, j;
  int gdf_header_flag;
  char type[STRLEN], global_type[STRLEN];
  char *c;
  short *points;
  int n_read;
  short label_value;
  float scale_x, scale_y;
  int source_x, source_y;
  char alt_compare[STRLEN];
  int empty_label_flag;
  int ascii_short_flag;
  int row;
  char *translate_start;
  int internal_structures_flag = FALSE;
  CMAoutlineField *of;
  int num_entries;
  char entry_name[STRLEN];

  if (!fgets(line, STRLEN, fp) && ferror(fp)) {
    ErrorPrintf(ERROR_BADFILE, "read_otl_file(): could not read file");
  }
  if (strncmp(line, "GDF FILE VERSION", 15) != 0) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "otl slice %s does not appear to be a GDF file", slice));
  }

  main_header_flag = FALSE;
  while (!main_header_flag) {
    if (feof(fp)) {
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF () in otl file %d", slice));
    }
    if (!fgets(line, STRLEN, fp) && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "read_otl_file(): could not read file");
    }
    if (strncmp(line, "START MAIN HEADER", 17) == 0) main_header_flag = TRUE;
  }

  n_cols = -1;
  type[0] = '\0';
  global_type[0] = '\0';

  while (main_header_flag) {
    if (feof(fp)) {
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF (in main header) in otl file %d", slice));
    }
    if (!fgets(line, STRLEN, fp) && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "read_otl_file(): could not read file");
    }
    if (strncmp(line, "END MAIN HEADER", 15) == 0) main_header_flag = FALSE;
    if (strncmp(line, "ONUM ", 5) == 0) sscanf(line, "%*s %d", &n_outlines);
    if (strncmp(line, "COL_NUM", 7) == 0) sscanf(line, "%*s %d", &n_cols);
    if (strncmp(line, "TYPE", 4) == 0) {
      strcpy(global_type, &(line[5]));
      c = strrchr(global_type, '\n');
      if (c != NULL) *c = '\0';
    }
  }

  if (n_outlines == -1) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "bad or undefined ONUM in otl file %d", slice));
  }

  of = CMAoutlineFieldAlloc(2 * mri->width, 2 * mri->height);
  if (of == NULL) return (ERROR_NOMEMORY);

  for (i = 0; i < n_outlines; i++) {
    if (feof(fp)) {
      CMAfreeOutlineField(&of);
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF (ready for next outline) in otl file %d", slice));
    }

    gdf_header_flag = FALSE;

    while (!gdf_header_flag) {
      if (feof(fp)) {
        CMAfreeOutlineField(&of);
        errno = 0;
        ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF (searching for gdf header) in otl file %d", slice));
      }
      if (!fgets(line, STRLEN, fp) && ferror(fp)) {
        ErrorPrintf(ERROR_BADFILE, "read_otl_file(): could not read file");
      }
      if (strncmp(line, "START GDF HEADER", 16) == 0) gdf_header_flag = TRUE;
    }

    n_rows = -1;
    seed_x = seed_y = -1;
    label[0] = '\0';
    type[0] = '\0';

    empty_label_flag = 0;

    while (gdf_header_flag) {
      if (feof(fp)) {
        CMAfreeOutlineField(&of);
        errno = 0;
        ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF (in gdf header) in otl file %d", slice));
      }
      if (!fgets(line, STRLEN, fp) && ferror(fp)) {
        ErrorPrintf(ERROR_BADFILE, "read_otl_file(): could not read file");
      }
      if (strncmp(line, "END GDF HEADER", 14) == 0) gdf_header_flag = FALSE;
      // getting rows, cols, type
      if (strncmp(line, "ROW_NUM", 7) == 0) sscanf(line, "%*s %d", &n_rows);
      if (strncmp(line, "COL_NUM", 7) == 0) sscanf(line, "%*s %d", &n_cols);
      if (strncmp(line, "TYPE", 4) == 0) {
        strcpy(type, &(line[5]));
        c = strrchr(type, '\n');
        if (c != NULL) *c = '\0';
      }
      if (strncmp(line, "SEED", 4) == 0) sscanf(line, "%*s %d %d", &seed_x, &seed_y);
      if (strncmp(line, "LABEL", 5) == 0) {
        strcpy(label, &(line[6]));
        c = strrchr(label, '\n');
        if (c != NULL) *c = '\0';

        /* exterior -> cortex, if desired */
        if (translate_label_flag) {
          translate_start = strstr(label, "Exterior");
          if (translate_start != NULL)
            sprintf(translate_start, "Cortex");
          else {
            translate_start = strstr(label, "exterior");
            if (translate_start != NULL) sprintf(translate_start, "cortex");
          }
        }

        /* warning if there's an "exterior" or
           "cortex" after any other label */
        if (strstr(label, "Exterior") == 0 || strstr(label, "exterior") == 0 || strstr(label, "Cortex") == 0 ||
            strstr(label, "cortex") == 0) {
          if (internal_structures_flag)
            printf(
                "WARNING: label \"%s\" following "
                "non-exterior labels in slice %d\n",
                label,
                slice);
          internal_structures_flag = FALSE;
        }
        else
          internal_structures_flag = TRUE;
      }
    }

    if (n_rows < 0) {
      CMAfreeOutlineField(&of);
      errno = 0;
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "bad or undefined ROW_NUM in otl file %d", slice));
    }

    if (n_cols != 2) {
      CMAfreeOutlineField(&of);
      errno = 0;
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "bad or undefined COL_NUM in otl file %d", slice));
    }

    if (label[0] == '\0') {
      empty_label_flag = 1;
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "empty LABEL in otl file %d (outline %d)", slice, i);
    }

    if (seed_x < 0 || seed_x >= 512 || seed_y < 0 || seed_y >= 512) {
      CMAfreeOutlineField(&of);
      errno = 0;
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "bad or undefined SEED in otl file %d", slice));
    }

    if (type[0] == '\0') strcpy(type, global_type);

    if (strcmp(type, "ascii short") == 0)
      ascii_short_flag = TRUE;
    else if (strcmp(type, "short") == 0)
      ascii_short_flag = FALSE;
    else if (type[0] == '\0') {
      CMAfreeOutlineField(&of);
      errno = 0;
      ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "undefined TYPE in otl file %d", slice));
    }
    else {
      CMAfreeOutlineField(&of);
      errno = 0;
      ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "unsupported TYPE \"%s\" in otl file %d", type, slice));
    }

    do {
      if (!fgets(line, STRLEN, fp) && ferror(fp)) {
        ErrorPrintf(ERROR_BADFILE, "read_otl_file(): could not read file");
      }
      if (feof(fp)) {
        CMAfreeOutlineField(&of);
        errno = 0;
        ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF (searching for points) in otl file %d", slice));
      }
    } while (strncmp(line, "START POINTS", 12) != 0);

    points = (short *)malloc(2 * n_rows * sizeof(short));
    if (points == NULL) {
      CMAfreeOutlineField(&of);
      errno = 0;
      ErrorReturn(ERROR_NOMEMORY, (ERROR_NOMEMORY, "error allocating memory for points in otl file %d", slice));
    }

    if (ascii_short_flag) {
      for (row = 0; row < n_rows; row++) {
        if (!fgets(line, STRLEN, fp) && ferror(fp)) {
          ErrorPrintf(ERROR_BADFILE, "read_otl_file(): could not read file");
        }
        if (feof(fp)) {
          free(points);
          CMAfreeOutlineField(&of);
          errno = 0;
          ErrorReturn(ERROR_BADFILE,
                      (ERROR_BADFILE,
                       "premature end of file while reading "
                       "points from otl file %d",
                       slice));
        }
        sscanf(line, "%hd %hd", &(points[2 * row]), &(points[2 * row + 1]));
      }
    }
    else {
      n_read = fread(points, 2, n_rows * 2, fp);
      if (n_read != n_rows * 2) {
        free(points);
        CMAfreeOutlineField(&of);
        errno = 0;
        ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "error reading points from otl file %d", slice));
      }
#if (BYTE_ORDER == LITTLE_ENDIAN)
#if defined(SunOS)
      swab((const char *)points, (char *)points, 2 * n_rows * sizeof(short));
#else
      {
	std::vector<short> tmp(2*n_rows);
	// Note:
	// void swab(const void *from, void *to, ssize_t n);
	// void *memcpy(void *dest, const void *src, size_t n);
	// Because consistency is the hobgoblin of small minds...
	swab(points, tmp.data(), 2 * n_rows * sizeof(short));
	memcpy(points, tmp.data(), 2 * n_rows * sizeof(short));
      }
#endif
#endif
    }

    if (!fgets(line, STRLEN, fp) && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "read_otl_file(): could not read file");
    }
    if (strncmp(line, "END POINTS", 10) != 0) {
      free(points);
      CMAfreeOutlineField(&of);
      errno = 0;
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE,
                   "error after points (\"END POINTS\" expected "
                   "but not found) in otl file %d",
                   slice));
    }

    if (!empty_label_flag) {
      strcpy(label_to_compare, label);
      for (j = 0; label_to_compare[j] != '\0'; j++) {
        if (label_to_compare[j] == '\n') label_to_compare[j] = '\0';
        if (label_to_compare[j] == ' ') label_to_compare[j] = '-';
      }

      /* --- strip 'Left' and 'Right'; --- */
      strcpy(alt_compare, label_to_compare);
      StrLower(alt_compare);
      if (strncmp(alt_compare, "left-", 5) == 0)
        strcpy(alt_compare, &(label_to_compare[5]));
      else if (strncmp(alt_compare, "right-", 6) == 0)
        strcpy(alt_compare, &(label_to_compare[6]));

      /* --- strip leading and trailing spaces (now dashes) --- */

      /* leading */
      for (j = 0; label_to_compare[j] == '-'; j++)
        ;
      if (label_to_compare[j] != '\0') {
        for (c = label_to_compare; *(c + j) != '\0'; c++) *c = *(c + j);
        *c = *(c + j);
      }

      /* trailing */
      for (j = strlen(label_to_compare) - 1; j >= 0 && label_to_compare[j] == '-'; j--)
        ;
      if (j < 0) /* all dashes! */
      {
        /* for now, let this fall through to an unknown label */
      }
      else /* j is the index of the last non-dash character */
      {
        label_to_compare[j + 1] = '\0';
      }

      label_value = -1;
      CTABgetNumberOfTotalEntries(ctab, &num_entries);
      for (j = 0; j < num_entries; j++) {
        if (NO_ERROR == CTABcopyName(ctab, j, entry_name, sizeof(entry_name))) {
          //	  printf("%s compared to  (%s, %s)\n", entry_name, label_to_compare, alt_compare);
          if (strcmp(entry_name, label_to_compare) == 0 || strcmp(entry_name, alt_compare) == 0) {
            //	      printf("FOUND \n");
            label_value = j;
          }
        }
      }

      if (label_value == -1) {
        register_unknown_label(label);
      }
      else {
        for (j = 0; j < n_rows; j++) cma_field[points[2 * j]][points[2 * j + 1]] = label_value;

        CMAclaimPoints(of, label_value, points, n_rows, seed_x, seed_y);
      }
    }

    free(points);
  }

  CMAassignLabels(of);

  if (zero_outlines_flag) CMAzeroOutlines(of);

  scale_x = 512 / mri->width;
  scale_y = 512 / mri->height;

  for (i = 0; i < mri->width; i++)
    for (j = 0; j < mri->height; j++) {
      source_x = (int)floor(scale_x * (float)i);
      source_y = (int)floor(scale_y * (float)j);
      if (source_x < 0) source_x = 0;
      if (source_x > 511) source_x = 511;
      if (source_y < 0) source_y = 0;
      if (source_y > 511) source_y = 511;
      MRISvox(mri, i, j, slice - 1) = of->fill_field[source_y][source_x];
    }

  CMAfreeOutlineField(&of);

  return (NO_ERROR);

} /* end read_otl_file() */

int list_labels_in_otl_file(FILE *fp)
{
  char line[STRLEN];
  int main_header_flag;
  int n_outlines = -1;
  int n_rows, n_cols;
  char type[STRLEN], global_type[STRLEN];
  char *c;
  int i, gdf_header_flag;
  int seed_x, seed_y;
  char label[STRLEN];
  int ascii_short_flag;
  short *points;
  int row;
  int n_read;
  int empty_label_flag;

  if (!fgets(line, STRLEN, fp) && ferror(fp)) {
    ErrorPrintf(ERROR_BADFILE, "list_labels_in_otl_file(): could not read file");
  }
  if (strncmp(line, "GDF FILE VERSION", 15) != 0) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "otl slice does not appear to be a GDF file"));
  }

  main_header_flag = FALSE;
  while (!main_header_flag) {
    if (feof(fp)) {
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF () in otl file"));
    }
    if (!fgets(line, STRLEN, fp) && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "list_labels_in_otl_file(): could not read file");
    }
    if (strncmp(line, "START MAIN HEADER", 17) == 0) main_header_flag = TRUE;
  }

  n_cols = -1;
  type[0] = '\0';
  global_type[0] = '\0';

  while (main_header_flag) {
    if (feof(fp)) {
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF (in main header) in otl file"));
    }
    if (!fgets(line, STRLEN, fp) && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "list_labels_in_otl_file(): could not read file");
    }
    if (strncmp(line, "END MAIN HEADER", 15) == 0) main_header_flag = FALSE;
    if (strncmp(line, "ONUM ", 5) == 0) sscanf(line, "%*s %d", &n_outlines);
    if (strncmp(line, "COL_NUM", 7) == 0) sscanf(line, "%*s %d", &n_cols);
    if (strncmp(line, "TYPE", 4) == 0) {
      strcpy(global_type, &(line[5]));
      c = strrchr(global_type, '\n');
      if (c != NULL) *c = '\0';
    }
  }

  if (n_outlines == -1) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "bad or undefined ONUM in otl file %d"));
  }

  for (i = 0; i < n_outlines; i++) {
    if (feof(fp)) {
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF (ready for next outline) in otl file"));
    }

    gdf_header_flag = FALSE;

    while (!gdf_header_flag) {
      if (feof(fp)) {
        errno = 0;
        ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF (searching for gdf header) in otl file"));
      }
      if (!fgets(line, STRLEN, fp) && ferror(fp)){
        ErrorPrintf(ERROR_BADFILE, "list_labels_in_otl_file(): could not read file");
      }
      if (strncmp(line, "START GDF HEADER", 16) == 0) gdf_header_flag = TRUE;
    }

    n_rows = -1;
    seed_x = seed_y = -1;
    label[0] = '\0';
    type[0] = '\0';

    empty_label_flag = 0;

    while (gdf_header_flag) {
      if (feof(fp)) {
        errno = 0;
        ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF (in gdf header) in otl file"));
      }
      if (!fgets(line, STRLEN, fp) && ferror(fp)) {
        ErrorPrintf(ERROR_BADFILE, "list_labels_in_otl_file(): could not read file"); 
      }
      if (strncmp(line, "END GDF HEADER", 14) == 0) gdf_header_flag = FALSE;
      // getting rows, cols, type
      if (strncmp(line, "ROW_NUM", 7) == 0) sscanf(line, "%*s %d", &n_rows);
      if (strncmp(line, "COL_NUM", 7) == 0) sscanf(line, "%*s %d", &n_cols);
      if (strncmp(line, "TYPE", 4) == 0) {
        strcpy(type, &(line[5]));
        c = strrchr(type, '\n');
        if (c != NULL) *c = '\0';
      }
      if (strncmp(line, "SEED", 4) == 0) sscanf(line, "%*s %d %d", &seed_x, &seed_y);
      if (strncmp(line, "LABEL", 5) == 0) {
        strcpy(label, &(line[6]));
        c = strrchr(label, '\n');
        if (c != NULL) *c = '\0';

        printf("%d: %s\n", i, label);
      }
    }

    if (n_rows < 0) {
      errno = 0;
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "bad or undefined ROW_NUM in otl file"));
    }

    if (n_cols != 2) {
      errno = 0;
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "bad or undefined COL_NUM in otl file"));
    }

    if (label[0] == '\0') {
      empty_label_flag = 1;
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "empty LABEL in otl file (outline %d)", i);
    }

    if (seed_x < 0 || seed_x >= 512 || seed_y < 0 || seed_y >= 512) {
      errno = 0;
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "bad or undefined SEED in otl file"));
    }

    if (type[0] == '\0') strcpy(type, global_type);

    if (strcmp(type, "ascii short") == 0)
      ascii_short_flag = TRUE;
    else if (strcmp(type, "short") == 0)
      ascii_short_flag = FALSE;
    else if (type[0] == '\0') {
      errno = 0;
      ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "undefined TYPE in otl file"));
    }
    else {
      errno = 0;
      ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "unsupported TYPE \"%s\" in otl file", type));
    }

    do {
      if (!fgets(line, STRLEN, fp) && ferror(fp)) {
        ErrorPrintf(ERROR_BADFILE, "list_labels_in_otl_file(): could not read file");
      }
      if (feof(fp)) {
        errno = 0;
        ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "premature EOF (searching for points) in otl file"));
      }
    } while (strncmp(line, "START POINTS", 12) != 0);

    points = (short *)malloc(2 * n_rows * sizeof(short));
    if (points == NULL) {
      errno = 0;
      ErrorReturn(ERROR_NOMEMORY, (ERROR_NOMEMORY, "error allocating memory for points in otl file"));
    }

    if (ascii_short_flag) {
      for (row = 0; row < n_rows; row++) {
        if (!fgets(line, STRLEN, fp) && ferror(fp)) {
          ErrorPrintf(ERROR_BADFILE, "list_labels_in_otl_file(): could not read file");
        }
        if (feof(fp)) {
          free(points);
          errno = 0;
          ErrorReturn(ERROR_BADFILE,
                      (ERROR_BADFILE,
                       "premature end of file while reading "
                       "points from otl file"));
        }
        sscanf(line, "%hd %hd", &(points[2 * row]), &(points[2 * row + 1]));
      }
    }
    else {
      n_read = fread(points, 2, n_rows * 2, fp);
      if (n_read != n_rows * 2) {
        free(points);
        errno = 0;
        ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "error reading points from otl file"));
      }
#if (BYTE_ORDER == LITTLE_ENDIAN)
#if defined(SunOS)
      swab((const char *)points, (char *)points, 2 * n_rows * sizeof(short));
#else
      {
	std::vector<short> tmp(2*n_rows);
	// Note:
	// void swab(const void *from, void *to, ssize_t n);
	// void *memcpy(void *dest, const void *src, size_t n);
	// Because consistency is the hobgoblin of small minds...
	swab(points, tmp.data(), 2 * n_rows * sizeof(short));
	memcpy(points, tmp.data(), 2 * n_rows * sizeof(short));
      }
#endif
#endif
    }

    if (!fgets(line, STRLEN, fp) && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "list_labels_in_otl_file(): could not read file");
    }
    if (strncmp(line, "END POINTS", 10) != 0) {
      free(points);
      errno = 0;
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE,
                   "error after points (\"END POINTS\" expected "
                   "but not found) in otl file"));
    }

    free(points);
  }

  return (NO_ERROR);

} /* end list_labels_in_otl_file() */

MRI *MRIreadOtl(const char *fname, int width, int height, int slices, const char *color_file_name, int flags)
{
  char stem[STRLEN];
  int i;
  MRI *mri;
  char *c;
  int one_file_exists;
  char first_name[STRLEN], last_name[STRLEN];
  FILE *fp;
  COLOR_TABLE *ctab;
  int read_volume_flag, fill_flag, translate_labels_flag, zero_outlines_flag;

  /* ----- set local flags ----- */

  read_volume_flag = FALSE;
  fill_flag = FALSE;
  translate_labels_flag = FALSE;
  zero_outlines_flag = FALSE;

  if (flags & READ_OTL_READ_VOLUME_FLAG) read_volume_flag = TRUE;

  if (flags & READ_OTL_FILL_FLAG) fill_flag = TRUE;

  if (flags & READ_OTL_TRANSLATE_LABELS_FLAG) translate_labels_flag = TRUE;

  if (flags & READ_OTL_ZERO_OUTLINES_FLAG) zero_outlines_flag = TRUE;

  /* ----- reset our unknown label list ----- */

  clear_unknown_labels();

  /* ----- defaults to width and height ----- */
  if (width <= 0) width = 512;
  if (height <= 0) height = 512;

  /* ----- strip the stem of the otl file name ----- */
  strcpy(stem, fname);
  c = strrchr(stem, '.');
  if (c == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "MRIreadOtl(): bad file name: %s", fname));
  }

  for (c--; c >= stem && isdigit(*c); c--)
    ;

  if (c < stem) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "MRIreadOtl(): bad file name: %s", fname));
  }

  c++;

  one_file_exists = FALSE;
  for (i = 1; i <= slices; i++) {
    sprintf(c, "%d.otl", i);
    if (FileExists(stem)) one_file_exists = TRUE;
    if (i == 1) strcpy(first_name, stem);
  }

  strcpy(last_name, stem);

  if (!one_file_exists) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "MRIreadOtl(): couldn't find any file between %s and %s", first_name, last_name));
  }

  if (!read_volume_flag) {
    mri = MRIallocHeader(width, height, slices, MRI_SHORT, 1);
    if (mri == NULL) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_NOMEMORY, "MRIreadOtl(): error allocating MRI structure"));
    }
    return (mri);
  }
  mri = MRIalloc(width, height, slices, MRI_SHORT);
  if (mri == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_NOMEMORY, "MRIreadOtl(): error allocating MRI structure"));
  }

  if ((ctab = CTABreadASCII(color_file_name)) == NULL) {
    MRIfree(&mri);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "MRIreadOtl(): error reading color file %s", color_file_name));
  }

  one_file_exists = FALSE;
  for (i = 1; i <= slices; i++) {
    sprintf(c, "%d.otl", i);
    if ((fp = fopen(stem, "r")) != NULL) {
      if (read_otl_file(fp, mri, i, ctab, fill_flag, translate_labels_flag, zero_outlines_flag) != NO_ERROR) {
        MRIfree(&mri);
        return (NULL);
      }
      one_file_exists = TRUE;
    }

    exec_progress_callback(i - 1, slices, 0, 1);
  }

  CTABfree(&ctab);

  if (!one_file_exists) {
    MRIfree(&mri);
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "MRIreadOtl(): found at least one file "
                 "between %s and %s but couldn't open it!",
                 first_name,
                 last_name));
  }

  if (n_unknown_labels == 0) {
    printf("no unknown labels\n");
  }
  else {
    printf("unknown labels:\n");
    print_unknown_labels("  ");
  }
  clear_unknown_labels();

  // default direction cosine set here to be CORONAL
  // no orientation info is given and thus set to coronal
  setDirectionCosine(mri, MRI_CORONAL);

  return (mri);

} /* end MRIreadOtl() */

#define XIMG_PIXEL_DATA_OFFSET 8432
#define XIMG_IMAGE_HEADER_OFFSET 2308

static MRI *ximgRead(const char *fname, int read_volume)
{
  char fname_format[STRLEN];
  char fname_dir[STRLEN];
  char fname_base[STRLEN];
  char *c;
  MRI *mri = NULL;
  int im_init;
  int im_low, im_high;
  char fname_use[STRLEN];
  char temp_string[STRLEN];
  FILE *fp;
  int width, height;
  int pixel_data_offset;
  int image_header_offset;
  float tl_r, tl_a, tl_s;
  float tr_r, tr_a, tr_s;
  float br_r, br_a, br_s;
  float c_r, c_a, c_s;
  float n_r, n_a, n_s;
  float xlength, ylength, zlength;
  int i, y;
  MRI *header;
  float xfov, yfov, zfov;
  float nlength;

  printf("XIMG: using XIMG header corrections\n");

  /* ----- check the first (passed) file ----- */
  if (!FileExists(fname)) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error opening file %s", fname));
  }

  /* ----- split the file name into name and directory ----- */
  const char *cc = strrchr(fname, '/');
  if (cc == NULL) {
    fname_dir[0] = '\0';
    strcpy(fname_base, fname);
  }
  else {
    strncpy(fname_dir, fname, (cc - fname + 1));
    fname_dir[cc - fname + 1] = '\0';
    strcpy(fname_base, cc + 1);
  }

  /* ----- derive the file name format (for sprintf) ----- */
  if (strncmp(fname_base, "I.", 2) == 0) {
    im_init = atoi(&fname_base[2]);
    sprintf(fname_format, "I.%%03d");
  }
  else if (strlen(fname_base) >= 3) /* avoid core dumps below... */
  {
    c = &fname_base[strlen(fname_base) - 3];
    if (strcmp(c, ".MR") == 0) {
      *c = '\0';
      for (c--; isdigit(*c) && c >= fname_base; c--)
        ;
      c++;
      im_init = atoi(c);
      *c = '\0';
      int req = snprintf(fname_format, STRLEN, "%s%%d.MR", fname_base);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    else {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADPARM, "genesisRead(): can't determine file name format for %s", fname));
    }
  }
  else {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "genesisRead(): can't determine file name format for %s", fname));
  }

  int req = snprintf(temp_string, STRLEN, "%s", fname_format);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  req = snprintf(fname_format, STRLEN, "%s%s", fname_dir, temp_string);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  /* ----- find the low and high files ----- */
  im_low = im_init;
  do {
    im_low--;
    sprintf(fname_use, fname_format, im_low);
  } while (FileExists(fname_use));
  im_low++;

  im_high = im_init;
  do {
    im_high++;
    sprintf(fname_use, fname_format, im_high);
  } while (FileExists(fname_use));
  im_high--;

  /* ----- allocate the mri structure ----- */
  header = MRIallocHeader(1, 1, 1, MRI_SHORT, 1);

  header->depth = im_high - im_low + 1;
  header->imnr0 = 1;
  header->imnr1 = header->depth;

  /* ----- get the header information from the first file ----- */
  sprintf(fname_use, fname_format, im_low);
  if ((fp = fopen(fname_use, "r")) == NULL) {
    MRIfree(&header);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error opening file %s\n", fname_use));
  }

  fseek(fp, 4, SEEK_SET);
  if (fread(&pixel_data_offset, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  pixel_data_offset = orderIntBytes(pixel_data_offset);
  printf("XIMG: pixel data offset is %d, ", pixel_data_offset);
  if (pixel_data_offset != XIMG_PIXEL_DATA_OFFSET) pixel_data_offset = XIMG_PIXEL_DATA_OFFSET;
  printf("using offset %d\n", pixel_data_offset);

  fseek(fp, 8, SEEK_SET);
  if (fread(&width, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  width = orderIntBytes(width);
  if (fread(&height, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  height = orderIntBytes(height);

  fseek(fp, 148, SEEK_SET);
  if (fread(&image_header_offset, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  image_header_offset = orderIntBytes(image_header_offset);
  printf("XIMG: image header offset is %d, ", image_header_offset);
  if (image_header_offset != XIMG_IMAGE_HEADER_OFFSET) image_header_offset = XIMG_IMAGE_HEADER_OFFSET;
  printf("using offset %d\n", image_header_offset);

  header->width = width;
  header->height = height;

  strcpy(header->fname, fname);

  fseek(fp, image_header_offset + 26 + 2, SEEK_SET);
  if (fread(&(header->thick), 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  header->thick = orderFloatBytes(header->thick);
  header->zsize = header->thick;

  fseek(fp, image_header_offset + 50 + 2, SEEK_SET);
  if (fread(&(header->xsize), 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  header->xsize = orderFloatBytes(header->xsize);
  if (fread(&(header->ysize), 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  header->ysize = orderFloatBytes(header->ysize);
  header->ps = header->xsize;

/* all in micro-seconds */
#define MICROSECONDS_PER_MILLISECOND 1e3
  fseek(fp, image_header_offset + 194 + 6, SEEK_SET);
  header->tr = freadInt(fp) / MICROSECONDS_PER_MILLISECOND;
  fseek(fp, image_header_offset + 198, SEEK_SET);
  header->ti = freadInt(fp) / MICROSECONDS_PER_MILLISECOND;
  fseek(fp, image_header_offset + 202, SEEK_SET);
  header->te = freadInt(fp) / MICROSECONDS_PER_MILLISECOND;
  fseek(fp, image_header_offset + 254, SEEK_SET);
  header->flip_angle = RADIANS(freadShort(fp)); /* was in degrees */

  fseek(fp, image_header_offset + 130 + 6, SEEK_SET);
  if (fread(&c_r, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  c_r = orderFloatBytes(c_r);
  if (fread(&c_a, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  c_a = orderFloatBytes(c_a);
  if (fread(&c_s, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  c_s = orderFloatBytes(c_s);
  if (fread(&n_r, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  n_r = orderFloatBytes(n_r);
  if (fread(&n_a, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  n_a = orderFloatBytes(n_a);
  if (fread(&n_s, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  n_s = orderFloatBytes(n_s);
  if (fread(&tl_r, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tl_r = orderFloatBytes(tl_r);
  if (fread(&tl_a, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tl_a = orderFloatBytes(tl_a);
  if (fread(&tl_s, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tl_s = orderFloatBytes(tl_s);
  if (fread(&tr_r, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tr_r = orderFloatBytes(tr_r);
  if (fread(&tr_a, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tr_a = orderFloatBytes(tr_a);
  if (fread(&tr_s, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  tr_s = orderFloatBytes(tr_s);
  if (fread(&br_r, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  br_r = orderFloatBytes(br_r);
  if (fread(&br_a, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  br_a = orderFloatBytes(br_a);
  if (fread(&br_s, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "genesisRead(): could not read file");
  }
  br_s = orderFloatBytes(br_s);

  nlength = sqrt(n_r * n_r + n_a * n_a + n_s * n_s);
  n_r = n_r / nlength;
  n_a = n_a / nlength;
  n_s = n_s / nlength;

  if (getenv("KILLIANY_SWAP") != NULL) {
    printf("WARNING - swapping normal direction!\n");
    n_a *= -1;
  }

  header->x_r = (tr_r - tl_r);
  header->x_a = (tr_a - tl_a);
  header->x_s = (tr_s - tl_s);
  header->y_r = (br_r - tr_r);
  header->y_a = (br_a - tr_a);
  header->y_s = (br_s - tr_s);

  /* --- normalize -- the normal vector from the file
     should have length 1, but just in case... --- */
  xlength = sqrt(header->x_r * header->x_r + header->x_a * header->x_a + header->x_s * header->x_s);
  ylength = sqrt(header->y_r * header->y_r + header->y_a * header->y_a + header->y_s * header->y_s);
  zlength = sqrt(n_r * n_r + n_a * n_a + n_s * n_s);

  header->x_r = header->x_r / xlength;
  header->x_a = header->x_a / xlength;
  header->x_s = header->x_s / xlength;
  header->y_r = header->y_r / ylength;
  header->y_a = header->y_a / ylength;
  header->y_s = header->y_s / ylength;
  header->z_r = n_r / zlength;
  header->z_a = n_a / zlength;
  header->z_s = n_s / zlength;

  header->c_r = (tl_r + br_r) / 2.0 + n_r * header->zsize * (header->depth - 1.0) / 2.0;
  header->c_a = (tl_a + br_a) / 2.0 + n_a * header->zsize * (header->depth - 1.0) / 2.0;
  header->c_s = (tl_s + br_s) / 2.0 + n_s * header->zsize * (header->depth - 1.0) / 2.0;

  header->ras_good_flag = 1;

  header->xend = header->xsize * (double)header->width / 2.0;
  header->xstart = -header->xend;
  header->yend = header->ysize * (double)header->height / 2.0;
  header->ystart = -header->yend;
  header->zend = header->zsize * (double)header->depth / 2.0;
  header->zstart = -header->zend;

  xfov = header->xend - header->xstart;
  yfov = header->yend - header->ystart;
  zfov = header->zend - header->zstart;

  header->fov = (xfov > yfov ? (xfov > zfov ? xfov : zfov) : (yfov > zfov ? yfov : zfov));

  fclose(fp);

  if (read_volume)
    mri = MRIalloc(header->width, header->height, header->depth, header->type);
  else
    mri = MRIallocHeader(header->width, header->height, header->depth, header->type, 1);

  MRIcopyHeader(header, mri);
  MRIfree(&header);

  /* ----- read the volume if required ----- */
  if (read_volume) {
    for (i = im_low; i <= im_high; i++) {
      sprintf(fname_use, fname_format, i);
      if ((fp = fopen(fname_use, "r")) == NULL) {
        MRIfree(&mri);
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error opening file %s", fname_use));
      }

      /*
        fseek(fp, 4, SEEK_SET);
        fread(&pixel_data_offset, 4, 1, fp);
        pixel_data_offset = orderIntBytes(pixel_data_offset);
      */
      pixel_data_offset = XIMG_PIXEL_DATA_OFFSET;
      fseek(fp, pixel_data_offset, SEEK_SET);

      for (y = 0; y < mri->height; y++) {
        if ((int)fread(mri->slices[i - im_low][y], 2, mri->width, fp) != mri->width) {
          fclose(fp);
          MRIfree(&mri);
          errno = 0;
          ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(): error reading from file file %s", fname_use));
        }
#if (BYTE_ORDER == LITTLE_ENDIAN)
	{
	  std::vector<short> tmp(2*mri->width);
	  // Note:
	  // void swab(const void *from, void *to, ssize_t n);
	  // void *memcpy(void *dest, const void *src, size_t n);
	  // Because consistency is the hobgoblin of small minds...
	  swab(mri->slices[i - im_low][y], tmp.data(), (size_t)(2 * mri->width));
	  memcpy(mri->slices[i - im_low][y], tmp.data(), (size_t)(2 * mri->width));
	}
#endif
      }

      fclose(fp);

      exec_progress_callback(i - im_low, im_high - im_low + 1, 0, 1);
    }
  }

  return (mri);

} /* end ximgRead() */

/*-----------------------------------------------------------
  MRISreadCurvAsMRI() - reads freesurfer surface curv format
  as an MRI.
  -----------------------------------------------------------*/
static MRI *MRISreadCurvAsMRI(const char *curvfile, int read_volume)
{
  int magno, k, vnum, fnum, vals_per_vertex;
  float curv;
  FILE *fp;
  MRI *curvmri;

  if (!IDisCurv(curvfile)) return (NULL);

  fp = fopen(curvfile, "r");
  fread3(&magno, fp);

  vnum = freadInt(fp);
  fnum = freadInt(fp);
  vals_per_vertex = freadInt(fp);
  if (vals_per_vertex != 1) {
    fclose(fp);
    printf("ERROR: MRISreadCurvAsMRI: %s, vals/vertex %d unsupported\n", curvfile, vals_per_vertex);
    return (NULL);
  }

  if (!read_volume) {
    curvmri = MRIallocHeader(vnum, 1, 1, MRI_FLOAT, 1);
    curvmri->nframes = 1;
    fclose(fp);
    return (curvmri);
  }

  curvmri = MRIalloc(vnum, 1, 1, MRI_FLOAT);
  for (k = 0; k < vnum; k++) {
    curv = freadFloat(fp);
    MRIsetVoxVal(curvmri, k, 0, 0, 0, curv);
  }
  fclose(fp);

  return (curvmri);
}

/*------------------------------------------------------------------
  nifti1Read() - note: there is also an niiRead(). Make sure to edit
  both. NIFTI1 is defined to be the two-file NIFTI standard, ie, there
  must be a .img and .hdr file. (see is_nifti1(char *fname))
  Automatically detects whether an input is Ico7 and
  reshapes.
   -----------------------------------------------------------------*/
static MRI *nifti1Read(const char *fname, int read_volume)
{
  char hdr_fname[STRLEN];
  char img_fname[STRLEN];
  char fname_stem[STRLEN];
  char *dot;
  FILE *fp;
  MRI *mri, *mritmp;
  struct nifti_1_header hdr;
  int nslices;
  int fs_type;
  float time_units_factor, space_units_factor;
  int swapped_flag;
  int n_read, i, j, k, t;
  int bytes_per_voxel, time_units, space_units;
  int ncols, IsIco7 = 0;

  strcpy(fname_stem, fname);
  dot = strrchr(fname_stem, '.');
  if (dot != NULL)
    if (strcmp(dot, ".img") == 0 || strcmp(dot, ".hdr") == 0) *dot = '\0';

  int req = snprintf(hdr_fname, STRLEN, "%s.hdr", fname_stem);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  req = snprintf(img_fname, STRLEN, "%s.img", fname_stem);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  fp = fopen(hdr_fname, "r");
  if (fp == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error opening file %s", hdr_fname));
  }

  if (fread(&hdr, sizeof(hdr), 1, fp) != 1) {
    fclose(fp);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error reading header from %s", hdr_fname));
  }

  fclose(fp);

  swapped_flag = FALSE;
  if (hdr.dim[0] < 1 || hdr.dim[0] > 7) {
    swapped_flag = TRUE;
    swap_nifti_1_header(&hdr);
    if (hdr.dim[0] < 1 || hdr.dim[0] > 7) {
      ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): bad number of dimensions (%hd) in %s", hdr.dim[0], hdr_fname));
    }
  }

  if (memcmp(hdr.magic, NIFTI1_MAGIC, 4) != 0)
    ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): bad magic number in %s", hdr_fname));

  if (hdr.dim[0] != 2 && hdr.dim[0] != 3 && hdr.dim[0] != 4)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "nifti1Read(): %hd dimensions in %s; unsupported", hdr.dim[0], hdr_fname));

  if (hdr.datatype == DT_NONE || hdr.datatype == DT_UNKNOWN)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "nifti1Read(): unknown or no data type in %s; bailing out", hdr_fname));

  space_units = XYZT_TO_SPACE(hdr.xyzt_units);
  if (space_units == NIFTI_UNITS_METER)
    space_units_factor = 1000.0;
  else if (space_units == NIFTI_UNITS_MM)
    space_units_factor = 1.0;
  else if (space_units == NIFTI_UNITS_MICRON)
    space_units_factor = 0.001;
  else if (space_units == NIFTI_UNITS_UNKNOWN) {
    space_units_factor = 1.0;
  }
  else {
    ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): unknown space units %d in %s", space_units, hdr_fname));
  }

  time_units = XYZT_TO_TIME(hdr.xyzt_units);
  if (time_units == NIFTI_UNITS_SEC)
    time_units_factor = 1000.0;
  else if (time_units == NIFTI_UNITS_MSEC)
    time_units_factor = 1.0;
  else if (time_units == NIFTI_UNITS_USEC)
    time_units_factor = 0.001;
  else {
    if (hdr.dim[4] > 1) {
      ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): unknown time units %d in %s", time_units, hdr_fname));
    }
    else
      time_units_factor = 0;
  }

  if (hdr.slice_code != 0 && DIM_INFO_TO_SLICE_DIM(hdr.dim_info) != 0 && hdr.slice_duration > 0.0) {
    if (hdr.slice_code != NIFTI_SLICE_SEQ_INC && hdr.slice_code != NIFTI_SLICE_SEQ_DEC &&
        hdr.slice_code != NIFTI_SLICE_ALT_INC && hdr.slice_code != NIFTI_SLICE_ALT_DEC &&
        hdr.slice_code != NIFTI_SLICE_ALT_INC2 && hdr.slice_code != NIFTI_SLICE_ALT_DEC2) {
      ErrorReturn(
          NULL,
          (ERROR_UNSUPPORTED, "nifti1Read(): unsupported slice timing pattern %d in %s", hdr.slice_code, hdr_fname));
    }
  }

  if (hdr.dim[0] < 4)
    nslices = 1;
  else
    nslices = hdr.dim[4];

  // check whether data needs to be scaled
  bool scaledata = (hdr.scl_slope != 0) && !((hdr.scl_slope == 1) && (hdr.scl_inter == 0));

  if (!scaledata) {
    // voxel values are unscaled -- we use the file's data type
    if (hdr.datatype == DT_UNSIGNED_CHAR) {
      fs_type = MRI_UCHAR;
      bytes_per_voxel = 1;
    }
    else if (hdr.datatype == DT_SIGNED_SHORT) {
      fs_type = MRI_SHORT;
      bytes_per_voxel = 2;
    }
    else if (hdr.datatype == DT_UINT16) {
      // This will not always work ...
      printf("INFO: This is an unsigned short.\n");    // I'll try to read it, but\n");
      //printf("      it might not work if there are values over 32k\n");
      fs_type = MRI_USHRT;
      bytes_per_voxel = 2;
    }
    else if (hdr.datatype == DT_SIGNED_INT) {
      fs_type = MRI_INT;
      bytes_per_voxel = 4;
    }
    else if (hdr.datatype == DT_FLOAT) {
      fs_type = MRI_FLOAT;
      bytes_per_voxel = 4;
    }
    else {
      ErrorReturn(NULL,
                  (ERROR_UNSUPPORTED,
                   "nifti1Read(): unsupported datatype %d "
                   "(with scl_slope = 0) in %s",
                   hdr.datatype,
                   hdr_fname));
    }
  }
  else  // we must scale the voxel values
  {
    if (hdr.datatype != DT_UNSIGNED_CHAR && hdr.datatype != DT_SIGNED_SHORT && hdr.datatype != DT_SIGNED_INT &&
        hdr.datatype != DT_FLOAT && hdr.datatype != DT_DOUBLE && hdr.datatype != DT_INT8 && hdr.datatype != DT_UINT16 &&
        hdr.datatype != DT_UINT32) {
      ErrorReturn(NULL,
                  (ERROR_UNSUPPORTED,
                   "nifti1Read(): unsupported datatype %d "
                   "(with scl_slope != 0) in %s",
                   hdr.datatype,
                   hdr_fname));
    }
    fs_type = MRI_FLOAT;
    bytes_per_voxel = 0; /* set below -- this line is to
                                              avoid the compiler warning */
  }

  // Check whether dim[1] is less than 0. This can happen when FreeSurfer
  // writes a nifti volume where the number of columns is more than shortmax.
  // In this case, dim[1] is set to -1, and glmin is set to the number of cols.
  if (hdr.dim[1] > 0)
    ncols = hdr.dim[1];
  else
    ncols = hdr.glmin;

  if (ncols * hdr.dim[2] * hdr.dim[3] == 163842) IsIco7 = 1;

  if (read_volume)
    mri = MRIallocSequence(ncols, hdr.dim[2], hdr.dim[3], fs_type, nslices);
  else {
    if (!IsIco7)
      mri = MRIallocHeader(ncols, hdr.dim[2], hdr.dim[3], fs_type, nslices);
    else
      mri = MRIallocHeader(163842, 1, 1, fs_type, nslices);
    mri->nframes = nslices;
  }
  if (mri == NULL) return (NULL);

  mri->xsize = hdr.pixdim[1];
  mri->ysize = hdr.pixdim[2];
  mri->zsize = hdr.pixdim[3];
  // Keep in msec as NIFTI_UNITS_MSEC is specified
  if (hdr.dim[0] == 4) mri->tr = hdr.pixdim[4];

  // Set the vox2ras matrix
  if (hdr.sform_code != 0) {
    // First, use the sform, if that is ok. Using the sform
    // first makes it more compatible with FSL.
    // fprintf(stderr, "INFO: using NIfTI-1 sform \n");
    if (niftiSformToMri(mri, &hdr) != NO_ERROR) {
      MRIfree(&mri);
      return (NULL);
    }
    mri->ras_good_flag = 1;
  }
  else if (hdr.qform_code != 0) {
    // Then, try the qform, if that is ok
    fprintf(stderr, "INFO: using NIfTI-1 qform \n");
    if (niftiQformToMri(mri, &hdr) != NO_ERROR) {
      MRIfree(&mri);
      return (NULL);
    }
    mri->ras_good_flag = 1;
  }
  else {
    // Should probably just die here.
    printf("WARNING: neither NIfTI-1 qform or sform are valid\n");
    printf("WARNING: your volume will probably be incorrectly oriented\n");
    mri->x_r = -1.0;
    mri->x_a = 0.0;
    mri->x_s = 0.0;
    mri->y_r = 0.0;
    mri->y_a = 1.0;
    mri->y_s = 0.0;
    mri->z_r = 0.0;
    mri->z_a = 0.0;
    mri->z_s = 1.0;
    mri->c_r = mri->xsize * mri->width / 2.0;
    mri->c_a = mri->ysize * mri->height / 2.0;
    mri->c_s = mri->zsize * mri->depth / 2.0;
    mri->ras_good_flag = 0;
  }

  mri->xsize = mri->xsize * space_units_factor;
  mri->ysize = mri->ysize * space_units_factor;
  mri->zsize = mri->zsize * space_units_factor;
  mri->c_r = mri->c_r * space_units_factor;
  mri->c_a = mri->c_a * space_units_factor;
  mri->c_s = mri->c_s * space_units_factor;
  if (hdr.dim[0] == 4) mri->tr = mri->tr * time_units_factor;

  if (!read_volume) return (mri);

  fp = fopen(img_fname, "r");
  if (fp == NULL) {
    MRIfree(&mri);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error opening file %s", img_fname));
  }

  if (!scaledata)  // no voxel value scaling needed
  {
    void *buf;

    for (t = 0; t < mri->nframes; t++)
      for (k = 0; k < mri->depth; k++) {
        for (j = 0; j < mri->height; j++) {
          buf = &MRIseq_vox(mri, 0, j, k, t);

          n_read = fread(buf, bytes_per_voxel, mri->width, fp);
          if (n_read != mri->width) {
            fclose(fp);
            MRIfree(&mri);
            errno = 0;
            ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error reading from %s", img_fname));
          }

          if (swapped_flag) {
            if (bytes_per_voxel == 2) byteswapbufshort(buf, bytes_per_voxel * mri->width);
            if (bytes_per_voxel == 4) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
          }
        }
        exec_progress_callback(k, mri->depth, t, mri->nframes);
      }
  }
  else  // voxel value scaling needed
  {
    if (hdr.datatype == DT_UNSIGNED_CHAR) {
      unsigned char *buf;
      bytes_per_voxel = 1;
      buf = (unsigned char *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = fread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              fclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error reading from %s", img_fname));
            }
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_SIGNED_SHORT) {
      short *buf;
      bytes_per_voxel = 2;
      buf = (short *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = fread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              fclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error reading from %s", img_fname));
            }
            if (swapped_flag) byteswapbufshort(buf, bytes_per_voxel * mri->width);
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_SIGNED_INT) {
      int *buf;
      bytes_per_voxel = 4;
      buf = (int *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = fread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              fclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error reading from %s", img_fname));
            }
            if (swapped_flag) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_FLOAT) {
      float *buf;
      bytes_per_voxel = 4;
      buf = (float *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = fread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              fclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error reading from %s", img_fname));
            }
            if (swapped_flag) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_DOUBLE) {
      double *buf;
      unsigned char *cbuf, ccbuf[8];
      bytes_per_voxel = 8;
      buf = (double *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = fread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              fclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error reading from %s", img_fname));
            }
            if (swapped_flag) {
              for (i = 0; i < mri->width; i++) {
                cbuf = (unsigned char *)&buf[i];
                memmove(ccbuf, cbuf, 8);
                cbuf[0] = ccbuf[7];
                cbuf[1] = ccbuf[6];
                cbuf[2] = ccbuf[5];
                cbuf[3] = ccbuf[4];
                cbuf[4] = ccbuf[3];
                cbuf[5] = ccbuf[2];
                cbuf[6] = ccbuf[1];
                cbuf[7] = ccbuf[0];
              }
            }
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_INT8) {
      char *buf;
      bytes_per_voxel = 1;
      buf = (char *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = fread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              fclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error reading from %s", img_fname));
            }
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_UINT16) {
      unsigned short *buf;
      bytes_per_voxel = 2;
      buf = (unsigned short *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = fread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              fclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error reading from %s", img_fname));
            }
            if (swapped_flag) byteswapbufshort(buf, bytes_per_voxel * mri->width);
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_UINT32) {
      unsigned int *buf;
      bytes_per_voxel = 4;
      buf = (unsigned int *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = fread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              fclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "nifti1Read(): error reading from %s", img_fname));
            }
            if (swapped_flag) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }
  }

  fclose(fp);

  // Check for ico7 surface
  if (IsIco7) {
    // printf("niiRead: reshaping\n");
    mritmp = mri_reshape(mri, 163842, 1, 1, mri->nframes);
    MRIfree(&mri);
    mri = mritmp;
  }

  return (mri);

} /* end nifti1Read() */

/*------------------------------------------------------------------
  nifti1Write() - note: there is also an niiWrite(). Make sure to
  edit both. Automatically detects whether an input is Ico7
  and reshapes.
  -----------------------------------------------------------------*/
static int nifti1Write(MRI *mri0, const char *fname)
{
  FILE *fp;
  int j, k, t;
  BUFTYPE *buf;
  struct nifti_1_header hdr;
  char fname_stem[STRLEN];
  char hdr_fname[STRLEN];
  char img_fname[STRLEN];
  char *dot;
  int error;
  int shortmax;
  MRI *mri;
  int FreeMRI = 0;

  // Check for ico7 surface
  if (mri0->width == 163842 && mri0->height == 1 && mri0->depth == 1) {
    // printf("nifit1Write: reshaping\n");
    mri = mri_reshape(mri0, 27307, 1, 6, mri0->nframes);
    FreeMRI = 1;
  }
  else
    mri = mri0;

  shortmax = (int)(pow(2.0, 15.0));
  if (0 && mri->width > shortmax) {
    printf("NIFTI FORMAT WARNING: ncols %d in input exceeds %d.\n", mri->width, shortmax);
    printf("So I'm going to put the true ncols in glmin and set dim[1]=-1.\n");
    printf("This should be ok within FreeSurfer, but you will not be\n");
    printf("able to use this volume with other software.\n");
    // exit(1);
  }
  if (mri->height > shortmax) {
    printf("NIFTI FORMAT ERROR: nrows %d in volume exceeds %d\n", mri->height, shortmax);
    exit(1);
  }
  if (mri->depth > shortmax) {
    printf("NIFTI FORMAT ERROR: nslices %d in volume exceeds %d\n", mri->depth, shortmax);
    exit(1);
  }
  if (mri->nframes > shortmax) {
    printf("NIFTI FORMAT ERROR:  nframes %d in volume exceeds %d\n", mri->nframes, shortmax);
    exit(1);
  }

  memset(&hdr, 0x00, sizeof(hdr));

  hdr.sizeof_hdr = 348;
  hdr.dim_info = 0;

  for (t = 0; t < 8; t++) {
    hdr.dim[t] = 1;
    hdr.pixdim[t] = 1;
  }  // needed for afni
  if (mri->nframes == 1)
    hdr.dim[0] = 3;
  else
    hdr.dim[0] = 4;

  if (mri->width < shortmax)
    hdr.dim[1] = mri->width;
  else {
    // number of columns too big, put in glmin
    hdr.dim[1] = -1;
    hdr.glmin = mri->width;
  }
  hdr.dim[2] = mri->height;
  hdr.dim[3] = mri->depth;
  hdr.dim[4] = mri->nframes;
  hdr.pixdim[1] = mri->xsize;
  hdr.pixdim[2] = mri->ysize;
  hdr.pixdim[3] = mri->zsize;
  hdr.pixdim[4] = mri->tr / 1000.0;  // see also xyzt_units

  if (mri->type == MRI_UCHAR) {
    hdr.datatype = DT_UNSIGNED_CHAR;
    hdr.bitpix = 8;
  }
  else if (mri->type == MRI_INT) {
    hdr.datatype = DT_SIGNED_INT;
    hdr.bitpix = 32;
  }
  else if (mri->type == MRI_LONG) {
    hdr.datatype = DT_SIGNED_INT;
    hdr.bitpix = 32;
  }
  else if (mri->type == MRI_FLOAT) {
    hdr.datatype = DT_FLOAT;
    hdr.bitpix = 32;
  }
  else if (mri->type == MRI_SHORT) {
    hdr.datatype = DT_SIGNED_SHORT;
    hdr.bitpix = 16;
  }
  else if (mri->type == MRI_USHRT) {
    hdr.datatype = DT_UINT16;
    hdr.bitpix = 16;
  }
  else if (mri->type == MRI_BITMAP) {
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "nifti1Write(): data type MRI_BITMAP unsupported"));
  }
  else if (mri->type == MRI_TENSOR) {
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "nifti1Write(): data type MRI_TENSOR unsupported"));
  }
  else {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "nifti1Write(): unknown data type %d", mri->type));
  }

  hdr.intent_code = NIFTI_INTENT_NONE;
  hdr.intent_name[0] = '\0';
  hdr.vox_offset = 0;
  hdr.scl_slope = 0.0;
  hdr.slice_code = 0;
  hdr.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;  // This may be wrong
  hdr.cal_max = 0.0;
  hdr.cal_min = 0.0;
  hdr.toffset = 0;
  sprintf(hdr.descrip, "FreeSurfer %s", __DATE__);

  /* set the nifti header qform values */
  error = mriToNiftiQform(mri, &hdr);
  if (error != NO_ERROR) return (error);

  /* set the nifti header sform values */
  // This just copies the vox2ras into the sform
  mriToNiftiSform(mri, &hdr);

  memmove(hdr.magic, NIFTI1_MAGIC, 4);

  strcpy(fname_stem, fname);
  dot = strrchr(fname_stem, '.');
  if (dot != NULL)
    if (strcmp(dot, ".img") == 0 || strcmp(dot, ".hdr") == 0) *dot = '\0';

  int needed = snprintf(hdr_fname, STRLEN, "%s.hdr", fname_stem);
  if( needed >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation of hdr_fname" << std::endl;
  }
  needed = snprintf(img_fname, STRLEN, "%s.img", fname_stem);
  if( needed >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation of img_fname" << std::endl;
  }

  fp = fopen(hdr_fname, "w");
  if (fp == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "nifti1Write(): error opening file %s", hdr_fname));
  }

  if (fwrite(&hdr, sizeof(hdr), 1, fp) != 1) {
    fclose(fp);
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "nifti1Write(): error writing header to %s", hdr_fname));
  }

  fclose(fp);

  fp = fopen(img_fname, "w");
  if (fp == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "nifti1Write(): error opening file %s", img_fname));
  }

  for (t = 0; t < mri->nframes; t++)
    for (k = 0; k < mri->depth; k++) {
      for (j = 0; j < mri->height; j++) {
        buf = &MRIseq_vox(mri, 0, j, k, t);
        if ((int)fwrite(buf, hdr.bitpix / 8, mri->width, fp) != mri->width) {
          fclose(fp);
          errno = 0;
          ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "nifti1Write(): error writing data to %s", img_fname));
        }
      }
      exec_progress_callback(k, mri->depth, t, mri->nframes);
    }

  fclose(fp);
  if (FreeMRI) MRIfree(&mri);

  return (NO_ERROR);

} /* end nifti1Write() */

/*------------------------------------------------------------------
  niiRead() - note: there is also an nifti1Read(). Make sure to
  edit both. Automatically detects whether an input is Ico7
  and reshapes.
  -----------------------------------------------------------------*/
static MRI *niiRead(const char *fname, int read_volume)
{
  znzFile fp;
  MRI *mri, *mritmp;
  struct nifti_1_header hdr;
  int nslices;
  int fs_type;
  float time_units_factor, space_units_factor;
  int swapped_flag;
  int n_read, i, j, k, t;
  int bytes_per_voxel, time_units, space_units;
  int use_compression, fnamelen;
  int ncols, IsIco7 = 0;

  use_compression = 0;
  fnamelen = strlen(fname);
  if (fname[fnamelen - 1] == 'z') use_compression = 1;
  if (Gdiag_no > 0) printf("niiRead: use_compression = %d\n", use_compression);

  fp = znzopen(fname, "r", use_compression);
  if (fp == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error opening file %s", fname));
  }

  if (znzread(&hdr, sizeof(hdr), 1, fp) != 1) {
    znzclose(fp);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error reading header from %s", fname));
  }

  znzclose(fp);

  swapped_flag = FALSE;
  if (hdr.dim[0] < 1 || hdr.dim[0] > 7) {
    swapped_flag = TRUE;
    swap_nifti_1_header(&hdr);
    if (hdr.dim[0] < 1 || hdr.dim[0] > 7) {
      ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): bad number of dimensions (%hd) in %s", hdr.dim[0], fname));
    }
  }

  if (memcmp(hdr.magic, NII_MAGIC, 4) != 0) {
    ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): bad magic number in %s", fname));
  }

  //  if (hdr.dim[0] != 2 && hdr.dim[0] != 3 && hdr.dim[0] != 4){
  if (hdr.dim[0] < 1 || hdr.dim[0] > 5) {
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "niiRead(): %hd dimensions in %s; unsupported", hdr.dim[0], fname));
  }

  if (hdr.datatype == DT_NONE || hdr.datatype == DT_UNKNOWN) {
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "niiRead(): unknown or no data type in %s; bailing out", fname));
  }
  if (hdr.dim[4] == 0) {
    printf("WARNING: hdr.dim[4] = 0 (nframes), setting to 1\n");
    hdr.dim[4] = 1;
  }

  space_units = XYZT_TO_SPACE(hdr.xyzt_units);
  if (space_units == NIFTI_UNITS_METER)
    space_units_factor = 1000.0;
  else if (space_units == NIFTI_UNITS_MM)
    space_units_factor = 1.0;
  else if (space_units == NIFTI_UNITS_MICRON)
    space_units_factor = 0.001;
  else if (space_units == NIFTI_UNITS_UNKNOWN) {
    space_units_factor = 1.0;
  }
  else
    ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): unknown space units %d in %s", space_units, fname));

  time_units = XYZT_TO_TIME(hdr.xyzt_units);
  if (time_units == NIFTI_UNITS_SEC)
    time_units_factor = 1000.0;
  else if (time_units == NIFTI_UNITS_MSEC)
    time_units_factor = 1.0;
  else if (time_units == NIFTI_UNITS_USEC)
    time_units_factor = 0.001;
  else {
    // This can be a tricky situation because the time units may not mean
    // anything even with multiple frames (eg, BO mag and phase).
    if (hdr.dim[4] > 1) printf("WARNING: niiRead(): unknown time units %d in %s\n", time_units, fname);
    time_units_factor = 0;
  }
  // printf("hdr.xyzt_units = %d, time_units = %d, %g, %g\n",
  // hdr.xyzt_units,time_units,hdr.pixdim[4],time_units_factor);

  if (hdr.slice_code != 0 && DIM_INFO_TO_SLICE_DIM(hdr.dim_info) != 0 && hdr.slice_duration > 0.0) {
    if (hdr.slice_code != NIFTI_SLICE_SEQ_INC && hdr.slice_code != NIFTI_SLICE_SEQ_DEC &&
        hdr.slice_code != NIFTI_SLICE_ALT_INC && hdr.slice_code != NIFTI_SLICE_ALT_DEC &&
        hdr.slice_code != NIFTI_SLICE_ALT_INC2 && hdr.slice_code != NIFTI_SLICE_ALT_DEC2) {
      ErrorReturn(
          NULL, (ERROR_UNSUPPORTED, "niiRead(): unsupported slice timing pattern %d in %s", hdr.slice_code, fname));
    }
  }

  if (hdr.dim[0] < 4)
    nslices = 1;
  else
    nslices = hdr.dim[4];

  // put extra dims in frames
  if (hdr.dim[0] > 4 && hdr.dim[5] > 0) nslices *= hdr.dim[5];
  if (Gdiag_no > 0)
    printf("niiRead(): hdr.dim %d %d %d %d %d %d\n",
           hdr.dim[0],
           hdr.dim[1],
           hdr.dim[2],
           hdr.dim[3],
           hdr.dim[4],
           hdr.dim[5]);

  // check whether data needs to be scaled
  bool scaledata = (hdr.scl_slope != 0) && !((hdr.scl_slope == 1) && (hdr.scl_inter == 0));

  if (!scaledata) {
    // voxel values are unscaled -- we use the file's data type
    if (hdr.datatype == DT_UNSIGNED_CHAR) {
      fs_type = MRI_UCHAR;
      bytes_per_voxel = 1;
    }
    else if (hdr.datatype == DT_SIGNED_SHORT) {
      fs_type = MRI_SHORT;
      bytes_per_voxel = 2;
    }
    else if (hdr.datatype == DT_UINT16) {
      // This will not always work ...
      printf("INFO: This is an unsigined short.\n");    // I'll try to read it, but\n");
      //printf("      it might not work if there are values over 32k\n");
      fs_type = MRI_USHRT;
      bytes_per_voxel = 2;
    }
    else if (hdr.datatype == DT_SIGNED_INT) {
      fs_type = MRI_INT;
      bytes_per_voxel = 4;
    }
    else if (hdr.datatype == DT_FLOAT) {
      fs_type = MRI_FLOAT;
      bytes_per_voxel = 4;
    }
    else if (hdr.datatype == DT_DOUBLE) {
      fs_type = MRI_FLOAT;
      bytes_per_voxel = 8;
      printf("niiRead(): detected input as 64 bit double, reading in as 32 bit float\n");
    }
    else {
      ErrorReturn(
          NULL,
          (ERROR_UNSUPPORTED, "niiRead(): unsupported datatype %d (with scl_slope = 0) in %s", hdr.datatype, fname));
    }
  }
  else {
    // we must scale the voxel values
    if (hdr.datatype != DT_UNSIGNED_CHAR && hdr.datatype != DT_SIGNED_SHORT && hdr.datatype != DT_SIGNED_INT &&
        hdr.datatype != DT_FLOAT && hdr.datatype != DT_DOUBLE && hdr.datatype != DT_INT8 && hdr.datatype != DT_UINT16 &&
        hdr.datatype != DT_UINT32) {
      ErrorReturn(
          NULL,
          (ERROR_UNSUPPORTED, "niiRead(): unsupported datatype %d (with scl_slope != 0) in %s", hdr.datatype, fname));
    }
    fs_type = MRI_FLOAT;
    bytes_per_voxel = 0; /* set below -- avoid the compiler warning */
  }

  // Check whether dim[1] is less than 0. This can happen when FreeSurfer
  // writes a nifti volume where the number of columns is more than shortmax.
  // In this case, dim[1] is set to -1, and glmin is set to the number of cols.
  if (hdr.dim[1] > 0)
    ncols = hdr.dim[1];
  else
    ncols = hdr.glmin;

  if (ncols * hdr.dim[2] * hdr.dim[3] == 163842) IsIco7 = 1;

  if (read_volume)
    mri = MRIallocSequence(ncols, hdr.dim[2], hdr.dim[3], fs_type, nslices);
  else {
    if (!IsIco7)
      mri = MRIallocHeader(ncols, hdr.dim[2], hdr.dim[3], fs_type, nslices);
    else
      mri = MRIallocHeader(163842, 1, 1, fs_type, nslices);
    mri->nframes = nslices;
  }
  if (mri == NULL) return (NULL);

  mri->xsize = hdr.pixdim[1];
  mri->ysize = hdr.pixdim[2];
  mri->zsize = hdr.pixdim[3];
  mri->tr = hdr.pixdim[4];

  // Set the vox2ras matrix
  if (hdr.sform_code != 0) {
    // First, use the sform, if that is ok. Using the sform
    // first makes it more compatible with FSL.
    // fprintf(stderr,"INFO: using NIfTI-1 sform \n");
    if (niftiSformToMri(mri, &hdr) != NO_ERROR) {
      MRIfree(&mri);
      return (NULL);
    }
    mri->ras_good_flag = 1;
  }
  else if (hdr.qform_code != 0) {
    // Then, try the qform, if that is ok
    fprintf(stderr, "INFO: using NIfTI-1 qform \n");
    if (niftiQformToMri(mri, &hdr) != NO_ERROR) {
      MRIfree(&mri);
      return (NULL);
    }
    mri->ras_good_flag = 1;
  }
  else {
    // Should probably just die here.
    fprintf(stderr, "WARNING: neither NIfTI-1 qform or sform are valid\n");
    fprintf(stderr, "WARNING: your volume will probably be incorrectly oriented\n");
    mri->x_r = -1.0;
    mri->x_a = 0.0;
    mri->x_s = 0.0;
    mri->y_r = 0.0;
    mri->y_a = 1.0;
    mri->y_s = 0.0;
    mri->z_r = 0.0;
    mri->z_a = 0.0;
    mri->z_s = 1.0;
    mri->c_r = mri->xsize * mri->width / 2.0;
    mri->c_a = mri->ysize * mri->height / 2.0;
    mri->c_s = mri->zsize * mri->depth / 2.0;
    mri->ras_good_flag = 0;
  }

  mri->xsize = mri->xsize * space_units_factor;
  mri->ysize = mri->ysize * space_units_factor;
  mri->zsize = mri->zsize * space_units_factor;
  mri->c_r = mri->c_r * space_units_factor;
  mri->c_a = mri->c_a * space_units_factor;
  mri->c_s = mri->c_s * space_units_factor;
  mri->tr = mri->tr * time_units_factor;

  mri->fov = mri->xsize * mri->width;
  mri->xstart = -(mri->xsize * mri->width) / 2.0;
  mri->xend = (mri->xsize * mri->width) / 2.0;
  mri->ystart = -(mri->ysize * mri->height) / 2.0;
  mri->yend = (mri->ysize * mri->height) / 2.0;
  mri->zstart = -(mri->zsize * mri->depth) / 2.0;
  mri->zend = (mri->zsize * mri->depth) / 2.0;

  if (Gdiag_no > 0) {
    printf("nifti header ---------------------------------\n");
    niiPrintHdr(stdout, &hdr);
    printf("-----------------------------------------\n");
  }

  if (!read_volume) return (mri);

  fp = znzopen(fname, "r", use_compression);
  if (fp == NULL) {
    MRIfree(&mri);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error opening file %s", fname));
  }

  if (znzseek(fp, (long)(hdr.vox_offset), SEEK_SET) == -1) {
    znzclose(fp);
    MRIfree(&mri);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error finding voxel data in %s", fname));
  }

  if (!scaledata) {
    // no voxel value scaling needed
    void *buf;
    float *fbuf;
    double *dbuf;
    int nn;
    fbuf = (float *)calloc(mri->width, sizeof(float));
    dbuf = (double *)calloc(mri->width, sizeof(double));

    for (t = 0; t < mri->nframes; t++) {
      for (k = 0; k < mri->depth; k++) {
        for (j = 0; j < mri->height; j++) {
          buf = &MRIseq_vox(mri, 0, j, k, t);

          if (hdr.datatype != DT_DOUBLE)
            n_read = znzread(buf, bytes_per_voxel, mri->width, fp);
          else
            n_read = znzread(dbuf, bytes_per_voxel, mri->width, fp);
          if (n_read != mri->width) {
            printf("ERROR: Read %d, expected %d\n", n_read, mri->width);
            znzclose(fp);
            MRIfree(&mri);
            errno = 0;
            ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error reading from %s", fname));
          }

          if (swapped_flag) {
            if (bytes_per_voxel == 2) byteswapbufshort(buf, bytes_per_voxel * mri->width);
            if (bytes_per_voxel == 4) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
            if (bytes_per_voxel == 8) byteswapbuffloat(dbuf, bytes_per_voxel * mri->width);
          }
          if (hdr.datatype == DT_DOUBLE) {
            for (nn = 0; nn < mri->width; nn++) fbuf[nn] = (float)dbuf[nn];
            memcpy(buf, fbuf, 4 * mri->width);
          }
        }
        exec_progress_callback(k, mri->depth, t, mri->nframes);
      }
    }
    free(fbuf);
    free(dbuf);
  }
  else {
    // voxel value scaling needed
    if (hdr.datatype == DT_UNSIGNED_CHAR) {
      unsigned char *buf;
      bytes_per_voxel = 1;
      buf = (unsigned char *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = znzread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              znzclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error reading from %s", fname));
            }
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_SIGNED_SHORT) {
      short *buf;
      bytes_per_voxel = 2;
      buf = (short *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = znzread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              znzclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error reading from %s", fname));
            }
            if (swapped_flag) byteswapbufshort(buf, bytes_per_voxel * mri->width);
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_SIGNED_INT) {
      int *buf;
      bytes_per_voxel = 4;
      buf = (int *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = znzread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              znzclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error reading from %s", fname));
            }
            if (swapped_flag) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_FLOAT) {
      float *buf;
      bytes_per_voxel = 4;
      buf = (float *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = znzread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              znzclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error reading from %s", fname));
            }
            if (swapped_flag) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_DOUBLE) {
      double *buf;
      unsigned char *cbuf, ccbuf[8];
      bytes_per_voxel = 8;
      buf = (double *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = znzread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              znzclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error reading from %s", fname));
            }
            if (swapped_flag) {
              for (i = 0; i < mri->width; i++) {
                cbuf = (unsigned char *)&buf[i];
                memmove(ccbuf, cbuf, 8);
                cbuf[0] = ccbuf[7];
                cbuf[1] = ccbuf[6];
                cbuf[2] = ccbuf[5];
                cbuf[3] = ccbuf[4];
                cbuf[4] = ccbuf[3];
                cbuf[5] = ccbuf[2];
                cbuf[6] = ccbuf[1];
                cbuf[7] = ccbuf[0];
              }
            }
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_INT8) {
      char *buf;
      bytes_per_voxel = 1;
      buf = (char *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = znzread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              znzclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error reading from %s", fname));
            }
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_UINT16) {
      unsigned short *buf;
      bytes_per_voxel = 2;
      buf = (unsigned short *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = znzread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              znzclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error reading from %s", fname));
            }
            if (swapped_flag) byteswapbufshort(buf, bytes_per_voxel * mri->width);
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }

    if (hdr.datatype == DT_UINT32) {
      unsigned int *buf;
      bytes_per_voxel = 4;
      buf = (unsigned int *)malloc(mri->width * bytes_per_voxel);
      for (t = 0; t < mri->nframes; t++)
        for (k = 0; k < mri->depth; k++) {
          for (j = 0; j < mri->height; j++) {
            n_read = znzread(buf, bytes_per_voxel, mri->width, fp);
            if (n_read != mri->width) {
              free(buf);
              znzclose(fp);
              MRIfree(&mri);
              errno = 0;
              ErrorReturn(NULL, (ERROR_BADFILE, "niiRead(): error reading from %s", fname));
            }
            if (swapped_flag) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
            for (i = 0; i < mri->width; i++)
              MRIFseq_vox(mri, i, j, k, t) = hdr.scl_slope * (float)(buf[i]) + hdr.scl_inter;
          }
          exec_progress_callback(k, mri->depth, t, mri->nframes);
        }
      free(buf);
    }
  }
  znzclose(fp);

  // Check for ico7 surface
  if (IsIco7) {
    //   printf("niiRead: reshaping\n");
    mritmp = mri_reshape(mri, 163842, 1, 1, mri->nframes);
    MRIfree(&mri);
    mri = mritmp;
  }

  return (mri);

} /* end niiRead() */



/*------------------------------------------------------------------
  niiReadFromMriFsStruct() - note: there are also nifti1Read() and niiRead(). Make sure to
  edit all. Automatically detects whether an input is Ico7
  and reshapes.
  This function is used with DICOMRead3().
  DICOMRead3() read/parse dicom files using dcm2niix_fswrapper.
  MRIFSSTRUCT holds the parsing output from dcm2niix, which contains nifti header, 
  image data, acqusition parameters, and bvecs. See nii_dicom_batch.h for details.
  -----------------------------------------------------------------*/
static MRI *niiReadFromMriFsStruct(MRIFSSTRUCT *mrifsStruct)
{
  MRI *mri, *mritmp;
  int nslices;
  int fs_type;
  float time_units_factor, space_units_factor;
  int swapped_flag;
  int bytes_per_voxel, time_units, space_units;
  int ncols, IsIco7 = 0;

  struct nifti_1_header *hdr = &(mrifsStruct->hdr0);
  const unsigned char *img = mrifsStruct->imgM;
  struct TDICOMdata *tdicomData = &mrifsStruct->tdicomData;

  swapped_flag = FALSE;
  if (hdr->dim[0] < 1 || hdr->dim[0] > 7) {
    swapped_flag = TRUE;
    swap_nifti_1_header(hdr);
    if (hdr->dim[0] < 1 || hdr->dim[0] > 7) {
      ErrorReturn(NULL, (ERROR_BADFILE, "niiReadFromMriFsStruct(): bad number of dimensions (%hd)", hdr->dim[0]));
    }
  }

  if (memcmp(hdr->magic, NII_MAGIC, 4) != 0) {
    ErrorReturn(NULL, (ERROR_BADFILE, "niiReadFromMriFsStruct(): bad magic number"));
  }

  //  if (hdr.dim[0] != 2 && hdr.dim[0] != 3 && hdr.dim[0] != 4){
  if (hdr->dim[0] < 1 || hdr->dim[0] > 5) {
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "niiReadFromMriFsStruct(): %hd dimensions; unsupported", hdr->dim[0]));
  }

  if (hdr->datatype == DT_NONE || hdr->datatype == DT_UNKNOWN) {
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "niiReadFromMriFsStruct(): unknown or no data type; bailing out"));
  }
  if (hdr->dim[4] == 0) {
    printf("WARNING: hdr->dim[4] = 0 (nframes), setting to 1\n");
    hdr->dim[4] = 1;
  }

  space_units = XYZT_TO_SPACE(hdr->xyzt_units);
  if (space_units == NIFTI_UNITS_METER)
    space_units_factor = 1000.0;
  else if (space_units == NIFTI_UNITS_MM)
    space_units_factor = 1.0;
  else if (space_units == NIFTI_UNITS_MICRON)
    space_units_factor = 0.001;
  else if (space_units == NIFTI_UNITS_UNKNOWN) {
    space_units_factor = 1.0;
  }
  else
    ErrorReturn(NULL, (ERROR_BADFILE, "niiReadFromMriFsStruct(): unknown space units %d", space_units));

  time_units = XYZT_TO_TIME(hdr->xyzt_units);
  if (time_units == NIFTI_UNITS_SEC)
    time_units_factor = 1000.0;
  else if (time_units == NIFTI_UNITS_MSEC)
    time_units_factor = 1.0;
  else if (time_units == NIFTI_UNITS_USEC)
    time_units_factor = 0.001;
  else {
    // This can be a tricky situation because the time units may not mean
    // anything even with multiple frames (eg, BO mag and phase).
    if (hdr->dim[4] > 1) printf("WARNING: niiRead(): unknown time units %d\n", time_units);
    time_units_factor = 0;
  }
  // printf("hdr.xyzt_units = %d, time_units = %d, %g, %g\n",
  // hdr.xyzt_units,time_units,hdr.pixdim[4],time_units_factor);

  if (hdr->slice_code != 0 && DIM_INFO_TO_SLICE_DIM(hdr->dim_info) != 0 && hdr->slice_duration > 0.0) {
    if (hdr->slice_code != NIFTI_SLICE_SEQ_INC && hdr->slice_code != NIFTI_SLICE_SEQ_DEC &&
        hdr->slice_code != NIFTI_SLICE_ALT_INC && hdr->slice_code != NIFTI_SLICE_ALT_DEC &&
        hdr->slice_code != NIFTI_SLICE_ALT_INC2 && hdr->slice_code != NIFTI_SLICE_ALT_DEC2) {
      ErrorReturn(
          NULL, (ERROR_UNSUPPORTED, "niiReadFromMriFsStruct(): unsupported slice timing pattern %d", hdr->slice_code));
    }
  }

  if (hdr->dim[0] < 4)
    nslices = 1;
  else
    nslices = hdr->dim[4];

  // put extra dims in frames
  if (hdr->dim[0] > 4 && hdr->dim[5] > 0) nslices *= hdr->dim[5];
  if (Gdiag_no > 0)
    printf("niiReadFromMriFsStruct(): hdr->dim %d %d %d %d %d %d\n",
           hdr->dim[0],
           hdr->dim[1],
           hdr->dim[2],
           hdr->dim[3],
           hdr->dim[4],
           hdr->dim[5]);

  // check whether data needs to be scaled
  bool scaledata = (hdr->scl_slope != 0) && !((hdr->scl_slope == 1) && (hdr->scl_inter == 0));

  // example values from tag (0028, 1052) and (0028, 1053):
  // (0028,1052) DS [-4096]                                  #   6, 1 RescaleIntercept
  // (0028,1053) DS [2]                                      #   2, 1 RescaleSlope
  if (scaledata)
  {
    printf("\nRescale Slope (scl_slope) = %f, Rescale Intercept (scl_inter) = %f\n", hdr->scl_slope, hdr->scl_inter);
    printf("Voxel values need to be scaled!\n");
    if (getenv("NII_RESCALE_OFF"))
    {
      printf("Environment Variable NII_RESCALE_OFF set. Turn off voxel values scaling ...\n\n");
      scaledata = false;
    }
  }

  if (!scaledata) {
    // voxel values are unscaled -- we use the file's data type
    if (hdr->datatype == DT_UNSIGNED_CHAR) {
      fs_type = MRI_UCHAR;
      bytes_per_voxel = 1;
    }
    else if (hdr->datatype == DT_SIGNED_SHORT) {
      fs_type = MRI_SHORT;
      bytes_per_voxel = 2;
    }
    else if (hdr->datatype == DT_UINT16) {
      printf("INFO: This is an unsigned short.\n");
      fs_type = MRI_USHRT;
      bytes_per_voxel = 2;
    }
    else if (hdr->datatype == DT_SIGNED_INT) {
      fs_type = MRI_INT;
      bytes_per_voxel = 4;
    }
    else if (hdr->datatype == DT_FLOAT) {
      fs_type = MRI_FLOAT;
      bytes_per_voxel = 4;
    }
    else if (hdr->datatype == DT_DOUBLE) {
      fs_type = MRI_FLOAT;
      bytes_per_voxel = 8;
      printf("niiReadFromMriFsStruct(): detected input as 64 bit double, reading in as 32 bit float\n");
    }
#if 0    // MRI_RBG not support in mghWrite
    else if (hdr->datatype == DT_RGB) {
      fs_type = MRI_UCHAR;    //MRI_RGB;
      bytes_per_voxel = 3;
      printf("niiReadFromMriFsStruct(): DT_RGB, MRI_RGB\n");
    }
#endif
    else {
      ErrorReturn(
          NULL,
          (ERROR_UNSUPPORTED, "niiReadFromMriFsStruct(): unsupported datatype %d (with scl_slope = 0)", hdr->datatype)); 
    }
  }
  else {
    // we must scale the voxel values
    if (hdr->datatype != DT_UNSIGNED_CHAR && hdr->datatype != DT_SIGNED_SHORT && hdr->datatype != DT_SIGNED_INT &&
        hdr->datatype != DT_FLOAT && hdr->datatype != DT_DOUBLE && hdr->datatype != DT_INT8 && hdr->datatype != DT_UINT16 &&
        hdr->datatype != DT_UINT32) {
      ErrorReturn(
          NULL,
          (ERROR_UNSUPPORTED, "niiReadFromMriFsStruct(): unsupported datatype %d (with scl_slope != 0)", hdr->datatype));
    }
    fs_type = MRI_FLOAT;
    bytes_per_voxel = 0; /* set below -- avoid the compiler warning */
  }

  // Check whether dim[1] is less than 0. This can happen when FreeSurfer
  // writes a nifti volume where the number of columns is more than shortmax.
  // In this case, dim[1] is set to -1, and glmin is set to the number of cols.
  if (hdr->dim[1] > 0)
    ncols = hdr->dim[1];
  else
    ncols = hdr->glmin;

  if (ncols * hdr->dim[2] * hdr->dim[3] == 163842) IsIco7 = 1;

  // process image
  int read_volume = 1;  // ????should it be 1 or 0????
  if (read_volume)
    mri = MRIallocSequence(ncols, hdr->dim[2], hdr->dim[3], fs_type, nslices);
  else {
    if (!IsIco7)
      mri = MRIallocHeader(ncols, hdr->dim[2], hdr->dim[3], fs_type, nslices);
    else
      mri = MRIallocHeader(163842, 1, 1, fs_type, nslices);
    mri->nframes = nslices;
  }
  if (mri == NULL) return (NULL);

  mri->xsize = hdr->pixdim[1];
  mri->ysize = hdr->pixdim[2];
  mri->zsize = hdr->pixdim[3];
  mri->tr = hdr->pixdim[4];

  // set MRI values from struct TDICOMdata
  mri->te = tdicomData->TE;
  mri->ti = tdicomData->TI;
  mri->flip_angle = M_PI * tdicomData->flipAngle / 180.0; 
  mri->FieldStrength = tdicomData->fieldStrength;
  if (tdicomData->phaseEncodingRC == 'R')
    mri->pedir = strcpyalloc("ROW");
  else if (tdicomData->phaseEncodingRC == 'C')
    mri->pedir = strcpyalloc("COL");

  /* start setting of bvals and bvecs */
  // set bvals and bvecs
  struct TDTI* tdti = mrifsStruct->tdti;
  int numDti = mrifsStruct->numDti;
  if (numDti > 0)
  {
    int IsDWI = 1, IsPhilipsDWI = 0;
    if (tdicomData->manufacturer == kMANUFACTURER_PHILIPS)
      IsPhilipsDWI = 1;
    printf("IsDWI = %d, IsPhilipsDWI = %d\n", IsDWI, IsPhilipsDWI);

    mri->bvals = MatrixAlloc(numDti, 1, MATRIX_REAL);
    mri->bvecs = MatrixAlloc(numDti, 3, MATRIX_REAL);
    if (DIAG_VERBOSE_ON)
    {
      for (int i = 0; i < numDti; i++)
	printf("bvals: %f\n", tdti[i].V[0]);  //mri->bvals->rptr[i][1] = tdti[i].V[v];

      for (int i = 0; i < numDti; i++)
	printf("bvecs: %16.14f %16.14f %16.14f\n", tdti[i].V[1], tdti[i].V[2], tdti[i].V[3]);
    }

    // mri->bvals, mri->bvecs start at index 1; tdti[].V[] starts at index 0
    for (int i = 0; i < numDti; i++)
    {
      mri->bvals->rptr[i+1][1] = tdti[i].V[0];
      mri->bvecs->rptr[i+1][1] = tdti[i].V[1];
      mri->bvecs->rptr[i+1][2] = tdti[i].V[2];
      mri->bvecs->rptr[i+1][3] = tdti[i].V[3];

#if 0
      for (int v = 0; v < 4; v++)
      {
	if (v == 0)
	  mri->bvals->rptr[i+1][1] = tdti[i].V[v];
	else  // v = 1, 2, 3
	  mri->bvecs->rptr[i+1][v] = tdti[i].V[v];
      }
#endif
    }

    // assign bvec_space
    mri->bvec_space = BVEC_SPACE_VOXEL;
    if (IsPhilipsDWI)
      mri->bvec_space = BVEC_SPACE_SCANNER;
  }
  /* end setting of bvals and bvecs */

  // Set the vox2ras matrix
  if (hdr->sform_code != 0) {
    // First, use the sform, if that is ok. Using the sform
    // first makes it more compatible with FSL.
    // fprintf(stderr,"INFO: using NIfTI-1 sform \n");
    if (niftiSformToMri(mri, hdr) != NO_ERROR) {
      MRIfree(&mri);
      return (NULL);
    }
    mri->ras_good_flag = 1;
  }
  else if (hdr->qform_code != 0) {
    // Then, try the qform, if that is ok
    fprintf(stderr, "INFO: using NIfTI-1 qform \n");
    if (niftiQformToMri(mri, hdr) != NO_ERROR) {
      MRIfree(&mri);
      return (NULL);
    }
    mri->ras_good_flag = 1;
  }
  else {
    // Should probably just die here.
    fprintf(stderr, "WARNING: neither NIfTI-1 qform or sform are valid\n");
    fprintf(stderr, "WARNING: your volume will probably be incorrectly oriented\n");
    mri->x_r = -1.0;
    mri->x_a = 0.0;
    mri->x_s = 0.0;
    mri->y_r = 0.0;
    mri->y_a = 1.0;
    mri->y_s = 0.0;
    mri->z_r = 0.0;
    mri->z_a = 0.0;
    mri->z_s = 1.0;
    mri->c_r = mri->xsize * mri->width / 2.0;
    mri->c_a = mri->ysize * mri->height / 2.0;
    mri->c_s = mri->zsize * mri->depth / 2.0;
    mri->ras_good_flag = 0;
  }

  mri->xsize = mri->xsize * space_units_factor;
  mri->ysize = mri->ysize * space_units_factor;
  mri->zsize = mri->zsize * space_units_factor;
  mri->c_r = mri->c_r * space_units_factor;
  mri->c_a = mri->c_a * space_units_factor;
  mri->c_s = mri->c_s * space_units_factor;
  mri->tr = mri->tr * time_units_factor;

  mri->fov = mri->xsize * mri->width;
  mri->xstart = -(mri->xsize * mri->width) / 2.0;
  mri->xend = (mri->xsize * mri->width) / 2.0;
  mri->ystart = -(mri->ysize * mri->height) / 2.0;
  mri->yend = (mri->ysize * mri->height) / 2.0;
  mri->zstart = -(mri->zsize * mri->depth) / 2.0;
  mri->zend = (mri->zsize * mri->depth) / 2.0;

  if (Gdiag_no > 0) {
    printf("nifti header ---------------------------------\n");
    niiPrintHdr(stdout, hdr);
    printf("-----------------------------------------\n");
  }

  if (!read_volume) return (mri);

  int j, k, t;
  if (!scaledata) {
    // no voxel value scaling needed
    void *buf;
    float *fbuf = (float *)calloc(mri->width, sizeof(float));
    double *dbuf = (double *)calloc(mri->width, sizeof(double));
    int nn;

    int bytesRead = 0;
    for (t = 0; t < mri->nframes; t++) {
      for (k = 0; k < mri->depth; k++) {
        for (j = 0; j < mri->height; j++) {
          buf = &MRIseq_vox(mri, 0, j, k, t);

          if (hdr->datatype != DT_DOUBLE)
	  {
	    memcpy(buf, img+bytesRead, bytes_per_voxel*mri->width);
	  }
          else
	  {
            memcpy(dbuf, img+bytesRead, bytes_per_voxel*mri->width);
          }
          bytesRead += bytes_per_voxel*mri->width;

          if (swapped_flag) {
            if (bytes_per_voxel == 2) byteswapbufshort(buf, bytes_per_voxel * mri->width);
            if (bytes_per_voxel == 4) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
            if (bytes_per_voxel == 8) byteswapbuffloat(dbuf, bytes_per_voxel * mri->width);
          }
          if (hdr->datatype == DT_DOUBLE) {
            for (nn = 0; nn < mri->width; nn++) fbuf[nn] = (float)dbuf[nn];
            memcpy(buf, fbuf, 4 * mri->width);
          }
        }  // height
        exec_progress_callback(k, mri->depth, t, mri->nframes);
      }  // depth
    }  // nframes
    free(fbuf);
    free(dbuf);
  } 
  else {
      int bytesRead = 0;

      printf("Voxel values scaling ... (Set environment variable 'NII_RESCALE_OFF' to turn off scaling)\n\n");

      /* these are set earlier:
       * fs_type = MRI_FLOAT;
       * bytes_per_voxel = 0;
       */
      if (hdr->datatype == DT_UNSIGNED_CHAR)
	bytes_per_voxel = 1;
      else if (hdr->datatype == DT_SIGNED_SHORT)
	bytes_per_voxel = 2;
      else if (hdr->datatype == DT_SIGNED_INT)
	bytes_per_voxel = 4;
      else if (hdr->datatype == DT_FLOAT)
	bytes_per_voxel = 4;
      else if (hdr->datatype == DT_DOUBLE)
	bytes_per_voxel = 8;
      else if (hdr->datatype == DT_UINT8)
	bytes_per_voxel = 1;
      else if (hdr->datatype == DT_UINT16)
	bytes_per_voxel = 2;
      else if (hdr->datatype == DT_UINT32)
	bytes_per_voxel = 4;

    // voxel value scaling needed
    unsigned char *buf = (unsigned char *)malloc(mri->width * bytes_per_voxel);
    int nn;

    for (t = 0; t < mri->nframes; t++) {
      for (k = 0; k < mri->depth; k++) {
	for (j = 0; j < mri->height; j++) {
          memcpy(buf, img+bytesRead, bytes_per_voxel*mri->width);
          bytesRead += bytes_per_voxel*mri->width;

          if (swapped_flag) {
            if (bytes_per_voxel == 2) byteswapbufshort(buf, bytes_per_voxel * mri->width);
            if (bytes_per_voxel == 4) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
            if (bytes_per_voxel == 8) byteswapbuffloat(buf, bytes_per_voxel * mri->width);
          }

          // additional logic to scale voxel values
          for (nn = 0; nn < mri->width; nn++)
	  {
            if (hdr->datatype == DT_UNSIGNED_CHAR) {
              MRIFseq_vox(mri, nn, j, k, t) = hdr->scl_slope * ((unsigned char*)buf)[nn] + hdr->scl_inter;
            }
            else if (hdr->datatype == DT_SIGNED_SHORT) {
              MRIFseq_vox(mri, nn, j, k, t) = hdr->scl_slope * ((short*)buf)[nn] + hdr->scl_inter;
            }
            else if (hdr->datatype == DT_UINT16) {
              MRIFseq_vox(mri, nn, j, k, t) = hdr->scl_slope * ((unsigned short*)buf)[nn] + hdr->scl_inter;
            }
            else if (hdr->datatype == DT_SIGNED_INT) {
              MRIFseq_vox(mri, nn, j, k, t) = hdr->scl_slope * ((int*)buf)[nn] + hdr->scl_inter;
            }
            else if (hdr->datatype == DT_FLOAT) {
              MRIFseq_vox(mri, nn, j, k, t) = hdr->scl_slope * ((float*)buf)[nn] + hdr->scl_inter;
            }
            else if (hdr->datatype == DT_DOUBLE) {
              MRIFseq_vox(mri, nn, j, k, t) = hdr->scl_slope * ((double*)buf)[nn] + hdr->scl_inter;
            }
            else if (hdr->datatype == DT_UINT32) {
              MRIFseq_vox(mri, nn, j, k, t) = hdr->scl_slope * ((unsigned int*)buf)[nn] + hdr->scl_inter;
            }
            else if (hdr->datatype == DT_INT8) {
              MRIFseq_vox(mri, nn, j, k, t) = hdr->scl_slope * ((char*)buf)[nn] + hdr->scl_inter;
            }
          }
        } // height
        exec_progress_callback(k, mri->depth, t, mri->nframes);
      } // depth
    } //nframes

    //free(fbuf);
    //free(dbuf);

#if 0 // original implementation
    int i;

    // voxel value scaling needed
    unsigned char *buf = (unsigned char *)malloc(mri->width * bytes_per_voxel);
    for (t = 0; t < mri->nframes; t++) {
      for (k = 0; k < mri->depth; k++) {
	for (j = 0; j < mri->height; j++) {
	  memcpy(buf, img+bytesRead, bytes_per_voxel*mri->width);
          bytesRead += bytes_per_voxel*mri->width;

	  if (swapped_flag) {
	    if (hdr->datatype == DT_SIGNED_SHORT)
	      byteswapbufshort(buf, bytes_per_voxel * mri->width);
	    else if (hdr->datatype == DT_SIGNED_INT)
	      byteswapbuffloat(buf, bytes_per_voxel * mri->width);
	    else if (hdr->datatype == DT_FLOAT)
	      byteswapbuffloat(buf, bytes_per_voxel * mri->width);
	    else if (hdr->datatype == DT_DOUBLE) {
	      unsigned char *cbuf, ccbuf[8];
	      for (i = 0; i < mri->width; i++) {
                cbuf = (unsigned char *)&buf[i];
                memmove(ccbuf, cbuf, 8);
                cbuf[0] = ccbuf[7];
                cbuf[1] = ccbuf[6];
                cbuf[2] = ccbuf[5];
                cbuf[3] = ccbuf[4];
                cbuf[4] = ccbuf[3];
                cbuf[5] = ccbuf[2];
                cbuf[6] = ccbuf[1];
                cbuf[7] = ccbuf[0];
	      }
	    }
	    else if (hdr->datatype == DT_UINT16)
	      byteswapbufshort(buf, bytes_per_voxel * mri->width);
	    else if (hdr->datatype == DT_UINT32)
	      byteswapbuffloat(buf, bytes_per_voxel * mri->width);
	  }

          for (i = 0; i < mri->width; i++)
            MRIFseq_vox(mri, i, j, k, t) = hdr->scl_slope * (float)(buf[i]) + hdr->scl_inter;
        }
        exec_progress_callback(k, mri->depth, t, mri->nframes);
      }
    }
    free(buf);
#endif
  }

  // Check for ico7 surface
  if (IsIco7) {
    //   printf("niiRead: reshaping\n");
    mritmp = mri_reshape(mri, 163842, 1, 1, mri->nframes);
    MRIfree(&mri);
    mri = mritmp;
  }

  return (mri);

} /* end niiReadFromMriFsStruct() */

/*------------------------------------------------------------------
  niiWrite() - note: there is also an nifti1Write(). Make sure to
  edit both. Automatically detects whether an input is Ico7
  and reshapes.
  -----------------------------------------------------------------*/
static int niiWrite(MRI *mri0, const char *fname)
{
  znzFile fp;
  int j, k, t;
  BUFTYPE *buf;
  char *chbuf;
  struct nifti_1_header hdr;
  int error, shortmax, use_compression, fnamelen, nfill;
  MRI *mri = NULL;
  int FreeMRI = 0;

  // printf("In niiWrite()\n");

  use_compression = 0;
  fnamelen = strlen(fname);
  if (fname[fnamelen - 1] == 'z') use_compression = 1;
  if (Gdiag_no > 0) printf("niiWrite: use_compression = %d\n", use_compression);

  // Check for ico7 surface
  if (mri0->width == 163842 && mri0->height == 1 && mri0->depth == 1) {
    // printf("niiWrite: reshaping\n");
    mri = mri_reshape(mri0, 27307, 1, 6, mri0->nframes);
    FreeMRI = 1;
  }
  else
    mri = mri0;

  shortmax = (int)(pow(2.0, 15.0));
  if (0 && mri->width > shortmax) {
    printf("NIFTI FORMAT WARNING: ncols %d in input exceeds %d.\n", mri->width, shortmax);
    printf("So I'm going to put the true ncols in glmin and set dim[1]=-1.\n");
    printf("This should be ok within FreeSurfer, but you will not be\n");
    printf("able to use this volume with other software.\n");
    // This usually happens with surfaces.
    // exit(1);
  }
  if (mri->height > shortmax) {
    printf("NIFTI FORMAT ERROR: nrows %d in volume exceeds %d\n", mri->height, shortmax);
    exit(1);
  }
  if (mri->depth > shortmax) {
    printf("NIFTI FORMAT ERROR: nslices %d in volume exceeds %d\n", mri->depth, shortmax);
    exit(1);
  }
  if (mri->nframes > shortmax) {
    printf("NIFTI FORMAT ERROR:  nframes %d in volume exceeds %d\n", mri->nframes, shortmax);
    exit(1);
  }

  memset(&hdr, 0x00, sizeof(hdr));

  hdr.sizeof_hdr = 348;
  hdr.dim_info = 0;

  for (t = 0; t < 8; t++) {
    hdr.dim[t] = 1;
    hdr.pixdim[t] = 1;
  }  // for afni
  if (mri->nframes == 1)
    hdr.dim[0] = 3;
  else
    hdr.dim[0] = 4;

  if (mri->width < shortmax)
    hdr.dim[1] = mri->width;
  else {
    // number of columns too big, put in glmin
    hdr.dim[1] = -1;
    hdr.glmin = mri->width;
  }
  hdr.dim[2] = mri->height;
  hdr.dim[3] = mri->depth;
  hdr.dim[4] = mri->nframes;
  hdr.pixdim[1] = mri->xsize;
  hdr.pixdim[2] = mri->ysize;
  hdr.pixdim[3] = mri->zsize;
  hdr.pixdim[4] = mri->tr / 1000.0;  // see also xyzt_units

  // printf("In niiWrite(): after init\n");

  if (mri->type == MRI_UCHAR) {
    hdr.datatype = DT_UNSIGNED_CHAR;
    hdr.bitpix = 8;
  }
  else if (mri->type == MRI_INT) {
    hdr.datatype = DT_SIGNED_INT;
    hdr.bitpix = 32;
  }
  else if (mri->type == MRI_LONG) {
    hdr.datatype = DT_SIGNED_INT;
    hdr.bitpix = 32;
  }
  else if (mri->type == MRI_FLOAT) {
    hdr.datatype = DT_FLOAT;
    hdr.bitpix = 32;
  }
  else if (mri->type == MRI_SHORT) {
    hdr.datatype = DT_SIGNED_SHORT;
    hdr.bitpix = 16;
  }
  else if (mri->type == MRI_USHRT) {
    hdr.datatype = DT_UINT16;
    hdr.bitpix = 16;
  }
  else if (mri->type == MRI_BITMAP) {
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "niiWrite(): data type MRI_BITMAP unsupported"));
  }
  else if (mri->type == MRI_TENSOR) {
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "niiWrite(): data type MRI_TENSOR unsupported"));
  }
  else {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "niiWrite(): unknown data type %d", mri->type));
  }

  hdr.intent_code = NIFTI_INTENT_NONE;
  hdr.intent_name[0] = '\0';
  hdr.vox_offset = 352;  // 352 is the min, dont use sizeof(hdr); See below
  hdr.scl_slope = 0.0;
  hdr.slice_code = 0;
  hdr.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;
  hdr.cal_max = 0.0;
  hdr.cal_min = 0.0;
  hdr.toffset = 0;
  sprintf(hdr.descrip, "FreeSurfer %s", __DATE__);

  /* set the nifti header qform values */
  error = mriToNiftiQform(mri, &hdr);
  if (error != NO_ERROR) return (error);

  // printf("In niiWrite():before sform\n");

  /* set the nifti header sform values */
  // This just copies the vox2ras into the sform
  mriToNiftiSform(mri, &hdr);

  memmove(hdr.magic, NII_MAGIC, 4);

  fp = znzopen(fname, "w", use_compression);
  if (fp == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "niiWrite(): error opening file %s", fname));
  }

  // White the header
  if (znzwrite(&hdr, sizeof(hdr), 1, fp) != 1) {
    znzclose(fp);
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "niiWrite(): error writing header to %s", fname));
  }

  // Fill in space to the voxel offset
  nfill = (int)hdr.vox_offset - sizeof(hdr);
  chbuf = (char *)calloc(nfill, sizeof(char));
  if ((int)znzwrite(chbuf, sizeof(char), nfill, fp) != nfill) {
    znzclose(fp);
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "niiWrite(): error writing data to %s", fname));
  }
  free(chbuf);

  // printf("In niiWrite():before dumping: %d, %d, %d, %d\n", mri->nframes,mri->depth,mri->width,mri->height );
  // Now dump the pixel data
  for (t = 0; t < mri->nframes; t++)
    for (k = 0; k < mri->depth; k++) {
      for (j = 0; j < mri->height; j++) {
        // printf("%d,%d,%d\n",t, k, j);
        buf = &MRIseq_vox(mri, 0, j, k, t);
        if ((int)znzwrite(buf, hdr.bitpix / 8, mri->width, fp) != mri->width) {
          znzclose(fp);
          errno = 0;
          ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "niiWrite(): error writing data to %s", fname));
        }
      }
      exec_progress_callback(k, mri->depth, t, mri->nframes);
    }

  znzclose(fp);

  if (FreeMRI) MRIfree(&mri);

  return (NO_ERROR);

} /* end niiWrite() */


/*------------------------------------------------------------------
  mriNrrdRead
  -----------------------------------------------------------------*/
static MRI *mriNrrdRead(const char *fname, int read_volume)
{
  if (!nrrdSanity()) {
    fprintf(stderr, "\n");
    fprintf(stderr, "!!! nrrd sanity check FAILED: fix and re-compile\n");
    char *err = biffGet(NRRD);
    fprintf(stderr, "%s\n", err);
    free(err);
    return NULL;
  }

  /* create a nrrd; at this point this is just an empty container */
  Nrrd *nin = nrrdNew();
  NrrdIoState *nio = nrrdIoStateNew();

  /* read in the nrrd from file */
  if (nrrdLoad(nin, fname, nio)) {
    char *err = biffGetDone(NRRD);
    fprintf(stderr, "mriNrrdRead: trouble reading \"%s\":\n%s", fname, err);
    free(err);
    return NULL;
  }

  // if it has more than 3 dimensions, then maybe its diffusion data
  if (((nin->dim != 3) && (nin->dim != 4)) || (nin->spaceDim != 3)) {
    /* say something about the array */
    printf("mriNrrdRead: \"%s\" is a %d-dimensional nrrd of type %d (%s)\n",
           fname,
           nin->dim,
           nin->type,
           airEnumStr(nrrdType, nin->type));
    printf("mriNrrdRead: the array contains %d elements, %d bytes in size\n",
           (int)nrrdElementNumber(nin),
           (int)nrrdElementSize(nin));
    if (nin->content) printf("mriNrrdRead: content: %s\n", nin->content);

    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "Nrrd input of diffusion data not supported!"));
    return NULL;
  }

  /* print out the key/value pairs present */
  int kvn = nrrdKeyValueSize(nin);
  if (kvn) {
    int kvi;
    for (kvi = 0; kvi < kvn; kvi++) {
      char *val, *key;
      nrrdKeyValueIndex(nin, &key, &val, kvi);
      printf("mriNrrdRead: key:value %d = %s:%s\n", kvi, key, val);
      free(key);
      free(val);
      key = val = NULL;
    }
  }

  // Get the component type.
  int type = MRI_UCHAR;
  switch (nin->type) {
    case nrrdTypeChar:
    case nrrdTypeUChar:
      type = MRI_UCHAR;
      break;
    case nrrdTypeShort:
    case nrrdTypeUShort:
      type = MRI_SHORT;
      break;
    case nrrdTypeInt:
    case nrrdTypeUInt:
      type = MRI_INT;
      break;
    case nrrdTypeLLong:
    case nrrdTypeULLong:
      type = MRI_LONG;
      break;
    case nrrdTypeFloat:
      type = MRI_FLOAT;
      break;
    default:
      printf("mriNrrdRead: Unsupported type: %d (%s)\n", nin->type, airEnumStr(nrrdType, nin->type));
      return NULL;
  }

  // alloc mri struct with the correct dimensions.
  int width = nin->axis[0].size;
  int height = nin->axis[1].size;
  int depth = nin->axis[2].size;
  int nframes = 1;                                 // default nin->dim = 3
  if (nin->dim == 4) nframes = nin->axis[3].size;  // multiple frames found
  MRI *mri = MRIallocSequence(width, height, depth, type, nframes);
  if (NULL == mri) {
    printf("mriNrrdRead: Couldn't allocate MRI of size %d %d %d %d\n", width, height, depth, nframes);
    return NULL;
  }

  // Copy all the pixel data.
  int x, y, z, f;
  if (type == MRI_UCHAR) {
    for (f = 0; f < nframes; f++)
      for (z = 0; z < depth; z++)
        for (y = 0; y < height; y++)
          for (x = 0; x < width; x++) {
            int index = x + (y * width) + (z * width * height) + (f * width * height * depth);
            unsigned char *_uc = (unsigned char *)nin->data;
            MRIseq_vox(mri, x, y, z, f) = (BUFTYPE)_uc[index];
          }
  }
  else if (type == MRI_SHORT) {
    for (f = 0; f < nframes; f++)
      for (z = 0; z < depth; z++)
        for (y = 0; y < height; y++)
          for (x = 0; x < width; x++) {
            int index = x + (y * width) + (z * width * height) + (f * width * height * depth);
            short *_s = (short *)nin->data;
            MRISseq_vox(mri, x, y, z, f) = (short)_s[index];
          }
  }
  else if (type == MRI_INT) {
    for (f = 0; f < nframes; f++)
      for (z = 0; z < depth; z++)
        for (y = 0; y < height; y++)
          for (x = 0; x < width; x++) {
            int index = x + (y * width) + (z * width * height) + (f * width * height * depth);
            int *_i = (int *)nin->data;
            MRIIseq_vox(mri, x, y, z, f) = (int)_i[index];
          }
  }
  else if (type == MRI_LONG) {
    for (f = 0; f < nframes; f++)
      for (z = 0; z < depth; z++)
        for (y = 0; y < height; y++)
          for (x = 0; x < width; x++) {
            int index = x + (y * width) + (z * width * height) + (f * width * height * depth);
            long *_l = (long *)nin->data;
            MRILseq_vox(mri, x, y, z, f) = (long)_l[index];
          }
  }
  else if (type == MRI_FLOAT) {
    for (f = 0; f < nframes; f++)
      for (z = 0; z < depth; z++)
        for (y = 0; y < height; y++)
          for (x = 0; x < width; x++) {
            int index = x + (y * width) + (z * width * height) + (f * width * height * depth);
            float *_f = (float *)nin->data;
            MRIFseq_vox(mri, x, y, z, f) = (float)_f[index];
          }
  }
  else {
    printf("mriNrrdRead: Unsupported type=%d\n", type);
    return NULL;
  }

  // get and set the origin
  mri->c_r = (float)nin->spaceOrigin[0];
  mri->c_a = (float)nin->spaceOrigin[1];
  mri->c_s = (float)nin->spaceOrigin[2];
  mri->ras_good_flag = 1;

  // get and set the spacing
  mri->xsize = (float)nin->axis[0].spaceDirection[0];
  mri->ysize = (float)nin->axis[1].spaceDirection[1];
  mri->zsize = (float)nin->axis[2].spaceDirection[2];

  // set orientation string
  switch (nin->space) {
    case nrrdSpaceRightAnteriorSuperior: {
      MRIorientationStringToDircos(mri, "RAS");
      break;
    }
    case nrrdSpaceLeftAnteriorSuperior: {
      MRIorientationStringToDircos(mri, "LAS");
      break;
    }
    case nrrdSpaceLeftPosteriorSuperior: {
      MRIorientationStringToDircos(mri, "LPS");
      break;
    }
  }

  return mri;
}

/*------------------------------------------------------------------
  mriNrrdWrite
  -----------------------------------------------------------------*/
static int mriNrrdWrite(MRI *mri, const char *fname)
{
  ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "mriNrrdWrite(): Nrrd output not supported!"));
}

/*------------------------------------------------------------------
  itkMorphWrite()
  -----------------------------------------------------------------*/
static int itkMorphWrite(MRI *mri, const char *fname)
{
  znzFile fp;
  int j, k, t;
  BUFTYPE *buf;
  char *chbuf;
  struct nifti_1_header hdr;
  int error, shortmax, use_compression, fnamelen, nfill;

  use_compression = 0;
  fnamelen = strlen(fname);
  if (fname[fnamelen - 1] == 'z') use_compression = 1;
  if (Gdiag_no > 0) printf("itkMorphWrite: use_compression = %d\n", use_compression);

  shortmax = (int)(pow(2.0, 15.0));
  if (0 && mri->width > shortmax) {
    printf("NIFTI FORMAT ERROR: nrows %d in volume exceeds %d\n", mri->height, shortmax);
    exit(1);
  }
  if (mri->height > shortmax) {
    printf("NIFTI FORMAT ERROR: nrows %d in volume exceeds %d\n", mri->height, shortmax);
    exit(1);
  }
  if (mri->depth > shortmax) {
    printf("NIFTI FORMAT ERROR: nslices %d in volume exceeds %d\n", mri->depth, shortmax);
    exit(1);
  }
  if (mri->nframes > shortmax) {
    printf("NIFTI FORMAT ERROR:  nframes %d in volume exceeds %d\n", mri->nframes, shortmax);
    exit(1);
  }

  memset(&hdr, 0x00, sizeof(hdr));

  hdr.sizeof_hdr = 348;
  hdr.dim_info = 0;

  for (t = 0; t < 8; t++) {
    hdr.dim[t] = 1;
    hdr.pixdim[t] = 0;
  }
  hdr.dim[0] = 5;

  hdr.dim[1] = mri->width;
  hdr.dim[2] = mri->height;
  hdr.dim[3] = mri->depth;
  hdr.dim[4] = 1;
  hdr.dim[5] = mri->nframes;
  hdr.pixdim[1] = mri->xsize;
  hdr.pixdim[2] = mri->ysize;
  hdr.pixdim[3] = mri->zsize;

  if (mri->type == MRI_FLOAT) {
    hdr.datatype = DT_FLOAT;
    hdr.bitpix = 32;
  }
  else {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "itkMorphWrite(): unknown data type %d", mri->type));
  }

  hdr.intent_code = NIFTI_INTENT_VECTOR;
  hdr.intent_name[0] = '\0';
  hdr.vox_offset = 352;  // 352 is the min, dont use sizeof(hdr); See below
  hdr.regular = 'r';
  hdr.scl_slope = 0.0;
  hdr.slice_code = 0;
  hdr.xyzt_units = NIFTI_UNITS_MM;
  hdr.cal_max = 0.0;
  hdr.cal_min = 0.0;
  hdr.toffset = 0;
  sprintf(hdr.descrip, "FreeSurfer %s", __DATE__);

  /* set the nifti header qform values */
  error = mriToNiftiQform(mri, &hdr);
  if (error != NO_ERROR) return (error);
  hdr.qform_code = NIFTI_XFORM_ALIGNED_ANAT;

  /* set the nifti header sform values */
  // This just copies the vox2ras into the sform
  mriToNiftiSform(mri, &hdr);

  memmove(hdr.magic, NII_MAGIC, 4);

  fp = znzopen(fname, "w", use_compression);
  if (fp == NULL) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "itkMorphWrite(): error opening file %s", fname));
  }

  // White the header
  if (znzwrite(&hdr, sizeof(hdr), 1, fp) != 1) {
    znzclose(fp);
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "itkMorphWrite(): error writing header to %s", fname));
  }

  // Fill in space to the voxel offset
  nfill = (int)hdr.vox_offset - sizeof(hdr);
  chbuf = (char *)calloc(nfill, sizeof(char));
  if ((int)znzwrite(chbuf, sizeof(char), nfill, fp) != nfill) {
    znzclose(fp);
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "itkMorphWrite(): error writing data to %s", fname));
  }
  free(chbuf);

  // Now dump the pixel data
  for (t = 0; t < mri->nframes; t++)
    for (k = 0; k < mri->depth; k++) {
      for (j = 0; j < mri->height; j++) {
        buf = &MRIseq_vox(mri, 0, j, k, t);
        if ((int)znzwrite(buf, hdr.bitpix / 8, mri->width, fp) != mri->width) {
          znzclose(fp);
          errno = 0;
          ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "itkMorphWrite(): error writing data to %s", fname));
        }
      }
      exec_progress_callback(k, mri->depth, t, mri->nframes);
    }

  znzclose(fp);

  return (NO_ERROR);

} /* end itkMorphWrite() */

static int niftiSformToMri(MRI *mri, struct nifti_1_header *hdr)
{
  int err;
  MATRIX *sform;
  sform = MatrixConstVal(0, 4, 4, NULL);
  /*
  R = srow_x[0] * c + srow_x[1] * r + srow_x[2] * s + srow_x[3]
  A = srow_y[0] * c + srow_y[1] * r + srow_y[2] * s + srow_y[3]
  S = srow_z[0] * c + srow_z[1] * r + srow_z[2] * s + srow_z[3]
  */
  sform->rptr[1][1] = hdr->srow_x[0];
  sform->rptr[1][2] = hdr->srow_x[1];
  sform->rptr[1][3] = hdr->srow_x[2];
  sform->rptr[2][1] = hdr->srow_y[0];
  sform->rptr[2][2] = hdr->srow_y[1];
  sform->rptr[2][3] = hdr->srow_y[2];
  sform->rptr[3][1] = hdr->srow_z[0];
  sform->rptr[3][2] = hdr->srow_z[1];
  sform->rptr[3][3] = hdr->srow_z[2];
  sform->rptr[1][4] = hdr->srow_x[3];
  sform->rptr[2][4] = hdr->srow_y[3];
  sform->rptr[3][4] = hdr->srow_z[3];
  // MatrixPrint(stdout,sform);
  err = MRIsetVox2RASFromMatrix(mri, sform);
  MatrixFree(&sform);
  return (err);
} /* end niftiSformToMri() */

static int niftiQformToMri(MRI *mri, struct nifti_1_header *hdr)
{
  float a, b, c, d;
  float r11, r12, r13;
  float r21, r22, r23;
  float r31, r32, r33;

  b = hdr->quatern_b;
  c = hdr->quatern_c;
  d = hdr->quatern_d;

  /* following quatern_to_mat44() */

  a = 1.0 - (b * b + c * c + d * d);
  if (a < 1.0e-7) {
    a = 1.0 / sqrt(b * b + c * c + d * d);
    b *= a;
    c *= a;
    d *= a;
    a = 0.0;
  }
  else
    a = sqrt(a);

  r11 = a * a + b * b - c * c - d * d;
  r12 = 2.0 * b * c - 2.0 * a * d;
  r13 = 2.0 * b * d + 2.0 * a * c;
  r21 = 2.0 * b * c + 2.0 * a * d;
  r22 = a * a + c * c - b * b - d * d;
  r23 = 2.0 * c * d - 2.0 * a * b;
  r31 = 2.0 * b * d - 2 * a * c;
  r32 = 2.0 * c * d + 2 * a * b;
  r33 = a * a + d * d - c * c - b * b;

  if (hdr->pixdim[0] < 0.0) {
    r13 = -r13;
    r23 = -r23;
    r33 = -r33;
  }

  mri->x_r = r11;
  mri->y_r = r12;
  mri->z_r = r13;
  mri->x_a = r21;
  mri->y_a = r22;
  mri->z_a = r23;
  mri->x_s = r31;
  mri->y_s = r32;
  mri->z_s = r33;

  /**/
  /* c_ras */
  /*

  [ r ]   [ mri->xsize * mri->x_r
  mri->ysize * mri->y_r    mri->zsize * mri->z_r  r_offset ]   [ i ]
  [ a ] = [ mri->xsize * mri->x_a
  mri->ysize * mri->y_a    mri->zsize * mri->z_a  a_offset ] * [ j ]
  [ s ]   [ mri->xsize * mri->x_s
  mri->ysize * mri->y_s    mri->zsize * mri->z_s  s_offset ]   [ k ]
  [ 1 ]   [            0
  0                        0              1     ]   [ 1 ]


  */

  mri->c_r = (mri->xsize * mri->x_r) * (mri->width / 2.0) + (mri->ysize * mri->y_r) * (mri->height / 2.0) +
             (mri->zsize * mri->z_r) * (mri->depth / 2.0) + hdr->qoffset_x;

  mri->c_a = (mri->xsize * mri->x_a) * (mri->width / 2.0) + (mri->ysize * mri->y_a) * (mri->height / 2.0) +
             (mri->zsize * mri->z_a) * (mri->depth / 2.0) + hdr->qoffset_y;

  mri->c_s = (mri->xsize * mri->x_s) * (mri->width / 2.0) + (mri->ysize * mri->y_s) * (mri->height / 2.0) +
             (mri->zsize * mri->z_s) * (mri->depth / 2.0) + hdr->qoffset_z;

  return (NO_ERROR);

} /* end niftiQformToMri() */

/*---------------------------------------------------------------------
  mriToNiftiSform() - just copies the vox2ras to sform and sets the
  xform code to scanner anat.
  ---------------------------------------------------------------------*/
static int mriToNiftiSform(MRI *mri, struct nifti_1_header *hdr)
{
  MATRIX *vox2ras;
  int c;
  vox2ras = MRIxfmCRS2XYZ(mri, 0);
  for (c = 0; c < 4; c++) {
    hdr->srow_x[c] = vox2ras->rptr[1][c + 1];
    hdr->srow_y[c] = vox2ras->rptr[2][c + 1];
    hdr->srow_z[c] = vox2ras->rptr[3][c + 1];
  }
  hdr->sform_code = NIFTI_XFORM_SCANNER_ANAT;
  MatrixFree(&vox2ras);
  return (0);
}

/*---------------------------------------------------------------------*/
static int mriToNiftiQform(MRI *mri, struct nifti_1_header *hdr)
{
  MATRIX *i_to_r;
  float r11, r12, r13;
  float r21, r22, r23;
  float r31, r32, r33;
  float qfac = -100000;
  float a, b, c, d;
  float xd, yd, zd;
  float r_det;

  /*

  nifti1.h (x, y, z are r, a, s, respectively):

  [ x ]   [ R11 R12 R13 ] [        pixdim[1] * i ]   [ qoffset_x ]
  [ y ] = [ R21 R22 R23 ] [        pixdim[2] * j ] + [ qoffset_y ]
  [ z ]   [ R31 R32 R33 ] [ qfac * pixdim[3] * k ]   [ qoffset_z ]

  mri.c extract_r_to_i():

  [ r ]   [ mri->xsize * mri->x_r
  mri->ysize * mri->y_r    mri->zsize * mri->z_r  r_offset ]   [ i ]
  [ a ] = [ mri->xsize * mri->x_a
  mri->ysize * mri->y_a    mri->zsize * mri->z_a  a_offset ] * [ j ]
  [ s ]   [ mri->xsize * mri->x_s
  mri->ysize * mri->y_s    mri->zsize * mri->z_s  s_offset ]   [ k ]
  [ 1 ]   [            0
  0                        0              1     ]   [ 1 ]

  where [ras]_offset are calculated by satisfying (r, a, s)
  = (mri->c_r, mri->c_a, mri->c_s)
  when (i, j, k) = (mri->width/2, mri->height/2, mri->depth/2)

  so our mapping is simple:

  (done outside this function)
  pixdim[1] mri->xsize
  pixdim[2] mri->ysize
  pixdim[3] mri->zsize

  R11 = mri->x_r  R12 = mri->y_r  R13 = mri->z_r
  R21 = mri->x_a  R22 = mri->y_a  R23 = mri->z_a
  R31 = mri->x_s  R32 = mri->y_s  R33 = mri->z_s

  qoffset_x = r_offset
  qoffset_y = a_offset
  qoffset_z = s_offset

  qfac (pixdim[0]) = 1

  unless det(R) == -1, in which case we substitute:

  R13 = -mri->z_r
  R23 = -mri->z_a
  R33 = -mri->z_s

  qfac (pixdim[0]) = -1

  we follow mat44_to_quatern() from nifti1_io.c

  for now, assume orthonormality in mri->[xyz]_[ras]

  */

  r11 = mri->x_r;
  r12 = mri->y_r;
  r13 = mri->z_r;
  r21 = mri->x_a;
  r22 = mri->y_a;
  r23 = mri->z_a;
  r31 = mri->x_s;
  r32 = mri->y_s;
  r33 = mri->z_s;

  r_det = +r11 * (r22 * r33 - r32 * r23) - r12 * (r21 * r33 - r31 * r23) + r13 * (r21 * r32 - r31 * r22);

  if (r_det == 0.0) {
    printf(
        "WARNING: bad orientation matrix (determinant = 0) "
        "in nifti1 file ...\n");
    printf(" ... continuing.\n");
  }
  else if (r_det > 0.0)
    qfac = 1.0;
  else {
    r13 = -r13;
    r23 = -r23;
    r33 = -r33;
    qfac = -1.0;
  }

  /* following mat44_to_quatern() */

  a = r11 + r22 + r33 + 1.0;

  if (a > 0.5) {
    a = 0.5 * sqrt(a);
    b = 0.25 * (r32 - r23) / a;
    c = 0.25 * (r13 - r31) / a;
    d = 0.25 * (r21 - r12) / a;
  }
  else {
    xd = 1.0 + r11 - (r22 + r33);
    yd = 1.0 + r22 - (r11 + r33);
    zd = 1.0 + r33 - (r11 + r22);

    if (xd > 1.0) {
      b = 0.5 * sqrt(xd);
      c = 0.25 * (r12 + r21) / b;
      d = 0.25 * (r13 + r31) / b;
      a = 0.25 * (r32 - r23) / b;
    }
    else if (yd > 1.0) {
      c = 0.5 * sqrt(yd);
      b = 0.25 * (r12 + r21) / c;
      d = 0.25 * (r23 + r32) / c;
      a = 0.25 * (r13 - r31) / c;
    }
    else {
      d = 0.5 * sqrt(zd);
      b = 0.25 * (r13 + r31) / d;
      c = 0.25 * (r23 + r32) / d;
      a = 0.25 * (r21 - r12) / d;
    }

    if (a < 0.0) {
      a = -a;
      b = -b;
      c = -c;
      d = -d;
    }
  }

  hdr->qform_code = NIFTI_XFORM_SCANNER_ANAT;

  hdr->pixdim[0] = qfac;

  hdr->quatern_b = b;
  hdr->quatern_c = c;
  hdr->quatern_d = d;

  i_to_r = extract_i_to_r(mri);
  if (i_to_r == NULL) return (ERROR_BADPARM);

  hdr->qoffset_x = *MATRIX_RELT(i_to_r, 1, 4);
  hdr->qoffset_y = *MATRIX_RELT(i_to_r, 2, 4);
  hdr->qoffset_z = *MATRIX_RELT(i_to_r, 3, 4);

  MatrixFree(&i_to_r);

  return (NO_ERROR);

} /* end mriToNiftiQform() */

static void swap_nifti_1_header(struct nifti_1_header *hdr)
{
  int i;

  hdr->sizeof_hdr = swapInt(hdr->sizeof_hdr);

  for (i = 0; i < 8; i++) hdr->dim[i] = swapShort(hdr->dim[i]);

  hdr->intent_p1 = swapFloat(hdr->intent_p1);
  hdr->intent_p2 = swapFloat(hdr->intent_p2);
  hdr->intent_p3 = swapFloat(hdr->intent_p3);
  hdr->intent_code = swapShort(hdr->intent_code);
  hdr->datatype = swapShort(hdr->datatype);
  hdr->bitpix = swapShort(hdr->bitpix);
  hdr->slice_start = swapShort(hdr->slice_start);

  for (i = 0; i < 8; i++) hdr->pixdim[i] = swapFloat(hdr->pixdim[i]);

  hdr->vox_offset = swapFloat(hdr->vox_offset);
  hdr->scl_slope = swapFloat(hdr->scl_slope);
  hdr->scl_inter = swapFloat(hdr->scl_inter);
  hdr->slice_end = swapShort(hdr->slice_end);
  hdr->cal_max = swapFloat(hdr->cal_max);
  hdr->cal_min = swapFloat(hdr->cal_min);
  hdr->slice_duration = swapFloat(hdr->slice_duration);
  hdr->toffset = swapFloat(hdr->toffset);
  hdr->qform_code = swapShort(hdr->qform_code);
  hdr->sform_code = swapShort(hdr->sform_code);
  hdr->quatern_b = swapFloat(hdr->quatern_b);
  hdr->quatern_c = swapFloat(hdr->quatern_c);
  hdr->quatern_d = swapFloat(hdr->quatern_d);
  hdr->qoffset_x = swapFloat(hdr->qoffset_x);
  hdr->qoffset_y = swapFloat(hdr->qoffset_y);
  hdr->qoffset_z = swapFloat(hdr->qoffset_z);

  for (i = 0; i < 4; i++) hdr->srow_x[i] = swapFloat(hdr->srow_x[i]);

  for (i = 0; i < 4; i++) hdr->srow_y[i] = swapFloat(hdr->srow_y[i]);

  for (i = 0; i < 4; i++) hdr->srow_z[i] = swapFloat(hdr->srow_z[i]);

  return;

} /* end swap_nifti_1_header */

MRI *MRIreadGeRoi(const char *fname, int n_slices)
{
  MRI *mri;
  int i;
  char prefix[STRLEN], postfix[STRLEN];
  int n_digits;
  FILE *fp;
  int width, height;
  std::string fname_use;
  int read_one_flag;
  int pixel_data_offset;
  int y;

  if ((fp = fopen(fname, "r")) == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "MRIreadGeRoi(): error opening file %s", fname));
  }

  fseek(fp, 8, SEEK_SET);
  if (fread(&width, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "MRIreadGeRoi(): could not read file");
  }
  width = orderIntBytes(width);
  if (fread(&height, 4, 1, fp) != 1) {
     ErrorPrintf(ERROR_BADFILE, "MRIreadGeRoi(): could not read file");
  }
  height = orderIntBytes(height);

  fclose(fp);

  for (i = strlen(fname); i >= 0 && fname[i] != '/'; i--)
    ;
  i++;

  n_digits = 0;
  for (; fname[i] != '\0' && n_digits < 3; i++) {
    if (isdigit(fname[i]))
      n_digits++;
    else
      n_digits = 0;
  }

  if (n_digits < 3) {
    errno = 0;
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MRIreadGeRoi(): bad GE file name (couldn't find "
                 "three consecutive digits in the base file name)"));
  }

  strcpy(prefix, fname);
  prefix[i - 3] = '\0';

  strcpy(postfix, &(fname[i]));

  printf("%s---%s\n", prefix, postfix);

  mri = MRIalloc(width, height, n_slices, MRI_SHORT);

  if (mri == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_NOMEMORY, "MRIreadGeRoi(): couldn't allocate MRI structure"));
  }

  MRIclear(mri);

  read_one_flag = FALSE;

  for (i = 0; i < n_slices; i++) {
    std::stringstream tmp;
    tmp << prefix << std::setw(3) << std::setfill('0') << i << postfix;
    fname_use = tmp.str();
    if ((fp = fopen(fname_use.c_str(), "r")) != NULL) {
      fseek(fp, 4, SEEK_SET);
      if (fread(&pixel_data_offset, 4, 1, fp) != 1) {
         ErrorPrintf(ERROR_BADFILE, "MRIreadGeRoi(): could not read file");
      }
      pixel_data_offset = orderIntBytes(pixel_data_offset);
      fseek(fp, pixel_data_offset, SEEK_SET);

      for (y = 0; y < mri->height; y++) {
        if ((int)fread(mri->slices[i][y], 2, mri->width, fp) != mri->width) {
          fclose(fp);
          MRIfree(&mri);
          errno = 0;
          ErrorReturn(NULL, (ERROR_BADFILE, "MRIreadGeRoi(): error reading from file file %s", fname_use.c_str()));
        }
#if (BYTE_ORDER == LITTLE_ENDIAN)
	{
	  std::vector<short> tmp(2*mri->width);
	  // Note:
	  // void swab(const void *from, void *to, ssize_t n);
	  // void *memcpy(void *dest, const void *src, size_t n);
	  // Because consistency is the hobgoblin of small minds...
	  swab(mri->slices[i][y], tmp.data(), (size_t)(2 * mri->width));
	  memcpy(mri->slices[i][y], tmp.data(), (size_t)(2 * mri->width));
	}
#endif
      }

      fclose(fp);
      read_one_flag = TRUE;
    }

    exec_progress_callback(0, n_slices, 0, 1);
  }

  if (!read_one_flag) {
    MRIfree(&mri);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "MRIreadGeRoi(): didn't read any ROI files"));
  }

  return (mri);

} /* end MRIreadGeRoi() */

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

static int data_size[] = {1, 4, 4, 4, 2};

static MRI *sdtRead(const char *fname, int read_volume)
{
  char header_fname[STR_LEN];
  char line[STR_LEN];
  char *colon, *dot;
  FILE *fp;
  MRI *mri;
  int ndim = -1, data_type = -1;
  int dim[4];
  float xsize = 1.0, ysize = 1.0, zsize = 1.0, dummy_size;
  int orientation = MRI_CORONAL;

  dim[0] = -1;
  dim[1] = -1;
  dim[2] = -1;
  dim[3] = -1;

  /* form the header file name */
  strcpy(header_fname, fname);
  if ((dot = strrchr(header_fname, '.')))
    sprintf(dot + 1, "spr");
  else
    strcat(header_fname, ".spr");

  /* open the header */
  if ((fp = fopen(header_fname, "r")) == NULL) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): could not open header file %s\n", fname, header_fname));
  }

  while (1)  // !feof(fp))
  {
    if (!fgets(line, STR_LEN, fp) && ferror(fp)){
      ErrorPrintf(ERROR_BADFILE, "stdRead(%s): could not read file", fname);
    }
    if (feof(fp))  // wow.  we read too many.  get out
      break;

    if ((colon = strchr(line, ':'))) {
      *colon = '\0';
      colon++;
      if (strcmp(line, "numDim") == 0) {
        sscanf(colon, "%d", &ndim);
        if (ndim < 3 || ndim > 4) {
          fclose(fp);
          errno = 0;
          ErrorReturn(NULL,
                      (ERROR_UNSUPPORTED,
                       "sdtRead(%s): only 3 or 4 dimensions "
                       "supported (numDim = %d)\n",
                       fname,
                       ndim));
        }
      }
      else if (strcmp(line, "dim") == 0) {
        if (ndim == -1) {
          fclose(fp);
          errno = 0;
          ErrorReturn(NULL,
                      (ERROR_BADFILE, "sdtRead(%s): 'dim' before 'numDim' in header file %s\n", fname, header_fname));
        }
        if (ndim == 3) {
          sscanf(colon, "%d %d %d", &dim[0], &dim[1], &dim[2]);
          dim[3] = 1;
        }
        else {
          sscanf(colon, "%d %d %d %d", &dim[0], &dim[1], &dim[2], &dim[3]);
          if (dim[3] != 1) {
            fclose(fp);
            errno = 0;
            ErrorReturn(NULL,
                        (ERROR_UNSUPPORTED,
                         "sdtRead(%s): nframes != 1 unsupported "
                         "for sdt (dim(4) = %d)\n",
                         fname,
                         dim[3]));
          }
        }
      }
      else if (strcmp(line, "dataType") == 0) {
        while (isspace((int)*colon)) colon++;
        if (strncmp(colon, "BYTE", 4) == 0)
          data_type = MRI_UCHAR;
        else if (strncmp(colon, "WORD", 4) == 0)
          data_type = MRI_SHORT;
        else if (strncmp(colon, "LWORD", 5) == 0)
          data_type = MRI_INT;
        else if (strncmp(colon, "REAL", 4) == 0)
          data_type = MRI_FLOAT;
        else if (strncmp(colon, "COMPLEX", 7) == 0) {
          fclose(fp);
          errno = 0;
          ErrorReturn(NULL, (ERROR_UNSUPPORTED, "sdtRead(%s): unsupported data type '%s'\n", fname, colon));
        }
        else {
          fclose(fp);
          errno = 0;
          ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): unknown data type '%s'\n", fname, colon));
        }
      }
      else if (strcmp(line, "interval") == 0) {
        if (ndim == 3)
          sscanf(colon, "%f %f %f", &xsize, &ysize, &zsize);
        else
          sscanf(colon, "%f %f %f %f", &xsize, &ysize, &zsize, &dummy_size);
        xsize *= 10.0;
        ysize *= 10.0;
        zsize *= 10.0;
      }
      else if (strcmp(line, "sdtOrient") == 0) {
        while (isspace((int)*colon)) colon++;
        if (strncmp(colon, "sag", 3) == 0)
          orientation = MRI_SAGITTAL;
        else if (strncmp(colon, "ax", 2) == 0)
          orientation = MRI_HORIZONTAL;
        else if (strncmp(colon, "cor", 3) == 0)
          orientation = MRI_CORONAL;
        else {
          fclose(fp);
          errno = 0;
          ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): unknown orientation %s\n", fname, colon));
        }
      }
      else {
      }
    }
  }

  fclose(fp);

  if (data_type == -1) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): data type undefined\n", fname));
  }
  if (dim[0] == -1 || dim[1] == -1 || dim[2] == -1) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): one or more dimensions undefined\n", fname));
  }

  if (read_volume) {
    if ((fp = fopen(fname, "r")) == NULL) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): error opening data file %s\n", fname, fname));
    }

    mri = MRIreadRaw(fp, dim[0], dim[1], dim[2], data_type);

    if (mri == NULL) return (NULL);

    fclose(fp);
  }
  else {
    mri = MRIallocHeader(dim[0], dim[1], dim[2], data_type, dim[3]);
    if (mri == NULL) return (NULL);
  }

  mri->xsize = xsize;
  mri->ysize = ysize;
  mri->zsize = zsize;

  setDirectionCosine(mri, orientation);

  mri->thick = mri->zsize;
  mri->xend = mri->xsize * mri->width / 2.;
  mri->yend = mri->ysize * mri->height / 2.;
  mri->zend = mri->zsize * mri->depth / 2.;
  mri->xstart = -mri->xend;
  mri->ystart = -mri->yend;
  mri->zstart = -mri->zend;

  mri->imnr0 = 1;
  mri->imnr1 = dim[2];

  mri->ps = 1.0 /*0.001*/;
  mri->tr = 0;
  mri->te = 0;
  mri->ti = 0;

  strcpy(mri->fname, fname);

  return (mri);

} /* end sdtRead() */

MRI *MRIreadRaw(FILE *fp, int width, int height, int depth, int type)
{
  MRI *mri;
  BUFTYPE *buf;
  int slice, pixels;
  int i;

  mri = MRIalloc(width, height, depth, type);
  if (!mri) return (NULL);

  pixels = width * height;
  buf = (BUFTYPE *)calloc(pixels, data_size[type]);

  /* every width x height pixels should be another slice */
  for (slice = 0; slice < depth; slice++) {
    if ((int)fread(buf, data_size[type], pixels, fp) != pixels) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADFILE, "%s: could not read %dth slice (%d)", Progname, slice, pixels));
    }
    if (type == 0) local_buffer_to_image(buf, mri, slice, 0);
    if (type == 1) {
      for (i = 0; i < pixels; i++) ((int *)buf)[i] = orderIntBytes(((int *)buf)[i]);
      int_local_buffer_to_image((int *)buf, mri, slice, 0);
    }
    if (type == 2) {
      for (i = 0; i < pixels; i++) ((long32 *)buf)[i] = orderLong32Bytes(((long32 *)buf)[i]);
      long32_local_buffer_to_image((long32 *)buf, mri, slice, 0);
    }
    if (type == 3) {
      for (i = 0; i < pixels; i++) ((float *)buf)[i] = orderFloatBytes(((float *)buf)[i]);
      float_local_buffer_to_image((float *)buf, mri, slice, 0);
    }
    if (type == 4) {
      for (i = 0; i < pixels; i++) ((short *)buf)[i] = orderShortBytes(((short *)buf)[i]);
      short_local_buffer_to_image((short *)buf, mri, slice, 0);
    }

    exec_progress_callback(slice, depth, 0, 1);
  }

  free(buf);
  return (mri);
}

static void int_local_buffer_to_image(int *buf, MRI *mri, int slice, int frame)
{
  int y, width, height;
  int *pslice;

  width = mri->width;
  height = mri->height;
  for (y = 0; y < height; y++) {
    pslice = &MRIIseq_vox(mri, 0, y, slice, frame);
    memmove(pslice, buf, width * sizeof(int));
    buf += width;
  }
}


static void long32_local_buffer_to_image(long32 *buf, MRI *mri, int slice, int frame)
{
  int y, width, height;
  long32 *pslice;

  width = mri->width;
  height = mri->height;
  for (y = 0; y < height; y++) {
    pslice = &MRILseq_vox(mri, 0, y, slice, frame);
    memmove(pslice, buf, width * sizeof(long));
    buf += width;
  }
}

static void float_local_buffer_to_image(float *buf, MRI *mri, int slice, int frame)
{
  int y, width, height;
  float *pslice;

  width = mri->width;
  height = mri->height;
  for (y = 0; y < height; y++) {
    pslice = &MRIFseq_vox(mri, 0, y, slice, frame);
    memmove(pslice, buf, width * sizeof(float));
    buf += width;
  }
}

static void short_local_buffer_to_image(short *buf, MRI *mri, int slice, int frame)
{
  int y, width, height;
  short *pslice;

  width = mri->width;
  height = mri->height;
  for (y = 0; y < height; y++) {
    pslice = &MRISseq_vox(mri, 0, y, slice, frame);
    memmove(pslice, buf, width * sizeof(short));
    buf += width;
  }
}

static void local_buffer_to_image(BUFTYPE *buf, MRI *mri, int slice, int frame)
{
  int y, width, height;
  BUFTYPE *pslice;

  width = mri->width;
  height = mri->height;
  for (y = 0; y < height; y++) {
    pslice = &MRIseq_vox(mri, 0, y, slice, frame);
    memmove(pslice, buf, width * sizeof(BUFTYPE));
    buf += width;
  }
}

int znzTAGwriteMRIframes(znzFile fp, MRI *mri)
{
  long long len = 0, fstart, fend, here;
  int fno, i;
  MRI_FRAME *frame;
  char *buf;

  // write some extra space so that we have enough room (can't seek in zz files)
  len = 10 * mri->nframes * sizeof(MRI_FRAME);
  znzTAGwriteStart(fp, TAG_MRI_FRAME, &fstart, len);
  here = znztell(fp);
  for (fno = 0; fno < mri->nframes; fno++) {
    frame = &mri->frames[fno];
    znzwriteInt(frame->type, fp);
    znzwriteFloat(frame->TE, fp);
    znzwriteFloat(frame->TR, fp);
    znzwriteFloat(frame->flip, fp);
    znzwriteFloat(frame->TI, fp);
    znzwriteFloat(frame->TD, fp);
    znzwriteFloat(frame->TM, fp);
    znzwriteInt(frame->sequence_type, fp);
    znzwriteFloat(frame->echo_spacing, fp);
    znzwriteFloat(frame->echo_train_len, fp);
    for (i = 0; i < 3; i++) znzwriteFloat(frame->read_dir[i], fp);
    for (i = 0; i < 3; i++) znzwriteFloat(frame->pe_dir[i], fp);
    for (i = 0; i < 3; i++) znzwriteFloat(frame->slice_dir[i], fp);
    znzwriteInt(frame->label, fp);
    znzwrite(frame->name, sizeof(char), STRLEN, fp);
    znzwriteInt(frame->dof, fp);
    if (frame->m_ras2vox && frame->m_ras2vox->rows > 0)
      znzWriteMatrix(fp, frame->m_ras2vox, 0);
    else {
      MATRIX *m = MatrixAlloc(4, 4, MATRIX_REAL);
      znzWriteMatrix(fp, m, 0);
      MatrixFree(&m);
    }
    znzwriteFloat(frame->thresh, fp);
    znzwriteInt(frame->units, fp);
    if (frame->type == FRAME_TYPE_DIFFUSION_AUGMENTED)  // also store diffusion info
    {
      znzwriteDouble(frame->DX, fp);
      znzwriteDouble(frame->DY, fp);
      znzwriteDouble(frame->DZ, fp);

      znzwriteDouble(frame->DR, fp);
      znzwriteDouble(frame->DP, fp);
      znzwriteDouble(frame->DS, fp);
      znzwriteDouble(frame->bvalue, fp);
      znzwriteDouble(frame->TM, fp);

      znzwriteLong(frame->D1_ramp, fp);
      znzwriteLong(frame->D1_flat, fp);
      znzwriteDouble(frame->D1_amp, fp);

      znzwriteLong(frame->D2_ramp, fp);
      znzwriteLong(frame->D2_flat, fp);
      znzwriteDouble(frame->D2_amp, fp);

      znzwriteLong(frame->D3_ramp, fp);
      znzwriteLong(frame->D3_flat, fp);
      znzwriteDouble(frame->D3_amp, fp);

      znzwriteLong(frame->D4_ramp, fp);
      znzwriteLong(frame->D4_flat, fp);
      znzwriteDouble(frame->D4_amp, fp);
    }
  }
  fend = znztell(fp);
  len -= (fend - here);  // unused space
  if (len > 0) {
    buf = (char *)calloc(len, sizeof(char));
    znzwrite(buf, len, sizeof(char), fp);
    free(buf);
  }
  znzTAGwriteEnd(fp, fend);

  return (NO_ERROR);
}
int znzTAGreadMRIframes(znzFile fp, MRI *mri, long len)
{
  int fno, i;
  long long fstart, fend;
  MRI_FRAME *frame;
  char *buf;

  fstart = znztell(fp);
  for (fno = 0; fno < mri->nframes; fno++) {
    frame = &mri->frames[fno];
    frame->type = znzreadInt(fp);
    frame->TE = znzreadFloat(fp);
    frame->TR = znzreadFloat(fp);
    frame->flip = znzreadFloat(fp);
    frame->TI = znzreadFloat(fp);
    frame->TD = znzreadFloat(fp);
    frame->TM = znzreadFloat(fp);
    frame->sequence_type = znzreadInt(fp);
    frame->echo_spacing = znzreadFloat(fp);
    frame->echo_train_len = znzreadFloat(fp);
    for (i = 0; i < 3; i++) frame->read_dir[i] = znzreadFloat(fp);
    for (i = 0; i < 3; i++) frame->pe_dir[i] = znzreadFloat(fp);
    for (i = 0; i < 3; i++) frame->slice_dir[i] = znzreadFloat(fp);
    frame->label = znzreadInt(fp);
    znzread(frame->name, sizeof(char), STRLEN, fp);
    frame->dof = znzreadInt(fp);
    frame->m_ras2vox = znzReadMatrix(fp);

    frame->thresh = znzreadFloat(fp);
    frame->units = znzreadInt(fp);
    if (frame->type == FRAME_TYPE_DIFFUSION_AUGMENTED) {
      frame->DX = znzreadDouble(fp);
      frame->DY = znzreadDouble(fp);
      frame->DZ = znzreadDouble(fp);

      frame->DR = znzreadDouble(fp);
      frame->DP = znzreadDouble(fp);
      frame->DS = znzreadDouble(fp);
      frame->bvalue = znzreadDouble(fp);
      frame->TM = znzreadDouble(fp);

      frame->D1_ramp = znzreadLong(fp);
      frame->D1_flat = znzreadLong(fp);
      frame->D1_amp = znzreadDouble(fp);

      frame->D2_ramp = znzreadLong(fp);
      frame->D2_flat = znzreadLong(fp);
      frame->D2_amp = znzreadDouble(fp);

      frame->D3_ramp = znzreadLong(fp);
      frame->D3_flat = znzreadLong(fp);
      frame->D3_amp = znzreadDouble(fp);

      frame->D4_ramp = znzreadLong(fp);
      frame->D4_flat = znzreadLong(fp);
      frame->D4_amp = znzreadDouble(fp);
    }
  }

  fend = znztell(fp);
  len -= (fend - fstart);
  if (len > 0) {
    buf = (char *)calloc(len, sizeof(char));
    znzread(buf, len, sizeof(char), fp);
    free(buf);
  }
  return (NO_ERROR);
}

#define UNUSED_SPACE_SIZE 256
#define USED_SPACE_SIZE (3 * sizeof(float) + 4 * 3 * sizeof(float))

#define MGH_VERSION 1

// declare function pointer
// static int (*myclose)(FILE *stream);

static MRI *mghRead(const char *fname, int read_volume, int frame)
{
  MRI *mri;
  znzFile fp;
  int start_frame, end_frame, width, height, depth, nframes, type, x, y, z, bpv, dof, bytes, version, ival,
      unused_space_size, good_ras_flag, i;
  BUFTYPE *buf;
  char unused_buf[UNUSED_SPACE_SIZE + 1];
  float fval, xsize, ysize, zsize, x_r, x_a, x_s, y_r, y_a, y_s, z_r, z_a, z_s, c_r, c_a, c_s, xfov, yfov, zfov;
  short sval;
  //  int tag_data_size;
  const char *ext;
  int gzipped = 0;
  int nread;
  int tag;

  ext = strrchr(fname, '.');
  int valid_ext = 0;
  if (ext) {
    ++ext;
    // if mgz, then it is compressed
    if (!stricmp(ext, "mgz") || strstr(fname, "mgh.gz")) {
      gzipped = 1;
      valid_ext = 1;
    }
    else if (!stricmp(ext, "mgh")) {
      valid_ext = 1;
    }
  }

  if (valid_ext) {
    fp = znzopen(fname, "rb", gzipped);
    if (znz_isnull(fp)) {
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADPARM, "mghRead(%s, %d): could not open file", fname, frame));
    }
  }
  else {
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "mghRead(%s, %d): could not open file.\n"
                 "Filename extension must be .mgh, .mgh.gz or .mgz",
                 fname,
                 frame));
  }

  /* keep the compiler quiet */
  xsize = ysize = zsize = 0;
  x_r = x_a = x_s = 0;
  y_r = y_a = y_s = 0;
  z_r = z_a = z_s = 0;
  c_r = c_a = c_s = 0;

  nread = znzreadIntEx(&version, fp);
  if (!nread) ErrorReturn(NULL, (ERROR_BADPARM, "mghRead(%s, %d): read error", fname, frame));

  width = znzreadInt(fp);
  height = znzreadInt(fp);
  depth = znzreadInt(fp);
  nframes = znzreadInt(fp);
  type = znzreadInt(fp);
  dof = znzreadInt(fp);

  unused_space_size = UNUSED_SPACE_SIZE - sizeof(short);

  good_ras_flag = znzreadShort(fp);
  if (good_ras_flag > 0) { /* has RAS and voxel size info */
    unused_space_size -= USED_SPACE_SIZE;
    xsize = znzreadFloat(fp);
    ysize = znzreadFloat(fp);
    zsize = znzreadFloat(fp);

    x_r = znzreadFloat(fp);
    x_a = znzreadFloat(fp);
    x_s = znzreadFloat(fp);
    y_r = znzreadFloat(fp);
    y_a = znzreadFloat(fp);
    y_s = znzreadFloat(fp);

    z_r = znzreadFloat(fp);
    z_a = znzreadFloat(fp);
    z_s = znzreadFloat(fp);
    c_r = znzreadFloat(fp);
    c_a = znzreadFloat(fp);
    c_s = znzreadFloat(fp);
  }
  /* so stuff can be added to the header in the future */
  znzread(unused_buf, sizeof(char), unused_space_size, fp);

  switch (type) {
    default:
    case MRI_FLOAT:
      bpv = sizeof(float);
      break;
    case MRI_UCHAR:
      bpv = sizeof(char);
      break;
    case MRI_SHORT:
      bpv = sizeof(short);
      break;
    case MRI_USHRT:
      bpv = sizeof(unsigned short);
      break;
    case MRI_INT:
      bpv = sizeof(int);
      break;
    case MRI_TENSOR:
      bpv = sizeof(float);
      nframes = 9;
      break;
  }
  bytes = width * height * bpv; /* bytes per slice */
  if (!read_volume) {
    mri = MRIallocHeader(width, height, depth, type, nframes);
    mri->dof = dof;
    mri->nframes = nframes;
    if (gzipped) {  // pipe cannot seek
      long count, total_bytes;
      uchar buf[STRLEN];

      total_bytes = (long)mri->nframes * width * height * depth * bpv;
      for (count = 0; count < total_bytes - STRLEN; count += STRLEN) znzread(buf, STRLEN, 1, fp);
      znzread(buf, total_bytes - count, 1, fp);
    }
    else
      znzseek(fp, (long)mri->nframes * width * height * depth * bpv, SEEK_CUR);
  }
  else {
    if (frame >= 0) {
      start_frame = end_frame = frame;
      if (gzipped) {  // pipe cannot seek
        long count;
        for (count = 0; count < (long)frame * width * height * depth * bpv; count++) znzgetc(fp);
      }
      else
        znzseek(fp, (long)frame * width * height * depth * bpv, SEEK_CUR);
      nframes = 1;
    }
    else { /* hack - # of frames < -1 means to only read in that
              many frames. Otherwise I would have had to change the whole
              MRIread interface and that was too much of a pain. Sorry.
           */
      if (frame < -1) nframes = frame * -1;

      start_frame = 0;
      end_frame = nframes - 1;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "read %d frames\n", nframes);
    }
    buf = (BUFTYPE *)calloc(bytes, sizeof(BUFTYPE));
    mri = MRIallocSequence(width, height, depth, type, nframes);
    mri->dof = dof;

    struct timespec begin, end;
    if (getenv("FS_MGZIO_TIMING"))
    {
      printf("INFO: Environment variable FS_MGZIO_TIMING set\n");
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
    }

    int USEVOXELBUF = 0;
    if (mri->ischunked && getenv("FS_MGZIO_USEVOXELBUF"))
    {
      USEVOXELBUF = 1;
      printf("INFO: Environment variable FS_MGZIO_USEVOXELBUF set\n");
      printf("mghRead(): read voxel buffer ...\n");
      if (buf) free(buf);

      int bytes_to_read = bytes * depth * (end_frame-start_frame+1);
      //BUFTYPE *bufsrc = (BUFTYPE*)calloc(bytes_to_read, sizeof(BUFTYPE));
      BUFTYPE *bufsrc = &MRIseq_vox(mri, 0, 0, 0, start_frame);
      if ((int)znzread(bufsrc, sizeof(BUFTYPE), bytes_to_read, fp) != bytes_to_read) {
        znzclose(fp);
        free(bufsrc);
        ErrorReturn(NULL, (ERROR_BADFILE, "mghRead (read voxel buffer): could not read %d bytes", bytes_to_read));
      }

#if (BYTE_ORDER == LITTLE_ENDIAN)
      if (bpv == 2) byteswapbufshort(bufsrc, bytes_to_read);
      if (bpv == 4) byteswapbuffloat(bufsrc, bytes_to_read);
      if (bpv == 8) byteswapbuffloat(bufsrc, bytes_to_read);
#endif

      //void *voxelbuf = &MRIseq_vox(mri, 0, 0, 0, start_frame);
      //memcpy(voxelbuf, bufsrc, bytes_to_read);

      //if (bufsrc) free(bufsrc);
    }
    else  // copy voxel point by point
    {
      for (frame = start_frame; frame <= end_frame; frame++) {
        for (z = 0; z < depth; z++) {
          if ((int)znzread(buf, sizeof(char), bytes, fp) != bytes) {
            // fclose(fp) ;
            znzclose(fp);
            free(buf);
            ErrorReturn(NULL, (ERROR_BADFILE, "mghRead(%s): could not read %d bytes at slice %d", fname, bytes, z));
          }
          switch (type) {
            case MRI_INT:
              for (i = y = 0; y < height; y++) {
                for (x = 0; x < width; x++, i++) {
                  ival = orderIntBytes(((int *)buf)[i]);
                  MRIIseq_vox(mri, x, y, z, frame - start_frame) = ival;
                }
              }
              break;
            case MRI_SHORT:
              for (i = y = 0; y < height; y++) {
                for (x = 0; x < width; x++, i++) {
                  sval = orderShortBytes(((short *)buf)[i]);
                  MRISseq_vox(mri, x, y, z, frame - start_frame) = sval;
                }
              }
              break;
            case MRI_USHRT:
              for (i = y = 0; y < height; y++) {
                for (x = 0; x < width; x++, i++) {
                  unsigned short usval = orderUShortBytes(((unsigned short *)buf)[i]);
                  MRIUSseq_vox(mri, x, y, z, frame - start_frame) = usval;
                }
              }
              break;
            case MRI_TENSOR:
            case MRI_FLOAT:
              for (i = y = 0; y < height; y++) {
                for (x = 0; x < width; x++, i++) {
                  fval = orderFloatBytes(((float *)buf)[i]);
                  MRIFseq_vox(mri, x, y, z, frame - start_frame) = fval;
                }
              }
              break;
            case MRI_UCHAR:
              local_buffer_to_image(buf, mri, z, frame - start_frame);
              break;
            default:
              errno = 0;
              ErrorReturn(NULL, (ERROR_UNSUPPORTED, "mghRead: unsupported type %d", mri->type));
              break;
          } // switch
          exec_progress_callback(z, depth, frame - start_frame, end_frame - start_frame + 1);
        } // depth
      } // frame
      if (buf) free(buf);
    } // end of copying voxel point by point

    if (getenv("FS_MGZIO_TIMING"))
    {
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
      printf("Total time (mghRead) = %ld.%09ld seconds%s\n", 
             (end.tv_nsec < begin.tv_nsec) ? (end.tv_sec - 1 - begin.tv_sec) : (end.tv_sec - begin.tv_sec), 
             (end.tv_nsec < begin.tv_nsec) ? (1000000000 + end.tv_nsec - begin.tv_nsec) : (end.tv_nsec - begin.tv_nsec),
             (USEVOXELBUF) ? " (USEVOXELBUF)" : "");
    }
  }

  if (good_ras_flag > 0) {
    mri->xsize = xsize;
    mri->ysize = ysize;
    mri->zsize = zsize;

    mri->ps = mri->xsize;
    mri->thick = mri->zsize;

    mri->x_r = x_r;
    mri->x_a = x_a;
    mri->x_s = x_s;

    mri->y_r = y_r;
    mri->y_a = y_a;
    mri->y_s = y_s;

    mri->z_r = z_r;
    mri->z_a = z_a;
    mri->z_s = z_s;

    mri->c_r = c_r;
    mri->c_a = c_a;
    mri->c_s = c_s;
    if (good_ras_flag > 0) mri->ras_good_flag = 1;
  }
  else {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr,
              "-----------------------------------------------------------------\n"
              "Could not find the direction cosine information.\n"
              "Will use the CORONAL orientation.\n"
              "If not suitable, please provide the information in %s.\n"
              "-----------------------------------------------------------------\n",
              fname);
    setDirectionCosine(mri, MRI_CORONAL);
  }
  // read TR, Flip, TE, TI, FOV
  if (znzreadFloatEx(&(mri->tr), fp)) {
    if (znzreadFloatEx(&fval, fp)) {
      mri->flip_angle = fval;
      // flip_angle is double. I cannot use the same trick.
      if (znzreadFloatEx(&(mri->te), fp)) {
        if (znzreadFloatEx(&(mri->ti), fp)) {
          if (znzreadFloatEx(&(mri->fov), fp)) {
            ;
	  }
	}
      }
    }
  }
  // tag reading
  if (getenv("FS_SKIP_TAGS") == NULL) {
    long long len;

    while (1) {
      tag = znzTAGreadStart(fp, &len);
      // printf("tag %d\n",tag);
      if (tag == 0) break;

      switch (tag) {
        case TAG_MRI_FRAME:
          if (znzTAGreadMRIframes(fp, mri, len) != NO_ERROR)
            fprintf(stderr, "couldn't read frame structure from file\n");
          break;

        case TAG_OLD_COLORTABLE:
          // if reading colortable fails, it will print its own error message
          mri->ct = znzCTABreadFromBinary(fp);
          break;

        case TAG_OLD_MGH_XFORM:
        case TAG_MGH_XFORM: {
          char *fnamedir;
          char tmpstr[1000];

          // First, try a path relative to fname (not the abs path)
          fnamedir = fio_dirname(fname);
          sprintf(tmpstr, "%s/transforms/talairach.xfm", fnamedir);
          free(fnamedir);
          fnamedir = NULL;
          znzgets(mri->transform_fname, len + 1, fp);
          // If this file exists, copy it to transform_fname
          if (FileExists(tmpstr)) strcpy(mri->transform_fname, tmpstr);
          if (FileExists(mri->transform_fname)) {
            // copied from corRead()
            if (input_transform_file(mri->transform_fname, &(mri->transform)) == NO_ERROR) {
              mri->linear_transform = get_linear_transform_ptr(&mri->transform);
              mri->inverse_linear_transform = get_inverse_linear_transform_ptr(&mri->transform);
              mri->free_transform = 1;
              if (DIAG_VERBOSE_ON) fprintf(stderr, "INFO: loaded talairach xform : %s\n", mri->transform_fname);
            }
            else {
              errno = 0;
              ErrorPrintf(ERROR_BAD_FILE, "error loading transform from %s", mri->transform_fname);
              mri->linear_transform = NULL;
              mri->inverse_linear_transform = NULL;
              mri->free_transform = 1;
              (mri->transform_fname)[0] = '\0';
            }
          }
          break;
        }
        case TAG_CMDLINE:
          if (mri->ncmds >= MAX_CMDS)
            ErrorExit(ERROR_NOMEMORY, "mghRead(%s): too many commands (%d) in file", fname, mri->ncmds);
          mri->cmdlines[mri->ncmds] = (char *)calloc(len + 1, sizeof(char));
          znzread(mri->cmdlines[mri->ncmds], sizeof(char), len, fp);
          mri->cmdlines[mri->ncmds][len] = 0;
          mri->ncmds++;
          break;

        case TAG_AUTO_ALIGN:
          mri->AutoAlign = znzReadAutoAlignMatrix(fp);
          break;

        case TAG_ORIG_RAS2VOX:
          mri->origRas2Vox = znzReadMatrix(fp);
          break;

        case TAG_PEDIR:
          mri->pedir = (char *)calloc(len + 1, sizeof(char));
          znzread(mri->pedir, sizeof(char), len, fp);
          break;

        case TAG_FIELDSTRENGTH:
          // znzreadFloatEx(&(mri->FieldStrength), fp); // Performs byte swap not in znzTAGwrite()
          znzTAGreadFloat(&(mri->FieldStrength), fp);
          break;

        default:
          znzTAGskip(fp, tag, (long long)len);
          break;
      }
    }
  }

  // fclose(fp) ;
  znzclose(fp);

  // xstart, xend, ystart, yend, zstart, zend are not stored
  mri->xstart = -mri->width / 2. * mri->xsize;
  mri->xend = mri->width / 2. * mri->xsize;
  mri->ystart = -mri->height / 2. * mri->ysize;
  mri->yend = mri->height / 2. * mri->ysize;
  mri->zstart = -mri->depth / 2. * mri->zsize;
  mri->zend = mri->depth / 2. * mri->zsize;
  xfov = mri->xend - mri->xstart;
  yfov = mri->yend - mri->ystart;
  zfov = mri->zend - mri->zstart;
  mri->fov = (xfov > yfov ? (xfov > zfov ? xfov : zfov) : (yfov > zfov ? yfov : zfov));

  strcpy(mri->fname, fname);
  return (mri);
}

static int mghWrite(MRI *mri, const char *fname, int frame)
{
  znzFile fp;
  int ival, start_frame, end_frame, x, y, z, width, height, depth, unused_space_size, flen;
  char buf[UNUSED_SPACE_SIZE + 1];
  float fval;
  short sval;
  int gzipped = 0;
  const char *ext;

  if (frame >= 0)
    start_frame = end_frame = frame;
  else {
    start_frame = 0;
    end_frame = mri->nframes - 1;
  }
  ////////////////////////////////////////////////////////////
  ext = strrchr(fname, '.');
  int valid_ext = 0;
  if (ext) {
    ++ext;
    // if mgz, then it is compressed
    if (!stricmp(ext, "mgz") || strstr(fname, "mgh.gz")) {
      gzipped = 1;
      valid_ext = 1;
    }
    else if (!stricmp(ext, "mgh")) {
      valid_ext = 1;
    }
  }
  if (valid_ext) {
    fp = znzopen(fname, "wb", gzipped);
    if (znz_isnull(fp)) {
      errno = 0;
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "mghWrite(%s, %d): could not open file", fname, frame));
    }
  }
  else {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "mghWrite: filename '%s' "
                 "needs to have an extension of .mgh or .mgz",
                 fname));
  }

  /* WARNING - adding or removing anything before nframes will
     cause mghAppend to fail.
  */
  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  // printf("(w,h,d) = (%d,%d,%d)\n", width, height, depth);
  znzwriteInt(MGH_VERSION, fp);
  znzwriteInt(mri->width, fp);
  znzwriteInt(mri->height, fp);
  znzwriteInt(mri->depth, fp);
  znzwriteInt(mri->nframes, fp);
  znzwriteInt(mri->type, fp);
  znzwriteInt(mri->dof, fp);

  unused_space_size = UNUSED_SPACE_SIZE - USED_SPACE_SIZE - sizeof(short);

  /* write RAS and voxel size info */
  znzwriteShort(mri->ras_good_flag ? 1 : -1, fp);
  znzwriteFloat(mri->xsize, fp);
  znzwriteFloat(mri->ysize, fp);
  znzwriteFloat(mri->zsize, fp);

  znzwriteFloat(mri->x_r, fp);
  znzwriteFloat(mri->x_a, fp);
  znzwriteFloat(mri->x_s, fp);

  znzwriteFloat(mri->y_r, fp);
  znzwriteFloat(mri->y_a, fp);
  znzwriteFloat(mri->y_s, fp);

  znzwriteFloat(mri->z_r, fp);
  znzwriteFloat(mri->z_a, fp);
  znzwriteFloat(mri->z_s, fp);

  znzwriteFloat(mri->c_r, fp);
  znzwriteFloat(mri->c_a, fp);
  znzwriteFloat(mri->c_s, fp);

  /* so stuff can be added to the header in the future */
  memset(buf, 0, UNUSED_SPACE_SIZE * sizeof(char));
  znzwrite(buf, sizeof(char), unused_space_size, fp);

  struct timespec begin, end;
  if (getenv("FS_MGZIO_TIMING"))
  {
    printf("INFO: Environment variable FS_MGZIO_TIMING set\n");
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
  }

  int USEVOXELBUF = 0;
  if (mri->ischunked && getenv("FS_MGZIO_USEVOXELBUF"))
  {
    USEVOXELBUF = 1;
    printf("INFO: Environment variable FS_MGZIO_USEVOXELBUF set\n");
    printf("mghWrite(): output voxel bufer ...\n");
    int bytes_per_voxel = MRIsizeof(mri->type);
    BUFTYPE *voxelbuf = &MRIseq_vox(mri, 0, 0, 0, start_frame);
    int bytes_to_write = bytes_per_voxel * width * height * depth * (end_frame-start_frame+1);

#if 0
    printf("<DEBUG> voxelbuf = %p - %p (%p - %p) %d (frame=%d, depth=%d, height=%d) (%d %d)\n", 
           voxelbuf, voxelbuf+(MRIsizeof(mri->type)*(width-1)), 
	   &MRISseq_vox(mri, 0, y, z, frame), &MRISseq_vox(mri, width-1, z, y, frame),
	   bytes_to_write, frame, y, z, (int)mri->bytes_per_vox, width);
#endif

#if (BYTE_ORDER == LITTLE_ENDIAN)
    if (bytes_per_voxel == 2) byteswapbufshort(voxelbuf, bytes_to_write);
    if (bytes_per_voxel == 4) byteswapbuffloat(voxelbuf, bytes_to_write);
    if (bytes_per_voxel == 8) byteswapbuffloat(voxelbuf, bytes_to_write);
#endif

    if ((int)znzwrite(voxelbuf, sizeof(BUFTYPE), bytes_to_write, fp) != bytes_to_write) 
    {
      errno = 0;
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "mghWrite: (output voxel buffer) could not write %d bytes", bytes_to_write));
    }
  }
  else
  {
    for (frame = start_frame; frame <= end_frame; frame++) {
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          switch (mri->type) {
            case MRI_SHORT:
              for (x = 0; x < width; x++) {
                if (z == 74 && y == 16 && x == 53) DiagBreak();
                sval = MRISseq_vox(mri, x, y, z, frame);
                znzwriteShort(sval, fp);
              }
              break;
            case MRI_USHRT:
              for (x = 0; x < width; x++) {
                if (z == 74 && y == 16 && x == 53) DiagBreak();
                unsigned short usval = MRIUSseq_vox(mri, x, y, z, frame);
                znzwriteUShort(usval, fp);
              }
              break;
            case MRI_INT:
              for (x = 0; x < width; x++) {
                if (z == 74 && y == 16 && x == 53) DiagBreak();
                ival = MRIIseq_vox(mri, x, y, z, frame);
                znzwriteInt(ival, fp);
              }
              break;
            case MRI_FLOAT:
              for (x = 0; x < width; x++) {
                if (z == 74 && y == 16 && x == 53) DiagBreak();
                // printf("mghWrite: MRI_FLOAT: curr (x, y, z, frame) = (%d, %d, %d, %d)\n", x, y, z, frame);
                fval = MRIFseq_vox(mri, x, y, z, frame);
                // if(x==10 && y == 0 && z == 0 && frame == 67)
                // printf("MRIIO: %g\n",fval);
                znzwriteFloat(fval, fp);
              }
              break;
            case MRI_UCHAR:
              if ((int)znzwrite(&MRIseq_vox(mri, 0, y, z, frame), sizeof(BUFTYPE), width, fp) != width) {
                errno = 0;
                ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "mghWrite: could not write %d bytes to %s", width, fname));
              }
              break;
            default:
              errno = 0;
              ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "mghWrite: unsupported type %d", mri->type));
              break;
	      } // switch
	} // height
        exec_progress_callback(z, depth, frame - start_frame, end_frame - start_frame + 1);
      } // depth
    } // frame
  } // end of writing voxel point by point

  if (getenv("FS_MGZIO_TIMING"))
  {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    printf("Total time (mghWrite) = %ld.%09ld seconds%s\n", 
           (end.tv_nsec < begin.tv_nsec) ? (end.tv_sec - 1 - begin.tv_sec) : (end.tv_sec - begin.tv_sec), 
           (end.tv_nsec < begin.tv_nsec) ? (1000000000 + end.tv_nsec - begin.tv_nsec) : (end.tv_nsec - begin.tv_nsec),
           (USEVOXELBUF) ? " (USEVOXELBUF)" : "");
  }

  znzwriteFloat(mri->tr, fp);
  znzwriteFloat(mri->flip_angle, fp);
  znzwriteFloat(mri->te, fp);
  znzwriteFloat(mri->ti, fp);
  znzwriteFloat(mri->fov, fp);

  // if mri->transform_fname has non-zero length
  // I write a tag with strlength and write it
  // I increase the tag_datasize with this amount
  if ((flen = strlen(mri->transform_fname)) > 0) {
    znzTAGwrite(fp, TAG_MGH_XFORM, mri->transform_fname, flen + 1);
  }
  // If we have any saved tag data, write it.
  if (NULL != mri->tag_data) {
    // Int is 32 bit on 32 bit and 64 bit os and thus it is safer
    znzwriteInt(mri->tag_data_size, fp);
    znzwrite(mri->tag_data, mri->tag_data_size, 1, fp);
  }

  if (mri->AutoAlign) znzWriteMatrix(fp, mri->AutoAlign, TAG_AUTO_ALIGN);
  if (mri->pedir)
    znzTAGwrite(fp, TAG_PEDIR, mri->pedir, strlen(mri->pedir) + 1);
  else
    znzTAGwrite(fp, TAG_PEDIR, (void *)"UNKNOWN", strlen("UNKNOWN"));
  if (mri->origRas2Vox)
  {
    printf("saving original ras2vox\n") ;
    znzWriteMatrix(fp, mri->origRas2Vox, TAG_ORIG_RAS2VOX);
  }

  znzTAGwrite(fp, TAG_FIELDSTRENGTH, (void *)(&mri->FieldStrength), sizeof(mri->FieldStrength));

  znzTAGwriteMRIframes(fp, mri);

  if (mri->ct) {
    znzwriteInt(TAG_OLD_COLORTABLE, fp);
    znzCTABwriteIntoBinary(mri->ct, fp);
  }

  // write other tags
  for (int i = 0; i < mri->ncmds; i++) znzTAGwrite(fp, TAG_CMDLINE, mri->cmdlines[i], strlen(mri->cmdlines[i]) + 1);

  // fclose(fp) ;
  znzclose(fp);

  return (NO_ERROR);
}

/*!
\fn MRI *MRIreorder4(MRI *mri, int order[4])
\brief Can reorders all 4 dimensions. Just copies old header to new.
\param mri - input
\param order[4] - new order of old dims
*/
MRI *MRIreorder4(MRI *mri, int order[4])
{
  MRI *result;
  int n, olddims[4], newdims[4];
  int dold[4], dnew[4];
  int c0, r0, s0, f0;
  int c1, r1, s1, f1;
  double v;

  olddims[0] = mri->width;
  olddims[1] = mri->height;
  olddims[2] = mri->depth;
  olddims[3] = mri->nframes;

  for (n = 0; n < 4; n++) newdims[n] = olddims[order[n] - 1];
  // for(n=0; n<4; n++)
  //  printf("%d %d  %d  %d\n",n,order[n],olddims[n],newdims[n]);
  result = MRIallocSequence(newdims[0], newdims[1], newdims[2], mri->type, newdims[3]);
  if (result == NULL) return (NULL);
  MRIcopyHeader(mri, result);

  for (c0 = 0; c0 < mri->width; c0++) {
    for (r0 = 0; r0 < mri->height; r0++) {
      for (s0 = 0; s0 < mri->depth; s0++) {
        for (f0 = 0; f0 < mri->nframes; f0++) {
          dold[0] = c0;
          dold[1] = r0;
          dold[2] = s0;
          dold[3] = f0;
          for (n = 0; n < 4; n++) dnew[n] = dold[order[n] - 1];
          c1 = dnew[0];
          r1 = dnew[1];
          s1 = dnew[2];
          f1 = dnew[3];
          v = MRIgetVoxVal(mri, c0, r0, s0, f0);
          MRIsetVoxVal(result, c1, r1, s1, f1, v);
        }
      }
    }
    exec_progress_callback(c0, mri->width, 0, 1);
  }

  return (result);
}

/*!
\fn MRIreorderVox2RAS(MRI *mri_src, MRI *mri_dst, int xdim, int ydim, int zdim)
\brief Changes the vox2ras matrix (and so xsize, x_r, ysize, etc) so that
it is appropriate for the given reordering of dimensions as performed by
MRIreorder(). xdim,ydim,zdim indicate which dimension the given dimension
moves to. So xdim=2 indicates that the cols of the source become the rows
of the destination. The dims can be signed.
*/
int MRIreorderVox2RAS(MRI *mri_src, MRI *mri_dst, int xdim, int ydim, int zdim)
{
  MATRIX *Q, *Msrc, *Mdst;

  // Source vox2ras
  Msrc = MRIxfmCRS2XYZ(mri_src, 0);

  // Q maps CRSdst to CRSsrc, ie, CRSsrc = Q*CRSdst
  Q = MatrixConstVal(0, 4, 4, NULL);
  Q->rptr[1][abs(xdim)] = ISIGN(xdim);
  Q->rptr[2][abs(ydim)] = ISIGN(ydim);
  Q->rptr[3][abs(zdim)] = ISIGN(zdim);
  Q->rptr[4][4] = 1;
  if (xdim < 0) Q->rptr[1][4] = mri_src->width - 1;
  if (ydim < 0) Q->rptr[2][4] = mri_src->height - 1;
  if (zdim < 0) Q->rptr[3][4] = mri_src->depth - 1;

  // Dst vox2ras = Msrc*Q
  Mdst = MatrixMultiply(Msrc, Q, NULL);
  if (Gdiag > 0) {
    printf("Msrc ------------------------\n");
    MatrixPrint(stdout, Msrc);
    printf("Q ------------------------\n");
    MatrixPrint(stdout, Q);
    printf("Qinv ------------------------\n");
    MatrixPrint(stdout, MatrixInverse(Q, NULL));
    printf("Mdst ------------------------\n");
    MatrixPrint(stdout, Mdst);
  }

  // Set destination vox2ras (does not change vox size!)
  MRIsetVox2RASFromMatrix(mri_dst, Mdst);

  // Free at last
  MatrixFree(&Msrc);
  MatrixFree(&Mdst);
  MatrixFree(&Q);
  return (0);
}

/*!
\fn MRIreorder(MRI *mri_src, MRI *mri_dst, int xdim, int ydim, int zdim)
\brief Reorders the dimensions of a volume.  xdim,ydim,zdim indicate
which dimension the given dimension moves to. So xdim=2 indicates that
the cols of the source become the rows of the destination. The dims
can be signed indicating a reversal. xdim, ydim, zdim values are only
+/- 1,2,3 only. The vox2ras matrix is updated so that the new and old
volumes will share RAS space (and can be alligned with a header
registration. Handles multiple frames.
*/
MRI *MRIreorder(MRI *mri_src, MRI *mri_dst, int xdim, int ydim, int zdim)
{
  int width, height, depth, xs, ys, zs, xd, yd, zd, x, y, z, f;
  int srcdims[3], dstdims[3];
  float srcsizes[3], dstsizes[3];

  if(Gdiag > 0){
    printf("MRIreorder() -----------\n");
    printf("xdim=%d ydim=%d zdim=%d\n", xdim, ydim, zdim);
  }
  /* check that the source ras coordinates are good and
     that each direction is used once and only once */
  if (abs(xdim) * abs(ydim) * abs(zdim) != 6 || abs(xdim) + abs(ydim) + abs(zdim) != 6) {
    printf("ERROR: replicated/incorrect dimension number\n");
    return (NULL);
  }

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  srcdims[0] = width;
  srcdims[1] = height;
  srcdims[2] = depth;
  srcsizes[0] = mri_src->xsize;
  srcsizes[1] = mri_src->ysize;
  srcsizes[2] = mri_src->zsize;

  dstdims[abs(xdim) - 1] = srcdims[0];
  dstdims[abs(ydim) - 1] = srcdims[1];
  dstdims[abs(zdim) - 1] = srcdims[2];
  dstsizes[abs(xdim) - 1] = srcsizes[0];
  dstsizes[abs(ydim) - 1] = srcsizes[1];
  dstsizes[abs(zdim) - 1] = srcsizes[2];

  if (!mri_dst) {
    // mri_dst = MRIclone(mri_src, NULL) ;
    mri_dst = MRIallocSequence(dstdims[0], dstdims[1], dstdims[2], mri_src->type, mri_src->nframes);
    MRIcopyHeader(mri_src, mri_dst);
    mri_dst->xsize = dstsizes[0];
    mri_dst->ysize = dstsizes[1];
    mri_dst->zsize = dstsizes[2];
  }
  MRIreorderVox2RAS(mri_src, mri_dst, xdim, ydim, zdim);

  if(Gdiag > 0){
    printf("src %d %d %d, %f %f %f\n",
	   mri_src->width,
	   mri_src->height,
	   mri_src->depth,
	   mri_src->xsize,
	   mri_src->ysize,
	   mri_src->zsize);
    printf("dst %d %d %d, %f %f %f\n",
	   mri_dst->width,
	   mri_dst->height,
	   mri_dst->depth,
	   mri_dst->xsize,
	   mri_dst->ysize,
	   mri_dst->zsize);
  }

  xd = yd = zd = 0;

  // XDIM=1, YDIM=2, ZDIM=3
  // xdim tells original x-axis/cols goes to where
  switch (abs(xdim)) {
    default:
    case XDIM:
      if (mri_dst->width != mri_src->width) {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "MRIreorder: incorrect dst width"));
      }
      break;
    case YDIM:
      if (mri_dst->height != mri_src->width) {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "MRIreorder: incorrect dst width"));
      }
      break;
    case ZDIM:
      if (mri_dst->depth != mri_src->width) {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "MRIreorder: incorrect dst width"));
      }
      break;
  }
  // ydim tells original y-axis/rows goes to where
  switch (abs(ydim)) {
    default:
    case XDIM:  // map srcrows to dstcols
      if (mri_src->height != mri_dst->width) {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "MRIreorder: incorrect dst width != src height"));
      }
      break;
    case YDIM:  // map srcrows to dstrows
      if (mri_src->height != mri_dst->height) {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "MRIreorder: incorrect dst height != src height"));
      }
      break;
    case ZDIM:  // map srcrows to dstrows
      if (mri_src->height != mri_dst->depth) {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "MRIreorder: incorrect dst depth != height"));
      }
      break;
  }
  // zdim tells original z-axis/slices goes to where
  switch (abs(zdim)) {
    default:
    case XDIM:
      if (mri_src->depth != mri_dst->width) {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "MRIreorder: incorrect dst depth"));
      }
      break;
    case YDIM:
      if (mri_src->depth != mri_dst->height) {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "MRIreorder: incorrect dst depth"));
      }
      break;
    case ZDIM:
      if (mri_src->depth != mri_dst->depth) {
        errno = 0;
        ErrorReturn(NULL, (ERROR_BADPARM, "MRIreorder: incorrect dst depth"));
      }
      break;
  }

  for (f = 0; f < mri_src->nframes; f++) {
    for (zs = 0; zs < depth; zs++) {
      if (zdim < 0)
        z = depth - zs - 1;
      else
        z = zs;
      switch (abs(zdim)) {
        case XDIM:
          xd = z;
          break;
        case YDIM:
          yd = z;
          break;
        case ZDIM:
        default:
          zd = z;
          break;
      }
      for (ys = 0; ys < height; ys++) {
        if (ydim < 0)
          y = height - ys - 1;
        else
          y = ys;
        switch (abs(ydim)) {
          case XDIM:
            xd = y;
            break;
          case YDIM:
            yd = y;
            break;
          case ZDIM:
          default:
            zd = y;
            break;
        }
        for (xs = 0; xs < width; xs++) {
          if (xdim < 0)
            x = width - xs - 1;
          else
            x = xs;
          switch (abs(xdim)) {
            case XDIM:
              xd = x;
              break;
            case YDIM:
              yd = x;
              break;
            case ZDIM:
            default:
              zd = x;
              break;
          }
          switch (mri_src->type) {
            case MRI_SHORT:
              MRISseq_vox(mri_dst, xd, yd, zd, f) = MRISseq_vox(mri_src, xs, ys, zs, f);
              break;
            case MRI_USHRT:
              MRIUSseq_vox(mri_dst, xd, yd, zd, f) = MRIUSseq_vox(mri_src, xs, ys, zs, f);
              break;
            case MRI_FLOAT:
              MRIFseq_vox(mri_dst, xd, yd, zd, f) = MRIFseq_vox(mri_src, xs, ys, zs, f);
              break;
            case MRI_INT:
              MRIIseq_vox(mri_dst, xd, yd, zd, f) = MRIIseq_vox(mri_src, xs, ys, zs, f);
              break;
            case MRI_UCHAR:
              MRIseq_vox(mri_dst, xd, yd, zd, f) = MRIseq_vox(mri_src, xs, ys, zs, f);
              break;
            default:
              errno = 0;
              ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIreorder: unsupported voxel format %d", mri_src->type));
              break;
          }
        }
      }

      exec_progress_callback(zs, depth, f, mri_src->nframes);
    }
  } /* end loop over frames */
  return (mri_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Write the MRI header information to the file
  COR-.info in the directory specified by 'fpref'
  ------------------------------------------------------*/
int MRIwriteInfo(MRI *mri, const char *fpref)
{
  FILE *fp;
  char fname[STRLEN];
  int slice_direction;
  sprintf(fname, "%s/%s", fpref, INFO_FNAME);
  fp = fopen(fname, "w");
  if (fp == NULL) {
    errno = 0;
    ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MRIwriteInfo(%s): could not open %s.\n", fpref, fname));
  }

  fprintf(fp, "%s %d\n", "imnr0", mri->imnr0);
  fprintf(fp, "%s %d\n", "imnr1", mri->imnr1);
  slice_direction = getSliceDirection(mri);
  fprintf(fp, "%s %d\n", "ptype", slice_direction == MRI_CORONAL ? 2 : slice_direction == MRI_HORIZONTAL ? 0 : 1);
  fprintf(fp, "%s %d\n", "x", mri->width);
  fprintf(fp, "%s %d\n", "y", mri->height);
  fprintf(fp, "%s %f\n", "fov", mri->fov / MM_PER_METER);
  fprintf(fp, "%s %f\n", "thick", mri->ps / MM_PER_METER);
  fprintf(fp, "%s %f\n", "psiz", mri->ps / MM_PER_METER);
  fprintf(fp, "%s %f\n", "locatn", mri->location);             /* locatn */
  fprintf(fp, "%s %f\n", "strtx", mri->xstart / MM_PER_METER); /* strtx */
  fprintf(fp, "%s %f\n", "endx", mri->xend / MM_PER_METER);    /* endx */
  fprintf(fp, "%s %f\n", "strty", mri->ystart / MM_PER_METER); /* strty */
  fprintf(fp, "%s %f\n", "endy", mri->yend / MM_PER_METER);    /* endy */
  fprintf(fp, "%s %f\n", "strtz", mri->zstart / MM_PER_METER); /* strtz */
  fprintf(fp, "%s %f\n", "endz", mri->zend / MM_PER_METER);    /* endz */
  fprintf(fp, "%s %f\n", "tr", mri->tr);
  fprintf(fp, "%s %f\n", "te", mri->te);
  fprintf(fp, "%s %f\n", "ti", mri->ti);
  if (mri->linear_transform) {
    char fname[STRLEN];

    /*
       this won't work for relative paths which are not the same for the
       destination directory as for the the source directory.
    */
    sprintf(fname, "%s", mri->transform_fname);
    fprintf(fp, "xform %s\n", fname);

  }

  fprintf(fp, "%s %d\n", "ras_good_flag", mri->ras_good_flag);
  fprintf(fp, "%s %f %f %f\n", "x_ras", mri->x_r, mri->x_a, mri->x_s);
  fprintf(fp, "%s %f %f %f\n", "y_ras", mri->y_r, mri->y_a, mri->y_s);
  fprintf(fp, "%s %f %f %f\n", "z_ras", mri->z_r, mri->z_a, mri->z_s);
  fprintf(fp, "%s %f %f %f\n", "c_ras", mri->c_r, mri->c_a, mri->c_s);

  fclose(fp);

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Write an MRI header and a set of data files to
  the directory specified by 'fpref'
  ------------------------------------------------------*/
int MRIappend(MRI *mri, const char *fpref)
{
  int type, frame;
  char fname[STRLEN];

  MRIunpackFileName(fpref, &frame, &type, fname);
  if (type == MRI_MGH_FILE)
    return (mghAppend(mri, fname, frame));
  else {
    errno = 0;
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "MRIappend(%s): file type not supported", fname));
  }

  return (NO_ERROR);
}

static int mghAppend(MRI *mri, const char *fname, int frame)
{
  FILE *fp;
  int start_frame, end_frame, x, y, z, width, height, depth, nframes;

  if (frame >= 0)
    start_frame = end_frame = frame;
  else {
    start_frame = 0;
    end_frame = mri->nframes - 1;
  }
  fp = fopen(fname, "rb");
  if (!fp) /* doesn't exist */
    return (mghWrite(mri, fname, frame));
  fclose(fp);
  fp = fopen(fname, "r+b");
  if (!fp) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "mghAppend(%s, %d): could not open file", fname, frame));
  }

  /* WARNING - this is dependent on the order of writing in mghWrite */
  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  fseek(fp, 4 * sizeof(int), SEEK_SET);
  nframes = freadInt(fp);
  fseek(fp, 4 * sizeof(int), SEEK_SET);
  fwriteInt(nframes + end_frame - start_frame + 1, fp);
  fseek(fp, 0, SEEK_END);

  for (frame = start_frame; frame <= end_frame; frame++) {
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        switch (mri->type) {
          case MRI_FLOAT:
            for (x = 0; x < width; x++) {
              fwriteFloat(MRIFseq_vox(mri, x, y, z, frame), fp);
            }
            break;
          case MRI_UCHAR:
            if ((int)fwrite(&MRIseq_vox(mri, 0, y, z, frame), sizeof(BUFTYPE), width, fp) != width) {
              errno = 0;
              ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "mghAppend: could not write %d bytes to %s", width, fname));
            }
            break;
          default:
            errno = 0;
            ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "mghAppend: unsupported type %d", mri->type));
            break;
        }
      }
    }
  }

  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRIunpackFileName(const char *inFname, int *pframe, int *ptype, char *outFname)
{
  char *number = NULL, *at = NULL, buf[STRLEN];
  struct stat stat_buf;

  strcpy(outFname, inFname);
  if (MRIIO_Strip_Pound)
    number = strrchr(outFname, '#');
  else
    number = NULL;

  at = strrchr(outFname, '@');

  if (at) *at = '\0';

  if (number) /* '#' in filename indicates frame # */
  {
    if (sscanf(number + 1, "%d", pframe) < 1) *pframe = -1;
    *number = 0;
  }
  else
    *pframe = -1;

  if (at) {
    at = StrUpper(strcpy(buf, at + 1));
    if (!strcmp(at, "MNC"))
      *ptype = MRI_MINC_FILE;
    else if (!strcmp(at, "MINC"))
      *ptype = MRI_MINC_FILE;
    else if (!strcmp(at, "BRIK"))
      *ptype = BRIK_FILE;
    else if (!strcmp(at, "SIEMENS"))
      *ptype = SIEMENS_FILE;
    else if (!strcmp(at, "MGH"))
      *ptype = MRI_MGH_FILE;
    else if (!strcmp(at, "MR"))
      *ptype = GENESIS_FILE;
    else if (!strcmp(at, "GE"))
      *ptype = GE_LX_FILE;
    else if (!strcmp(at, "IMG"))
      *ptype = MRI_ANALYZE_FILE;
    else if (!strcmp(at, "COR"))
      *ptype = MRI_CORONAL_SLICE_DIRECTORY;
    else if (!strcmp(at, "BSHORT"))
      *ptype = BSHORT_FILE;
    else if (!strcmp(at, "SDT"))
      *ptype = SDT_FILE;
    else {
      errno = 0;
      ErrorExit(ERROR_UNSUPPORTED, "unknown file type %s", at);
    }
  }
  else /* no '@' found */
  {
    *ptype = -1;

    if (is_genesis(outFname))
      *ptype = GENESIS_FILE;
    else if (is_ge_lx(outFname))
      *ptype = GE_LX_FILE;
    else if (is_brik(outFname))
      *ptype = BRIK_FILE;
    else if (is_siemens(outFname))
      *ptype = SIEMENS_FILE;
    else if (is_analyze(outFname))
      *ptype = MRI_ANALYZE_FILE;
    else if (is_signa(outFname))
      *ptype = SIGNA_FILE;
    else if (is_sdt(outFname))
      *ptype = SDT_FILE;
    else if (is_mgh(outFname))
      *ptype = MRI_MGH_FILE;
    else if (is_mnc(outFname))
      *ptype = MRI_MINC_FILE;
    else if (is_bshort(outFname))
      *ptype = BSHORT_FILE;
    else {
      if (stat(outFname, &stat_buf) < 0) {
        errno = 0;
        ErrorReturn(ERROR_BADFILE, (ERROR_BAD_FILE, "can't stat file %s", outFname));
      }
      if (S_ISDIR(stat_buf.st_mode)) *ptype = MRI_CORONAL_SLICE_DIRECTORY;
    }

    if (*ptype == -1) {
      errno = 0;
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "unrecognized file type for file %s", outFname));
    }
  }

  return (NO_ERROR);
}

/*---------------------------------------------------------------
  MRIwriteAnyFormat() - saves the data in the given mri structure to
  "any" format, where "any" is defined as anything writable by
  MRIwrite and wfile format. If wfile format is used, then mriframe
  must be a valid frame number. The val field of the surf is
  preserved. If another format is used, then, if mriframe = -1, then
  all frames are stored, otherwise the given frame is stored. If
  fmt is NULL, then it will attempt to infer the format from the
  name. surf only needs to be supplied when format is wfile. Legal
  formats include: wfile, paint, w, bshort, bfloat, COR, analyze,
  analyze4d, spm.
  ---------------------------------------------------------------*/
int MRIwriteAnyFormat(MRI *mri, const char *fileid, const char *fmt, int mriframe, MRIS *surf)
{
  int fmtid, err, n, r, c, s;
  float *v = NULL, f;
  MRI *mritmp = NULL;

  if (fmt != NULL && (!strcmp(fmt, "paint") || !strcmp(fmt, "w") || !strcmp(fmt, "wfile"))) {
    /* Save as a wfile */
    if (surf == NULL) {
      printf("ERROR: MRIwriteAnyFormat: need surf with paint format\n");
      return (1);
    }
    if (mriframe >= mri->nframes) {
      printf(
          "ERROR: MRIwriteAnyFormat: frame (%d) exceeds number of\n"
          "       frames\n",
          mriframe >= mri->nframes);
      return (1);
    }

    /* Copy current surf values into v (temp storage) */
    v = (float *)calloc(surf->nvertices, sizeof(float));
    for (n = 0; n < surf->nvertices; n++) v[n] = surf->vertices[n].val;

    /* Copy the mri values into the surf values */
    err = MRIScopyMRI(surf, mri, mriframe, "val");
    if (err) {
      printf("ERROR: MRIwriteAnyFormat: could not copy MRI to MRIS\n");
      return (1);
    }

    /* Write the surf values */
    err = MRISwriteValues(surf, fileid);

    /* Copy v back into surf values */
    for (n = 0; n < surf->nvertices; n++) surf->vertices[n].val = v[n];
    free(v);

    /* Check for errors from write of values */
    if (err) {
      printf("ERROR: MRIwriteAnyFormat: MRISwriteValues\n");
      return (1);
    }

    return (0);
  }

  /* Copy the desired frame if necessary */
  if (mriframe > -1) {
    if (mriframe >= mri->nframes) {
      printf(
          "ERROR: MRIwriteAnyFormat: frame (%d) exceeds number of\n"
          "       frames\n",
          mriframe >= mri->nframes);
      return (1);
    }
    mritmp = MRIallocSequence(mri->width, mri->height, mri->depth, mri->type, 1);
    for (c = 0; c < mri->width; c++) {
      for (r = 0; r < mri->height; r++) {
        for (s = 0; s < mri->depth; s++) {
          f = MRIgetVoxVal(mri, c, r, s, mriframe);
          MRIsetVoxVal(mritmp, c, r, s, 0, f);
        }
      }
    }
  }
  else
    mritmp = mri;

  /*------------ Save using MRIwrite or MRIwriteType ---------*/
  if (fmt != NULL) {
    /* Save as the given format */
    fmtid = string_to_type(fmt);
    if (fmtid == MRI_VOLUME_TYPE_UNKNOWN) {
      printf("ERROR: format string %s unrecognized\n", fmt);
      return (1);
    }
    err = MRIwriteType(mritmp, fileid, fmtid);
    if (err) {
      printf("ERROR: MRIwriteAnyFormat: could not write to %s\n", fileid);
      return (1);
    }
  }
  else {
    /* Try to infer the type and save (format is NULL) */
    err = MRIwrite(mritmp, fileid);
    if (err) {
      printf("ERROR: MRIwriteAnyFormat: could not write to %s\n", fileid);
      return (1);
    }
  }

  if (mri != mritmp) MRIfree(&mritmp);

  return (0);
}

/* EOF */

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* ----------- Obsolete functions below. --------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

#include "gca.h"
static MRI *readGCA(const char *fname, int start_frame, int end_frame)
{
  GCA *gca;
  MRI *mri;

  gca = GCAread(fname);
  if (!gca) return (NULL);
  printf("reading frame %d of gca\n", start_frame);
  switch (start_frame) {
    default:
    case -1: {
      MRI *mri_means, *mri_labels, *mri_pvals, *mri_pvals_resampled;

      mri = MRIallocSequence(gca->width, gca->height, gca->depth, MRI_FLOAT, 3);
      mri_means = MRIallocSequence(gca->width, gca->height, gca->depth, MRI_FLOAT, 1);
      mri_labels = MRIallocSequence(gca->width, gca->height, gca->depth, MRI_FLOAT, 1);

      MRIsetResolution(mri, gca->xsize, gca->ysize, gca->zsize);
      MRIsetResolution(mri_labels, gca->xsize, gca->ysize, gca->zsize);
      MRIsetResolution(mri_means, gca->xsize, gca->ysize, gca->zsize);
      GCAcopyDCToMRI(gca, mri);
      GCAcopyDCToMRI(gca, mri_means);
      GCAcopyDCToMRI(gca, mri_labels);
      GCAbuildMostLikelyVolume(gca, mri_means);
      GCAbuildMostLikelyLabelVolume(gca, mri_labels);
      mri_pvals = GCAbuildMostLikelyLabelProbabilityVolume(gca);

      MRIcopyFrames(mri_means, mri, 0, 0, 0);
      MRIcopyFrames(mri_labels, mri, 0, 0, 1);

      mri_pvals_resampled = MRIresample(mri_pvals, mri_means, SAMPLE_NEAREST);
      MRIcopyFrames(mri_pvals_resampled, mri, 0, 0, 2);
      MRIfree(&mri_means);
      MRIfree(&mri_labels);
      MRIfree(&mri_pvals);
      MRIfree(&mri_pvals_resampled);
      break;
    }
    case 0:
      mri = MRIallocSequence(gca->node_width, gca->node_height, gca->node_depth, MRI_FLOAT, 1);
      MRIsetResolution(
          mri, gca->xsize * gca->node_spacing, gca->ysize * gca->node_spacing, gca->zsize * gca->node_spacing);
      GCAcopyDCToMRI(gca, mri);
      GCAbuildMostLikelyVolume(gca, mri);
      break;
    case 1:
      mri = MRIallocSequence(gca->prior_width, gca->prior_height, gca->prior_depth, MRI_SHORT, 1);
      MRIsetResolution(
          mri, gca->xsize * gca->prior_spacing, gca->ysize * gca->prior_spacing, gca->zsize * gca->prior_spacing);
      GCAcopyDCToMRI(gca, mri);
      GCAbuildMostLikelyLabelVolume(gca, mri);
      break;
    case 2:
      printf("interpreting as probability volume\n");
      mri = GCAbuildMostLikelyLabelProbabilityVolume(gca);
      break;
  }
  GCAfree(&gca);
  return (mri);
}

MRI *MRIremoveNaNs(MRI *mri_src, MRI * mri_dst)
{
  if (mri_dst != mri_src) mri_dst = MRIcopy(mri_src, mri_dst);

  int x;
  int nans = 0;
  static int first = 1;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(shown_reproducible) shared(mri_dst) reduction(+ : nans)
#endif
  for (x = 0; x < mri_dst->width; x++) {
    ROMP_PFLB_begin
    
    int const height  = mri_dst->height;
    int const depth   = mri_dst->depth;
    int const nframes = mri_dst->nframes;

    int y, z, f;
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        for (f = 0; f < nframes; f++) {
          float val = MRIgetVoxVal(mri_dst, x, y, z, f);
          if (!std::isfinite(val)) {
            nans++;
	    if (getenv("FS_LEAVE_NANS") == NULL)
	      MRIsetVoxVal(mri_dst, x, y, z, f, 0);
	    if (first)
	    {
  printf("NaN found at voxel (%d, %d, %d, %d)\n", x, y, z, f);
  first = 0;
}

          }
        }
      }
    }
    exec_progress_callback(x, mri_dst->width, 0, 1);
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

  if (nans > 0) ErrorPrintf(ERROR_BADPARM, "WARNING: %d NaNs found in volume %s...\n", nans, mri_src->fname);
  return (mri_dst);
}

int MRIaddCommandLine(MRI *mri, const std::string& cmdline)
{
  if (mri->ncmds >= MAX_CMDS)
    fs::error() << "can't add cmd to mri since max cmds (" << mri->ncmds <<  ") has been reached";

  int i = mri->ncmds++;
  mri->cmdlines[i] = (char *)calloc(cmdline.size() + 1, sizeof(char));
  strcpy(mri->cmdlines[i], cmdline.c_str());
  return NO_ERROR;
}

/*------------------------------------------------------------------
  niiPrintHdr() - this dumps (most of) the nifti header to the given
  stream.
  ------------------------------------------------------------------*/
static int niiPrintHdr(FILE *fp, struct nifti_1_header *hdr)
{
  int n;
  fprintf(fp, "sizeof_hdr %d \n", hdr->sizeof_hdr);
  fprintf(fp, "data_type  %s \n", hdr->data_type);
  fprintf(fp, "db_name    %s \n", hdr->db_name);
  fprintf(fp, "extents    %d \n", hdr->extents);
  fprintf(fp, "session_error %d \n", hdr->session_error);
  fprintf(fp, "regular    %c \n", hdr->regular);
  fprintf(fp, "dim_info   %c \n", hdr->dim_info);
  for (n = 0; n < 8; n++) fprintf(fp, "dim[%d] %d\n", n, hdr->dim[n]);
  fprintf(fp, "intent_p1  %f \n", hdr->intent_p1);
  fprintf(fp, "intent_p2  %f \n", hdr->intent_p2);
  fprintf(fp, "intent_p3  %f \n", hdr->intent_p3);
  fprintf(fp, "intent_code  %d \n", hdr->intent_code);
  fprintf(fp, "datatype      %d \n", hdr->datatype);
  fprintf(fp, "bitpix        %d \n", hdr->bitpix);
  fprintf(fp, "slice_start   %d \n", hdr->slice_start);
  for (n = 0; n < 8; n++) fprintf(fp, "pixdim[%d] %f\n", n, hdr->pixdim[n]);
  fprintf(fp, "vox_offset    %f\n", hdr->vox_offset);
  fprintf(fp, "scl_slope     %f\n", hdr->scl_slope);
  fprintf(fp, "scl_inter     %f\n", hdr->scl_inter);
  fprintf(fp, "slice_end     %d\n", hdr->slice_end);
  fprintf(fp, "slice_code    %c\n", hdr->slice_code);
  fprintf(fp, "xyzt_units    %c\n", hdr->xyzt_units);
  fprintf(fp, "cal_max       %f\n", hdr->cal_max);
  fprintf(fp, "cal_min       %f\n", hdr->cal_min);
  fprintf(fp, "slice_duration %f\n", hdr->slice_duration);
  fprintf(fp, "toffset        %f \n", hdr->toffset);
  fprintf(fp, "glmax          %d \n", hdr->glmax);
  fprintf(fp, "glmin          %d \n", hdr->glmin);
  fprintf(fp, "descrip        %s \n", hdr->descrip);
  fprintf(fp, "aux_file       %s \n", hdr->aux_file);
  fprintf(fp, "qform_code     %d \n", hdr->qform_code);
  fprintf(fp, "sform_code     %d \n", hdr->sform_code);
  // There's more, but I ran out of steam ...
  // fprintf(fp,"          %d \n",hdr->);
  return (0);
}
