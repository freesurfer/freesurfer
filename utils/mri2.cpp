/**
 * @brief more routines for loading, saving, and operating on MRI structures
 *
 */
/*
 * Original Author: Douglas N. Greve
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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "bfileio.h"
#include "cma.h"
#include "corio.h"
#include "diag.h"
#include "error.h"
#include "fio.h"
#include "fmriutils.h"
#include "mri.h"
#include "mriBSpline.h"
#include "mrisurf.h"
#include "proto.h"
#include "region.h"
#include "sig.h"
#include "stats.h"
#include "mrimorph.h"
#include "mri_identify.h"
#include "mri2.h"

//#define MRI2_TIMERS

#include "affine.h"
#include "romp_support.h"

#ifndef FSIGN
# define FSIGN(f) (((f) < 0) ? -1 : 1)
#endif

/* overwrite generic nint to speed up execution
  Make sure NOT use calls like nint(f++) in this file!
*/
#define nint(f) (f < 0 ? ((int)(f - 0.5)) : ((int)(f + 0.5)))

#define VERBOSE_MODE
#undef VERBOSE_MODE

/*-------------------------------------------------------------
  mri_load_bvolume() -- same as bf_ldvolume() but returns an
  MRI structure instead of a BF_DATA. See bfileio.h.
  -------------------------------------------------------------*/
MRI *mri_load_bvolume(char *bfstem)
{
  BF_DATA *bfvol;
  MRI *vol;
  int r, c, s, f;
  float val;

  /* first load as a BF_DATA sturcture */
  bfvol = bf_ldvolume(bfstem);
  if (bfvol == NULL) return (NULL);

  /* allocate the MRI */
  vol = MRIallocSequence(bfvol->ncols, bfvol->nrows, bfvol->nslcs, MRI_FLOAT, bfvol->nfrms);
  if (vol == NULL) {
    bf_freebfd(&bfvol);
    fprintf(stderr, "mri_load_bvolume(): could not alloc vol\n");
    return (NULL);
  }

  /* copy data from the BF_DATA struct to the ARRAY4D*/
  for (r = 0; r < bfvol->nrows; r++) {
    for (c = 0; c < bfvol->ncols; c++) {
      for (s = 0; s < bfvol->nslcs; s++) {
        for (f = 0; f < bfvol->nfrms; f++) {
          val = BF_GETVAL(bfvol, r, c, s, f);
          MRIFseq_vox(vol, c, r, s, f) = val;
        }
      }
    }
  }

  bf_freebfd(&bfvol);
  return (vol);
}

/*-------------------------------------------------------------
  mri_save_as_bvolume() - same as bf_svvolume() but takes an MRI
  sturucture as input. See also bfileio.h.
  -------------------------------------------------------------*/
int mri_save_as_bvolume(MRI *vol, char *stem, int svendian, int svtype)
{
  BF_DATA *bfvol;
  int r, c, s, f;
  float val;

  /* allocate a temporary BF_DATA struct */
  bfvol = bf_allocbfd(vol->height, vol->width, vol->depth, vol->nframes);
  if (bfvol == NULL) return (1);

  /* copy data from ARRAY4D to BF_DATA */
  for (r = 0; r < bfvol->nrows; r++) {
    for (c = 0; c < bfvol->ncols; c++) {
      for (s = 0; s < bfvol->nslcs; s++) {
        for (f = 0; f < bfvol->nfrms; f++) {
          val = MRIFseq_vox(vol, c, r, s, f);
          BF_SETVAL(val, bfvol, r, c, s, f);
        }
      }
    }
  }

  /* save the BF_DATA volume */
  bf_svvolume(bfvol, stem, svendian, svtype);

  bf_freebfd(&bfvol);
  return (0);
}
/*-------------------------------------------------------------
  mri_load_bvolume_frame() -- loads a single frame as an MRI.
  -------------------------------------------------------------*/
MRI *mri_load_bvolume_frame(char *bfstem, int frameno)
{
  BF_DATA *bfvol;
  MRI *vol;
  int r, c, s;
  float val;

  /* first load as a BF_DATA sturcture */
  bfvol = bf_ldvolume(bfstem);
  if (bfvol == NULL) return (NULL);

  if (frameno >= bfvol->nfrms) {
    fprintf(stderr,
            "ERROR: mri_load_bvolume_frame(): frameno = %d, exceeds "
            "number of frames in %s = %d\n",
            frameno,
            bfstem,
            bfvol->nfrms);
    bf_freebfd(&bfvol);
    return (NULL);
  }

  /* allocate the MRI */
  vol = MRIallocSequence(bfvol->ncols, bfvol->nrows, bfvol->nslcs, MRI_FLOAT, 1);
  if (vol == NULL) {
    bf_freebfd(&bfvol);
    fprintf(stderr, "mri_load_bvolume_frame(): could not alloc vol\n");
    return (NULL);
  }

  /* copy data from the BF_DATA struct to the ARRAY4D*/
  for (r = 0; r < bfvol->nrows; r++) {
    for (c = 0; c < bfvol->ncols; c++) {
      for (s = 0; s < bfvol->nslcs; s++) {
        val = BF_GETVAL(bfvol, r, c, s, frameno);
        MRIFseq_vox(vol, c, r, s, 0) = val;
      }
    }
  }

  bf_freebfd(&bfvol);
  return (vol);
}
/*-------------------------------------------------------------
  mri_save_as_cor() - basically the same as sv_cor() (corio.c) but
  takes an MRI structure as input.  The MRi has much more
  flexibility than the COR structure which is limited to 256^3, one frame,
  and type uchar.  If any MRI dimension exceeds 256, it is truncated;
  if it is less than 256, the extra is filled with zeros. If the MRI
  has more than one frame, the frame can be set by the frame argument.
  If the values in the MRI are beyond the range of 0-255 (ie, that
  of unsigned char), the MRI can be rescaled to 0-255 by setting
  the rescale flag to non-zero. See also corio.h.
  -------------------------------------------------------------*/
int mri_save_as_cor(MRI *vol, char *cordir, int frame, int rescale)
{
  unsigned char **COR;
  int r, c, s;
  int rmax, cmax, smax;
  float val;

  if (frame >= vol->nframes) {
    fprintf(stderr, "mri_save_as_cor(): frame = %d, must be <= %d\n", frame, vol->nframes);
    return (1);
  }

  /* make sure the output directory is writable */
  if (!cordir_iswritable(cordir)) return (1);

  /* allocate a temporary COR volume */
  COR = alloc_cor();
  if (COR == NULL) return (1);

  /* rescale to 0-255 (range of uchar) */
  if (rescale) mri_rescale(vol, 0, 255, vol);

  /* make sure maximum subscript does not
     exceed either range */
  if (vol->height < 256)
    rmax = vol->height;
  else
    rmax = 256;
  if (vol->width < 256)
    cmax = vol->width;
  else
    cmax = 256;
  if (vol->depth < 256)
    smax = vol->depth;
  else
    smax = 256;

  /* copy data from ARRAY4D to COR */
  for (r = 0; r < rmax; r++) {
    for (c = 0; c < cmax; c++) {
      for (s = 0; s < smax; s++) {
        val = MRIFseq_vox(vol, c, r, s, 0);
        CORVAL(COR, r, c, s) = (unsigned char)val;
      }
    }
  }

  /* save the COR volume */
  sv_cor(COR, cordir);

  free_cor(&COR);
  return (0);
}
/*------------------------------------------------------------
  mri_rescale() -- rescales array to be min <= val <= max. Uses
  outvol if non-null, otherwise allocates an output volume. Can
  be done in-place. Rescales across all frames.
  ------------------------------------------------------------*/
MRI *mri_rescale(MRI *vol, float min, float max, MRI *outvol)
{
  int r, c, s, f;
  float val, volmin, volmax, range;
  // float volrange;
  MRI *tmpvol;

  if (outvol != NULL)
    tmpvol = outvol;
  else {
    tmpvol = MRIallocSequence(vol->width, vol->height, vol->depth, MRI_FLOAT, vol->nframes);
    if (tmpvol == NULL) return (NULL);
  }

  /* find the minimum and maximum */
  volmin = MRIFseq_vox(vol, 0, 0, 0, 0);
  volmax = MRIFseq_vox(vol, 0, 0, 0, 0);
  for (r = 0; r < vol->height; r++) {
    for (c = 0; c < vol->width; c++) {
      for (s = 0; s < vol->depth; s++) {
        for (f = 0; f < vol->nframes; f++) {
          val = MRIFseq_vox(vol, c, r, s, f);
          if (volmin > val) volmin = val;
          if (volmax < val) volmax = val;
        }
      }
    }
  }
  // volrange = volmax - volmin;
  range = max - min;

  /* rescale to the new range */
  for (r = 0; r < vol->height; r++) {
    for (c = 0; c < vol->width; c++) {
      for (s = 0; s < vol->depth; s++) {
        for (f = 0; f < vol->nframes; f++) {
          val = MRIFseq_vox(vol, c, r, s, f);
          val = range * val + min;
          MRIFseq_vox(tmpvol, c, r, s, f) = val;
        }
      }
    }
  }

  printf("volmin = %g, volmax = %g\n", volmin, volmax);

  return (tmpvol);
}
/*------------------------------------------------------------
  mri_minmax() -- gets min and max values of volume
  ------------------------------------------------------------*/
int mri_minmax(MRI *vol, float *min, float *max)
{
  int r, c, s, f;
  float val;

  *min = MRIFseq_vox(vol, 0, 0, 0, 0);
  *max = MRIFseq_vox(vol, 0, 0, 0, 0);
  for (r = 0; r < vol->height; r++) {
    for (c = 0; c < vol->width; c++) {
      for (s = 0; s < vol->depth; s++) {
        for (f = 0; f < vol->nframes; f++) {
          val = MRIgetVoxVal(vol, c, r, s, f);  // MRIFseq_vox(vol,c,r,s,f);
          if (*min > val) *min = val;
          if (*max < val) *max = val;
        }
      }
    }
  }
  // printf("volmin = %g, volmax = %g\n",*min,*max);
  return (0);
}
/*--------------------------------------------------------
  mri_framepower() -- raises the value in each frame to that
  indicated by the framepower vector.  This is a pre-processing
  step for averaging across voxels when the data contain both
  averages and standard deviations. The stddevs need to be
  squared before spatial averaging and then square-rooted before
  saving again.
  -----------------------------------------------------------*/
int mri_framepower(MRI *vol, float *framepower)
{
  int r, c, s, f;
  float val;

  for (f = 0; f < vol->nframes; f++) {
    /* power = 1 -- don't do anything */
    if (fabs(framepower[f] - 1.0) < .00001) continue;

    /* power = 0.5 -- use sqrt() */
    if (fabs(framepower[f] - 0.5) < .00001) {
      for (r = 0; r < vol->height; r++) {
        for (c = 0; c < vol->width; c++) {
          for (s = 0; s < vol->depth; s++) {
            val = MRIFseq_vox(vol, c, r, s, f);
            MRIFseq_vox(vol, c, r, s, f) = sqrt(val);
          }
        }
      }
      continue;
    }

    /* power = 2 -- use val*val */
    if (fabs(framepower[f] - 0.5) < .00001) {
      for (r = 0; r < vol->height; r++) {
        for (c = 0; c < vol->width; c++) {
          for (s = 0; s < vol->depth; s++) {
            val = MRIFseq_vox(vol, c, r, s, f);
            MRIFseq_vox(vol, c, r, s, f) = val * val;
          }
        }
      }
      continue;
    }

    /* generic: use pow() -- least efficient */
    for (r = 0; r < vol->height; r++) {
      for (c = 0; c < vol->width; c++) {
        for (s = 0; s < vol->depth; s++) {
          val = MRIFseq_vox(vol, c, r, s, f);
          MRIFseq_vox(vol, c, r, s, f) = pow(val, framepower[f]);
        }
      }
    }

  } /* end loop over frames */

  return (0);
}
/*-------------------------------------------------------------------
  mri_binarize() - converts each element to either 0 or 1 depending
  upon whether the value of the element is above or below the
  threshold. If the invert flag is set to 1, then the binarization is
  inverted. Can be done in-place (ie, vol=volbin). If volbin is NULL,
  then a new MRI vol is allocated and returned. tail is either
  positive, negative, or absolute (NULL=positive). nover is the
  number of voxels in vol that meet the threshold criteria
  (ie, the number of 1's in volbin).
  -----------------------------------------------------------------*/
MRI *mri_binarize(MRI *vol, float thresh, const char *tail, int invert, MRI *volbin, int *nover)
{
  int r, c, s, f, tailcode;
  float val;
  int b, onval, offval;
  MRI *voltmp;

  if (tail == NULL) tail = "positive";

  /* check the first 3 letters of tail */
  if (!strncasecmp(tail, "positive", 3))
    tailcode = 1;
  else if (!strncasecmp(tail, "negative", 3))
    tailcode = 2;
  else if (!strncasecmp(tail, "absolute", 3))
    tailcode = 3;
  else {
    fprintf(stderr, "mri_binarize: tail = %s unrecoginzed\n", tail);
    return (NULL);
  }

  if (volbin == NULL) {
    voltmp = MRIallocSequence(vol->width, vol->height, vol->depth, MRI_FLOAT, vol->nframes);
    if (voltmp == NULL) return (NULL);
    MRIcopyHeader(vol, voltmp);
    MRIcopyPulseParameters(vol, voltmp);
  }
  else
    voltmp = volbin;

  if (!invert) {
    onval = 1;
    offval = 0;
    printf("NOT INVERTING\n");
  }
  else {
    onval = 0;
    offval = 1;
    printf("INVERTING\n");
  }

  *nover = 0;
  for (r = 0; r < vol->height; r++) {
    for (c = 0; c < vol->width; c++) {
      for (s = 0; s < vol->depth; s++) {
        for (f = 0; f < vol->nframes; f++) {
          // val = MRIFseq_vox(vol,c,r,s,f);
          val = MRIgetVoxVal(vol, c, r, s, f);
          switch (tailcode) {
            case 2:
              val = -val;
              break;
            case 3:
              val = fabs(val);
              break;
          }
          if (val > thresh)
            b = onval;
          else
            b = offval;
          if (b) (*nover)++;
          // MRIFseq_vox(voltmp,c,r,s,f) = b;
          MRIsetVoxVal(voltmp, c, r, s, f, b);
        }
      }
    }
  }

  return (voltmp);
}
/*--------------------------------------------------------
  mri_load_cor_as_float()
  --------------------------------------------------------*/
MRI *mri_load_cor_as_float(char *cordir)
{
  MRI *ucvol;
  MRI *vol;
  int r, c, s;
  float val;

  /* read in the cor as unsigned char */
  ucvol = MRIread(cordir);
  if (ucvol == NULL) return (NULL);

  /* allocate a float volume */
  vol = MRIallocSequence(256, 256, 256, MRI_FLOAT, 1);
  if (vol == NULL) return (NULL);

  for (r = 0; r < vol->height; r++) {
    for (c = 0; c < vol->width; c++) {
      for (s = 0; s < vol->depth; s++) {
        val = (float)(MRIseq_vox(ucvol, c, r, s, 0));
        MRIFseq_vox(vol, c, r, s, 0) = val;
      }
    }
  }

  MRIfree(&ucvol);

  return (vol);
}
/* ---------------------------------------- */
/* not tested */
MRI *mri_load_wfile(char *wfile)
{
  FILE *fp;
  int i, ilat, num, vtx, nvertices;
  int *vtxnum;
  float *wval;
  MRI *w;

  fp = fopen(wfile, "r");
  if (fp == NULL) {
    fprintf(stderr, "ERROR: Progname: mri_load_wfile():\n");
    fprintf(stderr, "Could not open %s\n", wfile);
    fprintf(stderr, "(%s,%d,%s)\n", __FILE__, __LINE__, __DATE__);
    return (NULL);
  }

  fread2(&ilat, fp);
  fread3(&num, fp);

  vtxnum = (int *)calloc(sizeof(int), num);
  wval = (float *)calloc(sizeof(float), num);

  for (i = 0; i < num; i++) {
    fread3(&vtxnum[i], fp);
    wval[i] = freadFloat(fp);
  }
  fclose(fp);

  nvertices = vtxnum[num - 1] + 1;

  w = MRIallocSequence(nvertices, 1, 1, MRI_FLOAT, 1);
  for (i = 0; i < num; i++) {
    vtx = vtxnum[i];
    MRIFseq_vox(w, vtx, 0, 0, 0) = wval[i];
  }

  free(vtxnum);
  free(wval);
  return (w);
}
/*------------------------------------------------------------
  mri_sizeof() - returns the size of the data type of the MRI
  volume (in number of bytes).
  ------------------------------------------------------------*/
size_t mri_sizeof(MRI *vol)
{
  size_t bytes = 0;

  switch (vol->type) {
    case MRI_UCHAR:
      bytes = sizeof(BUFTYPE);
      break;
    case MRI_SHORT:
      bytes = sizeof(short);
      break;
    case MRI_USHRT:
      bytes = sizeof(unsigned short);
      break;
    case MRI_FLOAT:
      bytes = sizeof(float);
      break;
    case MRI_INT:
      bytes = sizeof(int);
      break;
    case MRI_LONG:
      bytes = sizeof(long);
      break;
  }

  return (bytes);
}
/*------------------------------------------------------------
  mri_reshape() --
  ------------------------------------------------------------*/
MRI *mri_reshape(MRI *vol, int ncols, int nrows, int nslices, int nframes)
{
  MRI *outvol;
  int r, c, s, f, nv1, nv2;
  int r2, c2, s2, f2;

  if (vol->nframes == 0) vol->nframes = 1;

  nv1 = vol->width * vol->height * vol->depth * vol->nframes;
  nv2 = ncols * nrows * nslices * nframes;

  if (nv1 != nv2) {
    printf("ERROR: mri_reshape: number of elements cannot change\n");
    printf("  nv1 = %d, nv1 = %d\n", nv1, nv2);
    return (NULL);
  }

  outvol = MRIallocSequence(ncols, nrows, nslices, vol->type, nframes);
  if (outvol == NULL) return (NULL);

  MRIcopyHeader(vol, outvol); /* does not change dimensions */
  MRIcopyPulseParameters(vol, outvol);

  // printf("vol1: %d %d %d %d %d\n",vol->width,vol->height,
  // vol->depth,vol->nframes,mri_sizeof(vol));
  // printf("vol2: %d %d %d %d %d\n",outvol->width,outvol->height,
  // outvol->depth,outvol->nframes,mri_sizeof(outvol));

  c2 = 0;
  r2 = 0;
  s2 = 0;
  f2 = 0;
  for (f = 0; f < vol->nframes; f++) {
    for (s = 0; s < vol->depth; s++) {
      for (r = 0; r < vol->height; r++) {
        for (c = 0; c < vol->width; c++) {
          switch (vol->type) {
            case MRI_UCHAR:
              MRIseq_vox(outvol, c2, r2, s2, f2) = MRIseq_vox(vol, c, r, s, f);
              break;
            case MRI_SHORT:
              MRISseq_vox(outvol, c2, r2, s2, f2) = MRISseq_vox(vol, c, r, s, f);
              break;
            case MRI_USHRT:
              MRIUSseq_vox(outvol, c2, r2, s2, f2) = MRIUSseq_vox(vol, c, r, s, f);
              break;
            case MRI_FLOAT:
              MRIFseq_vox(outvol, c2, r2, s2, f2) = MRIFseq_vox(vol, c, r, s, f);
              break;
            case MRI_INT:
              MRIIseq_vox(outvol, c2, r2, s2, f2) = MRIIseq_vox(vol, c, r, s, f);
              break;
            case MRI_LONG:
              MRILseq_vox(outvol, c2, r2, s2, f2) = MRILseq_vox(vol, c, r, s, f);
              break;
          }

          c2++;
          if (c2 == ncols) {
            c2 = 0;
            r2++;
            if (r2 == nrows) {
              r2 = 0;
              s2++;
              if (s2 == nslices) {
                s2 = 0;
                f2++;
              }
            }
          }
        }
      }
    }
  }

  return (outvol);
}

/*!
  \fn MRI *MRIreshape1d(MRI *src, MRI *trg)
  \brief Reshapes the MRI structure to have nvox columns,
  1 row, and 1 slice; does not change the number of frames.
 */
MRI *MRIreshape1d(MRI *src, MRI *trg)
{
  int ncols, nrows, nslices, nvox, nframes;
  ncols = src->width;
  nrows = src->height;
  nslices = src->depth;
  nvox = ncols * nrows * nslices;
  nframes = src->nframes;

  if (trg == NULL) {
    trg = MRIallocSequence(nvox, 1, 1, src->type, nframes);
    if (trg == NULL) return (NULL);
    MRIcopyHeader(src, trg);
  }
  trg = mri_reshape(src, nvox, 1, 1, nframes);
  return (trg);
}

/*!
  \fn int MRIvol2VolLTA::ReadLTAorTarg(char *fname)
  \brief Reads in either the target volume or an LTA; it figures it
  out for itself, which can be handy when running this class.
 */
int MRIvol2VolLTA::ReadLTAorTarg(char *fname)
{
  int mritype = mri_identify(fname);
  if(mritype != MRI_VOLUME_TYPE_UNKNOWN){
    printf("Reading %s as MRI\n",fname);
    targ = MRIread(fname);
    if(targ==NULL) return(1);
    return(0);
  }
  printf("Reading %s as LTA\n",fname);
  lta = LTAread(fname);
  if(lta==NULL) return(1);
  return(0);
}

/*!
  \fn MRI *MRIvol2VolLTA::vol2vol(MRI *outvol)
  \brief Converts one volume to another. See the notes for the class.
 */
MRI *MRIvol2VolLTA::vol2vol(MRI *outvol)
{
  extern double vg_isEqual_Threshold;
  vg_isEqual_Threshold = 10e-4;

  if(mov==NULL){
    printf("ERROR: MRIvol2VolLTA(): mov volume is NULL\n");
    return(NULL);
  }
  if(targ == NULL && lta == NULL){
    printf("ERROR: MRIvol2VolLTA(): both lta and targ volume are NULL\n");
    return(NULL);
  }

  LTA *ltacopy=NULL;
  if(lta == NULL) ltacopy = TransformRegDat2LTA(targ, mov, NULL);
  else            ltacopy = LTAcopy(lta,NULL);

  // Check whether the source and dest VGs are the same. If they are,
  // check if the registration is the identity. If not, print a warning
  // that we can't tell which direction to go
  if(vg_isEqual(&ltacopy->xforms[0].src, &ltacopy->xforms[0].dst)){
    // If they are the same, check whether the registration is the identity
    // in which case the direction is not important.
    int c,r,IsIdentity=1;
    double val;
    for(r=1; r<=4; r++){
      for(c=1; c<=4; c++){
	val = ltacopy->xforms[0].m_L->rptr[r][c];
	if(r==c && fabs(val-1.0) > 10e-4) IsIdentity = 0;
	if(r!=c && fabs(val)     > 10e-4) IsIdentity = 0;
      }
    }
    if(! IsIdentity){
      // Only print out a warning if if they are the same and the reg is not identity.
      printf("\nINFO: MRISvol2VolLTA(): LTA src and dst vg's are the same and reg is not identity.\n");
      printf("  Make sure you have the direction correct!\n\n");
    }
  }

  // Could check destination against targ if targ != NULL 

  // LTA needs to go from target to mov, so invert if needed.
  VOL_GEOM movvg;
  getVolGeom(mov, &movvg);
  if(vg_isEqual(&ltacopy->xforms[0].src, &movvg)){
    printf("MRIvol2VolLTA(): inverting LTA\n");
    LTAinvert(ltacopy,ltacopy);
  }

  // LTA must be in vox2vox space to use MRIvol2VolVSM()
  if(ltacopy->type != LINEAR_VOX_TO_VOX){
    printf("MRIvol2VolLTA(): changing type to LINEAR_VOX_TO_VOX\n");
    LTAchangeType(ltacopy, LINEAR_VOX_TO_VOX);
  }

  VOL_GEOM *outvg = &(ltacopy->xforms[0].src);
  if(outvol==NULL){
    outvol = MRIallocFromVolGeom(outvg, mov->type, mov->nframes, 0);
    if(outvol==NULL) return(NULL);
    MRIcopyPulseParameters(mov, outvol);
  }
  if(outvol->width != outvg->width || outvol->height != outvg->height ||
     outvol->depth != outvg->depth || outvol->nframes != mov->nframes){
    printf("ERROR: MRIvol2VolLTA(): dimension mismatch\n");
    return(NULL);
  }

  printf("MRIvol2VolLTA(): applying matrix %d %g\n",InterpCode, sinchw);
  MatrixPrint(stdout,ltacopy->xforms[0].m_L);
  int err = MRIvol2VolVSM(mov, outvol, ltacopy->xforms[0].m_L, InterpCode, sinchw, vsm);
  fflush(stdout);
  if(err){
    MRIfree(&outvol);
    return(NULL);
  }
  printf("MRIvol2VolLTA(): done\n");
  return(outvol);
}



/*---------------------------------------------------------------
  MRIvol2Vol() - samples the values of one volume into that of
  another. Handles multiple frames. Can do nearest-neighbor,
  trilinear, and sinc interpolation (sinc may be flaky).
  Use this function instead of vol2vol_linear() in resample.c.

  Vt2s is the 4x4 matrix which converts CRS in the target to
  CRS in the source (ie, it is a vox2vox). If it is NULL, then
  the vox2vox is computed from the src and targ vox2ras matrices,
  assuming the src and targ share the same RAS.

  InterpCode is either: SAMPLE_NEAREST, SAMPLE_TRILINEAR, or
  SAMPLE_SINC.

  param is a generic parameter. For sinc, param is the hw parameter,
  otherwise, it currently has no meaning.
  ---------------------------------------------------------------*/
int MRIvol2Vol(MRI *src, MRI *targ, MATRIX *Vt2s, int InterpCode, float param)
{
  int ct, show_progress_thread;
  int tid = 0;
  float *valvects[_MAX_FS_THREADS];
  int sinchw;
  MATRIX *V2Rsrc = NULL, *invV2Rsrc = NULL, *V2Rtarg = NULL;
  int FreeMats = 0;
  MRI_BSPLINE *bspline = NULL;
  int (*nintfunc)( double );

  /*
    This is a little bit of a hack for the case where there is only
    one slice. If the source and target are aligned by half a voxel
    off, then nint() will never map a target voxel to a valide index
    in the source, and the output will always be 0. nint2() has very
    slightly different behavior that will allow this case to work
    while only mildly affecting more generic cases.
   */
  nintfunc = &nint;
  if(src->width == 1 || src->height == 1 || src->depth == 1)
    nintfunc = &nint2;

#ifdef VERBOSE_MODE

  printf("%s: Begin\n", __FUNCTION__);

  printf("Sources sizes are w=%i h=%i d=%i f=%i\n", src->width, src->height, src->depth, src->nframes);
  printf("src type is %i\n", src->type);

  printf("Target sizes are w=%i h=%i d=%i f=%i\n", src->width, src->height, src->depth, src->nframes);
  printf("targ type is %i\n", targ->type);

  Timer tTotal;
#endif

  if (src->nframes != targ->nframes) {
    printf(
        "ERROR: MRIvol2vol: source and target have different number "
        "of frames\n");
    return (1);
  }

  // Compute vox2vox matrix based on vox2ras of src and target.
  // Assumes that src and targ have same RAS space.
  if (Vt2s == NULL) {
    V2Rsrc = MRIxfmCRS2XYZ(src, 0);
    invV2Rsrc = MatrixInverse(V2Rsrc, NULL);
    V2Rtarg = MRIxfmCRS2XYZ(targ, 0);
    Vt2s = MatrixMultiply(invV2Rsrc, V2Rtarg, NULL);
    FreeMats = 1;
  }
  if (Gdiag_no > 0) {
    printf("MRIvol2Vol: Vt2s Matrix (%d)\n", FreeMats);
    MatrixPrint(stdout, Vt2s);
  }

  sinchw = nint(param);

#ifdef VERBOSE_MODE
  Timer tSample;
#endif

  if (InterpCode == SAMPLE_CUBIC_BSPLINE) bspline = MRItoBSpline(src, NULL, 3);

#ifdef HAVE_OPENMP
  if (omp_get_max_threads() == 1)
    show_progress_thread = 0;
  else
    show_progress_thread = omp_get_max_threads() - 1;  // avoid master thread

  for (tid = 0; tid < _MAX_FS_THREADS; tid++) {
    valvects[tid] = (float *)calloc(sizeof(float), src->nframes);
  }
#else
  show_progress_thread = 0;
  valvects[0] = (float *)calloc(sizeof(float), src->nframes);
#endif

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) shared(show_progress_thread, targ, bspline, src, Vt2s, InterpCode)
#endif
  for (ct = 0; ct < targ->width; ct++) {
    ROMP_PFLB_begin
    
    int rt, st, f;
    int ics, irs, iss;
    float fcs, frs, fss, *valvect;
    double rval;

#ifdef HAVE_OPENMP
    int tid = omp_get_thread_num();
    valvect = valvects[tid];
#else
    valvect = valvects[0];
#endif

    for (rt = 0; rt < targ->height; rt++) {
      for (st = 0; st < targ->depth; st++) {
        /* Column in source corresponding to CRS in Target */
        fcs = Vt2s->rptr[1][1] * ct + Vt2s->rptr[1][2] * rt + Vt2s->rptr[1][3] * st + Vt2s->rptr[1][4];
        ics = nintfunc(fcs);
        if (ics < 0 || ics >= src->width) continue;

        /* Row in source corresponding to CRS in Target */
        frs = Vt2s->rptr[2][1] * ct + Vt2s->rptr[2][2] * rt + Vt2s->rptr[2][3] * st + Vt2s->rptr[2][4];
        irs = nintfunc(frs);
        if (irs < 0 || irs >= src->height) continue;

        /* Slice in source corresponding to CRS in Target */
        fss = Vt2s->rptr[3][1] * ct + Vt2s->rptr[3][2] * rt + Vt2s->rptr[3][3] * st + Vt2s->rptr[3][4];
        iss = nintfunc(fss);
        if (iss < 0 || iss >= src->depth) continue;

        /* Assign output volume values */
        if (InterpCode == SAMPLE_TRILINEAR)
          MRIsampleSeqVolume(src, fcs, frs, fss, valvect, 0, src->nframes - 1);
        else {
          for (f = 0; f < src->nframes; f++) {
            switch (InterpCode) {
              case SAMPLE_NEAREST:
                valvect[f] = MRIgetVoxVal(src, ics, irs, iss, f);
                break;
              case SAMPLE_CUBIC_BSPLINE:
                MRIsampleBSpline(bspline, fcs, frs, fss, f, &rval);
                valvect[f] = rval;
                break;
              case SAMPLE_SINC: /* no multi-frame */
                MRIsincSampleVolume(src, fcs, frs, fss, sinchw, &rval);
                valvect[f] = rval;
                break;
              default:
                printf("ERROR: MRIvol2vol: interpolation method %i unknown\n", InterpCode);
                exit(1);
            }
          }
        }

        for (f = 0; f < src->nframes; f++) MRIsetVoxVal(targ, ct, rt, st, f, valvect[f]);

      } /* target col */
    }   /* target row */
    if (tid == show_progress_thread) exec_progress_callback(ct, targ->width, 0, 1);
    ROMP_PFLB_end
  } /* target slice */
  ROMP_PF_end
  
#ifdef HAVE_OPENMP
  for (tid = 0; tid < _MAX_FS_THREADS; tid++) free(valvects[tid]);
#else
  free(valvects[0]);
#endif

#ifdef VERBOSE_MODE
  int tSampleTime = tSample.milliseconds();
#endif

  if (FreeMats) {
    MatrixFree(&V2Rsrc);
    MatrixFree(&invV2Rsrc);
    MatrixFree(&V2Rtarg);
    MatrixFree(&Vt2s);
  }

  if (bspline) MRIfreeBSpline(&bspline);

#ifdef VERBOSE_MODE
  printf("Timings ------------\n");
  printf("  tSample : %d ms\n", tSampleTime);
  printf("Total     : %d ms\n", tTotal.milliseconds());
  printf("%s: Done\n", __FUNCTION__);
#endif

  return (0);
}
int MRIvol2VolR(MRI *src, MRI *targ, MATRIX *Vt2s, int InterpCode, float param, MATRIX *RRot)
{
  int ct, rt, st, f;
  int ics, irs, iss;
  float fcs, frs, fss;
  float *valvect;
  int sinchw;
  double rval;
  MATRIX *V2Rsrc = NULL, *invV2Rsrc = NULL, *V2Rtarg = NULL;
  int FreeMats = 0;

  MATRIX *RRotT, *Tensor, *RotTensor;
  VECTOR *Vector, *RotVector;
  RRotT = MatrixIdentity(3, NULL);
  Tensor = MatrixIdentity(3, NULL);
  RotTensor = MatrixIdentity(3, NULL);
  RRotT = MatrixTranspose(RRot, NULL);
  Vector = VectorAlloc(3, MATRIX_REAL);
  RotVector = VectorAlloc(3, MATRIX_REAL);

  int nframes = targ->nframes;
  printf("Number of frames: %d\n", nframes);
  if (src->nframes != nframes) {
    printf(
        "ERROR: MRIvol2volR: source and target have different number "
        "of frames\n");
    return (1);
  }
  if (nframes != 3 && nframes != 9) {
    printf("MRIvol2VolR: Wrong number of frames. Exiting!\n");
    return 1;
  }

  // Compute vox2vox matrix based on vox2ras of src and target.
  // Assumes that src and targ have same RAS space.
  if (Vt2s == NULL) {
    V2Rsrc = MRIxfmCRS2XYZ(src, 0);
    invV2Rsrc = MatrixInverse(V2Rsrc, NULL);
    V2Rtarg = MRIxfmCRS2XYZ(targ, 0);
    Vt2s = MatrixMultiply(invV2Rsrc, V2Rtarg, NULL);
    FreeMats = 1;
  }
  if (Gdiag_no > 0) {
    printf("MRIvol2VolR: Vt2s Matrix (%d)\n", FreeMats);
    MatrixPrint(stdout, Vt2s);
  }

  sinchw = nint(param);
  valvect = (float *)calloc(sizeof(float), src->nframes);

  for (ct = 0; ct < targ->width; ct++) {
    for (rt = 0; rt < targ->height; rt++) {
      for (st = 0; st < targ->depth; st++) {
        /* Column in source corresponding to CRS in Target */
        fcs = Vt2s->rptr[1][1] * ct + Vt2s->rptr[1][2] * rt + Vt2s->rptr[1][3] * st + Vt2s->rptr[1][4];
        ics = nint(fcs);
        if (ics < 0 || ics >= src->width) continue;

        /* Row in source corresponding to CRS in Target */
        frs = Vt2s->rptr[2][1] * ct + Vt2s->rptr[2][2] * rt + Vt2s->rptr[2][3] * st + Vt2s->rptr[2][4];
        irs = nint(frs);
        if (irs < 0 || irs >= src->height) continue;

        /* Slice in source corresponding to CRS in Target */
        fss = Vt2s->rptr[3][1] * ct + Vt2s->rptr[3][2] * rt + Vt2s->rptr[3][3] * st + Vt2s->rptr[3][4];
        iss = nint(fss);
        if (iss < 0 || iss >= src->depth) continue;

        /* Assign output volume values */
        if (InterpCode == SAMPLE_TRILINEAR)
          MRIsampleSeqVolume(src, fcs, frs, fss, valvect, 0, src->nframes - 1);
        else {
          for (f = 0; f < src->nframes; f++) {
            switch (InterpCode) {
              case SAMPLE_NEAREST:
                valvect[f] = MRIgetVoxVal(src, ics, irs, iss, f);
                break;
              case SAMPLE_SINC: /* no multi-frame */
                MRIsincSampleVolume(src, fcs, frs, fss, sinchw, &rval);
                valvect[f] = rval;
                break;
            }
          }
        }

        // Two scenarios: either a 3-frame volume (eg.: eigvecs) or
        // a 9-frame volume (eg.: tensors)

        if (nframes == 3) {
          int col;
          for (col = 1; col <= 3; col++) VECTOR_ELT(Vector, col) = valvect[col - 1];

          RotVector = MatrixMultiply(RRot, Vector, NULL);
          // printf("RotVector:\n");
          // MatrixPrint(stdout,RotVector);
          // printf("Vec Loop\n");
          for (col = 1; col <= 3; col++) {
            MRIsetVoxVal(targ, ct, rt, st, col - 1, VECTOR_ELT(RotVector, col));
          }
        }
        else  // if(nframes == 9)
        {
          int row, col;
          for (row = 1; row <= 3; row++)
            for (col = 1; col <= 3; col++) Tensor->rptr[row][col] = valvect[(row - 1) * 3 + col - 1];
          RotTensor = MatrixMultiply(RRot, MatrixMultiply(Tensor, RRotT, NULL), NULL);
          // printf("RotTensor:\n");
          // MatrixPrint(stdout,RotTensor);
          //	    printf("Tensor Loop\n");
          for (row = 1; row <= 3; row++)
            for (col = 1; col <= 3; col++) {
              MRIsetVoxVal(targ, ct, rt, st, (row - 1) * 3 + col - 1, RotTensor->rptr[row][col]);
            }
        }
      } /* target col */
    }   /* target row */
  }     /* target slice */

  free(valvect);
  if (FreeMats) {
    MatrixFree(&V2Rsrc);
    MatrixFree(&invV2Rsrc);
    MatrixFree(&V2Rtarg);
    MatrixFree(&Vt2s);
  }

  return (0);
}
/*-------------------------------------------------------*/
/*!
  \fn int MRIvol2VolTkReg(MRI *mov, MRI *targ, MATRIX *Rtkreg,
      int InterpCode, float param)
  \brief Applies a tkregister matrix to a volume. Computes
    the vox2vox matrix then calls MRIvol2Vol().
  \param mov - source volume
  \param targ - target volume. Must be fully alloced.
  \param Rtkreg - tkreg ras2ras transform from target to source
    (assumes identity if NULL)
  \param InterCode - SAMPLE_NEAREST, SAMPLE_TRILINEAR, SAMPLE_SINC.
    Sinc probably does not work.
 */
int MRIvol2VolTkReg(MRI *mov, MRI *targ, MATRIX *Rtkreg, int InterpCode, float param)
{
  MATRIX *vox2vox = NULL;
  MATRIX *Tmov, *invTmov, *Ttarg;
  int err;

  if (Rtkreg != NULL) {
    // TkReg Vox2RAS matrices
    Tmov = MRIxfmCRS2XYZtkreg(mov);
    invTmov = MatrixInverse(Tmov, NULL);
    Ttarg = MRIxfmCRS2XYZtkreg(targ);
    // vox2vox = invTmov*R*Ttarg
    vox2vox = MatrixMultiply(invTmov, Rtkreg, vox2vox);
    MatrixMultiply(vox2vox, Ttarg, vox2vox);
  }
  else
    vox2vox = NULL;

  // resample
  err = MRIvol2Vol(mov, targ, vox2vox, InterpCode, param);

  if (vox2vox) {
    MatrixFree(&vox2vox);
    MatrixFree(&Tmov);
    MatrixFree(&invTmov);
    MatrixFree(&Ttarg);
  }

  return (err);
}

/*-----------------------------------------------------------*/
/*!
  \fn MRI *MRIvol2VolTLKernel(MRI *src, MRI *targ, MATRIX *Vt2s)
  \brief Computes the trilinear interpolation kernel at each voxel.
  \param src - source volume
  \param targ - target volume
  \param Vt2s - vox2vox transform from target to source (can be NULL)
 */
MRI *MRIvol2VolTLKernel(MRI *src, MRI *targ, MATRIX *Vt2s)
{
  int ct, rt, st, f;
  int ics, irs, iss;
  float fcs, frs, fss;
  double *kvect = NULL;
  MATRIX *V2Rsrc = NULL, *invV2Rsrc = NULL, *V2Rtarg = NULL;
  int FreeMats = 0;
  MRI *kernel;

  kernel = MRIallocSequence(targ->width, targ->height, targ->depth, MRI_FLOAT, 8);
  if (kernel == NULL) return (NULL);
  MRIcopyHeader(targ, kernel);

  // Compute vox2vox matrix based on vox2ras of src and target.
  // Assumes that src and targ have same RAS space.
  if (Vt2s == NULL) {
    V2Rsrc = MRIxfmCRS2XYZ(src, 0);
    invV2Rsrc = MatrixInverse(V2Rsrc, NULL);
    V2Rtarg = MRIxfmCRS2XYZ(targ, 0);
    Vt2s = MatrixMultiply(invV2Rsrc, V2Rtarg, NULL);
    FreeMats = 1;
  }

  for (ct = 0; ct < targ->width; ct++) {
    for (rt = 0; rt < targ->height; rt++) {
      for (st = 0; st < targ->depth; st++) {
        /* Column in source corresponding to CRS in Target */
        fcs = Vt2s->rptr[1][1] * ct + Vt2s->rptr[1][2] * rt + Vt2s->rptr[1][3] * st + Vt2s->rptr[1][4];
        ics = nint(fcs);
        if (ics < 0 || ics >= src->width) continue;

        /* Row in source corresponding to CRS in Target */
        frs = Vt2s->rptr[2][1] * ct + Vt2s->rptr[2][2] * rt + Vt2s->rptr[2][3] * st + Vt2s->rptr[2][4];
        irs = nint(frs);
        if (irs < 0 || irs >= src->height) continue;

        /* Slice in source corresponding to CRS in Target */
        fss = Vt2s->rptr[3][1] * ct + Vt2s->rptr[3][2] * rt + Vt2s->rptr[3][3] * st + Vt2s->rptr[3][4];
        iss = nint(fss);
        if (iss < 0 || iss >= src->depth) continue;

        kvect = MRItrilinKernel(src, fcs, frs, fss, kvect);
        for (f = 0; f < 8; f++) MRIsetVoxVal(targ, ct, rt, st, f, kvect[f]);

      } /* target col */
    }   /* target row */
  }     /* target slice */

  free(kvect);
  if (FreeMats) {
    MatrixFree(&V2Rsrc);
    MatrixFree(&invV2Rsrc);
    MatrixFree(&V2Rtarg);
    MatrixFree(&Vt2s);
  }

  return (0);
}
/*
  \fn MRI *MRImaskAndUpsample(MRI *src, MRI *mask, int UpsampleFactor, int DoConserve, LTA **src2out)

  \brief Masks and upsamples source volume and creates an LTA that
  maps from source voxel to output voxel. mask=NULL, the mask is
  generated from the source. If UpsampleFactor <= 1, upsampling is not
  done. If DoConserve=1, then the upsampled volume is divided by
  UpsampleFactor^3 and so conserves the sum of all voxel
  intensities. When masking is done, a bounding box around the
  non-zero voxels in the mask is created.

 */
MRI *MRImaskAndUpsample(MRI *src, MRI *mask, int UpsampleFactor, int nPad, int DoConserve, LTA **src2out)
{
  MRI *srcmask, *srcus;
  MRI_REGION *region;

  if (mask)
    region = REGIONgetBoundingBox(mask, nPad);
  else
    region = REGIONgetBoundingBox(src, nPad);
  if (Gdiag_no > 0) {
    printf(" mask-and-upsample bounding box %g ",
           ((float)src->width * src->height * src->depth) / (region->dx * region->dy * region->dz));
    REGIONprint(stdout, region);
  }

  srcmask = MRIextractRegion(src, NULL, region);
  if (srcmask == NULL) return (NULL);
  free(region);

  if (UpsampleFactor > 1) {
    if (DoConserve)
      srcus = MRIupsampleNConserve(srcmask, NULL, UpsampleFactor);
    else
      srcus = MRIupsampleN(srcmask, NULL, UpsampleFactor);
  }
  else
    srcus = srcmask;

  *src2out = TransformRegDat2LTA(src, srcus, NULL);  // src2srcus

  if (UpsampleFactor > 1) MRIfree(&srcmask);
  return (srcus);
}
/*---------------------------------------------------------------
  MRIdimMismatch() - checks whether the dimensions on two volumes
  are inconsistent. Returns:
    0 no mismatch
    1 columns mismatch
    2 rows mismatch
    3 slices mismatch
    4 frames mismatch - note: frameflag must = 1 to check frames
  ---------------------------------------------------------------*/
int MRIdimMismatch(const MRI *v1, const MRI *v2, int frameflag)
{
  if (v1->width != v2->width) return (1);
  if (v1->height != v2->height) return (2);
  if (v1->depth != v2->depth) return (3);
  if (frameflag && v1->nframes != v2->nframes) return (4);
  return (0);
}

/*---------------------------------------------------------------
 \fn int MRIfdr2vwth(MRI **vollist, int nvols, int *framelist, double fdr, int signid,
                int log10flag, MRI **masklist, double *vwth, MRI **ovollist)

  MRIfdr2vwth() - computes the voxel-wise threshold needed to realize
  the given False Discovery Rate (FDR) based on the values in the
  given frame. Optionally thresholds the the MRI volume. The input,
  mask, frame, and output are provided as a list. When multiple inputs
  are provided, then a single FDR threshold is computed over all
  values from all inputs. Some masks in the masklist can be null.
  Some outputs inthe outputlist can be null.

  frame - 0-based frame number of input volume to use
  fdr - false dicovery rate, between 0 and 1, eg: .05
  signid -
      0 = use all values regardless of sign
     +1 = use only positive values
     -1 = use only negative values
     If a vertex does not meet the sign criteria, its val2 is 0
  log10flag - interpret vol values as -log10(p)
  mask - binary mask volume. If the mask value at a voxel is 1,
     then its the input voxel value will be used to compute the
     threshold (if it also meets the sign criteria). If the mask is
     0, then ovol will be set to 0 at that location. Pass NULL if
     not using a mask.
  vwth - voxel-wise threshold between 0 and 1. If log10flag is set,
     then vwth = -log10(vwth). ovol voxels with p values
     GREATER than vwth will be 0. Note that this is the same
     as requiring -log10(p) > -log10(vwth).
  ovol - output volume. Pass NULL if no output volume is needed.

  So, for the ovol to be set to something non-zero, the vol must
  meet the sign, mask, and threshold criteria. If vol meets all
  the criteria, then ovol=vol (ie, no log10 transforms). vol
  itself is not altered, unless it is also passed as the ovol.
  Note: only frame=0 in the ovol is used.

  Return Values:
    0 - everything is OK
    1 - something went wrong

  Ref: http://www.sph.umich.edu/~nichols/FDR/FDR.m
  Thresholding of Statistical Maps in Functional Neuroimaging Using
  the False Discovery Rate.  Christopher R. Genovese, Nicole A. Lazar,
  Thomas E. Nichols (2002).  NeuroImage 15:870-878.

  *----------------------------------------------------*/
int MRIfdr2vwth(MRI **vollist,
                int nvols,
                int *framelist,
                double fdr,
                int signid,
                int log10flag,
                MRI **masklist,
                double *vwth,
                MRI **ovollist)
{
  MRI *vol, *mask, *ovol;
  double *p = NULL, val = 0.0, valnull = 0.0, maskval;
  int Nv, np, c, r, s, frame, nthvol;

  Nv = 0;
  for (nthvol = 0; nthvol < nvols; nthvol++) {
    vol = vollist[nthvol];
    frame = framelist[nthvol];
    ovol = NULL;
    mask = NULL;
    if (masklist) mask = masklist[nthvol];
    if (ovollist) ovol = ovollist[nthvol];

    if (vol->nframes <= frame) {
      printf("ERROR: MRIfdr2vwth: frame = %d, must be < nframes = %d\n", frame, vol->nframes);
      return (1);
    }
    if (vol->type != MRI_FLOAT) {
      printf("ERROR: MRIfdr2vwth: input volume is not of type MRI_FLOAT\n");
      return (1);
    }
    if (ovollist != NULL && ovol != NULL) {
      if (ovol->type != MRI_FLOAT) {
        printf("ERROR: MRIfdr2vwth: output volume is not of type MRI_FLOAT\n");
        return (1);
      }
      if (MRIdimMismatch(vol, ovol, 0)) {
        printf("ERROR: MRIfdr2vwth: output/input dimension mismatch\n");
        return (1);
      }
    }
    if (mask != NULL) {
      if (MRIdimMismatch(vol, mask, 0)) {
        printf("ERROR: MRIfdr2vwth: mask/input dimension mismatch\n");
        return (1);
      }
    }
    Nv += (vol->width * vol->height * vol->depth);
  }

  // Package all the p-values into a vector
  p = (double *)calloc(Nv, sizeof(double));
  np = 0;
  for (nthvol = 0; nthvol < nvols; nthvol++) {
    vol = vollist[nthvol];
    frame = framelist[nthvol];
    mask = NULL;
    if (masklist) mask = masklist[nthvol];
    for (c = 0; c < vol->width; c++) {
      for (r = 0; r < vol->height; r++) {
        for (s = 0; s < vol->depth; s++) {
          if (mask) {
            // Must be in the mask if using a mask
            maskval = MRIgetVoxVal(mask, c, r, s, 0);
            if (maskval < 0.5) continue;
          }
          val = MRIFseq_vox(vol, c, r, s, frame);

          // Check the sign
          if (signid == -1 && val > 0) continue;
          if (signid == +1 && val < 0) continue;

          // Get value, convert from log10 if needed
          val = fabs(val);
          if (log10flag) val = pow(10, -val);

          p[np] = val;
          np++;
        }
      }
    }
  }
  printf("MRIfdr2vwth: np = %d, nv = %d\n", np, Nv);

  // Check that something met the match criteria,
  // otherwise return an error
  if (np == 0) {
    printf("ERROR: MRIfdr2vwth: no matching voxels found\n");
    free(p);
    return (1);
  }

  // Compute the voxel-wise threshold using FDR
  *vwth = fdr2vwth(p, np, fdr);
  printf("MRIfdr2vwth: vwth = %lf, log10(vwhth) = %lf\n", *vwth, -log10(*vwth));
  free(p);

  if (ovollist == NULL) {
    // return here if no output
    if (log10flag) *vwth = -log10(*vwth);
    return (0);
  }

  if (log10flag)
    valnull = 0;
  else
    valnull = 1;

  // Perform the thresholding
  for (nthvol = 0; nthvol < nvols; nthvol++) {
    vol = vollist[nthvol];
    frame = framelist[nthvol];
    ovol = NULL;
    mask = NULL;
    if (masklist) mask = masklist[nthvol];
    if (ovollist) ovol = ovollist[nthvol];
    if (ovol == NULL) continue;

    for (c = 0; c < vol->width; c++) {
      for (r = 0; r < vol->height; r++) {
        for (s = 0; s < vol->depth; s++) {
          if (mask) {
            maskval = MRIgetVoxVal(mask, c, r, s, 0);
            if (maskval < 0.5) {
              // Set to null if out of mask
              MRIFseq_vox(ovol, c, r, s, 0) = valnull;
              continue;
            }
          }

          val = MRIFseq_vox(vol, c, r, s, frame);

          if (signid == -1 && val > 0) {
            // Set to null if wrong sign
            MRIFseq_vox(ovol, c, r, s, 0) = valnull;
            continue;
          }
          if (signid == +1 && val < 0) {
            // Set to null if wrong sign
            MRIFseq_vox(ovol, c, r, s, 0) = valnull;
            continue;
          }

          val = fabs(val);
          if (log10flag) val = pow(10, -val);

          if (val > *vwth) {
            // Set to null if greather than thresh
            MRIFseq_vox(ovol, c, r, s, 0) = valnull;
            continue;
          }

          // Otherwise, it meets all criteria, so
          // pass the original value through
          MRIFseq_vox(ovol, c, r, s, 0) = MRIFseq_vox(vol, c, r, s, frame);
        }
      }
    }
  }

  if (log10flag) *vwth = -log10(*vwth);

  return (0);
}

/*-------------------------------------------------------------------
  MRIcovarianceMatrix() - computes the cross-frame (temporal)
  covariance matrix. Returns M=D*D/Nv' where D is an Nt-by-Nv data
  set, where Nt is the number of frames/timepoints and Nv is the
  number of voxels/vertices which equals ncols*nrows*nslices in the
  mri strucutre. If mask is non-null, then voxels/vertices whose value
  in the mask are less than 0.5 are excluded from the computation of
  the covariance matrix, and Nv becomes the number points in the mask.
  ------------------------------------------------------------------*/
MATRIX *MRIcovarianceMatrix(MRI *mri, MRI *mask)
{
  int UseMask = 0, nmask, f1, f2;
  int r, c, s;
  double sum, v1, v2;
  MATRIX *M;

  // Handle masking
  if (mask != NULL) {
    // count number of points in the mask
    nmask = 0;
    for (c = 0; c < mri->width; c++) {
      for (r = 0; r < mri->height; r++) {
        for (s = 0; s < mri->depth; s++) {
          if (MRIgetVoxVal(mask, c, r, s, 0) > 0.5) nmask++;
        }
      }
    }
    // printf("Number of voxels in the mask %d\n",nmask);
    if (nmask == 0) {
      printf("ERROR: no voxels in mask\n");
      return (NULL);
    }
    UseMask = 1;
  }
  else {
    // Otherwise use all voxels/vertices
    nmask = mri->width * mri->height * mri->depth;
    UseMask = 0;
  }

  // Allocate the covariance matrix
  M = MatrixAlloc(mri->nframes, mri->nframes, MATRIX_REAL);

  // Compute the covariance matrix
  for (f1 = 0; f1 < mri->nframes; f1++) {
    // printf("f1 = %d\n",f1);
    for (f2 = f1; f2 < mri->nframes; f2++) {
      sum = 0;
      for (c = 0; c < mri->width; c++) {
        for (r = 0; r < mri->height; r++) {
          for (s = 0; s < mri->depth; s++) {
            if (UseMask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
            v1 = MRIgetVoxVal(mri, c, r, s, f1);
            v2 = MRIgetVoxVal(mri, c, r, s, f2);
            sum += (v1 * v2);
            // printf("%g %g %g\n",v1,v2,sum);
          }  // s
        }    // r
      }      // s
      M->rptr[f1 + 1][f2 + 1] = sum / nmask;
      M->rptr[f2 + 1][f1 + 1] = sum / nmask;
    }  // f1
  }    // f2

  // MatrixPrint(stdout,M);
  return (M);
}

/*--------------------------------------------------------------------
  MRIpca() - computes the SVD/PCA of the input volume D.  pU is the
  matrix of temporal EVs, pS is the vector singular values, and pV arg
  the spatial EVs (reshaped back into a volume).  Note: this uses the
  SVD code from matrix.c, which is not very accurate. If mask is
  non-null, uses only the voxels for which mask>0.5. When complete,
  D = U*S*V';
  -------------------------------------------------------------------*/
int MRIpca(MRI *D, MATRIX **pU, VECTOR **pS, MRI **pV, MRI *mask)
{
  int dim, dim_real, nvoxels, c, r, s, f, UseMask, nmask = 0;
  MATRIX *M, *VV, *UinvS, *Fd, *Fv;
  VECTOR *S2;
  double *sum2, v;

  nvoxels = D->width * D->height * D->depth;
  dim = MIN(nvoxels, D->nframes);

  // Count the number of voxels in the mask
  if (mask != NULL) {
    UseMask = 1;
    for (c = 0; c < D->width; c++)
      for (r = 0; r < D->height; r++)
        for (s = 0; s < D->depth; s++)
          if (UseMask && MRIgetVoxVal(mask, c, r, s, 0) > 0.5) nmask++;
  }
  else {
    UseMask = 0;
    nmask = nvoxels;
  }

  M = MRIcovarianceMatrix(D, mask);
  if (M == NULL) return (1);
  // MatrixWriteTxt("cvm.dat",M);

  // Compute the SVD of the Temporal Cov Matrix
  S2 = RVectorAlloc(D->nframes, MATRIX_REAL);
  *pU = MatrixCopy(M, NULL);      // It's done in-place so make a copy
  VV = MatrixSVD(*pU, S2, NULL);  // M = U*S2*VV, VV = U';
  // MatrixWriteTxt("s2.dat",S2);

  *pS = RVectorAlloc(D->nframes, MATRIX_REAL);
  for (c = 1; c <= dim; c++) {
    S2->rptr[1][c] *= nmask;  // Needs this
    (*pS)->rptr[1][c] = sqrt(S2->rptr[1][c]);
    // Exclude dims for which the singular value is very small
    if (S2->rptr[1][c] < FLT_EPSILON) break;
  }
  dim_real = c - 1;
  if(dim_real < 1) dim_real = 1; // keep at least one

  // Compute U*inv(S) (note square root to go from S2 to S)
  UinvS = MatrixAlloc(D->nframes, dim_real, MATRIX_REAL);
  for (c = 1; c <= dim_real; c++) {
    for (f = 1; f <= D->nframes; f++) UinvS->rptr[f][c] = (double)(*pU)->rptr[f][c] / sqrt(S2->rptr[1][c]);
  }
  // printf("UinvS %d -------------------------\n",dim_real);
  // MatrixPrint(stdout,UinvS);

  // Allocate the spatial EVs
  *pV = MRIallocSequence(D->width, D->height, D->depth, MRI_FLOAT, dim_real);
  MRIcopyHeader(D, *pV);

  // Compute V = D'*U*inv(S)
  sum2 = (double *)calloc(dim_real, sizeof(double));
  Fd = MatrixAlloc(1, D->nframes, MATRIX_REAL);
  Fv = MatrixAlloc(1, dim_real, MATRIX_REAL);
  for (c = 0; c < D->width; c++) {
    for (r = 0; r < D->height; r++) {
      for (s = 0; s < D->depth; s++) {
        if (UseMask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
        for (f = 0; f < D->nframes; f++) Fd->rptr[1][f + 1] = MRIgetVoxVal(D, c, r, s, f);
        MatrixMultiply(Fd, UinvS, Fv);
        for (f = 0; f < dim_real; f++) {
          MRIsetVoxVal(*pV, c, r, s, f, Fv->rptr[1][f + 1]);
          sum2[f] += (Fv->rptr[1][f + 1]) * (Fv->rptr[1][f + 1]);
        }
      }
    }
  }

  /* Normalize to a unit vector */
  for (c = 0; c < D->width; c++) {
    for (r = 0; r < D->height; r++) {
      for (s = 0; s < D->depth; s++) {
        if (UseMask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
        for (f = 0; f < dim_real; f++) {
          v = MRIgetVoxVal(*pV, c, r, s, f);
          MRIsetVoxVal(*pV, c, r, s, f, v / sqrt(sum2[f]));
        }
      }
    }
  }

  free(sum2);
  MatrixFree(&M);
  MatrixFree(&VV);
  MatrixFree(&UinvS);
  MatrixFree(&Fd);
  MatrixFree(&Fv);
  MatrixFree(&S2);

  return (0);
}

/*--------------------------------------------------------------
  PrintPCAStats() - prints pca summary statistics to a stream.
  The summary stats are: (1) nth EV (2) var spanned by nth,
  (3) var spanned by 1-nth EVs, (3) percent var spanned by nth,
  (4) perent var spanned by 1-nth EVs,
  --------------------------------------------------------------*/
int PrintPCAStats(FILE *fp, MATRIX *Spca)
{
  int n;
  double totvar, v, vsum;

  totvar = 0.0;
  for (n = 1; n <= Spca->cols; n++) totvar += (Spca->rptr[1][n] * Spca->rptr[1][n]);

  vsum = 0.0;
  for (n = 1; n <= Spca->cols; n++) {
    v = (Spca->rptr[1][n] * Spca->rptr[1][n]);
    vsum += v;
    fprintf(fp, "%3d   %8.2f %9.2f   %6.2f  %6.2f \n", n, v, vsum, 100 * v / totvar, 100 * vsum / totvar);
  }
  return (0);
}

/*--------------------------------------------------------------
  WritePCAStats() - writes pca summary statistics to a file.
  see PrintPCAStats() for more info.
  --------------------------------------------------------------*/
int WritePCAStats(char *fname, MATRIX *Spca)
{
  FILE *fp;
  fp = fopen(fname, "w");
  if (fp == NULL) {
    printf("ERROR: opening %s\n", fname);
    return (1);
  }
  PrintPCAStats(fp, Spca);
  return (0);
}
/*---------------------------------------------------------------
  MRIsqrt() - computes sqrt(fabs(v)). Calls MRIsquraRoot().
  ---------------------------------------------------------------*/
MRI *MRIsqrt(MRI *invol, MRI *outvol)
{
  outvol = MRIsquareRoot(invol, NULL, outvol);
  return (outvol);
}
/*------------------------------------------------------*/
/*!
 \fn MRI *MRImax(MRI *mri1, MRI *mri2, MRI *out)
 \brief Computes the voxel-by-voxel max between the input MRIs.
 \param MRI *mri1 - first input
 \param MRI *mri2 - second input
 \param MRI *out  - output (can be NULL)
 \return MRI * - pointer to output MRI
 Both inputs must be of the same size, but they can be of
 different data types. If the output is not specified, then
 it gets the data type of mri1. If it is specified, it can
 be any data type.
*/
MRI *MRImax(MRI *mri1, MRI *mri2, MRI *out)
{
  int cols = mri1->width;
  int rows = mri1->height;
  int slices = mri1->depth;
  int frames = mri1->nframes;

  if (!out) {
    out = MRIallocSequence(cols, rows, slices, mri1->type, frames);
    if (!out) {
      printf("ERROR: MRImax: could not alloc output\n");
      return nullptr;
    }
    MRIcopyHeader(mri1, out);  // ordinarily would need to change nframes
  }

  // check dimensions
  if (out->width != cols || out->height != rows || out->depth != slices || out->nframes != frames) {
    printf("ERROR: MRImax: dimension mismatch\n");
    return nullptr;
  }

  for (unsigned int f = 0; f < frames; f++) {
    for (unsigned int c = 0; c < cols; c++) {
      for (unsigned int r = 0; r < rows; r++) {
        for (unsigned int s = 0; s < slices; s++) {
          double v1 = MRIgetVoxVal(mri1, c, r, s, f);
          double v2 = MRIgetVoxVal(mri2, c, r, s, f);
          MRIsetVoxVal(out, c, r, s, f, std::max(v1, v2));
        }
      }
    }
  }

  return (out);
}

/*---------------------------------------------------------------
  MRImaxAbsDiff() - finds the voxel where the two volumes differ
  the most.
  ---------------------------------------------------------------*/
double MRImaxAbsDiff(MRI *vol1, MRI *vol2, int *cmax, int *rmax, int *smax, int *fmax)
{
  int c, r, s, f;
  double v1, v2, maxdiff;

  maxdiff = 0.0;
  for (c = 0; c < vol1->width; c++) {
    for (r = 0; r < vol1->height; r++) {
      for (s = 0; s < vol1->depth; s++) {
        for (f = 0; f < vol1->nframes; f++) {
          v1 = MRIgetVoxVal(vol1, c, r, s, f);
          v2 = MRIgetVoxVal(vol2, c, r, s, f);
          if (maxdiff < fabs(v1 - v2)) {
            maxdiff = fabs(v1 - v2);
            *cmax = c;
            *rmax = r;
            *smax = s;
            *fmax = f;
          }
        }
      }
    }
  }
  return (maxdiff);
}
/* --------------------------------------------------------------- */
MRI *MRImultiplyConst(MRI *src, double vconst, MRI *dst)
{
  int r, c, s, f;
  double v;

  if (dst == NULL) {
    dst = MRIallocSequence(src->width, src->height, src->depth, MRI_FLOAT, src->nframes);
    MRIcopyHeader(src, dst);
  }

  for (c = 0; c < src->width; c++) {
    for (r = 0; r < src->height; r++) {
      for (s = 0; s < src->depth; s++) {
        for (f = 0; f < src->nframes; f++) {
          v = MRIgetVoxVal(src, c, r, s, f);
          MRIsetVoxVal(dst, c, r, s, f, v * vconst);
        }
      }
    }
  }

  return (dst);
}
/* --------------------------------------------------------------- */
MRI *MRIaddConst(MRI *src, double vconst, MRI *dst)
{
  int r, c, s, f;
  double v;

  if (dst == NULL) {
    dst = MRIallocSequence(src->width, src->height, src->depth, MRI_FLOAT, src->nframes);
    MRIcopyHeader(src, dst);
  }

  for (c = 0; c < src->width; c++) {
    for (r = 0; r < src->height; r++) {
      for (s = 0; s < src->depth; s++) {
        for (f = 0; f < src->nframes; f++) {
          v = MRIgetVoxVal(src, c, r, s, f);
          MRIsetVoxVal(dst, c, r, s, f, v + vconst);
        }
      }
    }
  }

  return (dst);
}

/*--------------------------------------------------------------------
  MRIframeBinarize() - creates a binary mask of voxels for which the
  abs(mri->val) of all frames are > thresh. If input mask is not null,
  then that mask is pruned, ie, if a voxel was not in the input mask,
  it will not be in the return mask. If it was in the input mask but
  it does not meet the frame criteria in mri, then it will be set to
  0. Note: if mask is not null, it's values will be changed.
  --------------------------------------------------------------------*/
MRI *MRIframeBinarize(MRI *mri, double thresh, MRI *mask)
{
  int c, r, s, f, n, premask;
  double val, m;
  premask = 1;
  if (!mask) {
    mask = MRIcloneBySpace(mri, MRI_FLOAT, 1);
    MRIclear(mask);
    premask = 0;
  }

  for (c = 0; c < mri->width; c++) {
    for (r = 0; r < mri->height; r++) {
      for (s = 0; s < mri->depth; s++) {
        if (premask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        n = 0;
        for (f = 0; f < mri->nframes; f++) {
          val = MRIgetVoxVal(mri, c, r, s, f);
          if (fabs(val) > thresh) n++;
        }
        if (n == mri->nframes)
          MRIsetVoxVal(mask, c, r, s, 0, 1);
        else
          MRIsetVoxVal(mask, c, r, s, 0, 0);
      }
    }
  }
  return (mask);
}
/*!
  \fn MRI *MRIexp(MRI *mri, double a, double b, MRI *mask, MRI *out)
  \brief Computes a*exp(b*mri). If a mask is supplied, then values
         outside the mask (ie, mask < 0.5) are set to 0 (if out=NULL)
         or the previous value of out.
 */
MRI *MRIexp(MRI *mri, double a, double b, MRI *mask, MRI *out)
{
  int c, r, s, f;
  double val, valout, m;
  int err;

  if (out == NULL) {
    out = MRIcloneBySpace(mri, MRI_FLOAT, -1);
    if (out == NULL) {
      printf("ERROR: MRIexp: could not alloc\n");
      return (NULL);
    }
  }
  else {
    err = MRIdimMismatch(mri, out, 1);
    if (err) {
      printf("ERROR: MRIexp(): output dimension mismatch (%d)\n", err);
      return (NULL);
    }
    if (out->type != MRI_FLOAT) {
      printf("ERROR: MRIexp(): structure passed is not MRI_FLOAT\n");
      return (NULL);
    }
  }

  for (c = 0; c < mri->width; c++) {
    for (r = 0; r < mri->height; r++) {
      for (s = 0; s < mri->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        for (f = 0; f < mri->nframes; f++) {
          val = MRIgetVoxVal(mri, c, r, s, f);
          valout = a * exp(b * val);
          MRIsetVoxVal(out, c, r, s, f, valout);
        }
      }
    }
  }
  return (out);
}

/*!
  \fn MRI *MRIsum(MRI *mri1,MRI *mri2, double a,double b, MRI *mask,MRI *out)
  \brief Computes a*mri1 + b*mri2. If a mask is supplied, then values
         outside the mask (ie, mask < 0.5) are set to 0 (if out=NULL)
         or the previous value of out.
 */
MRI *MRIsum(MRI *mri1, MRI *mri2, double a, double b, MRI *mask, MRI *out)
{
  int c, r, s, f;
  double val1, val2, valout, m;
  int err;

  err = MRIdimMismatch(mri1, mri2, 1);
  if (err) {
    printf("ERROR: MRIsum(): input dimension mismatch (%d)\n", err);
    return (NULL);
  }

  if (out == NULL) {
    out = MRIcloneBySpace(mri1, MRI_FLOAT, -1);
    if (out == NULL) {
      printf("ERROR: MRIsum: could not alloc\n");
      return (NULL);
    }
  }
  else {
    err = MRIdimMismatch(mri1, out, 1);
    if (err) {
      printf("ERROR: MRIsum(): output dimension mismatch (%d)\n", err);
      return (NULL);
    }
  }

  for (c = 0; c < mri1->width; c++) {
    for (r = 0; r < mri1->height; r++) {
      for (s = 0; s < mri1->depth; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        for (f = 0; f < mri1->nframes; f++) {
          val1 = MRIgetVoxVal(mri1, c, r, s, f);
          val2 = MRIgetVoxVal(mri2, c, r, s, f);
          valout = a * val1 + b * val2;
          MRIsetVoxVal(out, c, r, s, f, valout);
        }
      }
    }
  }
  return (out);
}
/*!
  \fn MRI *MRIvote(MRI *in, MRI *mask, MRI *vote)
  \brief Select the most frequently occuring value measured
     across frames in each voxel.
  \param vote - has 2 frames:
     (1) Most freqently occuring value
     (2) Fraction of occurances (ie noccurances/nframes)
  Note: input will be sorted in asc order
 */
MRI *MRIvote(MRI *in, MRI *mask, MRI *vote)
{
  int c, r, s, f, f0, ncols, nrows, nslices, nframes;
  float m;
  double vmax, v, v0;
  int runlen, runlenmax;
  MRI *sorted;

  printf("MRIvote: sorting\n");
  sorted = MRIsort(in, mask, in);  // this sorts the input
  if (sorted == NULL) return (NULL);
  printf("MRIvote: done sorting\n");

  ncols = in->width;
  nrows = in->height;
  nslices = in->depth;
  nframes = in->nframes;

  if (vote == NULL) {
    vote = MRIallocSequence(ncols, nrows, nslices, in->type, 2);
    if (vote == NULL) {
      printf("ERROR: MRIvote: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(in, vote);
    vote->nframes = 2;
  }
  if (in->type != vote->type) {
    printf("ERROR: MRIvote: type mismatch\n");
    return (NULL);
  }
  if (vote->width != ncols || vote->height != nrows || vote->depth != nslices || vote->nframes != 2) {
    printf("ERROR: MRIvote: dimension mismatch\n");
    return (NULL);
  }

  for (c = 0; c < ncols; c++) {
    for (r = 0; r < nrows; r++) {
      for (s = 0; s < nslices; s++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        vmax = 0;
        runlenmax = 0;
        v0 = MRIgetVoxVal(sorted, c, r, s, 0);  // value at start of run
        f0 = 0;                                 // frame at start of run
        f = 1;
        while (f < nframes) {
          v = MRIgetVoxVal(sorted, c, r, s, f);
          if (v0 != v) {
            // new value is different than that of run start
            runlen = f - f0;  // runlength for v0
            if (runlenmax < runlen) {
              runlenmax = runlen;
              vmax = v0;
            }
            v0 = v;
            f0 = f;
          }
          f++;
        }
        // Need to do this one more time in case last value
        // has the longest run
        runlen = f - f0;
        if (runlenmax < runlen) {
          runlenmax = runlen;
          vmax = v0;
          v0 = v;
          f0 = f;
        }
        MRIsetVoxVal(vote, c, r, s, 0, vmax);
        MRIsetVoxVal(vote, c, r, s, 1, (double)runlenmax / nframes);

      }  // slices
    }    // rows
  }      // cols

  return (vote);
}

/*
\fn int MRImostFreqNeighbor(MRI *mri, int c, int r, int s, int f, int delta)
\brief Returns the most frequent value from within a neighborhood of the
column, row, and slice at the given frame. The neighborhood size is within
delta of the center voxel and includes the center voxel. Ie, the number
of voxels searched is (2*delta + 1)^3
*/
int MRImostFreqNeighbor(MRI *mri, int c, int r, int s, int f, int delta)
{
  int dc, dr, ds, cn, rn, sn, nbr;
  int nlist, *list, nmax, nside;

  nside = 2 * delta + 1;
  list = (int *)calloc(sizeof(int), nside * nside * nside);

  nlist = 0;
  for (dc = -delta; dc <= delta; dc++) {
    cn = c + dc;
    if (cn < 0 || cn >= mri->width) continue;
    for (dr = -delta; dr <= delta; dr++) {
      rn = r + dr;
      if (rn < 0 || rn >= mri->height) continue;
      for (ds = -delta; ds <= delta; ds++) {
        sn = s + ds;
        if (sn < 0 || sn >= mri->depth) continue;
        list[nlist] = MRIgetVoxVal(mri, cn, rn, sn, f);
        nlist++;
      }
    }
  }
  nbr = most_frequent_int_list(list, nlist, &nmax);
  free(list);
  return (nbr);
}

/* MRImakeVox2VoxReg() - takes a target volume, a movable volume, a
   registration method, and an optional registration filename, and
   will make an mriTransform in which A is the target non-comformed
   index space, ARAS is the target tkRegRAS space, B is the movable
   tkRegRAS space (registered), and B is the movable non-conformed
   index space. If *transform is NULL, this will make a new
   mriTransform and return it, otherwise it will set the values in the
   *transform passed in.

   In the use case of a functional volume "being registered onto" an
   anatomical volume, the anatomical volume should be the targ
   parameter, and the functional should be the mov parameter.
 */
int MRImakeVox2VoxReg(MRI *targ, MRI *mov, int regtype, char *regname, mriTransformRef *transform)
{
  int retcode = 0;
  char *cur_char, *base_end;
  int err;
  MATRIX *targ_idx_to_tkregras = NULL;
  MATRIX *mov_idx_to_tkregras = NULL;
  struct stat file_info;
  fMRI_REG *reg_info = NULL;
  char regpath[1000];
  char fullregname[1000];
  MATRIX *targ_tkregras_to_mov_tkregras = NULL;
  Trns_tErr trnscode;

  if (NULL == targ) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRImakeVox2VoxReg: targ was NULL"));

  if (NULL == mov) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRImakeVox2VoxReg: mov was NULL"));

  if (regtype < VOX2VOXREGTYPE_FILE || regtype > VOX2VOXREGTYPE_IDENTITY)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRImakeVox2VoxReg: invalid reg type"));

  if (VOX2VOXREGTYPE_FILE == regtype && NULL == regname)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "MRImakeVox2VoxReg: reg type was FILE but "
                 "regname was NULL"));

  /* We need to build three matrices to fill in the three transforms
     in the mriTransform object:

     targ (A) ----------> tkRegRAS unregistered (ARAS)
                             |
                             |
                             V
      mov (B) ----------> tkRegRAS registered (BRAS)

     A->ARAS and B->BRAS take each volume from native index space to
     tkRegRAS space, an RAS space in which the center of the world is
     the center of the volume. We use MRIxfmXRS2XYZtkreg to get these.

     For the ARAS->BRAS matrix, we're either going to use a matrix
     from a registration file or one generated by MRItkRegMtx.
  */

  /* Create the targ A->RAS matrix. */
  targ_idx_to_tkregras = MRIxfmCRS2XYZtkreg(targ);
  if (NULL == targ_idx_to_tkregras) {
    printf("ERROR: MRImakeVox2VoxReg: Couldn't create targ_idx_to_tkregras\n");
    goto error;
  }

  /* Create the mov B->RAS matrix. */
  mov_idx_to_tkregras = MRIxfmCRS2XYZtkreg(mov);
  if (NULL == mov_idx_to_tkregras) {
    printf("ERROR: MRImakeVox2VoxReg: Couldn't create mov_idx_to_tkregras\n");
    goto error;
  }

  /* Now we build the ARAS->BRAS matrix, which is the
     registration. Switch on our registration type. */
  switch (regtype) {
    case VOX2VOXREGTYPE_FILE:
    case VOX2VOXREGTYPE_FIND:

      /* If we're reading a file, copy the file from the input or
      generate one from our data file location. */
      if (VOX2VOXREGTYPE_FILE == regtype) {
	int written = snprintf(fullregname, 1000-1, "%s", regname);
	if( written == (1000-1)) {
	  std::cerr << __FUNCTION__ << ": Truncation writing fullregname" << std::endl;
	}
      }
      else if (VOX2VOXREGTYPE_FIND == regtype) {
        /* Copy the movable volume name and find the last / in the
           file name. From there, copy in "register.dat" for our file
           name. */
	int written = snprintf(regpath, 1000-1, "%s", mov->fname);
	if( written == (1000-1)) {
	  std::cerr << __FUNCTION__ << ": Truncation writing regpath" << std::endl;
	}
        cur_char = regpath;
        base_end = regpath;
        while (NULL != cur_char && '\0' != *cur_char) {
          if ('/' == *cur_char) base_end = cur_char;
          cur_char++;
        }
        *base_end = '\0';
        written = snprintf(fullregname, sizeof(fullregname), "%s/%s", regpath, "register.dat");
	if( written == sizeof(fullregname) ) {
	  std::cerr << __FUNCTION__ << ": Truncation writing fullregname (with regpath)" << std::endl;
	}
      }

      /* Check that the file exists. */
      err = stat(fullregname, &file_info);
      if (0 != err) {
        printf("ERROR: MRImakeVox2VoxReg: Couldn't find registration %s\n", fullregname);
        goto error;
      }

      /* Check if it's a regular file. */
      if (!S_ISREG(file_info.st_mode)) {
        printf("ERROR: MRImakeVox2VoxReg: Couldn't open registration %s\n", fullregname);
        goto error;
      }

      /* Read the registration */
      reg_info = StatReadRegistration(fullregname);
      if (NULL == reg_info) {
        printf("ERROR: MRImakeVox2VoxReg: Problem reading registration %s\n", fullregname);
        goto error;
      }

      /* Copy the registration matrix. */
      targ_tkregras_to_mov_tkregras = MatrixCopy(reg_info->mri2fmri, NULL);

      break;

    case VOX2VOXREGTYPE_IDENTITY:

      /* Use MRItkRegMtx to generate an identity registration between
      the two volumes. */
      targ_tkregras_to_mov_tkregras = MRItkRegMtx(targ, mov, NULL);

      break;
  }

  /* Now look at *transform and create a new one if it doesn't
     exist. */
  if (*transform == NULL) {
    trnscode = Trns_New(transform);
    if (Trns_tErr_NoErr != trnscode) {
      printf("ERROR: MRImakeVox2VoxReg: Error creating mriTransform\n");
      goto error;
    }
  }

#if 0
  printf("targ_idx_to_tkregras:\n");
  MatrixPrint(stdout,targ_idx_to_tkregras);
  printf("mov_idx_to_tkregras:\n");
  MatrixPrint(stdout,mov_idx_to_tkregras);
  printf("targ_tkregras_to_mov_tkregras:\n");
  MatrixPrint(stdout,targ_tkregras_to_mov_tkregras);
#endif

  /* Copy the matrices in. */
  trnscode = Trns_CopyAtoRAS(*transform, targ_idx_to_tkregras);
  if (Trns_tErr_NoErr != trnscode) {
    printf("ERROR: MRImakeVox2VoxReg: Error copying A to RAS\n");
    goto error;
  }
  trnscode = Trns_CopyBtoRAS(*transform, mov_idx_to_tkregras);
  if (Trns_tErr_NoErr != trnscode) {
    printf("ERROR: MRImakeVox2VoxReg: Error copying B to RAS\n");
    goto error;
  }
  trnscode = Trns_CopyARAStoBRAS(*transform, targ_tkregras_to_mov_tkregras);
  if (Trns_tErr_NoErr != trnscode) {
    printf("ERROR: MRImakeVox2VoxReg: Error copying ARAS to BRAS\n");
    goto error;
  }

  goto cleanup;
error:

  if (0 == retcode) retcode = -1;

cleanup:

  if (NULL != targ_idx_to_tkregras) MatrixFree(&targ_idx_to_tkregras);

  if (NULL != mov_idx_to_tkregras) MatrixFree(&mov_idx_to_tkregras);

  if (NULL != reg_info) StatFreeRegistration(&reg_info);

  if (NULL != targ_tkregras_to_mov_tkregras) MatrixFree(&targ_tkregras_to_mov_tkregras);

  return (retcode);
}
/*!
  \fn double MRIsum2All(MRI *mri)
  \brief squares and then sums all the voxels.
*/
double MRIsum2All(MRI *mri)
{
  int c, r, s, f;
  double sum2all, val;

  sum2all = 0;
  for (c = 0; c < mri->width; c++) {
    for (r = 0; r < mri->height; r++) {
      for (s = 0; s < mri->depth; s++) {
        for (f = 0; f < mri->nframes; f++) {
          val = MRIgetVoxVal(mri, c, r, s, f);
          sum2all += (val * val);
        }
      }
    }
  }
  return (sum2all);
}
/*!
  \fn MRI *MRIsquare(MRI *in, MRI *mask, MRI *out)
  \brief Squares the value at each voxel. Values outside
  of the mask are set to 0. mask can be NULL.
*/
MRI *MRIsquare(MRI *in, MRI *mask, MRI *out)
{
  int c, r, s, f;
  double val, mval;

  if (out == NULL) out = MRIclone(in, NULL);

  mval = 1;
  for (c = 0; c < in->width; c++) {
    for (r = 0; r < in->height; r++) {
      for (s = 0; s < in->depth; s++) {
        if (mask) mval = MRIgetVoxVal(mask, c, r, s, 0);
        if (mask) {
          val = MRIgetVoxVal(mask, c, r, s, 0);
          if (val < 0.5) continue;
        }
        for (f = 0; f < in->nframes; f++) {
          if (mval > 0.5)
            val = MRIgetVoxVal(in, c, r, s, f);
          else
            val = 0.0;
          MRIsetVoxVal(out, c, r, s, f, val * val);
        }
      }
    }
  }
  return (out);
}
/*!
  \fn MRI *MRIsquareRoot(MRI *in, MRI *mask, MRI *out)
  \brief Square root of the fabs(value) at each voxel.
  Values outside of the mask are set to 0.  mask can be
  NULL.
*/
MRI *MRIsquareRoot(MRI *in, MRI *mask, MRI *out)
{
  int c, r, s, f;
  double val, mval;

  if (out == NULL) {
    out = MRIallocSequence(in->width, in->height, in->depth, MRI_FLOAT, in->nframes);
    MRIcopyHeader(in, out);
    out->type = MRI_FLOAT;
  }

  mval = 1;
  for (c = 0; c < in->width; c++) {
    for (r = 0; r < in->height; r++) {
      for (s = 0; s < in->depth; s++) {
        if (mask) mval = MRIgetVoxVal(mask, c, r, s, 0);
        for (f = 0; f < in->nframes; f++) {
          if (mval > 0.5) {
            val = MRIgetVoxVal(in, c, r, s, f);
          }
          else
            val = 0.0;
          MRIsetVoxVal(out, c, r, s, f, sqrt(fabs(val)));
        }
      }
    }
  }
  return (out);
}
/*!
  \fn MRI *MRIsqr(MRI *in, MRI *out)
  \brief Square of the value at each voxel.
*/
MRI *MRIsqr(MRI *in, MRI *out)
{
  int c, r, s, f;
  double val;

  if (out == NULL) {
    out = MRIallocSequence(in->width, in->height, in->depth, MRI_FLOAT, in->nframes);
    MRIcopyHeader(in, out);
  }

  for (c = 0; c < in->width; c++) {
    for (r = 0; r < in->height; r++) {
      for (s = 0; s < in->depth; s++) {
        for (f = 0; f < in->nframes; f++) {
          val = MRIgetVoxVal(in, c, r, s, f);
          MRIsetVoxVal(out, c, r, s, f, val * val);
        }
      }
    }
  }
  return (out);
}
/*---------------------------------------------------------------*/
/*!
  \fn MRI *MRIchecker(MRI *mri, MRI *checker)
  \brief Creates a checkerboard pattern. Adjacent columns differ
  by 1. Adjacent rows differ by 2. Adjacent slices differ by a
  sign. Adjacent frames differ by a factor of 10. This is for
  visualizing the voxel sampling. Set fthresh = .1, fmax=3.5,
  and use linear blending.
*/
MRI *MRIchecker(MRI *mri, MRI *checker)
{
  int c, r, s, f;
  double cval = 0, rval = 0, sval = 0, fval = 0;

  if (checker == NULL) {
    checker = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, mri->nframes);
    if (checker == NULL) return (NULL);
    MRIcopyHeader(mri, checker);
    checker->type = MRI_FLOAT;
  }

  // Fill even cols with
  for (c = 0; c < mri->width; c++) {
    if (c % 2 == 0)
      cval = 1;
    else
      cval = 2;
    for (r = 0; r < mri->height; r++) {
      if (r % 2 == 0)
        rval = cval;
      else
        rval = cval + 2;
      for (s = 0; s < mri->depth; s++) {
        if (s % 2 == 0)
          sval = +rval;
        else
          sval = -rval;
        for (f = 0; f < mri->nframes; f++) {
          if (f % 2 == 0)
            fval = sval;
          else
            fval = 10 * sval;
          // printf("%d %d %d %d  %g\n",c,r,s,f,fval);
          MRIsetVoxVal(checker, c, r, s, f, fval);
        }
      }
    }
  }
  return (checker);
}

/*---------------------------------------------------------------*/
/*!
  \fn MRI *MRIgrid(MRI *mri, int dc, int dr, int ds, float val, MRI *grid)
  \brief Creates a grid pattern.
*/
MRI *MRIgrid(MRI *mri, int dc, int dr, int ds, float val, MRI *grid)
{
  int c, r, s, f;

  if (grid == NULL) {
    grid = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, mri->nframes);
    if (grid == NULL) return (NULL);
    MRIcopyHeader(mri, grid);
    grid->type = MRI_FLOAT;
  }

  // Make sure everyone is 0
  MRIconst(mri->width, mri->height, mri->depth, mri->nframes, 0, grid);

  for (c = 0; c < mri->width; c += dc) {
    for (r = 0; r < mri->height; r += dr) {
      for (s = 0; s < mri->depth; s += ds) {
        for (f = 0; f < mri->nframes; f++) {
          MRIsetVoxVal(grid, c, r, s, f, val);
        }
      }
    }
  }
  return (grid);
}

/*-----------------------------------------------------------*/
/*!
  \fn MRI *MRIvol2VolDelta(MRI *mov, MRI *targ, MATRIX *Rt2s)
  \brief Computes the amount by which each voxel moves
  \param mov - source volume
  \param targ - target volume
  \param Rt2s - ras2ras transform from target to source (can be NULL, but why?)
 */
MRI *MRIvol2VolDelta(MRI *mov, MRI *targ, MATRIX *Rt2s)
{
  int ct, rt, st;
  double dx, dy, dz;
  MATRIX *targCRS, *targRAS = NULL, *movRAS = NULL, *targVox2RAS, *targVox2movRAS;
  int FreeMats = 0;
  MRI *delta;

  delta = MRIallocSequence(targ->width, targ->height, targ->depth, MRI_FLOAT, 3);
  if (delta == NULL) return (NULL);
  MRIcopyHeader(targ, delta);

  // Compute ras2ras matrix based on vox2ras of mov and target.
  // Assumes that mov and targ have same RAS space.
  if (Rt2s == NULL) {
    Rt2s = MRItkRegMtx(targ, mov, NULL);
    FreeMats = 1;
  }
  targVox2RAS = MRIxfmCRS2XYZtkreg(targ);
  targVox2movRAS = MatrixMultiply(Rt2s, targVox2RAS, NULL);

  targCRS = MatrixAlloc(4, 1, MATRIX_REAL);
  targCRS->rptr[4][1] = 1;

  for (ct = 0; ct < targ->width; ct++) {
    for (rt = 0; rt < targ->height; rt++) {
      for (st = 0; st < targ->depth; st++) {
        targCRS->rptr[1][1] = ct;
        targCRS->rptr[2][1] = rt;
        targCRS->rptr[3][1] = st;
        targRAS = MatrixMultiply(targVox2RAS, targCRS, targRAS);
        movRAS = MatrixMultiply(targVox2movRAS, targCRS, movRAS);

        dx = targRAS->rptr[1][1] - movRAS->rptr[1][1];
        dy = targRAS->rptr[2][1] - movRAS->rptr[2][1];
        dz = targRAS->rptr[3][1] - movRAS->rptr[3][1];

        MRIsetVoxVal(delta, ct, rt, st, 0, dx);
        MRIsetVoxVal(delta, ct, rt, st, 1, dy);
        MRIsetVoxVal(delta, ct, rt, st, 2, dz);

      } /* target col */
    }   /* target row */
  }     /* target slice */

  if (FreeMats) MatrixFree(&Rt2s);
  MatrixFree(&targCRS);
  MatrixFree(&targRAS);
  MatrixFree(&movRAS);
  MatrixFree(&targVox2RAS);
  MatrixFree(&targVox2movRAS);

  return (delta);
}

/*-----------------------------------------------------------*/
/*!
  \fn MRI *MRIcrop(MRI *mri, int c1, int r1, int s1, int c2, int r2, int s2)
  \brief Crop volume, keep correct goemetry.
  \param mri - volume to be cropped
  \param crs1 - col, row, slice to start cropping (inclusive)
  \param crs2 - col, row, slice to end cropping (inclusive)
 */
MRI *MRIcrop(MRI *mri, int c1, int r1, int s1, int c2, int r2, int s2)
{
  int c, r, s, f, Nc, Nr, Ns;
  MRI *crop;
  MATRIX *Vox2RAS, *crs, *P0;
  double v;

  if (c1 < 0 || c1 >= mri->width || r1 < 0 || r1 >= mri->height || s1 < 0 || s1 >= mri->depth) {
    printf("MRIcrop(): start point %d %d %d out of range\n", c1, r1, s1);
    return (NULL);
  }

  if (c2 < 0 || c2 >= mri->width || r2 < 0 || r2 >= mri->height || s2 < 0 || s2 >= mri->depth) {
    printf("MRIcrop(): end point %d %d %d out of range\n", c2, r2, s2);
    return (NULL);
  }

  // Size of cropped volume. +1 to make inclusive.
  Nc = c2 - c1 + 1;
  Nr = r2 - r1 + 1;
  Ns = s2 - s1 + 1;

  crop = MRIallocSequence(Nc, Nr, Ns, mri->type, mri->nframes);
  MRIcopyHeader(mri, crop);

  // Compute location of 1st vox in cropped volume
  Vox2RAS = MRIxfmCRS2XYZ(mri, 0);
  crs = MatrixAlloc(4, 1, MATRIX_REAL);
  crs->rptr[1][1] = c1;
  crs->rptr[2][1] = r1;
  crs->rptr[3][1] = s1;
  crs->rptr[4][1] = 1;
  P0 = MatrixMultiply(Vox2RAS, crs, NULL);

  // Update header geometry for cropped
  MRIp0ToCRAS(crop, P0->rptr[1][1], P0->rptr[2][1], P0->rptr[3][1]);

  // Fill the value
  for (c = 0; c < crop->width; c++) {
    for (r = 0; r < crop->height; r++) {
      for (s = 0; s < crop->depth; s++) {
        for (f = 0; f < crop->nframes; f++) {
          v = MRIgetVoxVal(mri, c + c1, r + r1, s + s1, f);
          MRIsetVoxVal(crop, c, r, s, f, v);
        }
      }
    }
  }

  MatrixFree(&Vox2RAS);
  MatrixFree(&crs);
  MatrixFree(&P0);

  return (crop);
}

/*-----------------------------------------------------------*/
/*!
  \fn MRI *MRIuncrop(MRI *mri,MRI *crop,
  int c1,int r1,int s1, int c2, int r2, int s2)
  \brief Uncrops volume, keeps correct goemetry.
  \param mri - template for full, uncropped volume
  \param crs1 - col, row, slice in template at
  which cropping started (inclusive)
  \param crs2 - col, row, slice in template at
  which cropping ended   (inclusive)
 */
MRI *MRIuncrop(MRI *mri, MRI *crop, int c1, int r1, int s1, int c2, int r2, int s2)
{
  int c, r, s, f;
  MRI *uncrop;
  double v;

  if (c1 < 0 || c1 >= mri->width || r1 < 0 || r1 >= mri->height || s1 < 0 || s1 >= mri->depth) {
    printf("MRIuncrop(): start point %d %d %d out of range\n", c1, r1, s1);
    return (NULL);
  }

  if (c2 < 0 || c2 >= mri->width || r2 < 0 || r2 >= mri->height || s2 < 0 || r2 >= mri->depth) {
    printf("MRIuncrop(): end point %d %d %d out of range\n", c2, r2, s2);
    return (NULL);
  }

  uncrop = MRIcloneBySpace(mri, crop->type, crop->nframes);

  // Fill the values
  for (c = 0; c < crop->width; c++) {
    for (r = 0; r < crop->height; r++) {
      for (s = 0; s < crop->depth; s++) {
        for (f = 0; f < crop->nframes; f++) {
          v = MRIgetVoxVal(crop, c, r, s, f);
          MRIsetVoxVal(uncrop, c + c1, r + r1, s + s1, f, v);
        }
      }
    }
  }

  return (uncrop);
}
/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIreverseSlices(MRI *in, MRI *out)
  \brief Reverses the order of the slices and updates the
    vox2ras matrix.
*/
MRI *MRIreverseSlices(MRI *in, MRI *out)
{
  int c, r, s, f;
  MATRIX *M, *invM, *Sin, *Sout;
  double v;

  if (in == out) {
    printf("ERROR: MRIreverseSlices(): cannot be done in-place\n");
    return (NULL);
  }

  out = MRIcopy(in, out);
  if (out == NULL) return (NULL);

  // vox2ras for the input
  Sin = MRIxfmCRS2XYZ(in, 0);

  // M converts inCRS to outCRS
  M = MatrixAlloc(4, 4, MATRIX_REAL);
  M->rptr[1][1] = 1.0;
  M->rptr[2][2] = 1.0;
  // for reversal: sliceout = (Nslices-1) - slicein
  M->rptr[3][3] = -1.0;
  M->rptr[3][4] = in->depth - 1.0;
  M->rptr[4][4] = 1.0;
  invM = MatrixInverse(M, NULL);
  if (invM == NULL) {
    printf("ERROR: inverting M\n");
    return (NULL);
  }

  // vox2ras for the output
  Sout = MatrixMultiply(Sin, invM, NULL);
  MRIsetVoxelToRasXform(out, Sout);
  MatrixFree(&Sin);
  MatrixFree(&M);
  MatrixFree(&invM);
  MatrixFree(&Sout);

  for (s = 0; s < out->depth; s++) {
    for (c = 0; c < out->width; c++) {
      for (r = 0; r < out->height; r++) {
        for (f = 0; f < out->nframes; f++) {
          v = MRIgetVoxVal(in, c, r, s, f);
          MRIsetVoxVal(out, c, r, (in->depth - 1.0) - s, f, v);
        }
      }
    }
  }

  return (out);
}
/*----------------------------------------------------------------*/
MRI *MRIcutEndSlices(MRI *mri, int ncut)
{
  MRI *out;
  int nslices;
  int c, r, s, f, scut;
  double v;

  nslices = mri->depth - 2 * ncut;
  if (nslices <= 0) {
    printf("ERROR: MRIcutEndSlices(): ncut = %d, input only has %d \n", ncut, mri->depth);
    return (NULL);
  }

  out = MRIallocSequence(mri->width, mri->height, nslices, mri->type, mri->nframes);
  MRIcopyHeader(mri, out);

  for (c = 0; c < mri->width; c++) {
    for (r = 0; r < mri->height; r++) {
      scut = 0;
      for (s = ncut; s < mri->depth - ncut; s++) {
        for (f = 0; f < mri->nframes; f++) {
          v = MRIgetVoxVal(mri, c, r, s, f);
          MRIsetVoxVal(out, c, r, scut, f, v);
        }
        scut++;
      }
    }
  }
  return (out);
}

/* ----------------------------------------------------------*/
/*!
  \fn int *MRIsegIdListExclude0(MRI *seg, int *pnlist, int frame)
  \brief Returns a list of the unique segmentation ids in the volume
   excluding segid=0 (if there). The number in the list is
   *pnlist. The volume need not be an int or char, but it is probably
   what it will be.
*/
int *MRIsegIdListExclude0(MRI *seg, int *pnlist, int frame)
{
  int *segidlist, *segidlist2, n, m, has0;

  segidlist = MRIsegIdList(seg, pnlist, frame);
  if (segidlist == NULL) return (NULL);
  if (*pnlist == 0) return (NULL);  // not sure this can happen

  has0 = 0;
  for (n = 0; n < *pnlist; n++)
    if (segidlist[n] == 0) has0 = 1;

  if (!has0) return (segidlist);

  segidlist2 = (int *)calloc(sizeof(int), *pnlist - 1);
  m = 0;
  for (n = 0; n < *pnlist; n++) {
    if (segidlist[n] != 0) {
      segidlist2[m] = segidlist[n];
      m++;
    }
  }
  *pnlist = *pnlist - 1;
  free(segidlist);
  return (segidlist2);
}
/* ----------------------------------------------------------*/
/*!
  \fn int *MRIsegIdListNot0(MRI *seg, int *nsegs, int frame)
  \brief Returns a list of the unique segmentation ids in the volume
   excluding segid=0. The number in the list is *nsegs. The volume
   need not be an int or char, but it is probably what it will be.
*/
int *MRIsegIdListNot0(MRI *seg, int *nsegs, int frame)
{
  int *segidlist0, *segidlist, msegs, nthseg;
  segidlist0 = MRIsegIdList(seg, nsegs, frame);
  // remove 0 from the list
  segidlist = (int *)calloc(sizeof(int), *nsegs);
  msegs = 0;
  for (nthseg = 0; nthseg < *nsegs; nthseg++) {
    if (segidlist0[nthseg] != 0) {
      segidlist[msegs] = segidlist0[nthseg];
      msegs++;
    }
  }
  *nsegs = msegs;
  return (segidlist);
}
/* ----------------------------------------------------------*/
/*!
  \fn int *MRIsegIdList(MRI *seg, int *nlist, int frame)
  \brief Returns a list of the unique segmentation ids in
   the volume. The number in the list is *nlist. The volume need not
   be an int or char, but it is probably what it will be.
*/
int *MRIsegIdList(MRI *seg, int *nlist, int frame)
{
  int nvoxels, r, c, s, nth;
  int *tmplist = NULL;
  int *segidlist = NULL;

  nvoxels = seg->width * seg->height * seg->depth;
  tmplist = (int *)calloc(sizeof(int), nvoxels);

  // First, load all voxels into a list
  nth = 0;
  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        tmplist[nth] = (int)MRIgetVoxVal(seg, c, r, s, frame);
        nth++;
      }
    }
  }

  segidlist = unqiue_int_list(tmplist, nvoxels, nlist);
  free(tmplist);
  // for(nth=0; nth < *nlist; nth++)
  // printf("%3d %3d\n",nth,segidlist[nth]);
  return (segidlist);
}

/* ----------------------------------------------------------*/
/*!
  \fn double *MRIsegDice(MRI *seg1, MRI *seg2, int *nsegs, int **segidlist)
  \brief Computes dice coefficient for each segmentation. seg1 and seg2
  should have the same number of segs. Note: to get the name of
  the seg, CTABcopyName(ctab,segidlist[k],tmpstr,sizeof(tmpstr));
*/
double *MRIsegDice(MRI *seg1, MRI *seg2, int *nsegs, int **segidlist)
{
  int k, c, r, s, id1, id2, k1 = 0, k2 = 0;
  int nsegid1, *segidlist1;
  int nsegid2, *segidlist2;
  int *n1, *n2, *n12;
  double *dice;
  *nsegs = -1;

  // Extract a unique, sorted list of the ids
  segidlist1 = MRIsegIdList(seg1, &nsegid1, 0);
  segidlist2 = MRIsegIdList(seg1, &nsegid2, 0);

  if (nsegid1 != nsegid2) {
    printf("ERROR: MRIsegDice(): nsegs do not match %d %d\n", nsegid1, nsegid2);
    return (NULL);
  }
  printf("MRIsegDice(): found %d segs\n", nsegid1);
  *nsegs = nsegid1;

  n1 = (int *)calloc(nsegid1, sizeof(int));
  n2 = (int *)calloc(nsegid1, sizeof(int));
  n12 = (int *)calloc(nsegid1, sizeof(int));

  for (c = 0; c < seg1->width; c++) {
    for (r = 0; r < seg1->height; r++) {
      for (s = 0; s < seg1->depth; s++) {
        // segid for 1st seg vol
        id1 = MRIgetVoxVal(seg1, c, r, s, 0);
        for (k = 0; k < nsegid1; k++) {
          if (id1 == segidlist1[k]) {
            k1 = k;
            break;
          }
        }
        // segid for 2nd seg vol
        id2 = MRIgetVoxVal(seg2, c, r, s, 0);
        for (k = 0; k < nsegid2; k++) {
          if (id2 == segidlist2[k]) {
            k2 = k;
            break;
          }
        }
        n1[k1]++;
        n2[k2]++;
        if (id1 == id2) n12[k1]++;
      }
    }
  }

  dice = (double *)calloc(nsegid1, sizeof(double));
  for (k = 0; k < nsegid1; k++) dice[k] = (double)n12[k] / ((n1[k] + n2[k]) / 2.0);

  free(n1);
  free(n2);
  free(n12);
  free(segidlist2);
  *segidlist = segidlist1;

  return (dice);
}

/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIsegDiff(MRI *old, MRI *curr, int *DiffFlag)
  \brief Determines differences between old and curr segmentation.
  Voxels that are different will take the value of the curr
  segmentation. Voxels that are NOT different will take the
  value of VOXEL_UNCHANGED (defined as 256 in cma.h and LUT.txt).
  This allows the diff to be loaded as a segmentation in tkmedit.
  If a difference was detected DiffFlag is set to 1, otherwise 0.
  Note that if there is no difference, all voxels will have the
  value of VOXEL_UNCHANGED. See also MRIsegMergeDiff().
*/
MRI *MRIsegDiff(MRI *old, MRI *curr, int *DiffFlag)
{
  MRI *diff;
  int c, r, s;
  int vold, vnew, vdiff;

  diff = MRIallocSequence(curr->width, curr->height, curr->depth, MRI_INT, 1);
  MRIcopyHeader(curr, diff);

  *DiffFlag = 0;
  for (c = 0; c < curr->width; c++) {
    for (r = 0; r < curr->height; r++) {
      for (s = 0; s < curr->depth; s++) {
        vold = MRIgetVoxVal(old, c, r, s, 0);
        vnew = MRIgetVoxVal(curr, c, r, s, 0);
        if (vold == vnew)
          vdiff = VOXEL_UNCHANGED;
        else {
          vdiff = vnew;
          *DiffFlag = 1;
        }
        MRIsetVoxVal(diff, c, r, s, 0, vdiff);
      }
    }
  }

  return (diff);
}

/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIsegMergeDiff(MRI *old, MRI *diff)
  \brief Merges a segmentation with a "diff" segmentation to create a
  new segmentation. Voxels that have a diff value of VOXEL_UNCHANGED
  (cma.h) take their value from the "old" segmentation. Voxels that
  have a diff value other than VOXEL_UNCHANGED take their value from
  the "diff" segmentation.  See also MRIsegDiff().
*/
MRI *MRIsegMergeDiff(MRI *old, MRI *diff)
{
  MRI *curr;
  int c, r, s;
  int vold, vdiff, vnew;

  curr = MRIallocSequence(old->width, old->height, old->depth, MRI_INT, 1);
  MRIcopyHeader(old, curr);

  for (c = 0; c < curr->width; c++) {
    for (r = 0; r < curr->height; r++) {
      for (s = 0; s < curr->depth; s++) {
        vold = MRIgetVoxVal(old, c, r, s, 0);
        vdiff = MRIgetVoxVal(diff, c, r, s, 0);
        if (vdiff == VOXEL_UNCHANGED)
          vnew = vold;
        else
          vnew = MRIgetVoxVal(diff, c, r, s, 0);
        MRIsetVoxVal(curr, c, r, s, 0, vnew);
      }
    }
  }

  return (curr);
}

/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIhalfCosBias(MRI *in, double alpha, MRI *out)
  \brief Applies intensity bias using a half cosine. Each
    slice is rescaled separately. The scale factor is given
    by f = (alpha*cos(pi*s/(Ns-1))+(2-alpha))/2. Alpha is
    a control parameter. When 0, there is no bias. 1 is full
    bias.
*/
MRI *MRIhalfCosBias(MRI *in, double alpha, MRI *out)
{
  int c, r, s, f;
  double v, w;

  out = MRIcopy(in, out);

  for (s = 0; s < out->depth; s++) {
    w = (alpha * cos(M_PI * s / (out->depth - 1.0)) + (2.0 - alpha)) / 2.0;
    // printf("%2d %g\n",s,w);
    for (c = 0; c < out->width; c++) {
      for (r = 0; r < out->height; r++) {
        for (f = 0; f < out->nframes; f++) {
          v = MRIgetVoxVal(in, c, r, s, f);
          MRIsetVoxVal(out, c, r, s, f, w * v);
        }
      }
    }
  }

  return (out);
}



int MRIvol2VolVSM(MRI *src, MRI *targ, MATRIX *Vt2s, int InterpCode, float param, MRI *vsm, int pedir)
{
  int ct, rt, st, f;
  int ics, irs, iss, cvsm, rvsm;
  float fcs, frs, fss;
  float *valvect, dvsm;
  int sinchw;
  double rval, v;
  MATRIX *V2Rsrc = NULL, *invV2Rsrc = NULL, *V2Rtarg = NULL;
  MATRIX *crsT = NULL, *crsS = NULL;
  int FreeMats = 0;

  if(DIAG_VERBOSE_ON) {
    printf("Using MRIvol2VolVSM\n");
    printf("MRIvol2VolVSM interp=%d, param=%g, pedir=%d\n",InterpCode,param,pedir);
  }

  if(src->nframes != targ->nframes) {
    printf("ERROR: MRIvol2volVSM(): source and target have different number of frames\n"); 
    return (1);
  }
  if(vsm){
    if(abs(pedir) != 1 && abs(pedir) != 2 && abs(pedir) != 3){
      printf("ERROR: MRIvol2volVSM: pedir=%d, must be +/-1, +/-2, +/-3\n",pedir);
      exit(1);      
    }
  }

  // Compute vox2vox matrix based on vox2ras of src and target.
  // Assumes that src and targ have same RAS space.
  if (Vt2s == NULL) {
    V2Rsrc = MRIxfmCRS2XYZ(src, 0);
    invV2Rsrc = MatrixInverse(V2Rsrc, NULL);
    V2Rtarg = MRIxfmCRS2XYZ(targ, 0);
    Vt2s = MatrixMultiply(invV2Rsrc, V2Rtarg, NULL);
    FreeMats = 1;
  }
  if (Gdiag_no > 0) {
    printf("MRIvol2Vol: Vt2s Matrix (%d)\n", FreeMats);
    MatrixPrint(stdout, Vt2s);
  }

  sinchw = nint(param);
  valvect = (float *)calloc(sizeof(float), src->nframes);

  MRI_BSPLINE *bspline = NULL;
  if (InterpCode == SAMPLE_CUBIC_BSPLINE) bspline = MRItoBSpline(src, NULL, 3);

  crsT = MatrixAlloc(4, 1, MATRIX_REAL);
  crsT->rptr[4][1] = 1;
  crsS = MatrixAlloc(4, 1, MATRIX_REAL);
  for (ct = 0; ct < targ->width; ct++) {
    for (rt = 0; rt < targ->height; rt++) {
      for (st = 0; st < targ->depth; st++) {
        // Compute CRS in VSM space
        crsT->rptr[1][1] = ct;
        crsT->rptr[2][1] = rt;
        crsT->rptr[3][1] = st;
        crsS = MatrixMultiply(Vt2s, crsT, crsS);

        fcs = crsS->rptr[1][1];
        frs = crsS->rptr[2][1];
        fss = crsS->rptr[3][1];
        ics = nint(fcs);
        irs = nint(frs);
        iss = nint(fss);

        if (ics < 0 || ics >= src->width) continue;
        if (irs < 0 || irs >= src->height) continue;
        if (iss < 0 || iss >= src->depth) continue;

        if (vsm) {
          /* Compute the voxel shift (converts from vsm
             space to mov space). This does a 3d interp to
             get vsm, not sure if really want a 2d*/
          // Dont sample outside the BO mask
          cvsm = floor(fcs);
          rvsm = floor(frs);
          if (cvsm < 0 || cvsm + 1 >= src->width) continue;
          if (rvsm < 0 || rvsm + 1 >= src->height) continue;
          v = MRIgetVoxVal(vsm, cvsm, rvsm, iss, 0);
          if (fabs(v) < FLT_MIN) continue;
          v = MRIgetVoxVal(vsm, cvsm + 1, rvsm, iss, 0);
          if (fabs(v) < FLT_MIN) continue;
          v = MRIgetVoxVal(vsm, cvsm, rvsm + 1, iss, 0);
          if (fabs(v) < FLT_MIN) continue;
          v = MRIgetVoxVal(vsm, cvsm + 1, rvsm + 1, iss, 0);
          if (fabs(v) < FLT_MIN) continue;
          MRIsampleSeqVolume(vsm, fcs, frs, fss, &dvsm, 0, 0);
          if(dvsm == 0) continue;
	  if(abs(pedir) == 1){
	    fcs += (dvsm*FSIGN(pedir));
	    ics =  nint(fcs);
	    if(ics < 0 || ics >= src->width) continue;
	  }
	  if(abs(pedir) == 2){
	    frs += (dvsm*FSIGN(pedir));
	    irs =  nint(frs);
	    if(irs < 0 || irs >= src->height) continue;
	  }
	  if(abs(pedir) == 3){
	    if(dvsm == 0) continue;
	    fss += (dvsm*FSIGN(pedir));
	    iss = nint(fss);
	    if(iss < 0 || iss >= src->depth) continue;
	  }
        }

        /* Assign output volume values */
        if (InterpCode == SAMPLE_TRILINEAR)
          MRIsampleSeqVolume(src, fcs, frs, fss, valvect, 0, src->nframes - 1);
        else {
          for (f = 0; f < src->nframes; f++) {
            switch (InterpCode) {
              case SAMPLE_NEAREST:
                valvect[f] = MRIgetVoxVal(src, ics, irs, iss, f);
                break;
              case SAMPLE_CUBIC_BSPLINE:
                MRIsampleBSpline(bspline, fcs, frs, fss, f, &rval);
                valvect[f] = rval;
                break;
              case SAMPLE_SINC: /* no multi-frame */
                MRIsincSampleVolume(src, fcs, frs, fss, sinchw, &rval);
                valvect[f] = rval;
                break;
              default:
                printf("ERROR: MRIvol2volVSM: interpolation method %i unknown\n", InterpCode);
                exit(1);
            }
          }
        }

        for (f = 0; f < src->nframes; f++) MRIsetVoxVal(targ, ct, rt, st, f, valvect[f]);

      } /* target col */
    }   /* target row */
  }     /* target slice */

  free(valvect);
  MatrixFree(&crsS);
  MatrixFree(&crsT);
  if (bspline) MRIfreeBSpline(&bspline);
  if (FreeMats) {
    MatrixFree(&V2Rsrc);
    MatrixFree(&invV2Rsrc);
    MatrixFree(&V2Rtarg);
    MatrixFree(&Vt2s);
  }

  return (0);
}


void MRIConvertSurfaceVertexCoordinates(MRIS* mris, MRI* vol)
{
  int const nvertices = mris->nvertices;

  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  int vno;
  for (vno = 0; vno < nvertices; vno++) {

    VERTEX* const v = &mris->vertices[vno];

    double vx, vy, vz;
    MRIsurfaceRASToVoxel( vol,
                          v->x, v->y, v->z,
                          &vx, &vy, &vz);

    MRISsetXYZ(mris, vno, vx,vy,vz);
  }
}



/*---------------------------------------------------------------*/

MRI *MRIvol2surfVSM(const MRI *SrcVol,
                    const MATRIX *Rtk,
                    const MRI_SURFACE *TrgSurf,
                    const MRI *vsm,
                    int InterpMethod,
                    MRI *SrcHitVol,
                    float ProjFrac,
                    int ProjType,
                    int nskip,
                    MRI *TrgVol, int pedir)
{
  MATRIX *ras2vox, *vox2ras;
  AffineVector Scrs, Txyz;
  AffineMatrix ras2voxAffine;
  int irow, icol, islc; /* integer row, col, slc in source */
  int cvsm, rvsm;
  float frow, fcol, fslc; /* float row, col, slc in source */
  float srcval, *valvect, shift;
  int frm, vtx, nhits, err;
  double rval, val;
  float Tx, Ty, Tz;
  const VERTEX *v;

#ifdef MRI2_TIMERS
  Timer tLoop;
#endif

  if (vsm) {
    err = MRIdimMismatch(vsm, SrcVol, 0);
    if (err) {
      printf("ERROR: MRIvol2surfVSM: vsm dimension mismatch %d\n", err);
      exit(1);
    }
    if(abs(pedir) != 1 && abs(pedir) != 2 && abs(pedir) != 3){
      printf("ERROR: MRIvol2surfVSM: pedir=%d, must be +/-1, +/-2, +/-3\n",pedir);
      exit(1);      
    }
  }
  if(DIAG_VERBOSE_ON)  printf("MRIvol2surfVSM interp=%d, nskip=%d, pedir=%d\n",InterpMethod,nskip,pedir);

  vox2ras = MRIxfmCRS2XYZtkreg(SrcVol);
  ras2vox = MatrixInverse(vox2ras, NULL);
  if (Rtk != NULL) MatrixMultiply(ras2vox, Rtk, ras2vox);
  MatrixFree(&vox2ras);
  // ras2vox now converts surfacs RAS to SrcVol vox

  /* allocate a "volume" to hold the output */
  if (TrgVol == NULL) {
    TrgVol = MRIallocSequence(TrgSurf->nvertices, 1, 1, MRI_FLOAT, SrcVol->nframes);
    if (TrgVol == NULL) return (NULL);
    MRIcopyHeader(SrcVol, TrgVol);
  }
  else {
    if (TrgVol->width != TrgSurf->nvertices || TrgVol->nframes != SrcVol->nframes) {
      printf("ERROR: MRIvol2surfVSM: dimension mismatch (%d,%d), or (%d,%d)\n",
             TrgVol->width,
             TrgSurf->nvertices,
             TrgVol->nframes,
             SrcVol->nframes);
      return (NULL);
    }
    // make sure all values are zero
    MRIconst(TrgVol->width, TrgVol->height, TrgVol->depth, 1, 0, TrgVol);
  }
  // Dims here are meaningless, but setting to 1 means "volume" will be
  // number of vertices.
  TrgVol->xsize = 1;
  TrgVol->ysize = 1;
  TrgVol->zsize = 1;

  /* Zero the source hit volume */
  if (SrcHitVol != NULL) MRIconst(SrcHitVol->width, SrcHitVol->height, SrcHitVol->depth, 1, 0, SrcHitVol);

  srcval = 0;
  valvect = (float *)calloc(sizeof(float), SrcVol->nframes);
  nhits = 0;

  SetAffineMatrix(&ras2voxAffine, ras2vox);

  MRI_BSPLINE *bspline = NULL;
  if (InterpMethod == SAMPLE_CUBIC_BSPLINE) bspline = MRItoBSpline(SrcVol, NULL, 3);

/*--- loop through each vertex ---*/
#ifdef MRI2_TIMERS
  tLoop.reset();
  unsigned int skipped = 0;
#endif
  for (vtx = 0; vtx < TrgSurf->nvertices; vtx += nskip) {
    v = &TrgSurf->vertices[vtx];
    if (v->ripflag) {
#ifdef MRI2_TIMERS
      skipped++;
#endif
      continue;
    }

    if (ProjFrac != 0.0) {
      if (ProjType == 0)
        ProjNormDist(&Tx, &Ty, &Tz, TrgSurf, vtx, ProjFrac);
      else
        ProjNormFracThick(&Tx, &Ty, &Tz, TrgSurf, vtx, ProjFrac);
    }
    else {
      Tx = v->x;
      Ty = v->y;
      Tz = v->z;
    }

    /* Load the Target xyz vector */
    SetAffineVector(&Txyz, Tx, Ty, Tz);
    /* Compute the corresponding Source col-row-slc vector */
    AffineMV(&Scrs, &ras2voxAffine, &Txyz);
    GetAffineVector(&Scrs, &fcol, &frow, &fslc);

    icol = nint(fcol);
    irow = nint(frow);
    islc = nint(fslc);

    /* check that the point is in the bounds of the volume */
    if (irow < 0 || irow >= SrcVol->height || icol < 0 || icol >= SrcVol->width || islc < 0 || islc >= SrcVol->depth)
      continue;

    if (vsm) {
      /* Compute the voxel shift (converts from vsm
         space to mov space). This does a 3d interp to
         get vsm, not sure if really want a 2d*/
      // Dont sample outside the BO mask
      cvsm = floor(fcol);
      rvsm = floor(frow);
      if (cvsm < 0 || cvsm + 1 >= vsm->width) continue;
      if (rvsm < 0 || rvsm + 1 >= vsm->height) continue;
      val = MRIgetVoxVal(vsm, cvsm, rvsm, islc, 0);
      if (fabs(val) < FLT_MIN) continue;
      val = MRIgetVoxVal(vsm, cvsm + 1, rvsm, islc, 0);
      if (fabs(val) < FLT_MIN) continue;
      val = MRIgetVoxVal(vsm, cvsm, rvsm + 1, islc, 0);
      if (fabs(val) < FLT_MIN) continue;
      val = MRIgetVoxVal(vsm, cvsm + 1, rvsm + 1, islc, 0);
      if (fabs(val) < FLT_MIN) continue;
      MRIsampleSeqVolume(vsm, fcol, frow, fslc, &shift, 0, 0);
      if(shift == 0) continue;
      if(abs(pedir) == 1){
	fcol += (shift*FSIGN(pedir));
	icol =  nint(fcol);
	if(icol < 0 || icol >= SrcVol->width) continue;
      }
      if(abs(pedir) == 2){
	frow += (shift*FSIGN(pedir));
	irow =  nint(frow);
	if(irow < 0 || irow >= SrcVol->height) continue;
      }
      if(abs(pedir) == 3){
	if(shift == 0) continue;
	fslc += (shift*FSIGN(pedir));
	islc = nint(fslc);
	if(islc < 0 || islc >= SrcVol->depth) continue;
      }
    }

#if 0
    if (Gdiag_no == vtx)
    {
      printf("diag -----------------------------\n");
      printf("vtx = %d  %g %g %g\n",vtx,Txyz->rptr[1][1],
             Txyz->rptr[2][1],Txyz->rptr[3][1]);
      printf("fCRS  %g %g %g\n",Scrs->rptr[1][1],
             Scrs->rptr[2][1],Scrs->rptr[3][1]);
      printf("CRS  %d %d %d\n",icol,irow,islc);
    }
#endif

    /* only gets here if it is in bounds */
    nhits++;

    /* Assign output volume values */
    if (InterpMethod == SAMPLE_TRILINEAR) {
      MRIsampleSeqVolume(SrcVol, fcol, frow, fslc, valvect, 0, SrcVol->nframes - 1);
      if (Gdiag_no == vtx) printf("val = %f\n", valvect[0]);
      for (frm = 0; frm < SrcVol->nframes; frm++) MRIFseq_vox(TrgVol, vtx, 0, 0, frm) = valvect[frm];
    }
    else {
      for (frm = 0; frm < SrcVol->nframes; frm++) {
        switch (InterpMethod) {
          case SAMPLE_NEAREST:
            srcval = MRIgetVoxVal(SrcVol, icol, irow, islc, frm);
            break;
          case SAMPLE_CUBIC_BSPLINE:
            MRIsampleBSpline(bspline, fcol, frow, fslc, frm, &rval);
            srcval = rval;
            break;
          case SAMPLE_SINC: /* no multi-frame */
            MRIsincSampleVolume(SrcVol, fcol, frow, fslc, 5, &rval);
            srcval = rval;
            break;
          default:
            printf("ERROR: MRIvol2surfVSM: interpolation method %i unknown\n", InterpMethod);
            exit(1);
        }  // switch
        MRIFseq_vox(TrgVol, vtx, 0, 0, frm) = srcval;
        if (Gdiag_no == vtx) printf("val[%d] = %f\n", frm, srcval);
      }  // for
    }    // else
    if (SrcHitVol != NULL) MRIFseq_vox(SrcHitVol, icol, irow, islc, 0)++;
  }
#ifdef MRI2_TIMERS
  printf("%s: Main Loop complete in %d ms (%6u %6u)\n", __FUNCTION__, tLoop.milliseconds(), skipped, nhits);
#endif

  MatrixFree(&ras2vox);
  free(valvect);
  if (bspline) MRIfreeBSpline(&bspline);

  // printf("vol2surf_linear: nhits = %d/%d\n",nhits,TrgSurf->nvertices);

  return (TrgVol);
}

int MRIvol2VolTkRegVSM(MRI *mov, MRI *targ, MATRIX *Rtkreg, int InterpCode, float param, MRI *vsm, int pedir)
{
  MATRIX *vox2vox = NULL;
  MATRIX *Tmov, *invTmov, *Ttarg;
  int err;

  if (Rtkreg != NULL) {
    // TkReg Vox2RAS matrices
    Tmov = MRIxfmCRS2XYZtkreg(mov);
    invTmov = MatrixInverse(Tmov, NULL);
    Ttarg = MRIxfmCRS2XYZtkreg(targ);
    // vox2vox = invTmov*R*Ttarg
    vox2vox = MatrixMultiply(invTmov, Rtkreg, vox2vox);
    MatrixMultiply(vox2vox, Ttarg, vox2vox);
  }
  else
    vox2vox = NULL;

  // resample
  err = MRIvol2VolVSM(mov, targ, vox2vox, InterpCode, param, vsm, pedir);

  if (vox2vox) {
    MatrixFree(&vox2vox);
    MatrixFree(&Tmov);
    MatrixFree(&invTmov);
    MatrixFree(&Ttarg);
  }

  return (err);
}
/*
  \fn MRI *MRIvol2VolFill(MRI *src, LTA *lta, int DoConserve, MRI *outfill)
  \brief Maps one volume to another by summing (DoConserve=1) or
  averaving (DoConserve=0) all of the source voxels that land in a given
  target voxel. If DoConserve=1 is used, then the sum of all the voxels
  in the source equals that in the target.  The target geometry is
  taken from the LTA. The LTA can be source-to-target or the reverse.
  The direction is automatically determined and the inverse used if
  necessary. Uses nearest neighbor interpolation. This function should
  generally be used when the target has a much larger voxel size than
  the source. It is possible that no source voxels map to a given
  target voxel (but probably not likely if targ is lowres).
 */
MRI *MRIvol2VolFill(MRI *src, MRI *mask, LTA *lta, int UpsampleFactor, int DoConserve, MRI *outfill)
{
  int c, r, s, ct, rt, st, f, nhits, nPad = 2;
  MATRIX *crssrc, *crstarg, *v2v, *vmusinv;
  double vsrc, vout;
  LTA *ltatmp, *src2srcmus = NULL;
  MRI *hitmap = NULL;
  VOL_GEOM vgtarg;
  MRI *srcmus = NULL;

  if (lta->num_xforms > 1) {
    printf("ERROR: MRIvol2VolFill(): LTA can only have one xform\n");
    return (NULL);
  }
  if (!LTAmriIsSource(lta, src) && !LTAmriIsTarget(lta, src)) {
    printf("ERROR: MRIvol2VolFill(): src MRI is neither source nor target in LTA\n");
    return (NULL);
  }
  if (Gdiag_no > 0) printf("MRIvol2VolFill(): USF=%d, DoConserve=%d\n", UpsampleFactor, DoConserve);

  ltatmp = LTAcopy(lta, NULL);
  if (ltatmp->type != LINEAR_VOX_TO_VOX) ltatmp = LTAchangeType(ltatmp, LINEAR_VOX_TO_VOX);

  // Extract Vox2Vox
  if (LTAmriIsSource(ltatmp, src)) {
    if (Gdiag_no > 0) printf("MRIvol2VolFill(): not using inverse\n");
    v2v = ltatmp->xforms[0].m_L;
    vgtarg = ltatmp->xforms[lta->num_xforms - 1].dst;
  }
  else {
    // Invert the matrix if the LTA goes in the wrong direction
    if (Gdiag_no > 0) printf("MRIvol2VolFill(): using inverse\n");
    LTAfillInverse(ltatmp);
    v2v = ltatmp->inv_xforms[0].m_L;
    vgtarg = ltatmp->inv_xforms[lta->num_xforms - 1].dst;
  }

  if (UpsampleFactor > 0) {
    srcmus = MRImaskAndUpsample(src, mask, UpsampleFactor, nPad, DoConserve, &src2srcmus);
    // Recompute vox2vox
    vmusinv = MatrixInverse(src2srcmus->xforms[0].m_L, NULL);
    v2v = MatrixMultiply(v2v, vmusinv, v2v);
    MatrixFree(&vmusinv);
    LTAfree(&src2srcmus);
  }
  else
    srcmus = src;

  if (outfill == NULL) {
    outfill = MRIallocSequence(vgtarg.width, vgtarg.height, vgtarg.depth, MRI_FLOAT, src->nframes);
    useVolGeomToMRI(&vgtarg, outfill);
    outfill->tr = src->tr;
    outfill->te = src->te;
    outfill->flip_angle = src->flip_angle;
    outfill->ti = src->ti;
  }
  MRIclear(outfill);

  if (!DoConserve) {
    // This keeps track of the number of source voxels land in each target voxel
    hitmap = MRIallocSequence(vgtarg.width, vgtarg.height, vgtarg.depth, MRI_FLOAT, 1);
    useVolGeomToMRI(&vgtarg, hitmap);
    hitmap->tr = src->tr;
    hitmap->te = src->te;
    hitmap->flip_angle = src->flip_angle;
    hitmap->ti = src->ti;
  }

  crssrc = MatrixAlloc(4, 1, MATRIX_REAL);
  crssrc->rptr[4][1] = 1;
  crstarg = MatrixAlloc(4, 1, MATRIX_REAL);

  // Go through the source volume voxels
  for (c = 0; c < srcmus->width; c++) {
    for (r = 0; r < srcmus->height; r++) {
      for (s = 0; s < srcmus->depth; s++) {
        // Map source voxel to target voxel
        crssrc->rptr[1][1] = c;
        crssrc->rptr[2][1] = r;
        crssrc->rptr[3][1] = s;
        crstarg = MatrixMultiply(v2v, crssrc, crstarg);
        ct = nint(crstarg->rptr[1][1]);
        rt = nint(crstarg->rptr[2][1]);
        st = nint(crstarg->rptr[3][1]);
        if (ct < 0 || ct >= vgtarg.width) continue;
        if (rt < 0 || rt >= vgtarg.height) continue;
        if (st < 0 || st >= vgtarg.depth) continue;

        if (!DoConserve) {
          // Keep track of hits if not conserving
          nhits = MRIgetVoxVal(hitmap, ct, rt, st, 0);
          MRIsetVoxVal(hitmap, ct, rt, st, 0, nhits + 1);
        }

        for (f = 0; f < src->nframes; f++) {
          vsrc = MRIgetVoxVal(srcmus, c, r, s, f);
          vout = MRIgetVoxVal(outfill, ct, rt, st, f);
          MRIsetVoxVal(outfill, ct, rt, st, f, vsrc + vout);
        }
      }
    }
  }

  if (!DoConserve) {
    // Divide by number of hits if not conserving
    for (ct = 0; ct < vgtarg.width; ct++) {
      for (rt = 0; rt < vgtarg.height; rt++) {
        for (st = 0; st < vgtarg.depth; st++) {
          nhits = MRIgetVoxVal(hitmap, ct, rt, st, 0);
          if (nhits == 0) continue;
          for (f = 0; f < src->nframes; f++) {
            vout = MRIgetVoxVal(outfill, ct, rt, st, f);
            MRIsetVoxVal(outfill, ct, rt, st, f, vout / nhits);
          }
        }
      }
    }
    // MRIwrite(hitmap,"hitmap.mgh");
    MRIfree(&hitmap);
  }

  LTAfree(&ltatmp);
  if (srcmus != src) MRIfree(&srcmus);
  return (outfill);
}

/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIsegBoundary(MRI *seg)
  \brief Creates a segmentation boundary volume in which voxels
  whose neighbors are the same seg are set to 0.
*/
MRI *MRIsegBoundary(MRI *seg)
{
  MRI *boundary;
  int c, r, s, dc, dr, ds;
  int cseg, nseg, b;

  boundary = MRIclone(seg, NULL);

  for (c = 1; c < seg->width - 1; c++) {
    for (r = 1; r < seg->height - 1; r++) {
      for (s = 1; s < seg->depth - 1; s++) {
        cseg = (int)MRIgetVoxVal(seg, c, r, s, 0);
        if (cseg == 0) {
          MRIsetVoxVal(boundary, c, r, s, 0, 0);
          continue;
        }
        b = 0;
        for (dc = -1; dc < 2; dc++) {
          for (dr = -1; dr < 2; dr++) {
            for (ds = -1; ds < 2; ds++) {
              nseg = (int)MRIgetVoxVal(seg, c + dc, r + dr, s + ds, 0);
              if (cseg != nseg) {
                b = 1;  // It is a boundary voxel
                dc = 2;
                dr = 2;
                ds = 2;  // break from all three loops
              }
            }
          }
        }
        if (b == 0)
          MRIsetVoxVal(boundary, c, r, s, 0, 0);
        else
          MRIsetVoxVal(boundary, c, r, s, 0, cseg);
      }
    }
  }

  return (boundary);
}

/*!
  \fn MRI *MRIcrs(MRI *in, MRI *out)
  \brief Creates a 1-frame volume with the value
  at frame 0 equal to the slice number.
*/
MRI *MRIsliceNo(MRI *in, MRI *out)
{
  int c, r, s;
  if (out == NULL) {
    out = MRIalloc(in->width, in->height, in->depth, MRI_FLOAT);
    MRIcopyHeader(in, out);
  }

  for (c = 0; c < in->width; c++) {
    for (r = 0; r < in->height; r++) {
      for (s = 0; s < in->depth; s++) {
        MRIsetVoxVal(out, c, r, s, 0, s);
      }
    }
  }

  return (out);
}

/*!
  \fn MRI *MRIcrs(MRI *in, MRI *out)
  \brief Creates a 1-frame volume with the value
  at frame 0 equal to the voxel index, eg,
  sliceno*nrows*ncols + rowno*ncols + col.
*/
MRI *MRIindexNo(MRI *in, MRI *out)
{
  int c, r, s, index;
  if (out == NULL) {
    out = MRIalloc(in->width, in->height, in->depth, MRI_FLOAT);
    MRIcopyHeader(in, out);
  }

  index = 0;
  for (s = 0; s < in->depth; s++) {
    for (r = 0; r < in->height; r++) {
      for (c = 0; c < in->width; c++) {
        MRIsetVoxVal(out, c, r, s, 0, index);
        index++;
      }
    }
  }

  return (out);
}

/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIcrs(MRI *in, MRI *out)
  \brief Creates a 3-frame volume with the value
  at frame 0 equal to the column number, frame 1
  the row, and frame 2 the slice.
*/
MRI *MRIcrs(MRI *in, MRI *out)
{
  int c, r, s;

  if (out == NULL) {
    out = MRIallocSequence(in->width, in->height, in->depth, MRI_FLOAT, 3);
    MRIcopyHeader(in, out);
  }

  for (s = 0; s < in->depth; s++) {
    for (r = 0; r < in->height; r++) {
      for (c = 0; c < in->width; c++) {
        MRIsetVoxVal(out, c, r, s, 0, c);
        MRIsetVoxVal(out, c, r, s, 1, r);
        MRIsetVoxVal(out, c, r, s, 2, s);
      }
    }
  }

  return (out);
}

/*---------------------------------------------------------
  MRIsegStats() - computes statistics within a given
  segmentation. Returns the number of voxels in the
  segmentation.
  ---------------------------------------------------------*/
int MRIsegStats(MRI *seg, int segid, MRI *mri, int frame, float *min, float *max, float *range, float *mean, float *std)
{
  int id, nvoxels, r, c, s;
  double val, sum, sum2;

  *min = 0;
  *max = 0;
  sum = 0;
  sum2 = 0;
  nvoxels = 0;
  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        id = (int)MRIgetVoxVal(seg, c, r, s, 0);
        if (id != segid) {
          continue;
        }
        val = MRIgetVoxVal(mri, c, r, s, frame);
        nvoxels++;
        if (nvoxels == 1) {
          *min = val;
          *max = val;
        }
        if (*min > val) {
          *min = val;
        }
        if (*max < val) {
          *max = val;
        }
        sum += val;
        sum2 += (val * val);
      }
    }
  }

  *range = *max - *min;

  if (nvoxels != 0) {
    *mean = sum / nvoxels;
  }
  else {
    *mean = 0.0;
  }

  if (nvoxels > 1)
    *std = sqrt(((nvoxels) * (*mean) * (*mean) - 2 * (*mean) * sum + sum2) / (nvoxels - 1));
  else {
    *std = 0.0;
  }

  return (nvoxels);
}
/*------------------------------------------------------------*/
/*!
  \fn int MRIsegStatsRobust(MRI *seg, int segid, MRI *mri,int frame,
                      float *min, float *max, float *range,
                      float *mean, float *std, float Pct)
  \brief Computes stats based on the the middle 100-2*Pct values, ie,
         it trims Pct off the ends.
*/
int MRIsegStatsRobust(
    MRI *seg, int segid, MRI *mri, int frame, float *min, float *max, float *range, float *mean, float *std, float Pct)
{
  int id, nvoxels, r, c, s, k, m;
  double val, sum, sum2;
  float *vlist;

  *min = 0;
  *max = 0;
  *range = 0;
  *mean = 0;
  *std = 0;

  // Count number of voxels
  nvoxels = 0;
  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        id = (int)MRIgetVoxVal(seg, c, r, s, 0);
        if (id != segid) continue;
        nvoxels++;
      }
    }
  }
  if (nvoxels == 0) return (nvoxels);

  // Load voxels into an array
  vlist = (float *)calloc(sizeof(float), nvoxels);
  nvoxels = 0;
  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        id = (int)MRIgetVoxVal(seg, c, r, s, 0);
        if (id != segid) continue;
        vlist[nvoxels] = MRIgetVoxVal(mri, c, r, s, frame);
        nvoxels++;
      }
    }
  }
  // Sort the array
  qsort((void *)vlist, nvoxels, sizeof(float), compare_floats);

  // Compute stats excluding Pct of the values from each end
  sum = 0;
  sum2 = 0;
  m = 0;
  // printf("Robust Indices: %d %d\n",(int)nint(Pct*nvoxels/100.0),(int)nint((100-Pct)*nvoxels/100.0));
  for (k = 0; k < nvoxels; k++) {
    if (k < Pct * nvoxels / 100.0) continue;
    if (k > (100 - Pct) * nvoxels / 100.0) continue;
    val = vlist[k];
    if (m == 0) {
      *min = val;
      *max = val;
    }
    if (*min > val) *min = val;
    if (*max < val) *max = val;
    sum += val;
    sum2 += (val * val);
    m = m + 1;
  }

  *range = *max - *min;
  *mean = sum / m;
  if (m > 1)
    *std = sqrt(((m) * (*mean) * (*mean) - 2 * (*mean) * sum + sum2) / (m - 1));
  else
    *std = 0.0;

  free(vlist);
  vlist = NULL;
  return (m);
}
/*---------------------------------------------------------
  MRIsegFrameAvg() - computes the average time course withing the
  given segmentation. Returns the number of voxels in the
  segmentation. favg must be preallocated to number of
  frames. favg = (double *) calloc(sizeof(double),mri->nframes);
  ---------------------------------------------------------*/
int MRIsegFrameAvg(MRI *seg, int segid, MRI *mri, double *favg)
{
  int id, nvoxels, r, c, s, f;
  double val;

  /* zero it out */
  for (f = 0; f < mri->nframes; f++) {
    favg[f] = 0;
  }

  nvoxels = 0;
  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        id = (int)MRIgetVoxVal(seg, c, r, s, 0);
        if (id != segid) {
          continue;
        }
        for (f = 0; f < mri->nframes; f++) {
          val = MRIgetVoxVal(mri, c, r, s, f);
          favg[f] += val;
        }
        nvoxels++;
      }
    }
  }

  if (nvoxels != 0)
    for (f = 0; f < mri->nframes; f++) {
      favg[f] /= nvoxels;
    }

  return (nvoxels);
}

MRI *MRImask_with_T2_and_aparc_aseg(
    MRI *mri_src, MRI *mri_dst, MRI *mri_T2, MRI *mri_aparc_aseg, float T2_thresh, int mm_from_exterior)
{
  int x, y, z, nremoved, i;
  MRI *mri_bright, *mri_mask, *mri_tmp = NULL;

  mri_mask = MRIbinarize(mri_T2, NULL, T2_thresh, 255, 0);
  mri_bright = MRIcopy(mri_mask, NULL);

  if (mri_aparc_aseg)  // use T2 and aparc+aseg to remove non-brain stuff
  {
    MRIbinarize(mri_aparc_aseg, mri_aparc_aseg, 1, 0, 255);
    MRIdilate(mri_aparc_aseg, mri_aparc_aseg);
    MRInot(mri_aparc_aseg, mri_aparc_aseg);  // background now on, foreground off
    GetLargestCC6(mri_aparc_aseg);           // remove disconnected background components
    MRIand(mri_mask, mri_aparc_aseg, mri_mask, 1);
    MRIopenN(mri_mask, mri_mask, 3);  // third order open will remove thin chains of bright T2 that are in the interior
  }
  else  // just use T2
  {
    GetLargestCC6(mri_mask);
  }

  MRInot(mri_mask, mri_mask);  // 0 now means mask it out and 1 means retain it

  for (i = nremoved = 0; i < mm_from_exterior; i++) {
    mri_tmp = MRIcopy(mri_mask, mri_tmp);
    for (x = 0; x < mri_mask->width; x++)
      for (y = 0; y < mri_mask->height; y++)
        for (z = 0; z < mri_mask->depth; z++) {
          if (x == Gx && y == Gy && z == Gz) DiagBreak();
          if (MRIgetVoxVal(mri_mask, x, y, z, 0) == 0)  // already in the mask
            continue;
          if (MRIgetVoxVal(mri_aparc_aseg, x, y, z, 0) == 0)  // too close to brain
            continue;
          if (MRIgetVoxVal(mri_T2, x, y, z, 0) >= T2_thresh)  // bright in the T2
          {
            if (MRIneighborsOff(mri_mask, x, y, z, 1) > 0)  // touching the existing mask
            {
              if (x == Gx && y == Gy && z == Gz) DiagBreak();
              MRIsetVoxVal(mri_tmp, x, y, z, 0, 0);  // add this voxel to the mask
              nremoved++;
            }
          }
        }

    printf("%d T2-bright exterior voxels removed\n", nremoved);
    MRIcopy(mri_tmp, mri_mask);
  }

  mri_dst = MRImask(mri_dst, mri_mask, mri_dst, 0, 0);  // if mask == 0, then set dst as 0
  MRIfree(&mri_bright);
  MRIfree(&mri_mask);
  MRIfree(&mri_tmp);
  return (mri_dst);
}

/* ----------------------------------------------------------*/
/*!
  \fn int *MRIsegmentationList(MRI *seg, int *pListLength)
  \brief Extracts a list of unique segmentation IDs from the volume.
    Includes 0.
*/
int *MRIsegmentationList(MRI *seg, int *pListLength)
{
  int c, r, s, n, nvox;
  int *list, *voxlist;

  nvox = seg->width * seg->height * seg->depth;
  voxlist = (int *)calloc(nvox, sizeof(int));
  n = 0;
  for (s = 0; s < seg->depth; s++) {
    for (c = 0; c < seg->width; c++) {
      for (r = 0; r < seg->height; r++) {
        voxlist[n] = MRIgetVoxVal(seg, c, r, s, 0);
        n++;
      }
    }
  }
  list = unqiue_int_list(voxlist, nvox, pListLength);
  printf("MRIsegmentationList(): found %d unique segmentations\n", *pListLength);
  if (Gdiag_no > 0)
    for (n = 0; n < *pListLength; n++) printf("%2d %5d\n", n, list[n]);
  free(voxlist);
  return (list);
}
/* ----------------------------------------------------------*/
/*!
  \fn MATRIX *BuildGTM0(MRI *seg, MRI *mask,double cFWHM, double rFWHM,
                  double sFWHM, MATRIX *X)
  \brief Builds the Geometric Transfer Matrix (GTM) without taking
  into account the partial volume fraction within a voxel (thus the 0
  in the name).
  \param seg - segmentation. segid=0 will be ignored
  \param mask - binary brain mask
  \param {c,r,s}FWHM - PSF FWHM in mm for each dimension
  \param X - GTM nmask by nsegs. Segs are sorted in numerical order.
     The order of the voxels is row-fastest, then col, then slice,
     making it consistent with matlab. This order must be used when
     creating the y-matrix.
*/
MATRIX *BuildGTM0(MRI *seg, MRI *mask, double cFWHM, double rFWHM, double sFWHM, MATRIX *X)
{
  int c, r, s, nmask, nsegs, nthseg, mthseg, segid, *segidlist, has0;
  double cStd, rStd, sStd, val;
  MRI *roimask = NULL, *roimasksm = NULL;

  cStd = cFWHM / sqrt(log(256.0));
  rStd = rFWHM / sqrt(log(256.0));
  sStd = sFWHM / sqrt(log(256.0));

  segidlist = MRIsegIdList(seg, &nsegs, 0);

  has0 = 0;
  for (nthseg = 0; nthseg < nsegs; nthseg++)
    if (segidlist[nthseg] == 0) has0 = 1;

  // Count number of voxels in the mask
  nmask = 0;
  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
        nmask++;
      }
    }
  }
  if (Gdiag_no > 0) printf("BuildGTM0(): nmask = %d, nsegs = %d\n", nmask, nsegs);

  if (X == NULL) X = MatrixAlloc(nmask, nsegs - has0, MATRIX_REAL);
  if (X->rows != nmask || X->cols != nsegs - has0) {
    printf("ERROR: BuildGTM0(): X dim mismatch\n");
    return (NULL);
  }

  roimask = MRIconst(seg->width, seg->height, seg->depth, 1, 0.0, NULL);
  MRIcopyHeader(seg, roimask);

  mthseg = 0;
  for (nthseg = 0; nthseg < nsegs; nthseg++) {
    segid = segidlist[nthseg];
    if (segid == 0) continue;

    if (Gdiag_no > 0) {
      printf("BuildGTM0(): #@# %3d/%d %3d ---\n", mthseg, nsegs - has0, segid);
      fflush(stdout);
    }

    // Create a mask of the seg
    for (c = 0; c < seg->width; c++) {
      for (r = 0; r < seg->height; r++) {
        for (s = 0; s < seg->depth; s++) {
          if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
          if (MRIgetVoxVal(seg, c, r, s, 0) == segid)
            val = 1;
          else
            val = 0;
          MRIsetVoxVal(roimask, c, r, s, 0, val);
        }
      }
    }
    // Smooth the mask
    roimasksm = MRIgaussianSmoothNI(roimask, cStd, rStd, sStd, roimasksm);
    // Fill X
    // Creating X in this order makes it consistent with matlab
    // Note: y must be ordered in the same way.
    nmask = 0;
    for (s = 0; s < seg->depth; s++) {
      for (c = 0; c < seg->width; c++) {
        for (r = 0; r < seg->height; r++) {
          if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
          X->rptr[nmask + 1][mthseg + 1] = MRIgetVoxVal(roimasksm, c, r, s, 0);
          nmask++;
        }
      }
    }
    mthseg++;
  }  // seg

  MRIfree(&roimask);
  MRIfree(&roimasksm);
  free(segidlist);

  return (X);
}

/*!
  \fn MRI *MRIfisherTransform(MRI *rho, MRI *out)
  \brief Computes fisher transform. The input should be a correlation
  coefficient with values between -1 and +1
*/
MRI *MRIfisherTransform(MRI *rho, MRI *mask, MRI *out)
{
  int c, r, s, f;
  double v, ft;

  if (out == NULL) {
    out = MRIallocSequence(rho->width, rho->height, rho->depth, MRI_FLOAT, rho->nframes);
    if (out == NULL) return (NULL);
    MRIcopyHeader(rho, out);
  }

  for (s = 0; s < rho->depth; s++) {
    for (r = 0; r < rho->height; r++) {
      for (c = 0; c < rho->width; c++) {
        if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) {
          for (f = 0; f < rho->nframes; f++) MRIsetVoxVal(out, c, r, s, f, 0.0);
          continue;
        }
        for (f = 0; f < rho->nframes; f++) {
          v = MRIgetVoxVal(rho, c, r, s, f);
          ft = .5 * log((1 + v) / (1 - v));
          MRIsetVoxVal(out, c, r, s, f, ft);
        }
      }
    }
  }
  return (out);
}

/*!
  \fn MRI *MRIbinarizeMatch(MRI *seg, int match, int frame, MRI *out)
  \brief Binarizes a volume based on the voxels values that match the match value.
*/
MRI *MRIbinarizeMatch(MRI *seg, int *MatchList, int nList, int frame, MRI *out)
{
  int c, r, s, m, n;

  if (out == NULL) {
    out = MRIalloc(seg->width, seg->height, seg->depth, MRI_INT);
    MRIcopyHeader(seg, out);
  }
  MRIclear(out);

  for (s = 0; s < seg->depth; s++) {
    for (c = 0; c < seg->width; c++) {
      for (r = 0; r < seg->height; r++) {
        m = MRIgetVoxVal(seg, c, r, s, frame);
        for (n = 0; n < nList; n++) {
          if (m == MatchList[n]) {
            MRIsetVoxVal(out, c, r, s, 0, 1);
            break;
          }
        }
      }
    }
  }
  return (out);
}
/*
  \fn MRI *MRIhiresSeg(MRI *aseg, MRIS *lhw, MRIS *lhp, MRIS *rhw, MRIS *rhp, int USF, LTA **aseg2hrseg)
  \brief Creates a high-resolution (upsampled) segmentation given the
  aseg and surfaces. The USF is the upsampling factor. The result is
  upsampled and the FoV is reduced to the bare minimum so the number
  of voxels in a dimension will not necessarily be USF times the
  original number. The subcortical structures are the upsampled
  versions of the low-res aseg (ie, no new information is
  created). However, cortex benefits from the upsampling. aseg2hrseg
  is the transform between the aseg space and the hires space (they
  share a scanner RAS space). If aseg=NULL, then the VOL_GEOM from lhw
  is used. If USF=-1, then no FoV reduction is done. The surfaces can
  be NULL.  This is really meant for the aseg to be the aseg.mgz with
  cortex = 3,42 and cerebralwm = 2,41.
*/
MRI *MRIhiresSeg(MRI *aseg, MRIS *lhw, MRIS *lhp, MRIS *rhw, MRIS *rhp, int USF, LTA **aseg2hrseg)
{
  MRI *asegus, *lhwvol, *rhwvol, *lhpvol, *rhpvol, *seg;
  int c, r, s, asegv, lhwv, rhwv, lhpv, rhpv, segv, lhRibbon, rhRibbon, Ribbon, SubCort;
  int nPad = 2;
  int HasXCSF, UnknownFill;

  asegus = NULL;
  if (aseg != NULL) {
    if (aseg->type == MRI_UCHAR) {
      printf("ERROR: MRIhiresSeg(): aseg cannot be uchar\n");
      return (NULL);
    }
    if (USF > 0)
      asegus = MRImaskAndUpsample(aseg, aseg, USF, nPad, 0, aseg2hrseg);
    else {
      asegus = aseg;
      *aseg2hrseg = TransformRegDat2LTA(aseg, aseg, NULL);  // Identity
    }
  }
  else {
    // aseg might be null for testing purposes
    if (lhw->vg.valid == 0) {
      printf("ERROR: MRIhiresSeg(): if aseg=NULL, then surface volume geometry must be valid\n");
      return (NULL);
    }
    printf("Info: MRIhiresSeg(): aseg is NULL\n");
    aseg = MRIalloc(lhw->vg.width, lhw->vg.height, lhw->vg.depth, MRI_UCHAR);
    useVolGeomToMRI(&lhw->vg, aseg);
    asegus = MRIupsampleN(aseg, NULL, abs(USF));
    MRIfree(&aseg);
  }
  seg = MRIcopy(asegus, NULL);
  MRIcopyHeader(asegus, seg);

  // Check whether the seg has an extracerebral CSF segmentation
  HasXCSF = MRIcountMatches(aseg, CSF_ExtraCerebral, 0, aseg);
  if (HasXCSF > 0)
    UnknownFill = CSF_ExtraCerebral;
  else
    UnknownFill = 0;
  printf("  MRIhiresSeg(): filling unknowns with %d\n", UnknownFill);
  fflush(stdout);

  if (lhw) {
    // printf("lhw -------------\n");
    lhwvol = MRIcopy(asegus, NULL);
    MRIcopyHeader(asegus, lhwvol);
    MRISfillInterior(lhw, 0, lhwvol);
  }
  if (lhp) {
    // printf("lhp -------------\n");
    lhpvol = MRIcopy(asegus, NULL);
    MRIcopyHeader(asegus, lhpvol);
    MRISfillInterior(lhp, 0, lhpvol);
  }
  if (rhw) {
    // printf("rhw -------------\n");
    rhwvol = MRIcopy(asegus, NULL);
    MRIcopyHeader(asegus, rhwvol);
    MRISfillInterior(rhw, 0, rhwvol);
  }
  if (rhp) {
    // printf("rhp -------------\n");
    rhpvol = MRIcopy(asegus, NULL);
    MRIcopyHeader(asegus, rhpvol);
    MRISfillInterior(rhp, 0, rhpvol);
  }

  // stop compiler warnings
  lhwv = rhwv = lhpv = rhpv = 0;
  segv = 0;
  asegv = 0;

  if (Gdiag_no > 0) printf("Starting seg fill\n");
  for (c = 0; c < asegus->width; c++) {
    for (r = 0; r < asegus->height; r++) {
      for (s = 0; s < asegus->depth; s++) {
        if (aseg) asegv = MRIgetVoxVal(asegus, c, r, s, 0);
        if (lhw) lhwv = MRIgetVoxVal(lhwvol, c, r, s, 0);
        if (rhw) rhwv = MRIgetVoxVal(rhwvol, c, r, s, 0);
        if (lhp) lhpv = MRIgetVoxVal(lhpvol, c, r, s, 0);
        if (rhp) rhpv = MRIgetVoxVal(rhpvol, c, r, s, 0);

        // Check if voxel is in the "true" ribbon
        lhRibbon = 0;
        if ((lhpv && !lhwv)) lhRibbon = 1;
        rhRibbon = 0;
        if ((rhpv && !rhwv)) rhRibbon = 1;
        Ribbon = lhRibbon || rhRibbon;

        /* SubCort=1 if aseg says voxel is neither cortex nor cerebral
        WM nor background nor extracerebral CSF. Note that SubCort=1
        for other structures as well (head eyes, etc).  */
        SubCort = 0;
        if (asegv != Left_Cerebral_Cortex && asegv != Left_Cerebral_White_Matter && asegv != Right_Cerebral_Cortex &&
            asegv != Right_Cerebral_White_Matter && asegv != 0 && asegv != CSF_ExtraCerebral && aseg)
          SubCort = 1;

        if (SubCort)
          segv = asegv;  // subcort rules
        else if (Ribbon) {
          // Voxel is in the "true" ribbon but not in a subcortical structure
          // Aseg could say it is in xcCSF, but override that with ribbon
          // What if something other than aseg.mgz was passed? Seg might not
          // have X_Left_Cerebral_Cortex label. This is ok for computing PVF
          // but not if you are expecting the labels to be correct.
          if (lhRibbon) segv = Left_Cerebral_Cortex;
          if (rhRibbon) segv = Right_Cerebral_Cortex;
        }
        else {
          /* To get here, it cannot be in the true ribbon so, if the
             aseg says it is CorticalGM, that is wrong. The question
             is, what is right? Probably WM xCSF , but ideally, this should
             be set to the seg of nearby voxels.*/
          if (asegv == Left_Cerebral_Cortex) {
            if (lhpv)
              segv = Left_Cerebral_White_Matter;
            else
              segv = UnknownFill;
          }
          else if (asegv == Right_Cerebral_Cortex) {
            if (rhpv)
              segv = Right_Cerebral_White_Matter;
            else
              segv = UnknownFill;
          }
          else {
            // To get here aseg can only be CerebralWM, CSF_ExtraCerebral, Head_ExtraCerebral
            // or something else outside of the brain.
            if (asegv != Left_Cerebral_White_Matter && asegv != Right_Cerebral_White_Matter &&
                asegv != CSF_ExtraCerebral && asegv != Head_ExtraCerebral && asegv != 0 && aseg)
              if (Gdiag > 0)
                printf("WARNING: MRIhiresSeg(): voxel %d %d %d is %d, expecting WM, xcCSF, or head\n", c, r, s, asegv);
            segv = asegv;
          }
        }
        MRIsetVoxVal(seg, c, r, s, 0, segv);
      }
    }
  }
  if (asegus != aseg) MRIfree(&asegus);
  if (lhw) MRIfree(&lhwvol);
  if (rhw) MRIfree(&rhwvol);
  if (lhp) MRIfree(&lhpvol);
  if (rhp) MRIfree(&rhpvol);
  return (seg);
}
/*
  \fn MRI *MRIpartialVolumeFraction(LTA *seg2vol, MRI *seg, double resmm, COLOR_TABLE *ct, MRI *pvf)
  \brief Creates PVF maps for each tissue type in the color table. seg can be the aseg but
  is often a highres seg created from the aseg and surfaces (see MRIhiresSeg())
  USF is the upsample factor. If using a seg from MRIhiresSeg(), which has its
  own USF, then USF here can be set to 1. See also MRIpartialVolumeFractionAS().
  The return is an array of MRIs, one for each tissue type. The output volume
  is that of the dst volume geometry in seg2vol.
 */
MRI *MRIpartialVolumeFraction(LTA *seg2vol, MRI *seg, double resmm, COLOR_TABLE *ct, MRI *pvf)
{
  MRI *ttseg;
  int nTT, nsegs, *segidlist;
  VOL_GEOM *vg;

  if (ct->ctabTissueType == NULL) {
    printf("ERROR: MRIpartialVolumeFraction(): color table does not have tissue type ctab\n");
    return (NULL);
  }
  nTT = ct->ctabTissueType->nentries - 1;  // -1 to exclude background

  if (!LTAmriIsSource(seg2vol, seg) && !LTAmriIsTarget(seg2vol, seg)) {
    printf("ERROR: MRIpartialVolumeFraction(): seg MRI is neither source nor target in LTA\n");
    return (NULL);
  }
  vg = &(seg2vol->xforms[0].dst);

  if (pvf == NULL) {
    pvf = MRIallocSequence(vg->width, vg->height, vg->depth, MRI_FLOAT, nTT);
    if (pvf == NULL) {
      printf("ERROR: MRIpartialVolumeFraction(): could not alloc\n");
      return (NULL);
    }
    useVolGeomToMRI(vg, pvf);
    MRIcopyPulseParameters(seg, pvf);
  }
  if (pvf->width != vg->width || pvf->height != vg->height || pvf->depth != vg->depth || pvf->nframes != nTT) {
    printf("ERROR: MRIpartialVolumeFraction(): dimension mismatch\n");
    return (NULL);
  }

  // Create a tissue type segmentation from the seg
  ttseg = MRIseg2TissueType(seg, ct, NULL);
  if (ttseg == NULL) return (NULL);

  segidlist = MRIsegIdListNot0(ttseg, &nsegs, 0);  // nsegs=nTT
  pvf = MRIseg2SegPVF(ttseg, seg2vol, resmm, segidlist, nTT, NULL, 1, NULL, pvf);
  MRIseg2SegPVF(NULL, NULL, 0, NULL, 0, NULL, -1, NULL, NULL);  // clear cache
  free(segidlist);
  MRIfree(&ttseg);
  return (pvf);
  /*---------------------------------------------*/
  // output volume geometry
  // vmult = (seg->xsize*seg->ysize*seg->zsize)/(vol->xsize*vol->ysize*vol->zsize);
  // Go through each tissue type
  // for(tt = 0; tt < nTT; tt++){
  // binarize tissue type map
  // ttbin = MRIbinarizeMatch(ttseg, tt+1, 0, ttbin);
  // compute pvf based on number of seg voxels that fall into output vol vox
  // pvf[tt] = MRIvol2VolFill(ttbin, NULL, seg2vol, 1, 0, pvf[tt]);//USF=1 always here
  // if(pvf[tt]==NULL) return(NULL);
  // Better to turn off conserving in vol2volFill than to scale. The simple scaling
  // below creates a situation in which voxels in the middle of WM do not have
  // a PVF=1 because the number of highres voxels that land in a lowres voxel
  // is not constant.
  // Scale factor for mapping to a different voxel size
  // MRImultiplyConst(pvf[tt], vmult, pvf[tt]);
  //}
  // MRIfree(&ttseg);
  // MRIfree(&ttbin);
  // return(pvf);
}

/*
  \fn MRI **MRIpartialVolumeFractionAS(LTA *aseg2vol, MRI *aseg, MRIS *lhw, MRIS *lhp,
                                 MRIS *rhw, MRIS *rhp, int USF, COLOR_TABLE *ct, **pvf)
  \brief Creates PVF maps for each tissue type in the color table
  given the aseg and surfaces.  The return is an array of MRIs, one
  for each tissue type. The output volume is that of the dst volume
  geometry in seg2vol. aseg is usually the aseg.mgz. USF is the
  upsample factor, usually set to 2 or 3. This function calls
  MRIhiresSeg() then calls MRIpartialVolumeFractionAS().
 */
MRI *MRIpartialVolumeFractionAS(LTA *aseg2vol,
                                MRI *aseg,
                                MRIS *lhw,
                                MRIS *lhp,
                                MRIS *rhw,
                                MRIS *rhp,
                                int USF,
                                double resmm,
                                COLOR_TABLE *ct,
                                MRI *pvf)
{
  MRI *hrseg;
  LTA *aseg2hrseg, *hrseg2aseg, *hrseg2vol, *ltaArray[2];

  // Create a high resolution segmentation
  hrseg = MRIhiresSeg(aseg, lhw, lhp, rhw, rhp, USF, &aseg2hrseg);
  if (hrseg == NULL) return (NULL);
  hrseg2aseg = LTAinvert(aseg2hrseg, NULL);

  // Compute transform from high res to output volume
  ltaArray[0] = hrseg2aseg;
  ltaArray[1] = aseg2vol;
  hrseg2vol = LTAconcat(ltaArray, 2, 1);  // figures out inversions
  if (hrseg2vol == NULL) return (NULL);

  pvf = MRIpartialVolumeFraction(hrseg2vol, hrseg, resmm, ct, pvf);  // USF=1 here

  MRIfree(&hrseg);
  LTAfree(&hrseg2aseg);
  LTAfree(&hrseg2vol);
  return (pvf);
}

/*
  \fn int MRIcountMatches(MRI *seg, int MatchVal, int frame, MRI *mask)
  \brief Count the number of times there is a match in seg at the given frame
  within the mask (mask can be null).
 */
int MRIcountMatches(const MRI *seg, const int MatchVal, const int frame, const MRI *mask)
{
  int nMatches = 0;
  int c;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction(+ : nMatches)
#endif
  for (c = 0; c < seg->width; c++) {
    ROMP_PFLB_begin
    
    int r, s;
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
        if (MRIgetVoxVal(seg, c, r, s, frame) == MatchVal) nMatches++;
      }
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (nMatches);
}

/*
  \fn MRI *MRIaddExtraCerebralCSF(MRI *seg, int nDil, MRI *out)
  \brief Adds the CSF_ExtraCerebral segmentation to seg by dilating the segmentation
  then assigning CSF_ExtraCerebral to the new voxels. Dilates by nDil. If CSF_ExtraCerebral
  already exists in the seg, then nothing is done. If nDil<=0, then all voxels outside
  the segmentation are used. Note: this makes no attempt whatsoever to do a true segmentation
  of extra-cerebral CSF!!! The seg (or out if non-NULL) must not be UCHAR.
*/
MRI *MRIaddExtraCerebralCSF(MRI *seg, int nDil, MRI *out)
{
  int c, r, s, n, nxcsf;
  MRI *mask = NULL;

  if (seg->type == MRI_UCHAR && out != NULL && out->type == MRI_UCHAR) {
    printf("ERROR: MRIaddExtraCerebralCSF(): passed seg/out is of type UCHAR\n");
    return (NULL);
  }
  out = MRIcopy(seg, out);
  if (out == NULL) return (NULL);
  MRIcopyHeader(seg, out);
  MRIcopyPulseParameters(seg, out);

  // Check whether CSF_ExtraCerebral already exists
  n = MRIcountMatches(seg, CSF_ExtraCerebral, 0, NULL);
  if (n > 0) {
    if (Gdiag_no > 0)
      printf("MRIaddExtraCerebralCSF(): %d CSF_ExtraCerebral voxels already exist, not adding any more\n", n);
    return (out);
  }

  if (Gdiag_no > 0) printf("MRIaddExtraCerebralCSF(): nDil = %d %d\n", nDil, CSF_ExtraCerebral);
  if (nDil >= 1) {
    mask = MRIdilate(seg, NULL);
    for (n = 1; n < nDil; n++) MRIdilate(mask, mask);
  }

  nxcsf = 0;
  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        if (mask && MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;
        if (MRIgetVoxVal(seg, c, r, s, 0) > 0.5) continue;
        MRIsetVoxVal(out, c, r, s, 0, CSF_ExtraCerebral);
        nxcsf++;
      }  // slice
    }    // row
  }      // col
  if (Gdiag_no > 0) printf("MRIaddExtraCerebralCSF(): Added %d CSF_ExtraCerebral voxels\n", nxcsf);
  if (mask) MRIfree(&mask);
  return (out);
}
/*
\fn COLOR_TABLE *CTABpruneCTab(const COLOR_TABLE *ct0, MRI *seg)
\brief Creates a new CTAB with only the segments in seg
*/
COLOR_TABLE *CTABpruneCTab(const COLOR_TABLE *ct0, MRI *seg)
{
  COLOR_TABLE *ct;
  int *segidlist, nsegs, segid, n;

  segidlist = MRIsegIdList(seg, &nsegs, 0);  // list of segs and nsegs

  ct = CTABalloc(segidlist[nsegs - 1] + 1);
  for (n = 0; n < ct->nentries; n++) {  // start by setting all to NULL
    free(ct->entries[n]);
    ct->entries[n] = NULL;
  }
  strcpy(ct->TissueTypeSchema,ct0->TissueTypeSchema);

  for (n = 0; n < nsegs; n++) {
    segid = segidlist[n];
    if (ct0->entries[segid] == NULL) {
      printf("ERROR: CTABpruneCTab(): ctab does not have segid %d\n", segid);
      return (NULL);
    }
    ct->entries[segid] = (CTE *)calloc(sizeof(CTE), 1);
    memcpy(ct->entries[segid], ct0->entries[segid], sizeof(CTE));
  }

  if(ct0->ctabTissueType) 
    ct->ctabTissueType = CTABdeepCopy(ct0->ctabTissueType);

  free(segidlist);
  return (ct);
}

/*
  \fn MRI *MRIannot2CorticalSeg(MRI *seg, MRIS *lhw, MRIS *lhp, MRIS *rhw, MRIS *rhp, LTA *anat2seg, MRI *ctxseg)
  \brief Creates a segmentation of the cortical labels
  (X_Cerebral_Cortex) found in seg based upon the annotation of the
  nearest cortical vertex. For unknown areas, the segmentation is
  given CSF_ExtraCerebral if that segno already exists in seg,
  otherwise it is give 0. The annotion is expected to be in the pial
  surfaces. anat2seg is the registration between the
  surface/anatomical space and the segmentation space. If they share a
  space, then just use NULL. The surface space is obtained from
  lhw->vg. This function basically is what is done when creating
  aparc+aseg.mgz. It is recommended that MRISsetPialUnknownToWhite()
  be run on the pial surfaces before using this function. The output
  will be the same size as seg.
 */
MRI *MRIannot2CorticalSeg(MRI *seg, MRIS *lhw, MRIS *lhp, MRIS *rhw, MRIS *rhp, LTA *anat2seg, MRI *ctxseg)
{
  MATRIX *AnatVox2SurfRAS, *SegVox2SurfRAS;
  float hashres = 16;
  MHT *lhw_hash = NULL, *rhw_hash = NULL, *lhp_hash = NULL, *rhp_hash = NULL;
  LTA *lta;
  MRI *anat;
  int c, nunknown;
  int HasXCSF, UnknownFill;

  if (lhw->vg.valid != 1) {
    printf("ERROR: MRIannot2CorticalSeg(): lhw does not have a valid geometry\n");
    return (NULL);
  }

  if (ctxseg == NULL) {
    ctxseg = MRIallocSequence(seg->width, seg->height, seg->depth, MRI_INT, 1);
    MRIcopyHeader(seg, ctxseg);
    MRIcopyPulseParameters(seg, ctxseg);
  }
  if (MRIdimMismatch(seg, ctxseg, 0)) {
    printf("ERROR: MRIannot2CorticalSeg(): dimension mismatch\n");
    return (NULL);
  }

  // Create an MRI for the anatomical voume the surfaces were generated from
  anat = MRIallocFromVolGeom(&(lhw->vg), MRI_INT, 1, 1);

  // Compute an LTA that maps from the anatomical to the segmentation
  if (anat2seg == NULL)
    lta = TransformRegDat2LTA(anat, seg, NULL);
  else {
    if (LTAmriIsTarget(anat2seg, seg))
      lta = LTAcopy(anat2seg, NULL);
    else
      lta = LTAinvert(anat2seg, NULL);
  }
  if (lta->type != LINEAR_VOX_TO_VOX) LTAchangeType(lta, LINEAR_VOX_TO_VOX);
  LTAfillInverse(lta);

  // Anatomical Vox to Surface RAS
  AnatVox2SurfRAS = MRIxfmCRS2XYZtkreg(anat);
  // Segmentation Vox to Surface RAS
  SegVox2SurfRAS = MatrixMultiplyD(AnatVox2SurfRAS, lta->inv_xforms[0].m_L, NULL);

  // Create the hash for faster service
  lhw_hash = MHTcreateVertexTable_Resolution(lhw, CURRENT_VERTICES, hashres);
  lhp_hash = MHTcreateVertexTable_Resolution(lhp, CURRENT_VERTICES, hashres);
  rhw_hash = MHTcreateVertexTable_Resolution(rhw, CURRENT_VERTICES, hashres);
  rhp_hash = MHTcreateVertexTable_Resolution(rhp, CURRENT_VERTICES, hashres);

  // Check whether the seg has an extracerebral CSF segmentation
  HasXCSF = MRIcountMatches(seg, CSF_ExtraCerebral, 0, seg);
  if (HasXCSF > 0)
    UnknownFill = CSF_ExtraCerebral;
  else
    UnknownFill = 0;

  printf("  MRIannot2CorticalSeg(): looping over volume\n");
  fflush(stdout);
  nunknown = 0;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  printf("     nthreads = %d\n", omp_get_max_threads());
  #pragma omp parallel for if_ROMP(experimental) reduction(+ : nunknown)
#endif
  for (c = 0; c < seg->width; c++) {
    ROMP_PFLB_begin
    
    int r, s, asegv, annot, annotid, vtxno, wvtxno, pvtxno, segv;
    // int wmval;
    struct { float x,y,z; } vtx;
    MATRIX *RAS = NULL, *CRS = NULL;
    float wdw, pdw;
    MRIS *wsurf, *psurf;
    MHT *whash = NULL, *phash = NULL;
    CRS = MatrixAlloc(4, 1, MATRIX_REAL);
    CRS->rptr[4][1] = 1;
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        asegv = MRIgetVoxVal(seg, c, r, s, 0);

        if (asegv != Left_Cerebral_Cortex && asegv != Right_Cerebral_Cortex) {
          // If not in the ribbon, just copy seg to output and continue
          MRIsetVoxVal(ctxseg, c, r, s, 0, asegv);
          continue;
        }

        if (asegv == Left_Cerebral_Cortex) {
          wsurf = lhw;
          whash = lhw_hash;
          phash = lhp_hash;
          psurf = lhp;
          // wmval = Left_Cerebral_White_Matter;
        }
        else {
          wsurf = rhw;
          whash = rhw_hash;
          phash = rhp_hash;
          psurf = rhp;
          // wmval = Right_Cerebral_White_Matter;
        }

        // Compute location of voxel in surface space
        CRS->rptr[1][1] = c;
        CRS->rptr[2][1] = r;
        CRS->rptr[3][1] = s;
        RAS = MatrixMultiply(SegVox2SurfRAS, CRS, RAS);
        vtx.x = RAS->rptr[1][1];
        vtx.y = RAS->rptr[2][1];
        vtx.z = RAS->rptr[3][1];

        // Find closest white surface vertex and compute distance
        wvtxno = MHTfindClosestVertexNoXYZ(whash, wsurf, vtx.x, vtx.y, vtx.z, &wdw);
        if (wvtxno < 0) wvtxno = MRISfindClosestVertex(wsurf, vtx.x, vtx.y, vtx.z, &wdw, CURRENT_VERTICES);

        // Find closest pial surface vertex and compute distance
        pvtxno = MHTfindClosestVertexNoXYZ(phash, psurf, vtx.x, vtx.y, vtx.z, &pdw);
        if (pvtxno < 0) pvtxno = MRISfindClosestVertex(psurf, vtx.x, vtx.y, vtx.z, &pdw, CURRENT_VERTICES);

        // Use the vertex that is closest
        if (wdw <= pdw)
          vtxno = wvtxno;
        else
          vtxno = pvtxno;

        // From the surface annotation, get the annotation number
        annot = psurf->vertices[vtxno].annotation;
        // Convert annotation number to an entry number
        CTABfindAnnotation(psurf->ct, annot, &annotid);
        // Set segmentation to entry number + idbase
        if (annotid != -1)
          segv = annotid + psurf->ct->idbase;
        else {
          segv = UnknownFill;  // no annotation present
          nunknown++;
        }
        MRIsetVoxVal(ctxseg, c, r, s, 0, segv);
      }
    }
    MatrixFree(&CRS);
    if (RAS) MatrixFree(&RAS);
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  printf("  MRIannot2CorticalSeg(): found %d unknown, filled with %d\n", nunknown, UnknownFill);
  fflush(stdout);

  MHTfree(&lhw_hash);
  MHTfree(&rhw_hash);
  MHTfree(&rhp_hash);
  MHTfree(&lhp_hash);
  MatrixFree(&SegVox2SurfRAS);
  MatrixFree(&AnatVox2SurfRAS);
  MRIfree(&anat);
  LTAfree(&lta);

  return (ctxseg);
}

/*
  \fn MRI *MRIannot2CerebralWMSeg(MRI *seg, MRIS *lhw, MRIS *rhw, double DistThresh, LTA *anat2seg, MRI *wmseg)

  \brief Creates a segmentation of the cerebral WM
  (X_Cerebral_White_Matter) found in seg based upon the annotation of
  the nearest white vertex within DistThresh mm of the cortex. If it
  is beyond the distance threshold then 5001 or 5002 is used. If
  DistThresh is negative, then no distance threshold is used.
  anat2seg is the registration between the surface/anatomical space
  and the segmentation space. If they share a RAS space, then just use
  NULL. The surface space is obtained from lhw->vg. The output will be
  the same size as seg. Note that a voxel in seg must be labeled as
  X_Cerebral_White_Matter for it to receive a new label. This means
  that if you want hypointensities or CC to be re-labeled, you must
  change the hypointensities label to that of X_Cerebral_White_Matter
  (see MRIunsegmentWM()). This function basically is what is done when
  creating wmparc.mgz.
 */
MRI *MRIannot2CerebralWMSeg(MRI *seg, MRIS *lhw, MRIS *rhw, double DistThresh, LTA *anat2seg, MRI *wmseg)
{
  MATRIX *AnatVox2SurfRAS, *SegVox2SurfRAS;
  float hashres = 16;
  MHT *lhw_hash = NULL, *rhw_hash = NULL;
  LTA *lta;
  MRI *anat;
  int c;

  if (lhw->vg.valid != 1) {
    printf("ERROR: MRIannot2CerebralWMSeg(): lhw does not have a valid geometry\n");
    return (NULL);
  }

  if (wmseg == NULL) {
    wmseg = MRIallocSequence(seg->width, seg->height, seg->depth, MRI_INT, 1);
    MRIcopyHeader(seg, wmseg);
    MRIcopyPulseParameters(seg, wmseg);
  }
  if (MRIdimMismatch(seg, wmseg, 0)) {
    printf("ERROR: MRIannot2CerebralWMSeg(): dimension mismatch\n");
    return (NULL);
  }

  // Create an MRI for the anatomical voume the surfaces were generated from
  anat = MRIallocFromVolGeom(&(lhw->vg), MRI_INT, 1, 1);

  // Compute an LTA that maps from the anatomical to the segmentation
  if (anat2seg == NULL)
    lta = TransformRegDat2LTA(anat, seg, NULL);
  else {
    if (LTAmriIsTarget(anat2seg, seg))
      lta = LTAcopy(anat2seg, NULL);
    else
      lta = LTAinvert(anat2seg, NULL);
  }
  if (lta->type != LINEAR_VOX_TO_VOX) LTAchangeType(lta, LINEAR_VOX_TO_VOX);
  LTAfillInverse(lta);

  // Anatomical Vox to Surface RAS
  AnatVox2SurfRAS = MRIxfmCRS2XYZtkreg(anat);
  // Segmentation Vox to Surface RAS
  SegVox2SurfRAS = MatrixMultiplyD(AnatVox2SurfRAS, lta->inv_xforms[0].m_L, NULL);

  // Create the hash for faster service
  lhw_hash = MHTcreateVertexTable_Resolution(lhw, CURRENT_VERTICES, hashres);
  rhw_hash = MHTcreateVertexTable_Resolution(rhw, CURRENT_VERTICES, hashres);

  printf("  MRIannot2CerebralWMSeg(): looping over volume\n");
  fflush(stdout);
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  printf("     nthreads = %d\n", omp_get_max_threads());
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (c = 0; c < seg->width; c++) {
    ROMP_PFLB_begin
    
    int r, s, asegv, annot, annotid, wvtxno, segv, wmunknown;
    struct { float x,y,z; } vtx;
    MATRIX *RAS = NULL, *CRS = NULL;
    float wdw;
    MRIS *wsurf;
    MHT *whash = NULL;
    CRS = MatrixAlloc(4, 1, MATRIX_REAL);
    CRS->rptr[4][1] = 1;
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        asegv = MRIgetVoxVal(seg, c, r, s, 0);

        if (asegv != Left_Cerebral_White_Matter && asegv != Right_Cerebral_White_Matter) {
          // If not in WM, just copy seg to output and continue
          MRIsetVoxVal(wmseg, c, r, s, 0, asegv);
          continue;
        }

        if (asegv == Left_Cerebral_White_Matter) {
          wsurf = lhw;
          whash = lhw_hash;
          wmunknown = 5001;
        }
        else {
          wsurf = rhw;
          whash = rhw_hash;
          wmunknown = 5002;
        }

        // Compute location of voxel in surface space
        CRS->rptr[1][1] = c;
        CRS->rptr[2][1] = r;
        CRS->rptr[3][1] = s;
        RAS = MatrixMultiply(SegVox2SurfRAS, CRS, RAS);
        vtx.x = RAS->rptr[1][1];
        vtx.y = RAS->rptr[2][1];
        vtx.z = RAS->rptr[3][1];

        // Find closest white surface vertex and compute distance
        wvtxno = MHTfindClosestVertexNoXYZ(whash, wsurf, vtx.x, vtx.y, vtx.z, &wdw);
        if (wvtxno < 0) wvtxno = MRISfindClosestVertex(wsurf, vtx.x, vtx.y, vtx.z, &wdw, CURRENT_VERTICES);

        if (wdw <= DistThresh || DistThresh < 0) {
          // From the surface annotation, get the annotation number
          annot = wsurf->vertices[wvtxno].annotation;
          // Convert annotation number to an entry number
          CTABfindAnnotation(wsurf->ct, annot, &annotid);
          // Set segmentation to entry number + idbase
          if (annotid != -1)
            segv = annotid + wsurf->ct->idbase;
          else
            segv = wmunknown;  // no annotation present
        }
        else
          segv = wmunknown;
        MRIsetVoxVal(wmseg, c, r, s, 0, segv);
      }
    }
    MatrixFree(&CRS);
    if (RAS) MatrixFree(&RAS);
    ROMP_PFLB_end
  }
  ROMP_PF_end

  MHTfree(&lhw_hash);
  MHTfree(&rhw_hash);
  MatrixFree(&SegVox2SurfRAS);
  MatrixFree(&AnatVox2SurfRAS);
  MRIfree(&anat);
  LTAfree(&lta);

  return (wmseg);
}
/*
  \fn MRI *MRIunsegmentWM(MRI *seg, MRIS *lhw, MRIS *rhw, int *segidlist, int nlist, LTA *anat2seg, MRI *wmseg)
  \brief Changes a voxel segmentation to Left_Cerebral_White_Matter or Right_Cerebral_White_Matter
  depending on which surface it is closest to (no distance restriction). A voxel is relabeled if its
  segid in seg is in the segidlist. Can be done in place. anat2seg is an LTA that maps from the seg
  space the surface anatomical space. If NULL, then it assumes that the surface VOL_GEOM and the
  seg share a scanner RAS space. This function can be used to relabel hypointensities and CC.
 */
MRI *MRIunsegmentWM(MRI *seg, MRIS *lhw, MRIS *rhw, int *segidlist, int nlist, LTA *anat2seg, MRI *wmseg)
{
  MATRIX *AnatVox2SurfRAS, *SegVox2SurfRAS;
  float hashres = 16;
  MHT *lhw_hash = NULL, *rhw_hash = NULL;
  LTA *lta;
  MRI *anat;
  int c;

  if (lhw->vg.valid != 1) {
    printf("ERROR: MRIunsegmentWM(): lhw does not have a valid geometry\n");
    return (NULL);
  }

  if (wmseg == NULL) {
    wmseg = MRIallocSequence(seg->width, seg->height, seg->depth, MRI_INT, 1);
    MRIcopyHeader(seg, wmseg);
    MRIcopyPulseParameters(seg, wmseg);
  }
  if (MRIdimMismatch(seg, wmseg, 0)) {
    printf("ERROR: MRIunsegmentWM(): dimension mismatch\n");
    return (NULL);
  }

  // Create an MRI for the anatomical voume the surfaces were generated from
  anat = MRIallocFromVolGeom(&(lhw->vg), MRI_INT, 1, 1);

  // Compute an LTA that maps from the anatomical to the segmentation
  if (anat2seg == NULL)
    lta = TransformRegDat2LTA(anat, seg, NULL);
  else {
    if (LTAmriIsTarget(anat2seg, seg))
      lta = LTAcopy(anat2seg, NULL);
    else
      lta = LTAinvert(anat2seg, NULL);
  }
  if (lta->type != LINEAR_VOX_TO_VOX) LTAchangeType(lta, LINEAR_VOX_TO_VOX);
  LTAfillInverse(lta);

  // Anatomical Vox to Surface RAS
  AnatVox2SurfRAS = MRIxfmCRS2XYZtkreg(anat);
  // Segmentation Vox to Surface RAS
  SegVox2SurfRAS = MatrixMultiplyD(AnatVox2SurfRAS, lta->inv_xforms[0].m_L, NULL);

  // Create the hash for faster service
  lhw_hash = MHTcreateVertexTable_Resolution(lhw, CURRENT_VERTICES, hashres);
  rhw_hash = MHTcreateVertexTable_Resolution(rhw, CURRENT_VERTICES, hashres);

  printf("  MRIunsegmentWM(): looping over volume\n");
  fflush(stdout);
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  printf("     nthreads = %d\n", omp_get_max_threads());
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (c = 0; c < seg->width; c++) {
    ROMP_PFLB_begin
    
    int n, r, s, asegv, segv, hit, lhvtxno, rhvtxno;
    struct { float x,y,z; } vtx;
    MATRIX *RAS = NULL, *CRS = NULL;
    float lhd, rhd;
    CRS = MatrixAlloc(4, 1, MATRIX_REAL);
    CRS->rptr[4][1] = 1;
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        asegv = MRIgetVoxVal(seg, c, r, s, 0);

        hit = 0;
        for (n = 0; n < nlist; n++) {
          if (asegv == segidlist[n]) {
            hit = 1;
            break;
          }
        }
        if (hit == 0) {
          MRIsetVoxVal(wmseg, c, r, s, 0, asegv);
          continue;
        }

        CRS->rptr[1][1] = c;
        CRS->rptr[2][1] = r;
        CRS->rptr[3][1] = s;
        RAS = MatrixMultiply(SegVox2SurfRAS, CRS, RAS);
        vtx.x = RAS->rptr[1][1];
        vtx.y = RAS->rptr[2][1];
        vtx.z = RAS->rptr[3][1];

        lhvtxno = MHTfindClosestVertexNoXYZ(lhw_hash, lhw, vtx.x, vtx.y, vtx.z, &lhd);
        if (lhvtxno < 0) lhvtxno = MRISfindClosestVertex(lhw, vtx.x, vtx.y, vtx.z, &lhd, CURRENT_VERTICES);

        rhvtxno = MHTfindClosestVertexNoXYZ(rhw_hash, rhw, vtx.x, vtx.y, vtx.z, &rhd);
        if (rhvtxno < 0) rhvtxno = MRISfindClosestVertex(rhw, vtx.x, vtx.y, vtx.z, &rhd, CURRENT_VERTICES);

        if (lhd < rhd)
          segv = Left_Cerebral_White_Matter;
        else
          segv = Right_Cerebral_White_Matter;
        MRIsetVoxVal(wmseg, c, r, s, 0, segv);
      }
    }
    MatrixFree(&CRS);
    if (RAS) MatrixFree(&RAS);
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

  MHTfree(&lhw_hash);
  MHTfree(&rhw_hash);
  MatrixFree(&SegVox2SurfRAS);
  MatrixFree(&AnatVox2SurfRAS);
  MRIfree(&anat);
  LTAfree(&lta);

  return (wmseg);
}

/*
  \fn MRI *MRIrelabelHypoHemi(MRI *seg, MRIS *lhw, MRIS *rhw, LTA *anat2seg, MRI *wmseg)
  \brief Finds voxels labeled as WM_hypointensities to
  Left_WM_hypointensities or Right_WM_hypointensities depending on
  proximity to lh or rh white surface (lhw, rhw).  anat2seg is an LTA
  that maps from the seg space the surface anatomical space. If NULL,
  then it assumes that the surface VOL_GEOM and the seg share a
  scanner RAS space. See also MRIunsegmentWM().
 */
MRI *MRIrelabelHypoHemi(MRI *seg, MRIS *lhw, MRIS *rhw, LTA *anat2seg, MRI *wmseg)
{
  MATRIX *AnatVox2SurfRAS, *SegVox2SurfRAS;
  float hashres = 16;
  MHT *lhw_hash = NULL, *rhw_hash = NULL;
  LTA *lta;
  MRI *anat;
  int c;

  if (lhw->vg.valid != 1) {
    printf("ERROR: MRIrelabelHypoHemi(): lhw does not have a valid geometry\n");
    return (NULL);
  }

  if (wmseg == NULL) {
    wmseg = MRIallocSequence(seg->width, seg->height, seg->depth, MRI_INT, 1);
    MRIcopyHeader(seg, wmseg);
    MRIcopyPulseParameters(seg, wmseg);
  }
  if (MRIdimMismatch(seg, wmseg, 0)) {
    printf("ERROR: MRIrelabelHypoHemi(): dimension mismatch\n");
    return (NULL);
  }

  // Create an MRI for the anatomical volume the surfaces were generated from
  anat = MRIallocFromVolGeom(&(lhw->vg), MRI_INT, 1, 1);

  // Compute an LTA that maps from the anatomical to the segmentation
  if (anat2seg == NULL)
    lta = TransformRegDat2LTA(anat, seg, NULL);
  else {
    if (LTAmriIsTarget(anat2seg, seg))
      lta = LTAcopy(anat2seg, NULL);
    else
      lta = LTAinvert(anat2seg, NULL);
  }
  if (lta->type != LINEAR_VOX_TO_VOX) LTAchangeType(lta, LINEAR_VOX_TO_VOX);
  LTAfillInverse(lta);

  // Anatomical Vox to Surface RAS
  AnatVox2SurfRAS = MRIxfmCRS2XYZtkreg(anat);
  // Segmentation Vox to Surface RAS
  SegVox2SurfRAS = MatrixMultiplyD(AnatVox2SurfRAS, lta->inv_xforms[0].m_L, NULL);

  // Create the hash for faster service
  lhw_hash = MHTcreateVertexTable_Resolution(lhw, CURRENT_VERTICES, hashres);
  rhw_hash = MHTcreateVertexTable_Resolution(rhw, CURRENT_VERTICES, hashres);

  printf("  MRIrelabelHypoHemi(): looping over volume\n");
  fflush(stdout);
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  printf("     nthreads = %d\n", omp_get_max_threads());
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (c = 0; c < seg->width; c++) {
    ROMP_PFLB_begin
    
    int r, s, asegv, segv, lhvtxno, rhvtxno;
    struct { float x,y,z; } vtx;
    MATRIX *RAS = NULL, *CRS = NULL;
    float lhd, rhd;
    CRS = MatrixAlloc(4, 1, MATRIX_REAL);
    CRS->rptr[4][1] = 1;
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        asegv = MRIgetVoxVal(seg, c, r, s, 0);
        if (asegv != WM_hypointensities) {
          MRIsetVoxVal(wmseg, c, r, s, 0, asegv);
          continue;
        }

        CRS->rptr[1][1] = c;
        CRS->rptr[2][1] = r;
        CRS->rptr[3][1] = s;
        RAS = MatrixMultiply(SegVox2SurfRAS, CRS, RAS);
        vtx.x = RAS->rptr[1][1];
        vtx.y = RAS->rptr[2][1];
        vtx.z = RAS->rptr[3][1];

        lhvtxno = MHTfindClosestVertexNoXYZ(lhw_hash, lhw, vtx.x, vtx.y, vtx.z, &lhd);
        if (lhvtxno < 0) lhvtxno = MRISfindClosestVertex(lhw, vtx.x, vtx.y, vtx.z, &lhd, CURRENT_VERTICES);

        rhvtxno = MHTfindClosestVertexNoXYZ(rhw_hash, rhw, vtx.x, vtx.y, vtx.z, &rhd);
        if (rhvtxno < 0) rhvtxno = MRISfindClosestVertex(rhw, vtx.x, vtx.y, vtx.z, &rhd, CURRENT_VERTICES);

        if (lhd < rhd)
          segv = Left_WM_hypointensities;
        else
          segv = Right_WM_hypointensities;
        MRIsetVoxVal(wmseg, c, r, s, 0, segv);
      }
    }
    MatrixFree(&CRS);
    if (RAS) MatrixFree(&RAS);
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

  MHTfree(&lhw_hash);
  MHTfree(&rhw_hash);
  MatrixFree(&SegVox2SurfRAS);
  MatrixFree(&AnatVox2SurfRAS);
  MRIfree(&anat);
  LTAfree(&lta);

  return (wmseg);
}

/*
  \fn MRI *MRIunsegmentCortex(MRI *seg, int lhmin, int lhmax, int rhmin, int rhmax, MRI *out)
  \brief Replaces voxels in seg that have lhmin <= segid <= lhmax with
  Left_Cerebral_Cortex and rhmin <= segid <= rhmax with
  Right_Cerebral_Cortex. If lhmax or rhmax are less than 0, then no
  upper limit is used This is used prior to running MRIhiresSeg().
  It's a long story. Can be done in-place.  Ideally, seg is something
  like aparc+aseg.mgz */
MRI *MRIunsegmentCortex(MRI *seg, int lhmin, int lhmax, int rhmin, int rhmax, MRI *out)
{
  int c;

  if (out == NULL) {
    out = MRIallocSequence(seg->width, seg->height, seg->depth, MRI_INT, 1);
    if (out == NULL) return (NULL);
    MRIcopyHeader(seg, out);
    MRIcopyPulseParameters(seg, out);
  }
  if (MRIdimMismatch(seg, out, 0)) {
    printf("ERROR: MRIunsegmentCortex() dim mismatch between seg and out\n");
    return (NULL);
  }

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (c = 0; c < seg->width; c++) {
    ROMP_PFLB_begin
    
    int r, s, segid;
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        segid = MRIgetVoxVal(seg, c, r, s, 0);
        if (segid >= lhmin && (lhmax < 0 || segid <= lhmax))
          MRIsetVoxVal(out, c, r, s, 0, Left_Cerebral_Cortex);
        else if (segid >= rhmin && (rhmax < 0 || segid <= rhmax))
          MRIsetVoxVal(out, c, r, s, 0, Right_Cerebral_Cortex);
        else
          MRIsetVoxVal(out, c, r, s, 0, segid);
      }
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (out);
}

/*
  \fn MRI *MRIrelabelNonWMHypos(MRI *seg0, int *segidlist, int nsegs, int *outsegidlist)
  \brief Relabels non-wm hypointenties (80,81,82) based on
  proximity/most frequent nearest neighbor to labels listed in
  segidlist. The most frequent nearest neighbor to a hypo is
  determined.  If the most freq is in segidlist, then the hypo is
  relabled according to the corresponding segid in outsegidlist.  If
  none of the nearest neighbors are in the segidlist, the hypo is not
  relabeled.
 */
MRI *MRIrelabelNonWMHypos(MRI *seg0, int *segidlist, int nsegs, int *outsegidlist)
{
  int c, r, s, c0, r0, s0, k, n, loop, nchanged, mfsegid, nbrmax, nthnbr, nhits, segid, hit;
  int *clist, *rlist, *slist, *hitlist, nbrlist[3 * 3 * 3], nchangedtot;
  MRI *seg, *newseg;

  seg = MRIcopy(seg0, NULL);
  newseg = MRIcopy(seg0, NULL);

  // Get a count of non-wm-hypos
  nhits = 0;
  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        segid = MRIgetVoxVal(seg, c, r, s, 0);
        if (segid != 80 && segid != 81 && segid != 82) continue;
        nhits++;
      }
    }
  }
  printf("MRIrelabelNonWMHypos(): found %d non-WM-hypointensities\n", nhits);

  // Get a list of cols, rows, and slices of the non-wm-hypos
  clist = (int *)calloc(nhits, sizeof(int));
  rlist = (int *)calloc(nhits, sizeof(int));
  slist = (int *)calloc(nhits, sizeof(int));
  nhits = 0;
  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        segid = MRIgetVoxVal(seg, c, r, s, 0);
        if (segid != 80 && segid != 81 && segid != 82) continue;
        clist[nhits] = c;
        rlist[nhits] = r;
        slist[nhits] = s;
        nhits++;
      }
    }
  }

  // Loop dilating the segs adjacent to hypos until there are no changes
  hitlist = (int *)calloc(nhits, sizeof(int));
  nchangedtot = 0;
  loop = 0;
  nchanged = 1;
  while (nchanged != 0) {
    printf("loop %2d ", loop);
    fflush(stdout);
    loop++;
    nchanged = 0;
    // go through the hypo list
    for (k = 0; k < nhits; k++) {
      if (hitlist[k]) continue;
      c0 = clist[k];
      r0 = rlist[k];
      s0 = slist[k];
      // Get a list of neihbors in 3x3x3 neighborhood around the hypo
      nthnbr = 0;
      for (c = c0 - 1; c <= c0 + 1; c++) {
        for (r = r0 - 1; r <= r0 + 1; r++) {
          for (s = s0 - 1; s <= s0 + 1; s++) {
            if (c < 0 || c >= seg->width) continue;
            if (r < 0 || r >= seg->height) continue;
            if (s < 0 || s >= seg->depth) continue;
            if (c == c0 && r == r0 && s == s0) continue;
            segid = MRIgetVoxVal(seg, c, r, s, 0);
            /* Require that nbr be in segidlist. Can have a situation
            where the most freq nbr may not be in the segidlist but it
            gets labeled to the most freq neighbor that is in the
            segidlist*/
            hit = 0;
            for (n = 0; n < nsegs; n++) {
              if (segid == segidlist[n]) {
                hit = 1;
                break;
              }
            }
            if (hit == 0) continue;
            nbrlist[nthnbr] = segid;
            nthnbr++;
          }
        }
      }
      if (nthnbr == 0) continue;
      mfsegid = most_frequent_int_list(nbrlist, nthnbr, &nbrmax);
      fflush(stdout);
      for (n = 0; n < nsegs; n++) {
        if (mfsegid == segidlist[n]) {
          // relabel hypo in seg as most freq adjacent structure (dilation)
          MRIsetVoxVal(seg, c0, r0, s0, 0, segidlist[n]);
          // relabel hypo in the output based on the output segidlist
          MRIsetVoxVal(newseg, c0, r0, s0, 0, outsegidlist[n]);
          hitlist[k] = 1;
          nchanged++;
          nchangedtot++;
          break;
        }
      }  // nsegs
    }    // k
    printf("  nchanged %4d\n", nchanged);
    fflush(stdout);
  }  // while

  printf("MRIrelabelNonWMHypos(): relabeled %d non-WM-hypointensities\n", nchangedtot);
  fflush(stdout);

  free(clist);
  free(rlist);
  free(slist);
  free(hitlist);
  MRIfree(&seg);
  return (newseg);
}

/*!
\fn MRI *CTABcount2MRI(COLOR_TABLE *ct, MRI *seg)
\brief Creates an MRI structure with number of columns equal to the
number of non-null entries in the table. The value is set the to count
for that entry in the ctab.
*/
MRI *CTABcount2MRI(COLOR_TABLE *ct, MRI *seg)
{
  int n, ntot;
  MRI *mri;
  float voxsize;

  voxsize = seg->xsize * seg->ysize * seg->zsize;

  ntot = 0;
  for (n = 1; n < ct->nentries; n++)
    if (ct->entries[n]) ntot++;
  mri = MRIalloc(ntot, 1, 1, MRI_FLOAT);

  ntot = 0;
  for (n = 1; n < ct->nentries; n++)
    if (ct->entries[n]) {
      MRIsetVoxVal(mri, ntot, 0, 0, 0, voxsize * ct->entries[n]->count);
      ntot++;
    }

  return (mri);
}

/*!
  MRI *MRIreorientLIA2RAS(MRI *mriA, MRI *mriB)
  Reorient a volume that is LIA to be RAS
 */
MRI *MRIreorientLIA2RAS(MRI *mriA, MRI *mriB)
{
  MATRIX *vox2rasA, *MdcA, *MdcB, *crs0, *P0B, *DB, *vox2rasB;
  int r;
  char ostr[4];

  MRIdircosToOrientationString(mriA, ostr);
  if(strcmp(ostr,"LIA") != 0){
    printf("ERROR: MRIreorientLIA2RAS(): input ostring is %s, must be LIA\n",ostr);
    return(NULL);
  }
  if(mriB){
    MRIdircosToOrientationString(mriB, ostr);
    if(strcmp(ostr,"RAS") != 0){
      printf("ERROR: MRIreorientLIA2RAS(): output ostring is %s, must be RAS\n",ostr);
      return(NULL);
    }
    if(mriA->width != mriB->width){
      printf("ERROR: MRIreorientLIA2RAS(): input/output width mismatch %d %d\n",mriA->width,mriB->width);
      return(NULL);
    }
    if(mriA->height != mriB->depth){
      printf("ERROR: MRIreorientLIA2RAS(): input/output height mismatch %d %d\n",mriA->height,mriB->depth);
      return(NULL);
    }
    if(mriA->depth != mriB->height){
      printf("ERROR: MRIreorientLIA2RAS(): input/output depth mismatch %d %d\n",mriA->depth,mriB->height);
      return(NULL);
    }
  }

  vox2rasA = MRIxfmCRS2XYZ(mriA,0);

  // Create the new MdcB by swaping and negating apporpriately
  MdcA = MRImatrixOfDirectionCosines(mriA, NULL);
  MdcB = MatrixAlloc(4,4,MATRIX_REAL);
  for(r=1; r<=3; r++) MdcB->rptr[r][1] = -MdcA->rptr[r][1];
  for(r=1; r<=3; r++) MdcB->rptr[r][2] = +MdcA->rptr[r][3];
  for(r=1; r<=3; r++) MdcB->rptr[r][3] = -MdcA->rptr[r][2];

  // The origin (P0B) of the new volume will land at (Nx,0,Nz) of original
  crs0 = MatrixAlloc(4,1,MATRIX_REAL);
  crs0->rptr[1][1] = mriA->width-1;
  crs0->rptr[2][1] = mriA->height-1;
  crs0->rptr[3][1] = 0;
  crs0->rptr[4][1] = 1;
  P0B = MatrixMultiply(vox2rasA,crs0,NULL);

  // Matrix of voxel sizes for the new volume
  DB = MatrixIdentity(4,NULL);
  DB->rptr[1][1] = mriA->xsize;
  DB->rptr[2][2] = mriA->zsize;
  DB->rptr[3][3] = mriA->ysize;

  //vox2rasB = Mdc*DB, then fill in P0
  vox2rasB = MatrixMultiply(MdcB,DB,NULL);
  for(r=1; r<=3; r++) vox2rasB->rptr[r][4] = P0B->rptr[r][1];

  if(mriB == NULL) 
    mriB = MRIalloc(mriA->width,mriA->depth,mriA->height,mriA->type);

  // Set the MRI structure geometry from the vox2ras for B
  MRIsetVox2RASFromMatrix(mriB, vox2rasB);

  // Now do the resampling (might want to take this out so that caller can choose)
  MRIvol2Vol(mriA,mriB,NULL,SAMPLE_NEAREST,0);

  MatrixFree(&vox2rasA);
  MatrixFree(&MdcA);
  MatrixFree(&MdcB);
  MatrixFree(&crs0);
  MatrixFree(&P0B);
  MatrixFree(&DB);
  MatrixFree(&vox2rasB);

  return(mriB);
}

/*!
  \fn MATRIX *MRIvol2mat(MRI *vol, MRI *mask, int transposeFlag, MATRIX *M)
  \brief Converts a volume into a matrix. With transposeFlag=0, the
  output matrix will be nframes-by-nvoxels. If transposeFlag=1, then
  the transpose is returned. If mask is non-null, then nvoxels is the
  number of voxels > 0.5 in the mask. The order, from slowest to
  fastest, that the voxels are packed into the matrix is slice,
  column, row. This order is compatible with matlab. This function is
  compatible with MRImat2vol(). See also GTMvol2mat().
 */
MATRIX *MRIvol2mat(MRI *vol, MRI *mask, int transposeFlag, MATRIX *M)
{
  int nvox, c, r ,s, f, nrows, ncols, nthvox;

  if(mask && MRIdimMismatch(vol, mask, 0)) {
    printf("ERROR: MRIvol2mat(): mask and vol dimension mismatch\n");
    return (NULL);
  }

  if(mask) nvox = MRIcountAboveThreshold(mask,0.5);
  else     nvox = vol->width * vol->height * vol->depth;

  if(transposeFlag==0){
    nrows = vol->nframes;
    ncols = nvox;
  } else {
    nrows = nvox;
    ncols = vol->nframes;
  }

  if(M==NULL) M = MatrixAlloc(nrows,ncols,MATRIX_REAL);
  if(M->rows != nrows || M->cols != ncols){
    printf("ERROR: MRIvol2mat(): dimension mismatch expecting (%d,%d), got (%d,%d)\n",
	   nrows,ncols,M->rows,M->cols);
    return(NULL);
  }

  // col, row, slice order is compatible with MRImat2vol and matlab
  nthvox = 0;
  for(s=0; s < vol->depth; s++){
    for(c=0; c < vol->width; c++){
      for(r=0; r < vol->height; r++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	// Not the most efficient to parallelize here, but can't do above
        #ifdef HAVE_OPENMP
        #pragma omp parallel for 
        #endif
	for(f = 0; f < vol->nframes; f++){
	  int Mc, Mr;
	  if(transposeFlag==0){
	    Mr = f+1; 
	    Mc = nthvox+1;
	  }
	  else {
	    Mr = nthvox+1;
	    Mc = f+1; 
	  }
	  M->rptr[Mr][Mc] = MRIgetVoxVal(vol,c,r,s,f);
	} // loop over frame
	nthvox++;
      }
    }
  }

  return(M);
}

/*!
  \fn MRI *MRImat2vol(MATRIX *M, MRI *mask, int transposeFlag, MRI *vol)
  \brief Converts a matrix into a volume. With transposeFlag=0, the
  matrix should have size nframes-by-nvoxels (or nvoxel-by-nframes if
  transposeFlag=1).  If mask is non-null, then nvoxels is the number
  of voxels > 0.5 in the mask. The order, from slowest to fastest,
  that the voxels are unpacked from the matrix is slice, column,
  row. This order is compatible with matlab. This function is
  compatible with MRIvol2mat(). See also GTMmat2vol(). Either mask or
  vol must be non-NULL, otherwise, there is no way to know how big the
  volume should be. vol does not need to be zeroed.
 */
MRI *MRImat2vol(MATRIX *M, MRI *mask, int transposeFlag, MRI *vol)
{
  int nvox, c, r ,s, f, nthvox;
  int nframes;

  if(transposeFlag==0) nframes = M->rows;
  else                 nframes = M->cols;

  if(mask == NULL && vol == NULL){
    printf("ERROR: MRImat2vol(): both mask and vol are NULL\n");
    return(NULL);
  }
  if(vol == NULL){
    vol = MRIallocSequence(mask->width, mask->height, mask->depth, MRI_FLOAT, nframes);
    MRIcopyHeader(mask, vol);
    MRIcopyPulseParameters(mask, vol);
  }
  if(mask && MRIdimMismatch(vol, mask, 0)) {
    printf("ERROR: MRImat2vol(): mask and vol dimension mismatch\n");
    return (NULL);
  }
  if(vol->nframes != nframes){
    printf("ERROR: MRImat2vol(): vol and matrix frame dimension mismatch\n");
    return (NULL);
  }
  if(mask) nvox = MRIcountAboveThreshold(mask,0.5);
  else     nvox = vol->width * vol->height * vol->depth;
  if( (transposeFlag==0 && M->cols != nvox) || (transposeFlag==1 && M->rows != nvox)){
    printf("ERROR: MRImat2vol(): vol and matrix vox dimension mismatch\n");
    printf("   transposeFlag=%d, rows = %d, cols = %d, nvox = %d\n",transposeFlag,M->rows,M->cols,nvox);
    return (NULL);
  }

  // col, row, slice order is compatible with MRIvol2mat and matlab
  nthvox = 0;
  for(s=0; s < vol->depth; s++){
    for(c=0; c < vol->width; c++){
      for(r=0; r < vol->height; r++){
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5){
	  for(f = 0; f < vol->nframes; f++) MRIsetVoxVal(vol,c,r,s,f,0);
	  continue;
	}
	// Not the most efficient to parallelize here, but can't do above
        #ifdef HAVE_OPENMP
        #pragma omp parallel for 
        #endif
	for(f = 0; f < vol->nframes; f++){
	  int Mr, Mc;
	  if(transposeFlag==0){
	    Mr = f+1; 
	    Mc = nthvox+1;
	  }
	  else {
	    Mr = nthvox+1;
	    Mc = f+1; 
	  }
	  MRIsetVoxVal(vol,c,r,s,f,M->rptr[Mr][Mc]);
	}// loop over frame
	nthvox++;
      }
    }
  }

  return(vol);
}
/*!
  \fn MRI *MRImergeSegs(MRI *seg, int *seglist, int nsegs, int NewSegId, MRI *newseg)
  \brief Merges multiple segmentations into one. Can be done in-place.
  \parameter seg - original segmentation
  \parameter seglist - list of segmentation IDs to merge
  \parameter nsegs - length of list
  \parameter NewSegId - replace values in list with NewSegId
  \parameter newseg - new segmentation (also passed as output)
*/
MRI *MRImergeSegs(MRI *seg, int *seglist, int nsegs, int NewSegId, MRI *newseg)
{
  int c,r,s,n,segid;

  if(newseg == NULL) newseg = MRIcopy(seg,NULL);

  for (c=0; c < seg->width; c++){
    for (r=0; r < seg->height; r++){
      for (s=0; s < seg->depth; s++){
	segid = MRIgetVoxVal(seg,c,r,s,0);
	MRIsetVoxVal(newseg,c,r,s,0, segid);
	for(n=0; n < nsegs; n++){
	  if(segid == seglist[n]){
	    MRIsetVoxVal(newseg,c,r,s,0,NewSegId);
	    break;
	  }
	}
      }
    }
  }
  return(newseg);
}
/*
  \fn MRI *MRImatchSegs(MRI *seg, int *seglist, int nsegs, int MaskId, MRI *mask)
  \brief Creates a binary mask of voxels that match any of the IDs in the
    segmentations list. Can be done in-place.
  \parameter seg - original segmentation
  \parameter seglist - list of segmentation IDs to merge
  \parameter nsegs - length of list
  \parameter MaskId - replace values in list with MaskId
  \parameter mask - new segmentation (also passed as output)
*/
MRI *MRImatchSegs(MRI *seg, int *seglist, int nsegs, int MaskId, MRI *mask)
{
  int c,r,s,n,segid;

  if(mask == NULL) mask = MRIcopy(seg,NULL);

  for (c=0; c < seg->width; c++){
    for (r=0; r < seg->height; r++){
      for (s=0; s < seg->depth; s++){
	MRIsetVoxVal(mask,c,r,s,0, 0);
	segid = MRIgetVoxVal(seg,c,r,s,0);
	for(n=0; n < nsegs; n++){
	  if(segid == seglist[n]){
	    MRIsetVoxVal(mask,c,r,s,0,MaskId);
	    break;
	  }
	}
      }
    }
  }
  return(mask);
}


/*!
  \fn HISTOGRAM *HISTOseg(MRI *seg, int segid, MRI *vol, double bmin, double bmax, double bdelta)
  \brief Creates a histogram from the intensities in vol from the voxels in the given
  segmentation. The caller supplies the min, max, and delta for the bins of the histogram.
  Can't include this in histo.c because of circular dependence.
 */
HISTOGRAM *HISTOseg(MRI *seg, int segid, MRI *vol, double bmin, double bmax, double bdelta)
{
  HISTOGRAM *h;
  double v,vsegid;
  int c,r,s,nbins,binno;

  nbins = round((bmax-bmin)/bdelta) + 1;
  h = HISTOinit(NULL, nbins, bmin, bmax);
  for(c=0; c < seg->width; c++){
    for(r=0; r < seg->height; r++){
      for(s=0; s < seg->depth; s++){
	vsegid  = MRIgetVoxVal(seg,c,r,s,0);
	if(vsegid!=segid) continue;
	v = MRIgetVoxVal(vol,c,r,s,0);
	binno = round(v-bmin)/bdelta;
	if(binno < 0)      binno=0;
	if(binno >= nbins) binno=nbins-1;
	h->counts[binno] ++;
      }
    }
  }
  return(h);
}

/*!
\fn MRI *MRIvolTopoFix(MRI *binseg0, int onval, MRI *binseg, int ndilationsmax, int nerodesmax, int ntopoerodesmax)
\brief Fixes the topology of the binary segmentation using dilations
and erosions. The segmentation must be binary with on value =
onval. Fills holes only. The new segmentation will have all the same
voxels as the input seg plus whatever need to be added to make the seg
topo correct. The method is deterministic when given all the same
inputs, but there may be differences with simple shifts or crops. Some
weird things can happen if the dilation hits the edge of the volume;
need to debug at some point, but make sure you have good padding.
*/
MRI *MRIvolTopoFix(MRI *binseg0, int onval, MRI *binseg, int ndilationsmax, int nerodesmax, int ntopoerodesmax)
{
  if(binseg0 == binseg){
    printf("ERROR: FillIt(): output=input\n");
    return(NULL);
  }
  // mri has to be binarized with on-value=onval

  // First check whether it is already topo correct
  int nvertices, nfaces, nedges,eno,n;
  MRIS *surf;
  surf = MRIStessellate(binseg0,onval,0);
  eno = MRIScomputeEulerNumber(surf, &nvertices, &nfaces, &nedges) ;
  MRISfree(&surf);
  if(eno == 2) return(binseg0);

  binseg = MRIcopy(binseg0,binseg);
  MRIcopyPulseParameters(binseg0,binseg);
  if(binseg0->ct) binseg->ct = CTABdeepCopy(binseg0->ct);

  // Get a list of the origninal voxels to protect them while eroding
  // I think this could be multi-threaded
  int hitedge=0;
  std::vector<std::vector<int>> crs;
  for(int c=0; c < binseg->width; c++){
    for(int r=0; r < binseg->height; r++){
      for(int s=0; s < binseg->depth; s++){
	double val = MRIgetVoxVal(binseg,c,r,s,0);
	if(val < 0.5) continue;
	std::vector<int> crs0;
	crs0.push_back(c);
	crs0.push_back(r);
	crs0.push_back(s);
	crs.push_back(crs0);
	if(c == 0 || c == binseg->width-1 ||
	   r == 0 || r == binseg->height-1 ||
	   s == 0 || s == binseg->depth-1) hitedge ++;
      }
    }
  }
  printf("Found %d voxels in binseg\n",(int)crs.size());
  if(hitedge) printf("WARNING: binseg has %d edge voxels. You may want "
		     "to add some padding if results are not good\n",hitedge);

  // Make a copy so can protect during topoerode
  MRI *mrisrc = MRIcopy(binseg, NULL);
  MRIcopyPulseParameters(binseg,mrisrc);

  // Dilate until the topo defects go away. Weird stuff can happen in
  // the erosion if dilation hits the edge.  Could just fill the
  // entire volume or use MRIvolTopoDilateOne() or wrap an ellipsoid
  // around the shape.
  printf("\nDilating ============================\n");
  n=0;
  while(1){
    n++;
    printf("dilate n=%d eno=%d\n",n,eno);
    MRI *mritmp = MRIdilate6(binseg, NULL);
    MRIsetEdges(mritmp, 0, mritmp); // edge hack, does not work well
    MRIfree(&binseg);
    binseg = mritmp;
    surf = MRIStessellate(binseg,onval,0);
    eno = MRIScomputeEulerNumber(surf, &nvertices, &nfaces, &nedges) ;
    MRISfree(&surf);
    if(eno == 2 || n > ndilationsmax){
      if(eno == 2)  printf("  no defects found, breaking from dilation\n");
      else{
	printf("  ERROR: number of dilations > ndilationsmax=%d\n",ndilationsmax);
	return(NULL);
      }
      break;
    }
  }

  // Erode until a topo defect appears (take the previous one)
  // This is probably redundant with TopoErode below, but it 
  // may speed things up.
  printf("\nEroding ============================\n");
  n=0;
  while(1){
    n++;
    printf("erode n=%d eno=%d\n",n,eno);
    MRI *mricopy = MRIcopy(binseg, NULL); // copy in case defect shows up here
    // Erode
    MRI *mritmp = MRIerodeNN(binseg,NULL,NEAREST_NEIGHBOR_FACE,1);
    // Protect the true seg voxels from erosion
    for(int k=1; k < crs.size(); k++)
      MRIsetVoxVal(mritmp,crs[k][0],crs[k][1],crs[k][2],0,onval);
    // Test whether there is a defect after the erosion
    surf = MRIStessellate(mritmp,onval,0);
    eno = MRIScomputeEulerNumber(surf, &nvertices, &nfaces, &nedges) ;
    MRISfree(&surf);
    if(eno != 2 || n > nerodesmax){
      if(eno != 2)  printf("  defects found, breaking from erode, and reverting to previous\n");
      else          printf("  WARNING: number of erodes > nerodesmax=%d\n",nerodesmax);
      MRIfree(&mritmp);
      MRIfree(&binseg);
      binseg = mricopy;
      break;
    }
    MRIfree(&mricopy);
    MRIfree(&binseg);
    binseg = mritmp;
  }

  printf("\nTopoEroding ============================\n");
  n=0;
  while(1){
    n++;
    int nhits = MRIvolTopoErodeOne(binseg, mrisrc);
    surf = MRIStessellate(binseg,onval,0);
    eno = MRIScomputeEulerNumber(surf, &nvertices, &nfaces, &nedges) ;
    MRISfree(&surf);
    printf("n=%d nhits = %d eno=%d\n",n,nhits,eno);
    if(nhits == 0 || n > ntopoerodesmax){
      if(n > ntopoerodesmax) printf("  WARNING: number of topoerodes > max=%d\n",ntopoerodesmax);
      break;
    }
  }

  MRIcopyPulseParameters(mrisrc,binseg);
  MRIfree(&mrisrc);

  return(binseg);
}

/*!
\fn int MRIvolTopoErodeOne(MRI *binvol, const MRI *keepmask)
\brief Erodes a layer of boundary voxels. A boundary voxel is only eroded if its erosion 
would not change the Euler characteristic. Binvol is a binary volume where on is >0.5.
Voxels in the keepmask are not eroded. The final result may depend on the order that
the voxels are searched (eg, rotating the input could change results).
*/
int MRIvolTopoErodeOne(MRI *binvol, const MRI *keepmask)
{
  // only consider voxels on the boundary of the input binvol
  MRI *mrie6 = MRIerodeNN(binvol,NULL,NEAREST_NEIGHBOR_FACE,1);

  int nhits=0;
  for(int c=1; c < binvol->width-1; c++){
    for(int r=1; r < binvol->height-1; r++){
      for(int s=1; s < binvol->depth-1; s++){
	if(keepmask && MRIgetVoxVal(keepmask,c,r,s,0)>0.5) continue;
	int v = MRIgetVoxVal(binvol,c,r,s,0);
	if(v < 0.5) continue; // voxel not set so can't erode
	v = MRIgetVoxVal(mrie6,c,r,s,0);
	if(v > 0.5) continue; // voxel not a boundary voxel
	int dEC = QuadEulerCharChange(binvol, NULL, c, r, s);
	if(dEC != 0) continue; // skip voxel because eroding it would change the EC
	MRIsetVoxVal(binvol,c,r,s,0,0);
	nhits++;
      }
    }
  }
  MRIfree(&mrie6);
  return(nhits);
}
/*!
\fn int MRIvolTopoDilateOne(MRI *binvol, int onval)
\brief Dilate a layer of boundary voxels. A boundary voxel is only
dilated if its dilation would not change the Euler
characteristic. Binvol is a binary volume where on is >0.5.  The final
result may depend on the order that the voxels are searched (eg,
rotating the input could change results).
*/
int MRIvolTopoDilateOne(MRI *binvol, int onval)
{
  // only consider voxels on the boundary of the input binvol
  MRI *mrid6 = MRIdilate6(binvol,NULL);

  int nhits=0,dECsum=0;
  for(int c=0; c < binvol->width; c++){
    for(int r=0; r < binvol->height; r++){
      for(int s=0; s < binvol->depth; s++){
	int v = MRIgetVoxVal(binvol,c,r,s,0);
	if(v > 0.5) continue; // voxel set so can't dilate
	v = MRIgetVoxVal(mrid6,c,r,s,0);
	if(v < 0.5) continue; // voxel not a boundary voxel
	int dEC = QuadEulerCharChange(binvol, NULL, c, r, s);
	dECsum += dEC;
	if(dEC == 0) continue; // skip voxel because dilating it would not change the EC
	MRIsetVoxVal(binvol,c,r,s,0,onval); //dilate
	nhits++;
      }
    }
  }
  MRIfree(&mrid6);
  return(nhits);
}
/*!
  \fn int QuadEulerCharChange(MRI *vol, MRI *mask, int c, int r, int s)
  \brief Determines how the Euler Characteristic of a hypothetical
  square (quad) mesh would change if the given voxel were to be turned
  on (assuming the mesh was created with it off). The idea is that the
  hypothetical quad mesh was tessellated on the given binary volume
  (vol). If the mesh has no topological defects, then it should have
  EC=2.  If the given voxel were to be added and the mesh regenerated,
  the mesh's EC would change by the amount returned by this function
  (computed without having to have ever generate a mesh). When reading
  the code below, that three things can happen when adding an element:
  (1) it does not overlap with an existing element, so it will change
  the EC. (2) It does overlap with an existing element and the element
  will disappear from the mesh when the new one is added (and so
  changes the EC), and (3) It does overlap with an existing element
  and the element will not disappear (and so does NOT change the EC).
  No bounds checking is done on (c,r,s); since the nearest neighbors
  are evaluated, (c,r,s) should be within [1:dimsize-2]. If mask is
  used, it should be the same size as vol. No neighbors are considered
  if they are ouside of the mask (this has not been tested). A voxel
  in vol or mask is considered "set" if its value is greather than
  0.5.  There is no check to determine whether the passed voxel is in
  the mask.  See also QuadEulerCharChangeTest().
*/
int QuadEulerCharChange(MRI *vol, MRI *mask, int c, int r, int s)
{
  // dont just return if this voxel is set this as it prevents the
  // ability to see what would happen if this voxel were to be unset.
  //if(MRIgetVoxVal(vol,c,r,s,0)>0.5) return(0); 
  int dc, dr, ds, dsum, nhits,debug=0;
  if(c==Gx && r==Gy && s==Gz) debug = 1;
  if(debug) printf("  debug vox %d %d %d\n",c,r,s);
  int deltaEC=0;
  for(dc = -1; dc <= 1; dc++){
    for(dr = -1; dr <= 1; dr++){
      for(ds = -1; ds <= 1; ds++){
	if(mask && MRIgetVoxVal(mask,c+dc,r+dr,s+ds,0) < 0.5) continue;
	dsum = fabs(dc) + fabs(dr) + fabs(ds);
	if(dsum==0) continue; // no need to do itself
	nhits = 0;
	// determine whether neighbor voxel is out of bounds (oob)
	int coob = 0;
	if(c+dc < 0 || c+dc >= vol->width)  coob=1;
	int roob = 0;
	if(r+dr < 0 || r+dr >= vol->height) roob=1;
	int soob = 0;
	if(s+ds < 0 || s+ds >= vol->depth)  soob=1;
	if(dsum==1){ //face neighbor
	  // look at single voxel that shares this face
	  if(!coob && !roob && !soob && MRIgetVoxVal(vol,c+dc,r+dr,s+ds,0) > 0.5){
	    // face is already part of the surface so will lose both.
	    // -1 means dont add new face and remove the face that is there
	    deltaEC--;
	    if(debug) printf("face %2d %2d %2d removed dEC=%d\n",dc,dr,ds,deltaEC);
	    continue;
	  }
	  // If it gets here, then the face can be added (either
	  // because the neighbor voxel is not set or it is OOB),
	  // which increases the EC by 1
	  deltaEC++;
	  if(debug) printf("face %2d %2d %2d added dEC=%d\n",dc,dr,ds,deltaEC);
	}
	if(dsum==2){ //edge
	  // Look at the three other voxels that share this edge
	  // One of the voxels is always at +(dc,dr,ds)
	  if(!coob && !roob && !soob &&  MRIgetVoxVal(vol,c+dc,r+dr,s+ds,0) > 0.5 ) nhits++;
	  // For the remaining two voxels ...
	  if(dc==0){
	    // One voxel is a +dr the other is at +ds
	    if(!roob && MRIgetVoxVal(vol, c, r+dr, s,    0) > 0.5 ) nhits++;
	    if(!soob && MRIgetVoxVal(vol, c, r   , s+ds, 0) > 0.5 ) nhits++;
	  }
	  if(dr==0){
	    // One voxel is a +dc the other is at +ds
	    if(!coob && MRIgetVoxVal(vol, c+dc, r, s,    0) > 0.5 ) nhits++;
	    if(!soob && MRIgetVoxVal(vol, c,    r, s+ds, 0) > 0.5 ) nhits++;
	  }
	  if(ds==0){
	    // One voxel is a +dc the other is at +dr
	    if(!coob &&  MRIgetVoxVal(vol, c+dc, r,    s, 0)  > 0.5 ) nhits++;
	    if(!roob &&  MRIgetVoxVal(vol, c,    r+dr, s, 0)  > 0.5 ) nhits++;
	  }
	  if(nhits == 0) {
	    // No other voxels claims this edge, so, if this voxel is added
	    // this edge causes the EC to decrease by 1
	    deltaEC--;
	    if(debug) printf("edge %2d %2d %2d added dEC=%d\n",dc,dr,ds,deltaEC);
	  }
	  if(nhits == 3) {	  
	    // All 3 other voxels claim this edge. If the voxel is added,
	    // then this edge will be lost, so EC increases by 1
	    deltaEC++;
	    if(debug) printf("edge %2d %2d %2d removed dEC=%d\n",dc,dr,ds,deltaEC);
	  }
	  // If nhits != 0 and nhits != 3 then the edge will correspond to an existing
	  // edge that will not be lost by adding this voxel. As a result, no change in EC
	}
	if(dsum==3){ //corner/vertex
	  // Look at the seven other voxels that share this corner
	  if(!soob &&                   MRIgetVoxVal(vol, c,    r,    s+ds, 0) > 0.5 ) nhits++;
	  if(!roob &&                   MRIgetVoxVal(vol, c,    r+dr, s,    0) > 0.5 ) nhits++;
	  if(!roob && !soob &&          MRIgetVoxVal(vol, c,    r+dr, s+ds, 0) > 0.5 ) nhits++;
	  if(!coob &&                   MRIgetVoxVal(vol, c+dc, r,    s,    0) > 0.5 ) nhits++;
	  if(!coob && !soob &&          MRIgetVoxVal(vol, c+dc, r,    s+ds, 0) > 0.5 ) nhits++;
	  if(!coob && !roob &&          MRIgetVoxVal(vol, c+dc, r+dr, s,    0) > 0.5 ) nhits++;
	  if(!coob && !roob && !soob && MRIgetVoxVal(vol, c+dc, r+dr, s+ds, 0) > 0.5 ) nhits++;
	  if(nhits == 0) {
	    // No other voxels claims this corner, so, if this voxel is added
	    // this corner causes the EC to increase by 1
	    deltaEC++;
	    if(debug) printf("vertex %2d %2d %2d added dEC=%d\n",dc,dr,ds,deltaEC);
	  }
	  if(nhits == 7) {	  
	    // All 7 other voxels claim this corner. If the voxel is added,
	    // then this corner will be lost, so EC decreases by 1
	    deltaEC--;
	    if(debug) printf("vertex %2d %2d %2d removed dEC=%d\n",dc,dr,ds,deltaEC);
	  }
	  // If nhits != 0 and nhits != 7 then the vertex will correspond to an existing
	  // vertex that will not be lost by adding this voxel. As a result, no change in EC
	}
	//printf("%2d %2d %2d   %d   %3d\n",dc,dr,ds,dsum,deltaEC);
      }
    }
  }
  if(debug) printf("  debug vox %d %d %d %d\n",c,r,s,deltaEC);
  return(deltaEC);
}

/*!
  \fn int QuadEulerCharChangeTest2(void)
  \brief Tests int QuadEulerCharChange() by creating some simple structures
  and seeing how the EC changes when voxels are added nearby. This is not
  an exhaustive test, but it covers a lot of territory. I wrote this not
  knowing that I already had a QuadEulerCharChangeTest(). When TurnOnTestVoxel
  is 1, then it will set the voxel being tested to 1 (in some cases). This should
  not have any effect on the output.
 */
int QuadEulerCharChangeTest2(int TurnOnTestVoxel)
{
  // First, create an MRI structure and create a "sphere" where all
  // voxels within a 3x3x3 are set *except* for the center voxel.
  // This structure will be topologically defective.
  MRI *mri;
  mri = MRIallocSequence(10,10,10,MRI_INT,1);
  for(int dc = -1; dc <= 1; dc++){
    for(int dr = -1; dr <= 1; dr++){
      for(int ds = -1; ds <= 1; ds++){
	int dsum = fabs(dc) + fabs(dr) + fabs(ds);
	if(dsum == 0) continue;
	MRIsetVoxVal(mri,5+dc,5+dr,5+ds,0,1);
      }
    }
  }

  int err = 0, c;

  // Test adding a cube that shares a single face  (nothing to do with sphere)
  // dEC = 4 - 8 + (5-1) = 0
  c = QuadEulerCharChange(mri, NULL, 4, 4, 3);
  if(TurnOnTestVoxel) MRIsetVoxVal(mri,4,4,3,0,1);
  if(c != 0) {
    printf("ERROR: add face dEC = %d exp 0\n",c);
    err = 2;
  }
  if(TurnOnTestVoxel) MRIsetVoxVal(mri,4,4,3,0,0);

  // Add edge neighbor (nothing to do with sphere)
  // dEC = 6 - 11 + 6 = +1
  if(TurnOnTestVoxel) MRIsetVoxVal(mri,4,3,3,0,1);
  c = QuadEulerCharChange(mri, NULL, 4, 3, 3);
  if(c != 1) {
    printf("ERROR: add edge dEC = %d exp 1\n",c);
    err = 3;
  }
  if(TurnOnTestVoxel) MRIsetVoxVal(mri,4,3,3,0,0);

  // Add vertex neighbor (nothing to do with sphere)
  // dEC = 7 - 12 + 6 = +1
  if(TurnOnTestVoxel) MRIsetVoxVal(mri,3,3,3,0,1);
  c = QuadEulerCharChange(mri, NULL, 3, 3, 3);
  if(c != 1) {
    printf("ERROR: add vertex dEC = %d exp 1\n",c);
    err = 4;
  }
  if(TurnOnTestVoxel) MRIsetVoxVal(mri,3,3,3,0,0);

  // Test how filling in the center voxel changes the EC
  // Adds no v, e, or f and all are  removed so dEC = -8 -(-12) + 6 = -2
  c = QuadEulerCharChange(mri, NULL, 5, 5, 5);
  if(c != -2) {
    printf("ERROR: fill sphere dEC = %d exp -2\n",c);
    err = 1;
  }

  // Unset a corner cube in the 3x3x3 and see what happens if set it
  // dEC = 1 - (3-3) + (3-3) = 1
  MRIsetVoxVal(mri,4,4,4,0,0);
  c = QuadEulerCharChange(mri, NULL, 4, 4, 4);
  if(c != 1){
    printf("remove/add corner cube dEC = %d exp 1\n",c);
    err = 5;
  }
  MRIsetVoxVal(mri,4,4,4,0,1);

  // Unset a cube in the middle of the side of the 3x3x3 and see what happens if set it
  // dEC = 0 -4 + (2-4) = +2
  MRIsetVoxVal(mri,5,5,4,0,0);
  c =QuadEulerCharChange(mri, NULL, 5, 5, 4);
  if(c != 2) {
    printf("ERROR: remove/add mid cube exp dEC = %d exp 2\n",c);
    err = 6;
  }
  MRIsetVoxVal(mri,5,5,4,0,1);

  // Unset a cube in the edge of the side of the 3x3x3 and see what happens if set it
  // dEC = 0 - (1-4) + (2-4) = 1
  MRIsetVoxVal(mri,4,5,4,0,0);
  c = QuadEulerCharChange(mri, NULL, 4, 5, 4);
  if(c != 1){
    printf("ERROR: remove/add edge cube exp dEC = %d exp 1\n",c);
    err = 7;
  }
  MRIsetVoxVal(mri,4,5,4,0,1);

  MRIfree(&mri);
  return(err);
}

/*!
  \fn int QuadEulerCharChangeTest(int ForceFail)
  \brief Test for QuadEulerCharChange(). The test is done by creating
  various configurations of on and off voxels inside of a 3x3x3
  volume. In each configuration, the EC with and without the center
  voxel set is known, so the change is known.  Each configuration is
  tested by permuting the dimensions in various ways since the change
  in EC should be invariant under these manipulations. There are a
  possible 2^26 possible configurations. Not all are tested:). If
  ForceFail is set to non-zero, then all tests will fail (but error
  messages are still printed).
 */

int QuadEulerCharChangeTest(int ForceFail)
{
  MRI *mri;
  int dc, dr, ds, dsum;
  int err=0;
  char testname[1000];

  // Set up a simple 3x3x3 volume
  mri = MRIalloc(3,3,3,MRI_INT);

  // Set each corner voxel one-by-one, expect EC to increase by 1  in each case
  for(dc = -1; dc <= 1; dc++){
    for(dr = -1; dr <= 1; dr++){
      for(ds = -1; ds <= 1; ds++){
	dsum = fabs(dc) + fabs(dr) + fabs(ds);
	if(dsum != 3) continue;
	MRIconst(3,3,3,1,0,mri); // set MRI to 0
	MRIsetVoxVal(mri,1+dc,1+dr,1+ds,0,1);
	sprintf(testname,"single-corner %d %d %d\n",dc,dr,ds);
	err = QuadEulerCharChangeCheckReorder(mri, testname, 1+ForceFail);
      }
    }
  }

  // Set all corner voxels.  In this case, 0 corners will be added, 6
  // faces added, 12 edges added so one expects the EC to change by
  // 0+6-12=-6
  MRIconst(3,3,3,1,0,mri); // set MRI to 0
  for(dc = -1; dc <= 1; dc++){
    for(dr = -1; dr <= 1; dr++){
      for(ds = -1; ds <= 1; ds++){
	dsum = fabs(dc) + fabs(dr) + fabs(ds);
	if(dsum != 3) continue;
	MRIsetVoxVal(mri,1+dc,1+dr,1+ds,0,1);
      }
    }
  }
  err = QuadEulerCharChangeCheckReorder(mri, "all-corners", -6+ForceFail);

  // Set each edge voxel one-by-one, expect EC to increase by 1 in each case
  for(dc = -1; dc <= 1; dc++){
    for(dr = -1; dr <= 1; dr++){
      for(ds = -1; ds <= 1; ds++){
	dsum = fabs(dc) + fabs(dr) + fabs(ds);
	if(dsum != 2) continue;
	MRIconst(3,3,3,1,0,mri); // set MRI to 0
	MRIsetVoxVal(mri,1+dc,1+dr,1+ds,0,1);
	sprintf(testname,"single-edge %d %d %d\n",dc,dr,ds);
	err = QuadEulerCharChangeCheckReorder(mri, testname, 1+ForceFail);
      }
    }
  }

  // Set all edge voxels. No vertices and no edges are added, but 6 faces are added, so
  // expect EC to change by 0+6-0=6. 
  MRIconst(3,3,3,1,0,mri); // set MRI to 0
  for(dc = -1; dc <= 1; dc++){
    for(dr = -1; dr <= 1; dr++){
      for(ds = -1; ds <= 1; ds++){
	dsum = fabs(dc) + fabs(dr) + fabs(ds);
	if(dsum != 2) continue;
	MRIsetVoxVal(mri,1+dc,1+dr,1+ds,0,1);
      }
    }
  }
  err = QuadEulerCharChangeCheckReorder(mri, "all-edges", +6+ForceFail);

  // Set each face voxel one-by-one, expect EC to not change
  for(dc = -1; dc <= 1; dc++){
    for(dr = -1; dr <= 1; dr++){
      for(ds = -1; ds <= 1; ds++){
	dsum = fabs(dc) + fabs(dr) + fabs(ds);
	if(dsum != 1) continue;
	MRIconst(3,3,3,1,0,mri); // set MRI to 0
	MRIsetVoxVal(mri,1+dc,1+dr,1+ds,0,1);
	sprintf(testname,"single-face %d %d %d\n",dc,dr,ds);
	err = QuadEulerCharChangeCheckReorder(mri, testname, 0+ForceFail);
      }
    }
  }

  // Set each all face voxels. This configuration (before adding
  // the center voxel) has an EC=+8. After adding the center, the
  // EC should become +2, so expecting a drop of 6
  MRIconst(3,3,3,1,0,mri); // set MRI to 0
  for(dc = -1; dc <= 1; dc++){
    for(dr = -1; dr <= 1; dr++){
      for(ds = -1; ds <= 1; ds++){
	dsum = fabs(dc) + fabs(dr) + fabs(ds);
	if(dsum != 1) continue;
	MRIsetVoxVal(mri,1+dc,1+dr,1+ds,0,1);
      }
    }
  }
  err = QuadEulerCharChangeCheckReorder(mri, "all-faces", -6+ForceFail);

  // Set everything but the center voxel. In this configuration,
  // the EC is (9*6+6)+(4*6+2*12+8+8)-(12*6+3*12+12)=+4. When
  // adding the center, the EC changes to +2, so the change should
  // be -2
  MRIconst(3,3,3,1,0,mri); // set MRI to 0
  for(dc = -1; dc <= 1; dc++){
    for(dr = -1; dr <= 1; dr++){
      for(ds = -1; ds <= 1; ds++){
	dsum = fabs(dc) + fabs(dr) + fabs(ds);
	if(dsum == 0) continue;
	MRIsetVoxVal(mri,1+dc,1+dr,1+ds,0,1);
      }
    }
  }
  err = QuadEulerCharChangeCheckReorder(mri, "all", -2+ForceFail);

  // Set a corner and its 3 face neighbors. EC=2 and the change
  // should be 0.
  MRIconst(3,3,3,1,0,mri); // set MRI to 0
  MRIsetVoxVal(mri,0,0,0,0,1); // corner
  MRIsetVoxVal(mri,0,0,1,0,1); // face neighbor of corner
  MRIsetVoxVal(mri,0,1,0,0,1); // face neighbor of corner
  MRIsetVoxVal(mri,0,1,1,0,1); // face neighbor of corner
  err = QuadEulerCharChangeCheckReorder(mri, "cf3", 0+ForceFail);

  // Set a corner, edge, and face neighbor that are 
  // not neighbors to each other. EC=6. 
  MRIconst(3,3,3,1,0,mri); // set MRI to 0
  MRIsetVoxVal(mri,0,0,0,0,1); // corner
  MRIsetVoxVal(mri,1,2,0,0,1); // edge
  MRIsetVoxVal(mri,1,1,2,0,1); // face 
  err = QuadEulerCharChangeCheckReorder(mri, "cef", -2+ForceFail);

  MRIfree(&mri);
  return(err);
}

/*!
  \fn int QuadEulerCharChangeCheckReorder(MRI *mri, char *testname, int decExpected)
  \brief Runs QuadEulerCharChange() on the given 3x3x3 mri after
  permuting the dimensions in several ways. It compares the change in
  EC when setting the center voxel against the passed expected change.
  If they dont agree, then it prints an error message and returns
  non-zero.
 */
int QuadEulerCharChangeCheckReorder(MRI *mri, const char *testname, int decExpected)
{
  int dec,reorder,err;
  MRI *mri2;
  err=0;
  for(reorder=1; reorder <= 7; reorder++){
    switch(reorder){
    case 1: mri2 = MRIreorder(mri, NULL, -1, +2, +3); break; // reverse x
    case 2: mri2 = MRIreorder(mri, NULL, +1, -2, +3); break; // reverse y
    case 3: mri2 = MRIreorder(mri, NULL, +1, +2, -3); break; // reverse z
    case 4: mri2 = MRIreorder(mri, NULL, +2, +1, +3); break; // swap xy
    case 5: mri2 = MRIreorder(mri, NULL, +3, +2, +1); break; // swap xz
    case 6: mri2 = MRIreorder(mri, NULL, +1, +3, +2); break; // swap yz
    case 7: mri2 = MRIreorder(mri, NULL, +3, +1, +2); break; // rotate xyz
    default: 
      printf("ERROR: QuadEulerCharChangeReorder(): reorder option %d\n",reorder);
      return(-1);
    }
    // Determine the change in EC when setting the center voxel
    dec = QuadEulerCharChange(mri2, NULL, 1, 1, 1); 
    if(dec != decExpected) {
      printf("ERROR: QuadEulerCharChangeReorder(): %s reorder=%d, dec=%d, expected %d\n",
	     testname,reorder,dec,decExpected);
      err=1;
    }
    MRIfree(&mri2);
  }
  return(err);
}


/*!
  \fn MRI *MRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm)
  \brief Segments mri_T1 into three labels: 0,
  BRIGHT_BORDER_LABEL=100=border, BRIGHT_LABEL=130=bright that are
  (mostly) outside of the mri_wm mask (wm < WM_MIN_VAL=5).
 */
MRI *MRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm)
{
  int     width, height, depth, x, y, z, nlabeled, nwhite,
          xk, yk, zk, xi, yi, zi;
  BUFTYPE *pwm, val, wm ;
  MRI     *mri_labeled, *mri_tmp ;
  int aMIN_WHITE = ((3*3*3-1)/2); // = 14, hidden parameter, voxelsize dep

  mri_labeled = MRIclone(mri_T1, NULL) ;
  width  = mri_T1->width ;
  height = mri_T1->height ;
  depth  = mri_T1->depth ;

  /* This section creates a binary volume of voxels that are:
     1. Outside of the wm.mgz mask (wm < WM_MIN_VAL)
     2. Have a value of > 125 in eg, brain.finalsurfs
     3. Have < 14 FEC neighbors that are in the wm mask.
        So basically that are not near the wm mask
	125 = hidden parameter
	14 = hidden parameter
   */
  for (z = 0 ; z < depth ; z++)  {
    for (y = 0 ; y < height ; y++)    {
      pwm = &MRIvox(mri_wm, 0, y, z) ;
      for (x = 0 ; x < width ; x++)      {
        val = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
        wm = *pwm++ ;
        /* not white matter and bright (e.g. eye sockets) */
	// WM_MIN_VAL = 5 as of 9/5/19
	// wm < WM_MIN_VAL means outside of the wm.mgz mask
	// val > 125 means "bright" (hidden parameter)
        if ((wm < WM_MIN_VAL) && (val > 125)) {
	  // Found a bright voxel outside of WM
	  // Count the number of of nearest nbrs that are in wm.mgz
	  // Nearest means face, edge, and corner
          nwhite = 0 ;
          for (xk = -1 ; xk <= 1 ; xk++)          {
            xi = mri_T1->xi[x+xk] ;
            for (yk = -1 ; yk <= 1 ; yk++)            {
              yi = mri_T1->yi[y+yk] ;
              for (zk = -1 ; zk <= 1 ; zk++)              {
                zi = mri_T1->zi[z+zk] ;
                if (MRIvox(mri_wm, xi, yi, zi) >= WM_MIN_VAL)
                  nwhite++ ;
              }
            }
          }
	  // Number of WM neighbors must be less than 14 for this
          if (nwhite < aMIN_WHITE){
	    // voxel to be labeled bright
            MRIvox(mri_labeled, x, y, z) = BRIGHT_LABEL ;
	  }
        }
      }
    }
  }
  // At this point mri_labeled is a binary volume with 0 or BRIGHT_LABEL

  // Within a bounding box of the above label, dilate voxels that are > 115
  // in the brain.finalsurfs (mri_T1). This operation expands the binaization
  // above. Hidden parameters: 115 and 10=number of  dilations

  // ATH: Commenting line this out since it's actually doing nothing and causing a
  // memory leak, since the destination volume is not supplied and the return value
  // is not stored. Fixing it appropriately will affect recon-all output.
  // MRIdilateThreshLabel(mri_labeled, mri_T1, NULL, BRIGHT_LABEL, 10,115);
  
  // One dilation followed by one erosion
  MRIclose(mri_labeled, mri_labeled) ;

  /* expand once more to all neighboring voxels that are bright. At
     worst we will erase one voxel of white matter. */
  mri_tmp = MRIdilateThreshLabel(mri_labeled, mri_T1, NULL, BRIGHT_LABEL,1,100);

  // The xor essentially creates a shell (border)
  // xor (0,0->0), (1,0->1), (0,1->1), (1,1)->0
  // "1" means val is between 1 and 255 
  // "0" means val < 1 or > 255
  MRIxor(mri_labeled, mri_tmp, mri_tmp, 1, 255) ;
  MRIreplaceValues(mri_tmp, mri_tmp, 1, BRIGHT_BORDER_LABEL) ;

  // union mean val = MAX(mri_tmp,mri_labeled)
  MRIunion(mri_tmp, mri_labeled, mri_labeled) ;
  nlabeled = MRIvoxelsInLabel(mri_labeled, BRIGHT_LABEL) ;
  printf("MRIfindBrightNonWM(): %d bright non-wm voxels segmented.\n", nlabeled) ;

  /* dilate outwards if exactly 0 */
  // 3 = number of dilation iterations
  // 0 = threshold (exactly 0)
  MRIdilateInvThreshLabel(mri_labeled, mri_T1, mri_labeled, BRIGHT_LABEL, 3, 0) ;

  MRIfree(&mri_tmp) ;
  return(mri_labeled) ;
}


/*!
  \fn MRI *MRIzconcat(MRI *mri1, MRI *mri2, int nskip, MRI *out)
  \brief Concatenates mri2 onto mri1 in the z (slice) direction. The
  first nskip slices are removed from mri2 before concatenation. The
  original app for this was to combine two hires suscept slabs (JonP)
  where the top slice of the bottom slab overlapped with the first
  slice of the top slab. The geometry is such that it agrees with the
  bottom slab (mri1). 
*/
MRI *MRIzconcat(MRI *mri1, MRI *mri2, int nskip, MRI *out)
{
  int nslices = mri1->depth + mri2->depth - nskip;
  if(nskip >= mri2->depth){
    printf("ERROR: MRIzconcat(): nskip=%d >= nslices in mri2 %d\n",nskip,mri2->depth);
    return(NULL);
  }
  if(out == NULL) {
    out = MRIallocSequence(mri1->width, mri1->height, nslices, mri1->type, mri1->nframes);
    MRIcopyHeader(mri1, out);
  }
  if(mri1->width != mri2->width){
    printf("ERROR: MRIzconcat(): mri1 and mri2 mismatch width %d %d\n",mri1->width,mri2->width);
    return (NULL);
  }
  if(mri1->height != mri2->height){
    printf("ERROR: MRIzconcat(): mri1 and mri2 mismatch height %d %d\n",mri1->height,mri2->height);
    return (NULL);
  }
  if(mri1->nframes != out->nframes) {
    printf("ERROR: MRIzconcat(): out nframes mismatch %d %d\n",mri1->nframes, out->nframes);
    return (NULL);
  }
  if(mri1->width != out->width) {
    printf("ERROR: MRIzconcat(): out width mismatch %d %d\n",mri1->width, out->width);
    return (NULL);
  }
  if(mri1->height != out->height) {
    printf("ERROR: MRIzconcat(): out height mismatch %d %d\n",mri1->height, out->height);
    return (NULL);
  }
  if(nslices != out->depth) {
    printf("ERROR: MRIzconcat(): out depth mismatch %d %d\n",nslices, out->depth);
    return (NULL);
  }
  MRIcopyPulseParameters(mri1, out);

  int c;
  for(c=0; c < mri1->width; c++){
    int r,s,f;
    for(r=0; r < mri1->width; r++){
      for(f=0; f < mri1->nframes; f++){
	int sout = 0;
	for(s=0; s < mri1->depth; s++){
	  // Just copy the first volume
	  double v = MRIgetVoxVal(mri1,c,r,s,f);
	  MRIsetVoxVal(out,c,r,sout,f,v);
	  sout++;
	}
	for(s=nskip; s < mri2->depth; s++){
	  // Now append the second in the slice direction
	  double v = MRIgetVoxVal(mri2,c,r,s,f);
	  MRIsetVoxVal(out,c,r,sout,f,v);
	  sout++;
	}
      } //f
    } //r
  } //c

  // Need to fix geometry because we want to simply extend mri1
  MATRIX *M = MRIxfmCRS2XYZ(mri1, 0);
  MRIsetVox2RASFromMatrix(out, M);
  MatrixFree(&M);

  return(out);
}

// See class definition for docs
int FixSubCortMassHA::FixSCM(void)
{
  int c;

  mask = MRIalloc(aseg->width,aseg->height,aseg->depth,MRI_INT);
  MRIcopyHeader(aseg, mask);
  MRIcopyPulseParameters(aseg, mask);

  // create a mask with WM, cortex, and background. 
  // Importantly, hippocampus is not in the mask
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
  for(c=0; c < aseg->width; c++){
    int r,s;
    for(r=0; r < aseg->height; r++){
      for(s=0; s < aseg->depth; s++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	if(segid == 0 || segid == 2 || segid == 3 || segid == 41 || segid == 42){
	  MRIsetVoxVal(mask,c,r,s,0,1);
	}
      }
    }
  }

  // Dilate the mask. This will dilate the mask into hippocampus
  printf("Dilating %d voxels in 3d\n",nDilate);
  for(int n=0; n < nDilate; n++) MRIdilate(mask,mask);

  // Now, remove amyg and ILV and any hippo that is not in the mask
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
  for(c=0; c < aseg->width; c++){
    int r,s;
    for(r=0; r < aseg->height; r++){
      for(s=0; s < aseg->depth; s++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	if(IS_AMYGDALA(segid) || IS_INF_LAT_VENT(segid)){
	  // If in amyg or ventricle, just zero the input regardless of mask
	  MRIsetVoxVal(subcorticalmass,c,r,s,0, 0);
	  continue;
	}
	if(MRIgetVoxVal(mask,c,r,s,0)>0.5) continue;
	if(IS_HIPPO(segid) && MRIgetVoxVal(mask,c,r,s,0)<0.5){
	  // If in hipp but not in the mask, zero
	  MRIsetVoxVal(subcorticalmass,c,r,s,0, 0);
	  continue;
	}
      }
    }
  }
  return(0);
}


/*!
  \fn int MRIfillTriangle(MRI *vol, double p1[3], double p2[3], double p3[3], double dL, double FillVal)
  \brief Fills the voxels that intersect with the give triangle. p?[3]
  are the corners of the triangle where p?[0]=col, p?[1]=row, p?[2] =
  slice. The triangle is subsampled in uniform barycentric coordinate
  with spacing 0<dL<1. dL must be sufficiently small to assure that
  all voxels that intersect with the triangle are filled. The volume
  will be filled with FillVal.
 */
int MRIfillTriangle(MRI *vol, double p1[3], double p2[3], double p3[3], double dL, double FillVal)
{
  double l1, l2, l3;
  double r[3];
  int k,nhits=0;

  for(l1=0; l1 < 1; l1 += dL){
    for(l2=0; l2 < 1; l2 += dL){
      l3 = 1.0 - (l1+l2);
      if(l3 < 0.0 || l3 > 1.0) continue;
      // location of barycentric point
      for(k=0; k<3; k++) r[k] = l1*p1[k] + l2*p2[k] + l3*p3[k];
      int OutOfBounds = MRIindexNotInVolume(vol,r[0],r[1],r[2]);
      if(OutOfBounds) continue;
      MRIsetVoxVal(vol,r[0],r[1],r[2],0,FillVal);
      nhits ++;
    }
  }
  return(nhits);
}

int BinarizeMRI::dump(MRI *invol, const MRI *mask, MRI *outvol)
{
  fprintf(m_debugfp,"BinarizeMRI: ");
  fprintf(m_debugfp,"BinType=%d ",BinType);
  if(BinType == 1) fprintf(m_debugfp," thmin=%g thmax=%g ",thmin,thmax);
  if(BinType == 2) {
    fprintf(m_debugfp,"match N=%d ",(int)matchlist.size());
    for(int n=0; n < matchlist.size(); n++) fprintf(m_debugfp,"%d ",matchlist[n]);
  }
  if(mask) fprintf(m_debugfp,"maskthmin=%g maskthmax=%g ",maskthmin,maskthmax);
  fprintf(m_debugfp,"zce=%d zre=%d zse=%d ",ZeroColEdges,ZeroRowEdges,ZeroSliceEdges);
  fprintf(m_debugfp,"on=%g off=%g ",OnVal,OffVal);
  fprintf(m_debugfp,"fstart=%d fend=%d ",fstart,fend);
  fprintf(m_debugfp,"invert=%d ",invert);
  fprintf(m_debugfp,"DoAbs=%d ",DoAbs);
  fprintf(m_debugfp,"mritype=%d ",mritype);
  fprintf(m_debugfp,"\n");
  if(BinType == 1 && matchlist.size() > 0){
    fprintf(m_debugfp,"WARNING: BinType=1 but matchlist not empty N=%d\n",(int)matchlist.size());
  }
  fflush(m_debugfp);
  return(0);
}


/*!
  \fn int BinarizeMRI::qualifies(const MRI *invol, const MRI *mask,  int c,  int r,  int s,  int f)
  \brief Tests whether a given voxel "qualifies", ie, it meets the
  threshold criteria and (1) is in the mask (if mask is passed), (2)
  is not on an edge (if edges are being zeroed). If you just want to
  test for the mask or edge, then you can pass f = -1.
 */
int BinarizeMRI::qualifies(const MRI *invol, const MRI *mask,  int c,  int r,  int s,  int f)
{
  int qdebug=0;
  //if(c==4 && r==4 && s==6) qdebug=1;
  if(mask){
    double mval = MRIgetVoxVal(mask,c,r,s,0);
    if(mval < maskthmin || mval > maskthmax){
      if(qdebug) printf("q: f=%d mval=%g\n",f,mval);
      return(0);
    }
  }
  if((ZeroColEdges   && (c == 0 || c == invol->width-1))  ||
     (ZeroRowEdges   && (r == 0 || r == invol->height-1)) ||
     (ZeroSliceEdges && (s == 0 || s == invol->depth-1)) ) {
    if(qdebug) printf("q: f=%d edge\n",f);
    return(0);
  }
  if(qdebug) printf("q: f=%d\n",f);
  if(f < 0) return(1); // just checking edge or mask conditions
  double sval = MRIgetVoxVal(invol,c,r,s,f);
  if(DoAbs) sval = fabs(sval);
  int Q = 1;
  if(BinType == 1){
    if(sval < thmin || sval > thmax) Q = 0;
  }
  else {
    Q = 0;
    for(int n=0; n < matchlist.size(); n++){
      if(sval == matchlist[n]){
	Q = 1;
	break;
      }
    }
  }
  if(qdebug){
    printf("  q: sval=%g Q=%d\n",sval,Q);
    fflush(stdout);
  }
  return(Q);
}
MRI *BinarizeMRI::binarize(MRI *invol, const MRI *mask, MRI *outvol)
{
  if(m_debug) dump(invol, mask, outvol);

  if(invol == NULL) {
    printf("ERROR: BinarizeMRI::binarize(): input is NULL\n");
    return(NULL);
  }
  if(BinType == 1 && matchlist.size() > 0){
    printf("WARNING: BinarizeMRI::binarize(): BinType=1 but matchlist not empty N=%d\n",(int)matchlist.size());
    dump(invol, mask, outvol);
    fflush(stdout);
  }
  // Compute the limits of the frames
  if(fstart < 0) fstart = 0;
  if(fend < 0)   fend = invol->nframes;
  if(fend > invol->nframes){
    printf("ERROR: BinarizeMRI::binarize(): fend=%d >= nframes=%d\n",fend,invol->nframes);
    dump(invol, mask, outvol);
    return (NULL);
  }
  if(mask){
    if(MRIdimMismatch(invol, mask, 0)) {
      printf("ERROR: BinarizeMRI::binarize(): input/mask dim mismatch\n");
      return (NULL);
    }
  }
  if(outvol == NULL) {
    int mritypetmp = mritype;
    if(mritype < 0) mritypetmp = MRI_FLOAT;
    outvol = MRIallocSequence(invol->width, invol->height, invol->depth, mritypetmp, fend-fstart);
    if(outvol == NULL) {
      printf("ERROR: BinarizeMRI::binarize(): could not alloc with type %d\n",mritypetmp);
      dump(invol, mask, outvol);
      return(NULL);
    }
    MRIcopyHeader(invol, outvol);
    MRIcopyPulseParameters(invol, outvol);
  } 
  else {
    if(MRIdimMismatch(invol, outvol, 0)) {
      printf("ERROR: BinarizeMRI::binarize(): input/output dim mismatch\n");
      dump(invol, mask, outvol);
      return (NULL);
    }
  }
  // could check output type against On and Off values

  // Swap on and off if inverting (restored below)
  double OnValtmp=OnVal,OffValtmp=OffVal;
  if(invert){
    OnVal = OffValtmp;
    OffVal = OnValtmp;
  }
  if(m_debug) dump(invol, mask, outvol);

  int nhitslocal = 0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : nhitslocal)
  #endif
  for(int c=0; c < invol->width; c++){
    for(int r=0; r < invol->height; r++){
      for(int s=0; s < invol->depth; s++){
	if(! qualifies(invol,mask,c,r,s,-1)){
	  // Skip all frames if:
	  // (1) a mask is specified and vox is not in the mask or 
	  // (2) edges are being zeroed and vox is on the edge
	  for(int f=fstart; f < fend; f++) MRIsetVoxVal(outvol,c,r,s,f,OffVal);
	  continue;
	}
	for(int f=fstart; f < fend; f++){
	  double valset;
	  if(qualifies(invol,NULL,c,r,s,f)) {
	    // can pass NULL as mask above because we already know
	    // about mask at this voxel from above
	    valset = OnVal;
	    nhitslocal++;
	  }
	  else valset = OffVal;
	  MRIsetVoxVal(outvol,c,r,s,f-fstart,valset);
	}
      }
    }
  }
  if(invert){ // restore values
    OnVal = OnValtmp;
    OffVal = OffValtmp;
    nhitslocal = invol->width*invol->height*invol->depth*(fend-fstart)-nhitslocal;
  }

  // This will be the final number of activated voxels in the output
  nhits = nhitslocal;
  if(m_debug){
    fprintf(m_debugfp,"BinarizeMRI: nhits = %d\n",nhits);
    fflush(m_debugfp);
  }
  return(outvol);
}




MRI *DEMorphBinVol::morph(MRI *binvol)
{
  check(binvol,NULL,NULL);
  if(mm_debug) dump(binvol);
  MRI *outvol = copyFrame(binvol,frame);
  if(outvol==NULL) return(NULL);

  nchangestot=0;
  for(int nthmorph=0; nthmorph < nmorph; nthmorph++){
    int nchanges=0;
    MRI *binvol0 = copyFrame(outvol,0);
    #ifdef HAVE_OPENMP
    #pragma omp parallel for reduction(+:nchanges)
    #endif
    for(int c=1; c < binvol->width-1; c++){
      for(int r=1; r < binvol->height-1; r++){
	for(int s=1; s < binvol->depth-1; s++){
	  int binval = MRIgetVoxVal(binvol0,c,r,s,frame);
	  // For dil, look at unhit voxels to determine whether to hit them
	  // For ero, look at hit voxels to determine whether to erase them
	  if((morphtype==1 && binval > 0.5) || (morphtype==2 && binval < 0.5)) continue; 
	  // Count the number of neighbors that are hit (dil) or unhit (ero)
	  int nnbrs = 0;
	  for(int dc=-1; dc <= 1; dc++){
	    for(int dr=-1; dr <= 1; dr++){
	      for(int ds=-1; ds <= 1; ds++){
		if((abs(dc)+abs(dr)+abs(ds))>topo) continue;
		binval = MRIgetVoxVal(binvol0,c+dc,r+dr,s+ds,0);
		if((morphtype==1 && binval > 0.5) || (morphtype==2 && binval < 0.5)) nnbrs++;
		if(nnbrs > nnbrsthresh) break;
	      }//dslice
	      if(nnbrs > nnbrsthresh) break;
	    }//drow
	    if(nnbrs > nnbrsthresh) break;
	  }//dcol
	  if(nnbrs > nnbrsthresh){
	    // If the number of tagged neighbors is greater than thresh, then hit/erase this voxel
	    if(morphtype == 1) MRIsetVoxVal(outvol,c,r,s,0,1);
	    else               MRIsetVoxVal(outvol,c,r,s,0,0);
	    nchanges++;
	  }
	} // slice
      } // row
    } // col
    MRIfree(&binvol0);
    nchangestot += nchanges;
    if(mm_debug) printf("DEMorphBinVol::morph(): pass=%d/%d, nchanges = %d, tot = %d\n",nthmorph+1,nmorph,nchanges,nchangestot);
  }// nthmorph

  return(outvol);
}

MRI *DEMorphBinVol::morph(MRI *invol, double thmin, double thmax, MRI *mask, MRI *outvol)
{
  int err = check(invol, mask, outvol);
  if(err) return(NULL);
  MRI *involuse = invol;
  int free_involuse=0;
  if(invol->nframes > 1){
    involuse = copyFrame(invol,frame);
    free_involuse=1;
  }
  BinarizeMRI bm;
  bm.BinType = 1;
  bm.thmin = thmin;
  bm.thmax = thmax;
  bm.m_debug = mm_debug;
  bm.m_debugfp = mm_debugfp;
  MRI *binvol = bm.binarize(involuse, mask, NULL);
  if(binvol == NULL) return(NULL);
  outvol = morph(binvol);
  MRIfree(&binvol);
  return(outvol);
}
MRI *DEMorphBinVol::morph(MRI *invol, std::vector<int> matchlist, MRI *mask, MRI *outvol)
{
  int err = check(invol, mask, outvol);
  if(err) return(NULL);
  MRI *involuse = invol;
  int free_involuse=0;
  if(invol->nframes > 1){
    involuse = copyFrame(invol,frame);
    free_involuse=1;
  }
  BinarizeMRI bm;
  bm.BinType = 2;
  bm.matchlist = matchlist;
  bm.m_debug = mm_debug;
  bm.m_debugfp = mm_debugfp;
  MRI *binvol = bm.binarize(involuse, mask, NULL);
  if(binvol == NULL) return(NULL);
  outvol = morph(binvol);
  MRIfree(&binvol);
  return(outvol);
}

int DEMorphBinVol::check(MRI *invol, MRI *ubermask, MRI *outvol)
{
  if(invol == NULL) {
    printf("ERROR: DEMorphBinVol::check(): input is NULL\n");
    return(1);
  }
  if(ubermask){
    if(MRIdimMismatch(invol, ubermask, 0)) {
      printf("ERROR: DEMorphBinVol(): input/ubermask dim mismatch\n");
      return(1);
    }
  }
  if(outvol){
    if(MRIdimMismatch(invol, outvol, 0)) {
      printf("ERROR: DEMorphBinVol::check(): input/output dim mismatch\n");
      return(1);
    }
  }
  if(frame < 0 || frame > invol->nframes){
    printf("ERROR: DEMorphBinVol::check(): frame=%d, out of range %d\n",frame,invol->nframes);
    //dump(invol, mask, outvol);
    return(1);
  }
  if(topo < 1 || topo > 3){
    printf("ERROR: DEMorphBinVol::check(): topo=%d, must be 1, 2, or 3\n",topo);
    return(1);
  }
  if(morphtype != 1 && morphtype != 2){
    printf("ERROR: DEMorphBinVol::check(): morphtype=%d, must be 1 or 2\n",morphtype);
    return(1);
  }
  return(0);
}
int DEMorphBinVol::dump(MRI *binvol)
{
  fprintf(mm_debugfp,"DEMorphBinVol: ");
  fprintf(mm_debugfp,"topo=%d ",topo);
  fprintf(mm_debugfp,"morphtype=%d ",morphtype);
  fprintf(mm_debugfp,"nnbrsthresh=%d ",nnbrsthresh);
  fprintf(mm_debugfp,"frame=%d ",frame);
  fprintf(mm_debugfp,"nmorph=%d ",nmorph);
  fprintf(mm_debugfp,"mritype=%d ",mm_mritype);
  fprintf(mm_debugfp,"\n");
  fflush(mm_debugfp);
  return(0);
}

MRI *DEMorphBinVol::copyFrame(MRI *invol, int frameno)
{
  if(frameno < 0 || frameno > invol->nframes){
    printf("ERROR: DEMorphBinVol::copyFrame(): frameno=%d, out of range %d\n",frameno,invol->nframes);
    return(NULL);
  }
  int mritypetmp = mm_mritype;
  if(mm_mritype < 0) mritypetmp = MRI_UCHAR;
  MRI *involfr = MRIallocSequence(invol->width, invol->height, invol->depth, mritypetmp, 1);
  if(involfr == NULL) {
    printf("ERROR: DEMorphBinVol::copyFrame(): could not alloc\n");
    return(NULL);
  }
  MRIcopyHeader(invol, involfr);
  MRIcopyPulseParameters(invol, involfr);
  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for(int c=0; c < invol->width; c++){
    for(int r=0; r < invol->height; r++){
      for(int s=0; s < invol->depth; s++){
	double inval = MRIgetVoxVal(invol,c,r,s,frameno);
	MRIsetVoxVal(involfr,c,r,s,0,inval);
      }
    }
  }
  return(involfr);
}


MRI *SCMstopMask::getmask(void)
{
  // See the class definition for documentation
  printf("SCMstopMask: DoBFS255=%d, DoFilled=%d, DoLV=%d, DoWM255=%d DoWMSA=%d WMSAErodeMM=%g\n",
    DoBFS255,DoFilled,DoLV,DoWM255,DoWMSA,WMSAErodeMM);
  if(!DoBFS255 && !DoFilled && !DoLV && !DoWM255&& !DoWMSA){
    printf("ERROR: SCMstopMask::getmask(): nothing to do\n");
    return(NULL);
  }
  if(!bfs && DoBFS255){
    printf("ERROR: SCMstopMask::getmask(): bfs is NULL\n");
    return(NULL);
  }
  if(!aseg && (DoLV || DoWMSA)){
    printf("ERROR: SCMstopMask::getmask(): aseg is NULL\n");
    return(NULL);
  }
  if(!filledauto && DoFilled){
    printf("ERROR: SCMstopMask::getmask(): filledauto is NULL\n");
    return(NULL);
  }
  if(!filled && DoFilled){
    printf("ERROR: SCMstopMask::getmask(): filled is NULL\n");
    return(NULL);
  }
  if(!wm && DoWM255){
    printf("ERROR: SCMstopMask::getmask(): wm is NULL\n");
    return(NULL);
  }

  MRI *ref=NULL;
  if(aseg) ref = aseg;
  else if(filled) ref=filled;
  else if(wm) ref=wm;
  else if(bfs) ref=bfs;
  if(aseg && MRIdimMismatch(ref, aseg, 0)) {
    printf("ERROR: SCMstopMask::getmask(): aseg dim mismatch\n");
    return(NULL);
  }
  if(filledauto && MRIdimMismatch(ref, filledauto, 0)) {
    printf("ERROR: SCMstopMask::getmask(): filledauto dim mismatch\n");
    return(NULL);
  }
  if(filled && MRIdimMismatch(ref, filled, 0)) {
    printf("ERROR: SCMstopMask::getmask(): filled dim mismatch\n");
    return(NULL);
  }
  if(wm && MRIdimMismatch(ref, wm, 0)) {
    printf("ERROR: SCMstopMask::getmask(): wm dim mismatch\n");
    return(NULL);
  }
  if(bfs && MRIdimMismatch(ref, bfs, 0)) {
    printf("ERROR: SCMstopMask::getmask(): bfs dim mismatch\n");
    return(NULL);
  }

  MRI *bin = MRIallocSequence(ref->width, ref->height, ref->depth, MRI_UCHAR, 1);
  if(bin == NULL) {
    printf("ERROR: SCMstopMask::getmask(): could not alloc\n");
    return(NULL);
  }
  MRIcopyHeader(ref, bin);
  MRIcopyPulseParameters(ref, bin);

  // Create a color table for the bin
  bin->ct = CTABalloc(7);
  CTE *cte;
  cte = (bin->ct->entries[0]);
  sprintf(cte->name,"Background"); cte->ri = 0; cte->gi = 0;cte->bi = 0;
  cte = (bin->ct->entries[1]); 
  sprintf(cte->name,"WMEdits"); cte->ri = 255; cte->gi = 0; cte->bi = 0;
  cte = (bin->ct->entries[2]); 
  sprintf(cte->name,"FilledEdits"); cte->ri = 0;cte->gi = 255;cte->bi = 0;
  cte = (bin->ct->entries[3]); 
  sprintf(cte->name,"BFSEdits"); cte->ri = 0;cte->gi = 0;cte->bi = 255;
  cte = (bin->ct->entries[4]); 
  sprintf(cte->name,"LatVent"); cte->ri = 0;cte->gi = 255;cte->bi = 255;
  cte = (bin->ct->entries[5]); 
  sprintf(cte->name,"WMSA"); cte->ri = 255; cte->gi = 0; cte->bi = 255;

  // Handle WMSA exclusion mask
  MRI *WMSAexclusionMask=NULL;
  if(DoWMSA){
    nWMSAErode = round(WMSAErodeMM/((bfs->xsize+bfs->ysize+bfs->zsize)/3));
    printf("nWMSAErode = %d\n",nWMSAErode);
    if(nWMSAErode>0){
      // To create the exclusion mask, dilate the cortex label by
      // nWMSAErode.  If a WMSA voxel falls into this mask, then
      // ignore it. This is needed because cortex is sometimes
      // incorrectly labeled as WMSAs
      DEMorphBinVol mbv;
      mbv.morphtype = 1;
      mbv.nmorph = nWMSAErode;
      mbv.topo = 1;
      mbv.mm_mritype = MRI_UCHAR;
      std::vector<int> matchlist = {Left_Cerebral_Cortex,Right_Cerebral_Cortex}; //3,42
      WMSAexclusionMask = mbv.morph(aseg, matchlist, NULL,NULL);
    }
  }

  // these have to be local for reduction
  int nbfsedits=0;
  int nfedits=0;
  int nlv=0;
  int nwmedits=0;
  int nwmsa=0;
  int nwmsaexcluded=0;
  int nhits = 0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : nhits,nbfsedits,nfedits,nlv,nwmedits,nwmsaexcluded,nwmsa)
  #endif
  for(int c=0; c < bin->width; c++){
    for(int r=0; r < bin->height; r++){
      for(int s=0; s < bin->depth; s++){
	int hit = 0;
	if(DoWM255){
	  double wmval = MRIgetVoxVal(wm,c,r,s,0);
	  if(wmval == 255){
	    MRIsetVoxVal(bin,c,r,s,0,1);
	    hit = 1;
	    nwmedits++;
	  }
	}
	if(!hit && filledauto && filled){
	  int fa = MRIgetVoxVal(filledauto,c,r,s,0);
	  int f  = MRIgetVoxVal(filled,c,r,s,0);
	  if(fa != f && f > 0.5){
	    MRIsetVoxVal(bin,c,r,s,0,2);
	    hit = 1;
	    nfedits++;
	  }
	}
	if(!hit && DoBFS255){
	  double bfsval = MRIgetVoxVal(bfs,c,r,s,0);
	  if(bfsval == 255){
	    MRIsetVoxVal(bin,c,r,s,0,3);
	    hit = 1;
	    nbfsedits++;
	  }
	}
	if(!hit && DoLV){
	  int asegid = MRIgetVoxVal(aseg,c,r,s,0);
	  if(asegid == 4 || asegid == 43 || asegid == 63 || asegid == 31){
	    MRIsetVoxVal(bin,c,r,s,0,4);
	    hit = 1;
	    nlv++;
	  }
	}
	if(!hit && DoWMSA){
	  int asegid = MRIgetVoxVal(aseg,c,r,s,0);
	  if(asegid == WM_hypointensities || asegid == Left_WM_hypointensities || asegid == Right_WM_hypointensities){
	    int m=0;
	    if(WMSAexclusionMask) {
	      m = MRIgetVoxVal(WMSAexclusionMask,c,r,s,0);
	      if(m) nwmsaexcluded++;
	    }
	    if(!m){
	      MRIsetVoxVal(bin,c,r,s,0,5);
	      hit = 1;
	      nwmsa++;
	    }
	  }
	}
	if(!hit) continue;
	nhits ++;
      } // slice
    } // row
  } // cold
  if(WMSAexclusionMask) MRIfree(&WMSAexclusionMask);
  printf("SCMstopMask::getmask(): nhits=%d nbfsedits=%d, nfedits=%d, nlv=%d, nwmedits=%d, nwmsa=%d, nwmsaex=%d\n",
	 nhits,nbfsedits,nfedits,nlv,nwmedits,nwmsa,nwmsaexcluded);
  nhitslist[0]=nhits;
  nhitslist[1]=nbfsedits;
  nhitslist[2]=nfedits;
  nhitslist[3]=nwmedits;
  nhitslist[4]=nwmsa;
  nhitslist[5]=nwmsaexcluded;
  return(bin);
}


/*!
\fn MRI *MRIapplyEdits(MRI *newauto, MRI *oldauto, MRI *manedit, MRI *outvol)
\brief Copys newauto to outvol but keeping voxels manedit voxels when
they disagree with oldauto. The basic stream assumed here is that the
first pass of a process with create an auto volume and copy it to the
man volume (eg, filled.auto.mgz and filled.mgz) then, if need be, the
man volume is edited. Now the auto and man differ. The next time the
stream is run, several things could happen. (1) the auto is not
regenerated because its dependencies have not changed in which case
nothing needs to change and this function is not run. (2) The program
that generates the auto is rerun because its dependencies have changed
but the auto itself does not change. In this case it generates
"newauto" the same as "oldauto", this function is run, and voxels that
were manually changed are then inserted into the new/oldauto again
resulting in the same man volume as before. However, if upstream
processing does force the newauto to change, then this function will
be able to figure out which voxels need to be copied from manedit to
the newauto by looking at the difference between the oldauto and the
manedit. Can be run in-place.
*/
MRI *MRIapplyEdits(MRI *newauto, MRI *oldauto, MRI *manedit, MRI *outvol)
{
  int c;
  if(outvol==NULL){
    outvol = MRIallocSequence(oldauto->width, oldauto->height, oldauto->depth, oldauto->type, oldauto->nframes);
    if(outvol==NULL) return(NULL);
    MRIcopyHeader(oldauto,outvol);
    MRIcopyPulseParameters(oldauto,outvol);
  }
  int nhits = 0;
  #pragma omp parallel for reduction(+ : nhits)
  for(c=0; c < oldauto->width; c++){
    for(int r=0; r < oldauto->height; r++){
      for(int s=0; s < oldauto->depth; s++){
	for(int f=0; f < oldauto->nframes; f++){
	  double vnew = MRIgetVoxVal(newauto,c,r,s,f);
	  double vold = MRIgetVoxVal(oldauto,c,r,s,f);
	  double vman = MRIgetVoxVal(manedit,c,r,s,f);
	  if(vold==vman) MRIsetVoxVal(outvol,c,r,s,f,vnew);
	  else{
	    MRIsetVoxVal(outvol,c,r,s,f,vman);
	    nhits++;
	  }
	}
      }
    }
  }
  printf("MRIapplyEdits(): copying %d from man to new\n",nhits);
  return(outvol);
}

/*!
  \fn int MRIfixEntoWM(MRI *invol, MRI *entowm, double lhVal, double rhVal)
  \brief This is a bit of a hack to fix the WM around the gyrus
  ambiens (entorhinal cortex). This WM is often so thin that it looks
  like GM. This creates a "bite" in the surface. The idea here is to
  pass a segmentation (entowm) where this area is labeled as {3,4}006
  and/or {3,4}201, where 3=lh and 4=rh as created by mri_entowm_seg.
  Voxels with this label in the entowm volume, are replaced with
  either lhVal or rhVal.  X006 is a large chunk of WM around
  ento. X201 is the gyrus ambiens.  Having these two allows a way to
  modulate the aggressiveness of the fix. Level=1 only fills X006
  (probably not useful); Level=2 fills only X201; Level=3 fills both.
  For example, whe using it to fix the wm.seg.mgz and/or the
  filled.mgz, it is ok to be aggressive as the surface will be refined
  from there. When using it to fix the brain.finalsurfs, then just
  using gyrus ambiens might be called for because the fix will be very
  heavy-handed. Eg, when fixing the wm.mgz or wm.seg.mgz or
  wm.asegedit.mgz, set lhVal=rhVal=255 (this is the same value that is
  used to fill the ventricles, etc).  For brain.finalsurf, use
  lhVal=rhVal=255. For filled, use lhVal=255 and rhVal=127. Currently,
  the entowm.mgz can be created using mri_entowm_seg (a DL seg
  routine). The network was trained from manually labeling this area.
 */
int MRIfixEntoWM(MRI *invol, const MRI *entowm, int Level, double lhVal, double rhVal, int ACJ)
{
  printf("MRIfixEntoWM(): %g %g Level=%d\n",lhVal,rhVal,Level);

  int lhga = 3201;
  int rhga = 4201;
  if(ACJ){
    lhga = 7030;
    rhga = 7031;
  }
  int nchanged=0;
  for(int c=0; c < invol->width; c++){
    for(int r=0; r < invol->height; r++){
      for(int s=0; s < invol->depth; s++){
	int i = MRIgetVoxVal(entowm,c,r,s,0);
	if((i==3006 && (Level==1 || Level==3)) || 
	   (i==lhga && (Level==2 || Level==3)) ){
	  MRIsetVoxVal(invol,c,r,s,0,lhVal);
	  nchanged++;
	  continue;
	}
	if((i==4006 && (Level==1 || Level==3)) || 
	   (i==rhga && (Level==2 || Level==3)) ){
	  MRIsetVoxVal(invol,c,r,s,0,rhVal);
	  nchanged++;
	  continue;
	}
      }
    }
  }
  printf("MRIfixEntoWM(): nchanged = %d\n",nchanged);
  return(nchanged);
}

/*!
  \fn MRI *LabelAmygalaCortalJunction(MRI *seg, int topo, MRI *out)
  \brief Labels the voxel on the boundary bet amyg and ctx. Needed
  to create a white matter strand for surface placement. Seg should
  be an aseg-like volume wiht amyg and ctx labeled before fixing
  with surfaces (likely the aseg.presurf). topo should be 1,2,3;
  likely 3. See also MRIfixEntoWM() for applying the mask.
 */
MRI *LabelAmygalaCortalJunction(MRI *seg, int topo, MRI *out)
{
  if(out==NULL){
    out = MRIallocSequence(seg->width,seg->height,seg->depth,MRI_INT,1);
    MRIcopyHeader(seg, out);
    MRIcopyPulseParameters(seg,out);
    // Create a color table
    out->ct = CTABalloc(7031+1);
    for(int n = 0; n < 7030; n++){
      // delete superfluous entries
      free(out->ct->entries[n]);
      out->ct->entries[n] = NULL;
    }
    for(int n = 0; n < 2; n++){
      // These should match the default color table
      CTE *cte;
      if(n == 0){
	cte = out->ct->entries[7030];
	sprintf(cte->name,"Left-Amygdala-Cortical-Junction");
	cte->ri = 255;
	cte->gi = 85;
	cte->bi = 255;
      } else {
	cte = out->ct->entries[7031];
	sprintf(cte->name,"Right-Amygdala-Cortical-Junction");
	cte->ri = 254;
	cte->gi = 85;
	cte->bi = 255;
      }
    }
  }

  int nhits = 0;
#ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : nhits)
#endif
  for(int c=1; c < seg->width-1; c++){
    for(int r=1; r < seg->height-1; r++){
      for(int s=1; s < seg->depth-1; s++){
	int id = MRIgetVoxVal(seg,c,r,s,0);
	int outlabel;
	if(id == 18) outlabel = 7030; // lh
	else if(id == 54) outlabel = 7031; // rh
	else continue;
	for(int dc = -1; dc < 2; dc++){
	  for(int dr = -1; dr < 2; dr++){
	    for(int ds = -1; ds < 2; ds++){
	      int dsum = abs(dc)+abs(dr)+abs(ds);
	      if(dsum > topo) continue;
	      int did = MRIgetVoxVal(seg,c+dc,r+dr,s+ds,0);
	      if(did != 3 && did != 42) continue;
	      MRIsetVoxVal(out,c+dc,r+dr,s+ds,0,outlabel);
	      nhits ++;
	    }
	  }
	}

      }
    }
  }
  printf("LabelAmygalaCortalJunction((): nhits = %d\n",nhits);
  // Could remove islands and/or fill in holes, but have to separate lh and rh
  //MRI *MRIremoveVolumeIslands(MRI *mask, double thresh, int nKeep, MRI *outvol)
  //tmpvol = MRIremoveVolumeHoles(outvol, 0.5, 1, fillval, NULL);
  return(out);
}

/*!
  \fn MRI *MRIshiftDim(MRI *src, int dim, int nshift, int wrap)
  \brief Shifts the source mri in the given dimension by nshift.
  dim={1,2,3}, nshift can be pos or negative. If wrap != 0, then
  voxels will be wrapped around. Does not change the geometry
  to account for the shift.
 */
MRI *MRIshiftDim(MRI *src, int dim, int nshift, int wrap)
{
  printf("MRIshiftDim() dim = %d nshift = %d  wrap = %d\n",dim,nshift,wrap);

  if(dim < 1 || dim > 3) {
    printf("MRIshiftDim() dim out of range\n");
    return(NULL);
  }
  MRI *out = MRIallocSequence(src->width,src->height,src->depth,src->type,src->nframes);
  if(!out) return(NULL);
  MRIcopyHeader(src,out);
  MRIcopyPulseParameters(src,out);
  if(src->ct) out->ct = CTABdeepCopy(src->ct);

  for(int c=0; c < src->width; c++){
    int c2 = c;
    if(dim==1){
      c2 = c + nshift;
      if(c2 < 0){
	if(!wrap) continue;
	c2 += src->width;
      }
      if(c2 >= src->width){
	if(!wrap) continue;
	c2 -= src->width;
      }
    }
    for(int r=0; r < src->height; r++){
      int r2 = r;
      if(dim==2){
	r2 = r + nshift;
	if(r2 < 0){
	  if(!wrap) continue;
	  r2 += src->width;
	}
	if(r2 >= src->width){
	  if(!wrap) continue;
	  r2 -= src->width;
	}
      }
      for(int s=0; s < src->depth; s++){
	int s2 = s;
	if(dim==3){
	  s2 = s + nshift;
	  if(s2 < 0){
	    if(!wrap) continue;
	    s2 += src->width;
	  }
	  if(s2 >= src->width){
	    if(!wrap) continue;
	    s2 -= src->width;
	  }
	}
	for(int f=0; f < src->nframes; f++){
	  double val = MRIgetVoxVal(src,c,r,s,f);
	  MRIsetVoxVal(out,c2,r2,s2,f,val);
	}
      } // slice
    } // row
  } // col

  return(out);
}
