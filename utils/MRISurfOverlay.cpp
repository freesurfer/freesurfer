#include "MRISurfOverlay.h"

#include "fio.h"
#include "mri_identify.h"
#include "gifti.h"

// constructor
MRISurfOverlay::MRISurfOverlay(MRIS *mris, int noverlay, OverlayInfoStruct *poverlayInfo)
{
  __nVertices = mris->nvertices;
  __nFaces    = mris->nfaces;

  __noverlay = noverlay;
  __overlayInfo = new OverlayInfoStruct[__noverlay];
  for (int n = 0; n < __noverlay; n++)
  {
    __overlayInfo[n].__foverlay = poverlayInfo[n].__foverlay;
    __overlayInfo[n].__type = poverlayInfo[n].__type;
    __overlayInfo[n].__format = poverlayInfo[n].__format;
    __overlayInfo[n].__stframe = poverlayInfo[n].__stframe;  // assume one frame for each curvature
    __overlayInfo[n].__numframe = poverlayInfo[n].__numframe;

    __overlayInfo[n].__stframe = n;  // assume one frame for each overlay

    const char *curv_fname = __overlayInfo[n].__foverlay;
    __overlayInfo[n].__shapedatatype = NULL;
    if (strstr(curv_fname, ".thickness"))
      __overlayInfo[n].__shapedatatype = "Thickness";
    else if (strstr(curv_fname, ".curv"))
      __overlayInfo[n].__shapedatatype = "CurvatureRadial";
    else if (strstr(curv_fname, ".sulc"))
      __overlayInfo[n].__shapedatatype = "SulcalDepth";
    else if (strstr(curv_fname, ".area"))
      __overlayInfo[n].__shapedatatype = "Area";
    else if (strstr(curv_fname, ".volume"))
      __overlayInfo[n].__shapedatatype = "Volume";
    else if (strstr(curv_fname, ".jacobian"))
      __overlayInfo[n].__shapedatatype = "Jacobian";
    else
      __overlayInfo[n].__shapedatatype = NULL; 

    __overlayInfo[n].__giftiIntent = getGIFTIIntent(n);

#if 0
    memcpy(__overlayInfo[n], poverlayInfo[n], sizeof(__overlayInfo[n]));
    __overlayInfo[n].__type = poverlayInfo[n]__type;
    __overlayInfo[n].__format = MRI_VOLUME_TYPE_UNKNOWN;
    __overlayInfo[n].__stframe = n;  // assume one frame for each overlay
    __overlayInfo[n].__numframe = 1;

    __overlayInfo[n].__nValsPerVertex = 1;
#endif
  }

  // handle multi overlay in one input file
  __nframes = __getFrameCount();
  std::vector<int> shape{__nVertices, 1, 1, __nframes};
  __overlaymri = new MRI(shape, MRI_FLOAT);
  __currFrame = 0;
} 


// destructor
MRISurfOverlay::~MRISurfOverlay()
{
  // do not free __overlaymri
  // it is returned to caller, the caller is resposible to delete it after they are done
}


/* static member method to return file type for the given file
 *
 *  The following file types are considered valid overlay files:
 *    MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE, ASCII_FILE, VTK_FILE
 * For file types other than those, return MRI_VOLUME_TYPE_UNKNOWN to callers.
 *
 * Notes: This class only handles MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE now.
 */
int MRISurfOverlay::getFileFormat(const char *foverlay)
{
  // check if overlay file type is valid
  int mritype = mri_identify(foverlay);

  int mristype = MRISfileNameType(foverlay);
  if (mristype == MRIS_ASCII_FILE)
    mritype = ASCII_FILE;

  if ((mritype != MRI_CURV_FILE &&     // it is NEW_VERSION_MAGIC_NUMBER if it has type MRI_CURV_FILE
       mritype != MRI_MGH_FILE  && mritype != GIFTI_FILE) &&
      (mritype != ASCII_FILE && mritype != VTK_FILE)) 
    mritype = MRI_VOLUME_TYPE_UNKNOWN;

  return mritype;
}


// read all overlay
int MRISurfOverlay::read(int read_volume, MRIS *mris)
{
  for (int n = 0; n < __noverlay; n++)
  {
    __currOverlay = n;
    __currFrame = n;
    __overlayInfo[n].__stframe = __currFrame;

    int error = __readOneOverlay(n, read_volume, mris);
    if (error != NO_ERROR)
      return error;
  }

  return NO_ERROR;
}


/*
 * The following file formats are supported:
 *   shape measurements (.curv, .sulc, .area, .thickness), stats (.mgh), .gii
 *   (MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE). 
 *   MRI_CURV_FILE is the new CURV format with MAGICNO. = 16777215.
 *
 *   Overlay files can also be in ASCII_FILE, VTK_FILE, and old CURV formats (read). 
 *
 * The overlay data has 1D morphometry data (vertex-wise measures) or other per-vertex information.
 * The data is read into MRI representation in this class for MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE.
 */
int MRISurfOverlay::__readOneOverlay(int nthOverlay, int read_volume, MRIS *mris)
{
  if (mris == NULL)
  {
    printf("ERROR MRISurfOverlay::readOneOverlay() - need MRIS\n");
    return ERROR_BADPARM;
  }

  // check if we support the file format
  int overlayFormat = __overlayInfo[__currFrame].__format;
  if (overlayFormat == MRI_VOLUME_TYPE_UNKNOWN)
    overlayFormat = getFileFormat(__overlayInfo[nthOverlay].__foverlay);

  if (overlayFormat == MRI_VOLUME_TYPE_UNKNOWN)
  {
    printf("ERROR MRISurfOverlay::read() - unsupported overlay input type\n");
    return ERROR_BADFILE;
  }

  __overlayInfo[__currFrame].__format = overlayFormat;

  if (overlayFormat == MRI_CURV_FILE)
  {
    int error = __readCurvatureAsMRI(__overlayInfo[__currFrame].__foverlay, read_volume);
    if (error != NO_ERROR)
      return error;

    __overlayInfo[__currOverlay].__numframe = 1;
    __copyOverlay2MRIS(mris);
  }
  else if (overlayFormat == MRI_MGH_FILE)
  {
    MRI *tempMRI = mghRead(__overlayInfo[__currFrame].__foverlay, read_volume, -1);
    if (tempMRI->width != __nVertices || tempMRI->height != 1 || tempMRI->depth != 1)
    {
      printf("[ERROR] MRISurfOverlay::readOneOverlay() - %s dimensions (%d x %d x %d) doesn't match number of surface vertices (%d)\n",
             __overlayInfo[nthOverlay].__foverlay, tempMRI->width, tempMRI->height, tempMRI->depth, __nVertices);
      return ERROR_BADFILE;
    }

    __overlayInfo[__currOverlay].__numframe = 1;

    // MRI_MGH_FILE can have multiple frames if it is functional statistical data
    if (tempMRI->nframes > 1)
    {
      printf("[INFO] MRISurfOverlay::readOneOverlay() - %s has multiple frames = %d.\n", __overlayInfo[nthOverlay].__foverlay, tempMRI->nframes);

      printf("[INFO] MRISurfOverlay::readOneOverlay() - Each frame will be treated as one overlay.\n");
      __overlayInfo[__currOverlay].__numframe = tempMRI->nframes;
    }

    // copy the data to MRI frames
    int stframe = __overlayInfo[__currOverlay].__stframe;
    int endframe = __overlayInfo[__currOverlay].__stframe + __overlayInfo[__currOverlay].__numframe;
    for (int f = stframe; f < endframe; f++) {
      for (int s = 0; s < __overlaymri->depth; s++) {
        for (int r = 0; r < __overlaymri->height; r++) {
          for (int c = 0; c < __overlaymri->width; c++) {
            float fval = MRIgetVoxVal(tempMRI, c, r, s, f-stframe);
            MRIsetVoxVal(__overlaymri, c, r, s, f, fval);
	  }
        }
      }
    }
    MRIfree(&tempMRI);

    __copyOverlay2MRIS(mris);
  }
  else if (overlayFormat == GIFTI_FILE)
  {
    // after read, 
    //   first SHAPE is saved in mris->curv;
    //   first <STATS> is saved in mris->val and mris->stat;
    //   all SHAPE and <STATS> data arrays are saved as multi-frame MRI
    int currFrame_saved = __currFrame;
    mrisReadGIFTIfile(__overlayInfo[__currFrame].__foverlay, mris, __overlaymri, &__currFrame);

    __overlayInfo[__currOverlay].__numframe = (__currFrame - currFrame_saved);
  }
  else if (overlayFormat == ASCII_FILE)
  {
    mrisReadAsciiCurvatureFile(mris, __overlayInfo[__currFrame].__foverlay, __overlaymri, __currFrame);
    __overlayInfo[__currOverlay].__numframe = 1;
  }
  else if (overlayFormat == VTK_FILE)
  {
    // ??? VTK might have multiple FIELD ???
    MRISreadVTK(mris, __overlayInfo[__currFrame].__foverlay, __overlaymri, __currFrame);
    __overlayInfo[__currOverlay].__numframe = 1;
  }
  else // assume it is in old curv format
  {
    int error = __readOldCurvature(__overlayInfo[__currFrame].__foverlay);
    if (error != NO_ERROR)
      return error;

    __overlayInfo[__currOverlay].__numframe = 1;
    __copyOverlay2MRIS(mris);
  }

  return NO_ERROR;
}


/*
 * private function to read old curvature format
 *
 * here is the file format:
 *   int vnum (nvertices)
 *   int fnum (nfaces)
 *   int curv x nvertices
 */
int MRISurfOverlay::__readOldCurvature(const char *fname)
{
  FILE *fp = fopen(fname, "r");
  if (fp == NULL) 
  {
    printf("ERROR MRISurfOverlay::__readOldCurvature(): could not open %s\n", fname);
    return ERROR_BADFILE;
  }

  int vnum, fnum;

  fread3(&vnum, fp);
  /*
   * if (vnum == NEW_VERSION_MAGIC_NUMBER) {
   *  // If the first 4 bytes int = NEW_VERSION_MAGIC_NUMBER, IDisCurv() returns TRUE; 
   *  // and mri_identify() identifies it as MRI_CURV_FILE, which should be handled earlier already.
   *  // this should be treated as an error if it reaches here
   *  fclose(fp);
   *  return (MRISreadNewCurvatureFile(mris, fname));
   * }
   */

  fread3(&fnum, fp);
  if (vnum != __nVertices) {
    fclose(fp);
    printf("ERROR MRISurfOverlay::__readOldCurvature(): incompatible vertex number in file %s", fname);
    return ERROR_BADFILE;
  }

  //float curvmin = 10000.0f;
  //float curvmax = -10000.0f; /* for compiler warnings */
  for (int k = 0; k < vnum; k++) {
    int i;
    fread2(&i, fp);
    float curv = i / 100.0;

#if 0
    if (k == 0) {
      curvmin = curvmax = curv;
    }
    if (curv > curvmax) {
      curvmax = curv;
    }
    if (curv < curvmin) {
      curvmin = curv;
    }

    outmris->vertices[k].curv = curv;
#endif
    MRIsetVoxVal(__overlaymri, k, 0, 0, __currFrame, curv);
  }

#if 0
  outmris->max_curv = curvmax;
  outmris->min_curv = curvmin;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "done. min=%2.3f max=%2.3f\n", curvmin, curvmax);
  }
#endif

  fclose(fp);

  __overlayInfo[__currOverlay].__numframe = 1;

  return NO_ERROR;
}


/*
 * private function to copy curvature data to MRIS structure
 * shape measurement - copy to MRIS->curv and MRIS->val;
 * stats data        - copy to MRIS->stat and MRIS->val.
 * ??? check ripflag ???
 */
int MRISurfOverlay::__copyOverlay2MRIS(MRIS *outmris)
{
  int frame = MRISgetReadFrame();
  if (__overlayInfo[__currFrame].__numframe <= frame) {
    printf("ERROR: attempted to read frame %d from %s\n", frame, __overlayInfo[__currFrame].__foverlay);
    printf("  but this file only has %d frames.\n", __overlayInfo[__currFrame].__numframe);
    return (ERROR_BADFILE);
  }

  int nv = __overlaymri->width * __overlaymri->height * __overlaymri->depth;
  if (nv != outmris->nvertices) {
    printf("ERROR: number of vertices in %s does not match surface (%d,%d)\n", 
           __overlayInfo[__currFrame].__foverlay, nv, outmris->nvertices);
    return 1;
  }

  bool shapeMeasurement = __isShapeMeasurement(__currFrame);
  bool statsData = __isStatsData(__currFrame);

  int vno = 0;
  float curvmin = 10000.0f;
  float curvmax = -10000.0f; /* for compiler warnings */
  for (int s = 0; s < __overlaymri->depth; s++) {
    for (int r = 0; r < __overlaymri->height; r++) {
      for (int c = 0; c < __overlaymri->width; c++) {
        float f = MRIgetVoxVal(__overlaymri, c, r, s, frame);

        if (shapeMeasurement)
	{ 
          if (s == 0 && r == 0 && c == 0) {
            curvmin = curvmax = f;
          }
          if (f > curvmax) {
            curvmax = f;
          }
          if (f < curvmin) {
            curvmin = f;
          }
          outmris->vertices[vno].curv = f;
          outmris->vertices[vno].val = f;
        }
        else if (statsData)
        {
          outmris->vertices[vno].stat = f;
          outmris->vertices[vno].val = f;
        }

        vno++;
      }
    }
  }

  if (shapeMeasurement)
  {
    outmris->max_curv = curvmax;
    outmris->min_curv = curvmin;
  }

  return NO_ERROR;
}


/* The method reads MRI_CURV_FILE into MRI structure.
 * 
 * The implementation is taken from MRISreadCurvAsMRI().
 * MRISreadCurvAsMRI() will be changed to use this method.
 */
int MRISurfOverlay::__readCurvatureAsMRI(const char *curvfile, int read_volume)
{
  if (!IDisCurv(curvfile))
    return ERROR_BADFILE;

  FILE *fp = fopen(curvfile, "r");

  int magno;
  fread3(&magno, fp);

  __nVertices = freadInt(fp);
  __nFaces = freadInt(fp);

  int nValsPerVertex = freadInt(fp);
  if (nValsPerVertex != 1) {
    fclose(fp);
    printf("ERROR: MRISurfOverlay::readCurvatureBinary(): %s, vals/vertex %d unsupported\n", curvfile, nValsPerVertex);
    return ERROR_BADFILE;
  }

  if (!read_volume) {
    fclose(fp);
    return NO_ERROR;
  }

  for (int k = 0; k < __nVertices; k++) {
    float curv = freadFloat(fp);
    MRIsetVoxVal(__overlaymri, k, 0, 0, __currFrame, curv);
  }

  __overlayInfo[__currOverlay].__numframe = 1;

  fclose(fp);

  return NO_ERROR;
}


/* - Output overlay data to disk.
 * - file formats supported are MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE, ASCII_FILE, VTK_FILE.
 * - ASCII_FILE, VTK_FILE, MRI_CURV_FILE can have only one frame of data. Use mris to output.
 * - MRI_MGH_FILE and GIFTI_FILE can have multi-frame data. mris is needed for GIFTI_FILE output.
 */
int MRISurfOverlay::write(const char *fout, MRIS *mris, bool mergegifti)
{
  int error = 0;

#ifdef __MRISURFOVERLAY_DEBUG
  // debug
  char dbgvol[1024] = {'\0'};  //"/space/papancha/2/users/yh887/fs_test/mris_convert-gifti/dbgvol.mgz";
  sprintf(dbgvol, "%s.dbg.mgz", fout);
  printf("[DEBUG] write debug __overlaymri volume as %s\n", dbgvol);
  mghWrite(__overlaymri, dbgvol, -1);
#endif


  if (!mergegifti && __noverlay > 1)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - more than one overlay to output"));

  MRI *outmri = __overlaymri;
  if (outmri == NULL && mris == NULL)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - no data available"));

  // check if we support the file format
  int  outtype = getFileFormat(fout);
  if (outtype == MRI_VOLUME_TYPE_UNKNOWN)
    outtype = MRI_CURV_FILE;    // write as MRI_CURV_FILE

  if ((outtype == ASCII_FILE || outtype == VTK_FILE || outtype == MRI_CURV_FILE) && __noverlay > 1)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - ASCII_FILE/VTK_FILE/MRI_CURV_FILE has more than one overlay"));

  // ASCII_FILE, VTK_FILE, MRI_CURV_FILE can only have contain one overlay data
  // ASCII_FILE, VTK_FILE only output from MRIS
  if (outtype == MRI_MGH_FILE)
  {
#if 0
    if (mris != NULL)
    {
      MRI *TempMRI = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT);
      if (TempMRI == NULL)
        return (ERROR_NOMEMORY);

      for (int vno = 0; vno < mris->nvertices; vno++) {
        VERTEX *v = &mris->vertices[vno];
        if (vno == Gdiag_no)
          DiagBreak();

        MRIsetVoxVal(TempMRI, vno, 0, 0, 0, v->curv);
      }

      MRIwrite(TempMRI, fout);
      MRIfree(&TempMRI);
    }
    else
#endif
    {
      error = mghWrite(outmri, fout, -1);
    }
  }
  else if (outtype == GIFTI_FILE)
  {
    if (mergegifti)
      MRISwriteGIFTICombined(mris, this, fout);
    else
    {
      /* ???__foverlay will be read again in MRISwriteGIFTI() why???
       * ...
       * if (intent_code == NIFTI_INTENT_SHAPE) {
       *   if (MRISreadCurvatureFile(mris, curv_fname)) {
       *     fprintf(stderr, "MRISwriteGIFTI: couldn't read %s\n", curv_fname);
       *     gifti_free_image(image);
       *     return ERROR_BADFILE;
       *   }
       *   ...
       * }
       */
      int giftiintent = getGIFTIIntent(0);
      error = MRISwriteGIFTI(mris, giftiintent, fout, __overlayInfo[0].__foverlay);
    }
  }
  else if (outtype == ASCII_FILE)
  {
    error = mrisWriteAsciiCurvatureFile(mris, (char*)fout);
  }
  else if (outtype == VTK_FILE)
  {
    MRISwriteVTK(mris, fout);
    MRISwriteCurvVTK(mris, fout);
    error = NO_ERROR;
  }
  else if (outtype == MRI_CURV_FILE)
  {
    error = __writeCurvFromMRIS(mris, fout);
#if 0
    if (inmris != NULL)
    {
      error = __writeCurvFromMRIS(mris, fout);
    }
    else
    {
      error = __writeCurvFromMRI(outmri, fout);
    }
#endif
  }

  return error;
}


// MRI_CURV_FILE can only have one frame of MRI data
int MRISurfOverlay::__writeCurvFromMRI(MRI *outmri, const char *fout)
{
    // check MRI dimensions
    if (__nVertices != outmri->width)
    {
      printf("ERROR number of vertices don't match.\n");
      return 1;
    }
    
    if (outmri->height != 1 || outmri->depth != 1 || outmri->nframes != 1)
    {
      printf("ERROR MRI row, slice, or frame dimension is greater than one\n");
      return 1;
    }

    // The following outputs the new CURV format with MAGICNO. = 16777215.
    // the logic was modified based on MRISwriteCurvature() binary output.
    FILE *fp = fopen(fout, "wb");
    if (fp == NULL) 
      ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::__writeCurvFromMRI() - could not open %s", fout));

    fwrite3(-1, fp); /* same old trick - mark it as new format */
    fwriteInt(__nVertices, fp);
    fwriteInt(__nFaces, fp);
    fwriteInt(1, fp); /* 1 value per vertex */
    
    // loop through MRI crs
    for (int s = 0; s < outmri->depth; s++)
    {
      for (int r = 0; r < outmri->height; r++)
      {
        for (int c = 0; c < outmri->width; c++)
	{
          float curv = MRIgetVoxVal(outmri, c, r, s, 0);
          fwriteFloat(curv, fp);
        }
      }
    }
    fclose(fp);

  return NO_ERROR;
}


int MRISurfOverlay::__writeCurvFromMRIS(MRIS *outmris, const char *fout)
{
  // output curv new format
  FILE *fp = fopen(fout, "wb");
  if (fp == NULL) 
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::__writeCurvFromMRIS(): could not open %s", fout));

  fwrite3(-1, fp); /* same old trick - mark it as new format */
  fwriteInt(outmris->nvertices, fp);
  fwriteInt(outmris->nfaces, fp);
  fwriteInt(1, fp); /* 1 value per vertex */

  for (int k = 0; k < outmris->nvertices; k++) {
    float curv = outmris->vertices[k].curv;
    fwriteFloat(curv, fp);
  }
  fclose(fp);

  return NO_ERROR;
}


int MRISurfOverlay::getGIFTIIntent(int nthOverlay)
{
  int giftiIntent = __overlayInfo[nthOverlay].__giftiIntent;

  if (giftiIntent != NIFTI_INTENT_NONE)
    return giftiIntent;

  int type = __overlayInfo[nthOverlay].__type;
  if (type == FS_MRISURFOVERLAY_SHAPE      ||
      type == FS_MRISURFOVERLAY_SHAPE_CURV || 
      type == FS_MRISURFOVERLAY_SHAPE_SULC ||
      type == FS_MRISURFOVERLAY_SHAPE_AREA ||
      type == FS_MRISURFOVERLAY_SHAPE_THICKNESS)
    giftiIntent = NIFTI_INTENT_SHAPE;
  else if (type == FS_MRISURFOVERLAY_STATS_CORREL)
    giftiIntent = NIFTI_INTENT_CORREL; 
  else if (type == FS_MRISURFOVERLAY_STATS_TTEST)
    giftiIntent = NIFTI_INTENT_TTEST;
  else if (type == FS_MRISURFOVERLAY_STATS_FTEST)
     giftiIntent = NIFTI_INTENT_FTEST;
  else if (type == FS_MRISURFOVERLAY_STATS_ZSCORE)
    giftiIntent = NIFTI_INTENT_ZSCORE; 
  else if (type == FS_MRISURFOVERLAY_STATS_CHISQ)
    giftiIntent = NIFTI_INTENT_CHISQ;
  else if (type == FS_MRISURFOVERLAY_STATS_BETA)
    giftiIntent = NIFTI_INTENT_BETA;
  else if (type == FS_MRISURFOVERLAY_STATS_BINOM)
    giftiIntent = NIFTI_INTENT_BINOM;
  else if (type == FS_MRISURFOVERLAY_STATS_GAMMA)
    giftiIntent = NIFTI_INTENT_GAMMA; 
  else if (type == FS_MRISURFOVERLAY_STATS_POISSON)
    giftiIntent = NIFTI_INTENT_POISSON;
  else if (type == FS_MRISURFOVERLAY_STATS_NORMAL)
    giftiIntent = NIFTI_INTENT_NORMAL; 
  else if (type == FS_MRISURFOVERLAY_STATS_FTEST_NONC)
    giftiIntent = NIFTI_INTENT_FTEST_NONC;
  else if (type == FS_MRISURFOVERLAY_STATS_CHISQ_NONC)
    giftiIntent = NIFTI_INTENT_CHISQ_NONC; 
  else if (type == FS_MRISURFOVERLAY_STATS_LOGISTIC)
    giftiIntent = NIFTI_INTENT_LOGISTIC;
  else if (type == FS_MRISURFOVERLAY_STATS_LAPLACE)
    giftiIntent = NIFTI_INTENT_LAPLACE; 
  else if (type == FS_MRISURFOVERLAY_STATS_UNIFORM)
    giftiIntent = NIFTI_INTENT_UNIFORM;
  else if (type == FS_MRISURFOVERLAY_STATS_TTEST_NONC)
    giftiIntent = NIFTI_INTENT_TTEST_NONC; 
  else if (type == FS_MRISURFOVERLAY_STATS_WEIBULL)
    giftiIntent = NIFTI_INTENT_WEIBULL;
  else if (type == FS_MRISURFOVERLAY_STATS_CHI)
    giftiIntent = NIFTI_INTENT_CHI; 
  else if (type == FS_MRISURFOVERLAY_STATS_INVGAUSS)
    giftiIntent = NIFTI_INTENT_INVGAUSS; 
  else if (type == FS_MRISURFOVERLAY_STATS_EXTVAL)
    giftiIntent = NIFTI_INTENT_EXTVAL;
  else if (type == FS_MRISURFOVERLAY_STATS_PVAL)
    giftiIntent = NIFTI_INTENT_PVAL; 
  else if (type == FS_MRISURFOVERLAY_STATS_LOGPVAL)
    giftiIntent = NIFTI_INTENT_LOGPVAL;
  else if (type == FS_MRISURFOVERLAY_STATS_LOG10PVAL)
    giftiIntent = NIFTI_INTENT_LOG10PVAL; 
  else if (type == FS_MRISURFOVERLAY_STATS_ESTIMATE)
    giftiIntent = NIFTI_INTENT_ESTIMATE;

  __overlayInfo[nthOverlay].__giftiIntent = giftiIntent;
  return giftiIntent;
}


bool MRISurfOverlay::__isShapeMeasurement(int nthOverlay)
{
  bool shapeMeasurement = false;
  int type = __overlayInfo[nthOverlay].__type;
  if (type == FS_MRISURFOVERLAY_SHAPE      ||
      type == FS_MRISURFOVERLAY_SHAPE_CURV || 
      type == FS_MRISURFOVERLAY_SHAPE_SULC ||
      type == FS_MRISURFOVERLAY_SHAPE_AREA ||
      type == FS_MRISURFOVERLAY_SHAPE_THICKNESS)
    shapeMeasurement = true;

  return shapeMeasurement;
}

bool MRISurfOverlay::__isStatsData(int nthOverlay)
{
  bool statsData = false;
  int type = __overlayInfo[nthOverlay].__type;
  if (type == FS_MRISURFOVERLAY_STATS_CORREL     ||
      type == FS_MRISURFOVERLAY_STATS_TTEST      ||
      type == FS_MRISURFOVERLAY_STATS_FTEST      ||
      type == FS_MRISURFOVERLAY_STATS_ZSCORE     ||
      type == FS_MRISURFOVERLAY_STATS_CHISQ      ||
      type == FS_MRISURFOVERLAY_STATS_BETA       ||
      type == FS_MRISURFOVERLAY_STATS_BINOM      ||
      type == FS_MRISURFOVERLAY_STATS_GAMMA      ||
      type == FS_MRISURFOVERLAY_STATS_POISSON    ||
      type == FS_MRISURFOVERLAY_STATS_NORMAL     ||
      type == FS_MRISURFOVERLAY_STATS_FTEST_NONC ||
      type == FS_MRISURFOVERLAY_STATS_CHISQ_NONC ||
      type == FS_MRISURFOVERLAY_STATS_LOGISTIC   ||
      type == FS_MRISURFOVERLAY_STATS_LAPLACE    ||
      type == FS_MRISURFOVERLAY_STATS_UNIFORM    ||
      type == FS_MRISURFOVERLAY_STATS_TTEST_NONC ||
      type == FS_MRISURFOVERLAY_STATS_WEIBULL    ||
      type == FS_MRISURFOVERLAY_STATS_CHI        ||
      type == FS_MRISURFOVERLAY_STATS_INVGAUSS   ||
      type == FS_MRISURFOVERLAY_STATS_EXTVAL     ||
      type == FS_MRISURFOVERLAY_STATS_PVAL       ||
      type == FS_MRISURFOVERLAY_STATS_LOGPVAL    ||
      type == FS_MRISURFOVERLAY_STATS_LOG10PVAL  ||
      type == FS_MRISURFOVERLAY_STATS_ESTIMATE)
    statsData = true;

  return statsData;
}

int MRISurfOverlay::__getFrameCount()
{
  int count = 0;
  for (int n = 0; n < __noverlay; n++)
  {
    int overlayFormat = getFileFormat(__overlayInfo[n].__foverlay);
    if (overlayFormat == MRI_MGH_FILE)
    {
      MRI *tempMRI = MRIreadHeader(__overlayInfo[n].__foverlay, MRI_MGH_FILE);
      if (tempMRI == NULL)
        return count;
  
      count += tempMRI->nframes;
      MRIfree(&tempMRI);
    }
    else if (overlayFormat == GIFTI_FILE)
    {
      // each intent is a frame in MRI
      int intents = getShapeStatIntentCount(__overlayInfo[n].__foverlay);
      count += intents;
    }
    else // MRI_CURV_FILE, ASCII_FILE, VTK_FILE, and old curv format
      count++;
  
    __overlayInfo[n].__format = overlayFormat;
  }
  
  return count;
}
