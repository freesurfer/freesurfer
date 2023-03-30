#include "MRISurfOverlay.h"

#include "fio.h"
#include "mri_identify.h"
#include "gifti.h"

// constructor
MRISurfOverlay::MRISurfOverlay(MRIS *mris, int nfcurvs, const char **fcurvs, int nfstats, const char **fstats)
{
  __nVertices = mris->nvertices;
  __nFaces    = mris->nfaces;

  __nfcurvature = nfcurvs;
  __fcurvatures = fcurvs;
  __nfstats     = nfstats;
  __fstats      = fstats;

  // handle multi overlay in one input file
  __currFrame = 0;
  __nframes = __getFrameCount();
  std::vector<int> shape{__nVertices, 1, 1, __nframes};
  __overlaymri = new MRI(shape, MRI_FLOAT);
} 


// constructor
MRISurfOverlay::MRISurfOverlay(const char *fgifti)
{
  __nVertices = 0;
  __nFaces    = 0;

  __nfcurvature = 1;
  __fcurvatures = new const char*[1];
  __fcurvatures[0] = fgifti;
  __nfstats     = 0;
  __fstats      = NULL;

  // handle multi overlay in one input file
  __currFrame = 0;
  __nframes = __getFrameCount(&__nVertices, &__nFaces);
  std::vector<int> shape{__nVertices, 1, 1, __nframes};
  __overlaymri = new MRI(shape, MRI_FLOAT);
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
  int n = 0;

  // read curvature files
  for (; n < __nfcurvature; n++)
  {
    OverlayInfoStruct overlayInfo;
    overlayInfo.__foverlay = new char[strlen(__fcurvatures[n])+1];
    strcpy(overlayInfo.__foverlay, __fcurvatures[n]);
    overlayInfo.__type = FS_MRISURFOVERLAY_SHAPE;
    overlayInfo.__giftiIntent = NIFTI_INTENT_SHAPE;    //getGIFTIIntent(overlayInfo.__type);

    // assign GIFTI ShapeDataType metadata
    const char *curv_fname = overlayInfo.__foverlay;
    memset(overlayInfo.__shapedatatype, 0, sizeof(overlayInfo.__shapedatatype));
    if (strstr(curv_fname, ".thickness"))
      strcpy(overlayInfo.__shapedatatype, "Thickness");
    else if (strstr(curv_fname, ".curv"))
      strcpy(overlayInfo.__shapedatatype, "CurvatureRadial");
    else if (strstr(curv_fname, ".sulc"))
      strcpy(overlayInfo.__shapedatatype, "SulcalDepth");
    else if (strstr(curv_fname, ".area"))
      strcpy(overlayInfo.__shapedatatype, "Area");
    else if (strstr(curv_fname, ".volume"))
      strcpy(overlayInfo.__shapedatatype, "Volume");
    else if (strstr(curv_fname, ".jacobian"))
      strcpy(overlayInfo.__shapedatatype, "Jacobian");

    overlayInfo.__format = getFileFormat(overlayInfo.__foverlay);
    overlayInfo.__stframe = __currFrame;

    // read curvature file
    int error = __readOneOverlay(&overlayInfo, read_volume, mris);
    if (error != NO_ERROR)
      return error;

    // __currFrame will be updated in __readOneOverlay() (mrisReadGIFTIfile() for GIFTI_FILE)
    //__currFrame += overlayInfo.__numframe;
  }

  // read stats files
  // ...

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
int MRISurfOverlay::__readOneOverlay(OverlayInfoStruct *overlayInfo, int read_volume, MRIS *mris)
{
  if (mris == NULL)
  {
    printf("ERROR MRISurfOverlay::readOneOverlay() - need MRIS\n");
    return ERROR_BADPARM;
  }

  // check if we support the file format
  int overlayFormat = overlayInfo->__format;
  if (overlayFormat == MRI_VOLUME_TYPE_UNKNOWN)
    overlayFormat = getFileFormat(overlayInfo->__foverlay);

  if (overlayFormat == MRI_VOLUME_TYPE_UNKNOWN)
  {
    printf("ERROR MRISurfOverlay::read() - unsupported overlay input type\n");
    return ERROR_BADFILE;
  }

  overlayInfo->__format = overlayFormat;

  if (overlayFormat == MRI_CURV_FILE)
  {
    int error = __readCurvatureAsMRI(overlayInfo->__foverlay, read_volume);
    if (error != NO_ERROR)
      return error;

    overlayInfo->__numframe = 1;

    // add to OverlayInfoStruct vector
    __overlayInfo.push_back(*overlayInfo);

    __copyOverlay2MRIS(overlayInfo, mris);

    __currFrame += overlayInfo->__numframe;
  }
  else if (overlayFormat == MRI_MGH_FILE)
  {
    MRI *tempMRI = mghRead(overlayInfo->__foverlay, read_volume, -1);
    if (tempMRI->width != __nVertices || tempMRI->height != 1 || tempMRI->depth != 1)
    {
      printf("[ERROR] MRISurfOverlay::readOneOverlay() - %s dimensions (%d x %d x %d) doesn't match number of surface vertices (%d)\n",
             overlayInfo->__foverlay, tempMRI->width, tempMRI->height, tempMRI->depth, __nVertices);
      return ERROR_BADFILE;
    }

    overlayInfo->__numframe = 1;

    // MRI_MGH_FILE can have multiple frames if it is functional statistical data
    if (tempMRI->nframes > 1)
    {
      printf("[INFO] MRISurfOverlay::readOneOverlay() - %s has multiple frames = %d.\n", overlayInfo->__foverlay, tempMRI->nframes);

      printf("[INFO] MRISurfOverlay::readOneOverlay() - Each frame will be treated as one overlay.\n");
      overlayInfo->__numframe = tempMRI->nframes;
    }

    // copy the data to MRI frames
    int stframe = overlayInfo->__stframe;
    int endframe = overlayInfo->__stframe + overlayInfo->__numframe;
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

    // add to OverlayInfoStruct vector
    __overlayInfo.push_back(*overlayInfo);

    __copyOverlay2MRIS(overlayInfo, mris);

    __currFrame += overlayInfo->__numframe;
  }
  else if (overlayFormat == GIFTI_FILE)
  {
    // after read, 
    //   first SHAPE is saved in mris->curv;
    //   first <STATS> is saved in mris->val and mris->stat;
    //   all SHAPE and <STATS> data arrays are saved as multi-frame MRI
    // __overlayInfo is updated in mrisReadGIFTIfile() while reading the intents
    mrisReadGIFTIfile(overlayInfo->__foverlay, mris, __overlaymri, &__currFrame, &__overlayInfo);
  }
  else if (overlayFormat == ASCII_FILE)
  {
    mrisReadAsciiCurvatureFile(mris, overlayInfo->__foverlay, __overlaymri, __currFrame);
    overlayInfo->__numframe = 1;

    // add to OverlayInfoStruct vector
    __overlayInfo.push_back(*overlayInfo);

    __currFrame += overlayInfo->__numframe;
  }
  else if (overlayFormat == VTK_FILE)
  {
    // ??? VTK might have multiple FIELD ???
    MRISreadVTK(mris, overlayInfo->__foverlay, __overlaymri, __currFrame);
    overlayInfo->__numframe = 1;

    // add to OverlayInfoStruct vector
    __overlayInfo.push_back(*overlayInfo);

    __currFrame += overlayInfo->__numframe;
  }
  else // assume it is in old curv format
  {
    int error = __readOldCurvature(overlayInfo->__foverlay);
    if (error != NO_ERROR)
      return error;

    overlayInfo->__numframe = 1;

    // add to OverlayInfoStruct vector
    __overlayInfo.push_back(*overlayInfo);

    __copyOverlay2MRIS(overlayInfo, mris);

    __currFrame += overlayInfo->__numframe;
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

  return NO_ERROR;
}


/*
 * private function to copy curvature data to MRIS structure
 * shape measurement - copy to MRIS->curv and MRIS->val;
 * stats data        - copy to MRIS->stat and MRIS->val.
 * ??? check ripflag ???
 */
int MRISurfOverlay::__copyOverlay2MRIS(OverlayInfoStruct *overlayInfo, MRIS *outmris)
{
  int frame = MRISgetReadFrame();
  if (overlayInfo->__numframe <= frame) {
    printf("ERROR: attempted to read frame %d from %s\n", frame, overlayInfo->__foverlay);
    printf("  but this file only has %d frames.\n", overlayInfo->__numframe);
    return (ERROR_BADFILE);
  }

  int nv = __overlaymri->width * __overlaymri->height * __overlaymri->depth;
  if (nv != outmris->nvertices) {
    printf("ERROR: number of vertices in %s does not match surface (%d,%d)\n", 
           overlayInfo->__foverlay, nv, outmris->nvertices);
    return 1;
  }

  bool shapeMeasurement = __isShapeMeasurement(overlayInfo->__type);
  bool statsData = __isStatsData(overlayInfo->__type);

  int stframe = overlayInfo->__stframe;
  int endframe = overlayInfo->__stframe + overlayInfo->__numframe;

  if (shapeMeasurement && overlayInfo->__numframe > 1)
  {
    printf("ERROR: curvature data has multiple frames\n");
    return (ERROR_BADFILE);
  }

  int vno = 0;
  float curvmin = 10000.0f;
  float curvmax = -10000.0f; /* for compiler warnings */
  for (int f = stframe; f < endframe; f++) {
    for (int s = 0; s < __overlaymri->depth; s++) {
      for (int r = 0; r < __overlaymri->height; r++) {
        for (int c = 0; c < __overlaymri->width; c++) {
          float fval = MRIgetVoxVal(__overlaymri, c, r, s, f);

          if (shapeMeasurement)
	  { 
            if (s == 0 && r == 0 && c == 0)
              curvmin = curvmax = fval;

            if (fval > curvmax)
              curvmax = fval;

            if (fval < curvmin)
              curvmin = fval;

            outmris->vertices[vno].curv = fval;
            outmris->vertices[vno].val = fval;
          }
          else if (statsData)
          {
            outmris->vertices[vno].stat = fval;
            outmris->vertices[vno].val = fval;
          }

          vno++;
        }
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


  int noverlay = __overlayInfo.size();
  if (!mergegifti && noverlay > 1)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - more than one overlay to output"));

  MRI *outmri = __overlaymri;
  if (outmri == NULL && mris == NULL)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - no data available"));

  // check if we support the file format
  int  outtype = getFileFormat(fout);
  if (outtype == MRI_VOLUME_TYPE_UNKNOWN)
    outtype = MRI_CURV_FILE;    // write as MRI_CURV_FILE

  if (noverlay > 1 &&
      (outtype == ASCII_FILE   || outtype == VTK_FILE || 
       outtype == MRI_MGH_FILE || outtype == MRI_CURV_FILE))
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - ASCII_FILE/VTK_FILE/MRI_MGH_FILE/MRI_CURV_FILE has more than one overlay"));

  // ASCII_FILE, VTK_FILE, MRI_MGH_FILE, MRI_CURV_FILE can only have contain one overlay data
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
      int giftiintent = getGIFTIIntent(__overlayInfo[0].__type);
      error = MRISwriteGIFTI(mris, giftiintent, fout, __overlayInfo[0].__foverlay);
      // ??? change it to use separateGIFTIDataArray() to write the first overlay ???
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
    error = __writeCurvFromMRI(outmri, fout);
    //error = __writeCurvFromMRIS(mris, fout);
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


int MRISurfOverlay::getGIFTIIntent(int overlaytype)
{
  int giftiIntent = NIFTI_INTENT_NONE;

  if (overlaytype == FS_MRISURFOVERLAY_SHAPE      ||
      overlaytype == FS_MRISURFOVERLAY_SHAPE_CURV || 
      overlaytype == FS_MRISURFOVERLAY_SHAPE_SULC ||
      overlaytype == FS_MRISURFOVERLAY_SHAPE_AREA ||
      overlaytype == FS_MRISURFOVERLAY_SHAPE_THICKNESS)
    giftiIntent = NIFTI_INTENT_SHAPE;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_CORREL)
    giftiIntent = NIFTI_INTENT_CORREL; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_TTEST)
    giftiIntent = NIFTI_INTENT_TTEST;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_FTEST)
     giftiIntent = NIFTI_INTENT_FTEST;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_ZSCORE)
    giftiIntent = NIFTI_INTENT_ZSCORE; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_CHISQ)
    giftiIntent = NIFTI_INTENT_CHISQ;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_BETA)
    giftiIntent = NIFTI_INTENT_BETA;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_BINOM)
    giftiIntent = NIFTI_INTENT_BINOM;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_GAMMA)
    giftiIntent = NIFTI_INTENT_GAMMA; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_POISSON)
    giftiIntent = NIFTI_INTENT_POISSON;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_NORMAL)
    giftiIntent = NIFTI_INTENT_NORMAL; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_FTEST_NONC)
    giftiIntent = NIFTI_INTENT_FTEST_NONC;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_CHISQ_NONC)
    giftiIntent = NIFTI_INTENT_CHISQ_NONC; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_LOGISTIC)
    giftiIntent = NIFTI_INTENT_LOGISTIC;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_LAPLACE)
    giftiIntent = NIFTI_INTENT_LAPLACE; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_UNIFORM)
    giftiIntent = NIFTI_INTENT_UNIFORM;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_TTEST_NONC)
    giftiIntent = NIFTI_INTENT_TTEST_NONC; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_WEIBULL)
    giftiIntent = NIFTI_INTENT_WEIBULL;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_CHI)
    giftiIntent = NIFTI_INTENT_CHI; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_INVGAUSS)
    giftiIntent = NIFTI_INTENT_INVGAUSS; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_EXTVAL)
    giftiIntent = NIFTI_INTENT_EXTVAL;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_PVAL)
    giftiIntent = NIFTI_INTENT_PVAL; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_LOGPVAL)
    giftiIntent = NIFTI_INTENT_LOGPVAL;
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_LOG10PVAL)
    giftiIntent = NIFTI_INTENT_LOG10PVAL; 
  else if (overlaytype == FS_MRISURFOVERLAY_STATS_ESTIMATE)
    giftiIntent = NIFTI_INTENT_ESTIMATE;

  return giftiIntent;
}


bool MRISurfOverlay::__isShapeMeasurement(int overlaytype)
{
  bool shapeMeasurement = false;
  if (overlaytype == FS_MRISURFOVERLAY_SHAPE      ||
      overlaytype == FS_MRISURFOVERLAY_SHAPE_CURV || 
      overlaytype == FS_MRISURFOVERLAY_SHAPE_SULC ||
      overlaytype == FS_MRISURFOVERLAY_SHAPE_AREA ||
      overlaytype == FS_MRISURFOVERLAY_SHAPE_THICKNESS)
    shapeMeasurement = true;

  return shapeMeasurement;
}

bool MRISurfOverlay::__isStatsData(int overlaytype)
{
  bool statsData = false;
  if (overlaytype == FS_MRISURFOVERLAY_STATS_CORREL     ||
      overlaytype == FS_MRISURFOVERLAY_STATS_TTEST      ||
      overlaytype == FS_MRISURFOVERLAY_STATS_FTEST      ||
      overlaytype == FS_MRISURFOVERLAY_STATS_ZSCORE     ||
      overlaytype == FS_MRISURFOVERLAY_STATS_CHISQ      ||
      overlaytype == FS_MRISURFOVERLAY_STATS_BETA       ||
      overlaytype == FS_MRISURFOVERLAY_STATS_BINOM      ||
      overlaytype == FS_MRISURFOVERLAY_STATS_GAMMA      ||
      overlaytype == FS_MRISURFOVERLAY_STATS_POISSON    ||
      overlaytype == FS_MRISURFOVERLAY_STATS_NORMAL     ||
      overlaytype == FS_MRISURFOVERLAY_STATS_FTEST_NONC ||
      overlaytype == FS_MRISURFOVERLAY_STATS_CHISQ_NONC ||
      overlaytype == FS_MRISURFOVERLAY_STATS_LOGISTIC   ||
      overlaytype == FS_MRISURFOVERLAY_STATS_LAPLACE    ||
      overlaytype == FS_MRISURFOVERLAY_STATS_UNIFORM    ||
      overlaytype == FS_MRISURFOVERLAY_STATS_TTEST_NONC ||
      overlaytype == FS_MRISURFOVERLAY_STATS_WEIBULL    ||
      overlaytype == FS_MRISURFOVERLAY_STATS_CHI        ||
      overlaytype == FS_MRISURFOVERLAY_STATS_INVGAUSS   ||
      overlaytype == FS_MRISURFOVERLAY_STATS_EXTVAL     ||
      overlaytype == FS_MRISURFOVERLAY_STATS_PVAL       ||
      overlaytype == FS_MRISURFOVERLAY_STATS_LOGPVAL    ||
      overlaytype == FS_MRISURFOVERLAY_STATS_LOG10PVAL  ||
      overlaytype == FS_MRISURFOVERLAY_STATS_ESTIMATE)
    statsData = true;

  return statsData;
}

int MRISurfOverlay::__getFrameCount(int *nVertices, int *nFaces)
{
  int count = 0;
  for (int n = 0; n < __nfcurvature; n++)
  {
    int overlayFormat = getFileFormat(__fcurvatures[n]);
    if (overlayFormat == MRI_MGH_FILE)
    {
      MRI *tempMRI = MRIreadHeader(__fcurvatures[n], MRI_MGH_FILE);
      if (tempMRI == NULL)
        return count;
  
      count += tempMRI->nframes;
      MRIfree(&tempMRI);
    }
    else if (overlayFormat == GIFTI_FILE)
    {
      // count SHAPE & <STATS> intents
      // each intent is a frame in MRI
      int intents = getShapeStatIntentCount(__fcurvatures[n], nVertices, nFaces);
      count += intents;
    }
    else // MRI_CURV_FILE, ASCII_FILE, VTK_FILE, and old curv format
      count++;
  }
  
#if 0
  for (int n = 0; n < __nfstats; n++)
  {
    int overlayFormat = getFileFormat(__fstats[n]);
    if (overlayFormat == MRI_MGH_FILE)
    {
      MRI *tempMRI = MRIreadHeader(__fstats[n], MRI_MGH_FILE);
      if (tempMRI == NULL)
        return count;
  
      count += tempMRI->nframes;
      MRIfree(&tempMRI);
    }
    else if (overlayFormat == GIFTI_FILE)
    {
      // count <STATS> intents
      // each intent is a frame in MRI
      int intents = getShapeStatIntentCount(__fstats[n]);
      count += intents;
    }
    else // MRI_CURV_FILE, ASCII_FILE, VTK_FILE, and old curv format
      count++;
  }
#endif

  return count;
}


MRIS* MRISurfOverlay::readGIFTICombined(const char *fgifti, MRIS *mris)
{
  return mrisReadGIFTIfile(fgifti, mris, __overlaymri, &__currFrame, &__overlayInfo);
}


int MRISurfOverlay::separateGIFTIDataArray(MRIS *mris, const char *outdir, const char *fout)
{
#if 0  // surface will be written out by the caller mris_convert
  int error = MRISwriteGIFTI(mris, NIFTI_INTENT_POINTSET, fout, NULL);
  if (error != NO_ERROR)
    return error;
#endif

  int error = NO_ERROR;

  // write out individual overlay in GIFTI.
  // Later, we can expand the logic to output in other format too.
  int noverlay = __overlayInfo.size();
  for (int n = 0; n < noverlay; n++)
  {
    int giftiintent = getGIFTIIntent(__overlayInfo[n].__type);
    int stframe  = __overlayInfo[n].__stframe;
    int endframe = __overlayInfo[n].__stframe + __overlayInfo[n].__numframe;
    const char *datatype = __overlayInfo[n].__shapedatatype;
    char curv_fname[1024] = {'\0'};
    if (__overlayInfo[n].__foverlay != NULL)
      sprintf(curv_fname, "%s/%s.gii", outdir, __overlayInfo[n].__foverlay);
    else
      sprintf(curv_fname, "%s/intent-dataarray-%d.gii", outdir, n);
    error =  MRISwriteGIFTI(mris, __overlaymri, stframe, endframe, giftiintent, curv_fname, curv_fname, datatype);
    if (error != NO_ERROR)
      break; 
  }

  return error;
}
