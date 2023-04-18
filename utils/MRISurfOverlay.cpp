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
int MRISurfOverlay::read(int read_volume, MRIS *mris, bool mergegifti, bool splitgifti)
{
  bool usemri = false;
  if (mergegifti || splitgifti)
    usemri = true;

  // read curvature files
  for (int n = 0; n < __nfcurvature; n++)
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

    // read curvature file
    int error = __readOneOverlay(&overlayInfo, read_volume, mris, usemri);
    if (error != NO_ERROR)
      return error;
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
int MRISurfOverlay::__readOneOverlay(OverlayInfoStruct *overlayInfo, int read_volume, MRIS *mris, bool usemri)
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

  overlayInfo->__format = overlayFormat;
#if 0  
  // file in curv format will have MRI_VOLUME_TYPE_UNKNOWN if it is not found
  // MRISreadCurvatureFile() will also try to find the file in surf/ directory
  // don't exit now
  if (overlayFormat == MRI_VOLUME_TYPE_UNKNOWN)
  {
    printf("ERROR MRISurfOverlay::read() - unsupported overlay input type\n");
    return ERROR_BADFILE;
  }
#endif

  if (!usemri)
  {
    MRISreadCurvatureFile(mris, overlayInfo->__foverlay);

    // add to OverlayInfoStruct vector
    __overlayInfo.push_back(*overlayInfo);
  }
  else
  {
    if (overlayFormat == GIFTI_FILE)
    {
      // after read, 
      //   first SHAPE is saved in mris->curv;
      //   first <STATS> is saved in mris->val and mris->stat;
      //   all SHAPE and <STATS> data arrays are saved as multi-frame MRI
      // __overlayInfo is updated in mrisReadGIFTIfile() while reading the intents
      MRISreadCurvatureFile(mris, overlayInfo->__foverlay, NULL, &__overlayInfo);
    }
    else
    {
      std::vector<int> shape{__nVertices, 1, 1, 1};
      overlayInfo->__overlaymri =  new MRI(shape, MRI_FLOAT);
      MRISreadCurvatureFile(mris, overlayInfo->__foverlay, overlayInfo->__overlaymri);

      // add to OverlayInfoStruct vector
      __overlayInfo.push_back(*overlayInfo);
    }
  }

  return NO_ERROR;
}


/* - Output overlay data to disk.
 * - file formats supported are MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE, ASCII_FILE, VTK_FILE.
 * - ASCII_FILE, VTK_FILE, MRI_CURV_FILE can have only one frame of data. Use mris to output.
 * - MRI_MGH_FILE and GIFTI_FILE can have multi-frame data. mris is needed for GIFTI_FILE output.
 */
int MRISurfOverlay::write(MRIS *mris, const char *fout, bool mergegifti, bool splitgifti, const char *outdir)
{
#ifdef __MRISURFOVERLAY_DEBUG
  // debug
  char dbgvol[1024] = {'\0'};
  sprintf(dbgvol, "%s.dbg.mgz", fout);
  printf("[DEBUG] write debug __overlaymri volume as %s\n", dbgvol);
  mghWrite(__overlaymri, dbgvol, -1);
#endif


  int noverlay = __overlayInfo.size();
  if (!mergegifti && noverlay > 1)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - more than one overlay to output"));

  // check if we support the file format
  int  outtype = getFileFormat(fout);
  if (outtype == MRI_VOLUME_TYPE_UNKNOWN)
    outtype = MRI_CURV_FILE;    // write as MRI_CURV_FILE

  if (noverlay > 1 &&
      (outtype == ASCII_FILE   || outtype == VTK_FILE || 
       outtype == MRI_MGH_FILE || outtype == MRI_CURV_FILE))
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - ASCII_FILE/VTK_FILE/MRI_MGH_FILE/MRI_CURV_FILE has more than one overlay"));

  int error = NO_ERROR;
  if (mergegifti)
    error = MRISwriteGIFTICombined(mris, &__overlayInfo, fout);
  else if (splitgifti)
    error = separateGIFTIDataArray(mris, outdir);
  else
    error = MRISwriteCurvature(mris, fout, __overlayInfo[0].__foverlay);

  return error;
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


// the combined GIFTI must contain NIFTI_INTENT_POINTSET and NIFTI_INTENT_TRIANGLE data array
MRIS* MRISurfOverlay::readGIFTICombined(const char *fgifti, MRIS *mris)
{
  return mrisReadGIFTIfile(fgifti, mris, &__overlayInfo);
}


// this function outputs only GIFTI intents for overlays
// output of NIFTI_INTENT_POINTSET and NIFTI_INTENT_TRIANGLE is handled separately
int MRISurfOverlay::separateGIFTIDataArray(MRIS *mris, const char *outdir, const char *fout)
{
  int error = NO_ERROR;

  // write out individual overlay in GIFTI.
  // Later, we can expand the logic to output in other format too.
  int noverlay = __overlayInfo.size();
  for (int n = 0; n < noverlay; n++)
  {
    int giftiintent = getGIFTIIntent(__overlayInfo[n].__type);
    const char *datatype = __overlayInfo[n].__shapedatatype;
    char curv_fname[1024] = {'\0'};
    if (__overlayInfo[n].__foverlay != NULL)
      sprintf(curv_fname, "%s/%s.gii", outdir, __overlayInfo[n].__foverlay);
    else
      sprintf(curv_fname, "%s/intent-dataarray-%d.gii", outdir, n);
    error =  MRISwriteGIFTI(mris, __overlayInfo[n].__overlaymri, giftiintent, curv_fname, curv_fname, datatype);
    if (error != NO_ERROR)
      break; 
  }

  return error;
}
