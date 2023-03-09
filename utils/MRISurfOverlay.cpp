#include "MRISurfOverlay.h"

#include "fio.h"
#include "mri_identify.h"
#include "gifti.h"

// constructor
MRISurfOverlay::MRISurfOverlay()
{
  memset(__foverlay, 0, sizeof(__foverlay));
  __type = FS_MRISURFOVERLAY_UNKNOWN;
  __format = MRI_VOLUME_TYPE_UNKNOWN;
  __stframe = 0;
  __nframes = 1;
  __overlaymri = NULL;

  __nVertices = 0;
  __nFaces = 0;
  __nValsPerVertex = 1;
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
 *    MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE, MRIS_ASCII_FILE, MRIS_VTK_FILE
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


/*
 * The following file formats are supported:
 *   shape measurements (.curv, .sulc, .area, .thickness), .mgh (stats), .gii
 *   (MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE). 
 *   MRI_CURV_FILE is the new CURV format with MAGICNO. = 16777215.
 *
 *   Overlay files can also be in MRIS_ASCII_FILE, MRIS_VTK_FILE, and old CURV formats (read). 
 *
 * The overlay data has 1D morphometry data (vertex-wise measures) or other per-vertex information.
 * The data is read into MRI representation in this class for MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE.
 */
MRI *MRISurfOverlay::read(const char *foverlay, int read_volume, MRIS *outmris)
{
  // check if we support the file format
  __format = getFileFormat(foverlay);
  if (__format == MRI_VOLUME_TYPE_UNKNOWN)
  {
    printf("ERROR MRISurfOverlay::read() - unsupported overlay input type\n");
    return NULL;
  }

  __overlaymri = NULL;
  memcpy(__foverlay, foverlay, sizeof(__foverlay));
  if (__format == MRI_CURV_FILE)
  {
    readCurvatureAsMRI(__foverlay, read_volume);
    if (outmris != NULL)
      __copyOverlay2MRIS(outmris);
  }
  else if (__format == MRI_MGH_FILE)
  {
    __overlaymri = mghRead(__foverlay, read_volume, -1);
    if (outmris != NULL)
      __copyOverlay2MRIS(outmris);
  }
  else if (__format == GIFTI_FILE)
  {
    if (outmris != NULL)
      mrisReadGIFTIfile(__foverlay, outmris);
    else
      __overlaymri = MRISreadGiftiAsMRI(__foverlay, read_volume);
  }
  else if (__format == ASCII_FILE)
  {
    mrisReadAsciiCurvatureFile(outmris, __foverlay);
  }
  else if (__format == VTK_FILE)
  {
    MRISreadVTK(outmris, __foverlay);
  }
  else // assume it is in old curv format
  {
    __readOldCurvature(outmris, __foverlay);
  }

  return __overlaymri;
}


/*
 * private function to read old curvature format
 *
 * here is the file format:
 *   int vnum (nvertices)
 *   int fnum (nfaces)
 *   int curv x nvertices
 */
int MRISurfOverlay::__readOldCurvature(MRIS *outmris, const char *fname)
{
  FILE *fp = fopen(fname, "r");
  if (fp == NULL) 
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::__readOldCurvature(): could not open %s", fname));

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
  if (vnum != outmris->nvertices) {
    fclose(fp);
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISurfOverlay::__readOldCurvature(): incompatible vertex "
                 "number in file %s",
                 fname));
  }

  float curvmin = 10000.0f;
  float curvmax = -10000.0f; /* for compiler warnings */
  for (int k = 0; k < vnum; k++) {
    int i;
    fread2(&i, fp);
    float curv = i / 100.0;

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
  }
  outmris->max_curv = curvmax;
  outmris->min_curv = curvmin;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "done. min=%2.3f max=%2.3f\n", curvmin, curvmax);
  }

  fclose(fp);

  return NO_ERROR;
}


/*
 * private function to copy curvature data to MRIS structure
 * ??? optional flag indicating the overlay data ???
 * ??? shape or stats ???
 */
int MRISurfOverlay::__copyOverlay2MRIS(MRIS *outmris)
{
  int frame = MRISgetReadFrame();
  if (__overlaymri->nframes <= frame) {
    printf("ERROR: attempted to read frame %d from %s\n", frame, __foverlay);
    printf("  but this file only has %d frames.\n", __overlaymri->nframes);
    return (ERROR_BADFILE);
  }

  int nv = __overlaymri->width * __overlaymri->height * __overlaymri->depth;
  if (nv != outmris->nvertices) {
    printf("ERROR: number of vertices in %s does not match surface (%d,%d)\n", __foverlay, nv, outmris->nvertices);
    return 1;
  }

  int vno = 0;
  float curvmin = 10000.0f;
  float curvmax = -10000.0f; /* for compiler warnings */
  for (int s = 0; s < __overlaymri->depth; s++) {
    for (int r = 0; r < __overlaymri->height; r++) {
      for (int c = 0; c < __overlaymri->width; c++) {
        float curv = MRIgetVoxVal(__overlaymri, c, r, s, frame);
        if (s == 0 && r == 0 && c == 0) {
          curvmin = curvmax = curv;
        }
        if (curv > curvmax) {
          curvmax = curv;
        }
        if (curv < curvmin) {
          curvmin = curv;
        }
        outmris->vertices[vno].curv = curv;
        vno++;
      }
    }
  }

  outmris->max_curv = curvmax;
  outmris->min_curv = curvmin;

  return NO_ERROR;
}


/* The method reads MRI_CURV_FILE into MRI structure.
 * 
 * The implementation is taken from MRISreadCurvAsMRI().
 * MRISreadCurvAsMRI() will be changed to use this method.
 */
MRI *MRISurfOverlay::readCurvatureAsMRI(const char *curvfile, int read_volume)
{
  if (!IDisCurv(curvfile))
    return (NULL);

  FILE *fp = fopen(curvfile, "r");

  int magno;
  fread3(&magno, fp);

  __nVertices = freadInt(fp);
  __nFaces = freadInt(fp);
  __nValsPerVertex = freadInt(fp);
  if (__nValsPerVertex != 1) {
    fclose(fp);
    printf("ERROR: MRISurfOverlay::readCurvatureBinary(): %s, vals/vertex %d unsupported\n", curvfile, __nValsPerVertex);
    return (NULL);
  }

  MRI *mri = NULL;
  if (!read_volume) {
    mri = MRIallocHeader(__nVertices, 1, 1, MRI_FLOAT, 1);
    mri->nframes = 1;
    fclose(fp);
    __overlaymri = mri;
    return (mri);
  }

  mri = MRIalloc(__nVertices, 1, 1, MRI_FLOAT);
  for (int k = 0; k < __nVertices; k++) {
    float curv = freadFloat(fp);
    MRIsetVoxVal(mri, k, 0, 0, 0, curv);
  }

  fclose(fp);

  __overlaymri = mri;

  return (mri);
}


/* Output overlay MRI representation to disk. Output inmris if it is given. 
 * file formats supported are MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE,
 * MRIS_ASCII_FILE, MRIS_VTK_FILE.
 */
int MRISurfOverlay::write(const char *fout, MRIS *inmris)
{
  int error = 0;

  MRI *outmri = __overlaymri;
  if (outmri == NULL && inmris == NULL)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - no data available"));

  // check if we support the file format
  int  outtype = getFileFormat(fout);
  if (outtype == MRI_VOLUME_TYPE_UNKNOWN)
    outtype = MRI_CURV_FILE;    // write as MRI_CURV_FILE

  if (outtype == MRI_MGH_FILE)
  {
    if (inmris != NULL)
    {
      MRI *TempMRI = MRIalloc(inmris->nvertices, 1, 1, MRI_FLOAT);
      if (TempMRI == NULL)
        return (ERROR_NOMEMORY);

      for (int vno = 0; vno < inmris->nvertices; vno++) {
        VERTEX *v = &inmris->vertices[vno];
        if (vno == Gdiag_no)
          DiagBreak();

        MRIsetVoxVal(TempMRI, vno, 0, 0, 0, v->curv);
      }

      MRIwrite(TempMRI, fout);
      MRIfree(&TempMRI);
    }
    else
    {
      error = mghWrite(outmri, fout, -1);
    }
  }
  else if (outtype == GIFTI_FILE)
  {
    if (inmris != NULL)
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
      error = MRISwriteGIFTI(inmris, NIFTI_INTENT_SHAPE, fout, __foverlay);
    }
    else
      error = mriWriteGifti(outmri, fout);
  }
  else if (outtype == ASCII_FILE)
  {
    error = mrisWriteAsciiCurvatureFile(inmris, (char*)fout);
  }
  else if (outtype == VTK_FILE)
  {
    MRISwriteVTK(inmris, fout);
    MRISwriteCurvVTK(inmris, fout);
    error = NO_ERROR;
  }
  else if (outtype == MRI_CURV_FILE)
  {
    if (inmris != NULL)
    {
      error = __writeCurvFromMRIS(inmris, fout);
    }
    else
    {
      error = __writeCurvFromMRI(outmri, fout);
    }
  }

  return error;
}


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
