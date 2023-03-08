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


// static member method to return file type for the given file
int MRISurfOverlay::getFileFormat(const char *foverlay)
{
  // check if we support the file format
  int mritype = mri_identify(foverlay);
  if (mritype != MRI_CURV_FILE &&     // it is NEW_VERSION_MAGIC_NUMBER if it has type MRI_CURV_FILE
      mritype != MRI_MGH_FILE  &&
      mritype != GIFTI_FILE) 
  {
    printf("ERROR unsupported overlay file type %d\n", mritype);
  }

  return mritype;
}


/*
 * The following file formats are supported:
 *   shape measurements (.curv, .sulc, .area, .thickness), .mgh (stats), .gii
 *   (MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE). 
 *   MRI_CURV_FILE is the new CURV format with MAGICNO. = 16777215.
 *
 *   Overlay files can also be in MRIS_ASCII_FILE, MRIS_VTK_FILE, and old CURV formats. 
 *   These formats are still handled by MRISreadCurvatureFile()/MRISwriteCurvature().
 *
 * The overlay data has 1D morphometry data (vertex-wise measures) or other per-vertex information.
 * The data is read into MRI representation in this class.
 */
MRI *MRISurfOverlay::read(const char *foverlay, int read_volume)
{
  // check if we support the file format
  __format = getFileFormat(foverlay);
  if (__format == MRI_VOLUME_TYPE_UNKNOWN)
    return NULL;

  memcpy(__foverlay, foverlay, sizeof(__foverlay));
  if (__format == MRI_CURV_FILE)
    readCurvatureBinary(__foverlay, read_volume);
  else if (__format == MRI_MGH_FILE)
    __overlaymri = mghRead(__foverlay, read_volume, -1);
  else if (__format == GIFTI_FILE)
    __overlaymri = MRISreadGiftiAsMRI(__foverlay, read_volume);

  return __overlaymri;
}


/* The method reads MRI_CURV_FILE into MRI structure.
 * 
 * The implementation is taken from MRISreadCurvAsMRI().
 * MRISreadCurvAsMRI() will be changed to use this method.
 */
MRI *MRISurfOverlay::readCurvatureBinary(const char *curvfile, int read_volume)
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


/* Output overlay MRI representation to disk. Output inmri if it is given. 
 * file formats supported are MRI_CURV_FILE, MRI_MGH_FILE, and GIFTI_FILE.
 * 
 * MRIS_ASCII_FILE, MRIS_VTK_FILE, and old CURV formats are handled in MRISwriteCurvature().
 */
int MRISurfOverlay::write(const char *fout, MRI *inmri)
{
  int error = 0;

  MRI *outmri = __overlaymri;
  if (inmri != NULL)
    outmri = inmri;

  if (outmri == NULL)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - empty MRI"));

  int  mritype  = mri_identify(fout);
  if (mritype == MRI_MGH_FILE)
    error = mghWrite(outmri, fout, -1);
  else if (mritype == GIFTI_FILE)
    error = mriWriteGifti(outmri, fout);
  else   // default, write as MRI_CURV_FILE
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
      ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISurfOverlay::write() - could not open %s", fout));

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
  }

  return error;
}


