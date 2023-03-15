#ifndef MRISURFOVERLAY_H
#define MRISURFOVERLAY_H

#include "mri.h"
#include "mrisurf.h"

// list of Freesurfer surface overlays
// the list should be mapped to INTENTS in GIFTI
#define FS_MRISURFOVERLAY_UNKNOWN                 -1
#define FS_MRISURFOVERLAY_SHAPE                  100
#define FS_MRISURFOVERLAY_SHAPE_CURV             101
#define FS_MRISURFOVERLAY_SHAPE_SULC             102
#define FS_MRISURFOVERLAY_SHAPE_AREA             103
#define FS_MRISURFOVERLAY_SHAPE_THICKNESS        104

#define FS_MRISURFOVERLAY_STATS                  200
#define FS_MRISURFOVERLAY_STATS_CORREL           201
#define FS_MRISURFOVERLAY_STATS_TTEST            202
#define FS_MRISURFOVERLAY_STATS_FTEST            203
#define FS_MRISURFOVERLAY_STATS_ZSCORE           204
#define FS_MRISURFOVERLAY_STATS_CHISQ            205
#define FS_MRISURFOVERLAY_STATS_BETA             206
#define FS_MRISURFOVERLAY_STATS_BINOM            207
#define FS_MRISURFOVERLAY_STATS_GAMMA            208
#define FS_MRISURFOVERLAY_STATS_POISSON          209
#define FS_MRISURFOVERLAY_STATS_NORMAL           210
#define FS_MRISURFOVERLAY_STATS_FTEST_NONC       211
#define FS_MRISURFOVERLAY_STATS_CHISQ_NONC       212
#define FS_MRISURFOVERLAY_STATS_LOGISTIC         213
#define FS_MRISURFOVERLAY_STATS_LAPLACE          214
#define FS_MRISURFOVERLAY_STATS_UNIFORM          215
#define FS_MRISURFOVERLAY_STATS_TTEST_NONC       216
#define FS_MRISURFOVERLAY_STATS_WEIBULL          217
#define FS_MRISURFOVERLAY_STATS_CHI              218
#define FS_MRISURFOVERLAY_STATS_INVGAUSS         219
#define FS_MRISURFOVERLAY_STATS_EXTVAL           220
#define FS_MRISURFOVERLAY_STATS_PVAL             221
#define FS_MRISURFOVERLAY_STATS_LOGPVAL          222
#define FS_MRISURFOVERLAY_STATS_LOG10PVAL        223
#define FS_MRISURFOVERLAY_STATS_ESTIMATE         224
// ...

struct OverlayInfoStruct
{
  //char __foverlay[1024]; // full path to overlay file
  const char *__foverlay;
  int  __type;           // FS_MRISURFOVERLAY*
  const char *__datatype;// CurvatureRadial, SulcalDepth, Thickness, Area, (Volume, Jacobian)
  int  __format;         // MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE, ASCII_FILE, VTK_FILE
  int  __stframe;        // starting frame, this can be used when we combine multiple overlays in one MRI
  int  __numframe;       // for functions time series*, nframes can be > 1

  int __nValsPerVertex = 1;  // number of values at each vertex, should be 1  
};


/* This class implements methods to read/write Freesurfer overlay files.
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
class MRISurfOverlay
{
public:
  MRISurfOverlay(MRIS *mris, int noverlay, OverlayInfoStruct *poverlayInfo);
  ~MRISurfOverlay();

  int read(int read_volume, MRIS *mris);
  int write(const char *fout, MRIS *inmris=NULL, bool mergegifti=false);  // MRISwriteCurvature()

  int getNumOverlay() { return __noverlay; }
  const char *getOverlayFilename(int nthOverlay) { return __overlayInfo[nthOverlay].__foverlay; }
  int getOverlayType(int nthOverlay) { return __overlayInfo[nthOverlay].__type; }
  int getGIFTIIntent(int nthOverlay);
  int getFirstFrameNo(int nthOverlay) { return __overlayInfo[nthOverlay].__stframe; }
  int getNumFrames(int nthOverlay) { return __overlayInfo[nthOverlay].__stframe + __overlayInfo[nthOverlay].__numframe; }
  const char *getDataType(int nthOverlay) { return __overlayInfo[nthOverlay].__datatype; }

  // return overlay data in multi-frame MRI
  MRI *getOverlayMRI() { return __overlaymri; }
   
  static int getFileFormat(const char *foverlay);

private:
  int __readOneOverlay(int nthOverlay, int read_volume, MRIS *mris);  // MRISreadCurvatureFile()
  int __copyOverlay2MRIS(MRIS *outmris);
  int __readCurvatureAsMRI(const char *curvfile, int read_volume);
  int __readOldCurvature(const char *fname);

  int __writeCurvFromMRI(MRI *outmri, const char *fout);
  int __writeCurvFromMRIS(MRIS *outmris, const char *fout);

  bool __isStatsData(int nthOverlay);
  bool __isShapeMeasurement(int nthOverlay);

  int __getFrameCount();

  int __nVertices;
  int __nFaces;              // # of triangles, we need this for MRI_CURV_FILE output

  int __noverlay;
  OverlayInfoStruct *__overlayInfo;

  int __currOverlay;
  int __currFrame;

  MRI  *__overlaymri;    // overlay data in MRI representation
};

#endif

