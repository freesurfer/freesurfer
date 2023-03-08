#ifndef MRISURFOVERLAY_H
#define MRISURFOVERLAY_H

#include "mri.h"

// list of Freesurfer surface overlays
// the list for FS_MRISURFOVERLAY_STATS* is not complete
// the list should be mapped to INTENTS in GIFTI
#define FS_MRISURFOVERLAY_UNKNOWN                 -1
#define FS_MRISURFOVERLAY_SHAPE                  100
#define FS_MRISURFOVERLAY_SHAPE_CURV             101
#define FS_MRISURFOVERLAY_SHAPE_SULC             102
#define FS_MRISURFOVERLAY_SHAPE_AREA             103
#define FS_MRISURFOVERLAY_SHAPE_THICKNESS        104
#define FS_MRISURFOVERLAY_STATS                  200
#define FS_MRISURFOVERLAY_STATS_GAMMA            201
#define FS_MRISURFOVERLAY_STATS_BETA             202
// ...


/* This class implements methods to read/write Freesurfer overlay files.
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
class MRISurfOverlay
{
public:
  MRISurfOverlay();
  ~MRISurfOverlay();

  MRI *read(const char *foverlay, int read_volume);
  int write(const char *fout, MRI *inmri=NULL);  // MRISwriteCurvature()

  MRI *readCurvatureBinary(const char *curvfile, int read_volume);

  static int getFileFormat(const char *foverlay);

private:
  char __foverlay[1024]; // full path to overlay file
  int  __type;           // not assigned yet, FS_MRISURFOVERLAY*
  int  __format;         // MRI_CURV_FILE, MRI_MGH_FILE, GIFTI_FILE
  int  __stframe;        // starting frame, this can be used when we combine multiple overlays in one MRI
  int  __nframes;        // for functions time series*, nframes can be > 1
  MRI  *__overlaymri;    // overlay data in MRI representation

  int __nVertices;
  int __nFaces;          // # of triangles, we need this for MRI_CURV_FILE output
  int __nValsPerVertex;  // number of values at each vertex, should be 1
};


#endif
