/**
 * @brief utils for adding tags (meta info) to mgz/h files
 *
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

#include "matrix.h"
#include "znzlib.h"

#ifndef TAGS_H
#define TAGS_H

#define TAG_OLD_COLORTABLE          1
#define TAG_OLD_USEREALRAS          2
#define TAG_CMDLINE                 3
#define TAG_USEREALRAS              4
#define TAG_COLORTABLE              5

#define TAG_GCAMORPH_GEOM           10
#define TAG_GCAMORPH_TYPE           11
#define TAG_GCAMORPH_LABELS         12
#define TAG_GCAMORPH_META           13
#define TAG_GCAMORPH_AFFINE         14

#define TAG_OLD_SURF_GEOM           20
#define TAG_SURF_GEOM               21

#define TAG_OLD_MGH_XFORM           30
#define TAG_MGH_XFORM               31
#define TAG_GROUP_AVG_SURFACE_AREA  32

#define TAG_AUTO_ALIGN              33

#define TAG_SCALAR_DOUBLE           40
#define TAG_PEDIR                   41
#define TAG_MRI_FRAME               42
#define TAG_FIELDSTRENGTH           43
#define TAG_ORIG_RAS2VOX            44


int TAGreadStart(FILE *fp, long long *plen) ;
int TAGwriteStart(FILE *fp, int tag, long long *phere, long long len) ;
int TAGwriteEnd(FILE *fp, long long there) ;
int TAGskip(FILE *fp, int tag, long long len) ;
int TAGmakeCommandLineString(int argc, char **argv, char *cmd_line) ;
int TAGwriteCommandLine(FILE *fp, char *cmd_line) ;
int TAGwrite(FILE *fp, int tag, void *buf, long long len) ;
int TAGwriteMatrix(FILE *fp, MATRIX *M);
MATRIX *TAGreadMatrix(FILE *fp);

/* zlib i/o support */
int znzTAGreadStart(znzFile fp, long long *plen, int tagwithzerolen=0) ;
int znzTAGwriteStart(znzFile fp, int tag, long long *phere, long long len) ;
int znzTAGwriteEnd(znzFile fp, long long there) ;
int znzTAGskip(znzFile fp, int tag, long long len) ;
int znzTAGwriteCommandLine(znzFile fp, char *cmd_line) ;
int znzTAGwrite(znzFile fp, int tag, void *buf, long long len) ;
int znzWriteMatrix(znzFile fp, MATRIX *M, int tag);
MATRIX *znzReadMatrix(znzFile fp);
int znzTAGreadFloat(float *pf, znzFile fp);

int znzWriteAutoAlignMatrix(znzFile fp, MATRIX *M);
MATRIX *znzReadAutoAlignMatrix(znzFile fp);

#endif
