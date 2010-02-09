/**
 * @file  tags.h
 * @brief utils for adding tags (meta info) to mgz/h files
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2010/02/09 17:52:29 $
 *    $Revision: 1.19 $
 *
 * Copyright (C) 2005-2010,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
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

#define TAG_OLD_SURF_GEOM           20
#define TAG_SURF_GEOM               21

#define TAG_OLD_MGH_XFORM           30
#define TAG_MGH_XFORM               31
#define TAG_GROUP_AVG_SURFACE_AREA  32

#define TAG_AUTO_ALIGN              33

#define TAG_SCALAR_DOUBLE           40
#define TAG_PEDIR                   41

int TAGreadStart(FILE *fp, long long *plen) ;
int TAGwriteStart(FILE *fp, int tag, long long *phere, long long len) ;
int TAGwriteEnd(FILE *fp, long long there) ;
int TAGskip(FILE *fp, int tag, long long len) ;
int TAGmakeCommandLineString(int argc, char **argv, char *cmd_line) ;
int TAGwriteCommandLine(FILE *fp, char *cmd_line) ;
int TAGwrite(FILE *fp, int tag, void *buf, long long len) ;
int TAGwriteAutoAlign(FILE *fp, MATRIX *M);
MATRIX *TAGreadAutoAlign(FILE *fp);

/* zlib i/o support */
int znzTAGreadStart(znzFile fp, long long *plen) ;
int znzTAGwriteStart(znzFile fp, int tag, long long *phere, long long len) ;
int znzTAGwriteEnd(znzFile fp, long long there) ;
int znzTAGskip(znzFile fp, int tag, long long len) ;
int znzTAGwriteCommandLine(znzFile fp, char *cmd_line) ;
int znzTAGwrite(znzFile fp, int tag, void *buf, long long len) ;
int znzTAGwriteAutoAlign(znzFile fp, MATRIX *M);
MATRIX *znzTAGreadAutoAlign(znzFile fp);

#endif
