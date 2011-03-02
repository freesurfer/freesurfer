/**
 * @file  mris_expand.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef _MRIS_EXPAND_H
#define _MRIS_EXPAND_H

#define MAXVERTICES 10000
#define MAXFACES 10000

#define TRUE 1
#define FALSE 0

#define SQR(x) ((x)*(x))

/* FUNCTION PROTOTYPES */

static void write_geometry(char *fname) ;
static void read_geometry(char *fname);
static void compute_normals(void);
static void normal_face(int f,float *n);
static void expand_geometry(float mm);

/* TYPE DEFINITIONS */

typedef struct ss_vertex_type_
{
  float x,y,z;
  float nx,ny,nz;
  float xb,yb,zb;
  float nxb,nyb,nzb;
  float ox,oy,oz;
  float dx,dy,dz;
  float mx,my,mz;
  float nc,onc,snc;
  float thickness;
  int vnum;
  int v[10];
  int fnum;
  int f[10];
}
ss_vertex_type;

#endif
