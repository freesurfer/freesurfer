
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <GL/glut.h>
#include <GL/gl.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "macros.h"
#include "oglutil.h"

#if 0
static char vcid[] = "$Id: oglutil.c,v 1.9 1997/12/15 20:41:33 fischl Exp $";
#endif

/*-------------------------------- CONSTANTS -----------------------------*/

#define GRAY               100

#define FIELDSIGN_POS      4    /* blue */
#define FIELDSIGN_NEG      5    /* yellow */
#define BORDER             6 
#define MARKED             7 

#define FOV                (256.0f*SCALE_FACTOR)


/*-------------------------------- PROTOTYPES ----------------------------*/

static void load_brain_coords(float x,float y, float z, float v[]) ;
void   OGLUsetLightingModel(float lite0, float lite1, float lite2, 
             float lite3, float newoffset) ;
static int mrisFindMaxExtents(MRI_SURFACE *mris) ;
static int ogluSetFOV(MRI_SURFACE *mris, double fov) ;

double oglu_fov = FOV ;

static int oglu_scale = 1 ;

int
OGLUinit(MRI_SURFACE *mris, long frame_xdim, long frame_ydim)
{
  int    error ;

  glViewport(0,0,frame_xdim,frame_ydim);
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);

  ogluSetFOV(mris, oglu_fov) ;
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  OGLUsetLightingModel(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f) ;

  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);

  /* now specify where the model is in space */
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  error = glGetError() ;
  if (error != GL_NO_ERROR)
  {
    const char *errstr ;

    errstr = gluErrorString(error) ;
    fprintf(stderr, "GL error %d: %s\n", error, errstr) ;
    exit(1) ;
  }
  return(NO_ERROR) ;
}

#define LIGHT0_BR  0.4 /* was 0.2 */
#define LIGHT1_BR  0.0 
#define LIGHT2_BR  0.6 /* was 0.3 */
#define LIGHT3_BR  0.2 /* was 0.1 */
#define OFFSET 0.25   /* was 0.15 */
void
OGLUsetLightingModel(float lite0, float lite1, float lite2, float lite3, 
                   float newoffset)
{
  float offset = OFFSET ;
  static GLfloat mat0_ambient[] =  { 0.0, 0.0, 0.0, 1.0 };
  static GLfloat mat0_diffuse[] =  { OFFSET, OFFSET, OFFSET, 1.0 };
  static GLfloat mat0_emission[] = { 0.0, 0.0, 0.0, 1.0 };
  static GLfloat light0_diffuse[] = { LIGHT0_BR, LIGHT0_BR, LIGHT0_BR, 1.0 };
  static GLfloat light1_diffuse[] = { LIGHT1_BR, LIGHT1_BR, LIGHT1_BR, 1.0 };
  static GLfloat light2_diffuse[] = { LIGHT2_BR, LIGHT2_BR, LIGHT2_BR, 1.0 };
  static GLfloat light3_diffuse[] = { LIGHT3_BR, LIGHT3_BR, LIGHT3_BR, 1.0 };
  static GLfloat light0_position[] = { 0.0, 0.0, 1.0, 0.0 };
  static GLfloat light1_position[] = { 0.0, 0.0,-1.0, 0.0 };
  static GLfloat light2_position[] = { 0.6, 0.6, 1.6, 0.0 };
  static GLfloat light3_position[] = {-1.0, 0.0, 0.0, 0.0 };
  static GLfloat lmodel_ambient[] =  { 0.0, 0.0, 0.0, 0.0 };

  if (lite0 < 0.0)     lite0 = light0_diffuse[0];
  if (lite1 < 0.0)     lite1 = light1_diffuse[0];
  if (lite2 < 0.0)     lite2 = light2_diffuse[0];
  if (lite3 < 0.0)     lite3 = light3_diffuse[0];
  if (newoffset < 0.0)
    newoffset = offset;
  else
    offset = newoffset ;

  /* material: change DIFFUSE,EMISSION (purpler: EMISSION=0.05*newoffset) */
  mat0_diffuse[0] = mat0_diffuse[1] = mat0_diffuse[2] = newoffset;
  mat0_emission[0] = mat0_emission[1] = mat0_emission[2] = 0.1*newoffset;

  /* lights: change DIFFUSE */
  light0_diffuse[0] = light0_diffuse[1] = light0_diffuse[2] = lite0;
  light1_diffuse[0] = light1_diffuse[1] = light1_diffuse[2] = lite1;
  light2_diffuse[0] = light2_diffuse[1] = light2_diffuse[2] = lite2;
  light3_diffuse[0] = light3_diffuse[1] = light3_diffuse[2] = lite3;

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
    glLoadIdentity();

    glMaterialfv(GL_FRONT, GL_AMBIENT, mat0_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat0_diffuse);
    glMaterialfv(GL_FRONT, GL_EMISSION, mat0_emission);

    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_diffuse);
    glLightfv(GL_LIGHT3, GL_DIFFUSE, light3_diffuse);

    /* might need to move */
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
    glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
    glLightfv(GL_LIGHT3, GL_POSITION, light3_position);

    /* turn off default global 0.2 ambient (inf. viewpnt, not 2-sided OK) */
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    glEnable(GL_LIGHT3);
  glPopMatrix();
}

int
OGLUcompile(MRI_SURFACE *mris, int *marked_vertices, int flags, float cslope)
{
  int          k, n, red, green, blue, error, mv, marked ;
  face_type    *f;
  VERTEX       *v, *vn ;
  float        v1[3], min_curv, max_curv, offset ;

/*  ogluSetFOV(mris) ;*/
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "compiling surface tessellation...") ;

  min_curv = mris->min_curv ;
  max_curv = mris->max_curv ;
  if (-min_curv > max_curv)
    max_curv = -min_curv ;

  for (k=0;k<mris->nfaces;k++) if (!mris->faces[k].ripflag)
  {
    f = &mris->faces[k];
    marked = 0 ;

    if (flags & MESH_FLAG)
    {
      float v2[3], v3[3] ;

      glBegin(GL_LINES);
      glColor3ub(255,255,255);
      for (n=0;n<4;n++)
      {
        v = &mris->vertices[f->v[n]];
        vn = n == 3 ? &mris->vertices[f->v[0]] : &mris->vertices[f->v[n+1]] ;
        load_brain_coords(v->nx,v->ny,v->nz,v1);
        glNormal3fv(v1);
        load_brain_coords(v->x,v->y,v->z,v2);
        load_brain_coords(vn->x, vn->y, vn->z,v3);
        glVertex3fv(v2);
        glVertex3fv(v3);
      }
      glEnd() ;
    }
#if 0
    for (n=0;n<4;n++)
      for (mv = 0 ; marked_vertices[mv] >= 0 ; mv++)
        if (marked_vertices[mv] == f->v[n])
          marked = 1 ;
#else
    for (n=0;n<4;n++)
      if (mris->vertices[f->v[n]].marked)
        marked = 1 ;
#endif

    glBegin(GL_QUADS) ;
    for (n=0;n<4;n++)
    {
      v = &mris->vertices[f->v[n]];

      if (v->ripflag)
        mv = 0 ;

      /* don't display negative flat stuff */
      if (flags & PATCH_FLAG && v->nz < 0) 
        continue ;

      if (marked)
        glColor3ub(0,255,255) ;      /* paint the marked vertex blue */
      else if (v->border)
        glColor3f(240,240,0.0);
      else   /* color it depending on curvature */
      {
#define MIN_GRAY  50
#define MAX_COLOR ((float)(255 - MIN_GRAY))

        red = green = blue = MIN_GRAY ;
        if (FZERO(max_curv))  /* no curvature info */
          red = green = blue = GRAY ;    /* paint it gray */

        /* curvatures seem to be sign-inverted??? */
        if (v->curv > 0)  /* color it red */
        {
          offset = MAX_COLOR*tanh(cslope * (v->curv/max_curv)) ;
          red = MIN_GRAY + nint(offset);
          green = MIN_GRAY - nint(offset) ;
          blue = MIN_GRAY - nint(offset) ;
        }
        else              /* color it green */
        {
          offset = MAX_COLOR*tanh(cslope*(-v->curv/max_curv)) ;
          green = MIN_GRAY + nint(offset);
          red = MIN_GRAY - nint(offset) ;
          blue = MIN_GRAY - nint(offset) ;
        }
        if (red > 255)
          red = 255 ;
        if (green > 255)
          green = 255 ;
        if (blue > 255)
          blue = 255 ;
        if (red < 0)
          red = 0 ;
        if (green < 0)
          green = 0 ;
        if (blue < 0)
          blue = 0 ;
        glColor3ub(red,green,blue);  /* specify the RGB color */
      }
      load_brain_coords(v->nx,v->ny,v->nz,v1);
      glNormal3fv(v1);                /* specify the normal for lighting */
      load_brain_coords(v->x,v->y,v->z,v1);
      glVertex3fv(v1);                /* specify the position of the vertex*/
    }
    glEnd() ;   /* done specifying this polygon */
  }



  if (flags & TP_FLAG) for (mv = 0 ; marked_vertices[mv] >= 0 ; mv++)
  {
    float v2[3], v3[3] ;
#define LINE_LEN  0.1*MAX(MAX(mris->xhi,mris->yhi),mris->zhi)

    if (marked_vertices[mv] < mris->nvertices)
    {
      v = &mris->vertices[marked_vertices[mv]] ;
      glLineWidth(2.0f);
      
      glBegin(GL_LINES);
      glColor3ub(0,255,255);
      load_brain_coords(v->e1x,v->e1y,v->e1z,v1);
      glNormal3fv(v1);
      load_brain_coords(v->x,v->y,v->z,v2);
      load_brain_coords(v->x+LINE_LEN*v->nx,v->y+LINE_LEN*v->ny,
                        v->z+LINE_LEN*v->nz,v3);
      glVertex3fv(v2);
      glVertex3fv(v3);
      glEnd() ;
      
#if 1
      glBegin(GL_LINES);
      glColor3ub(255,255,0);
      load_brain_coords(v->nx,v->ny,v->nz,v1);
      glNormal3fv(v1);
      load_brain_coords(v->x,v->y,v->z,v2);
      load_brain_coords(v->x+LINE_LEN*v->e1x,v->y+LINE_LEN*v->e1y,
                        v->z+LINE_LEN*v->e1z,v3);
      glVertex3fv(v2);
      glVertex3fv(v3);
      glEnd() ;
      
      glBegin(GL_LINES);
      glColor3ub(255,255,0);
      load_brain_coords(v->nx,v->ny,v->nz,v1);
      glNormal3fv(v1);
      load_brain_coords(v->x,v->y,v->z,v2);
      load_brain_coords(v->x+LINE_LEN*v->e2x,v->y+LINE_LEN*v->e2y,
                        v->z+LINE_LEN*v->e2z,v3);
      glVertex3fv(v2);
      glVertex3fv(v3);
      glEnd() ;
#endif
    }
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;

  error = glGetError() ;
  if (error != GL_NO_ERROR)
  {
    const char *errstr ;

    errstr = gluErrorString(error) ;
    fprintf(stderr, "GL error %d: %s\n", error, errstr) ;
    exit(1) ;
  }
  return(NO_ERROR) ;
}

static void
load_brain_coords(float x,float y, float z, float v[])
{
  v[0] = -x;
  v[1] = z;
  v[2] = y;
}
#if SHOW_AXES
  /* build axes */
  glNewList(2, GL_COMPILE) ;

  /* x axis */
  glBegin(GL_LINE) ;
  glColor3ub(meshr,0,0) ;
  load_brain_coords(mris->xlo-30,mris->yctr,mris->zctr,v1);
  glVertex3fv(v1);               
  load_brain_coords(mris->xhi+30,mris->yctr,mris->zctr,v1);
  glColor3ub(meshr,0,0) ;
  glVertex3fv(v1);               
  glEnd() ;

  /* y axis */
  glBegin(GL_LINE) ;
  load_brain_coords(mris->xctr,mris->ylo-30,mris->zctr,v1);
  glColor3ub(meshr,0,0) ;
  glVertex3fv(v1);               
  load_brain_coords(mris->xctr,mris->yhi+30,mris->zctr,v1);
  glColor3ub(meshr,0,0) ;
  glVertex3fv(v1);               
  glEnd() ;

  /* z axis */
  glBegin(GL_LINE) ;
  load_brain_coords(mris->xctr,mris->yctr,mris->zlo-30,v1);
  glColor3ub(meshr,0,0) ;
  glVertex3fv(v1);               
  load_brain_coords(mris->xctr,mris->yctr,mris->zhi+30,v1);
  glColor3ub(meshr,0,0) ;
  glVertex3fv(v1);               
  glEnd() ;

  glEndList() ;
#endif



static int
mrisFindMaxExtents(MRI_SURFACE *mris)
{
  int     vno, xlo, xhi, ylo, yhi, zlo, zhi, x, y, z ;
  VERTEX  *v ;

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo=  10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    x = v->x;
    y = v->y;
    z = v->z;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }

  mris->xlo = xlo ; mris->xhi = xhi ;
  mris->ylo = ylo ; mris->yhi = yhi ;
  mris->zlo = zlo ; mris->zhi = zhi ;
  return(NO_ERROR) ;
}

static int
ogluSetFOV(MRI_SURFACE *mris, double fov)
{
  double zfov, max_dim ;

  mrisFindMaxExtents(mris) ;
  oglu_fov = fov ;
  if (oglu_scale)
  {
    max_dim = MAX(mris->xhi, MAX(mris->yhi, mris->zhi)) ;
    max_dim = MAX(max_dim, MAX(-mris->xlo, MAX(-mris->ylo, -mris->zlo))) ;
    max_dim *= 2 ;
#if 0
    if (3*max_dim < FOV)
      oglu_fov = 3*max_dim ;
    else if (max_dim > FOV) /* try and keep it constant for a range of sizes */
    {
      oglu_fov = FOV*1.5 ;
      while (oglu_fov < max_dim)
        oglu_fov *= 2.0 ;
    }
    else
      oglu_fov = FOV ;
#else
    oglu_fov = max_dim*1.1 ;
#endif
    
#if 0
    if (mris->xlo < -oglu_fov)
      oglu_fov = -1.1f * mris->xlo ;
    if (mris->ylo < -oglu_fov)
      oglu_fov = -1.1f * mris->ylo ;
    if (mris->xhi > oglu_fov)
      oglu_fov = 1.1f * mris->xhi ;
    if (mris->yhi > oglu_fov)
      oglu_fov = 1.1f * mris->yhi ;
#endif
  }

  zfov = 10.0f*oglu_fov ; 
  glOrtho(-oglu_fov, oglu_fov, -oglu_fov, oglu_fov, -zfov, zfov);

  return(NO_ERROR) ;
}

int
OGLUnoscale(void)
{
  oglu_scale = 0 ;
  return(NO_ERROR) ;
}

int
OGLUsetFOV(int fov)
{
  OGLUnoscale() ;
  oglu_fov = (double)fov*SCALE_FACTOR ;
  return(NO_ERROR) ;
}

