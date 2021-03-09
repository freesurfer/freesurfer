/**
 * @brief OpenGL utils
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

#ifdef HAVE_OPENGL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#ifdef HAVE_APPLE_OPENGL_FRAMEWORK
#  include <OpenGL/glu.h>
#  include <OpenGL/gl.h>
#else
#  include <GL/glu.h>
#  include <GL/gl.h>
#endif
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "macros.h"
#include "TexFont.h"
#include "oglutil.h"

#if 0
#endif

/*-------------------------------- CONSTANTS -----------------------------*/

#define RGBcolor(R,G,B)  glColor3ub((GLubyte)(R),(GLubyte)(G),(GLubyte)(B))
#define GRAY               100

#define FIELDSIGN_POS      4    /* blue */
#define FIELDSIGN_NEG      5    /* yellow */
#define BORDER             6
#define MARKED             7

#define FOV                (256.0f*SCALE_FACTOR)
#define COLSCALEBAR_XPOS    1.7
#define COLSCALEBAR_YPOS   -1.3
#define COLSCALEBAR_WIDTH   0.15
#define COLSCALEBAR_HEIGHT  1.0

static float fov = FOV ;
static float colscalebar_xpos  = COLSCALEBAR_XPOS ;
static float colscalebar_ypos = COLSCALEBAR_YPOS ;
static float colscalebar_width = COLSCALEBAR_WIDTH ;
static float colscalebar_height = COLSCALEBAR_HEIGHT ;
static float sf=0.55;      /* initial scale factor */


/*-------------------------------- PROTOTYPES ----------------------------*/

static void draw_colscalebar(void) ;
static void draw_colscalebar_time(int flags) ;
static int  read_environment_variables(void) ;
static int set_color(float val,float curv, int flags) ;
static int set_stat_color(float f, float *rp, float *gp, float *bp,
                          float tmpoffset) ;
static int set_stat_color_time(float f, float *rp, float *gp, float *bp,
                               float tmpoffset) ;
static void load_brain_coords(float x,float y, float z, float v[]) ;
static int mrisFindMaxExtents(MRI_SURFACE *mris) ;
static int ogluSetFOV(MRI_SURFACE *mris, double fov) ;

double oglu_fov = FOV ;
static double oglu_coord_thickness = 1.0 ;   /* sigma of coord line */
static double oglu_coord_spacing = 18.0 ;    /* spacing between coord lines */

static int oglu_scale = 1 ;

static double cvfact = 1.5;
static double offset = 0.25 ;
static double fcurv = 0.0f ;

double  fthresh = 0.0;
double  pre_fthresh = 0.0 ;
double  fmid = 1.0 ;
double  fslope = 1.0 ;
double  time_fthresh = 0.9 ;

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

    errstr = (const char *)gluErrorString(error) ;
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
  static GLfloat mat0_ambient[] =
    {
      0.0, 0.0, 0.0, 1.0
    };
  static GLfloat mat0_diffuse[] =
    {
      OFFSET, OFFSET, OFFSET, 1.0
    };
  static GLfloat mat0_emission[] =
    {
      0.0, 0.0, 0.0, 1.0
    };
  static GLfloat light0_diffuse[] =
    {
      LIGHT0_BR, LIGHT0_BR, LIGHT0_BR, 1.0
    };
  static GLfloat light1_diffuse[] =
    {
      LIGHT1_BR, LIGHT1_BR, LIGHT1_BR, 1.0
    };
  static GLfloat light2_diffuse[] =
    {
      LIGHT2_BR, LIGHT2_BR, LIGHT2_BR, 1.0
    };
  static GLfloat light3_diffuse[] =
    {
      LIGHT3_BR, LIGHT3_BR, LIGHT3_BR, 1.0
    };
  static GLfloat light0_position[] =
    {
      0.0, 0.0, 1.0, 0.0
    };
  static GLfloat light1_position[] =
    {
      0.0, 0.0,-1.0, 0.0
    };
#if 0
  static GLfloat light2_position[] =
    {
      0.6, 0.6, 1.6, 0.0
    };
#else
  static GLfloat light2_position[] =
    {
      0.6, 0.6, 0.6, 0.0
    };
#endif
  static GLfloat light3_position[] =
    {
      -1.0, 0.0, 0.0, 0.0
    };
  static GLfloat lmodel_ambient[] =
    {
      0.0, 0.0, 0.0, 0.0
    };

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
  mat0_emission[0] = mat0_emission[1] = 
    mat0_emission[2] = (GLfloat)0.1*newoffset;

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

#define COORD_RED   255
#define COORD_BLUE  255
#define COORD_GREEN 255

static float min_gray = 0.2f ;
static float brightness = 255.0f ;

int
OGLUcompile(MRI_SURFACE *mris, int *marked_vertices, int flags, float cslope)
{
  int          k, n, red, green, blue, mv, marked, vno ;
  int          error;
  face_type    *f;
  VERTEX       *v, *vn ;
  float        v1[3], min_curv, max_curv, coord_coef = 0.0, color_val;

  TexFont *txf;
  char text[256] ;

  read_environment_variables() ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "compiling surface tessellation...") ;

  min_curv = mris->min_curv ;
  max_curv = mris->max_curv ;
  if (-min_curv > max_curv)
    max_curv = -min_curv ;

  if (FZERO(max_curv) && FZERO(min_curv))
    max_curv = 1 ;

  /** Pre-processing to paint MEG/EEG activation/temporal information **/
  if (flags & VAL_FLAG)
  {
    int find_flag = 0 ;
    /* find range of values */
    min_curv = 1000000.0f ;
    max_curv = -min_curv ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->val == 0.0)
        continue ;
      find_flag = 1 ;
      if (v->val > max_curv)
        max_curv = v->val ;
      if (v->val < min_curv)
        min_curv = v->val ;
    }
    if (!find_flag) min_curv = max_curv = 0.0 ;
    fprintf(stderr, "min val = %2.3f, max val = %2.3f\n",
            min_curv, max_curv) ;

    if ((flags & TIME_FLAG) || (flags & STAN_FLAG))
    {
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag || v->val == 0.0)
          continue ;
        v->val = (v->val - (min_curv)) / (max_curv-(min_curv)) ;
      }
    }
  }

  for (k=0;k<mris->nfaces;k++) if (!mris->faces[k].ripflag)
    {
      f = &mris->faces[k];
      marked = 0 ;

      if (flags & MESH_FLAG)
      {
        float v2[3], v3[3] ;

        glBegin(GL_LINES);
        glColor3ub(255,255,255);
        for (n=0;n<VERTICES_PER_FACE ;n++)
        {
          v = &mris->vertices[f->v[n]];
          vn = n == VERTICES_PER_FACE-1 ?
               &mris->vertices[f->v[0]] : &mris->vertices[f->v[n+1]] ;
          load_brain_coords(v->nx,v->ny,v->nz,v1);
          glNormal3fv(v1);
          load_brain_coords(v->x,v->y,v->z,v2);
          load_brain_coords(vn->x, vn->y, vn->z,v3);
          glVertex3fv(v2);
          glVertex3fv(v3);
        }
        glEnd() ;
      }
      for (n=0;n<VERTICES_PER_FACE;n++)
        if (mris->vertices[f->v[n]].marked)
          marked = mris->vertices[f->v[n]].marked ;

      glBegin(GL_TRIANGLES) ;
      for (n=0;n<VERTICES_PER_FACE;n++)
      {
        v = &mris->vertices[f->v[n]];

        if (v->ripflag)
          mv = 0 ;

        if ((flags & PATCH_FLAG) && (v->nz < 0))  /* invert negative vertices */
          v->nz *= -1 ;

        /* don't display negative flat stuff */
        if ((flags & PATCH_FLAG) && (v->nz < 0) && !(flags & NEG_FLAG))
          continue ;

        if (flags & COORD_FLAG && getenv("THIN_LINES") == NULL)
        {
          int    itheta, iphi ;
          float  phi_dist, theta_dist, dist, xc, yc, zc, radius, theta, phi ;

          radius = sqrt(SQR(v->cx)+SQR(v->cy)+SQR(v->cz)) ;
          itheta =
            nint((DEGREES(v->theta) / oglu_coord_spacing)) 
            * oglu_coord_spacing ;
          iphi =
            nint((DEGREES(v->phi)   / oglu_coord_spacing)) 
            * oglu_coord_spacing ;
          phi = RADIANS((double)iphi) ;
          theta = RADIANS((double)itheta) ;
          xc = radius * sin(phi) * cos(v->theta) ;
          yc = radius * sin(phi) * sin(v->theta) ;
          zc = radius * cos(phi) ;
          phi_dist = sqrt(SQR(xc-v->cx)+SQR(yc-v->cy)+SQR(zc-v->cz)) ;
          if (zc >126)
            DiagBreak() ;
          phi = RADIANS((double)iphi) ;
          theta = RADIANS((double)itheta) ;
          xc = radius * sin(v->phi) * cos(theta) ;
          yc = radius * sin(v->phi) * sin(theta) ;
          zc = radius * cos(v->phi) ;
          theta_dist = sqrt(SQR(xc-v->cx)+SQR(yc-v->cy)+SQR(zc-v->cz)) ;

          if (theta_dist < phi_dist)
            dist = theta_dist ;
          else
            dist = phi_dist ;

          if (k == 50322 && n == 2)
            DiagBreak() ;
          if (zc >126 && dist < oglu_coord_thickness/2)
            DiagBreak() ;
#if 0
          if (dist > oglu_coord_thickness)
            coord_coef = 0.0f ;
          else   /* ramp values down from 1 based on distance */
          {
            coord_coef = 1.0f - dist/oglu_coord_thickness ;
          }
#else
          coord_coef = exp(-(SQR(dist)) / (SQR(2*oglu_coord_thickness))) ;
#endif
        }
        else
          coord_coef = 0 ;

        if (marked)
        {
          switch (marked)
          {
          default:
          case MARK_RED:
            glColor3ub(255,0,0) ;
            break ;
            break ;
          case MARK_WHITE:
            glColor3ub(255,255,255) ;
            break ;
            break ;
          case MARK_GREEN:
            glColor3ub(0,255,0) ;
            break ;
            break ;
          case MARK_BLUE:
            glColor3ub(0,0,255) ;
            break ;
          case MARK_CYAN:
            glColor3ub(0,255,255) ;
            break ;
          case MARK_YELLOW:
            glColor3ub(255,255, 0) ;
            break ;
          case MARK_ORANGE:
            glColor3ub(255,128, 128) ;   /* orange */
            break ;
          case MARK_LIGHTGREEN:
            glColor3ub(200,255, 200) ;
            break ;
          case MARK_PURPLE:
            glColor3ub(255,0,255) ;
            break ;
          }
        }
        else if (v->border && !(flags & BW_FLAG) && !(flags & NOBORDER_FLAG))
          glColor3f(240,240,0.0);
        else   /* color it depending on curvature */
        {
#define DARK_GRAY    (brightness - (brightness / 3.0f))
#define BRIGHT_GRAY  (brightness + (brightness / 3.0f))
#define MIN_GRAY     min_gray
#define BRIGHTNESS   brightness

          red = green = blue = MIN_GRAY ;
          if (FZERO(max_curv))  /* no curvature info */
            red = green = blue = GRAY ;    /* paint it gray */

          if (f->v[n] == 62152)
            DiagBreak() ;

          if (flags & VAL_FLAG && (fabs(v->val) >= fthresh))
          {
            color_val = v->val ;
            set_color(v->val, v->curv, flags) ;
          }
          else
          {
            float  abs_color_val  ;

            if (v->curv < 0)
              color_val = -v->curv / min_curv ;
            else
              color_val = v->curv / max_curv ;

            color_val = tanh(cslope * color_val) ;
            abs_color_val = fabs(color_val) ;
            red =
              BRIGHTNESS * (MIN_GRAY * (1.0f - abs_color_val)+MAX(0,color_val));
            green =
              BRIGHTNESS * (MIN_GRAY * (1.0f - abs_color_val)+MAX(0,-color_val));
            blue =
              BRIGHTNESS * (MIN_GRAY * (1.0f - abs_color_val)) ;

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
            if (flags & BW_FLAG)
            {
              if (v->curv > 0)
                red = green = blue = BRIGHT_GRAY ;
              else
                red = green = blue = DARK_GRAY ;

            }
            red += (COORD_RED - red) * coord_coef ;
            blue += (COORD_BLUE - blue) * coord_coef ;
            green += (COORD_GREEN - green) * coord_coef ;
            glColor3ub(red,green,blue);  /* specify the RGB color */
          }
        }
        load_brain_coords(v->nx,v->ny,v->nz,v1);
        glNormal3fv(v1);                /* specify the normal for lighting */
        load_brain_coords(v->x,v->y,v->z,v1);
        glVertex3fv(v1);                /* specify the position of the vertex*/
      }
      glEnd() ;   /* done specifying this polygon */
    }


  /* draw canonical coordinate system */
  if (flags & COORD_FLAG && getenv("THIN_LINES") != NULL)
  {
    float  cheight = 15*0.05 ; /* dist edge above surface (0.05 in surfer.c) */
    FACE   *f ;
    float  nx,ny,nz ;
    float  cross_x[4],cross_y[4],cross_z[4];
    int    nintersections;
    VERTEX *va, *vb, *vf[VERTICES_PER_FACE] ;
    double avg_coord, phi_avg, theta_avg, coorda, coordb, coord_line ;
    double fcoords[VERTICES_PER_FACE], frac ;
    int    cno, fvno ;

    for (k = 0 ; k < mris->nfaces ; k++)
    {
      f = &mris->faces[k];
      if (f->ripflag)
        continue ;
      phi_avg = theta_avg = nx = ny = nz = 0.0f ;

      /* calculate average surface coordinate and average normal */
      for (fvno = 0 ; fvno < VERTICES_PER_FACE ; fvno++)
      {
        if (f->v[fvno] == 1210)
          DiagBreak() ;
        vf[fvno] = &mris->vertices[f->v[fvno]] ;
        phi_avg += vf[fvno]->phi ;
        theta_avg += vf[fvno]->theta ;
        nx += vf[fvno]->nx ;
        ny += vf[fvno]->ny ;
        nz += vf[fvno]->nz ;
      }
      phi_avg /= (float)VERTICES_PER_FACE ;
      theta_avg /= (float)VERTICES_PER_FACE ;
      nx /= (float)VERTICES_PER_FACE ;
      ny /= (float)VERTICES_PER_FACE ;
      nz /= (float)VERTICES_PER_FACE ;

      /* find intersection of face with nearest coord line for each coord */
      for (cno = 0 ; cno < 2 ; cno++)    /* once per coordinate */
      {
        if (!cno)    /* phi */
        {
          avg_coord = DEGREES(phi_avg) ;
          RGBcolor(255,230,0);
          RGBcolor(180,135,255);
          RGBcolor(150,255,255);
          RGBcolor(255,230, 0) ;
          RGBcolor(0,0, 255) ;      /* blue */
          RGBcolor(255,255, 0) ;
          glLineWidth(oglu_coord_thickness) ;
          for (fvno = 0 ; fvno < VERTICES_PER_FACE ; fvno++)
            fcoords[fvno] = DEGREES(vf[fvno]->phi) ;
        }
        else         /* theta */
        {
          avg_coord = DEGREES(theta_avg) ;
          RGBcolor(180,135,255);
          RGBcolor(150,255,255);
          RGBcolor(255,230, 0) ;
          RGBcolor(0,0, 255) ;      /* blue */
          RGBcolor(255,255, 0) ;    /* yellow */
          glLineWidth(oglu_coord_thickness) ;
          for (fvno = 0 ; fvno < VERTICES_PER_FACE ; fvno++)
            fcoords[fvno] = DEGREES(vf[fvno]->theta) ;
        }


        coord_line = nint(avg_coord/oglu_coord_spacing) * oglu_coord_spacing;

        /* find out if coordinate line passes between these two vertices */
        nintersections = 0 ;   /* assume no coordinate lines cross this face */
        for (fvno = 0 ; nintersections<25 && fvno < VERTICES_PER_FACE ; fvno++)
        {
          coorda = fcoords[fvno] ;
          coordb = fcoords[fvno < VERTICES_PER_FACE-1 ? fvno+1 : 0] ;
          va = vf[fvno] ;
          vb = vf[fvno < VERTICES_PER_FACE-1 ? fvno+1 : 0] ;

          if (coorda <= coord_line && coordb >= coord_line)
          {
            /* intersects between 0 and 1 */
            frac = (coord_line - coorda) / (coordb - coorda) ;
            cross_x[nintersections] = va->x + frac * (vb->x - va->x) ;
            cross_y[nintersections] = va->y + frac * (vb->y - va->y) ;
            cross_z[nintersections] = va->z + frac * (vb->z - va->z) ;
            nintersections++ ;
          }
          else if (coorda >= coord_line && coordb <= coord_line)
          {
            /* intersects between 1 and 0 */
            frac = (coord_line - coordb) / (coorda - coordb) ;
            cross_x[nintersections] = vb->x + frac * (va->x - vb->x) ;
            cross_y[nintersections] = vb->y + frac * (va->y - vb->y) ;
            cross_z[nintersections] = vb->z + frac * (va->z - vb->z) ;
            nintersections++ ;
          }
          if (nintersections >= 2)   /* draw the line */
          {
            glBegin(GL_LINES) ;

            load_brain_coords(nx,ny,nz,v1);
            glNormal3fv(v1);
            load_brain_coords(cross_x[0]+cheight*nx,cross_y[0]+cheight*ny,
                              cross_z[0]+cheight*nz,v1);
            glVertex3fv(v1);
            load_brain_coords(cross_x[1]+cheight*nx,cross_y[1]+cheight*ny,
                              cross_z[1]+cheight*nz,v1);
            glVertex3fv(v1);
            glEnd() ;
          }
        }
#if 0
        if (nintersections == 1 || nintersections == 3)
        {
          fprintf(stderr, "warning, f %d: %d intersections!\n",
                  k,nintersections);
          DiagBreak() ;
          RGBcolor(255,255,255);
          glBegin(GL_LINES) ;
          load_brain_coords(nx,ny,nz,v1);
          glNormal3fv(v1);
          load_brain_coords(vf[0]->x, vf[0]->y, vf[0]->z, v1);
          glVertex3fv(v1);
          load_brain_coords(vf[2]->x, vf[2]->y, vf[2]->z, v1);
          glVertex3fv(v1);
          glEnd() ;
        }
#endif
      }
    }
  }


#if 0
  if (flags & COORD_FLAG)    /* draw canonical coordinate system */
  {
    int itheta, iphi, vno ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;

      theta = DEGREES(v->theta) ;
      phi = DEGREES(v->phi) ;
      itheta = nint(theta) ;
      iphi = nint(phi) ;
      coord =
        ((((iphi/oglu_coord_spacing)   * oglu_coord_spacing) == iphi) ||
         (((itheta/oglu_coord_spacing) * oglu_coord_spacing) == itheta));
      if (!coord)
        continue ;

      glColor3ub(255,255,255) ;    /* paint the coordinate line white */
      load_brain_coords(v->nx,v->ny,v->nz,v1);
      glNormal3fv(v1);                /* specify the normal for lighting */
      load_brain_coords(v->x,v->y,v->z,v1);
      glVertex3fv(v1);                /* specify the position of the vertex*/
    }
  }
#endif

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
      }
    }


  /** DS Color Scale and writing ... **/


  /** Then draw the color scale (stolen from tksurfer) **/
#if 1
  if (flags & CSCALE_FLAG)
  {
    if ((flags & STAN_FLAG) || (flags & TIME_FLAG))
    {
      draw_colscalebar_time(flags) ;
    }
    else
    {
      draw_colscalebar() ;
    }
  }
#else
  for (i=0;i<100-1;i++)
  {
    stat = i/(100.0-1.0);
    set_color(stat,0.0, flags);
    glBegin(GL_POLYGON);
    v1[0] = FOV*SCALE_FACTOR*(COLSCALEBAR_XPOS-COLSCALEBAR_WIDTH/2);
    v1[1] = FOV*SCALE_FACTOR*(COLSCALEBAR_YPOS+COLSCALEBAR_HEIGHT*(i/(100-1.0)-0.5));
    v1[2] = FOV*SCALE_FACTOR*9.99;
    glVertex3fv(v1);
    v1[0] = FOV*SCALE_FACTOR*(COLSCALEBAR_XPOS+COLSCALEBAR_WIDTH/2);
    glVertex3fv(v1);
    v1[1] = FOV*SCALE_FACTOR*(COLSCALEBAR_YPOS+COLSCALEBAR_HEIGHT*((i+1)/(100-1.0)-0.5));
    glVertex3fv(v1);
    v1[0] = FOV*SCALE_FACTOR*(COLSCALEBAR_XPOS-COLSCALEBAR_WIDTH/2);
    glVertex3fv(v1);
    glEnd() ;
  }
#endif

  if (flags & LEG_FLAG)
  {
    /** All the following mess is just to write 2 numbers using a package
      found on the net DS. **/

    /** ** First load the texture corresponding to a given font, i.e.
      default in this case ** **/
    txf = txfLoadFont("/space/raid/2/users/inverse/build/mris2rgb/default.txf");
    if (txf == NULL)
    {
      fprintf(stderr, "Problem loading %s, %s\n",
              "default.txf", txfErrorString());
      exit(1);
    }

    glLoadIdentity() ;
    v1[0] = v1[1] = 0;
    v1[2] = 1.0 ;
    glNormal3fv(v1) ;

    set_color(0.01,0.0, flags);

    /** ** Texture computation ** **/
    txfEstablishTexture(txf, 0, GL_TRUE);

    glEnable(GL_TEXTURE_2D);
    glAlphaFunc(GL_GEQUAL, 0.0625);
    glEnable(GL_ALPHA_TEST);

    /** ** put here the text you want to render ... ** **/
    if (flags & TIME_FLAG)
    {
      sprintf(text, "%.0fms", min_curv) ;
    }
    else
    {
      sprintf(text, "%.2f", min_curv) ;
    }

    /** ** Text positioning and scaling ** **/
    glMatrixMode(GL_MODELVIEW) ;

    if (flags & STAN_FLAG)
    {
      glTranslatef(100.0, -137, 1.0);
    }
    else
    {
      glTranslatef(90.0, -137, 1.0);
    }
    glScalef(0.2, 0.2, 0.2);

    /** ** and finally render "text" ** **/
    txfRenderString(txf, text, strlen(text));


    /** ** Once again for the second number ** **/
    glLoadIdentity() ;
    glNormal3fv(v1) ;
    glMatrixMode(GL_MODELVIEW);

    if (flags & STAN_FLAG)
    {
      glTranslatef(100.0, -67, 1.0);
    }
    else
    {
      glTranslatef(90.0, -67, 1.0);
    }
    glScalef(0.2, 0.2, 0.2);

    if (flags & TIME_FLAG)
    {
      sprintf(text, "%.0fms", max_curv) ;
    }
    else
    {
      sprintf(text, "%.2f", max_curv) ;
    }
    txfRenderString(txf, text, strlen(text));

    /** Restore geometrical information **/
    glPopMatrix();

    /** End of the text rendering ... **/
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;

  error = glGetError() ;
  if (error != GL_NO_ERROR)
  {
    const char *errstr ;

    errstr = (const char*)gluErrorString(error) ;
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

  mris->xlo = xlo ;
  mris->xhi = xhi ;
  mris->ylo = ylo ;
  mris->yhi = yhi ;
  mris->zlo = zlo ;
  mris->zhi = zhi ;
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
    oglu_fov = max_dim*1.1 ;
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

int
OGLUsetCoordParms(double coord_thickness, int num_coord_lines)
{
  oglu_coord_thickness = coord_thickness ;
  oglu_coord_spacing = 360.0 / (double)num_coord_lines ;
  return(NO_ERROR) ;
}


static int
set_color(float val,float curv, int flags)
{
  short r,g,b;
  float fr,fg,fb,tmpoffset;

  if (curv<0)  tmpoffset = cvfact*offset;
  else         tmpoffset = offset;

  if (flags & VAL_FLAG)
  {
    if (flags & TIME_FLAG)
    {
      set_stat_color_time(val,&fr,&fg,&fb,tmpoffset);
    }
    else
    {
      set_stat_color(val,&fr,&fg,&fb,tmpoffset);
    }
    r=fr;
    g=fg;
    b=fb;
    r = (r<0)?0:(r>255)?255:r;
    g = (g<0)?0:(g>255)?255:g;
    b = (b<0)?0:(b>255)?255:b;
    RGBcolor(r,g,b);
  }
  return(NO_ERROR) ;
}

/** Table color used to paint activation (linest, sptf) **/
static int
set_stat_color(float f, float *rp, float *gp, float *bp, float tmpoffset)
{
  float r,g,b;
  float ftmp,c1,c2;

  if (fabs(f)>fthresh && fabs(f)<fmid)
  {
    ftmp = fabs(f);
    c1 = 1.0/(fmid-fthresh);
    if (fcurv!=1.0)
      c2 = (fmid-fthresh-fcurv*c1*SQR(fmid-fthresh))/
           ((1-fcurv)*(fmid-fthresh));
    else
      c2 = 0;
    ftmp = fcurv*c1*SQR(ftmp-fthresh)+c2*(1-fcurv)*(ftmp-fthresh)+fthresh;
    f = (f<0)?-ftmp:ftmp;
  }

  if (f>=0)
  {
    r = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
        ((f<fthresh)?0:(f<fmid)?(f-fthresh)/(fmid-fthresh):1);
    g = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
        ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
    b = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0);
  }
  else
  {
    f = -f;
    b = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
        ((f<fthresh)?0:(f<fmid)?(f-fthresh)/(fmid-fthresh):1);
    g = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
        ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
    r = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0);
  }
  r = r*255;
  g = g*255;
  b = b*255;

  *rp = r;
  *gp = g;
  *bp = b;
  return(NO_ERROR) ;
}

static int
read_environment_variables(void)
{
  char *cp ;

  cp = getenv("FCURV") ;
  if (cp)
    fcurv = atof(cp) ;

  cp = getenv("MIN_GRAY") ;
  if (cp)
    min_gray = atoi(cp) ;

  cp = getenv("BRIGHTNESS") ;
  if (cp)
    brightness = atoi(cp) ;

  return(NO_ERROR) ;
}

/****************************************************************************/


/** Color table used to paint temporal map **/
static int
set_stat_color_time(float f, float *rp, float *gp, float *bp, float tmpoffset)
{
  float r,g,b;
#if 0
  float fr,fg,fb;
  float ftmp,c1,c2;
  double blufact = 1.0;
#endif

  if (f>0.0)
  {
    f = 1.0-f ;
    b = tmpoffset/5.0 + ((f<0.25)?f:((f<0.50)?(0.25)*(1-(f-0.25)/(1-0.25)):0));
    g = tmpoffset/5.0 + ((f<0.25)?0:((f<0.50)?2*(f-0.25):(f<time_fthresh)?2*(0.50-0.25)*(1-(f-0.50)/(1-0.50)):0));
    r = ((f<0.50)?0:(f<time_fthresh)?(f-0.50)/(time_fthresh-0.50):10.0);
  }
  else
  {
    r = g = b = tmpoffset ;
  }

  r = r*255;
  g = g*255;
  b = b*255;

  *rp = r;
  *gp = g;
  *bp = b;
  return(NO_ERROR) ;
}

static void
draw_colscalebar(void)
{
  int i;
  float v[3], stat, maxval, v1[3] ;
  int NSEGMENTS = 100;

  /** Save geometry / set starting point ... **/
  maxval = fmid+1.0/fslope;
  glPushMatrix();
  glLoadIdentity() ;
  v1[0] = v1[1] = 0;
  v1[2] = 1.0 ;
  glNormal3fv(v1) ;

  for (i=0;i<NSEGMENTS-1;i++)
  {
    /*
        stat = fthresh+i*(maxval-fthresh)/(NSEGMENTS-1.0);
    */
    stat = -maxval+i*(2*maxval)/(NSEGMENTS-1.0);
    set_color(stat,0.0,VAL_FLAG /*REAL_VAL*/);
    glBegin(GL_POLYGON);
    v[0] = fov*sf*(colscalebar_xpos-colscalebar_width/2);
    v[1] = fov*sf*(colscalebar_ypos+colscalebar_height*(i/(NSEGMENTS-1.0)-0.5));
    v[2] = fov*sf*9.99;
    glVertex3fv(v);
    v[0] = fov*sf*(colscalebar_xpos+colscalebar_width/2);
    glVertex3fv(v);
    v[1] =
      fov*sf*(colscalebar_ypos+colscalebar_height*((i+1)/(NSEGMENTS-1.0)-0.5));
    glVertex3fv(v);
    v[0] = fov*sf*(colscalebar_xpos-colscalebar_width/2);
    glVertex3fv(v);
    glEnd();
  }
  glPopMatrix();
}


static void
draw_colscalebar_time(int flags)
{
  int i;
  float v[3], stat, maxval, v1[3] ;
  int NSEGMENTS = 100;

  /** Save geometry / set starting point ... **/
  maxval = fmid+1.0/fslope;
  glPushMatrix();
  glLoadIdentity() ;
  v1[0] = v1[1] = 0;
  v1[2] = 1.0 ;
  glNormal3fv(v1) ;

  for (i=0;i<NSEGMENTS-1;i++)
  {
    stat = fthresh+i*(maxval-fthresh)/(NSEGMENTS-1.0);
    set_color(stat, 0.0, flags);
    glBegin(GL_POLYGON);
    v[0] = fov*sf*(colscalebar_xpos-colscalebar_width/2);
    v[1] = fov*sf*(colscalebar_ypos+colscalebar_height*(i/(NSEGMENTS-1.0)-0.5));
    v[2] = fov*sf*9.99;
    glVertex3fv(v);
    v[0] = fov*sf*(colscalebar_xpos+colscalebar_width/2);
    glVertex3fv(v);
    v[1] =
      fov*sf*(colscalebar_ypos+colscalebar_height*((i+1)/(NSEGMENTS-1.0)-0.5));
    glVertex3fv(v);
    v[0] = fov*sf*(colscalebar_xpos-colscalebar_width/2);
    glVertex3fv(v);
    glEnd();
  }
  glPopMatrix();
}

#endif // HAVE_OPENGL





