
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
static char vcid[] = "$Id: oglutil.c,v 1.18 1998/09/22 02:45:38 fischl Exp $";
#endif

/*-------------------------------- CONSTANTS -----------------------------*/

#define RGBcolor(R,G,B)  glColor3ub((GLubyte)(R),(GLubyte)(G),(GLubyte)(B))
#define GRAY               100

#define FIELDSIGN_POS      4    /* blue */
#define FIELDSIGN_NEG      5    /* yellow */
#define BORDER             6 
#define MARKED             7 

#define FOV                (256.0f*SCALE_FACTOR)


/*-------------------------------- PROTOTYPES ----------------------------*/

static int  read_environment_variables(void) ;
static int set_color(float val,float curv, int flags) ;
static int set_stat_color(float f, float *rp, float *gp, float *bp, 
                          float tmpoffset) ;
static void load_brain_coords(float x,float y, float z, float v[]) ;
static int mrisFindMaxExtents(MRI_SURFACE *mris) ;
static int ogluSetFOV(MRI_SURFACE *mris, double fov) ;

double oglu_fov = FOV ;
static double oglu_coord_thickness = 1.0 ;   /* sigma of coord line */
static double oglu_coord_spacing = 18.0 ;    /* spacing between coord lines */

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
#if 0
  static GLfloat light2_position[] = { 0.6, 0.6, 1.6, 0.0 };
#else
  static GLfloat light2_position[] = { 0.6, 0.6, 0.6, 0.0 };
#endif
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

#define COORD_RED   255
#define COORD_BLUE  255
#define COORD_GREEN 255

static float min_gray = 0.2f ;
static float brightness = 255.0f ;

int
OGLUcompile(MRI_SURFACE *mris, int *marked_vertices, int flags, float cslope)
{
  int          k, n, red, green, blue, error, mv, marked, vno ;
  face_type    *f;
  VERTEX       *v, *vn ;
  float        v1[3], min_curv, max_curv, coord_coef = 0.0, color_val;
#if 0
  float        phi, theta ;
#endif

  read_environment_variables() ;
/*  ogluSetFOV(mris) ;*/
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "compiling surface tessellation...") ;

  min_curv = mris->min_curv ;
  max_curv = mris->max_curv ;
  if (-min_curv > max_curv)
    max_curv = -min_curv ;

  if (flags & VAL_FLAG)
  {
    /* find range of values */
    min_curv = 1000.0f ; max_curv = -min_curv ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (v->val > max_curv)
        max_curv = v->val ;
      if (v->val < min_curv)
        min_curv = v->val ;
    }
    fprintf(stderr, "min val = %2.3f, max val = %2.3f\n",
            min_curv, max_curv) ;
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
    for (n=0;n<4;n++)
      if (mris->vertices[f->v[n]].marked)
        marked = mris->vertices[f->v[n]].marked ;

    glBegin(GL_QUADS) ;
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
#if 1
      if (flags & COORD_FLAG)
      {
        int    itheta, iphi ;
        float  phi_dist, theta_dist, dist, xc, yc, zc, radius, theta, phi ;

        radius = sqrt(SQR(v->cx)+SQR(v->cy)+SQR(v->cz)) ;
        itheta = 
          nint((DEGREES(v->theta) / oglu_coord_spacing)) * oglu_coord_spacing ;
        iphi =   
          nint((DEGREES(v->phi)   / oglu_coord_spacing)) * oglu_coord_spacing ;
        phi = RADIANS((double)iphi) ; theta = RADIANS((double)itheta) ;
        xc = radius * sin(phi) * cos(v->theta) ;
        yc = radius * sin(phi) * sin(v->theta) ;
        zc = radius * cos(phi) ;
        phi_dist = sqrt(SQR(xc-v->cx)+SQR(yc-v->cy)+SQR(zc-v->cz)) ;
        if (zc >126)
          DiagBreak() ;
        phi = RADIANS((double)iphi) ; theta = RADIANS((double)itheta) ;
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
#endif
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
#define DARK_GRAY    (brightness - (brightness / 4.0f))
#define BRIGHT_GRAY  (brightness + (brightness / 4.0f))
#define MIN_GRAY     min_gray
#define BRIGHTNESS   brightness

        red = green = blue = MIN_GRAY ;
        if (FZERO(max_curv))  /* no curvature info */
          red = green = blue = GRAY ;    /* paint it gray */

        if (flags & VAL_FLAG)
        {
          color_val = v->val ;
          set_color(v->val, 0.0f, flags) ;
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
              red = green = blue = DARK_GRAY ;
            else
              red = green = blue = BRIGHT_GRAY ;
            
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


#if 0
  if (flags & COORD_FLAG)    /* draw canonical coordinate system */
  {
    float  cheight = 15*0.05 ; /* dist edge above surface (0.05 in surfer.c) */
    FACE   *f ;
    float  nx,ny,nz;
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
        nx += vf[fvno]->nx ; ny += vf[fvno]->ny ; nz += vf[fvno]->nz ;
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
          glLineWidth(1.0f);
          for (fvno = 0 ; fvno < VERTICES_PER_FACE ; fvno++)
            fcoords[fvno] = DEGREES(vf[fvno]->phi) ;
        }
        else         /* theta */
        {
          avg_coord = DEGREES(theta_avg) ;
          RGBcolor(180,135,255);
          RGBcolor(150,255,255);
          glLineWidth(1.0f);
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

#endif

#if 0
  if (flags & COORD_FLAG)    /* draw canonical coordinate system */
  {
    int itheta, iphi, vno ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;

      theta = DEGREES(v->theta) ; phi = DEGREES(v->phi) ;
      itheta = nint(theta) ; iphi = nint(phi) ;
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

int
OGLUsetCoordParms(double coord_thickness, double coord_spacing)
{
  oglu_coord_thickness = coord_thickness ;
  oglu_coord_spacing = coord_spacing ;
  return(NO_ERROR) ;
}

static double cvfact = 1.5;
/*static double fadef = 0.7;*/
static double fthresh = 0.0f ;
static double fslope = 30.0f ;
static double offset = 0.25 ;
static double fcurv = 0.0f ;
static double fmid = 0.01f ;

static int
set_color(float val,float curv, int flags)
{
  short r,g,b;
  float fr,fg,fb,tmpoffset;

  if (curv<0)  tmpoffset = cvfact*offset;
  else         tmpoffset = offset;

  if (flags & VAL_FLAG)
  {
    set_stat_color(val,&fr,&fg,&fb,tmpoffset);
    r=fr; g=fg; b=fb;
    r = (r<0)?0:(r>255)?255:r;
    g = (g<0)?0:(g>255)?255:g;
    b = (b<0)?0:(b>255)?255:b;
    RGBcolor(r,g,b);
  }

#if 0
  if (mode==GREEN_RED_CURV)
  {
    f = tanh(cslope*(curv-cmid));
    if (f>0) {
      r = 255 * (offset/blufact + 0.95*(1-offset/blufact)*fabs(f));
      g = 255 * (offset/blufact*(1 - fabs(f)));
    }
    else {
      r = 255 * (offset/blufact*(1 - fabs(f)));
      g = 255 * (offset/blufact + 0.95*(1-offset/blufact)*fabs(f));
    }
    b = 255 * (offset*blufact*(1 - fabs(f)));
  }

  if (mode==REAL_VAL)   /* single val positive or signed */
  {
    if (colscale==HEAT_SCALE)  /* stat */
    {
      set_stat_color(val,&fr,&fg,&fb,tmpoffset);
      r=fr; g=fg; b=fb;
    }
    else  /* positive */
    if (colscale==CYAN_TO_RED || colscale==BLU_GRE_RED || colscale==JUST_GRAY)
    {
      if (val<fthresh) 
      {
        r = g = 255 * (tmpoffset/blufact);
        b =     255 * (tmpoffset*blufact);
      }
      else 
      {
        if (fslope!=0)
          f = (tanh(fslope*fmid)+tanh(fslope*(val-fmid)))/(2-tanh(fslope*fmid));
        else
          f = (val<0)?0:((val>1)?1:val);
        set_positive_color(f,&fr,&fg,&fb,tmpoffset);
        r=fr; g=fg; b=fb;
      }
    }
    else /* signed */
    {
      if (fabs(val)<fthresh) 
      {
        r = g = 255 * (tmpoffset/blufact);
        b =     255 * (tmpoffset*blufact);
      }
      else 
      {
        if (fslope!=0)
        {
          if (fmid==0)
            f = tanh(fslope*(val));
          else
          {
            if (val<0)
              f = -(tanh(fslope*fmid) + tanh(fslope*(-val-fmid)))/
                   (2-tanh(fslope*fmid));
            else
              f = (tanh(fslope*fmid) + tanh(fslope*( val-fmid)))/
                  (2-tanh(fslope*fmid));
          }
        }
        else
          f = (val<-1)?-1:((val>1)?1:val);
        if (revphaseflag)
          f = -f;
        set_signed_color(f,&fr,&fg,&fb,tmpoffset);
        r=fr; g=fg; b=fb;
      }
    }
  }

  if (mode==FIELDSIGN_POS || mode==FIELDSIGN_NEG) {
    if (val<fthresh) {
      r = g = 255 * (tmpoffset/blufact);
      b =     255 * (tmpoffset*blufact);
    }
    else {
      f = (1.0 + tanh(fslope*(val-fmid)))/2.0;
      if (mode==FIELDSIGN_POS) {
        b = 255 * (tmpoffset + 0.95*(1-tmpoffset)*fabs(f));
        r = g = 255* (tmpoffset*(1 - fabs(f)));
      }
      else {
        b = 255 * (tmpoffset*(1 - fabs(f)));
        r = g = 255 * (tmpoffset + 0.95*(1-tmpoffset)*fabs(f));
      }
    }
  }

  if (mode==BORDER)  /* AMD 5/27/95 */
  {
    r = 255;
    g = 255;
    b = 0;
  }

  if (mode==MARKED)
  {
    r = 255;
    g = 255;
    b = 255;
  }
  r = (r<0)?0:(r>255)?255:r;
  g = (g<0)?0:(g>255)?255:g;
  b = (b<0)?0:(b>255)?255:b;
  RGBcolor(r,g,b);
#endif
  return(NO_ERROR) ;
}


static int
set_stat_color(float f, float *rp, float *gp, float *bp, float tmpoffset)
{
  float r,g,b;
  float ftmp,c1,c2;

#if 0
  if (invphaseflag)
     f = -f;
  if (truncphaseflag && f<0)
    f = 0;
  if (rectphaseflag)
     f = fabs(f);
#endif

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
  } else
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

#if 0
  if (colscale==HEAT_SCALE)
  {
    if (f>=0)
    {
      r = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fthresh)?0:(f<fmid)?(f-fthresh)/(fmid-fthresh):1);
      g = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
      b = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0);
    } else
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
  }
  else if (colscale==BLU_GRE_RED)
  {
/*
    if (f<0) f = -f;
    b = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
        ((f<fthresh)?0:(f<fmid)?(f-fthresh)/(fmid-fthresh):
         (f<fmid+0.25/fslope)?1-4*(f-fmid)*fslope:
         (f<fmid+0.75/fslope)?0:
         (f<fmid+1.00/fslope)?4*(f-(fmid+0.75/fslope))*fslope:1);
    g = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
        ((f<fmid)?0:(f<fmid+0.25/fslope)?4*(f-fmid)*fslope:
         (f<fmid+0.50/fslope)?1-4*(f-(fmid+0.25/fslope))*fslope:
         (f<fmid+0.75/fslope)?4*(f-(fmid+0.50/fslope))*fslope:1);
    r = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
         ((f<fmid+0.25/fslope)?0:(f<fmid+0.50/fslope)?4*(f-(fmid+0.25/fslope))
                                                                    *fslope:1);
*/
    if (f>=0)
    {
      r = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fthresh)?0:(f<fmid)?(f-fthresh)/(fmid-fthresh):1);
      g = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
      b = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
    } else
    {
      f = -f;
      b = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fthresh)?0:(f<fmid)?(f-fthresh)/(fmid-fthresh):1);
      g = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
      r = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
    }
    r = r*255;
    g = g*255;
    b = b*255;
  }
  else if (colscale==JUST_GRAY)
  {
    if (f<0) f = -f;
    r = g = b = f*255;
  }
#endif
  *rp = r;
  *gp = g;
  *bp = b;
  return(NO_ERROR) ;
}
#if 0
set_complexval_color(x,y,stat,curv)
float x,y,stat,curv;
{
  short sr,sg,sb;
  float f,a,r,g,b;
  float tmpoffset,fscale;
  float a_cycles, oa;

  if (statflag)
    f = stat;
  else
    f = sqrt(x*x+y*y);

  if (curv<0.0) tmpoffset = cvfact*offset;
  else          tmpoffset = offset;

  if (fabs(f)<fthresh)  /* trunc */
  {
    r = g = 255 * (tmpoffset/blufact);
    b =     255 * (tmpoffset*blufact);
  }
  else  /* use complex (or stat vals which ignore complex!!) */
  {
    if (!statflag)
    {
      if (fslope!=0)
        f = (1.0 + tanh(fslope*(f-fmid)))/2.0;
      else
        f = (f<0)?0:((f>1)?1:f);

      if (truncphaseflag) 
      {
        a = atan2(y,x)/(2*M_PI);
        if (revphaseflag)
          a = -a;
        if (invphaseflag)
          a += 0.5;
        a -= angle_offset;
        a = fmod(a,1.0);
        if (a<0) a += 1;
        if (a>0.5)
          f = 0;
      }
    }

    fscale = f;

    if (colscale==HEAT_SCALE || colscale==CYAN_TO_RED ||
        colscale==BLU_GRE_RED || colscale==JUST_GRAY)
    {
      if (statflag)
        set_stat_color(f,&r,&g,&b,tmpoffset);
      else
        set_positive_color(f,&r,&g,&b,tmpoffset);
    }

    else if (colscale==TWOCOND_GREEN_RED)
    {
      a = atan2(y,x)/(2*M_PI);
      if (revphaseflag)
        a = -a;
      if (invphaseflag)
        a += 0.5;
      a -= angle_offset;       /* pos num cancels delay */
      a = fmod(a,1.0);         /* -1 to 1 */
      if (a<0) a += 1;         /* make positive */
      r = g = b = 0;
      f = sin(a*2*M_PI);
      if (f>0.0)
        r = 1;
      else
        g = 1;
      f = fabs(f)*fscale;
      r = 255 * (tmpoffset/blufact*(1-f)+f*r);
      g = 255 * (tmpoffset/blufact*(1-f)+f*g);
      b = 255 * (tmpoffset*blufact*(1-f));
    }

    else if (colscale==COLOR_WHEEL)
    {
      a = atan2(y,x)/(2*M_PI);
      if (revphaseflag)
        a = -a;
      if (invphaseflag)
        a += 0.5;
      a_cycles = angle_cycles;
      oa = a;

      if (fmod(angle_cycles,1.0)==0.0)  /* integral cycles (eccentricity) */
      {
        a -= angle_offset;
        a = fmod(a,1.0);
        if (a<0) a += 1;
        a -= 0.333333;           /* center on blue (1/3)*/
        a = a_cycles*a;          /* allow multiple */
        a = fmod(a,1.0);
        if (a<0) a += 1;
        a = 3*a;
        r = g = b = 0;
        if (a<1.0)
        {
          r = 1-a;
          b = a;
        }
        else if (a<2.0)
        {
          b = 1-(a-1);
          g = a-1;
        }
        else
        {
          r = a-2;
          g = 1-(a-2);
        }
      }
      else /* non-integral cycles (polar angle) */
      {
        a -= angle_offset;
        a = fmod(a,1.0);
        if (a<0) a += 1;
        a -= 0.5;          /* center on blue (1/2) */
        a = a_cycles*a;
        r = g = b = 0;
        if (a<-0.33333)
        {
          r = 1;
        }
        else if (a<0.0)
        {
          r = 1-(a-(-0.33333))/(0.0-(-0.33333));
          b = (a-(-0.33333))/(0.0-(-0.33333));
        }
        else if (a<0.33333)
        {
          b = 1-(a)/(0.33333);
          g = (a)/(0.33333);
        }
        else
        {
          g = 1;
        }

        if (a>fadef*a_cycles/2)
        {
           f = 1-(a-fadef*a_cycles/2)/(a_cycles/2-fadef*a_cycles/2);
           r = (tmpoffset*(1-f)+f*r);
           g = (tmpoffset*(1-f)+f*g);
           b = (tmpoffset*(1-f)+f*b);
        }
        if (a<-fadef*a_cycles/2)
        {
           f = (a-(-a_cycles/2))/(a_cycles/2-fadef*a_cycles/2);
           r = (tmpoffset*(1-f)+f*r);
           g = (tmpoffset*(1-f)+f*g);
           b = (tmpoffset*(1-f)+f*b);
        }

      } /* end non-integral */
      r = (tmpoffset*(1-fscale)+fscale*r)*255;
      b = (tmpoffset*(1-fscale)+fscale*b)*255;
      g = (tmpoffset*(1-fscale)+fscale*g)*255;
    }  /* end color wheel */

    else if (colscale==RYGB_WHEEL)
    {
      a = atan2(y,x)/(2*M_PI);
      if (revphaseflag)
        a = -a;
      if (invphaseflag)
        a += 0.5;
      a_cycles = angle_cycles;
      oa = a;

      a -= angle_offset;
      a = fmod(a,1.0);
      if (a<0) a += 1;
      a -= 0.25;           /* center on blue (1/4)*/
      a = a_cycles*a;          /* allow multiple */
      a = fmod(a,1.0);
      if (a<0) a += 1;
      a = 4*a;
      r = g = b = 0;
      if (a<1.0)
      {
        r = 1.0;
        g = a;
      }
      else if (a<2.0)
      {
        r = 1-(a-1);
        g = 1.0;
      }
      else if (a<3.0)
      {
        g = 1-(a-2);
        b = a-2;
      }
      else
      {
        r = a-3;
        b = 1-(a-3);
      }
      r = (tmpoffset*(1-fscale)+fscale*r)*255;
      b = (tmpoffset*(1-fscale)+fscale*b)*255;
      g = (tmpoffset*(1-fscale)+fscale*g)*255;
    }  /* end RYGB wheel */

    if (phasecontourflag) {
      if (phasecontour_min < phasecontour_max) {
        if (oa>phasecontour_min&&oa<phasecontour_max) {
          if (phasecontourmodflag)
            r = g = b = (tmpoffset*(1-fscale)+fscale*1.0)*255;
          else
            r = g = b = phasecontour_bright;
        }
      }
      else { /* wrap */
        if (oa>phasecontour_min||oa<phasecontour_max) {
          if (phasecontourmodflag)
            r = g = b = (tmpoffset*(1-fscale)+fscale*1.0)*255;
          else
            r = g = b = phasecontour_bright;
        }
      }
    }
  }
  sr = (r<0)?0:(r>255)?255:r;
  sg = (g<0)?0:(g>255)?255:g;
  sb = (b<0)?0:(b>255)?255:b;
  RGBcolor(sr,sg,sb);
}
#endif
static int
read_environment_variables(void)
{
  char *cp ;

  cp = getenv("FTHRESH") ;
  if (cp)
    fthresh = atof(cp) ;

  cp = getenv("FSLOPE") ;
  if (cp)
    fslope = atof(cp) ;

  cp = getenv("FCURV") ;
  if (cp)
    fcurv = atof(cp) ;

  cp = getenv("FMID") ;
  if (cp)
    fmid = atof(cp) ;

  cp = getenv("MIN_GRAY") ;
  if (cp)
    min_gray = atoi(cp) ;

  cp = getenv("BRIGHTNESS") ;
  if (cp)
    brightness = atoi(cp) ;

  return(NO_ERROR) ;
}

