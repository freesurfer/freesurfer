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

static char vcid[] = "$Id: mris_show.c,v 1.18 1997/12/19 20:28:53 fischl Exp $";


/*-------------------------------- CONSTANTS -----------------------------*/

#define BIG_FOV                 300   /* for unfolded hemispheres */
#define SMALL_FOV               200

#define FRAME_SIZE         600

#define DELTA_ANGLE  18.0f
#define RIGHT_HEMISPHERE_ANGLE   90.0
#define LEFT_HEMISPHERE_ANGLE   -90.0

#define ORIG_SURFACE_LIST   1
#define ELLIPSOID_LIST      1

#ifndef GLUT_ESCAPE_KEY
#define GLUT_ESCAPE_KEY  27
#endif

#define MAX_MARKED       5000

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
/* static void rotate_brain(float angle, char c) ;*/

static void keyboard_handler(unsigned char key, int x, int y) ;
static void special_key_handler(int key, int x, int y) ;
static void reshape_handler(int width, int height) ;
static void mouse_handler(int button, int state, int x, int y) ;
#if 0
static void controlLights(int value) ;
#endif

static void display_handler(void) ;
static void home(MRI_SURFACE *mris) ;

/*-------------------------------- DATA ----------------------------*/

char *Progname ;
static char *surf_fname ;
static MRI_SURFACE  *mris, *mris_ellipsoid = NULL ;
static int nmarked = 0 ;
static int marked_vertices[MAX_MARKED+1] = { -1 } ;
static long frame_xdim = FRAME_SIZE;
static long frame_ydim = FRAME_SIZE;

static int talairach_flag = 0 ;
static int ellipsoid_flag = 0 ;

static float x_angle = 0.0f ;
static float y_angle = 0.0f ;
static float z_angle = 0.0f ;

static float delta_angle = 2*DELTA_ANGLE ;


static int current_list = ORIG_SURFACE_LIST ;

static INTEGRATION_PARMS  parms ;


static void findAreaExtremes(MRI_SURFACE *mris) ;
static char curvature_fname[100] = "" ;


static float cslope = 1.75f ;
static int patch_flag = 0 ;
static int compile_flags = 0 ;
static int mean_curvature_flag = 0 ;
static int gaussian_curvature_flag = 0 ;
static int fit_flag = 0 ;
static int fov = -1 ;
static int noscale = 1 ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[])
{
  char         **av, *in_fname, *out_fname, wname[200], fname[100], hemi[10],
               *cp, path[100], name[100] ;
  int          ac, nargs ;
  float        angle ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.projection = PROJECT_ELLIPSOID ;
  parms.tol = TOL ;
  parms.n_averages = N_AVERAGES ;
  parms.l_angle = L_ANGLE ;
  parms.l_area = L_AREA ;
  parms.niterations = NITERATIONS ;
  parms.write_iterations = WRITE_ITERATIONS ;
  parms.a = parms.b = parms.c = 0.0f ;  /* ellipsoid parameters */

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    usage_exit() ;

  surf_fname = in_fname = argv[1] ;
  out_fname = argv[2] ;

  if (patch_flag)   /* read in orig surface before reading in patch */
  {
    FileNamePath(surf_fname, path) ;
    cp = strrchr(surf_fname, '.') ;
    if (cp)
    {
      strncpy(hemi, cp-2, 2) ;
      hemi[2] = 0 ;
    }
    else
      strcpy(hemi, "lh") ;
    sprintf(fname, "%s/%s.orig", path, hemi) ;
    mris = MRISread(fname) ;
    MRISreadTriangleProperties(mris, fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    FileNameOnly(surf_fname, name) ;
    if (MRISreadPatch(mris, name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                Progname, surf_fname) ;
  }
  else
  {
    mris = MRISread(in_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
  }
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeSecondFundamentalForm(mris) ;
  if (mean_curvature_flag)
    MRISuseMeanCurvature(mris) ;
  if (gaussian_curvature_flag)
    MRISuseGaussianCurvature(mris) ;
  if (curvature_fname[0])
    MRISreadCurvatureFile(mris, curvature_fname) ;
  if (ellipsoid_flag)
  {
    MRISupdateEllipsoidSurface(mris) ;
    findAreaExtremes(mris) ;
  }
  if (talairach_flag)
    MRIStalairachTransform(mris, mris) ;
  MRIScenter(mris, mris) ;
  if (noscale)
    OGLUnoscale() ;
  if (fov > 0)
    OGLUsetFOV(fov) ;
  else   /* try and pick an appropriate one */
  {
    if (patch_flag)
    {
      int   ngood, vno ;
      float pct_vertices ;
      
      for (vno = ngood = 0 ; vno < mris->nvertices ; vno++)
        if (!mris->vertices[vno].ripflag)
          ngood++ ;
      pct_vertices = (float)ngood / (float)mris->nvertices ;
      if (pct_vertices > 0.75f)
        OGLUsetFOV(BIG_FOV) ;
      else if (pct_vertices < 0.4f)
        OGLUsetFOV(SMALL_FOV) ;
    }
  }
#if 0
  fprintf(stderr, "(%2.1f, %2.1f, %2.1f) --> (%2.1f, %2.1f, %2.1f), ctr "
          "(%2.1f, %2.1f, %2.1f)\n",
          mris->xlo, mris->ylo, mris->zlo, mris->xhi, mris->yhi, mris->zhi,
          mris->xctr, mris->yctr, mris->zctr);
#endif

  glutInit(&argc, argv) ;  /* pass -geometry and such to  GLUT */
  glutInitWindowSize(frame_xdim, frame_ydim) ;

/* 
   use double buffering, RBF color stuff, and a depth buffer for 
   hidden surface removal 
   */
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH) ;
  sprintf(wname, "%scortical surface reconstruction '%s'", 
          talairach_flag ? "talairached " : "", in_fname) ;
  glutCreateWindow(wname) ;          /* open the window */
  glutHideWindow() ;
  glutDisplayFunc(display_handler) ; /* specify a drawing callback function */
  OGLUinit(mris, frame_xdim, frame_ydim) ; /* specify lighting and such */
  if (fit_flag)
    oglu_fov = 0.0 ;

  /* now compile the surface tessellation */
  glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
  OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
  glEndList() ;

  glMatrixMode(GL_MODELVIEW);
  if (patch_flag)
  {
    angle = 90.0f ;
    glRotatef(angle, 1.0f, 0.0f, 0.0f) ;
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
  }
  else
  {
    if (mris->hemisphere == RIGHT_HEMISPHERE)
      angle = RIGHT_HEMISPHERE_ANGLE ;
    else
      angle = LEFT_HEMISPHERE_ANGLE ;
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
  }


  /* handle reshaping, keyboard and mouse events */
  glutKeyboardFunc(keyboard_handler) ;
  glutSpecialFunc(special_key_handler) ;
  glutMouseFunc(mouse_handler) ;
  glutReshapeFunc(reshape_handler) ;

  /* create a menu and attach it to the window */
#if 0
  glutCreateMenu(controlLights) ;
  glutAddMenuEntry("toggle light 1", 1) ;
  glutAddMenuEntry("toggle light 2", 2) ;
  glutAddMenuEntry("toggle light 3", 3) ;
  glutAddMenuEntry("toggle light 4", 4) ;
  glutAttachMenu(GLUT_RIGHT_BUTTON) ;
#endif

  glutShowWindow() ;

  glRotatef(x_angle, 1.0f, 0.0f, 0.0f) ;
  glRotatef(y_angle, 0.0f, 1.0f, 0.0f) ;
  glRotatef(z_angle, 0.0f, 0.0f, 1.0f) ;
  glutPostRedisplay() ;

  glutMainLoop() ;               /* enter event handling loop */
  MRISfree(&mris) ;

  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "cslope"))
  {
    sscanf(argv[2], "%f", &cslope) ;
    nargs = 1 ;
    fprintf(stderr, "using color slope compression = %2.4f\n", cslope) ;
  }
  else if (!stricmp(option, "tp"))
    compile_flags |= TP_FLAG ;
  else if (!stricmp(option, "mesh"))
    compile_flags |= MESH_FLAG ;
  else if (!stricmp(option, "noscale"))
    noscale = 1 ;
  else if (!stricmp(option, "scale"))
    noscale = 0 ;
  else if (!stricmp(option, "fov"))
  {
    fov = atoi(argv[2]) ;
    fprintf(stderr, "using fov = %d\n", fov) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
  {
  case 'M':
    mean_curvature_flag = 1 ;
    break ;
  case 'G':
    gaussian_curvature_flag = 1 ;
    break ;
  case 'P':
    patch_flag = 1 ;
    compile_flags |= PATCH_FLAG ;
    nargs = 0 ;
    break ;
  case 'C':
    strcpy(curvature_fname, argv[2]) ;
    nargs = 1 ;
    break ;
  case 'V':
    if (nmarked >= MAX_MARKED)
      ErrorExit(ERROR_NOMEMORY, "%s: too many vertices marked (%d)",
                Progname, nmarked) ;
    if (*argv[2] == '#' || *argv[2] == ':')
    {
      int  skip, vno, vindex ;

      sscanf(argv[2], "#%d", &skip) ;
      for (vindex = vno = 0 ; vindex < MAX_MARKED ; vindex++, vno += skip)
        marked_vertices[vindex] = vno ;
      nmarked = MAX_MARKED ;
    }
    else
      sscanf(argv[2], "%d", &marked_vertices[nmarked++]) ;
    marked_vertices[nmarked] = -1 ;  /* terminate list */
    nargs = 1 ;
    break ;
  case 'X':
    sscanf(argv[2], "%f", &x_angle) ;
    x_angle = RADIANS(x_angle) ;
    nargs = 1 ;
    break ;
  case 'F':
    fit_flag = 1 ;
    break ;
  case 'Y':
    sscanf(argv[2], "%f", &y_angle) ;
    y_angle = RADIANS(y_angle) ;
    nargs = 1 ;
    break ;
  case 'Z':
    sscanf(argv[2], "%f", &z_angle) ;
    z_angle = RADIANS(z_angle) ;
    nargs = 1 ;
    break ;
  case 'E':
  case 'S':
    ellipsoid_flag = 1 ;
    break ;
  case 'T':
    talairach_flag = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'N':
    sscanf(argv[2], "%d", &parms.niterations) ;
    nargs = 1 ;
    fprintf(stderr, "using niterations = %d\n", parms.niterations) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <input image file>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program will display an MRI surface reconstruction.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


#define SHOW_AXES  0

#define COMPILE_SURFACE 1
static void
display_handler(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ;
#if COMPILE_SURFACE
  glCallList(current_list) ;      /* render sphere display list */
#else
  OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
#endif

#if SHOW_AXES
  glCallList(2) ;      /* show axes */
#endif
  glutSwapBuffers() ;  /* swap rendered buffer into display */
}

static void
keyboard_handler(unsigned char key, int x, int y)
{
  int   redraw = 1 ;
  char  wname[100] ;
  float angle ;

  glMatrixMode(GL_MODELVIEW);
  switch (key)
  {
  case 'W':  /* write it out */
    if (mris->hemisphere == RIGHT_HEMISPHERE)
      sprintf(wname, "rh.surf") ;
    else
      sprintf(wname, "lh.surf") ;
    fprintf(stderr, "writing surface to %s\n", wname) ;
    MRISwrite(mris, wname) ;
    break ;
  case '7':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(-10.0f, 10.0f, 10.0f) ;
    break ;
  case '9':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(10.0f, 10.0f, 0.0f) ;
    break ;
  case '1':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(-10.0f, -10.0f, 0.0f) ;
    break ;
  case '3':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(10.0f, -10.0f, 0.0f) ;
    break ;
  case '8':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(0.0f, 10.0f, 0.0f) ;
    break ;
  case '2':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(0.0f, -10.0f, 0.0f) ;
    break ;
  case '4':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(-10.0f, 0.0f, 0.0f) ;
    break ;
  case '6':
    glMatrixMode(GL_PROJECTION);
    glTranslatef(10.0f, 0.0f, 0.0f) ;
    break ;
  case '+':
    delta_angle *= 2.0f ;
    redraw = 0 ;
    break ;
  case '/':
    cslope /= 1.5f ;
    fprintf(stderr, "cslope = %2.3f\n", cslope) ;
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
    break ;
  case '*':
    cslope *= 1.5f ;
    fprintf(stderr, "cslope = %2.3f\n", cslope) ;
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
    break ;
  case '-':
    delta_angle *= 0.5f ;
    redraw = 0 ;
    break ;
  case GLUT_ESCAPE_KEY:
  case 'q':
    exit(0) ;
    break ;
  case 'h':
    home(mris) ;
    break ;
#if 0
  case 'x':   /* local x axis is reversed */
    glRotatef(-delta_angle, 1.0f, 0.0f, 0.0f) ;
    break ;
  case 'X':
    glRotatef(delta_angle, 1.0f, 0.0f, 0.0f) ;
    break ;
  case 'y':   /* local y is actually z coordinate */
    glRotatef(delta_angle, 0.0f, 0.0f, 1.0f) ;
    break ;
  case 'Y':
    glRotatef(-delta_angle, 0.0f, 0.0f, 1.0f) ;
    break ;
  case 'z':   /* local z is actually y axis */
    glRotatef(delta_angle, 0.0f, 1.0f, 0.0f) ;
    break ;
  case 'Z':
    glRotatef(-delta_angle, 0.0f, 1.0f, 0.0f) ;
    break ;
#endif
  case 'D':
  case 'd':
  case 'r':
    redraw = 1 ;
    break ;
  case 'T':
  case 't':
    sprintf(wname, "talairached cortical surface reconstruction '%s'", 
            surf_fname);
    MRIStalairachTransform(mris, mris) ;
    MRIScenter(mris, mris) ;
#if COMPILE_SURFACE
    glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
#endif
    glutSetWindowTitle(wname) ;
    break ;
  case 'P':
  case 'p':
    mris_ellipsoid = MRISprojectOntoEllipsoid(mris, NULL, 0.0f, 0.0f, 0.0f);
    MRISfree(&mris) ;
    mris = mris_ellipsoid ;
#if COMPILE_SURFACE
    glDeleteLists(ORIG_SURFACE_LIST, 1) ;
    glNewList(ELLIPSOID_LIST, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
#endif
    current_list = ELLIPSOID_LIST ;
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    if (mris_ellipsoid->hemisphere == RIGHT_HEMISPHERE)
      angle = RIGHT_HEMISPHERE_ANGLE ;
    else
      angle = LEFT_HEMISPHERE_ANGLE ;
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
    break ;
#if 0
  case 'u':
    niter = parms.niterations ;
  case 'U':
    MRISunfold(mris, &parms) ;

#if COMPILE_SURFACE
    glDeleteLists(ELLIPSOID_LIST, 1) ;
    glNewList(ELLIPSOID_LIST, GL_COMPILE) ;
    OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
    glEndList() ;
#endif
    break ;
#endif
  default:
    fprintf(stderr, "unknown normal key=%d\n", key) ;
    redraw = 0 ;
    break ;
  }
  if (redraw)
    glutPostRedisplay() ;
}

#if 0
static void
rotate_brain(float angle, char c)
{
  int i,j,k;
  float m[4][4], m1[4][4], m2[4][4]; /* Matrix m,m1,m2; */
  float sa,ca;

  glGetFloatv(GL_MODELVIEW_MATRIX,(float *)(m)) ;

  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
      m1[i][j] = (i==j)?1.0:0.0;
  }
  angle = angle*PI/180;
  sa = sin(angle);
  ca = cos(angle);
  if (c=='y')
  {
    m1[0][0] = m1[2][2] = ca;
    m1[2][0] = -(m1[0][2] = sa);
  } 
  else if (c=='x') 
  {
    m1[1][1] = m1[2][2] = ca;
    m1[1][2] = -(m1[2][1] = sa);
  } 
  else if (c=='z') 
  {
    m1[0][0] = m1[1][1] = ca;
    m1[1][0] = -(m1[0][1] = sa);
  } 
  else 
  {
    printf("surfer: ### Illegal axis %c\n",c);return;
  }
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++) 
    {
      m2[i][j] = 0;
      for (k=0;k<4;k++)
        m2[i][j] += m[i][k]*m1[k][j];
    }
  }

  glLoadMatrixf((float *)m) ;
}
#endif
static void
reshape_handler(int width, int height)
{
  int new_dim ;

  new_dim = MIN(width, height) ;
  glViewport(0,0,frame_xdim = new_dim, frame_ydim = new_dim);
}

static void
mouse_handler(int button, int state, int x, int y)
{
  int kb_state, redraw = 1, width, height ;

  kb_state = glutGetModifiers() ;
  switch (button)
  {
  case GLUT_LEFT_BUTTON:
    width = glutGet(GLUT_WINDOW_WIDTH) ;
    height = glutGet(GLUT_WINDOW_HEIGHT) ;
    x -= width/2 ; y -= height/2 ;  /* put origin to center of window */
    if (kb_state & GLUT_ACTIVE_CTRL && state == GLUT_UP)/* zoom around click */
    {
      glMatrixMode(GL_PROJECTION) ;
      glTranslatef(-SCALE_FACTOR*(float)x, SCALE_FACTOR*(float)y, 0.0f) ;
      glScalef(2.0f, 2.0f, 2.0f) ;
      /*      glTranslatef(0.0f, SCALE_FACTOR*y, SCALE_FACTOR*x) ;*/
    }
    else
    {
      redraw = 0 ;
      fprintf(stderr, "(%d, %d)      \r", x, y) ;
    }
    break ;
  case GLUT_RIGHT_BUTTON:
    if (kb_state & GLUT_ACTIVE_CTRL && state == GLUT_UP)
      home(mris) ;
    break ;
  case GLUT_MIDDLE_BUTTON:
  default:
    redraw = 0 ;
    break ;
  }
  if (redraw)
    glutPostRedisplay() ;
}

static void
special_key_handler(int key, int x, int y)
{
  int   redraw = 1, kb_state ;
  float angle ;

  kb_state = glutGetModifiers() ;

  glMatrixMode(GL_MODELVIEW);
  if (kb_state & GLUT_ACTIVE_SHIFT)
    angle = 180.0f ;
  else
    angle = delta_angle ;

  switch (key)
  {
  case GLUT_KEY_F1:
    break ;
  case GLUT_KEY_LEFT:
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
    break ;
  case GLUT_KEY_RIGHT:
    glRotatef(-angle, 0.0f, 1.0f, 0.0f) ;
    break ;
  case GLUT_KEY_DOWN:
    glRotatef(-angle, 1.0f, 0.0f, 0.0f) ;
    break ;
  case GLUT_KEY_UP:
    glRotatef(angle, 1.0f, 0.0f, 0.0f) ;
    break ;
#define GLUT_KEY_DELETE (GLUT_KEY_INSERT+1)
  case GLUT_KEY_DELETE:
  case GLUT_KEY_PAGE_UP:
    glRotatef(angle, 0.0f, 0.0f, 1.0f) ;
    break ;
  case GLUT_KEY_PAGE_DOWN:
    glRotatef(-angle, 0.0f, 0.0f, 1.0f) ;
    break ;
  case GLUT_KEY_INSERT:
    break ;
  case GLUT_KEY_HOME:
    break ;
#if 0
  case GLUT_KEY_ESCAPE:
    exit(0) ;
    break ;
#endif
  default:
    fprintf(stderr, "unknown special key=%d\n", key) ;
    redraw = 0 ;
    break ;
  }
  if (redraw)
    glutPostRedisplay() ;
}


#if 0
#define NUM_LIGHTS 4
static void
controlLights(int value)
{
  static int light_enabled[NUM_LIGHTS] = {1, 1, 1, 1} ;

  light_enabled[value] = !light_enabled[value] ;

  switch (value)
  {
  case 1:
    if (light_enabled[value])
      glEnable(GL_LIGHT0);
    else
      glDisable(GL_LIGHT0) ;
    break ;
  case 2:
    if (light_enabled[value])
      glEnable(GL_LIGHT1);
    else
      glDisable(GL_LIGHT1) ;
    break ;
  case 3:
    if (light_enabled[value])
      glEnable(GL_LIGHT2) ;
    else
      glDisable(GL_LIGHT2) ;
    break ;
  case 4:
    if (light_enabled[value])
      glEnable(GL_LIGHT3);
    else
      glDisable(GL_LIGHT3) ;
    break ;
  default:
    break ;
  }

  glutPostRedisplay() ;
}
#endif

static void
home(MRI_SURFACE *mris)
{
  float angle ;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if (patch_flag)
  {
    angle = 90.0f ;
    glRotatef(angle, 1.0f, 0.0f, 0.0f) ;
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
  }
  else
  {
    if (mris->hemisphere == RIGHT_HEMISPHERE)
      angle = RIGHT_HEMISPHERE_ANGLE ;
    else
      angle = LEFT_HEMISPHERE_ANGLE ;
    glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
  }
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-oglu_fov, oglu_fov, -oglu_fov, oglu_fov, 
          -10.0f*oglu_fov, 10.0f*oglu_fov);
}
static void
findAreaExtremes(MRI_SURFACE *mris)
{
  int     fno, tno, fmax, tmax, fmin, tmin ;
  FACE    *face ;
  float   min_area, max_area;

  fmax = fmin = tmax = tmin = 0 ;
  min_area = 10000.0f ; max_area = -10000.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      if (face->area[tno] > max_area)
      {
        max_area = face->area[tno] ;
        fmax = fno ;
        tmax = tno ;
      }
      if (face->area[tno] < min_area)
      {
        min_area = face->area[tno] ;
        fmin = fno ;
        tmin = tno ;
      }
    }
  }
  fprintf(stderr, "min_area = %2.3f at f %d, t %d\n", min_area, fmin, tmin) ;
  fprintf(stderr, "max_area = %2.3f at f %d, t %d\n", max_area, fmax, tmax) ;
}

#if 0
static void
set_color(float val, float curv, int mode)
{
  short r,g,b;
  float f,fr,fg,fb,tmpoffset;

  if (curv<0)  tmpoffset = cvfact*offset;
  else         tmpoffset = offset;

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

  if (mode==GREEN_RED_VAL)   /* single val positive or signed */
  {
    if (colscale==HEAT_SCALE)
    {
      set_stat_color(val,&fr,&fg,&fb,tmpoffset);
      r=fr; g=fg; b=fb;
    }
    else
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
    else
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
  glColor3ub(r,g,b);
}
#endif
