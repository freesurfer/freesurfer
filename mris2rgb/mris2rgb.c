#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>

#define OPENGL_CODE  1
#include <GL/glut.h>
#include <GL/glx.h>

#if defined(Linux)|| defined(__sun__)
#include "GL/osmesa.h"
#endif

#include "rgb_image.h"

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "oglutil.h"
#include "tiff.h"
#include "tiffio.h"
#include "label.h"

static char vcid[] = "$Id: mris2rgb.c,v 1.13 1998/05/20 19:37:53 fischl Exp $";

/*-------------------------------- CONSTANTS -----------------------------*/

#define MAX_MARKED               20000
#define RIGHT_HEMISPHERE_ANGLE   90.0
#define LEFT_HEMISPHERE_ANGLE   -90.0
#define BIG_FOV                 300   /* for unfolded hemispheres */
#define SMALL_FOV               200

#define ORIG_SURFACE_LIST        1

#define FILE_CURVATURE           0
#define MEAN_CURVATURE           1
#define GAUSSIAN_CURVATURE       2
#define AREA_ERRORS              3

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static void save_rgb(char *fname, int width, int height, unsigned short *red, 
                     unsigned short *green, unsigned short *blue) ;

static void save_TIFF(char *fname, int width, int height, unsigned char *rgb);
static void grabPixelsTIFF(unsigned int width, unsigned int height, 
         unsigned char *rgb);

static void grabPixels(unsigned int width, unsigned int height,
                       unsigned short *red, unsigned short *green, 
                       unsigned short *blue) ;

static void clear_pixmaps(MRI_SURFACE *mris) ;
static void xinit(void) ;
static void xend(void) ;

#define MAX_CENTROIDS  100
static int centroid_files = 0 ;
static char *centroid_fnames[MAX_CENTROIDS] ;
static int centroid_colors[MAX_CENTROIDS] ;
static int mark_centroids(MRI_SURFACE *mris, char *centroid_fnames[], 
                          int *centroid_colors, int centroid_files) ;

/*-------------------------------- DATA ----------------------------*/

char *Progname ;
static MRI_SURFACE  *mris ;

static long frame_xdim = 600;
static long frame_ydim = 600;

static int talairach_flag = 0 ;
static int medial_flag = 0 ;
static int lateral_flag = 1 ;

static int normalize_flag = 0 ;
static float angle_offset = 0.0f ;
static float x_angle = 0.0f ;
static float y_angle = 0.0f ;
static float z_angle = 0.0f ;

static float scale = 1.0f ;

static int current_list = ORIG_SURFACE_LIST ;

static int curvature_flag = 0 ;

#if defined(Linux) || defined(__sun__)
static OSMesaContext context ;
static void *buffer = NULL ;

#else
static Display     *display = NULL ;
static XVisualInfo *visual ;
static Pixmap      pixmap ;
static GLXPixmap   glxpixmap ;
static GLXContext  context ;
static int         configuration[] =
{
  GLX_DOUBLEBUFFER, GLX_RGBA, GLX_DEPTH_SIZE, 16, GLX_RED_SIZE, 1, 
  GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, None
} ;
#endif
static char *curvature_fname = NULL ;
static char *coord_fname = NULL ;
static float cslope = 5.0f ;
static int   patch_flag = 0 ;

static int compile_flags = 0 ;

static int nmarked = 0 ;
static int marked_vertices[MAX_MARKED] = { -1 } ;
static int noscale = 1 ;
static int fov = -1 ;
static MRI_SURFACE_PARAMETERIZATION *mrisp = NULL ;
static int tiff_flag = 0;
static int param_no = -1 ;
static int normalize_param = 0 ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[])
{
  char            **av, *in_fname, *out_prefix, out_fname[100], name[100],
    path[100], *cp, hemi[100], fname[100], *surf_fname ;
  int             ac, nargs, size, ino, old_status, i ;
  float           angle = 0.0f ;
  unsigned short  *red=NULL, *green=NULL, *blue=NULL;
  unsigned char   *rgb=NULL;
  
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  
  if (argc < 3)
    usage_exit() ;
  
  out_prefix = argv[argc-1] ;
  
  size = frame_xdim * frame_ydim ;
  if (!tiff_flag)
  {
    red = (unsigned short *)calloc(size, sizeof(unsigned short)) ;
    green = (unsigned short *)calloc(size, sizeof(unsigned short)) ;
    blue = (unsigned short *)calloc(size, sizeof(unsigned short)) ;
  }
  else
    rgb = (unsigned char *)calloc(frame_xdim * frame_ydim * 3, 
                                  sizeof(unsigned char));
  
  xinit() ;             /* open connection with x server */
  
  for (ino = 1 ; ino < argc-1 ; ino++)
  {
    surf_fname = in_fname = argv[ino] ;
    FileNameOnly(surf_fname, name) ;
    if (patch_flag)  /* read the orig surface, then the patch file */
    {
      FileNamePath(in_fname, path) ;
      cp = strchr(name, '.') ;
      if (cp)
      {
        strncpy(hemi, cp-2, 2) ;
        hemi[2] = 0 ;
      }
      else
        strcpy(hemi, "lh") ;

      sprintf(fname, "%s/%s.orig", path, hemi) ;
      mris = MRISfastRead(fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, in_fname) ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "reading patch file %s...\n", in_fname) ;
      if (MRISreadPatch(mris, name) != NO_ERROR)
          ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                    Progname, surf_fname) ;
      
    }
    else   /* just read the surface normally */
    {
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "reading surface file %s...\n", in_fname) ;
      mris = MRISfastRead(surf_fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, in_fname) ;
    }

    if (coord_fname)
    {
      if (MRISreadCanonicalCoordinates(mris, coord_fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read canonical coordinate system",
                  Progname) ;
    }
    else
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

    for (i = 0 ; i < nmarked ; i++)
      mris->vertices[marked_vertices[i]].marked = 1 ;
    if (talairach_flag)
    {
      MRIStalairachTransform(mris, mris) ;
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
      
    }
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

    if (ino == 1)    /* do initialization first time through */
      OGLUinit(mris, frame_xdim, frame_ydim) ;/* specify lighting and such */
    
    mark_centroids(mris, centroid_fnames, centroid_colors, centroid_files) ;
    if (mrisp)
    {
      MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
      if (normalize_param)
      {
        MRISnormalizeFromParameterization(mrisp, mris, param_no) ;
      }
      else
        MRISfromParameterization(mrisp, mris, 0) ;
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    }
    if (curvature_fname)
      MRISreadCurvatureFile(mris, curvature_fname) ;

    switch (curvature_flag)
    {
    case MEAN_CURVATURE:
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseMeanCurvature(mris) ;
      break ;
    case GAUSSIAN_CURVATURE:
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseGaussianCurvature(mris) ;
      break ;
    case AREA_ERRORS:
      MRISstoreCurrentPositions(mris) ;
      old_status = mris->status ;
      mris->status = MRIS_SURFACE ;
      MRISreadVertexPositions(mris, "smoothwm") ;
      MRIScomputeMetricProperties(mris) ;
      MRISstoreMetricProperties(mris) ;
      MRISrestoreOldPositions(mris) ;
      mris->status = old_status ;
      MRIScomputeMetricProperties(mris) ;
      MRISuseAreaErrors(mris) ;
      break ;
    }

    if (normalize_flag)
      MRISnormalizeCurvature(mris) ;

    if (patch_flag)
    {
      angle = 90.0f ;
      glRotatef(angle, 1.0f, 0.0f, 0.0f) ;
      glRotatef(angle, 0.0f, 1.0f, 0.0f) ;
      glRotatef(x_angle, 1.0f, 0.0f, 0.0f) ;
      glRotatef(y_angle, 0.0f, 1.0f, 0.0f) ;
      glRotatef(z_angle, 0.0f, 0.0f, 1.0f) ;
        
      glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
      OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
      glEndList() ;
      glCallList(current_list) ;      /* render sphere display list */
      if (!tiff_flag)
      {
        grabPixels(frame_xdim, frame_ydim, red, green, blue) ;
        sprintf(out_fname, "%s/%s.%s.rgb", out_prefix, "lateral", name) ;
        fprintf(stderr, "writing rgb file %s\n", out_fname) ;
        save_rgb(out_fname, frame_xdim, frame_ydim, red, green, blue) ;
      }
      else
      {
        grabPixelsTIFF(frame_xdim, frame_ydim, rgb) ;
        sprintf(out_fname, "%s/%s.%s.tiff", out_prefix, "lateral", name) ;
        fprintf(stderr, "writing TIFF file %s\n", out_fname) ;
        save_TIFF(out_fname, frame_xdim, frame_ydim, rgb) ;
      }
    }
    else
    {
      if (lateral_flag)
      {
        if (mris->hemisphere == RIGHT_HEMISPHERE)
          angle = RIGHT_HEMISPHERE_ANGLE ;
        else
          angle = LEFT_HEMISPHERE_ANGLE ;
        angle += angle_offset ;
        glRotatef(angle, 0.0f, 1.0f, 0.0f) ;

        glRotatef(x_angle, 1.0f, 0.0f, 0.0f) ;
        glRotatef(y_angle, 0.0f, 1.0f, 0.0f) ;
        glRotatef(z_angle, 0.0f, 0.0f, 1.0f) ;
        
        glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
        OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
        glEndList() ;
        glCallList(current_list) ;      /* render sphere display list */
        if (!tiff_flag)
        {
          grabPixels(frame_xdim, frame_ydim, red, green, blue) ;
          sprintf(out_fname, "%s/%s.%s.rgb", out_prefix, "lateral", name) ;
          fprintf(stderr, "writing rgb file %s\n", out_fname) ;
          save_rgb(out_fname, frame_xdim, frame_ydim, red, green, blue) ;
        }
        else
        {
          grabPixelsTIFF(frame_xdim, frame_ydim, rgb) ;
          sprintf(out_fname, "%s/%s.%s.tiff", out_prefix, "lateral", name) ;
          fprintf(stderr, "writing TIFF file %s\n", out_fname) ;
          save_TIFF(out_fname, frame_xdim, frame_ydim, rgb) ;
        }
      }
      
      if (medial_flag)
      {
        if (lateral_flag)  /* clear old display */
        {
          clear_pixmaps(mris) ;
#ifdef IRIX
          /*          glRotatef(-angle, 0.0f, 1.0f, 0.0f) ;*/
#endif
        }

        if (mris->hemisphere == LEFT_HEMISPHERE)
          angle = RIGHT_HEMISPHERE_ANGLE ;
        else
          angle = LEFT_HEMISPHERE_ANGLE ;
        angle += angle_offset ;
        glRotatef(angle, 0.0f, 1.0f, 0.0f) ;

        
        glRotatef(x_angle, 1.0f, 0.0f, 0.0f) ;
        glRotatef(y_angle, 0.0f, 1.0f, 0.0f) ;
        glRotatef(z_angle, 0.0f, 0.0f, 1.0f) ;
        
        glNewList(ORIG_SURFACE_LIST, GL_COMPILE) ;
        OGLUcompile(mris, marked_vertices, compile_flags, cslope) ;
        glEndList() ;
        glCallList(current_list) ;      /* render sphere display list */
        if (!tiff_flag) 
        {
          grabPixels(frame_xdim, frame_ydim, red, green, blue) ;
          sprintf(out_fname, "%s/%s.%s.rgb", out_prefix, "medial", name) ;
          fprintf(stderr, "writing rgb file %s\n", out_fname) ;
          save_rgb(out_fname, frame_xdim, frame_ydim, red, green, blue) ;
        }
        else
        {
          grabPixelsTIFF(frame_xdim, frame_ydim, rgb) ;
          sprintf(out_fname, "%s/%s.%s.tiff", out_prefix, "medial", name) ;
          fprintf(stderr, "writing TIFF file %s\n", out_fname) ;
          save_TIFF(out_fname, frame_xdim, frame_ydim, rgb) ;
        }
      }
    }
    if (ino+1 < argc-1)
      clear_pixmaps(mris) ;

    MRISfree(&mris) ;
  }
  if (!tiff_flag)
  {
    free(red) ; free(blue) ; free(green) ;
  }
  else
    free(rgb);
      
  xend() ;
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
  else if (!stricmp(option, "param"))
  {
    mrisp = MRISPread(argv[2]) ;
    if (!mrisp)
      ErrorExit(ERROR_NOFILE, "%s: could not read parameterization file %s",
                Progname, argv[2]) ;
    fprintf(stderr, "reading parameterized curvature from %s\n", argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "nparam"))
  {
    char *cp ;
    cp = strchr(argv[2], '#') ;
    if (cp)   /* # explicitly given */
    {
      param_no = atoi(cp+1) ;
      *cp = 0 ;
    }
    else
      param_no = 0 ;
    mrisp = MRISPread(argv[2]) ;
    normalize_param = 1 ;
    if (!mrisp)
      ErrorExit(ERROR_NOFILE, "%s: could not read parameterization file %s",
                Progname, argv[2]) ;
    fprintf(stderr, "reading parameterized curvature from %s\n", argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!strcmp(option, "centroids"))
  {
    centroid_fnames[centroid_files] = argv[2] ; 
    centroid_colors[centroid_files] = atoi(argv[3]) ;
    nargs = 2 ;
    fprintf(stderr, "marking all vertices in file %s with %d\n", 
            centroid_fnames[centroid_files],centroid_colors[centroid_files]);
    centroid_files++ ;
  }
  else if (!stricmp(option, "cslope"))
  {
    sscanf(argv[2], "%f", &cslope) ;
    nargs = 1 ;
    fprintf(stderr, "using color slope compression = %2.4f\n", cslope) ;
  }
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
  else if (!stricmp(option, "mean"))
    curvature_flag = MEAN_CURVATURE ;
  else if (!stricmp(option, "tp"))
    compile_flags |= TP_FLAG ;
  else if (!stricmp(option, "bw"))
    compile_flags |= BW_FLAG ;
  else if (!stricmp(option, "neg"))
    compile_flags |= NEG_FLAG ;
  else if (!stricmp(option, "noborder"))
    compile_flags |= NOBORDER_FLAG ;
  else if (!stricmp(option, "cparms"))
  {
    double coord_thickness, coord_spacing ;

    coord_thickness = atof(argv[2]) ;
    coord_spacing = atof(argv[3]) ;
    fprintf(stderr, 
        "spacing coordinate lines %2.2f degrees apart, with thickness %2.1f\n",
            coord_spacing, coord_thickness) ;
    nargs = 2 ;
    OGLUsetCoordParms(coord_thickness, coord_spacing) ;
  }
  else if (!stricmp(option, "coord"))
  {
    compile_flags |= COORD_FLAG ;
    coord_fname = argv[2] ;
    fprintf(stderr, "reading coordinate locations from %s.\n", coord_fname);
    nargs = 1 ;
  }
  else if (!stricmp(option, "canon"))
  {
    coord_fname = argv[2] ;
    fprintf(stderr, "reading coordinate locations from %s.\n", coord_fname);
    nargs = 1 ;
  }
  else if (!stricmp(option, "mesh"))
    compile_flags |= MESH_FLAG ;
  else if (!stricmp(option, "tiff"))
    tiff_flag = 1;
  else switch (toupper(*option))
  {
  case 'L':
  {
    LABEL *area ;
    int   i ;

    fprintf(stderr, "reading label file %s\n",argv[2]) ;
    nargs = 1 ;
    area = LabelRead(NULL, argv[2]) ;
    if (nmarked+area->n_points >= MAX_MARKED)
      area->n_points = MAX_MARKED - nmarked ;
    for (i = 0 ; i < area->n_points ; i++)
      marked_vertices[nmarked++] = area->lv[i].vno ;
    LabelFree(&area) ;
  }
  case 'N':
    normalize_flag = 1 ;
    break ;
  case 'S':
    scale = atof(argv[2]) ;
    fprintf(stderr, "scaling window by %2.1f\n", scale) ;
    frame_xdim *= scale ; frame_ydim *= scale ;
    nargs = 1 ;
    break ;
  case 'G':
    curvature_flag = GAUSSIAN_CURVATURE ;
    break ;
  case 'P':
    patch_flag = 1 ;
    compile_flags |= PATCH_FLAG ;
    nargs = 0 ;
    break ;
  case 'V':
    if (nmarked >= MAX_MARKED-1)
      ErrorExit(ERROR_NOMEMORY, "%s: too many vertices marked (%d)",
                Progname, nmarked) ;
    sscanf(argv[2], "%d", &marked_vertices[nmarked++]) ;
    marked_vertices[nmarked] = -1 ;  /* terminate list */
    nargs = 1 ;
    break ;
  case 'E':
    curvature_flag = AREA_ERRORS ;
    break ;
  case 'M':
    medial_flag = 1 ;
    lateral_flag = 0 ;
    break ;
  case 'T':
    talairach_flag = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'A':
    sscanf(argv[2], "%f", &angle_offset) ;
    nargs = 1 ;
    break ;
  case 'X':
    sscanf(argv[2], "%f", &x_angle) ;
    nargs = 1 ;
    break ;
  case 'Y':
    sscanf(argv[2], "%f", &y_angle) ;
    nargs = 1 ;
    break ;
  case 'Z':
    sscanf(argv[2], "%f", &z_angle) ;
    nargs = 1 ;
    break ;
  case 'B':   /* output both medial and lateral views */
    lateral_flag = medial_flag = 1 ;
    break ;
  case 'C':
    curvature_fname =  argv[2] ;
    nargs = 1 ;
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
          "usage: %s [options] <input surface> ... <output directory>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program will convert a sequence of MR surface\n") ;
  fprintf(stderr,"reconstructions\nto SGI rgb files.\nvalid options are:\n\n");
  fprintf(stderr, "-b          output both medial and lateral views\n") ;
  fprintf(stderr, "-c <cslope> change the color slope (default = 5)\n") ;
  fprintf(stderr, "-medial     only output a medial view\n") ;
  fprintf(stderr, "-tiff       output a TIFF file rather than SGI rgb\n");
  
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


static void
xend(void)
{
#if defined(Linux) || defined(__sun__)
  OSMesaDestroyContext(context);
  if (buffer)
  {
    free(buffer) ;
    buffer = NULL ;
  }
#else
  glXDestroyContext(display, context) ;
  glXDestroyGLXPixmap(display, glxpixmap) ;
  XCloseDisplay(display) ;
#endif
}

static void
xinit(void)
{
#if defined(Linux) || defined(__sun__)
  context = OSMesaCreateContext( GL_RGBA, NULL );
  /* Allocate the image buffer */
  buffer = calloc(1, frame_xdim * frame_ydim * 4 );
  
  /* Bind the buffer to the context and make it current */
  OSMesaMakeCurrent( context, buffer,GL_UNSIGNED_BYTE,frame_xdim,frame_ydim);

#else
  display = XOpenDisplay(NULL) ;
  if (!display)
    ErrorExit(ERROR_BADPARM, "could not open display\n") ;
  if (!glXQueryExtension(display, NULL, NULL))
    ErrorExit(ERROR_UNSUPPORTED, "X server has no OpenGL GLX extension\n") ;
  visual = glXChooseVisual(display,DefaultScreen(display),&configuration[1]);
  if (!visual)
  {
    visual = 
      glXChooseVisual(display,DefaultScreen(display),&configuration[0]);
    if (!visual)
      ErrorExit(ERROR_UNSUPPORTED, "could not find an appropriate visual.") ;
  }
  context = glXCreateContext(display, visual, NULL, False) ;
  if (!context)
    ErrorExit(ERROR_UNSUPPORTED, "could not create a drawing context") ;
  
  pixmap = XCreatePixmap(display, RootWindow(display, visual->screen), 
                         frame_xdim, frame_ydim, visual->depth) ;
  glxpixmap = glXCreateGLXPixmap(display, visual, pixmap) ;
  glXMakeCurrent(display, glxpixmap, context) ;
#endif
}

static void
save_rgb(char *fname, int width, int height, unsigned short *red, 
         unsigned short *green, unsigned short *blue)
{
  RGB_IMAGE  *image ;
  int    y ;
  unsigned short *r, *g, *b ;

#ifdef IRIX
  image = iopen(fname,"w",RLE(1), 3, width, height, 3);
#else
  image = iopen(fname,"w",VERBATIM(1), 3, width, height, 3);
#endif
  for(y = 0 ; y < height; y++) 
  {
    r = red + y * width ;
    g = green + y * width ;
    b = blue + y * width ;

    /* fill rbuf, gbuf, and bbuf with pixel values */
    putrow(image, r, y, 0);    /* red row */
    putrow(image, g, y, 1);    /* green row */
    putrow(image, b, y, 2);    /* blue row */
  }
  iclose(image);
}

static void
save_TIFF(char *fname, int width, int height, unsigned char *rgb)
{
  TIFF *out = TIFFOpen(fname,"w");

  if (!out) 
    ErrorExit(ERROR_BADFILE,"Could not open file %s for writing",fname);

  TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 3); 
  TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, height);
  TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

  TIFFWriteEncodedStrip(out, 0, rgb, width*height*3);

  TIFFFlushData(out);
  TIFFClose(out);
}

static void
grabPixelsTIFF(unsigned int width, unsigned int height, unsigned char *rgb)
{
  GLint    swapbytes, lsbfirst, rowlength ;
  GLint    skiprows, skippixels, alignment ;
  unsigned char *tmp, *ttmp, *trgb;
  int i;

  tmp = (unsigned char *)malloc(sizeof(unsigned char)*width*height*3);
  if (!tmp)
    ErrorExit(ERROR_NOMEMORY, "Malloc failed in grabPixelsTIFF\n");

  glGetIntegerv(GL_UNPACK_SWAP_BYTES, &swapbytes) ;
  glGetIntegerv(GL_UNPACK_LSB_FIRST, &lsbfirst) ;
  glGetIntegerv(GL_UNPACK_ROW_LENGTH, &rowlength) ;
  glGetIntegerv(GL_UNPACK_SKIP_ROWS, &skiprows) ;
  glGetIntegerv(GL_UNPACK_SKIP_PIXELS, &skippixels) ;
  glGetIntegerv(GL_UNPACK_ALIGNMENT, &alignment) ;

  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE) ;
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0) ;
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0) ;
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0) ;
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1) ;

  glReadPixels(0, 0, width, height, GL_RGB,GL_UNSIGNED_BYTE, (GLvoid *)tmp);

  glPixelStorei(GL_UNPACK_SWAP_BYTES, swapbytes) ;
  glPixelStorei(GL_UNPACK_LSB_FIRST, lsbfirst) ;
  glPixelStorei(GL_UNPACK_ROW_LENGTH, rowlength) ;
  glPixelStorei(GL_UNPACK_SKIP_ROWS, skiprows) ;
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, skippixels) ;
  glPixelStorei(GL_UNPACK_ALIGNMENT, alignment) ;

  for(i=0;i<height;i++)
    {
      ttmp = tmp + (height-(i+1)) * (width*3);
      trgb = rgb + i*(width*3);
      memcpy(trgb, ttmp, width*3);
    }
  free(tmp);
}

static void
grabPixels(unsigned int width, unsigned int height, unsigned short *red, 
           unsigned short *green, unsigned short *blue)
{
  GLint    swapbytes, lsbfirst, rowlength ;
  GLint    skiprows, skippixels, alignment ;

  glGetIntegerv(GL_UNPACK_SWAP_BYTES, &swapbytes) ;
  glGetIntegerv(GL_UNPACK_LSB_FIRST, &lsbfirst) ;
  glGetIntegerv(GL_UNPACK_ROW_LENGTH, &rowlength) ;
  glGetIntegerv(GL_UNPACK_SKIP_ROWS, &skiprows) ;
  glGetIntegerv(GL_UNPACK_SKIP_PIXELS, &skippixels) ;
  glGetIntegerv(GL_UNPACK_ALIGNMENT, &alignment) ;

  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE) ;
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0) ;
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0) ;
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0) ;
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1) ;

  glReadPixels(0, 0, width, height,GL_RED,GL_UNSIGNED_SHORT, (GLvoid *)red);
  glReadPixels(0, 0,width,height,GL_GREEN,GL_UNSIGNED_SHORT,(GLvoid *)green);
  glReadPixels(0, 0, width, height,GL_BLUE,GL_UNSIGNED_SHORT,(GLvoid *)blue);

  glPixelStorei(GL_UNPACK_SWAP_BYTES, swapbytes) ;
  glPixelStorei(GL_UNPACK_LSB_FIRST, lsbfirst) ;
  glPixelStorei(GL_UNPACK_ROW_LENGTH, rowlength) ;
  glPixelStorei(GL_UNPACK_SKIP_ROWS, skiprows) ;
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, skippixels) ;
  glPixelStorei(GL_UNPACK_ALIGNMENT, alignment) ;
}
static void
clear_pixmaps(MRI_SURFACE *mris)
{
#if defined(Linux) || defined(__sun__)
  OSMesaDestroyContext(context);
  if (buffer)
  {
    free(buffer) ;
    buffer = NULL ;
  }
  xinit() ;   /* create new ones */
  OGLUinit(mris, frame_xdim, frame_ydim) ;  /* reinitialize stuff */
#else
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
#endif
  glLoadIdentity();
}
static int 
mark_centroids(MRI_SURFACE *mris, char *centroid_fnames[], 
               int *centroid_colors, int centroid_files)
{
  FILE    *fp ;
  int     cno, vno, color ;
  float   x, y, z ;
  VERTEX  *v ;
  char    line[100] ;

  for (cno = 0 ; cno < centroid_files ; cno++)
  {
    fp = fopen(centroid_fnames[cno], "r") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open centroid file %s\n",
                Progname, centroid_fnames[cno]) ;

    fgetl(line, 99, fp) ;
    if (sscanf(line, "%f  %f  %f\n", &x, &y, &z) != 3)
      ErrorExit(ERROR_BADFILE, "%s: could not scan area centroid from %s\n",
                Progname, line) ;
    vno = MRISfindClosestCanonicalVertex(mris, x, y, z) ;
    v = &mris->vertices[vno] ;
    color = centroid_colors[cno] ;
    v->marked = color ;
    color++ ;  /* mark the individuals in a different color than the group */
    while (fgetl(line, 99, fp))
    {
      if (sscanf(line, "%f  %f  %f\n", &x, &y, &z) != 3)
        ErrorExit(ERROR_BADFILE, 
                  "%s: could not scan area centroid from %s\n",
                  Progname, line) ;
      vno = MRISfindClosestCanonicalVertex(mris, x, y, z) ;
      v = &mris->vertices[vno] ;
      if (!v->marked)
        v->marked = color ;
    }

    fclose(fp) ;
  }

  return(NO_ERROR) ;
}

