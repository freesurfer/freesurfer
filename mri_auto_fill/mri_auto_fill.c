#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "transform.h"
#include "timer.h"

int main(int argc, char *argv[]) ;
static BUFTYPE findLabel(MRI *mri, int x0, int y0, int z0) ;
static int get_option(int argc, char *argv[]) ;
static MRI *MRIcombineHemispheres(MRI *mri_lh, MRI *mri_rh, MRI *mri_dst,
                                  MRI *mri_lh_prob, MRI *mri_rh_prob) ;

static int MRIfillVolume(MRI *mri) ;
#if 0
static MRI *MRImaskThreshold(MRI *mri_src, MRI *mri_mask, MRI *mri_dst,
                             float threshold, int out_label) ;
static int MRIgrowLabel(MRI *mri, MRI *mri_bg, int in_label, int out_label) ;
static int MRIturnOnFG(MRI *mri, MRI *mri_fg) ;
static int MRIturnOffBG(MRI *mri, MRI *mri_bg) ;
#endif
static MRI *MRIcheckHemisphereOverlap(MRI *mri_lh, MRI *mri_rh, MRI *mri_dst,
                                     MRI *mri_lh_prob, MRI *mri_rh_prob) ;

char *Progname ;

static void usage_exit(int code) ;

#define LH_FILLED_VOLUME 4
#define RH_FILLED_VOLUME 6

static float pct = 90.0f ;
static int fix = 1 ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs, i ;
  MRI    *mri, *mri_lh_template, *mri_lh_inverse_template, *mri_lh, *mri_rh,
         *mri_rh_template, *mri_rh_inverse_template ;
  char   *in_fname, *template_fname, *out_fname, *xform_fname, fname[100] ;
  M3D    *m3d ;
  struct  timeb start ;
  int     msec ;

  TimerStart(&start) ;

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

  if (argc < 4)
    usage_exit(1) ;

  in_fname = argv[1] ; xform_fname = argv[2] ;
  template_fname = argv[3] ; out_fname = argv[4] ;

  mri = MRIread(in_fname) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s.\n",
              Progname, in_fname) ;

  if (strchr(template_fname, '#') == NULL)
    sprintf(fname, "%s#%d", template_fname, LH_FILLED_VOLUME);
  else
    strcpy(fname, template_fname) ;
  mri_lh_template = MRIread(fname) ;
  if (!mri_lh_template)
    ErrorExit(ERROR_NOFILE, "%s: could not read lh template volume %s.\n",
              Progname, template_fname) ;

  if (strchr(template_fname, '#') == NULL)
    sprintf(fname, "%s#%d", template_fname, RH_FILLED_VOLUME);
  else
    strcpy(fname, template_fname) ;
  mri_rh_template = MRIread(fname) ;
  if (!mri_rh_template)
    ErrorExit(ERROR_NOFILE, "%s: could not read rh template volume %s.\n",
              Progname, template_fname) ;

  fprintf(stderr, "reading transform %s...", xform_fname) ;
  m3d = MRI3DreadSmall(xform_fname) ;
  if (!m3d)
    ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
              Progname, xform_fname) ;
  fprintf(stderr, "done.\n") ;
  fprintf(stderr, "applying inverse transform...\n") ;
  mri_lh_inverse_template = MRIapplyInverse3DMorph(mri_lh_template,m3d,NULL);
  MRIfree(&mri_lh_template) ;
  mri_rh_inverse_template = MRIapplyInverse3DMorph(mri_rh_template,m3d,NULL);
  MRIfree(&mri_rh_template) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    fprintf(stderr, "writing inverse image to lh and rh_inverse.mgh\n") ;
    MRIwrite(mri_lh_inverse_template, "lh_inverse.mgh") ;
    MRIwrite(mri_rh_inverse_template, "rh_inverse.mgh") ;
  }
  MRI3DmorphFree(&m3d) ;
  fprintf(stderr, "using inverse image to fill lh volume...\n") ;
  mri_lh = MRImaskThreshold(mri, mri_lh_inverse_template, NULL, 
                            pct,MRI_LEFT_HEMISPHERE);
  fprintf(stderr, "using inverse image to fill rh volume...\n") ;
  mri_rh = MRImaskThreshold(mri, mri_rh_inverse_template, NULL,
                            pct,MRI_RIGHT_HEMISPHERE);
  MRIfree(&mri) ;
  mri = MRIcombineHemispheres(mri_lh, mri_rh, NULL, mri_lh_inverse_template, 
                              mri_rh_inverse_template) ;
  fprintf(stderr, "done.\n") ;

  if (fix)
    MRIcheckHemisphereOverlap(mri_lh, mri_rh, mri, mri_lh_inverse_template,
                              mri_rh_inverse_template) ;
  MRIfillVolume(mri) ;
  MRIfree(&mri_lh) ; MRIfree(&mri_rh) ;
  MRIfree(&mri_lh_inverse_template) ; MRIfree(&mri_rh_inverse_template) ;
  if (MRIwrite(mri, out_fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not write output volume to %s.\n",
              Progname, out_fname) ;

  msec = TimerStop(&start) ;
  fprintf(stderr, "auto filling took %2.2f minutes\n",
          (float)msec/(1000.0f*60.0f));
  exit(0) ;
  return(0) ;
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
  switch (toupper(*option))
  {
  case 'F':
    fix = 1 ;
    fprintf(stderr, "correcting hemisphere overlap...\n") ;
    break ;
  case 'T':
  case 'P':
    pct = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using threshold = %2.1f%%\n", pct) ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code)
{
  printf("usage: %s <input volume> <transform> <template volume> "
         "<output volume>\n", 
         Progname) ;
  exit(code) ;
}
static MRI *
MRIcombineHemispheres(MRI *mri_lh, MRI *mri_rh, MRI *mri_dst,
                      MRI *mri_lh_prob, MRI *mri_rh_prob)
{
  int       x, y, z, width, height, depth ;
  BUFTYPE   *plh_prob, *prh_prob, *plh, *prh, *pdst,
            lh_prob, rh_prob, lh_label, rh_label, out_label ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_lh, NULL) ;

  width = mri_dst->width ; height = mri_dst->height ; depth = mri_dst->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      plh = &MRIvox(mri_lh, 0, y, z) ;
      prh = &MRIvox(mri_rh, 0, y, z) ;
      plh_prob = &MRIvox(mri_lh_prob, 0, y, z) ;
      prh_prob = &MRIvox(mri_rh_prob, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        lh_label = *plh++ ;
        rh_label = *prh++ ;
        lh_prob =  *plh_prob++ ;
        rh_prob =  *prh_prob++ ;
        if (lh_label && !rh_label)
          out_label = lh_label ;
        else if (rh_label && !lh_label)
          out_label = rh_label ;
        else if (rh_label && lh_label)
        {
          if (rh_prob > lh_prob)
            out_label = rh_label ;
          else
            out_label = lh_label ;
        }
        else
          out_label = 0 ;
        *pdst++ = out_label ;
      }
    }
  }
      
  return(mri_dst) ;
}
/* Talairach seed points for white matter in left and right hemispheres */
#define SEED_SEARCH_SIZE            9
#define MIN_NEIGHBORS               3
#define CORPUS_CALLOSUM_TAL_X       0.0
#define CORPUS_CALLOSUM_TAL_Y       0.0
#define CORPUS_CALLOSUM_TAL_Z       27.0

static int neighbors_on(MRI *mri, int x0, int y0, int z0, int label) ;

static int
MRIfillVolume(MRI *mri)
{
  BUFTYPE  val ;
  Real     wm_rh_tal_x = 29.0 ;
  Real     wm_rh_tal_y = -12.0 ;
  Real     wm_rh_tal_z = 28.0 ;
  Real     wm_lh_tal_x = -29.0 ;
  Real     wm_lh_tal_y = -12.0 ;
  Real     wm_lh_tal_z = 28.0 ;
  int      wm_lh_x, wm_lh_y, wm_lh_z, wm_rh_x, wm_rh_y, wm_rh_z,xi,yi,zi,
           xlim0, xlim1, ylim0, ylim1, zlim0, zlim1, x, y, z, xnew, ynew,znew;
  Real     xr, yr, zr, dist, min_dist, xd, yd, zd ;
  MRI      *mri_fg, *mri_bg, *mri_filled ;

  /* find bounding box for valid data */
  ylim0 = mri->height ; xlim0 = mri->width ;
  ylim1 = xlim1 = 0;
  for (z = 0 ; z < mri->depth ; z++)
  {
    for (y = 0 ; y < mri->height; y++)
    {
      for (x = 0 ; x < mri->width ; x++)
      {
        val = MRIvox(mri, x, y, z) ;
        if (val == MRI_RIGHT_HEMISPHERE || val == MRI_LEFT_HEMISPHERE) ;
        {
          if (y<ylim0) ylim0=y;
          if (y>ylim1) ylim1=y;
          if (x<xlim0) xlim0=x;
          if (x>xlim1) xlim1=x;
        }
      }
    }
  }

  /* find white matter seed point for the left hemisphere */
  MRItalairachToVoxel(mri,CORPUS_CALLOSUM_TAL_X+(SEED_SEARCH_SIZE+1),
                      CORPUS_CALLOSUM_TAL_Y, CORPUS_CALLOSUM_TAL_Z, 
                      &xr, &yr, &zr);
  MRItalairachToVoxel(mri, wm_lh_tal_x,wm_lh_tal_y,wm_lh_tal_z,&xr,&yr,&zr);
  wm_lh_x = nint(xr) ; wm_lh_y = nint(yr) ; wm_lh_z = nint(zr) ;
  if ((MRIvox(mri, wm_lh_x, wm_lh_y, wm_lh_z) != MRI_LEFT_HEMISPHERE) ||
      (neighbors_on(mri, wm_lh_x, wm_lh_y, wm_lh_z, MRI_LEFT_HEMISPHERE) 
       < MIN_NEIGHBORS))
  {
    xnew = ynew = znew = 0 ;
    min_dist = 10000.0f ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "searching for lh wm seed...") ;
    for (z = wm_lh_z-SEED_SEARCH_SIZE ; z <= wm_lh_z+SEED_SEARCH_SIZE ; z++)
    {
      zi = mri->zi[z] ;
      for (y = wm_lh_y-SEED_SEARCH_SIZE ; y <= wm_lh_y+SEED_SEARCH_SIZE ; y++)
      {
        yi = mri->yi[y] ;
        for (x = wm_lh_x-SEED_SEARCH_SIZE;x <= wm_lh_x+SEED_SEARCH_SIZE ; x++)
        {
          xi = mri->xi[x] ;
          if ((MRIvox(mri, xi, yi, zi) == MRI_LEFT_HEMISPHERE) &&
              neighbors_on(mri, xi, yi, zi, MRI_LEFT_HEMISPHERE) >= 
              MIN_NEIGHBORS)
          {
            xd = (xi - wm_lh_x) ; yd = (yi - wm_lh_y) ; zd = (zi - wm_lh_z) ;
            dist = xd*xd + yd*yd + zd*zd ;
            if (dist < min_dist)
            {
              xnew = xi ; ynew = yi ; znew = zi ;
              min_dist = dist ;
            }
          }
        }
      }
    }
    wm_lh_x = xnew ; wm_lh_y = ynew ; wm_lh_z = znew ; 
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "found at (%d, %d, %d)\n", xnew, ynew, znew) ;
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "lh seed point at (%d, %d, %d): %d neighbors on.\n",
            wm_lh_x, wm_lh_y, wm_lh_z, 
            neighbors_on(mri, wm_lh_x, wm_lh_y, wm_lh_z,MRI_LEFT_HEMISPHERE));

  /* find white matter seed point for the right hemisphere */
  MRItalairachToVoxel(mri,CORPUS_CALLOSUM_TAL_X+(SEED_SEARCH_SIZE+1),
                      CORPUS_CALLOSUM_TAL_Y, CORPUS_CALLOSUM_TAL_Z, 
                      &xr, &yr, &zr);
  MRItalairachToVoxel(mri, wm_rh_tal_x,wm_rh_tal_y,wm_rh_tal_z,&xr,&yr,&zr);
  wm_rh_x = nint(xr) ; wm_rh_y = nint(yr) ; wm_rh_z = nint(zr) ;
  if ((MRIvox(mri, wm_rh_x, wm_rh_y, wm_rh_z) == MRI_RIGHT_HEMISPHERE) ||
      (neighbors_on(mri, wm_rh_x, wm_rh_y, wm_rh_z, MRI_RIGHT_HEMISPHERE) 
       < MIN_NEIGHBORS))
  {
    xnew = ynew = znew = 0 ;
    min_dist = 10000.0f ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "searching for rh wm seed...") ;
    for (z = wm_rh_z-SEED_SEARCH_SIZE ; z <= wm_rh_z+SEED_SEARCH_SIZE ; z++)
    {
      zi = mri->zi[z] ;
      for (y = wm_rh_y-SEED_SEARCH_SIZE ; y <= wm_rh_y+SEED_SEARCH_SIZE ; y++)
      {
        yi = mri->yi[y] ;
        for (x = wm_rh_x-SEED_SEARCH_SIZE ;x <= wm_rh_x+SEED_SEARCH_SIZE; x++)
        {
          xi = mri->xi[x] ;
          if ((MRIvox(mri, xi, yi, zi) == MRI_RIGHT_HEMISPHERE) &&
              (neighbors_on(mri, xi, yi, zi, MRI_RIGHT_HEMISPHERE) 
               >= MIN_NEIGHBORS))
          {
            xd = (xi - wm_rh_x) ; yd = (yi - wm_rh_y) ; zd = (zi - wm_rh_z) ;
            dist = xd*xd + yd*yd + zd*zd ;
            if (dist < min_dist)
            {
              xnew = xi ; ynew = yi ; znew = zi ;
              min_dist = dist ;
            }
          }
        }
      }
    }
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "found at (%d, %d, %d)\n", xnew, ynew, znew) ;
    wm_rh_x = xnew ; wm_rh_y = ynew ; wm_rh_z = znew ; 

  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "rh seed point at (%d, %d, %d): %d neighbors on.\n",
            wm_rh_x, wm_rh_y, wm_rh_z, 
            neighbors_on(mri, wm_rh_x, wm_rh_y,wm_rh_z,MRI_RIGHT_HEMISPHERE));

  /* initialize the fill with the detected seed points */
  mri_fg = MRIclone(mri, NULL) ;
  MRIvox(mri_fg, wm_rh_x, wm_rh_y, wm_rh_z) = MRI_RIGHT_HEMISPHERE ;
  MRIvox(mri_fg, wm_lh_x, wm_lh_y, wm_lh_z) = MRI_LEFT_HEMISPHERE ;
  fprintf(stderr, "filling left hemisphere...\n") ;
  MRIgrowLabel(mri, mri_fg, MRI_LEFT_HEMISPHERE, MRI_LEFT_HEMISPHERE) ;
  MRIwrite(mri_fg, "lh.mgh") ;
  fprintf(stderr, "filling right hemisphere...\n") ;
  MRIgrowLabel(mri, mri_fg, MRI_RIGHT_HEMISPHERE, MRI_RIGHT_HEMISPHERE) ;
  MRIwrite(mri_fg, "rh.mgh") ;

  mri_bg = MRIclone(mri, NULL) ;
  MRIvox(mri_bg, 0, 0, 0) = 1 ;
  fprintf(stderr, "filling background...\n") ;
  MRIgrowLabel(mri, mri_bg, 0, 1) ;
  MRIwrite(mri_bg, "bg.mgh") ;
  MRIturnOnFG(mri, mri_fg) ;
  MRIturnOffBG(mri, mri_bg) ;
#if 0
  if (!Gdiag)
    fprintf(stderr, "filling volume: pass 1 of 3...") ;
  fill_brain(WM_MIN_VAL);    /* fill from 2 seeds in wm outwards */

  /*MRIwrite(mri_filled, "fill1.mnc") ;*/
  
  /* set im to initial outward fill */
  MRIcopy(mri_filled, mri) ;
  MRIclear(mri_filled) ;
  MRIvox(mri_filled,1,1,1) = 255;              /* seed for background */
  
  if (!Gdiag)
    fprintf(stderr, "\rfilling volume: pass 2 of 3...") ;
  fill_brain(-WM_MIN_VAL);/* fill in connected component of background */

  if (fill_holes_flag)              
    fill_holes();         /* fill in islands with fewer than 10 nbrs on */
  
  /* put complement into im (im==on means part of background) */
  MRIcopy(mri_filled, mri) ;
  MRIclear(mri_filled) ;

  MRIvox(mri_filled, wm_rh_x, wm_rh_y,wm_rh_z) = RIGHT_HEMISPHERE_WHITE_MATTER;
  MRIvox(mri_filled, wm_lh_x, wm_lh_y, wm_lh_z) = LEFT_HEMISPHERE_WHITE_MATTER;


  /* fill in background of complement (also sometimes called the foreground) */
  if (!Gdiag)
    fprintf(stderr, "\rfilling volume: pass 3 of 3...") ;
  fill_brain(-WM_MIN_VAL);   
  ylim0 = mri->height ; xlim0 = mri->width ;
  ylim1 = xlim1 = 0;
  for (z = 0 ; z < mri_filled->depth ; z++)
  {
    for (y = 0 ; y < mri_filled->height ; y++)
    {
      for (x = 0 ; x < mri_filled->width;x++)
      {
        if (MRIvox(mri_filled, x, y, z) >= WM_MIN_VAL)
        {
          if (y<ylim0) ylim0=y;
          if (y>ylim1) ylim1=y;
          if (x<xlim0) xlim0=x;
          if (x>xlim1) xlim1=x;
        }
      }
    }
  }
  
  if (fill_holes_flag)
    fill_holes();
  
  if (!Gdiag)
    fprintf(stderr, "done.\n") ;
  for (z = 0 ; z < mri_filled->depth ; z++)
    for (y = 0; y < mri_filled->height ; y++)
      for (x = 0; x < mri_filled->width ; x++)
      {
        if (z==0||z==mri_filled->depth-1) 
          MRIvox(mri_filled, x, y, z) = 0;
      }
#endif
  return(NO_ERROR) ;
}
static int
neighbors_on(MRI *mri, int x0, int y0, int z0, int label)
{
  int   nbrs = 0 ;

  if (MRIvox(mri,x0-1,y0,z0) == label)
    nbrs++ ;
  if (MRIvox(mri,x0+1,y0,z0) == label)
    nbrs++ ;
  if (MRIvox(mri,x0,y0+1,z0) == label)
    nbrs++ ;
  if (MRIvox(mri,x0,y0-1,z0) == label)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,z0+1) == label)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,z0-1) == label)
    nbrs++ ;
  return(nbrs) ;
}
#if 0

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static MRI *
MRImaskThreshold(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, float threshold,
                 int out_label)
{
  BUFTYPE   *pmask, *pdst, *psrc, out_val, mask, in_val ;
  int       width, height, depth, x, y, z, nchanged, noff, non ;
  float     nvox ;

  if (mri_mask->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, 
                       "MRI3Dthreshold: mask must be MRI_FLOAT")) ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ; 
  /* now apply the inverse morph to build an average wm representation
     of the input volume 
     */


  noff = non = nchanged = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pmask = &MRIvox(mri_mask, 0, y, z) ;
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        if (x == 131 && y == 97 && z == 127)
          DiagBreak() ;  /* marked as 255, should be 127 */
        out_val = 0 ;
        mask = *pmask++ ;   /* value from inverse morphed volume */
        in_val = *psrc++ ;
        if (mask < 100-threshold)        /* probably off */
          out_val = 0 ;
        else  if (mask > threshold)      /* probably on */
          out_val = out_label ;
        else                             /* not sure, use original fill val */
          out_val = in_val ;
        if (out_val != in_val)
        {
          if (out_val)
            non++ ;
          else
            noff++ ;
          nchanged++ ;
        }
        *pdst++ = out_val ;
      }
    }
  }

  nvox = (float)(width * height * depth) ;
  fprintf(stderr, "%8d of %8d voxels changed - %2.1f%%.\n",
          nchanged, (int)nvox, 100.0f*(float)nchanged/nvox) ;
  fprintf(stderr, "%8d of %8d voxels off     - %2.1f%%.\n",
          noff, (int)nvox,   100.0f*(float)noff/nvox) ;
  fprintf(stderr, "%8d of %8d voxels on      - %2.1f%%.\n",
          nchanged, (int)nvox, 100.0f*(float)non/nvox) ;
  return(mri_dst) ;
}

static BUFTYPE
findLabel(MRI *mri, int x, int y, int z)
{
  int   xi, yi, zi, xk, yk, zk, left_count, right_count, label,
        width, height, depth ;

  width = mri->width ; height = mri->height ; depth = mri->depth ;
  left_count = right_count = 0 ;
  for (zk = -1 ; zk <= 1 ; zk++)
  {
    zi = z + zk ;
    if (zi < 0 || zi >= depth)
      continue ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = y + yk ;
      if (yi < 0 || yi >= height)
        continue ;
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = x + xk ;
        if (xi < 0 || xi >= width)
          continue ;
        label = MRIvox(mri, xi, yi, zi) ;
        if (label == MRI_LEFT_HEMISPHERE)
          left_count++ ;
        else if (label == MRI_RIGHT_HEMISPHERE)
          right_count++ ;
      }
    }
  }
  return ((left_count > right_count) ? 
          MRI_LEFT_HEMISPHERE : 
          MRI_RIGHT_HEMISPHERE) ;
}

static int
MRIgrowLabel(MRI *mri, MRI *mri_filled, int in_label, int out_label)
{
  int      x, y, z, width, height, depth, xi, yi, zi, xk, yk, zk, nfilled,
           total_filled, xmin, xmax, ymin, ymax, zmin, zmax ;

  width = mri->width ; height = mri->height ; depth = mri->depth ;
  total_filled = 0 ;

  xmin = width ; ymin = height ; zmin = depth ;
  xmax = ymax = zmax = 0 ;
  if (in_label) for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == 131 && y == 97 && z == 127)
          DiagBreak() ;
        if (MRIvox(mri, x, y, z) == in_label)
        {
          if (x > xmax)
            xmax = x ;
          if (y > ymax)
            ymax = y ;
          if (z > zmax)
            zmax = z ;
          if (x < xmin)
            xmin = x ;
          if (y < ymin)
            ymin = y ;
          if (z < zmin)
            zmin = z ;
        }
      }
    }
  }
  else   /* filling bg - do outside first (hack, but much faster)*/
  {
    /* find bounding box for real data */
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
        if (x == 131 && y == 97 && z == 127)
          DiagBreak() ;
          if (MRIvox(mri, x, y, z))
          {
            if (x > xmax)
              xmax = x ;
            if (y > ymax)
              ymax = y ;
            if (z > zmax)
              zmax = z ;
            if (x < xmin)
              xmin = x ;
            if (y < ymin)
              ymin = y ;
            if (z < zmin)
              zmin = z ;
          }
        }
      }
    }

    /* fill everything outside it */
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == 131 && y == 97 && z == 127)
            DiagBreak() ;
          if (z <= zmin || z >= zmax || 
              y <= ymin || y >= ymax || 
              x <= xmin || x >= xmax)
          {
            total_filled++ ;
            MRIvox(mri_filled, x, y, z) = out_label ;
          }
        }
      }
    }
  }

#if 0
  xmin = ymin = zmin = 0 ;
  xmax = width-1 ; ymax = height-1 ; zmax = depth-1 ;
#endif

  do
  {
    nfilled = 0 ;
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        for (x = xmin ; x <= xmax ; x++)
        {
          if (MRIvox(mri_filled, x, y, z) == out_label)
          {
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = z+zk ;
              if (zi < 0 || zi >= depth)
                continue ;
              for (yk = -1 ; yk <= 1 ; yk++)
              {
                yi = y+yk ;
                if (yi < 0 || yi >= height)
                  continue ;
                for (xk = -1 ; xk <= 1 ; xk++)
                {
                  xi = x+xk ;
                  if (xi < 0 || xi >= width)
                    continue ;
                  if (xi == 131 && yi == 97 && zi == 127)
                    DiagBreak() ;
                  if ((MRIvox(mri, xi, yi, zi) == in_label) &&
                      (MRIvox(mri_filled, xi, yi, zi) != out_label))
                  {
                    nfilled++ ;
                    MRIvox(mri_filled, xi, yi, zi) = out_label ;
                  }
                }
              }
            }
          }
        }
      }
    }
    total_filled += nfilled ;
    fprintf(stderr, "%d voxels filled, total = %d.\n", nfilled, total_filled);
  } while (nfilled > 0) ;
  return(NO_ERROR) ;
}
static int
MRIturnOnFG(MRI *mri, MRI *mri_fg)
{
  int    x, y, z, width, height, depth ;

  width = mri->width ; height = mri->height ; depth = mri->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (MRIvox(mri_fg, x, y, z) > 0)
          MRIvox(mri, x, y, z) = MRIvox(mri_fg, x, y, z) ;
      }
    }
  }
  return(NO_ERROR) ;
}
/*
  turn off all voxels which are set in the bg image 
  */
static int
MRIturnOffBG(MRI *mri, MRI *mri_bg)
{
  int    x, y, z, width, height, depth ;

  width = mri->width ; height = mri->height ; depth = mri->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (MRIvox(mri_bg, x, y, z) > 0)
          MRIvox(mri, x, y, z) = 0 ;
      }
    }
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if ((MRIvox(mri_bg, x, y, z) == 0) &&
            (MRIvox(mri, x, y, z) == 0))
          MRIvox(mri, x, y, z) = findLabel(mri, x, y, z) ;
      }
    }
  }
  return(NO_ERROR) ;
}
#endif
static MRI *
MRIcheckHemisphereOverlap(MRI *mri_lh, MRI *mri_rh, MRI *mri_dst,
                                     MRI *mri_lh_prob, MRI *mri_rh_prob)
{
  int       x, y, z, width, height, depth ;
  BUFTYPE   *plh_prob, *prh_prob, *plh, *prh, *pdst, *psrc,
            lh_prob, rh_prob, lh_label, rh_label, label ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_lh, NULL) ;

  width = mri_dst->width ; height = mri_dst->height ; depth = mri_dst->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      plh = &MRIvox(mri_lh, 0, y, z) ;
      prh = &MRIvox(mri_rh, 0, y, z) ;
      plh_prob = &MRIvox(mri_lh_prob, 0, y, z) ;
      prh_prob = &MRIvox(mri_rh_prob, 0, y, z) ;
      psrc = &MRIvox(mri_dst, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        if (x == 123 && y == 70 && z == 97)
          DiagBreak() ;  /* marked as 255, should be 127 */
        if (x == 125 && y == 148 && z == 100)
          DiagBreak() ;  /* marked as 255, should be 127 */

        label = *psrc++ ;
        lh_label = *plh++ ;
        rh_label = *prh++ ;
        lh_prob =  *plh_prob++ ;
        rh_prob =  *prh_prob++ ;
        if (label)  /* white matter of one hemi or the other */
        {
#if 0
          if ((label == lh_label) && lh_prob < 10 && rh_prob > 2*lh_prob)
            label = rh_label ;
          else if ((label == rh_label) && rh_prob < 10 && lh_prob > 2*rh_prob)
            label = lh_label ;
#else
          if (lh_prob < 10 && rh_prob > 2*lh_prob)
            label = MRI_RIGHT_HEMISPHERE ;
          else if ((label == rh_label) && rh_prob < 10 && lh_prob > 2*rh_prob)
            label = MRI_LEFT_HEMISPHERE ;
#endif
        }
        *pdst++ = label ;
      }
    }
  }
  return(mri_dst) ;
}

