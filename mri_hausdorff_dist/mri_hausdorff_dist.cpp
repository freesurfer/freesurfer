/**
 * @brief Computes the modified (mean) Hausdorff distance
 *    
 *
 * Computes the modified (mean) Hausdorff distance in an arbitrary set of volumes.
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

#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "const.h"
#include "machine.h"
#include "cma.h"
#include "fio.h"
#include "utils.h"
#include "mri.h"
#include "gcamorph.h"
#include "minc.h"
#include "analyze.h"
#include "mri_identify.h"
#include "error.h"
#include "diag.h"
#include "version.h"
#include "mghendian.h"
#include "fio.h"
#include "cmdargs.h"
#include "macros.h"

static int  get_option(int argc, char *argv[]) ;
static void print_version(void) ;
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static double compute_hdist(MRI **mri, int nvolumes, int index, double *hdists, int which);
static int fromFile = 0;

const char *Progname ;

static int use_vox = 0 ;

#define MAX_VOLUMES 100

static int binarize = 0 ;
static float binarize_thresh = 0 ;
static int target_label = 1 ;
static double blur_sigma = 0 ;

#define MEAN_HDIST 0
#define MAX_HDIST  1

static int which = MEAN_HDIST ;

/***-------------------------------------------------------****/
int main(int argc, char *argv[]) 
{
  int   nargs, nvolumes = 0, n, ac, filecount ;
  char  *name = NULL, fname[STRLEN], *out_fname, **av, *list_fname, in_fname[STRLEN];
  MRI   *mri[MAX_VOLUMES], *mri_tmp ;
  double hdist, hdists[MAX_VOLUMES] ;
  FILE   *fp  ;

  nargs = handleVersionOption(argc, argv, "mri_hausdorff_dist");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) 
    usage_exit();

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  out_fname = argv[argc-1] ;

  if(fromFile)
    {
      list_fname = argv[1] ;
      if (! FileExists(list_fname))
	ErrorExit(Gerror, "%s: fopen(%s) for input volume filelist failed", Progname, list_fname) ;
      else 
	{ 
	  FILE           *fp ;
	  fprintf(stderr, "reading input volume filenames from %s...\n", list_fname) ;
	  fp = fopen(list_fname, "r") ;
	  if (!fp)
	    ErrorExit(Gerror, "Volumelist file cannot be opened\n") ;
	  
	  filecount = 0;
	  while (fscanf(fp,"%s",in_fname) != EOF)
	    {
	      filecount ++;
	    }
          fclose(fp);
	  nvolumes = filecount;
	  if (nvolumes < 2)
	    usage_exit() ;
	  if (nvolumes < 2)
	    ErrorExit(ERROR_BADPARM, "%s: must specify at least 2 input volumes and an output text file",
		      Progname);

	  fprintf(stderr, "processing %d volumes and writing output to %s\n", nvolumes, out_fname) ;
	  fflush(stderr) ;
	  /*Similar to original loop below*/
	  fp = fopen(list_fname, "r") ;
	  for (n = 0 ; n < nvolumes ; n++)
	    {
	      fscanf(fp,"%s",in_fname); 
	      fprintf(stderr, "reading input volume %d from %s\n", n+1, in_fname) ;
	      mri_tmp = MRIread(in_fname) ;
	      if (mri_tmp == NULL)
		ErrorExit(ERROR_BADPARM, "%s: could not read %dth input volume from %s", Progname,n,in_fname);
	      if (mri_tmp->depth == 1)
	      {
		MRI_REGION reg ;
		MRI *mri_tmp2 ;
		reg.x = reg.y = reg.z = 0 ;
		reg.dx = mri_tmp->width ;
		reg.dy = mri_tmp->height ;
		reg.dz = mri_tmp->depth ;
		mri_tmp2 = MRIextractRegionAndPad(mri_tmp, NULL, &reg, 1);
		MRIfree(&mri_tmp) ;
		mri_tmp = mri_tmp2 ;
	      }
	      if (use_vox)
		mri_tmp->xsize = mri_tmp->ysize = mri_tmp->zsize = 1 ;
		  
	      if (mri_tmp->type != MRI_FLOAT)
		{
		  MRI *m ;
		  m = MRIchangeType(mri_tmp, MRI_FLOAT, 0, 1, 1) ;
		  MRIfree(&mri_tmp) ; mri_tmp = m ;
		}
	      if (blur_sigma > 0)
		{
		  MRI *mri_kernel, *mri_smooth ; 
		  mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
		  mri_smooth = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel) ;
		  MRIfree(&mri_kernel) ; MRIfree(&mri_tmp) ;
		  mri_tmp = mri_smooth ;
		}
	      
	      sprintf(fname, "d%d.mgz", n) ;
#define USE_DISTANCE_TRANSFORM 1
#if USE_DISTANCE_TRANSFORM 
	      if (binarize)
		{
		  MRI *mri_tmp2 ;
		  mri_tmp2 = MRIbinarize(mri_tmp, NULL, binarize_thresh, 0, target_label) ;
		  MRIfree(&mri_tmp) ;
		  mri_tmp = mri_tmp2 ;
		}
	      mri[n] = MRIdistanceTransform(mri_tmp, NULL, target_label, 
					    mri_tmp->width+mri_tmp->height+mri_tmp->depth, 
					    DTRANS_MODE_SIGNED, NULL) ;
	      MRIfree(&mri_tmp) ;
	      //    MRIwrite(mri[n], fname) ;
#else
	      mri[n] = mri_tmp ;
#endif
	    }
	  fclose(fp);	 
	  /*End similar of the original loop*/
	}
    }
  else
   {
      nvolumes = argc-2 ;  // last argument is output file name
      if (nvolumes < 2)
	usage_exit() ;
      if (nvolumes < 2)
	ErrorExit(ERROR_BADPARM, "%s: must specify at least 2 input volumes and an output text file",
		  Progname);
      //      out_fname = argv[nvolumes+1] ;
      fprintf(stderr, "processing %d volumes and writing output to %s\n", nvolumes, out_fname) ;
      fflush(stderr) ;

      for (n = 0 ; n < nvolumes ; n++)
	{
	  name = argv[n+1] ;
	  fprintf(stderr, "reading input volume %d from %s\n", n+1, name) ;
	  mri_tmp = MRIread(name) ;
	  if (mri_tmp == NULL)
	    ErrorExit(ERROR_BADPARM, "%s: could not read %dth input volume from %s", Progname,n,name);
	  if (mri_tmp->depth == 1)
	  {
	    MRI_REGION reg ;
	    MRI *mri_tmp2 ;
	    reg.x = reg.y = reg.z = 0 ;
	    reg.dx = mri_tmp->width ;
	    reg.dy = mri_tmp->height ;
	    reg.dz = mri_tmp->depth ;
	    mri_tmp2 = MRIextractRegionAndPad(mri_tmp, NULL, &reg, 1);
	    MRIfree(&mri_tmp) ;
	    mri_tmp = mri_tmp2 ;
	  }
		  
	  if (use_vox)
	    mri_tmp->xsize = mri_tmp->ysize = mri_tmp->zsize = 1 ;
	  if (mri_tmp->type != MRI_FLOAT)
	    {
	      MRI *m ;
	      m = MRIchangeType(mri_tmp, MRI_FLOAT, 0, 1, 1) ;
	      MRIfree(&mri_tmp) ; mri_tmp = m ;
	    }
	  if (blur_sigma > 0)
	    {
	      MRI *mri_kernel, *mri_smooth ; 
	      mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
	      mri_smooth = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel) ;
	      MRIfree(&mri_kernel) ; MRIfree(&mri_tmp) ;
	      mri_tmp = mri_smooth ;
	    }

	  sprintf(fname, "d%d.mgz", n) ;
#define USE_DISTANCE_TRANSFORM 1
#if USE_DISTANCE_TRANSFORM 
	  if (binarize)
	    {
	      MRI *mri_tmp2 ;
	      mri_tmp2 = MRIbinarize(mri_tmp, NULL, binarize_thresh, 0, target_label) ;
	      MRIfree(&mri_tmp) ;
	      mri_tmp = mri_tmp2 ;
	    }
	  mri[n] = MRIdistanceTransform(mri_tmp, NULL, target_label, 
					mri_tmp->width+mri_tmp->height+mri_tmp->depth, 
					DTRANS_MODE_SIGNED, NULL) ;
	  MRIfree(&mri_tmp) ;
	  //    MRIwrite(mri[n], fname) ;
#else
	  mri[n] = mri_tmp ;
#endif
	}
    }


  fp = fopen(out_fname, "w") ;
  if (fp == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open %s", Progname, out_fname) ;

#if 1
  hdist = compute_hdist(mri, nvolumes, 0, hdists, which) ;
  if (hdist < 0)
    DiagBreak() ;
  for (n = 1 ; n < nvolumes ; n++)
  {
    if (hdists[n] < 0)
      DiagBreak() ;
    fprintf(fp, "%f\n", hdists[n]) ;
  }
#else
  for (n = 0 ; n < nvolumes ; n++)
  {
    hdist = compute_hdist(mri, nvolumes, n, hdists, which) ;
    fprintf(fp, "%f\n", hdist) ;
    fflush(fp) ;
  }
#endif
  fclose(fp) ;
  exit(0);

} /* end main() */

/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s <options> <vol1> <vol2> ... <output text file> \n",Progname) ;
  printf("\twhere options are:\n") ;
  printf("\t-b <thresh>        binarize input volumes with threshold <thresh>\n") ;
  printf("\t-F                 read volumes from an input file (first argument is the input filename)\n") ;
  printf("\t-g <sigma>         blur the input image with Gaussian <sigma>\n");
  printf("\t-max               compute the max of the min distances instead of the mean\n") ;
  printf("\t-l <label index>   use <label index> as the target label\n") ;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf(
    "\n"
    "Computes the mean of the min distances between point sets\n") ;

  exit(1) ;
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

#if USE_DISTANCE_TRANSFORM 
static double
compute_hdist(MRI **mri, int nvolumes, int index, double *hdists, int which)
{
  int   x, y, z, width, depth, height, nvox, n, xk, yk, zk, xi, yi, zi ;
  float d, d2, dist, dx, dy, dz ;
  MRI   *mri_src ;
  double hdists_sigma[MAX_VOLUMES], hdist, max_vox,xf, yf, zf, zval, max_hdist, max_hdists[MAX_VOLUMES];
  //  FILE   *fp ;
  static int i = 0 ;
  char fname[STRLEN] ;

  sprintf(fname, "hdists%d.txt", i++) ;
  //  fp = fopen(fname, "w") ;

  mri_src = mri[index] ;

  width = mri_src->width ; height = mri_src->height ;  depth = mri_src->depth ; 

  max_hdist = 0.0 ;
  max_vox = MAX(mri_src->xsize, MAX(mri_src->ysize, mri_src->zsize)) ;
  for (hdist = 0.0, n = 0 ; n < nvolumes ; n++)
  {
    if (n == index)
      continue ;
    for (max_hdists[n] = hdists_sigma[n] = hdists[n] = 0.0, nvox = x = 0 ; x < width ; x++)
      for (y = 0 ; y < height ; y++)
        for (z = 0 ; z < depth ; z++)
        {
          d = MRIgetVoxVal(mri_src, x, y, z, 0) ;
          if (fabs(d) > max_vox)
            continue ;

          // find locations on either side of 0
          for (xk = -1 ; xk <= 1 ; xk++)
          {
            xi = mri_src->xi[x+xk] ;
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_src->yi[y+yk] ;
              for (zk = -1 ; zk <= 1 ; zk++)
              {
                if (fabs(xk) + fabs(yk) + fabs(zk) != 1)
                  continue ;

                zi = mri_src->zi[z+zk] ;
                
                d2 = MRIgetVoxVal(mri_src, xi, yi, zi, 0) ;
                if (d * d2 > 0)
                  continue ;  // not on either side of 0

                // compute dist in mm to nbr voxel
                dx = (xi - x) * mri_src->xsize ;
                dy = (yi - y) * mri_src->ysize ;
                dz = (zi - z) * mri_src->zsize ;
                dist = sqrt(dx*dx + dy*dy + dz*dz);

                // convert to voxel coords of 0-crossing
                dx = (xi - x) ; dy = (yi - y) ; dz = (zi - z) ;
                xf = x + dx * fabs(d2)/dist ;
                yf = y + dy * fabs(d2)/dist ;
                zf = z + dz * fabs(d2)/dist ;
		//printf("mri_hausdorff: dist = %f\n", dist);
                MRIsampleVolume(mri[n], xf, yf, zf, &zval) ;
                if (zval < 0)
                  DiagBreak() ;
                zval = fabs(zval) ;
                //                fprintf(fp, "%d %d %d %d %2.3f\n", n, x, y, z, d) ;
                if (zval > 40)
                  DiagBreak() ;
                if (zval > max_hdists[n])
                  max_hdists[n] = zval ;
                hdists[n] += zval ;
                hdists_sigma[n] += zval*zval ;
                nvox++ ;
              }
            }
          }
        }
    if (which == MAX_HDIST)
      hdists[n] = max_hdists[n] ;
    else {
      if (nvox > 0)
	hdists[n] /= (double)nvox ;
      else hdists[n] = 0;
    }
    hdists_sigma[n] = sqrt(hdists_sigma[n]/(double)nvox - hdists[n]*hdists[n]);
    hdist += hdists[n] ;
  }

  if (which == MAX_HDIST)
    hdist = max_hdist ;
  else
    hdist /= (nvolumes-1) ;
  //  fclose(fp) ;
  return(hdist) ;
}

#else
#include "voxlist.h"
static double
compute_hdist(MRI **mri, int nvolumes, int index)
{
  int   i1, i2, width, depth, height, n ;
  float dist, min_dist ;
  double hdists[MAX_VOLUMES], hdists_sigma[MAX_VOLUMES], hdist ;
  FILE   *fp ;
  static int i = 0 ;
  char fname[STRLEN] ;
  VOXEL_LIST  *vl[MAX_VOLUMES], *vl1, *vl2 ;

  sprintf(fname, "hdists%d.txt", i++) ;
  fp = fopen(fname, "w") ;
  if (fp == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open %s", Progname, fname) ;

  width = mri[index]->width ; height = mri[index]->height ;  depth = mri[index]->depth ; 

  for (n = 0 ; n < nvolumes ; n++)
    vl[n] = VLSTcreate(mri[n], 1, 256, NULL, 0, 1) ;
  vl1 = vl[index] ;
  for (hdist = 0.0, n = 0 ; n < nvolumes ; n++)
  {
    if (n == index)
      continue ;
    vl2 = vl[n] ;
    for (i1 = 0; i1 < vl1->nvox ; i1++)
    {
      min_dist = width+height+depth ;
      for (i2 = 0; i2 < vl2->nvox ; i2++)
      {
        dist = sqrt(SQR(vl1->xi[i1]-vl2->xi[i2])+
                    SQR(vl1->yi[i1]-vl2->yi[i2])+
                    SQR(vl1->zi[i1]-vl2->zi[i2])) ;
        if (dist < min_dist)
          min_dist = dist ;
      }
      if (min_dist > 40)
        DiagBreak() ;
      hdists[n] += min_dist ;
      hdists_sigma[n] += min_dist*min_dist ;
    }
    if (vl1->nvox>0)
      hdists[n] /= (double)vl1->nvox ;
    else
      hdists[n] = 0 ;
    hdists_sigma[n] = sqrt(hdists_sigma[n]/(double)vl1->nvox - hdists[n]*hdists[n]);
    hdist += hdists[n] ;
  }

  for (n = 0 ; n < nvolumes ; n++)
    VLSTfree(&vl[n]) ;
  hdist /= (nvolumes-1) ;
  fclose(fp) ;
  return(hdist) ;
}
#endif
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "max"))
  {
    printf("using max of min distances\n") ;
    which = MAX_HDIST ;
  } else switch (toupper(*option)) {
  case 'G':
    blur_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("blurring input image with sigma = %2.2f\n", blur_sigma) ;
    break ;
  case 'L':
    target_label = atoi(argv[2]) ;
    printf("using %d (%s) as target label\n", target_label,
           cma_label_to_name(target_label)) ;
    nargs = 1 ;
    break ;
  case 'B':
    binarize = 1 ;
    binarize_thresh = atof(argv[2]) ;
    printf("binarizing input data with threshold %2.2f\n", binarize_thresh) ;
    nargs = 1 ;
    break ;
  case 'F':
    fromFile = 1 ;
    printf("read input volumes from file\n") ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'V':
    use_vox = 1 ;
    printf("ignoring voxel sizes to compute distances in voxel coords") ;
    break;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }
  return(nargs) ;
}
static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}
