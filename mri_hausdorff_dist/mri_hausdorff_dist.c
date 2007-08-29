/**
 * @file  mri_hausdorff_dist.c
 * @brief Computes the modified (mean) Hausdorff distance
 *    
 *
 * Computes the modified (mean) Hausdorff distance in an arbitrary set of volumes.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/08/29 13:47:36 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

char *MRI_HAUSDORFF_DIST_VERSION = "$Revision: 1.1 $";

#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "const.h"
#include "machine.h"
#include "fio.h"
#include "utils.h"
#include "mri.h"
#include "gcamorph.h"
#include "volume_io.h"
#include "analyze.h"
#include "mri_identify.h"
#include "error.h"
#include "diag.h"
#include "version.h"
#include "mghendian.h"
#include "fio.h"
#include "cmdargs.h"
#include "macros.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static double compute_hdist(MRI **mri, int nvolumes, int index) ;
static char vcid[] = "$Id: mri_hausdorff_dist.c,v 1.1 2007/08/29 13:47:36 fischl Exp $";

char *Progname ;

#define MAX_VOLUMES 100

/***-------------------------------------------------------****/
int main(int argc, char *argv[]) 
{
  int   nargs, nvolumes, n ;
  char  *name, fname[STRLEN], *out_fname;
  MRI   *mri[MAX_VOLUMES], *mri_tmp ;
  double hdist ;
  FILE   *fp  ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) 
    usage_exit();

  parse_commandline(argc, argv);
  check_options();

  nvolumes = argc-1 ;
  out_fname = argv[nvolumes] ;
  fprintf(stderr, "processing %d volumes and writing output to %s\n", nvolumes, out_fname) ;

  for (n = 0 ; n < nvolumes ; n++)
  {
    name = argv[n] ;
    fprintf(stderr, "reading input volume %d from %s\n", n+1, name) ;
    mri_tmp = MRIread(name) ;
    if (mri_tmp == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not read %dth input volume from %s", Progname,n,name);
    sprintf(fname, "d%d.mgz", n) ;
#define USE_DISTANCE_TRANSFORM 1
#if USE_DISTANCE_TRANSFORM 
    mri[n] = MRIdistanceTransform(mri_tmp, NULL, 1, 100, DTRANS_MODE_SIGNED) ;
    MRIfree(&mri_tmp) ;
    //    MRIwrite(mri[n], fname) ;
#else
    mri[n] = mri_tmp ;
#endif
  }


  fp = fopen(out_fname, "w") ;
  for (n = 0 ; n < nvolumes ; n++)
  {
    hdist = compute_hdist(mri, nvolumes, n) ;
    fprintf(fp, "%f\n", hdist) ;
    fflush(fp) ;
  }
  fclose(fp) ;
  exit(0);

} /* end main() */

/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  
      print_help() ;
    else {
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s fname1 <fname2> <options> \n",Progname) ;
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
/* --------------------------------------------- */
static void check_options(void) {
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

#if USE_DISTANCE_TRANSFORM 
static double
compute_hdist(MRI **mri, int nvolumes, int index)
{
  int   x, y, z, width, depth, height, nvox, n ;
  float d ;
  MRI   *mri_src ;
  double hdists[MAX_VOLUMES], hdists_sigma[MAX_VOLUMES], hdist ;
  FILE   *fp ;
  static int i = 0 ;
  char fname[STRLEN] ;

  sprintf(fname, "hdists%d.txt", i++) ;
  fp = fopen(fname, "w") ;

  mri_src = mri[index] ;

  width = mri_src->width ; height = mri_src->height ;  depth = mri_src->depth ; 

  for (hdist = 0.0, n = 0 ; n < nvolumes ; n++)
  {
    if (n == index)
      continue ;
    for (hdists_sigma[n] = hdists[n] = 0.0, nvox = x = 0 ; x < width ; x++)
      for (y = 0 ; y < height ; y++)
        for (z = 0 ; z < depth ; z++)
        {
          d = MRIgetVoxVal(mri_src, x, y, z, 0) ;
          if (FEQUAL(d, -0.5))
          {
            d = MRIgetVoxVal(mri[n], x, y, z, 0) ;
            fprintf(fp, "%d %d %d %d %2.3f\n", n, x, y, z, d) ;
            if (d > 40)
              DiagBreak() ;
            d = fabs(d-0.5) ;
            hdists[n] += d ;
            hdists_sigma[n] += d*d ;
            nvox++ ;
          }
        }
    hdists[n] /= (double)nvox ;
    hdists_sigma[n] = sqrt(hdists_sigma[n]/(double)nvox - hdists[n]*hdists[n]);
    hdist += hdists[n] ;
  }

  hdist /= (nvolumes-1) ;
  fclose(fp) ;
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
    hdists[n] /= (double)vl1->nvox ;
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
