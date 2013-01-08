/**
 * @file  mris_make_map_surfaces.c
 * @brief surface deformation that maximizes likelihood of underlying MRI data

 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2013/01/08 22:01:16 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2013 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>

#include "cma.h"
#include "diag.h"
#include "error.h"
#include "proto.h"
#include "mrisutils.h"
#include "mri.h"
#include "version.h"
#include "utils.h"
#include "const.h"
#include "macros.h"
#include "mrisurf.h"

static char vcid[] =
  "$Id: mris_make_map_surfaces.c,v 1.3 2013/01/08 22:01:16 nicks Exp $";
char *Progname ;
static char sdir[STRLEN] = "" ;
static double l_surf_repulse = 5.0 ;
static int nbrs = 2 ;

static char white_name[STRLEN]  = "white" ;
static char pial_name[STRLEN]  = "pial" ;
static MRI *mri_aseg = NULL ;

static MRI *mask_volume(MRI *mri_src, MRI *mri_aseg, MRI *mri_dst) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static INTEGRATION_PARMS parms ;

static int max_averages = 16 ;
static int min_averages = 2 ;

int
main(int argc, char *argv[])
{
  char           *output_suffix, *subject, fname[STRLEN], *cp, **av ;
  MRI            *mri ;
  MRIS           *mris_lh, *mris_rh, *mris ;
  int            ac, nargs, averages, i ;
  double         orig_dt ;

  char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string
  (argc, argv,
   "$Id: mris_make_map_surfaces.c,v 1.3 2013/01/08 22:01:16 nicks Exp $",
   "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
    (argc, argv,
     "$Id: mris_make_map_surfaces.c,v 1.3 2013/01/08 22:01:16 nicks Exp $",
     "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }

  argc -= nargs;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.l_map = 1 ;
  parms.n_averages = 16 ;
  parms.projection = NO_PROJECTION ;
  parms.tol = 1e-4 ;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.dt = 0.5f ;
  parms.l_curv = .1 ;
  parms.l_tspring = 1.0 ;
  parms.l_nspring = .5 ;
  parms.check_tol = 1 ;
  parms.write_iterations = 0 ;
  parms.l_surf_repulse = 0.0 ;
  parms.l_repulse = 5 ;
  parms.niterations = 1000 ;
  parms.resolution = 0.25 ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
  {
    usage_exit() ;
  }

  subject = argv[1] ;
  printf("deforming surfaces based on input volume %s\n", argv[2]) ;
  mri = MRIread(argv[2]) ;
  if (mri == NULL)
  {
    ErrorExit(ERROR_NOFILE, 
              "%s: could not load MRI volume from %s\n",
              Progname,argv[2]) ;
  }
  output_suffix = argv[3] ;

  if (!strlen(sdir))
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n",
                Progname) ;
    strcpy(sdir, cp) ;
  }

  sprintf(parms.base_name, "%s", output_suffix) ;
  sprintf(fname, "%s/%s/surf/lh.%s", sdir, subject, white_name) ;
  printf("reading input surface from %s\n", fname) ;
  mris_lh = MRISread(fname) ;
  sprintf(fname, "%s/%s/surf/rh.%s", sdir, subject, white_name) ;
  printf("reading input surface from %s\n", fname) ;
  mris_rh = MRISread(fname) ;
  if (nbrs > 1)
  {
    MRISsetNeighborhoodSize(mris_lh, nbrs) ;
    MRISsetNeighborhoodSize(mris_rh, nbrs) ;
  }
  MRIScomputeMetricProperties(mris_lh) ;
  MRIScomputeMetricProperties(mris_rh) ;
  MRISstoreMetricProperties(mris_lh) ;
  MRISstoreMetricProperties(mris_rh) ;
  MRISsaveVertexPositions(mris_lh, WHITE_VERTICES) ;
  MRISsaveVertexPositions(mris_rh, WHITE_VERTICES) ;
  printf("reading pial coordinates from %s\n", pial_name) ;
  MRISreadPialCoordinates(mris_lh, pial_name);
  MRISreadPialCoordinates(mris_rh, pial_name);
  mris = MRISconcat(mris_lh, mris_rh, NULL) ;
  if (nbrs > 1)
  {
    MRISsetNeighborhoodSize(mris, nbrs) ;
  }
  MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
  MRISstoreMetricProperties(mris) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;

  if (mri_aseg)
  {
    mask_volume(mri, mri_aseg, mri) ;
  }

  parms.mri_aseg = mri_aseg ;
  orig_dt = parms.dt ;
  for (i = 0 ; i < 3 ; i++)
  {
    for (averages = max_averages ; averages >= min_averages ; averages /= 2)
    {
      printf("----------- Positioning Surfaces with %d iterative averages ------------------\n", averages) ;
      parms.n_averages = averages ;
//    parms.dt = orig_dt * sqrt(averages) ;
      MRISpositionSurface(mris, mri, mri, &parms) ;
    }
  }

  sprintf(fname, "%s/%s/surf/both.%s", sdir, subject, output_suffix) ;
  printf("writing output to %s\n", fname) ;
  MRISwrite(mris, fname) ;

  return(0) ;
}

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "SDIR"))
  {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "ASEG"))
  {
    mri_aseg = MRIread(argv[2]) ;
    if (mri_aseg == NULL)
    {
      ErrorExit(ERROR_NOFILE, 
                "%s: could not load aseg from %s", 
                Progname, argv[2]) ;
    }
    printf("using aseg from %s to mask cerebellum and "
           "other non-cortical structures\n", argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "DT"))
  {
    parms.dt = atof(argv[2]) ;
    printf("setting integration dt=%2.3f\n", parms.dt) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "MAP"))
  {
    parms.l_map = atof(argv[2]) ;
    printf("set MAP coefficient to %2.3f\n", parms.l_map) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "MAP2D"))
  {
    parms.l_map2d = atof(argv[2]) ;
    printf("set MAP2D coefficient to %2.3f\n", parms.l_map2d) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "TSPRING"))
  {
    parms.l_tspring = atof(argv[2]) ;
    printf("set tspring coefficient to %2.3f\n", parms.l_tspring) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "NSPRING"))
  {
    parms.l_nspring = atof(argv[2]) ;
    printf("set nspring coefficient to %2.3f\n", parms.l_nspring) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "TOL"))
  {
    parms.tol = atof(argv[2]) ;
    printf("set tol to %f\n", parms.tol) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "SPRING"))
  {
    parms.l_spring = atof(argv[2]) ;
    printf("set spring coefficient to %2.3f\n", parms.l_spring) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  }
  else switch (toupper(*option))
    {
    case 'Q':
      parms.flags |= IPFLAG_NO_SELF_INT_TEST ;
      fprintf(stderr,
              "doing quick (no self-intersection) surface positioning.\n") ;
      break ;
    case 'A':
      min_averages = atoi(argv[2]) ;
      max_averages = atoi(argv[3]) ;
      printf("setting n_averages %d --> %d\n", min_averages, max_averages) ;
      nargs = 2 ;
      break ;
    case 'M':
      parms.integration_type = INTEGRATE_MOMENTUM ;
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "momentum = %2.2f\n", parms.momentum) ;
      break ;
    case 'R':
      l_surf_repulse = atof(argv[2]) ;
      fprintf(stderr, "l_surf_repulse = %2.3f\n", l_surf_repulse) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'W':
      sscanf(argv[2], "%d", &parms.write_iterations) ;
      nargs = 1 ;
      fprintf(stderr, "write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'N':
      sscanf(argv[2], "%d", &parms.niterations) ;
      nargs = 1 ;
      fprintf(stderr, "niterations = %d\n", parms.niterations) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      print_help() ;
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

#include "mris_make_map_surfaces.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_make_map_surfaces_help_xml,
                mris_make_map_surfaces_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static MRI *
mask_volume(MRI *mri_src, MRI *mri_aseg, MRI *mri_dst)
{
  int  x, y, z, l ;

  if (mri_dst == NULL)
  {
    mri_dst = MRIcopy(mri_src, NULL) ;
  }

  for (x = 0 ; x < mri_src->width ; x++)
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        l = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        switch (l)
        {
        case Left_Cerebellum_Cortex:
        case Right_Cerebellum_Cortex:
          if ((MRIlabelsInNbhd(mri_aseg, x, y, z, 1, Left_Cerebral_Cortex) == 0) &&
              (MRIlabelsInNbhd(mri_aseg, x, y, z, 1, Right_Cerebral_Cortex) == 0))
          {
            MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
          }
          break ;
        case Brain_Stem:
        case Left_Cerebellum_White_Matter:
        case Right_Cerebellum_White_Matter:
          MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
          break ;
        case Unknown:
          if ((MRIlabelsInNbhd(mri_aseg, x, y, z, 2, Left_Cerebral_Cortex) == 0) &&
              (MRIlabelsInNbhd(mri_aseg, x, y, z, 2, Right_Cerebral_Cortex) == 0))
          {
            MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
          }
          break ;
        default:
          break ;
        }
      }

  return(mri_dst) ;
}

