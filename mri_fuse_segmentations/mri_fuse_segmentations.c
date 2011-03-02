/**
 * @file  mri_fuse_segmentations.c
 * @brief fuse a set of segmentations (asegs)
 *
 * program to fuse a group of cross sectional segmentations into 
 * an initial estimate of a longitudinal one
 * See Sabuncu et al., MICCA 2009 (SNIP paper).
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:15 $
 *    $Revision: 1.8 $
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
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "cma.h"
#include "transform.h"


int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;

static char *norm_name = "norm.mgz" ;
static char *aseg_name = "aseg.mgz" ;
static char *aseg_noCC_name = "aseg.auto_noCCseg.mgz" ;
static char sdir[STRLEN] = "" ;

static MRI *MRIfuseSegmentations(MRI *mri_in,
                                 MRI *mri_fused,
                                 int nvols,
                                 LTA *xforms[],
                                 MRI *mri_asegs[],
                                 MRI *mri_asegs_noCC[],
                                 MRI *mri_norms[],
                                 double std) ;

#define MAX_VOLUMES 1000
static int use_identity = 0 ;
int
main(int argc, char *argv[])
{
  char        **av, fname[STRLEN], *cp, *out_fname, *subject ;
  int          ac, nargs, i, ntimepoints, msec, minutes, seconds ;
  struct timeb start ;
  MRI         *mri_in, *mri_fused;
  MRI         *mri_norms[MAX_VOLUMES];
  MRI         *mri_asegs[MAX_VOLUMES];
  MRI         *mri_asegs_noCC[MAX_VOLUMES];
  double      std = 3.0 ;
  LTA         *xforms[MAX_VOLUMES] ;

  /* rkt: check for and handle version tag */
  nargs = 
    handle_version_option
    (argc, argv,
     "$Id: mri_fuse_segmentations.c,v 1.8 2011/03/02 00:04:15 nicks Exp $",
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;

  subject = argv[1] ;
  out_fname = argv[argc-1] ;
  if (!strlen(sdir))
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  ntimepoints = argc-3 ;
  sprintf(fname, "%s/%s/mri/%s", sdir, subject, norm_name) ;
  printf("reading %s and processing %d timepoints, writing to %s\n",
         fname, ntimepoints,out_fname) ;

  mri_in = MRIread(fname) ;
  if (mri_in == NULL)
    ErrorExit(ERROR_NOFILE,
              "%s: could not read input volume %s", Progname, fname) ;
  for (i = 0 ; i < ntimepoints ; i++)
  {
    printf("loading data for subject %s ...\n", argv[i+2]) ;

    sprintf(fname, "%s/%s/mri/%s", sdir, argv[i+2], aseg_name) ;
    mri_asegs[i] = MRIread(fname) ;
    if (mri_asegs[i] == NULL)
      ErrorExit(ERROR_NOFILE,
                "%s: could not read input volume %s", Progname, fname) ;

    sprintf(fname, "%s/%s/mri/%s", sdir, argv[i+2], aseg_noCC_name) ;
    mri_asegs_noCC[i] = MRIread(fname) ;
    if (mri_asegs_noCC[i] == NULL)
      ErrorExit(ERROR_NOFILE,
                "%s: could not read input volume %s", Progname, fname) ;

    sprintf(fname, "%s/%s/mri/%s", sdir, argv[i+2], norm_name) ;
    mri_norms[i] = MRIread(fname) ;
    if (mri_norms[i] == NULL)
      ErrorExit(ERROR_NOFILE,
                "%s: could not read input volume %s", Progname, fname) ;

    if (use_identity)
      xforms[i] = LTAalloc(1, NULL) ;
    else
    {
      sprintf(fname, "%s/%s/mri/transforms/%s_to_%s.lta",
              sdir, subject, argv[i+2], subject) ;
      xforms[i] = LTAread(fname) ;
      if (xforms[i] == NULL)
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read input transform %s", Progname, fname) ;
    }
  }

  printf("fusing...\n");
  mri_fused = MRIfuseSegmentations(mri_in,
                                   NULL,
                                   ntimepoints,
                                   xforms,
                                   mri_asegs,
                                   mri_asegs_noCC,
                                   mri_norms,
                                   std) ;

  printf("writing fused segmentation to %s\n", out_fname) ;
  MRIwrite(mri_fused, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr,
          "fusing took %d minutes and %d seconds.\n", minutes, seconds) ;

  MRIfree(&mri_in);
  MRIfree(&mri_fused);
  for (i = 0 ; i < ntimepoints ; i++)
  {
    MRIfree(&mri_asegs[i]);
    MRIfree(&mri_asegs_noCC[i]);
    MRIfree(&mri_norms[i]);
    LTAfree(&xforms[i]);
  }

  exit(0) ;
  return(0) ;
}



static int get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "SDIR"))
  {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
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
  else if (!stricmp(option, "identity"))
  {
    use_identity = 1 ;
    printf("using identity matrix as xform\n") ;
  }
  else switch (toupper(*option))
    {
    case 'N':
      norm_name = argv[2] ;
      printf("using %s as norm volume\n", norm_name) ;
      nargs = 1 ;
      break ;
    case 'A':
      aseg_name = argv[2] ;
      printf("using %s as aseg volume\n", aseg_name) ;
      nargs = 1 ;
      break ;
    case 'C':
      aseg_noCC_name = argv[2] ;
      printf("using %s as aseg no-CC-label volume\n", aseg_noCC_name) ;
      nargs = 1 ;
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

static void usage_exit(int code)
{
  printf("usage: %s [options] <subject to segment> "
         "<subject 1> <subject 2> <...> <output volume>\n",
         Progname) ;
  printf("\noptions:\n");
  printf("  -a <filename>  - name of aseg file to use (default: aseg.mgz)\n") ;
  printf("  -c <filename>  - name of aseg file w/o CC labels "
         "(default: aseg.auto_noCCseg.mgz)\n") ;
  printf("  -n <filename>  - name of norm file to use (default: norm.mgz)\n") ;
  printf("  -SDIR SUBJECTS_DIR  - alternate SUBJECTS_DIR\n\n");
  printf("example:\n");
  printf("  %s tp3.long.longbase tp1 tp2 tp3 aseg.fused.mgz\n\n",Progname);
  exit(code) ;
}


static MRI *MRIfuseSegmentations(MRI *mri_in,
                                 MRI *mri_fused,
                                 int nvols,
                                 LTA *xforms[],
                                 MRI *mri_asegs[],
                                 MRI *mri_asegs_noCC[],
                                 MRI *mri_norms[],
                                 double sigma)
{
  int    x, y, z, i, label, min_label, max_label, label_counts[MAX_CMA_LABELS],
         total;
  MATRIX *m_vox2vox[MAX_VOLUMES], *m ;
  VECTOR *v1, *v2 ;
  float  xd, yd, zd ;
  double label_pvals[MAX_CMA_LABELS], p, val, oval, dif ;
  int cc_relabel_count = 0;

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = 1.0 ;
  VECTOR_ELT(v2, 4) = 1.0 ;

  for (i = 0 ; i < nvols ; i++)
  {
    switch (xforms[i]->type)
    {
    case LINEAR_RAS_TO_RAS:
      m = MRIrasXformToVoxelXform(mri_norms[i],
                                  mri_in,
                                  xforms[i]->xforms[0].m_L,
                                  NULL);
      break ;
    case LINEAR_VOX_TO_VOX:
      m = MatrixCopy(xforms[i]->xforms[0].m_L, NULL);
      break ;
    default:
      ErrorExit(ERROR_UNSUPPORTED,
                "unsupported transform type %d\n", xforms[i]->type) ;
      break ;
    }
    m_vox2vox[i] = MatrixInverse(m, NULL);
    MatrixFree(&m) ;
  }
  memset(label_pvals, 0, sizeof(label_pvals)) ;
  memset(label_counts, 0, sizeof(label_counts)) ;
  if (mri_fused == NULL)
    mri_fused = MRIclone(mri_in, NULL) ;

  for (x = 0 ; x < mri_fused->width ; x++)
  {
    V3_X(v1) = x ;
    for (y = 0 ; y < mri_fused->height ; y++)
    {
      V3_Y(v1) = y ;
      for (z = 0 ; z < mri_fused->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        V3_Z(v1) = z ;
        val = MRIgetVoxVal(mri_in, x, y, z, 0) ;
        min_label = MAX_CMA_LABELS ; 
        max_label = 0 ;
        total = 0 ;
        for (i = 0 ; i < nvols ; i++)
        {
          MatrixMultiply(m_vox2vox[i], v1, v2) ;
          xd = V3_X(v2) ;
          yd = V3_Y(v2) ;
          zd = V3_Z(v2) ;
          if (xd < 0 || yd < 0 || zd < 0 || 
              xd > mri_asegs[i]->width-1 || 
              yd > mri_asegs[i]->height-1 || 
              zd > mri_asegs[i]->depth-1)
            continue ; // out of the FOV
          MRIsampleVolume(mri_norms[i], xd, yd, zd, &oval) ;
          label = MRIgetVoxVal(mri_asegs[i], nint(xd), nint(yd), nint(zd), 0) ;
          switch (label)
          {
          case CC_Posterior:
          case CC_Mid_Posterior:
          case CC_Central:
          case CC_Mid_Anterior:
          case CC_Anterior:
          {
            // need to relabel these back to white matter, from noCC aseg
            // (downstream mri_cc will re-seg the CC)
            int label_noCC =
              MRIgetVoxVal(mri_asegs_noCC[i],nint(xd),nint(yd),nint(zd),0);
            if ((label_noCC == Left_Cerebral_White_Matter) ||
                (label_noCC == Right_Cerebral_White_Matter))
            {
              //printf("relabeling %d to %d\n",label, label_noCC);
              label = label_noCC;
              cc_relabel_count++;
            }
            break;
          }
          default:
            break;
          }
          if (label < min_label)
            min_label = label ;
          if (label > max_label)
            max_label = label ;
          dif = oval - val ;
          p = (1.0 / (sqrt(2*M_PI)*sigma)) *
              exp(-0.5 * (dif*dif) / (sigma*sigma)) ;
          label_pvals[label] += p ;
          label_counts[label]++ ;
          total++ ;
        }
        p = 0 ;
        for (label = min_label ; label <= max_label ; label++)
        {
          //          label_pvals[label] *= (double)label_counts[label]/(double)total ;
          label_pvals[label] /= (double)total ;  // prior is already done by sum
          if (label_pvals[label] > p)
          {
            p = label_pvals[label] ;
            MRIsetVoxVal(mri_fused, x, y, z, 0, label) ;
          }
          label_pvals[label] = 0 ; // reset for next time to avoid big memset
          label_counts[label] = 0 ;
        }
      }
    }
  }

  // confirm that no CC labels still exist in the fused volume
  for (x = 0 ; x < mri_fused->width ; x++)
  {
    for (y = 0 ; y < mri_fused->height ; y++)
    {
      for (z = 0 ; z < mri_fused->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;

        label = MRIgetVoxVal(mri_fused, x, y, z, 0) ;
        switch (label)
        {
        case CC_Posterior:
        case CC_Mid_Posterior:
        case CC_Central:
        case CC_Mid_Anterior:
        case CC_Anterior:
        {
          //printf("WARNING: found CC label: %d, xyz: %d %d %d \n", 
          //     label, x, y, z);
          // looks like there's a CC voxel that wasnt relabeled 
          // look at neighbors to determine who is more predominant: left or
          // right white matter voxels; predominance determines relabeling
          int leftwm=0;
          int rightwm=0;
          int xk, yk, zk, xi, yi, zi;
          // brute force 3x3 search:
          for (zk = -1 ; zk <= 1 ; zk++)
          {
            zi = mri_fused->zi[z+zk] ;
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_fused->yi[y+yk] ;
              for (xk = -1 ; xk <= 1 ; xk++)
              {
                xi = mri_fused->xi[x+xk] ;
                if (!zk && !yk && !xk)
                  continue ;
                int lbl = MRIgetVoxVal(mri_fused, xi, yi, zi, 0);
                if (lbl == Left_Cerebral_White_Matter) leftwm++;
                else if (lbl == Right_Cerebral_White_Matter) rightwm++;
              }
            }
          }
          //printf("Neighbors: leftwm=%d, rightwm=%d\n", leftwm, rightwm);
          if (leftwm > rightwm)
          {
            MRIsetVoxVal(mri_fused, x, y, z, 0, Left_Cerebral_White_Matter) ;
            cc_relabel_count++;
          }
          else if (rightwm > 0)
          {
            MRIsetVoxVal(mri_fused, x, y, z, 0, Right_Cerebral_White_Matter) ;
            cc_relabel_count++;
          }
          else
          {
            printf("ERROR: could not determine a label for relabeling\n");
            printf("CC label (%d) at xyz %d %d %d\n",label, x, y, z);
          }
          
          break;
        }
        default:
          break;
        }
      }
    }
  }
  printf("relabeled %d CC labels\n",cc_relabel_count);

  for (i = 0 ; i < nvols ; i++)
  {
    MatrixFree(&m_vox2vox[i]) ;
  }
  VectorFree(&v1);
  VectorFree(&v2);

  return(mri_fused) ;
}

