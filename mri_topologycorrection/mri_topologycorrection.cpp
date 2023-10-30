/**
 * @brief topology correction routines
 *
 * "A Hybrid Approach to the Skull-Stripping Problem in MRI",
 * Ségonne, F., Dale, A.M., Busa, E., Glessner, M., Salvolini, U.,
 * Hahn, H.K., Fischl, B.
 * (2004) NeuroImage, 22:1160-1075.
 */
/*
 * Original Author: F. Segonne
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
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "mrisurf.h"
#include "matrix.h"
#include "proto.h"
#include "stats.h"
#include "timer.h"
#include "const.h"
#include "mrishash.h"
#include "icosahedron.h"
#include "tritri.h"
#include "timer.h"
#include "diag.h"
#include "mri_topology.h"
#include "gca.h"
#include "mri_tess.h"
#include "mrisutils.h"
#include "cma.h"

#if 0
static int
check_volume(MRI *mri_save, MRI *mri_out, int target_label) {
  int   x, y, z, sval, oval ;

  for (x = 0 ; x < mri_save->width ; x++)
    for (y = 0 ; y < mri_save->height ; y++)
      for (z = 0 ; z < mri_save->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        sval = nint(MRIgetVoxVal(mri_save, x, y, z, 0)) ;
        oval = nint(MRIgetVoxVal(mri_out, x, y, z, 0)) ;
        if ((oval == target_label && sval == 0) ||
            (oval != target_label && sval > 0))
          DiagBreak() ;
      }
  return(NO_ERROR) ;
}
#endif

const char *Progname;

static int resegment_erased_voxels(MRI *mri_T1, MRI *mri_in, MRI *mri_out, int label) ;
static int build_label_histograms(MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histos) ;
static MRI_TOPOLOGY_PARMS parms;

static void Error(const char *string);
static int get_option(int argc, char *argv[]);

//--------------------------------------------
/*Error routine*/
static void Error(const char *string) {
  fprintf(stderr, "\nError %s\n",string) ;
  exit(0) ;
}

/*get_option routine*/
static int get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  if (!strcmp(option, "nothing")) {
    fprintf(stderr,"Mode:          NOTHING\n") ;
    nargs = 0 ;
  } else if (!strcmp(option, "final_surf")) {
    parms.final_surface_file=argv[2];
    fprintf(stderr,"Mode:          Writing out final surface into %s\n",parms.final_surface_file) ;
    nargs = 1 ;
  } else if (!strcmp(option, "verbose")) {
    parms.verbose_mode=1;
    fprintf(stderr,"Mode:          verbose mode on\n") ;
    nargs = 0 ;
  } else if (!strcmp(option, "initial_surf")) {
    parms.initial_surface_file=argv[2];
    fprintf(stderr,"Mode:          Writing out initial surface into %s\n",parms.initial_surface_file) ;
    nargs = 1 ;
  } else if (!strcmp(option, "foreground")) {
    parms.background_priority=0;
    fprintf(stderr,"Mode:          Making corrections on the background first\n") ;
    nargs = 0 ;
  } else if (!strcmp(option, "background")) {
    parms.background_priority=1;
    fprintf(stderr,"Mode:          Making corrections on the foreground first\n") ;
    nargs = 0 ;
  } else if (!strcmp(option, "only")) {
    parms.only=1;
    fprintf(stderr,"Mode:          Making corrections on the foreground/background only\n") ;
    nargs = 0 ;
  } else if (!strcmp(option, "tess")) {
    parms.tesselation_mode=atoi(argv[2]);
    fprintf(stderr,"Mode:          Tesselating with mode %d\n",parms.tesselation_mode) ;
    nargs = 1 ;
  } else if (!strcmp(option, "priors")) {
    parms.using_gca_maps=1;
    parms.transform_fname=argv[2];
    parms.gca_fname=argv[3];
    fprintf(stderr,"Mode:          Using gca information: "
            "\n               transform %s"
            "\n               gca %s\n"
            ,parms.transform_fname,parms.gca_fname) ;
    nargs = 2 ;
  } else if (!strcmp(option, "priormap")) {
    parms.prior_map_file=argv[2];
    fprintf(stderr,"Mode:          Loading prior map from %s\n",parms.prior_map_file) ;
    nargs = 1 ;
  } else if (!strcmp(option, "VOXEL")) {
    parms.mode=VOXEL_MODE;
    fprintf(stderr,"Mode:          Cost computed from the number of voxels\n") ;
    nargs = 0 ;
  } else if (!strcmp(option, "MAP")) {
    parms.mode=MAP_MODE;
    fprintf(stderr,"Mode:          Cost computed from map estimate\n") ;
    nargs = 0 ;
  } else if (!strcmp(option, "PROB")) {
    parms.mode=PROB_MODE;
    fprintf(stderr,"Mode:          Cost computed from the sum of the voxel prob\n") ;
    nargs = 0 ;
  } else if (!strcmp(option, "PROB_MAP")) {
    parms.mode=PROB_MAP_MODE;
    fprintf(stderr,"Mode:          Cost computed from map estimate\n") ;
    nargs = 0 ;
  } else if (!strcmp(option, "guess")) {
    parms.guess_initial_segmentation=1;
    fprintf(stderr,"Mode:          Guess the initial Segmentation\n") ;
    nargs = 0 ;
  } else if (!strcmp(option, "maps")) {
    parms.debugging_map_folder=argv[2];
    fprintf(stderr,"Mode:          Writing out the maps into folder %s\n",parms.debugging_map_folder) ;
    nargs = 1 ;
  } else if (!strcmp(option, "connectivity")) {
    parms.connectivity=atoi(argv[2]);
    fprintf(stderr,"Mode:          Connectivity %d\n",parms.connectivity) ;
    nargs = 1 ;
  } else if (!strcmp(option, "label")) {
    parms.labels[parms.nlabels]=atoi(argv[2]);
    parms.nlabels++;
    fprintf(stderr,"Mode:          Correction topology for label %d\n",atoi(argv[2])) ;
    nargs = 1 ;
  } else if (!strcmp(option, "beta")) {
    parms.beta=MAX(0.0,MIN(1.0,atof(argv[2])));
    fprintf(stderr,"Mode:          Mixing parameters set to %3.3f\n",parms.beta) ;
    nargs = 1 ;
  } else if (!strcmp(option, "alpha")) {
    parms.alpha=MAX(0.0,MIN(1.0,atof(argv[2])));
    fprintf(stderr,"Mode:          Mixing parameters set to %3.3f\n",parms.alpha) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'L':
      parms.labels[parms.nlabels]=atoi(argv[2]);
      parms.nlabels++;
      fprintf(stderr,"Mode:          Correction topology for label %d\n",atoi(argv[2])) ;
      nargs = 1 ;
      break;
    default:
      printf("Mode:          unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}


static void MRI_TOPOLOGY_PARMSdefault(MRI_TOPOLOGY_PARMS *parms) {
  parms->connectivity=1;
  parms->tesselation_mode=-1;
  parms->MarchingCubes=0;
  parms->nlabels=0;
  parms->mode=VOXEL_MODE;
  parms->using_gca_maps=0;
  parms->transform_fname=0;
  parms->gca_fname=0;
  parms->initial_surface_file=0;
  parms->final_surface_file=0;
  parms->debugging_map_folder=0;
  parms->background_priority=0;
  parms->only=0;
  parms->prior_map_file=0;
  parms->alpha=1.0f;
  parms->beta=1.0f;
  parms->guess_initial_segmentation=0;
  parms->generate_surface=0;
  parms->verbose_mode=0;
  parms->gca=0;
  parms->transform=0;
}

int main(int argc, char *argv[]) {
  MRIS *mris;
  char  *in_orig_fname=NULL, *in_seg_fname=NULL,*out_fname=NULL;
  MRI *mri_orig=NULL,*mri_seg=NULL,*mri_out=NULL;
  int nargs,n;
  char fname[512];

  Progname=argv[0];
  fprintf(stderr,"\n");
  MRI_TOPOLOGY_PARMSdefault(&parms);

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  if (parms.tesselation_mode==-1)
    parms.tesselation_mode=parms.connectivity;
  if (argc<4) {
    fprintf(stderr, "\nUsage: %s options input_orig_file input_segmented_file output_folder\n", Progname);
    exit(1);
  };

  in_orig_fname=argv[argc-3];
  in_seg_fname = argv[argc-2];
  out_fname = argv[argc-1];

  fprintf(stderr,"************************************************************"
          "\nThe input orig volume is %s"
          "\nThe input segmented volume is %s"
          "\nThe output volume is %s"
          "\nIf this is incorrect, please exit quickly the program (Ctl-C)\n",in_orig_fname,in_seg_fname,out_fname);
  for (n=0;n<parms.nlabels;n++)
    fprintf(stderr,"label = %d: %s \n",parms.labels[n],cma_label_to_name(parms.labels[n]));
  if (parms.using_gca_maps)
    fprintf(stderr,"mixing parameters: alpha=%1.3f , beta=%1.3f \n",parms.alpha,parms.beta);
  else {
    parms.beta=1.0f;
    parms.alpha=1.0f;
  }
  fprintf(stderr,"connectivity = %d\n",parms.connectivity);

  mri_orig=MRIread(in_orig_fname);
  if (!mri_orig && parms.using_gca_maps)
    Error("orig volume: orig volume could not be read\n");
  mri_seg=MRIread(in_seg_fname);
  if (!mri_seg)
    Error("segmented volume: segmented volume could not be read\n");


  //check euler characteristic of initial surface
  if (parms.initial_surface_file) {
    int i,j,k,val,euler,pnvertices,  pnfaces, pnedges;
    MRI *mri_tmp;
    mri_tmp=MRIclone(mri_seg,NULL);
    for (k=0;k<mri_seg->depth;k++)
      for (j=0;j<mri_seg->height;j++)
        for (i=0;i<mri_seg->width;i++)
          for (n=0;n<parms.nlabels;n++) {
            val=MRIgetVoxVal(mri_seg,i,j,k, 0);
            if (val==parms.labels[n]) {
              MRIsetVoxVal(mri_tmp,i,j,k,0,1);
              break;
            }
          }
    mris=MRIScreateSurfaceFromVolume(mri_tmp,1,parms.connectivity);
    euler=MRIScomputeEulerNumber(mris,&pnvertices,&pnfaces,&pnedges);
    fprintf(stderr,"\ninitial euler characteristic = %d, %d vertices, %d faces, %d edges"
            ,euler,pnvertices,pnfaces,pnedges);
    MRISwrite(mris,parms.initial_surface_file);
    MRISfree(&mris);
    MRIfree(&mri_tmp);
  }

  mri_out=MRIcorrectTopology(mri_orig,mri_seg,NULL,&parms);

  if (parms.nlabels == 1) {
    MRI *mri_tmp ;

    // turn off all voxels that are going to be on in the output
    MRImask(mri_seg, mri_out, mri_seg, 1, 0) ;
    /* whatever ones are left are now incorrect and should be labeled
      something else
    */
    resegment_erased_voxels(mri_orig, mri_seg, mri_seg, parms.labels[0]) ;
    MRIreplaceValues(mri_out, mri_out, 1, parms.labels[0]) ;
    mri_tmp = MRIcopy(mri_seg, NULL) ;
    MRIcopyLabel(mri_out, mri_tmp, parms.labels[0]) ;
    MRIfree(&mri_out) ;
    mri_out = mri_tmp ;
    //  check_volume(mri_save, mri_out, parms.labels[0]) ;
  }
  MRIwrite(mri_out,out_fname);

  ////TEMPORARY VALIDATION STUFF //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
#if 0
  //validation of the algo
  {
    FILE *f;
    MRIS *mristb[20],*mrisr;
    int n,i,j,k,depth,height,width,count,count2;
    int tab[20]={4,43,51,12,52,13,54,18,53,17,49,10,50,11};//,6,7,10,11,12,13,17,18,43,44,45,46,49,50,51,52,53,54};
    MRI *mri_val=MRIclone(parms.mri_seg,NULL);
    parms.nlabels=1;

    depth=parms.mri_seg->depth;
    height=parms.mri_seg->height;
    width=parms.mri_seg->width;
    for (n=0;n<14;n++) {
      MRIfree(&parms.mri_output);
      MRIfree(&parms.mri_bin);
      MRIfree(&parms.mri_dist);
      MRIfree(&parms.mri_fcost);
      MRIfree(&parms.mri_bcost);
      MRIfree(&parms.mri_fprior);
      MRIfree(&parms.mri_bprior);
      MRIfree(&parms.mri_labeled);
      segmentationFree(&parms.F_Bseg);
      segmentationFree(&parms.F_Rseg);
      segmentationFree(&parms.B_Bseg);
      segmentationFree(&parms.B_Rseg);
      CCSfree(&parms.F_Bccs);
      CCSfree(&parms.F_Rccs);
      CCSfree(&parms.B_Bccs);
      CCSfree(&parms.B_Rccs);

      parms.labels[0]=tab[n];
      MRIcorrectTopology(parms.mri_orig,parms.mri_seg,&parms.mri_output,mris
                         ,parms.labels,parms.nblabels,parms.f_c,parms);



      MRISwrite(*mris,"./tmp");
      mristb[n]=MRISread("./tmp");
#if 0
      count=0;
      count2=0;
      for (k=0;k<depth;k++)
        for (j=0;j<height;j++)
          for (i=0;i<width;i++) {
            if (MRIvox(parms.mri_seg,i,j,k)==parms.labels[0])
              count2++;
            if (MRIvox(parms.mri_output,i,j,k)==1) {
              MRIvox(mri_val,i,j,k)++;
              if (MRIvox(parms.mri_seg,i,j,k)!=parms.labels[0])
                count++;
            } else if (MRIvox(parms.mri_seg,i,j,k)==parms.labels[0])
              count++;
          }
      fprintf(stderr,"\n yeh %d %d %f \n",count,count2,100.*count/count2);
      sprintf(fname,"./label%d",tab[n]);
      f=fopen(fname,"a+");
      fprintf(f,"\n %d %d %f ",count,count2,(float)100.*count/count2);
      fclose(f);
#endif

#if 0
      sprintf(fname,"./surf%d",n);
      MRISwrite(mristb[n],fname);
      MRISsmoothSurface2(mristb[n],5,0.5,0);
      MRISsmoothSurface2(mristb[n],5,0.25,2);
      MRISsmoothSurface2(mristb[n],10,0.05,5);
      sprintf(fname,"./surfsmooth%d",n);
      mristb[n]->type=MRIS_TRIANGULAR_SURFACE;//MRIS_BINARY_QUADRANGLE_FILE;
      MRISwrite(mristb[n],fname);

      MRISsetNeighborhoodSizeAndDist(mristb[n],3) ;
      MRIScomputeMetricProperties(mristb[n]) ;
      MRIScomputeSecondFundamentalForm(mristb[n]) ;
      MRISuseMeanCurvature(mristb[n]);
      MRISaverageCurvatures(mristb[n],2) ;
      MRISnormalizeCurvature(mristb[n], NORM_MEAN) ;
      sprintf(fname,"./curv%d",n);
      MRISwriteCurvature(mristb[n],fname);
#endif
    }

#if 0
    mrisr=MRISconcatenateQuadSurfaces(n,mristb);
    mrisr->type=MRIS_TRIANGULAR_SURFACE;
    MRISwrite(mrisr,"./lh.ZURFACE");


    //    for(k=0;k<mrisr->nvertices;k++)
    // mrisr->vertices[k].curv=0.3;

    //MRISnormalizeCurvature(mrisr, NORM_MEAN) ;
    MRISwriteCurvature(mrisr,"./ZURFACE_CURVATURE");
    for (k=0;k<mrisr->nvertices;k++)
      mrisr->vertices[k].curv=mrisr->vertices[k].val;
    MRISwriteCurvature(mrisr,"./ZURFACE_VAL");
#endif

    n=0;
    count=0;
    for (k=0;k<depth;k++)
      for (j=0;j<height;j++)
        for (i=0;i<width;i++) {
          if (MRIgetVoxVal(mri_val,i,j,k,0)>=1) {
            n++;
            if (MRIsetVoxVal(mri_val,i,j,k,0)>1)
              count++;
          }
        }
    //    sprintf(fname,"./labeltotal");
    /// f=fopen(fname,"a+");
    //fprintf(f,"\n %s %d %d %f ",in_seg_fname,count,n,(float)100.*count/n);
    //fclose(f);




#if 0
    MRIwrite(mri_val,"/tmp/tmp");
#endif

    fprintf(stderr,"\n WE HAVE %d %d %f   \n",count,n,100.*count/n);

  }
#endif
  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  if (parms.final_surface_file) {
    int euler,pnvertices,  pnfaces, pnedges;
    mris=MRIScreateSurfaceFromVolume(mri_out,1,parms.connectivity);
    euler=MRIScomputeEulerNumber(mris,&pnvertices,&pnfaces,&pnedges);
    fprintf(stderr,"\nfinal euler characteristic = %d, %d vertices, %d faces, %d edges"
            ,euler,pnvertices,pnfaces,pnedges);
    sprintf(fname,"%s",parms.final_surface_file);
    MRISwrite(mris,fname);

#if 0
    MRISsmoothSurface(mris,7,0.2);
    strcat(fname,"_smooth");
    MRISwrite(mris,fname);
    if (parms.fit) {
      sprintf(fname,parms.surfname);
      strcat(fname,"_fit");
      MRISmatchSurfaceToLabel(parms.mris,parms.mri_output,1,NULL,NULL,parms.f_c);
      MRISwrite(parms.mris,fname);
    }
#endif
    MRISfree(&mris);
  }

  if (mri_out)
    MRIfree(&mri_out);
  if (mri_orig)
    MRIfree(&mri_orig);
  if (mri_seg)
    MRIfree(&mri_seg);
  fprintf(stderr,"\n");
  return NO_ERROR;
}

/*
 figure out what to do with voxels that were turned 'off' by the
 topology correction. This is a hack, but for now just make them
 the most likely of the nbr voxel labels.
*/
static int
resegment_erased_voxels(MRI *mri_T1, MRI *mri_in, MRI *mri_out, int target_label) {
  int       x, y, z, label_in, label_out, xi, yi, zi, xk, yk, zk, label, changed=0 ;
  HISTOGRAM *histos[MAX_CMA_LABEL+1] ;
  double    p, max_p, val ;

  build_label_histograms(mri_in, mri_T1, histos) ;
  for (x = 0 ; x < mri_in->width ; x++) {
    for (y = 0 ; y < mri_in->height ; y++) {
      for (z = 0 ; z < mri_in->depth ; z++) {
        label_in = nint(MRIgetVoxVal(mri_in, x, y, z, 0)) ;
        label_out = nint(MRIgetVoxVal(mri_out, x, y, z, 0)) ;
        if (label_in == target_label) {
          // find most likely nbr label
          max_p = 0 ;
          label_out = label_in ;
          for (xk = -1 ; xk <= 1 ; xk++) {
            xi = x + xk ;
            if (xi < 0 || xi >= mri_in->width)
              continue ;
            for (yk = -1 ; yk <= 1 ; yk++) {
              yi = y + yk ;
              if (yi < 0 || yi >= mri_in->height)
                continue ;
              for (zk = -1 ; zk <= 1 ; zk++) {
                zi = z + zk ;
                if (zi < 0 || zi >= mri_in->depth)
                  continue ;
                label = nint(MRIgetVoxVal(mri_in, xi, yi, zi, 0)) ;
                if (label == label_in)
                  continue ;  // would be topologically incorrect
                val = MRIgetVoxVal(mri_T1, xi, yi, zi, 0) ;
                p = HISTOvalToCount(histos[label], val) ;
                if (p > max_p) {
                  max_p = p ;
                  label_out = label ;
                }
              }
            }
          }
          changed++ ;
          MRIsetVoxVal(mri_out, x, y, z, 0, label_out) ;
        }
      }
    }
  }
  printf("%d voxels resegmented to be ML\n", changed) ;
  return(NO_ERROR) ;
}

static int
build_label_histograms(MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histos) {
  int    labels[MAX_LABEL+1], x, y, z, l ;
  MRI    *mri_xformed ;
  MATRIX *m_vox2vox ;

  memset(labels, 0, sizeof(labels)) ;

  m_vox2vox = MRIgetVoxelToVoxelXform(mri_intensities, mri_labels) ;
  mri_xformed = MRIclone(mri_labels, NULL) ;
  MRIlinearTransform(mri_intensities, mri_xformed, m_vox2vox) ;
  MatrixFree(&m_vox2vox) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_xformed, "x.mgz") ;

  for (x = 0 ; x < mri_labels->width ; x++)
    for (y = 0 ; y < mri_labels->height ; y++)
      for (z = 0 ; z < mri_labels->depth ; z++) {
        l = nint(MRIgetVoxVal(mri_labels, x, y, z, 0)) ;
        //    if (l == 0)
        //     continue ;
        if (labels[l] == 0)  // first time
        {
          char fname[STRLEN] ;
          histos[l] = MRIhistogramLabel(mri_xformed, mri_labels, l, 50) ;
          HISTOmakePDF(histos[l], histos[l]) ;
          sprintf(fname, "label%d.plt", l) ;
          if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
            HISTOplot(histos[l], fname) ;
        }

        labels[l] = 1 ;
      }

  MRIfree(&mri_xformed) ;
  return(NO_ERROR) ;
}
