/**
 * @brief program to add a template into an average surface
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[],INTEGRATION_PARMS *parms) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static const char *surf_dir = "surf" ;
static char *annot_name = NULL ;
static const char *surface_names[] =
  {
    "inflated",
    "smoothwm",
    "smoothwm"
  } ;

static const char *curvature_names[] =
  {
    "inflated.H",
    "sulc",
    NULL
  } ;

#define IMAGES_PER_SURFACE   3   /* mean, variance, and dof */
#define SURFACES         sizeof(curvature_names) / sizeof(curvature_names[0])
#define PARAM_IMAGES         (IMAGES_PER_SURFACE*SURFACES)

static int nbrs = 3 ;
static int navgs = 0 ;
static float scale = 1 ;
static int no_rot = 1 ;
static char subjects_dir[STRLEN] ;

/* VECTORIAL_REGISTRATION */
static void setParms(INTEGRATION_PARMS *parms);
static int multiframes = 0 ; /* to use multiframes */
static int base_default = 1 ;  /* using default vector fields */
static int atlas_size = 3;

#define MAX_OVERLAYS 1000
static int noverlays = 0 ;
static char *overlays[MAX_OVERLAYS] ;
static const char *overlay_dir = "label";

static int which_norm = NORM_MEAN;

int
main(int argc, char *argv[])
{
  char         **av, surf_fname[STRLEN], *template_fname, *hemi, *sphere_name,
  *cp, *subject, fname[STRLEN] ;
  int          ac, nargs, ino, sno, nbad = 0, failed, n,nfields;
  VERTEX *v;
  VALS_VP *vp;
  MRI_SURFACE  *mris ;
  MRI_SP       *mrisp, /* *mrisp_aligned,*/ *mrisp_template ;
  INTEGRATION_PARMS parms ;

  nargs = handleVersionOption(argc, argv, "mris_make_template");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  /* setting default values for vectorial registration */
  setParms(&parms);

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv,&parms) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5) usage_exit() ;

  /* multiframe registration */
  if (multiframes) parms.flags |= IP_USE_MULTIFRAMES;

  if (!strlen(subjects_dir))  /* not specified on command line*/
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n",
                Progname) ;
    strcpy(subjects_dir, cp) ;
  }
  hemi = argv[1] ;
  sphere_name = argv[2] ;
  template_fname = argv[argc-1] ;
  if (1 || !FileExists(template_fname))  /* first time - create it */
  {
    fprintf(stderr, "creating new parameterization...\n") ;
    if (multiframes)
    {
      mrisp_template = MRISPalloc(scale, atlas_size * IMAGES_PER_SURFACE );
      /*    if (no_rot)  /\* don't do rigid alignment *\/ */
      /*     mrisp_aligned = NULL ; */
      /*    else */
      /*     mrisp_aligned = MRISPalloc(scale, PARAM_FRAMES);  */
    }
    else
    {
      mrisp_template = MRISPalloc(scale, PARAM_IMAGES);
      /*    if (no_rot)  /\* don't do rigid alignment *\/ */
      /*     mrisp_aligned = NULL ; */
      /*    else */
      /*     mrisp_aligned = MRISPalloc(scale, PARAM_IMAGES);  */
    }

  }
  else
  {
    fprintf(stderr, "reading template parameterization from %s...\n",
            template_fname) ;
    /* mrisp_aligned = NULL ; */
    mrisp_template = MRISPread(template_fname) ;
    if (!mrisp_template)
      ErrorExit(ERROR_NOFILE, "%s: could not open template file %s",
                Progname, template_fname) ;
  }

  argv += 3 ;
  argc -= 3 ;
  for (ino = 0 ; ino < argc-1 ; ino++)
  {
    failed = 0 ;
    subject = argv[ino] ;
    fprintf(stderr, "\nprocessing subject %s (%d of %d)\n", subject,
            ino+1, argc-1) ;
    int req = snprintf(surf_fname, STRLEN, "%s/%s/%s/%s.%s",
		       subjects_dir, subject, surf_dir, hemi, sphere_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stderr, "reading spherical surface %s...\n", surf_fname) ;
    mris = MRISread(surf_fname) ;
    if (!mris)
    {
      nbad++ ;
      ErrorPrintf(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, surf_fname) ;
      exit(1) ;
    }
    if (annot_name)
    {
      if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
        ErrorExit(ERROR_BADPARM,
                  "%s: could not read annot file %s",
                  Progname, annot_name) ;
      MRISripMedialWall(mris) ;
    }

    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    MRISstoreMetricProperties(mris) ;

    if (Gdiag & DIAG_WRITE)
    {
      char *cp1 ;

      FileNameOnly(template_fname, fname) ;
      cp = strchr(fname, '.') ;
      if (cp)
      {
        cp1 = strrchr(fname, '.') ;
        if (cp1 && cp1 != cp)
          strncpy(parms.base_name, cp+1, cp1-cp-1) ;
        else
          strcpy(parms.base_name, cp+1) ;
      }
      else
        strcpy(parms.base_name, "template") ;
      int req = snprintf(fname, STRLEN, "%s.%s.out", hemi, parms.base_name);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      INTEGRATION_PARMS_openFp(&parms, fname, "w") ;
      printf("writing output to '%s'\n", fname) ;
    }

    /* multiframe registration */
    if (multiframes)
    {
      nfields=parms.nfields;

      for ( n = 0; n < mris->nvertices ; n++) /* allocate the VALS_VP
                                                                 structure */
      {
        v=&mris->vertices[n];
        vp=(VALS_VP *)calloc(1,sizeof(VALS_VP));
        vp->nvals=nfields;
        vp->orig_vals=(float*)malloc(nfields*sizeof(float)); /* before
                                                                blurring */
        vp->vals=(float*)malloc(nfields*sizeof(float));     /* values used by
                                                               MRISintegrate */
        v->vp=(void*)vp;
      }

      /* load the different fields */
      for (n = 0 ; n < parms.nfields ; n++)
      {
        if (parms.fields[n].name != NULL)
        {
          int req = snprintf(surf_fname, STRLEN, "%s/%s/%s/%s.%s", subjects_dir,
			     subject, overlay_dir, hemi, parms.fields[n].name) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          printf("reading overlay file %s...\n", surf_fname) ;
          if (MRISreadValues(mris, surf_fname) != NO_ERROR)
            ErrorExit(ERROR_BADPARM, "%s: could not read overlay file %s",
                      Progname, surf_fname) ;
          MRIScopyValuesToCurvature(mris) ;
        }
        else if (ReturnFieldName(parms.fields[n].field))
        {
          /* read in precomputed curvature file */
          int req = snprintf(surf_fname, STRLEN, "%s/%s/%s/%s.%s", subjects_dir,
			     subject, surf_dir, hemi, ReturnFieldName(parms.fields[n].field)) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          // fprintf(stderr,"\nreading field %d from %s(type=%d,frame=%d)\n",parms.fields[n].field,surf_fname,parms.fields[n].type,parms.fields[n].frame);
          if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR)
          {
            fprintf(stderr,"\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
            fprintf(stderr, "%s: could not read curvature file '%s'\n",
                    Progname, surf_fname) ;
            failed = 1;
            break;
          }
        }
        else
        {                       /* compute curvature of surface */
          int req = snprintf(surf_fname, STRLEN, "%s/%s/%s/%s.%s", subjects_dir,
			     subject, surf_dir, hemi, surface_names[parms.fields[n].field]) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          /*if(parms.fields[n].field==0)
           sprintf(fname, "inflated") ;
           else
           sprintf(fname, "smoothwm") ;*/
          //fprintf(stderr,"\ngenerating field %d(type=%d,frame=%d) (from %s)\n",parms.fields[n].field,parms.fields[n].type,parms.fields[n].frame,surf_fname);
          //     MRISsaveVertexPositions(mris, TMP_VERTICES) ;
          if (MRISreadVertexPositions(mris, surf_fname) != NO_ERROR)
          {
            fprintf(stderr,"\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
            ErrorPrintf(ERROR_NOFILE, "%s: could not read surface file %s",
                        Progname, surf_fname) ;
            fprintf(stderr,"setting up correlation coefficient to zero\n");
            parms.fields[n].l_corr=parms.fields[n].l_pcorr=0.0;
            failed=1;
            break;
          }

          if (nbrs > 1) MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
          MRIScomputeMetricProperties(mris) ;
          MRIScomputeSecondFundamentalForm(mris) ;
          MRISuseMeanCurvature(mris) ;
          MRISaverageCurvatures(mris, navgs) ;
          MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
        }
        /*    if(parms.fields[n].field!=SULC_CORR_FRAME)*/
        MRISnormalizeField(mris,parms.fields[n].type,
                           parms.fields[n].which_norm); /* normalize values */
        MRISsetCurvaturesToOrigValues(mris,n);
        MRISsetCurvaturesToValues(mris,n);
      }

      if (failed)
      {
        fprintf(stderr,"\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
        fprintf(stderr,"Subject %s Failed",subject);
        fprintf(stderr,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n");
        /* free cal structure*/
        for ( n = 0; n < mris->nvertices ; n++)
        {
          v=&mris->vertices[n];
          vp=(VALS_VP*)v->vp;
          free(vp->orig_vals);
          free(vp->vals);
          free(vp);
          v->vp=NULL;
        }
        /* free surface */
        MRISfree(&mris);
        /* go onto the next subject */
        continue;
      }
    }

    if (multiframes && (!no_rot))
    { /* rigid body alignment */
      parms.frame_no = 3 ;  /* don't use single field correlation functions */
      parms.l_corr = parms.l_pcorr = 0.0f ;

      parms.mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
      parms.mrisp_template = mrisp_template ;

      MRISrigidBodyAlignVectorGlobal(mris, &parms, 1.0, 64.0, 8) ;
      if (Gdiag & DIAG_WRITE) MRISwrite(mris, "sphere.rot.global") ;
      MRISrigidBodyAlignVectorLocal(mris, &parms) ;
      if (Gdiag & DIAG_WRITE) MRISwrite(mris, "sphere.rot.local") ;
      MRISPfree(&parms.mrisp) ;
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    };
    if ((!multiframes) && (!no_rot) && ino > 0)
    { /* rigid body alignment */
      int req = snprintf(surf_fname, STRLEN, "%s/%s/%s/%s.%s",
			 subjects_dir, subject, surf_dir, hemi, "sulc") ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR)
      {
        ErrorPrintf(Gerror, "%s: could not read curvature file '%s'\n",
                    Progname, surf_fname) ;
        nbad++ ;
        MRISfree(&mris) ;
        continue ;
      }
      parms.frame_no = 3 ; /* use sulc for rigid registration */
      parms.mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
      parms.mrisp_template = mrisp_template ;
      parms.l_corr = 1.0f ;

      MRISrigidBodyAlignGlobal(mris, &parms, 1.0, 64.0, 8) ;
      if (Gdiag & DIAG_WRITE)
        MRISwrite(mris, "sphere.rot.global") ;
      MRISrigidBodyAlignLocal(mris, &parms) ;
      if (Gdiag & DIAG_WRITE)
        MRISwrite(mris, "sphere.rot.local") ;
      MRISPfree(&parms.mrisp) ;
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    }

    if (multiframes)
    {
      for (n = 0; n < parms.nfields ; n++)
      {
        MRISsetOrigValuesToCurvatures(mris,n);
        MRISaverageCurvatures(mris, parms.fields[n].navgs) ;
        mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
        MRISPcombine(mrisp,
                     mrisp_template,
                     parms.fields[n].frame * IMAGES_PER_SURFACE) ;
        MRISPfree(&mrisp) ;
      }
      /* free the VALS_VP structure */
      for ( n = 0; n < mris->nvertices ; n++)
      {
        v=&mris->vertices[n];
        vp=(VALS_VP*)v->vp;
        free(vp->orig_vals);
        free(vp->vals);
        free(vp);
        v->vp=NULL;
      }
      MRISfree(&mris) ;
    }
    else
    {
      for (sno = 0; sno < SURFACES ; sno++)
      {
        if (curvature_names[sno])  /* read in precomputed curvature file */
        {
          int req = snprintf(surf_fname, STRLEN, "%s/%s/%s/%s.%s",
			     subjects_dir, subject, surf_dir, hemi, curvature_names[sno]) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR)
          {
            nbad++ ;
            ErrorPrintf(Gerror, "%s: could not read curvature file '%s'\n",
                        Progname, surf_fname) ;
            failed = 1 ;
            break ;
          }
          /* the two next lines were not in the original code */
          MRISaverageCurvatures(mris, navgs) ;
          MRISnormalizeCurvature(mris, which_norm) ;
        } else                       /* compute curvature of surface */
        {
          int req = snprintf(surf_fname, STRLEN, "%s/%s/%s/%s.%s",
			     subjects_dir, subject, surf_dir, hemi, surface_names[sno]) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          if (MRISreadVertexPositions(mris, surf_fname) != NO_ERROR)
          {
            ErrorPrintf(ERROR_NOFILE, "%s: could not read surface file %s",
                        Progname, surf_fname) ;
            nbad++ ;
            failed = 1 ;
            break ;
          }

          if (nbrs > 1)
            MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
          MRIScomputeMetricProperties(mris) ;
          MRIScomputeSecondFundamentalForm(mris) ;
          MRISuseMeanCurvature(mris) ;
          MRISaverageCurvatures(mris, navgs) ;
          MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
          MRISnormalizeCurvature(mris, which_norm) ;
        }
        fprintf(stderr, "computing parameterization for surface %s...\n",
                surf_fname);
        if (failed)
        {
          continue ;
          MRISfree(&mris) ;
        }
        mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
        MRISPcombine(mrisp, mrisp_template, sno*3) ;
        MRISPfree(&mrisp) ;
      }
      MRISfree(&mris) ;
    }
  }

#if 0
  if (mrisp_aligned)  /* new parameterization - use rigid alignment */
  {
    MRI_SP *mrisp_tmp ;

    if (Gdiag & DIAG_WRITE)
    {
      char *cp1 ;

      FileNameOnly(template_fname, fname) ;
      cp = strchr(fname, '.') ;
      if (cp)
      {
        cp1 = strrchr(fname, '.') ;
        if (cp1 && cp1 != cp)
          strncpy(parms.base_name, cp+1, cp1-cp-1) ;
        else
          strcpy(parms.base_name, cp+1) ;
      }
      else
        strcpy(parms.base_name, "template") ;
      sprintf(fname, "%s.%s.out", hemi, parms.base_name);
      parms.fp = fopen(fname, "w") ;
      printf("writing output to '%s'\n", fname) ;
    }
    for (ino = 0 ; ino < argc-1 ; ino++)
    {
      subject = argv[ino] ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms.fp, "processing subject %s\n", subject) ;
      fprintf(stderr, "processing subject %s\n", subject) ;
      sprintf(surf_fname, "%s/%s/%s/%s.%s",
              subjects_dir, subject, surf_dir, hemi, sphere_name) ;
      fprintf(stderr, "reading spherical surface %s...\n", surf_fname) ;
      mris = MRISread(surf_fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, surf_fname) ;
      MRIScomputeMetricProperties(mris) ;
      MRISstoreMetricProperties(mris) ;
      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
      sprintf(surf_fname, "%s/%s/%s/%s.%s",
              subjects_dir, subject, surf_dir, hemi, "sulc") ;
      if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read curvature file '%s'\n",
                  Progname, surf_fname) ;
      parms.frame_no = 3 ;
      parms.mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
      parms.mrisp_template = mrisp_template ;
      parms.l_corr = 1.0f ;

      MRISrigidBodyAlignGlobal(mris, &parms, 1.0, 32.0, 8) ;
      if (Gdiag & DIAG_WRITE)
        MRISwrite(mris, "sphere.rot.global") ;
      MRISrigidBodyAlignLocal(mris, &parms) ;
      if (Gdiag & DIAG_WRITE)
        MRISwrite(mris, "sphere.rot.local") ;
      MRISPfree(&parms.mrisp) ;

#if 0
      /* write out rotated surface */
      sprintf(surf_fname, "%s.rot", mris->fname) ;
      fprintf(stderr, "writing out rigidly aligned surface to '%s'\n",
              surf_fname) ;
      MRISwrite(mris, surf_fname) ;
#endif

      /* now generate new parameterization using the optimal alignment */
      for (sno = 0; sno < SURFACES ; sno++)
      {
        if (curvature_names[sno])  /* read in precomputed curvature file */
        {
          sprintf(surf_fname, "%s/%s/%s/%s.%s",
                  subjects_dir, subject, surf_dir, hemi, curvature_names[sno]) ;
          if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR)
            ErrorExit(Gerror, "%s: could not read curvature file '%s'\n",
                      Progname, surf_fname) ;
        } else                       /* compute curvature of surface */
        {
          sprintf(surf_fname, "%s/%s/%s/%s.%s",
                  subjects_dir, subject, surf_dir, hemi, surface_names[sno]) ;
          if (MRISreadVertexPositions(mris, surf_fname) != NO_ERROR)
            ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                      Progname, surf_fname) ;

          if (nbrs > 1)
            MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
          MRIScomputeMetricProperties(mris) ;
          MRIScomputeSecondFundamentalForm(mris) ;
          MRISuseMeanCurvature(mris) ;
          MRISaverageCurvatures(mris, navgs) ;
          MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
          MRISnormalizeCurvature(mris) ;
        }
        fprintf(stderr, "computing parameterization for surface %s...\n",
                surf_fname);
        mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
        MRISPcombine(mrisp, mrisp_aligned, sno*3) ;
        MRISPfree(&mrisp) ;
      }
      MRISfree(&mris) ;
    }

    if (Gdiag & DIAG_WRITE)
      fclose(parms.fp) ;

    mrisp_tmp = mrisp_aligned ;
    mrisp_aligned = mrisp_template ;
    mrisp_template = mrisp_tmp ;
    MRISPfree(&mrisp_aligned) ;
  }
#endif
  fprintf(stderr,
          "writing updated template with %d subjects to %s...\n",
          argc-1-nbad, template_fname) ;
  MRISPwrite(mrisp_template, template_fname) ;
  MRISPfree(&mrisp_template) ;
  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[],INTEGRATION_PARMS *parms)
{
  int  n,nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "help") || !stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "version") || !stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size = %d\n", nbrs) ;
  }
  else if (!stricmp(option, "sdir"))
  {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using SUBJECTS_DIR=%s\n", subjects_dir) ;
  }
  else if (!stricmp(option, "median"))
  {
    which_norm = NORM_MEDIAN ;
    printf("using median normalization\n") ;
  }
  else if (!stricmp(option, "nonorm"))
  {
    which_norm = NORM_NONE ;
    printf("not normalizing input data\n") ;
  }
  else if (!stricmp(option, "infname"))
  {
    char fname[STRLEN] ;
    surface_names[0] = argv[2] ;
    nargs = 1 ;
    printf("using %s instead of inflated\n", argv[2]) ;
    sprintf(fname, "%s.H", argv[2]) ;
    curvature_names[0] = (char *)calloc(strlen(fname)+1, sizeof(char)) ;
    strcpy(const_cast<char*>(curvature_names[0]), fname) ; // const_cast and strcpy
  }
  else if (!stricmp(option, "sulc"))
  {
    nargs = 1 ;
    printf("using %s instead of sulc\n", argv[2]) ;
    curvature_names[1] = (char *)calloc(strlen(argv[2])+1, sizeof(char)) ;
    strcpy(const_cast<char*>(curvature_names[1]), argv[2]) ; // const_cast and strcpy
  }
  else if (!stricmp(option, "norot"))
  {
    no_rot = 1 ;
    fprintf(stderr, "not aligning hemispheres before averaging.\n") ;
  }
  else if (!stricmp(option, "rot"))
  {
    no_rot = 0 ;
    fprintf(stderr, "rigidly aligning hemispheres before averaging.\n") ;
  }
  else if (!stricmp(option, "nodefault"))
  {
    base_default=0;
    atlas_size = 0;
    parms->nfields=0;
    for (n=0;n<3;n++)
      InitFieldLabel(&parms->fields[n]);
  }
  else if (!stricmp(option, "size"))
  {
    int size=atoi(argv[2]);
    if (size<=0)
    {
      fprintf(stderr,"ERROR: size should at least 1!\n");
      exit(1);
    }
    atlas_size = size;
    fprintf(stderr,"Setting up atlas size to %d\n",size);
    nargs=1;
  }
  else if (!stricmp(option, "annot"))
  {
    annot_name = argv[2] ;
    fprintf(stderr,"zeroing medial wall in %s\n", annot_name) ;
    nargs=1;
  }
  else if (!stricmp(option, "surf_dir"))
  {
    surf_dir = argv[2] ;
    fprintf(stderr,"using %s subdirectory, instead of 'surf'\n", surf_dir) ;
    nargs=1;
  }
  else if (!strcmp(option, "overlay"))
  {
    int navgs ;
    overlays[noverlays++] = argv[2] ;
    navgs = atoi(argv[3]) ;
    printf("reading overlay from %s...\n", argv[2]) ;
    multiframes = 1 ;
    n=parms->nfields++;
    SetFieldLabel(&parms->fields[n],
                  OVERLAY_FRAME,
                  atlas_size,
                  0.0,
                  0.0,
                  navgs,
                  which_norm);
    SetFieldName(&parms->fields[n], argv[2]) ;
    atlas_size++ ;
    nargs = 2 ;
  }
  else if (!strcmp(option, "distance"))
  {
    int navgs ;
    overlays[noverlays++] = argv[2] ;
    navgs = atoi(argv[3]) ;
    printf("reading overlay from %s...\n", argv[2]) ;
    multiframes = 1 ;
    n=parms->nfields++;
    SetFieldLabel(&parms->fields[n],
                  DISTANCE_TRANSFORM_FRAME,
                  atlas_size,
                  0.0,
                  0.0,
                  navgs,
                  NORM_MAX);
    SetFieldName(&parms->fields[n], argv[2]) ;
    atlas_size++ ;
    nargs = 2 ;
  }
  else if (!strcmp(option, "overlay-dir"))
  {
    overlay_dir = argv[2] ;
    printf("Setting overlay_dir to %s\n",overlay_dir);
    nargs = 1 ;
  }
  else if (!stricmp(option, "vector"))
  {
    fprintf(stderr,"\nMultiframe Mode:\n");
    fprintf(stderr,
            "Use <addframe> option to add extra-fields into average atlas\n");
    fprintf(stderr,"Use <size> option to set up the # of frames\n");
    fprintf(stderr,"field code:\n");
    for (n=0 ; n < NUMBER_OF_VECTORIAL_FIELDS ; n++)
      fprintf(stderr,"     field %d is '%s' (type = %d)\n",n,
              ReturnFieldName(n),IsDistanceField(n));

  }
  else if (!stricmp(option, "addframe"))
  {
    int which_field,where_in_atlas;

    if (multiframes==0) /* activate multiframes mode */
      multiframes = 1;

    which_field=atoi(argv[2]);
    where_in_atlas=atoi(argv[3]);

    fprintf(stderr,
            "adding field %d into average atlas at location %d\n",
            which_field,where_in_atlas) ;
    if ( base_default && (where_in_atlas<3) )
    {
      fprintf(stderr,"ERROR: Cannot add field %d into atlas\n",which_field);
      fprintf(stderr,"       By default, the first 3 fields are:\n");
      fprintf(stderr,"       -curvature of inflated surface \n");
      fprintf(stderr,"       -sulc of smoothwm surface \n");
      fprintf(stderr,"       -curvature of smoothwm surface \n");
      fprintf(stderr,"To cancel this option, use <nodefault> first\n");
      exit(1);
    }
    /* check if this field exist or not */
    for (n = 0 ; n < parms->nfields ; n++)
    {
      if (parms->fields[n].field==which_field)
      {
        fprintf(stderr,"ERROR: field already exists\n");
        exit(1);
      }
    }
    /* adding field into parms */
    n=parms->nfields++;
    SetFieldLabel(&parms->fields[n],
                  which_field,
                  where_in_atlas,
                  0.0,
                  0.0,
                  0,
                  which_norm);
    if (where_in_atlas >= atlas_size)
    {
      atlas_size = where_in_atlas+1;
      fprintf(stderr,
              "Atlas size increased to contain at least %d frames\n",
              atlas_size);
    }
    nargs = 2 ;
  }
  else switch (toupper(*option))
    {
    case 'W':
      Gdiag |= DIAG_WRITE ;
      if (isdigit((int)*argv[2]))
        nargs = 1 ;
      break ;
    case 'S':
      scale = atof(argv[2]) ;
      fprintf(stderr, "scaling parameterization by %2.1f\n", scale) ;
      nargs = 1 ;
      break ;
    case 'O':
      surface_names[1] = surface_names[2] = argv[2] ;
      nargs = 1 ;
      printf("using %s instead of smoothwm\n", argv[2]) ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      fprintf(stderr, "averaging curvature patterns %d times.\n", navgs) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
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
  printf("%s [options] <hemi> <surface name> <subject> <subject> ... "
         "<output name>\n", Progname) ;
  printf("\n");
  printf("Options:\n\n") ;
  printf(" -addframe which_field where_in_atlas \n");
  printf(" -vector : printf additional information for <addframe>\n");
  printf(" -norot  : not aligning hemispheres before averaging (default)\n");
  printf(" -rot    : rough rigid alignment of hemispheres before averaging\n");
  printf(" -annot  : zero medial wall\n");
  printf(" -overlay overlay naverages : read overlay from <overlay>\n");
  printf(" -overlay-dir dir           : use subject/<dir>/hemi.overlay\n");
  printf(" -s scale\n");
  printf(" -surf_dir dir : use 'dir' instead of 'surf' for subdirectory\n");
  printf(" -a N    : smooth curvature <N> iterations\n");
  printf(" -sdir SUBJECTS_DIR\n");
  printf("\n");
}

static void
print_help(void)
{
  print_usage() ;
  printf("This program will add a template into an average surface.\n");
  exit(1) ;
}

static void
print_version(void)
{
  printf( "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

/* set the correct set of vectorial parameters */
static void setParms(INTEGRATION_PARMS *parms)
{
  /* default template fields*/
  parms->nfields=3;

  SetFieldLabel(&parms->fields[0],
                INFLATED_CURV_CORR_FRAME,
                0,
                0.0,
                0.0,
                0,
                which_norm);
  /* only use sulc for rigid registration */
  SetFieldLabel(&parms->fields[1],SULC_CORR_FRAME,1,1.0,0.0,0,which_norm);
  SetFieldLabel(&parms->fields[2],CURVATURE_CORR_FRAME,2,0.0,0.0,0,which_norm);
}
