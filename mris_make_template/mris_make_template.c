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

static char vcid[] = "$Id: mris_make_template.c,v 1.12 2005/02/09 16:25:25 tosa Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static char *surface_names[] = 
{
  "inflated",
  "smoothwm",
  "smoothwm"
} ;

static char *curvature_names[] = 
{
  NULL,
  "sulc",
  NULL
} ;


#define IMAGES_PER_SURFACE   3   /* mean, variance, and dof */
#define SURFACES         sizeof(curvature_names) / sizeof(curvature_names[0])
#define PARAM_IMAGES         (IMAGES_PER_SURFACE*SURFACES)

static int nbrs = 1 ;
static int navgs = 0 ;
static float scale = 1 ;
static int no_rot = 0 ;
static char subjects_dir[STRLEN] ;

/* VECTORIAL_REGISTRATION */
static void setParms(INTEGRATION_PARMS *parms);
static char *FRAME_FIELD_NAMES[]=  /* order correspond to maccros defined in mrisurf.h */
	{
		NULL,
		"sulc",
		NULL, /* curvature directly computed */
		GRAYMID_NAME,
		T1MID_NAME,
		T2MID_NAME,
		PDMID_NAME,
		AMYGDALA_DIST_NAME,       
		HIPPOCAMPUS_DIST_NAME,        
		PALLIDUM_DIST_NAME,        
		PUTAMEN_DIST_NAME,          
		CAUDATE_DIST_NAME,          
		LAT_VENTRICLE_DIST_NAME,     
		INF_LAT_VENTRICLE_DIST_NAME,
	};

static int multiframes = 0;
#define NUMBER_OF_FRAMES NUMBER_OF_FIELDS_IN_VECTORIAL_REGISTRATION
#define PARAM_FRAMES  (IMAGES_PER_SURFACE*NUMBER_OF_FRAMES)

int
main(int argc, char *argv[])
{
  char         **av, surf_fname[STRLEN], *template_fname, *hemi, *sphere_name,
               *cp, *subject, fname[STRLEN] ;
  int          ac, nargs, ino, sno, nbad = 0, failed, n,failed_field[MNOFIV],ncorrs;
  VERTEX *v;
  VALS_VP *vp;
  MRI_SURFACE  *mris ;
  MRI_SP       *mrisp, /* *mrisp_aligned,*/ *mrisp_template ;
  INTEGRATION_PARMS parms ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_make_template.c,v 1.12 2005/02/09 16:25:25 tosa Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  memset(&parms, 0, sizeof(parms)) ;
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

  if (argc < 5)
    usage_exit() ;

	/* multiframe registration */
	if(multiframes){
		parms.flags |= IP_USE_MULTIFRAMES;
		setParms(&parms);
	}

  if (!strlen(subjects_dir))  /* not specified on command line*/
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment.\n",
                Progname) ;
    strcpy(subjects_dir, cp) ;
  }
  hemi = argv[1] ;
  sphere_name = argv[2] ;
  template_fname = argv[argc-1] ;
  if (1 || !FileExists(template_fname))  /* first time - create it */
  {
    fprintf(stderr, "creating new parameterization...\n") ;
    if(multiframes){
			mrisp_template = MRISPalloc(scale, PARAM_FRAMES); 
			/* 			if (no_rot)  /\* don't do rigid alignment *\/ */
			/* 				mrisp_aligned = NULL ; */
			/* 			else */
			/* 				mrisp_aligned = MRISPalloc(scale, PARAM_FRAMES);  */
		}else{
			mrisp_template = MRISPalloc(scale, PARAM_IMAGES); 
			/* 			if (no_rot)  /\* don't do rigid alignment *\/ */
			/* 				mrisp_aligned = NULL ; */
			/* 			else */
			/* 				mrisp_aligned = MRISPalloc(scale, PARAM_IMAGES);  */
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

  argv += 3 ; argc -= 3 ;
  for (ino = 0 ; ino < argc-1 ; ino++)
  {
		failed = 0 ;
		memset(failed_field,0,sizeof(int)*MNOFIV); /* for multiframes only */
    if(multiframes) setParms(&parms); /* reset parameter values to default */
    subject = argv[ino] ;
    fprintf(stderr, "processing subject %s (%d of %d)\n", subject,
						ino+1, argc-1) ;
    sprintf(surf_fname, "%s/%s/surf/%s.%s", 
            subjects_dir, subject, hemi, sphere_name) ;
    fprintf(stderr, "reading spherical surface %s...\n", surf_fname) ;
    mris = MRISread(surf_fname) ;
    if (!mris)
		{
			nbad++ ;
      ErrorPrintf(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_fname) ;
			continue ;
		}
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ; MRISstoreMetricProperties(mris) ;

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
		
		/* multiframe registration */
		if(multiframes){
			ncorrs=parms.ncorrs;
			/* allocate the VALS_VP structure */
			for( n = 0; n < mris->nvertices ; n++){
				v=&mris->vertices[n];
				vp=calloc(1,sizeof(VALS_VP));
				vp->nvals=ncorrs;
				vp->orig_vals=(float*)malloc(ncorrs*sizeof(float)); /* before blurring */
				vp->vals=(float*)malloc(ncorrs*sizeof(float));     /* values used by MRISintegrate */
				v->vp=(void*)vp;
			}
			
			/* load the different fields */
			for(n = 0 ; n < parms.ncorrs ; n++){
				if (FRAME_FIELD_NAMES[parms.corrfields[n]]){  /* read in precomputed curvature file */
					sprintf(surf_fname, "%s/%s/surf/%s.%s", subjects_dir, subject, hemi, FRAME_FIELD_NAMES[parms.corrfields[n]]) ;
					if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR){
						fprintf(stderr,"\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
						fprintf(stderr, "%s: could not read curvature file '%s'\n",Progname, surf_fname) ;
						fprintf(stderr,"setting up correlation coefficient to zero\n");
						parms.l_corrs[n]=parms.l_pcorrs[n]=0.0;
						failed_field[n]=1;
						continue;
					}
				}else{                       /* compute curvature of surface */
					sprintf(surf_fname, "%s/%s/surf/%s.%s", subjects_dir, subject, hemi, surface_names[parms.corrfields[n]]) ;
					/*					if(parms.corrfields[n]==0)
											sprintf(fname, "inflated") ;
											else
											sprintf(fname, "smoothwm") ;*/
					
					MRISsaveVertexPositions(mris, TMP_VERTICES) ;
					if (MRISreadVertexPositions(mris, surf_fname) != NO_ERROR){
						fprintf(stderr,"\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
						ErrorPrintf(ERROR_NOFILE, "%s: could not read surface file %s",
												Progname, surf_fname) ;
						fprintf(stderr,"setting up correlation coefficient to zero\n");
						parms.l_corrs[n]=parms.l_pcorrs[n]=0.0;
						failed_field[n]=1;
						continue;
					}
					MRISsetNeighborhoodSize(mris, -1) ;  /* back to max */
					MRIScomputeMetricProperties(mris) ;
					MRIScomputeSecondFundamentalForm(mris) ;
					MRISuseMeanCurvature(mris) ;
					MRISresetNeighborhoodSize(mris,1);/*only use nearest neighbor distances*/
					MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
				}
				MRISnormalizeField(mris,IS_DISTANCE_FIELD(parms.corrfields[n])); /* normalize values */
				MRISsetCurvaturesToOrigValues(mris,n);
				MRISsetCurvaturesToValues(mris,n);
			}
		}
		
		if (!no_rot && ino > 0)
		{
      if(multiframes){
				parms.frame_no = 3 ; 
				parms.l_corr = parms.l_pcorr = 0.0f ; /* don't use single field correlation functions */
				
				parms.mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
				parms.mrisp_template = mrisp_template ;
				
				MRISrigidBodyAlignVectorGlobal(mris, &parms, 1.0, 64.0, 8) ;
				if (Gdiag & DIAG_WRITE)
					MRISwrite(mris, "sphere.rot.global") ;
				MRISrigidBodyAlignVectorLocal(mris, &parms) ;
				if (Gdiag & DIAG_WRITE)
					MRISwrite(mris, "sphere.rot.local") ;
				MRISPfree(&parms.mrisp) ;
				MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
			}else{
				sprintf(surf_fname, "%s/%s/surf/%s.%s", 
								subjects_dir, subject, hemi, "sulc") ;
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
		}

				if(multiframes){
			for (n = 0; n < parms.ncorrs ; n++)
			{
				if(failed_field[n]) continue;
				
				//				fprintf(stderr, "computing parameterization for surface %s...\n",);
				
				MRISsetOrigValuesToCurvatures(mris,n);
				//				MRISnormalizeField(mris,IS_DISTANCE_FIELD(parms.corrfields[n]));
				mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
				MRISPcombine(mrisp, mrisp_template, n*3) ;
				MRISPfree(&mrisp) ;
			}
			/* free the VALS_VP structure */
			for( n = 0; n < mris->nvertices ; n++){
				v=&mris->vertices[n];
				vp=(VALS_VP*)v->vp;
				free(vp->orig_vals);
				free(vp->vals);
				free(vp);
				v->vp=NULL;
			}
			MRISfree(&mris) ;
		}else{
			for (sno = 0; sno < SURFACES ; sno++)
			{
				if (curvature_names[sno])  /* read in precomputed curvature file */
				{
					sprintf(surf_fname, "%s/%s/surf/%s.%s", 
									subjects_dir, subject, hemi, curvature_names[sno]) ;
					if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR)
					{
					nbad++ ;
          ErrorPrintf(Gerror, "%s: could not read curvature file '%s'\n",
											Progname, surf_fname) ;
					failed = 1 ;
					break ;
					}
				}
				else                       /* compute curvature of surface */
				{
					sprintf(surf_fname, "%s/%s/surf/%s.%s", 
									subjects_dir, subject, hemi, surface_names[sno]) ;
					if (MRISreadVertexPositions(mris, surf_fname) != NO_ERROR)
					{
						ErrorPrintf(ERROR_NOFILE, "%s: could not read surface file %s",
												Progname, surf_fname) ;
						nbad++ ;
						failed = 1 ;
						break ;
					}
					
					if (nbrs > 1)
						MRISsetNeighborhoodSize(mris, nbrs) ;
					MRIScomputeMetricProperties(mris) ;
					MRIScomputeSecondFundamentalForm(mris) ;
					MRISuseMeanCurvature(mris) ;
					MRISaverageCurvatures(mris, navgs) ;
					MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
					MRISnormalizeCurvature(mris) ;
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
      sprintf(surf_fname, "%s/%s/surf/%s.%s", 
              subjects_dir, subject, hemi, sphere_name) ;
      fprintf(stderr, "reading spherical surface %s...\n", surf_fname) ;
      mris = MRISread(surf_fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, surf_fname) ;
      MRIScomputeMetricProperties(mris) ; MRISstoreMetricProperties(mris) ;
      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
      sprintf(surf_fname, "%s/%s/surf/%s.%s", 
              subjects_dir, subject, hemi, "sulc") ;
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
          sprintf(surf_fname, "%s/%s/surf/%s.%s", 
                  subjects_dir, subject, hemi, curvature_names[sno]) ;
          if (MRISreadCurvatureFile(mris, surf_fname) != NO_ERROR)
            ErrorExit(Gerror, "%s: could not read curvature file '%s'\n",
                      Progname, surf_fname) ;
        }
        else                       /* compute curvature of surface */
        {
          sprintf(surf_fname, "%s/%s/surf/%s.%s", 
                  subjects_dir, subject, hemi, surface_names[sno]) ;
          if (MRISreadVertexPositions(mris, surf_fname) != NO_ERROR)
            ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                      Progname, surf_fname) ;
          
          if (nbrs > 1)
            MRISsetNeighborhoodSize(mris, nbrs) ;
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

    mrisp_tmp = mrisp_aligned ; mrisp_aligned = mrisp_template ;
    mrisp_template = mrisp_tmp ;
    MRISPfree(&mrisp_aligned) ;
  }
#endif
  fprintf(stderr, "writing updated template with %d subjects to %s...\n", argc-1-nbad, template_fname) ;
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
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
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
  else if (!stricmp(option, "norot"))
  {
    no_rot = 1 ;
    fprintf(stderr, "not aligning hemispheres before averaging.\n") ;
  }
	else if (!stricmp(option, "vector"))
  {
    multiframes = 1 ;
    fprintf(stderr, "generating vector atlas\n") ;
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
  fprintf(stderr, 
          "usage: %s [options] <hemi> <surface name> <subject> <subject> ... "
          "<output name>\n", Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program will add a template into an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

/* set the correct set of vectorial parameters */
static void setParms(INTEGRATION_PARMS *parms){
	parms->ncorrs=NUMBER_OF_FRAMES;
		
	parms->corrfields[0]=INFLATED_CURV_CORR_FRAME ; 
	parms->frames[0]=0;
	parms->l_corrs[0]=0.0f;parms->l_pcorrs[0]=0.0f;
	
	parms->corrfields[1]=SULC_CORR_FRAME;
	parms->frames[1]=1;
	parms->l_corrs[1]=1.0f;parms->l_pcorrs[1]=1.0f;
	
	parms->corrfields[2]=CURVATURE_CORR_FRAME;
	parms->frames[2]=2;
	parms->l_corrs[2]=1.0f;parms->l_pcorrs[2]=1.0f;
	
	parms->corrfields[3]=GRAYMID_CORR_FRAME;
	parms->frames[3]=3;
	parms->l_corrs[3]=1.0f;parms->l_pcorrs[3]=1.0f;
	
	parms->corrfields[4]=AMYGDALA_CORR_FRAME;
	parms->frames[4]=4;
	parms->l_corrs[4]=1.0f;parms->l_pcorrs[4]=1.0f;
	
	parms->corrfields[5]=HIPPOCAMPUS_CORR_FRAME;
	parms->frames[5]=5;
	parms->l_corrs[5]=1.0f;parms->l_pcorrs[5]=1.0f;
	
	parms->corrfields[6]=PALLIDUM_CORR_FRAME;
	parms->frames[6]=6;
	parms->l_corrs[6]=1.0f;parms->l_pcorrs[6]=1.0f;
	
	parms->corrfields[7]=PUTAMEN_CORR_FRAME;
	parms->frames[7]=7;
	parms->l_corrs[7]=1.0f;parms->l_pcorrs[7]=1.0f;
	
	parms->corrfields[8]=CAUDATE_CORR_FRAME;
	parms->frames[8]=8;
	parms->l_corrs[8]=1.0f;parms->l_pcorrs[8]=1.0f;
	
	parms->corrfields[9]=LAT_VENTRICLE_CORR_FRAME;
	parms->frames[9]=9;
	parms->l_corrs[9]=1.0f;parms->l_pcorrs[9]=1.0f;
	
	parms->corrfields[10]=INF_LAT_VENTRICLE_CORR_FRAME;
	parms->frames[10]=10;
	parms->l_corrs[10]=1.0f;parms->l_pcorrs[10]=1.0f;
}
