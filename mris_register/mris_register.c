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
#include "gcsa.h"

static char vcid[] = "$Id: mris_register.c,v 1.24 2005/02/05 23:37:59 segonne Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  compute_area_ratios(MRI_SURFACE *mris) ;
static double gcsaSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;

static int max_passes = 4 ;
static float min_degrees = 0.5 ;
static float max_degrees = 64.0 ;
static int   nangles = 8 ;
static int nbrs = 1 ;
static float scale = 1.0f ;

static int reverse_flag = 0 ;

static float dalpha = 0.0f ;
static float dbeta = 0.0f ;
static float dgamma = 0.0f ;

char *Progname ;
static char curvature_fname[STRLEN] = "" ;
static char *orig_name = "smoothwm" ;
static char *jacobian_fname = NULL ;

#define MAX_LABELS 100
static int nlabels = 0 ;
static LABEL *labels[MAX_LABELS] ;
static char  *label_names[MAX_LABELS] ;
static GCSA  *label_gcsa[MAX_LABELS] ;
static int   label_indices[MAX_LABELS] ;

/* multiframe registration */
static int multiframes = 0;
#define NUMBER_OF_FRAMES NUMBER_OF_FIELDS_IN_VECTORIAL_REGISTRATION
#define PARAM_FRAMES  (IMAGES_PER_SURFACE*NUMBER_OF_FRAMES)

static void initParms(void);
static void setParms(void);

static int use_defaults = 1 ;

static INTEGRATION_PARMS  parms ;

int
main(int argc, char *argv[])
{
  char         **av, *surf_fname, *template_fname, *out_fname, fname[STRLEN],*cp;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;
  MRI_SP       *mrisp_template ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_register.c,v 1.24 2005/02/05 23:37:59 segonne Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag = DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.projection = PROJECT_SPHERE ;
  parms.flags |= IP_USE_CURVATURE ;
  parms.tol = 1e-0*10 ;
  parms.min_averages = 0 ;
  parms.l_area = 0.0 ;
  parms.l_parea = 0.2f ;
  parms.l_dist = 0.1 ;
  parms.l_corr = 1.0f ;
  parms.l_nlarea = 1 ;
  parms.l_pcorr = 0.0f ;
  parms.niterations = 25 ;
  parms.n_averages = 256 ;
  parms.write_iterations = 100 ;
  parms.dt_increase = 1.01 /* DT_INCREASE */;
  parms.dt_decrease = 0.99 /* DT_DECREASE*/ ;
  parms.error_ratio = 1.03 /*ERROR_RATIO */;
  parms.dt_increase = 1.0 ;
  parms.dt_decrease = 1.0 ;
	parms.l_external = 1000 ;   /* in case manual label is specified */
  parms.error_ratio = 1.1 /*ERROR_RATIO */;
  parms.integration_type = INTEGRATE_ADAPTIVE ;
  parms.integration_type = INTEGRATE_MOMENTUM /*INTEGRATE_LINE_MINIMIZE*/ ;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
  parms.dt = 0.9 ;
  parms.momentum = 0.95 ;
  parms.desired_rms_height = -1.0 ;
  parms.nbhd_size = 0 ;
  parms.max_nbrs = 0 ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  surf_fname = argv[1] ;
  template_fname = argv[2] ;
  out_fname = argv[3] ;

  if (parms.base_name[0] == 0)
  {
    FileNameOnly(out_fname, fname) ;
    cp = strchr(fname, '.') ;
    if (cp)
      strcpy(parms.base_name, cp+1) ;
    else
      strcpy(parms.base_name, "sphere") ;
  }

  fprintf(stderr, "reading surface from %s...\n", surf_fname) ;
  mris = MRISread(surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_fname) ;

  if (!FZERO(dalpha) || !FZERO(dbeta) || !FZERO(dgamma))
    MRISrotate(mris, mris, RADIANS(dalpha), RADIANS(dbeta), 
               RADIANS(dgamma)) ;

  if (curvature_fname[0])
  {
    fprintf(stderr, "reading source curvature from %s\n",curvature_fname) ;
    MRISreadCurvatureFile(mris, curvature_fname) ;
  }
  fprintf(stderr, "reading template parameterization from %s...\n",
          template_fname) ;
  mrisp_template = MRISPread(template_fname) ;
  if (!mrisp_template)
    ErrorExit(ERROR_NOFILE, "%s: could not open template file %s",
                Progname, template_fname) ;
  if (use_defaults)
  {
    if (*IMAGEFseq_pix(mrisp_template->Ip, 0, 0, 2) <= 1.0)  /* 1st time */
    {
      parms.l_dist = 0.5 ; parms.l_corr = 1.0 ; parms.l_parea = 0.1 ;
    }
    else   /* subsequent alignments */
    {
      parms.l_dist = 0.1 ; parms.l_corr = 1.0 ; parms.l_parea = 0.2 ;
    }
		if(multiframes) parms.l_corr=0.0f;
#if 0
		// test XXXX
		parms.l_dist = 0.0 ; parms.l_corr = 0.0 ; parms.l_parea = 0.0 ;
		parms.l_area = 0.0 ;
		parms.l_parea = 0.0f ;
		parms.l_dist = 0.0 ;
		parms.l_corr = 0.0f ;
		parms.l_nlarea = 0.0f ;
		parms.l_pcorr = 0.0f ;
#endif
  }

  if (nbrs > 1)
    MRISsetNeighborhoodSize(mris, nbrs) ;
  MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
  if (reverse_flag)
    MRISreverse(mris, REVERSE_X) ; 
  mris->status = MRIS_PARAMETERIZED_SPHERE ;
  MRIScomputeMetricProperties(mris) ;
  if (!FZERO(parms.l_dist))
    MRISscaleDistances(mris, scale) ;
#if 0
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISzeroNegativeAreas(mris) ;
  MRISstoreMetricProperties(mris) ;
#endif
  MRISstoreMeanCurvature(mris) ;  /* use curvature from file */
  /*  MRISsetOriginalFileName(mris, orig_name) ;*/
  MRISreadOriginalProperties(mris, orig_name) ;
	if(multiframes)
		MRISvectorRegister(mris, mrisp_template, &parms, max_passes, min_degrees, max_degrees, nangles) ;
	else
		MRISregister(mris, mrisp_template, &parms, max_passes, min_degrees, max_degrees, nangles) ;
  fprintf(stderr, "writing registered surface to %s...\n", out_fname) ;
  MRISwrite(mris, out_fname) ;
  if (jacobian_fname)
  {
    MRIScomputeMetricProperties(mris) ;
    compute_area_ratios(mris) ;  /* will put results in v->curv */
#if 0
    MRISwriteArea(mris, jacobian_fname) ;
#else
    MRISwriteCurvature(mris, jacobian_fname) ;
#endif
  }

  MRISPfree(&mrisp_template) ;
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
  int    nargs = 0 ;
  char   *option ;
  float  f ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
	else if (!stricmp(option, "vector"))
  {
    multiframes = 1 ;
		setParms();
    fprintf(stderr, "using vectorial registration \n") ;
  }
	else if (!stricmp(option, "hippocampus"))
  {
		if(multiframes==0){
			multiframes = 1 ;
			initParms();
		}
		parms.l_corrs[5]=atof(argv[2]);
    fprintf(stderr, "using hippocampus distance map\n") ;
		nargs=1;
  }
	else if (!stricmp(option, "curvature"))
  {
		if(multiframes==0){
			multiframes = 1 ;
			initParms();
		}
		parms.l_corrs[2]=atof(argv[2]);
    fprintf(stderr, "using curvature map \n") ;
		nargs=1;
  }
	else if (!stricmp(option, "sulcus"))
  {
		if(multiframes==0){
			multiframes = 1 ;
			initParms();
		}
		parms.l_corrs[1]=atof(argv[2]);
    fprintf(stderr, "using sulcal depth map \n") ;
		nargs=1;
  }
	else if (!stricmp(option, "amygdala"))
  {
		if(multiframes==0){
			multiframes = 1 ;
			initParms();
		}
		parms.l_corrs[4]=atof(argv[2]);
    fprintf(stderr, "using amygdala distance map \n") ;
		nargs=1;
  }
  else if (!stricmp(option, "vnum") || !stricmp(option, "distances"))
  {
    parms.nbhd_size = atof(argv[2]) ;
    parms.max_nbrs = atof(argv[3]) ;
    nargs = 2 ;
    fprintf(stderr, "nbr size = %d, max neighbors = %d\n",
            parms.nbhd_size, parms.max_nbrs) ;
  }
  else if (!stricmp(option, "rotate"))
  {
    dalpha = atof(argv[2]) ;
    dbeta = atof(argv[3]) ;
    dgamma = atof(argv[4]) ;
    fprintf(stderr, "rotating brain by (%2.2f, %2.2f, %2.2f)\n",
            dalpha, dbeta, dgamma) ;
    nargs = 3 ;
  }
  else if (!stricmp(option, "reverse"))
  {
    reverse_flag = 1 ;
    fprintf(stderr, "mirror image reversing brain before morphing...\n") ;
  }
  else if (!stricmp(option, "min_degrees"))
  {
    min_degrees = atof(argv[2]) ;
    fprintf(stderr, "setting min angle for search to %2.2f degrees\n", min_degrees) ;
		nargs = 1 ;
  }
  else if (!stricmp(option, "max_degrees"))
  {
    max_degrees = atof(argv[2]) ;
    fprintf(stderr, "setting max angle for search to %2.2f degrees\n", max_degrees) ;
		nargs = 1 ;
  }
  else if (!stricmp(option, "nangles"))
  {
    nangles = atoi(argv[2]) ;
    fprintf(stderr, "setting # of angles/search per scale to %d\n", nangles) ;
		nargs = 1 ;
  }
  else if (!stricmp(option, "jacobian"))
  {
    jacobian_fname = argv[2] ;
    nargs = 1 ;
    printf("writing out jacobian of mapping to %s\n", jacobian_fname) ;
  }
  else if (!stricmp(option, "dist"))
  {
    sscanf(argv[2], "%f", &parms.l_dist) ;
    nargs = 1 ;
    use_defaults = 0 ;
    fprintf(stderr, "l_dist = %2.3f\n", parms.l_dist) ;
  }
  else if (!stricmp(option, "norot"))
  {
    fprintf(stderr, "disabling initial rigid alignment...\n") ;
    parms.flags |= IP_NO_RIGID_ALIGN ;
  }
  else if (!stricmp(option, "lm"))
  {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    fprintf(stderr, "integrating with line minimization\n") ;
  }
  else if (!stricmp(option, "search"))
  {
    parms.integration_type = INTEGRATE_LM_SEARCH ;
    fprintf(stderr, "integrating with binary search line minimization\n") ;
  }
  else if (!stricmp(option, "dt"))
  {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = .2*parms.dt ;
    nargs = 1 ;
    fprintf(stderr, "momentum with dt = %2.2f\n", parms.dt) ;
  }
  else if (!stricmp(option, "area"))
  {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "parea"))
  {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_parea) ;
    nargs = 1 ;
    fprintf(stderr, "using l_parea = %2.3f\n", parms.l_parea) ;
  }
  else if (!stricmp(option, "nlarea"))
  {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_nlarea) ;
    nargs = 1 ;
    fprintf(stderr, "using l_nlarea = %2.3f\n", parms.l_nlarea) ;
  }
  else if (!stricmp(option, "spring"))
  {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_spring) ;
    nargs = 1 ;
    fprintf(stderr, "using l_spring = %2.3f\n", parms.l_spring) ;
  }
  else if (!stricmp(option, "corr"))
  {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_corr) ;
    nargs = 1 ;
    fprintf(stderr, "using l_corr = %2.3f\n", parms.l_corr) ;
  }
  else if (!stricmp(option, "curv"))
  {
    parms.flags |= IP_USE_CURVATURE ;
    fprintf(stderr, "using smoothwm curvature for final alignment\n") ;
  }
  else if (!stricmp(option, "nocurv"))
  {
    parms.flags &= ~IP_USE_CURVATURE ;
    fprintf(stderr, "using smoothwm curvature for final alignment\n") ;
  }
  else if (!stricmp(option, "adaptive"))
  {
    parms.integration_type = INTEGRATE_ADAPTIVE ;
    fprintf(stderr, "using adaptive time step integration\n") ;
  }
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
  }
  else if (!stricmp(option, "tol"))
  {
    if (sscanf(argv[2], "%e", &f) < 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan tol from %s",
                Progname, argv[2]) ;
    parms.tol = (double)f ;
    nargs = 1 ;
    fprintf(stderr, "using tol = %2.2e\n", (float)parms.tol) ;
  }
  else if (!stricmp(option, "error_ratio"))
  {
    parms.error_ratio = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "error_ratio=%2.3f\n", parms.error_ratio) ;
  }
  else if (!stricmp(option, "dt_inc"))
  {
    parms.dt_increase = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt_increase=%2.3f\n", parms.dt_increase) ;
  }
  else if (!stricmp(option, "vnum"))
  {
    parms.nbhd_size = atof(argv[2]) ;
    parms.max_nbrs = atof(argv[3]) ;
    nargs = 2 ;
    fprintf(stderr, "nbr size = %d, max neighbors = %d\n",
            parms.nbhd_size, parms.max_nbrs) ;
  }
  else if (!stricmp(option, "dt_dec"))
  {
    parms.dt_decrease = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt_decrease=%2.3f\n", parms.dt_decrease) ;
  }
  else switch (toupper(*option))
  {
  case 'M':
    parms.integration_type = INTEGRATE_MOMENTUM ;
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "momentum = %2.2f\n", (float)parms.momentum) ;
    break ;
	case 'L':
		if (nlabels >= MAX_LABELS-1)
			ErrorExit(ERROR_NO_MEMORY, "%s: too many labels specified (%d max)", Progname, MAX_LABELS) ;
		nargs = 3 ;
		labels[nlabels] = LabelRead(NULL, argv[2]) ;
		if (labels[nlabels] == NULL)
			ErrorExit(ERROR_NOFILE, "%s: could not read label file %s", Progname, argv[2]) ;
		label_gcsa[nlabels] = GCSAread(argv[3]) ;
		if (label_gcsa[nlabels] == NULL)
			ErrorExit(ERROR_NOFILE, "%s: could not read GCSA file %s", Progname, argv[3]) ;
		label_names[nlabels] = argv[4] ;
		label_indices[nlabels] = CTABnameToAnnotation(label_gcsa[nlabels]->ct, argv[4]) ;

		if (label_indices[nlabels] < 0)
			ErrorExit(ERROR_NOFILE, "%s: could not map name %s to index", Progname, argv[3]) ;
		nlabels++ ;
		gMRISexternalSSE = gcsaSSE ;
		break ;
	case 'E':
		parms.l_external = atof(argv[2]) ;
		nargs = 1 ;
		printf("setting l_external = %2.1f\n", parms.l_external) ;
		break ;
  case 'C':
    strcpy(curvature_fname, argv[2]) ;
    nargs = 1 ;
    break ;
  case 'A':
    sscanf(argv[2], "%d", &parms.n_averages) ;
    nargs = 1 ;
    fprintf(stderr, "using n_averages = %d\n", parms.n_averages) ;
    break ;
  case 'S':
    scale = atof(argv[2]) ;
    fprintf(stderr, "scaling distances by %2.2f\n", scale) ;
    nargs = 1 ;
    break ;
  case 'N':
    sscanf(argv[2], "%d", &parms.niterations) ;
    nargs = 1 ;
    fprintf(stderr, "using niterations = %d\n", parms.niterations) ;
    break ;
  case 'W':
    Gdiag |= DIAG_WRITE ;
    sscanf(argv[2], "%d", &parms.write_iterations) ;
    nargs = 1 ;
    fprintf(stderr, "using write iterations = %d\n", parms.write_iterations) ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'O':
    orig_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "using %s for original properties...\n", orig_name) ;
    break ;
  case 'P':
    max_passes = atoi(argv[2]) ;
    fprintf(stderr, "limitting unfolding to %d passes\n", max_passes) ;
    nargs = 1 ;
    break ;
  case '?':
	case 'H':
  case 'U':
    print_help() ;
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
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
       "usage: %s [options] <input surface> <average surface> <output surface>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program register a surface with  an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
	fprintf(stderr, "\n\t-l <label file> <atlas (*.gcs)> <label name>\n"
					"\tthis option will specify a manual label to align with atlas label <label name>\n");
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static int
compute_area_ratios(MRI_SURFACE *mris)
{
  VERTEX  *v ;
  int     vno ;
  float   area_scale ;

  area_scale = mris->total_area / mris->orig_area  ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    v->curv = v->area / (v->origarea*area_scale) ;
  }

  return(NO_ERROR) ;
}

static double 
gcsaSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
	int       vno, ano, lno, vno_prior, n, found ;
	VERTEX    *v, *v_prior ;
	double    sse ;
	LABEL     *area ;
  CP_NODE   *cpn ;
  CP        *cp ;
	GCSA      *gcsa ;

	for (sse = 0.0, ano = 0 ; ano < nlabels ; ano++)
	{
		area = labels[ano] ;
		gcsa = label_gcsa[ano] ;
		for (lno = 0 ; lno < area->n_points ; lno++)
		{
			vno = area->lv[lno].vno ;
			if (vno < 0)
				continue ;
			v = &mris->vertices[vno] ;
			found = 0 ;
			v_prior = GCSAsourceToPriorVertex(gcsa, v) ;
			vno_prior = v_prior - gcsa->mris_priors->vertices ;
			cpn = &gcsa->cp_nodes[vno_prior] ;
			for (n = 0 ; n < cpn->nlabels ; n++)
			{
				cp = &cpn->cps[n] ;
				if (cpn->labels[n] == label_indices[ano])
				{
					found = 1 ;
					break ;
				}
			}
			if (found == 0)
				sse += parms->l_external ;
		}
	}

	return(sse) ;
}

void initParms(void){
	int n;
	parms.l_corr=parms.l_pcorr=0.0f; 
	parms.flags |= IP_USE_MULTIFRAMES; 
	parms.ncorrs=NUMBER_OF_FRAMES;
	parms.corrfields[0]=INFLATED_CURV_CORR_FRAME ; parms.frames[0]=0;
	parms.corrfields[1]=SULC_CORR_FRAME;parms.frames[1]=1;
	parms.corrfields[2]=CURVATURE_CORR_FRAME;parms.frames[2]=2;
	parms.corrfields[3]=GRAYMID_CORR_FRAME;parms.frames[3]=3;
	parms.corrfields[4]=AMYGDALA_CORR_FRAME;parms.frames[4]=4;
	parms.corrfields[5]=HIPPOCAMPUS_CORR_FRAME;parms.frames[5]=5;
	parms.corrfields[6]=PALLIDUM_CORR_FRAME;parms.frames[6]=6;
	parms.corrfields[7]=PUTAMEN_CORR_FRAME;parms.frames[7]=7;
	parms.corrfields[8]=CAUDATE_CORR_FRAME;parms.frames[8]=8;
	parms.corrfields[9]=LAT_VENTRICLE_CORR_FRAME;parms.frames[9]=9;
	parms.corrfields[10]=INF_LAT_VENTRICLE_CORR_FRAME;parms.frames[10]=10;
	for(n = 0 ; n < NUMBER_OF_FRAMES;n++){
		parms.frames[n]=n;
		parms.l_corrs[n]=parms.l_pcorrs[n]=0.0f;
	}
}

void setParms(void){
	parms.l_corr=parms.l_pcorr=0.0f; 
	parms.flags |= IP_USE_MULTIFRAMES; 
	parms.ncorrs=NUMBER_OF_FRAMES;
	parms.corrfields[0]=INFLATED_CURV_CORR_FRAME ; parms.frames[0]=0;parms.l_corrs[0]=0.0f;parms.l_pcorrs[0]=0.0f;
		parms.corrfields[1]=SULC_CORR_FRAME;parms.frames[1]=1;parms.l_corrs[1]=1.0f;parms.l_pcorrs[1]=0.0f;
		parms.corrfields[2]=CURVATURE_CORR_FRAME;parms.frames[2]=2;parms.l_corrs[2]=1.0f;parms.l_pcorrs[2]=0.0f;
		parms.corrfields[3]=GRAYMID_CORR_FRAME;parms.frames[3]=3;parms.l_corrs[3]=1.0f;parms.l_pcorrs[3]=0.0f;
		parms.corrfields[4]=AMYGDALA_CORR_FRAME;parms.frames[4]=4;parms.l_corrs[4]=10.0f;parms.l_pcorrs[4]=0.0f;    /* amygdala */
		parms.corrfields[5]=HIPPOCAMPUS_CORR_FRAME;parms.frames[5]=5;parms.l_corrs[5]=10.0f;parms.l_pcorrs[5]=0.0f; /* hippocampus */
		parms.corrfields[6]=PALLIDUM_CORR_FRAME;parms.frames[6]=6;parms.l_corrs[6]=1.0f;parms.l_pcorrs[6]=0.0f;
		parms.corrfields[7]=PUTAMEN_CORR_FRAME;parms.frames[7]=7;parms.l_corrs[7]=10.0f;parms.l_pcorrs[7]=0.0f;    /* putamen */
		parms.corrfields[8]=CAUDATE_CORR_FRAME;parms.frames[8]=8;parms.l_corrs[8]=10.0f;parms.l_pcorrs[8]=0.0f;    /* caudate */
		parms.corrfields[9]=LAT_VENTRICLE_CORR_FRAME;parms.frames[9]=9;parms.l_corrs[9]=1.0f;parms.l_pcorrs[9]=0.0f;
		parms.corrfields[10]=INF_LAT_VENTRICLE_CORR_FRAME;parms.frames[10]=10;parms.l_corrs[10]=1.0f;parms.l_pcorrs[10]=0.0f;
}
