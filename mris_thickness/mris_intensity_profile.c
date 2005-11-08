
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
#include "transform.h"

static char vcid[] = "$Id: mris_intensity_profile.c,v 1.5 2005/11/08 18:03:31 kteich Exp $";

int main(int argc, char *argv[]) ;

MRI *MRISmeasureCorticalIntensityProfiles(MRI_SURFACE *mris, MRI *mri, int nbhd_size, 
																					float max_thick, int normalize, int curv_thresh);
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;
static char pial_name[100] = "pial" ;
static char white_name[100] = WHITE_MATTER_NAME ;
#define MAX_LABELS 10000
static char *label_names[MAX_LABELS] ;
static int nlabels = 0 ;

static int nbhd_size = 2 ;
static float max_thick = 5.0 ;
static char *sdir = NULL ;
static int num_erode = 0 ;
static float thresh ;
static int normalize = 0 ;
static char *curv_fname = NULL ;
static int  curv_thresh = 0 ;
static int navgs = 0 ;
#define MIN_BORDER_DIST 20.0  // mm from border
#define MAX_SAMPLES 20

static int max_samples = MAX_SAMPLES ;

/* The following specifies the src and dst volumes of the input FSL/LTA transform */
MRI          *lta_src = 0;
MRI          *lta_dst = 0;
static int invert = 0 ;
static char *xform_fname = NULL;
static char *mean_outname = NULL;

int
main(int argc, char *argv[])
{
  char          **av, *out_fname, *sname, buf[STRLEN], *cp, fname[STRLEN], *hemi ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;
	LABEL         *label = NULL ;
	MRI           *mri, *mri_cor, *mri_profiles ;
	LTA           *lta ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_intensity_profile.c,v 1.5 2005/11/08 18:03:31 kteich Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

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

  sname = argv[1] ;
  hemi = argv[2] ;
  out_fname = argv[4] ;
	if (!sdir)
	{
		sdir = buf ;
		cp = getenv("SUBJECTS_DIR") ;
		if (!cp)
			ErrorExit(ERROR_BADPARM, 
								"%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
		strcpy(sdir, cp) ;
	}
  
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, pial_name) ;
  fprintf(stderr, "reading pial surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

	if (curv_fname)
	{
		if (MRISreadCurvatureFile(mris, curv_fname)  != NO_ERROR)
			ErrorExit(ERROR_NOFILE, "%s: could not read curv file %s", curv_fname) ;
	}

  fprintf(stderr, "reading intensity volume from %s...\n", argv[3]) ;
  mri = MRIread(argv[3]) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read intensity volume %s",
              Progname, argv[3]) ;

  if (MRISreadOriginalProperties(mris, white_name) != NO_ERROR)
    ErrorExit(Gerror, "%s: could not read white matter surface", Progname) ;
  fprintf(stderr, "measuring gray matter thickness...\n") ;

  //use the specified xform if it's given 
  if(xform_fname != NULL){
    LTA *lta = 0;
    int transform_type;
    
    printf("INFO: Applying transformation from file %s...\n", 
	   xform_fname);
    transform_type =  TransformFileNameType(xform_fname);

    if(transform_type == MNI_TRANSFORM_TYPE || 
       transform_type == TRANSFORM_ARRAY_TYPE ||
       transform_type == REGISTER_DAT ||
       transform_type == FSLREG_TYPE
       ){
      
      printf("Reading transform ...\n");
      lta = LTAreadEx(xform_fname) ;
      
      if (!lta)
	ErrorExit(ERROR_NOFILE, "%s: could not read transform file %s",
		  Progname, xform_fname) ;
      
      if (transform_type == FSLREG_TYPE){
	if(lta_src == 0 || lta_dst == 0){
	  fprintf(stderr, "ERROR: fslmat does not have information on the src and dst volumes\n");
	  fprintf(stderr, "ERROR: you must give options '-src' and '-dst' to specify the src and dst volume infos for the registration\n");
	}
	
	
	LTAmodifySrcDstGeom(lta, lta_src, lta_dst); // add src and dst information
	//The following is necessary to interpret FSLMAT correctly!!!
	LTAchangeType(lta, LINEAR_VOX_TO_VOX);
      }
      if (lta->xforms[0].src.valid == 0){
	if(lta_src == 0){
	  fprintf(stderr, "The transform does not have the valid src volume info.\n");
	  fprintf(stderr, "Either you give src volume info by option -src or\n");
	  fprintf(stderr, "make the transform to have the valid src info.\n");
	  ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
	}else{
	  LTAmodifySrcDstGeom(lta, lta_src, NULL); // add src information
	}
      }
      if (lta->xforms[0].dst.valid == 0){
	if(lta_dst == 0){
	  fprintf(stderr, "The transform does not have the valid dst volume info.\n");
	  fprintf(stderr, "Either you give src volume info by option -dst or\n");
	  fprintf(stderr, "make the transform to have the valid dst info.\n");
	  fprintf(stderr, "If the dst was average_305, then you can set\n");
	  fprintf(stderr, "environmental variable USE_AVERAGE305 true\n");
	  fprintf(stderr, "instead.\n");
	  ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
	}else{
	  LTAmodifySrcDstGeom(lta, NULL, lta_dst); // add  dst information
	}
      }

      if (invert){
	VOL_GEOM vgtmp;
	LT *lt;
	MATRIX *m_tmp = lta->xforms[0].m_L ;
	lta->xforms[0].m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
	MatrixFree(&m_tmp) ;
	lt = &lta->xforms[0];
	if (lt->dst.valid == 0 || lt->src.valid == 0){
	  fprintf(stderr, "WARNING:***************************************************************\n");
	  fprintf(stderr, "WARNING:dst volume infor is invalid.  Most likely produce wrong inverse.\n");
	  fprintf(stderr, "WARNING:***************************************************************\n");
	}
	copyVolGeom(&lt->dst, &vgtmp);
	copyVolGeom(&lt->src, &lt->dst);
	copyVolGeom(&vgtmp, &lt->src);
      }
      

      LTAchangeType(lta, LINEAR_RAS_TO_RAS);

      printf("---------------------------------\n");
      printf("INFO: Transform Matrix \n");
      MatrixPrint(stdout,lta->xforms[0].m_L);
      printf("---------------------------------\n");
      
      // convert surface into hires volume surface if needed
      //////////////////////////////////////////////////////////////////////
      MRISsurf2surfAll(mris, mri, lta);
      
      if (lta_src)
	MRIfree(&lta_src);
      if (lta_dst)
	MRIfree(&lta_dst);
      
    }else
      {
	fprintf(stderr, "unknown transform type in file %s\n", 
		xform_fname);
	exit(1);
      }
  }
  // now convert the surface to be in the volume ras coords
  else if (mriConformed(mri) == 0)
  {
    LINEAR_TRANSFORM *lt = 0;
    MATRIX           *m_L = 0;
    TRANSFORM        *xform ;
    
    printf("volume is not conformed - transforming surfaces....\n") ;
    sprintf(fname, "%s/%s/mri/T1/COR-256", sdir, sname) ;
    if (!FileExists(fname))
      sprintf(fname, "%s/%s/mri/T1.mgz", sdir, sname) ;
    fprintf(stderr, "conformed volume from %s...\n", fname) ;
    mri_cor = MRIread(fname) ;
    if (!mri_cor)
      ErrorExit(ERROR_NOFILE, "%s: could not read conformed volume from %s",
		Progname, fname) ;
    
    fprintf(stderr, "allocating identity RAS-to-RAS xform...\n") ;
    // allocate xform->xform 
    xform = TransformAlloc(MNI_TRANSFORM_TYPE, NULL); 
    if (!xform)
      ErrorExit(ERROR_NOFILE, "%s: could not allocate hires xform %s", Progname, argv[3]) ;
    lta = (LTA *) (xform->xform);
    lt = &lta->xforms[0];
    lt->sigma = 1.0f ;
    lt->x0 = lt->y0 = lt->z0 = 0 ;
    lta->type = LINEAR_RAS_TO_RAS;
    m_L = lt->m_L;
    MatrixIdentity(4, m_L);
    getVolGeom(mri, &lt->dst);
    getVolGeom(mri_cor, &lt->src);
    if (Gdiag & DIAG_WRITE)
      {
	MRI *mri_tmp;
	mri_tmp = MRIclone(mri, NULL);
	MRIapplyRASlinearTransform(mri_cor, mri_tmp, lt->m_L);
	MRIwrite(mri_tmp, "T1_xformed.mgz") ;
	MRIfree(&mri_tmp) ;
      }
    MRIfree(&mri_cor) ;
    ///////////////////////////////////////////////////////////////////////
    // convert surface into hires volume surface if needed
    //////////////////////////////////////////////////////////////////////
    MRISsurf2surfAll(mris, mri, lta);
  }
  
  if (nlabels > 0)
    {
      int l ;
      char label_name[STRLEN] ;
      LABEL *ltotal = NULL ;
      
      for (l = 0 ; l < nlabels ; l++)
	{
	  sprintf(label_name, "%s/%s/label/%s.%s.label", sdir, sname, hemi,label_names[l]) ;
	  
	  label = LabelRead(NULL, label_name) ;
	  if (!label)
	    ErrorExit(ERROR_NOFILE, "%s: could not read label file %s...\n", Progname,
		      label_name) ;
	  if (num_erode > 0)
	    {
	      printf("eroding label %d times, npoints went from %d ", num_erode,label->n_points) ;
	      LabelErode(label, mris, num_erode) ;
	      printf("to %d ", label->n_points) ;
	    }
	  ltotal = LabelCombine(label, ltotal) ;
	}
      if (nlabels == 0)
	ltotal = LabelInFOV(mris, mri, MIN_BORDER_DIST) ;
      
      LabelRipRestOfSurfaceWithThreshold(ltotal, mris, thresh) ;
    }
  
  mri_profiles = 
    MRISmeasureCorticalIntensityProfiles(mris, mri, nbhd_size, max_thick, normalize, curv_thresh) ;

  if (navgs > 0)
    {
      printf("smoothing profiles %d times\n", navgs) ;
      MRISsmoothFrames(mris, mri_profiles, navgs) ;
    }
  MRIwrite(mri_profiles, out_fname) ;

  if(mean_outname != NULL){
    MRI *mri_mean;
    int vtx, index;
    double sum_profile;
 
    printf("Compute the mean of the intensity profile ...\n");

    mri_mean = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 1);

    for(vtx = 0; vtx < mris->nvertices; vtx++){
      sum_profile = 0;
      
      if (mris->vertices[vtx].ripflag){
	MRIsetVoxVal(mri_mean, vtx, 0, 0, 0, sum_profile);
	continue ; 
      }

      for(index= 0; index < mri_profiles->nframes; index++)
      	sum_profile += MRIFseq_vox(mri_profiles, vtx, 0, 0, index);

      sum_profile /= (float)mri_profiles->nframes;

      MRIsetVoxVal(mri_mean, vtx, 0, 0, 0, sum_profile);
    }
    
    MRIScopyMRI(mris, mri_mean, 0, "curv");

    printf("writing average intensity profiles to curvature file %s...\n", mean_outname) ;
    MRISwriteCurvature(mris, mean_outname);

    MRIfree(&mri_mean);
  }

  printf("writing cortical intensity profiles to %s...\n", out_fname) ;
  MRIfree(&mri_profiles) ;
  
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
  else if (!stricmp(option, "pial"))
  {
    strcpy(pial_name, argv[2]) ;
    fprintf(stderr,  "reading pial surface from file named %s\n", pial_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "white"))
  {
    strcpy(white_name, argv[2]) ;
    fprintf(stderr,  "reading white surface from file named %s\n", white_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "erode"))
  {
    num_erode = atoi(argv[2]) ;
    fprintf(stderr,  "eroding label %d times\n", num_erode) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "sdir"))
  {
      sdir = argv[2] ;
      nargs = 1 ;
  }
  else if (!stricmp(option, "normalize"))
  {
    normalize = 1 ;
    printf("normalizing profiles to be same length\n") ;
  }
  else if (!stricmp(option, "nsamples") || !stricmp(option, "samples"))
  {
    max_samples = atoi(argv[2]) ;
    normalize = 1 ;
    nargs = 1 ;
    printf("normalizing profiles to have %d samples\n", max_samples) ;
  }
  else if (!stricmp(option, "max"))
  {
    max_thick = atof(argv[2]) ;
    fprintf(stderr,  "limiting maximum cortical thickness to %2.2f mm.\n",
            max_thick) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "mean"))
  {
    mean_outname = argv[2];
    fprintf(stderr,  "output the intensity profile mean to file %s.\n",
            mean_outname) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "xform") || 
	   !stricmp(option, "at")
	   ){
    xform_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "transform file name is %s\n", xform_fname);
  }
  else if (!stricmp(option, "ait")
	   ){
    xform_fname = argv[2];
    invert = 1;
    nargs = 1;
    fprintf(stderr, "Inversely apply the transform given by %s\n", xform_fname);
  }
  else if (!stricmp(option, "invert")){
    invert = 1;
    fprintf(stderr, "Inversely apply the given registration transform\n");
  }
  else if (!stricmp(option, "lta_src") ||
	   !stricmp(option, "src")
	   ){
    fprintf(stderr, "src volume for the given transform (given by -xform) is %s\n",argv[2]); 
    fprintf(stderr, "Reading the src volume...\n");
    lta_src = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!lta_src){
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  }
  else if (!stricmp(option, "lta_dst") ||
	   !stricmp(option, "dst")
	   ){
    fprintf(stderr, "dst volume for the transform (given by -xform) is %s\n",argv[2]); 
    fprintf(stderr, "Reading the dst volume...\n");
    lta_dst = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!lta_dst){
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  }
  else switch (toupper(*option))
  {
         case 'C':
	   curv_fname = argv[2] ;
	   curv_thresh = atoi(argv[3]) ;
	   printf("reading curvature from %s, and limiting calculations to %s regions\n",
		  curv_fname, curv_thresh > 0 ? "sulcal" : "gyral") ;
	   nargs = 2 ;
	   break ;
	case 'L':
	  label_names[nlabels] = argv[2] ;
	  nargs = 1 ;
	  printf("limiting profile calculation to label %s\n", label_names[nlabels]) ;
	  nlabels++ ;
	  break ;
  case 'N':
    nbhd_size = atoi(argv[2]) ;
    fprintf(stderr, "using neighborhood size=%d\n", nbhd_size) ;
    nargs = 1 ;
    break ;
  case 'A':
    navgs = atoi(argv[2]) ;
    nargs = 1;
    printf("smoothing profiles %d times across space\n", navgs) ;
    break ;
  case 'T':
    thresh = atof(argv[2]) ;
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
          "usage: %s [options] <subject name> <hemi> <volume> <output file>\n", 
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program computes the intensity profile of the cortical ribbon\n"
          "and writes the resulting measurement into a 'curvature' file\n"
          "<output file>.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-sdir %%s specifies the SUBJECTS_DIR \n") ;
  fprintf(stderr, "-white %%s specifies WHITE surface filename.\n") ;
  fprintf(stderr, "-pial %%s specifies PIAL surface filename.\n") ;
  fprintf(stderr, "-normalize normalized sampling w.r.t. thickness  \n") ;
  fprintf(stderr, "-mean %%s outputs the mean profile-integral to the specified file (output is in curv format).\n") ;
  fprintf(stderr, "-xform %%s specify the registration xform that maps the T1 volume (from which the surfaces were computed) to the input volume to be sampled.\n") ;
  fprintf(stderr, "-src %%s source volume used when computing the registration xform\n");
  fprintf(stderr, "-dst %%s destination volume used when computing the registration xform\n");
  fprintf(stderr, "-invert this flag asks to apply the registration xform inversely \n");
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


MRI *
MRISmeasureCorticalIntensityProfiles(MRI_SURFACE *mris, MRI *mri, int nbhd_size, 
																		 float max_thick,
																		 int normalize, int curv_thresh)
{
  int     vno, n, vlist[100000], vtotal, ns, i, vnum, min_n,
          pial_vno, nsamples ;
  VERTEX  *v, *vn, *vn2 ;
  float   d, dx, dy, dz, dist, min_dist, nx, ny, nz, dot, sample_dist, thick ;
	Real    x, y, z, xv, yv, zv, val ;
	MRI     *mri_profiles ;

	mri_profiles = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, max_samples) ;

  /* current vertex positions are gray matter, orig are white matter */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
	{
		if (!(vno % 25000))
			fprintf(stdout, "%d of %d vertices processed\n", vno,mris->nvertices) ;
		v = &mris->vertices[vno] ;
		if (v->ripflag)
			continue ;
		if (curv_thresh != 0 && curv_thresh*v->curv < 0)
			continue ;
		nx = v->nx ; ny = v->ny ; nz = v->nz ;
		if (vno == Gdiag_no)
			DiagBreak() ;
		dx = v->x - v->origx ; dy = v->y - v->origy ; dz = v->z - v->origz ; 
		pial_vno = vno ;
		min_dist = sqrt(dx*dx + dy*dy + dz*dz) ;
		v->marked = 1 ; vtotal = 1 ; vlist[0] = vno ;
		min_n = 0 ;
		for (ns = 1 ; ns <= nbhd_size ; ns++)
		{
			vnum = 0 ;  /* will be # of new neighbors added to list */
			for (i = 0 ; i < vtotal ; i++)
	    {
	      vn = &mris->vertices[vlist[i]] ;
	      if (vn->ripflag)
					continue ;
	      if (vn->marked && vn->marked < ns-1)
					continue ;
	      for (n = 0 ; n < vn->vnum ; n++)
				{
					vn2 = &mris->vertices[vn->v[n]] ;
					if (vn2->ripflag || vn2->marked)  /* already processed */
						continue ;
					vlist[vtotal+vnum++] = vn->v[n] ;
					vn2->marked = ns ;
					dx = vn2->x-v->origx ; dy = vn2->y-v->origy ; dz = vn2->z-v->origz ;
					dot = dx*nx + dy*ny + dz*nz ;
					if (dot < 0) /* must be outwards from surface */
						continue ;
					dot = vn2->nx*nx + vn2->ny*ny + vn2->nz*nz ;
					if (dot < 0) /* must be outwards from surface */
						continue ;
					dist = sqrt(dx*dx + dy*dy + dz*dz) ;
					if (dist < min_dist)
					{
						min_n = ns ;
						min_dist = dist ;
						if (min_n == nbhd_size && DIAG_VERBOSE_ON)
							fprintf(stdout, "%d --> %d = %2.3f\n",
											vno,vn->v[n], dist) ;
						pial_vno = vn->v[n] ;
					}
				}
	    }
			vtotal += vnum ;
		}

		// unmark stuff for next time
		for (n = 0 ; n < vtotal ; n++)
		{
			vn = &mris->vertices[vlist[n]] ;
			if (vn->ripflag)
				continue ;
			vn->marked = 0 ;
		}

		vn2 = &mris->vertices[pial_vno] ;
		// vector pointing from white to pial
		dx = vn2->x-v->origx ; dy = vn2->y-v->origy ; dz = vn2->z-v->origz ;
		thick = sqrt(dx*dx + dy*dy + dz*dz) ;
		dx /= thick ; dy /= thick ; dz /= thick ;
#define SAMPLE_DIST 0.25
		if (normalize)
		{
			sample_dist = thick / (max_samples-1);
			for (nsamples = 0, d = 0.0 ; nsamples < max_samples ; d += sample_dist, nsamples++)
			{
				x = v->origx + d*dx ; y = v->origy + d*dy ; z = v->origz + d*dz ;
				MRIsurfaceRASToVoxel(mri, x, y, z, &xv, &yv, &zv) ;
				MRIsampleVolume(mri, xv, yv, zv, &val) ;
				MRIFseq_vox(mri_profiles, vno, 0, 0, nsamples) = val ;
			}
		}
		else
		{
			sample_dist = SAMPLE_DIST ;
			//			for (nsamples = 0, d = 0.0 ; d <= max_thick ; d += sample_dist, nsamples++)
			for (nsamples = 0, d = 0.0 ;nsamples < max_samples; d += sample_dist, nsamples++)
			{
				if (d > thick)  // so each row has the same # of cols
				{
					MRIFseq_vox(mri_profiles, vno, 0, 0, nsamples) = -1 ;
				}
				else
				{
					x = v->origx + d*dx ; y = v->origy + d*dy ; z = v->origz + d*dz ;
					MRIsurfaceRASToVoxel(mri, x, y, z, &xv, &yv, &zv) ;
					MRIsampleVolume(mri, xv, yv, zv, &val) ;
					MRIFseq_vox(mri_profiles, vno, 0, 0, nsamples) = val ;
				}
			}
		}
	}

  return(mri_profiles) ;
}
