#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "fio.h"
#include "version.h"
#include "transform.h"

static char vcid[] = "$Id: mri_mask.c,v 1.8 2006/11/01 20:17:48 nicks Exp $";

void usage(int exit_val);

char *Progname;

static int  get_option(int argc, char *argv[]) ;

static int   InterpMethod = SAMPLE_NEAREST;

/* The following specifies the src and dst volumes 
   of the input FSL/LTA transform */
MRI          *lta_src = 0;
MRI          *lta_dst = 0;

static int invert = 0 ;
static char *xform_fname = NULL;

static float threshold = -1e10;

int main(int argc, char *argv[])
{
  char **av;
  MRI *mri_src, *mri_mask, *mri_dst ;
  int nargs, ac;
  int x, y, z;
  float value;

  LTA          *lta = 0;
  int          transform_type;

  nargs =
    handle_version_option 
    (
     argc, argv, 
     "$Id: mri_mask.c,v 1.8 2006/11/01 20:17:48 nicks Exp $", "$Name:  $"
     );
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs ;

  Progname = argv[0];

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

  if(argc != 4){
    printf("Incorrect number of arguments, argc = %d\n", argc);
    usage(1);
  }
  

  mri_src = MRIread(argv[1]) ;
  if (!mri_src)
    ErrorExit(ERROR_BADPARM, "%s: could not read source volume %s",
              Progname, argv[1]) ;
  mri_mask = MRIread(argv[2]) ;
  if (!mri_mask)
    ErrorExit(ERROR_BADPARM, "%s: could not read mask volume %s",
              Progname, argv[1]) ;

  /* Read LTA transform and apply it to mri_mask */
  if(xform_fname != NULL){
    MRI *mri_tmp;

    printf("Apply the given LTA xfrom to the mask volume\n");
    // read transform
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
          fprintf(stderr, 
                  "ERROR: fslmat does not have information on "
                  "the src and dst volumes\n");
          fprintf(stderr, 
                  "ERROR: you must give options '-lta_src' "
                  "and '-lta_dst' to specify the src and dst volume infos\n");
        }
        
        
        LTAmodifySrcDstGeom(lta, lta_src, lta_dst); 
        // add src and dst information
        LTAchangeType(lta, LINEAR_VOX_TO_VOX);
      }
      
      if (lta->xforms[0].src.valid == 0){
        if(lta_src == 0){
          fprintf(stderr, 
                  "The transform does not have the valid src volume info.\n");
          fprintf(stderr, 
                  "Either you give src volume info by option -lta_src or\n");
          fprintf(stderr, 
                  "make the transform to have the valid src info.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        }else{
          LTAmodifySrcDstGeom(lta, lta_src, NULL); // add src information
          //      getVolGeom(lta_src, &lt->src);
        }
      }
      if (lta->xforms[0].dst.valid == 0){
        if(lta_dst == 0){
          fprintf(stderr, 
                  "The transform does not have the valid dst volume info.\n");
          fprintf(stderr, 
                  "Either you give src volume info by option -lta_dst or\n");
          fprintf(stderr, 
                  "make the transform to have the valid dst info.\n");
          fprintf(stderr, 
                  "If the dst was average_305, then you can set\n");
          fprintf(stderr, 
                  "environmental variable USE_AVERAGE305 true\n");
          fprintf(stderr, 
                  "without giving the dst volume for RAS-to-RAS transform.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        }else{
          LTAmodifySrcDstGeom(lta, NULL, lta_dst); // add  dst information
        }
      }
    }
    else{
      ErrorExit(ERROR_BADPARM, 
                "transform is not of MNI, nor Register.dat, nor FSLMAT type");
    }
    
    if (invert){
      VOL_GEOM vgtmp;
      LT *lt;
      MATRIX *m_tmp = lta->xforms[0].m_L ;
      lta->xforms[0].m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
      MatrixFree(&m_tmp) ;
      lt = &lta->xforms[0];
      if (lt->dst.valid == 0 || lt->src.valid == 0){
        fprintf(stderr, 
                "WARNING:**************************************"
                "*************************\n");
        fprintf(stderr, 
                "WARNING:dst volume information is invalid.  "
                "Most likely produced wrong inverse.\n");
        fprintf(stderr, 
                "WARNING:**************************************"
                "*************************\n");
      }
      copyVolGeom(&lt->dst, &vgtmp);
      copyVolGeom(&lt->src, &lt->dst);
      copyVolGeom(&vgtmp, &lt->src);
    }
    
    //    LTAchangeType(lta, LINEAR_VOX_TO_VOX);
    mri_tmp = 
      MRIalloc(mri_src->width, 
               mri_src->height, 
               mri_src->depth, 
               mri_mask->type) ;
    MRIcopyHeader(mri_src, mri_tmp) ; 
 
    mri_tmp = LTAtransformInterp(mri_mask, mri_tmp, lta, InterpMethod);

    // mri_tmp = 
    //MRIlinearTransformInterp
    //  (
    //   mri_mask, mri_tmp, lta->xforms[0].m_L, InterpMethod
    //  );

    MRIfree(&mri_mask);

    mri_mask = mri_tmp;
    
    if (lta_src)
      MRIfree(&lta_src);
    if (lta_dst)
      MRIfree(&lta_dst);
    
    if(lta)
      LTAfree(&lta);
  }   /* if (xform_fname != NULL) */
  
  for (z = 0 ; z <mri_mask->depth ; z++)
    
    for (y = 0 ; y < mri_mask->height ; y++)
      
      for (x = 0 ; x < mri_mask->width ; x++){
        value = MRIgetVoxVal(mri_mask, x, y, z, 0);
        
        if(value <= threshold)
          MRIsetVoxVal(mri_mask,x,y,z,0,0);
        
      }
  
  mri_dst = MRImask(mri_src, mri_mask, NULL, 0, 0) ;
  if (!mri_dst)
    ErrorExit(Gerror, "%s: stripping failed", Progname) ;
  
  printf("writing masked volume to %s...\n", argv[3]) ;
  MRIwrite(mri_dst, argv[3]);
  
  MRIfree(&mri_src);
  MRIfree(&mri_mask);
  MRIfree(&mri_dst);

  exit(0);

}  /*  end main()  */

void usage(int exit_val)
{

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, 
          "Usage: %s [options] <in vol> <mask vol> <out vol>\n", 
          Progname);
  fprintf(fout, 
          "This program applies a mask volume (typically skull stripped).\n") ;

  fprintf(fout, "Options:\n") ;
  fprintf(fout, 
          "   -xform %%s    apply LTA transform to align mask "
          "to input volume\n");
  fprintf(fout, 
          "   -invert      reversely apply -xform \n");
  fprintf(fout, 
          "   -lta_src %%s  source volume for -xform "
          "(if not available from the xform file) \n");
  fprintf(fout, 
          "   -lta_dst %%s  target volume for -xform "
          "(if not available from the xform file) \n");
  fprintf(fout, 
          "   -T #         threshold mask volume at # "
          "(i.e., all values <= T is considered as zero) \n");

  fprintf(fout, "\n");

  exit(exit_val);

}  /*  end usage()  */
/*  EOF  */

/* --------------------------------------------- */
static void print_version(void)
{
  fprintf(stdout, "%s\n", vcid) ;
  exit(1) ;
}

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "help"))
    usage(1) ;
  else if (!stricmp(option, "version"))
    print_version() ;
  else if (!stricmp(option, "xform")){
    xform_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "transform file name is %s\n", xform_fname);
  }
  else if (!stricmp(option, "T")
           || !stricmp(option, "threshold")
           ){
    threshold = (float)atof(argv[2]);
    nargs = 1;
    fprintf(stderr, "threshold mask volume at %g\n", threshold);
  }
  else if (!stricmp(option, "invert")){
    invert = 1;
    fprintf(stderr, "Inversely apply the given LTA transform\n");
  }
  else if (!stricmp(option, "lta_src") ||
           !stricmp(option, "src")
           ){
    fprintf(stderr, "src volume for the given transform "
	    "(given by -xform) is %s\n",argv[2]); 
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
    fprintf(stderr, "dst volume for the transform "
	    "(given by -xform) is %s\n",argv[2]); 
    fprintf(stderr, "Reading the dst volume...\n");
    lta_dst = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!lta_dst){
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  }
  else{
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    usage(1) ;
    exit(1) ;
  }

  return(nargs) ;
}
