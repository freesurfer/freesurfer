/**
 * @brief program to find all edits made to a subject and write out
 *    a .mgz volume summarizing them.
 *
 * combines all the edits made to various volumes during the recon
 * stream and writes them into a single volume with different labels.
 *
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mri.h"
#include "colortab.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static char sdir[STRLEN] = "";

#define EDIT_WM_OFF           1
#define EDIT_WM_ON            2
#define EDIT_BRAIN_OFF        3
#define EDIT_BRAIN_ON         4
#define EDIT_BM_CHANGED       5 // brainmask.mgz changed
#define EDIT_FINALSURFS_OFF   6
#define EDIT_FINALSURFS_ON    7
#define EDIT_ASEG_CHANGED     8
#define CTAB_ENTRIES          EDIT_ASEG_CHANGED+1

int
main(int argc, char *argv[])
{
  char         **av, fname[STRLEN] ;
  int          ac, nargs, i ;
  char         *subject, *cp, mdir[STRLEN], *out_fname;
  const char *name ;
  int          r, g, b, nedits = 0 ;
  MRI          *mri=NULL,*mri_edits=NULL,*mri_aseg_auto=NULL,*mri_bm_auto=NULL;
  FILE         *ctfp;

  // default output file name:
  out_fname = strcpyalloc("edits.mgz");

  nargs = handleVersionOption(argc, argv, "mri_compile_edits");
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

  if (argc < 2)
    usage_exit(1) ;

  if (strlen(sdir) == 0)
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (cp == NULL)
      ErrorExit
        (ERROR_BADPARM,
         "%s: SUBJECTS_DIR must be defined in the env or on cmdline "
         "with -sdir", Progname) ;
    strcpy(sdir, cp) ;
  }

  subject = argv[1] ;
  if (argv[2]) out_fname = argv[2] ;
  printf("Compiling volume edits for subject %s...\n", subject) ;
  fflush(stdout);

  int req = snprintf(mdir, STRLEN, "%s/%s/mri", sdir, subject) ;  
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  /*
   * brain.mgz
   */
  req = snprintf(fname, STRLEN, "%s/brain.mgz", mdir) ; 
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  
  mri = MRIread(fname) ;
  if (mri)
  {
    if(NULL == mri_edits) mri_edits = MRIclone(mri, NULL) ;
    int edits = MRIsetVoxelsWithValue(mri,
                                      mri_edits,
                                      WM_EDITED_OFF_VAL,
                                      EDIT_BRAIN_OFF) ;
    if (edits) printf("Found %d OFF edits in brain.mgz\n",edits);
    nedits += edits;
    edits = MRIsetVoxelsWithValue(mri,
                                  mri_edits,
                                  WM_EDITED_ON_VAL,
                                  EDIT_BRAIN_ON) ;
    if (edits) printf("Found %d ON edits in brain.mgz\n",edits);
    nedits += edits;
    MRIfree(&mri) ;
    fflush(stdout);
  }

  /*
   * wm.mgz
   */
  req = snprintf(fname, STRLEN, "%s/wm.mgz", mdir) ;  
  if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri = MRIread(fname) ;
  if (mri)
  {
    if(NULL == mri_edits) mri_edits = MRIclone(mri, NULL) ;
    int edits = MRIsetVoxelsWithValue(mri,
                                      mri_edits,
                                      WM_EDITED_OFF_VAL,
                                      EDIT_WM_OFF) ;
    if (edits) printf("Found %d OFF edits in wm.mgz\n",edits);
    nedits += edits;
    edits = MRIsetVoxelsWithValue(mri, mri_edits,
                                  WM_EDITED_ON_VAL,
                                  EDIT_WM_ON) ;
    if (edits) printf("Found %d ON edits in wm.mgz\n",edits);
    nedits += edits;
    MRIfree(&mri) ;
    fflush(stdout);
  }

  /*
   * brainmask.mgz
   */
  req = snprintf(fname, STRLEN, "%s/brainmask.mgz", mdir) ; 
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri = MRIread(fname) ;
  if (mri)
  {
    if(NULL == mri_edits) mri_edits = MRIclone(mri, NULL) ;
    int req = snprintf(fname, STRLEN, "%s/brainmask.auto.mgz", mdir) ;  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    mri_bm_auto = MRIread(fname) ;
    if (mri_bm_auto)
    {
      int edits = MRIsetDifferentVoxelsWithValue(mri,
                                                 mri_bm_auto,
                                                 mri_edits,
                                                 EDIT_BM_CHANGED);
      nedits += edits;
      MRIfree(&mri) ;
      MRIfree(&mri_bm_auto) ;
      if (edits) printf("Found %d edits in brainmask.mgz\n",edits);
      fflush(stdout);
    }
  }

  /*
   * brain.finalsurfs.mgz
   */
  req = snprintf(fname, STRLEN, "%s/brain.finalsurfs.mgz", mdir) ;   
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri = MRIread(fname) ;
  if (mri)
  {
    if(NULL == mri_edits) mri_edits = MRIclone(mri, NULL) ;
    int edits = MRIsetVoxelsWithValue(mri,
                                      mri_edits,
                                      WM_EDITED_OFF_VAL,
                                      EDIT_FINALSURFS_OFF) ;
    if (edits) printf("Found %d OFF edits in brain.finalsurfs.mgz\n",edits);
    nedits += edits;
    edits = MRIsetVoxelsWithValue(mri,
                                  mri_edits,
                                  WM_EDITED_ON_VAL,
                                  EDIT_FINALSURFS_ON) ;
    if (edits) printf("Found %d ON edits in brain.finalsurfs.mgz\n",edits);
    nedits += edits;
    MRIfree(&mri) ;
    fflush(stdout);
  }

  /*
   * aseg.mgz
   */
  req = snprintf(fname, STRLEN, "%s/aseg.mgz", mdir) ;  
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri = MRIread(fname) ;
  if (mri)
  {
    if(NULL == mri_edits) mri_edits = MRIclone(mri, NULL) ;
    int req = snprintf(fname, STRLEN, "%s/aseg.auto.mgz", mdir) ;   
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    mri_aseg_auto = MRIread(fname) ;
    if (mri_aseg_auto)
    {
      int edits = MRIsetDifferentVoxelsWithValue(mri,
                                                 mri_aseg_auto,
                                                 mri_edits,
                                                 EDIT_ASEG_CHANGED);
      nedits += edits;
      MRIfree(&mri) ;
      MRIfree(&mri_aseg_auto) ;
      if (edits) printf("Found %d edits in aseg.mgz\n",edits);
      fflush(stdout);
    }
  }

  if (mri_edits)
  {
    mri_edits->ct = CTABalloc(CTAB_ENTRIES) ;
    strcpy (mri_edits->fname, "mri_compile_edits");
    for (i = 0 ; i < CTAB_ENTRIES ; i++)
    {
      switch (i)
      {
      case EDIT_WM_OFF:
        name = "wm-OFF" ;
        r = 0 ;
        g = 0 ;
        b = 255 ;
        break ;
      case EDIT_WM_ON:
        name = "wm-ON" ;
        r = 255 ;
        g = 0 ;
        b = 0 ;
        break ;
      case EDIT_BRAIN_OFF:
        name = "brain-OFF" ;
        r = 0 ;
        g = 255 ;
        b = 255 ;
        break ;
      case EDIT_BRAIN_ON:
        name = "brain-ON" ;
        r = 255 ;
        g = 255 ;
        b = 0 ;
        break ;
      case EDIT_BM_CHANGED:
        name = "brainmask-CHANGED" ;
        r = 0 ;
        g = 64 ;
        b = 255 ;
        break ;
      case EDIT_FINALSURFS_OFF:
        name = "brain.finalsurf-OFF" ;
        r = 0 ;
        g = 128 ;
        b = 255 ;
        break ;
      case EDIT_FINALSURFS_ON:
        name = "brain.finalsurfs-ON" ;
        r = 255 ;
        g = 128 ;
        b = 0 ;
        break ;
      case EDIT_ASEG_CHANGED:
        name = "aseg-CHANGED" ;
        r = 255 ;
        g = 255 ;
        b = 128 ;
        break ;
      default:
        name = "Unknown" ;
        r = 0 ;
        g = 0 ;
        b = 0 ;
        break ;
        continue ;
      }

      strcpy(mri_edits->ct->entries[i]->name, name) ;
      mri_edits->ct->entries[i]->ri = r ;
      mri_edits->ct->entries[i]->gi = g ;
      mri_edits->ct->entries[i]->bi = b ;
      mri_edits->ct->entries[i]->ai = 255;

      /* Now calculate the float versions. */
      mri_edits->ct->entries[i]->rf =
        (float)mri_edits->ct->entries[i]->ri / 255.0;
      mri_edits->ct->entries[i]->gf =
        (float)mri_edits->ct->entries[i]->gi / 255.0;
      mri_edits->ct->entries[i]->bf =
        (float)mri_edits->ct->entries[i]->bi / 255.0;
      mri_edits->ct->entries[i]->af =
        (float)mri_edits->ct->entries[i]->ai / 255.0;
    }
  }

  // recon-all -show-edits greps on this text, so dont change it
  if (nedits && mri_edits)
  {
    printf("%d mri_compile_edits_found, saving results to %s\n",
           nedits, out_fname) ;
    MRIwrite(mri_edits, out_fname);
    if (mri_edits->ct)
    {
      int req = snprintf(fname, STRLEN, "%s/mri_compile_edits_LUT", mdir) ; 
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("Colortable saved to %s :\n",fname);
      CTABprintASCII(mri_edits->ct,stdout);
      ctfp = fopen(fname,"w");
      CTABprintASCII(mri_edits->ct,ctfp);
      fclose(ctfp);
    }

    printf("Edits can be viewed with command:\n");
    if (mri_edits->ct)
    {
      printf("tkmedit %s T1.mgz -segmentation %s %s\n",
             subject,out_fname,fname);
    }
    else
    {
      printf("tkmedit %s T1.mgz -segmentation %s\n",subject,out_fname);
    }
  }
  else
  {
    printf("0 mri_compile_edits_found\n") ;
  }

  exit(0) ;
  return(0) ;
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
  if (!stricmp(option, "sdir"))
  {
    strcpy(sdir, argv[2]) ;
    nargs = 1;
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
  }
  else switch (toupper(*option))
  {
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


/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static void
usage_exit(int code)
{
  printf("usage: %s [options] <subject name> <output volume>\n",
         Progname) ;
  printf("program to create a single volume showing\n"
         "all the volumetric edits made to a subject\n") ;
  exit(code) ;
}





