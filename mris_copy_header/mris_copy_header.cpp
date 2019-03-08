#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "diag.h"
#include "error.h"
#include "utils.h"
#include "macros.h"
#include "histo.h"
#include "mri.h"
#include "mrisurf.h"

const char *Progname ;


int
main(int argc, char *argv[])
{
  MRI_SURFACE    *mris_in, *mris_template ;
  

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc < 4)
    ErrorExit(ERROR_BADPARM, "%s: <input surface> <template surface>  <output surface>",Progname);

  mris_in = MRISread(argv[1]) ;
  if (mris_in == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load source MRIS from %s\n", Progname, argv[1]) ;

  mris_template = MRISread(argv[2]) ;
  if (mris_template == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load template MRIS from %s\n", Progname, argv[2]) ;

  if (mris_template->mri_sras2vox)
    mris_in->mri_sras2vox = MRIcopy(mris_template->mri_sras2vox, NULL) ;
  if (mris_template->lta)
    mris_in->lta = LTAcopy(mris_template->lta, NULL) ;
  if (mris_template->SRASToTalSRAS_)
    mris_in->SRASToTalSRAS_ = MatrixCopy(mris_template->SRASToTalSRAS_,NULL) ;
  if (mris_template->TalSRASToSRAS_)
    mris_in->TalSRASToSRAS_ = MatrixCopy(mris_template->TalSRASToSRAS_,NULL) ;
  if (mris_template->m_sras2vox)
    mris_in->m_sras2vox = MatrixCopy(mris_template->m_sras2vox,NULL) ;
  strcpy(mris_in->subject_name, mris_template->subject_name) ;
  *(&mris_in->vg) = *(&mris_template->vg) ;
  mris_in->useRealRAS = mris_template->useRealRAS ;
  mris_in->ct = mris_template->ct ;

  printf("writing new surface to %s\n", argv[3]) ;
  MRISwrite(mris_in, argv[3]) ;
  return(0) ;
}
  

