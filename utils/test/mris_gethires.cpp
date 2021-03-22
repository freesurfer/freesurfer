/**
 * @brief ?
 *
 */
/*
 * Original Author: Y. Tosa
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

#include <iostream>
#include <iomanip>

extern "C"
{
#include "error.h"
#include "mri.h"
#include "mrisurf.h"

const char *Progname = "mris_gethires";
}

using namespace std;

int main(int argc, char *argv[])
{
  MRI *src = 0;
  if (argc < 2)
  {
    cerr << "Usage: mris_gethighres <surface> <highresvol> <outsurf> [<lowresvol>]" << endl;
    cerr << "   where lowresvol is only used when surface does not contain the volume info" << endl;
    return 0;
  }
  MRIS *mris = MRISread(argv[1]);
  if (!mris)
  {
    cerr << "could not load surface: " << argv[1] << endl;
    return 0;
  }
  if (mris->vg.valid == 0)
  {
    cerr << "Trying to read lowres volume to get info" << endl;
    src = MRIreadHeader(argv[4], MRI_VOLUME_TYPE_UNKNOWN);
    if (!src)
    {
      cerr << "could not read lowres volume : " << argv[2] << endl;
      return 0;
    }
    getVolGeom(src, &mris->vg);
  }
  MRI *mri_hires = MRIread(argv[2]);
  if (!mri_hires)
  {
    cerr << "could not highres volume : " << argv[2] << endl;
    return 0;
  }
  LINEAR_TRANSFORM *lt = 0;
  MATRIX *m_L = 0;
  fprintf(stderr, "allocating identity RAS-to-RAS xform...\n") ;
  // allocate hires_xform->xform
  TRANSFORM *hires_xform = TransformAlloc(MNI_TRANSFORM_TYPE, NULL);
  if (!hires_xform)
    ErrorExit(ERROR_NOFILE, "%s: could not allocate hires xform %s", Progname, argv[3]) ;
  LTA *hires_lta = (LTA *) (hires_xform->xform);
  lt = &hires_lta->xforms[0];
  lt->sigma = 1.0f ;
  lt->x0 = lt->y0 = lt->z0 = 0 ;
  hires_lta->type = LINEAR_RAS_TO_RAS;
  m_L = lt->m_L;
  MatrixIdentity(4, m_L);
  // geometry for src and dst are not specified yet
  MRISsurf2surf(mris, mri_hires, hires_lta);

  MRISwrite(mris, argv[3]);

  MRIfree(&mri_hires);
  if (src)
    MRIfree(&src);
  MRISfree(&mris);
  TransformFree(&hires_xform);
}
