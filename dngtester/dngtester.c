#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "mrisurf.h"
#include "utils.h"
#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "volume_io.h"
#include "region.h"
#include "machine.h"
#include "fio.h"
#include "mri_identify.h" 
#include "mrisurf.h" 
#include "fmriutils.h" 
#include "gca.h"
#include "gcsa.h"
#include "fsgdf.h"
#include "icosahedron.h"
#include "gca.h"
#include "gcamorph.h"
#include "DICOMRead.h"

char *Progname = "dngtester";
char *subject=NULL, *hemi=NULL, *surfname=NULL, *outsurfname=NULL;
char *SUBJECTS_DIR = NULL;
char tmpstr[2000];
MRI *mri, *mri2, *template;
int err;
char *fsh;
char *dcmfile, *dcmdir;

/*----------------------------------------*/
int main(int argc, char **argv)
{
  dcmfile = argv[1];
  mri = DICOMRead2(dcmfile,1);
  MRIwrite(mri,"tp1b.mgh");
  return(0);
}

/*------------------------------------------------------------------------
  MRIcorVox2RAS() - computes vox2ras for the standard COR volume, ie,
  256^3, 1mm^3. The RAS is in "tkregister" space (also known as  "surface"
  space). 
  ------------------------------------------------------------------------*/
MATRIX *MRIcorVox2RAS(MATRIX *vox2ras)
{
  if(vox2ras==NULL) vox2ras = MatrixConstVal(0,4,4,NULL);
  vox2ras->rptr[1][1] = -1;
  vox2ras->rptr[1][4] = +128;
  vox2ras->rptr[2][3] = +1;
  vox2ras->rptr[2][4] = -128;
  vox2ras->rptr[3][2] = -1;
  vox2ras->rptr[3][4] = +128;
  vox2ras->rptr[4][4] = +1;
  return(vox2ras);
}

