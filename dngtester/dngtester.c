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

MATRIX *MRIcorVox2RAS(MATRIX *vox2ras);

char *Progname = "dngtester";
char *subject=NULL, *hemi=NULL, *surfname=NULL, *outsurfname=NULL;
char *SUBJECTS_DIR = NULL;
MRI *mri, *mri2, *template;
char tmpstr[2000];
int err;
GCA  *gca ;
MATRIX *V2Rsrc,*V2Rtemplate,*V2V;
MRIS *mris;
GCSA *gcsa;
char *fsh;
GCA_MORPH *gcam;
float c,r,s;
char *dcmfile, *dcmdir;

/*----------------------------------------*/
int main(int argc, char **argv)
{
  char **FileNames;
  DICOMInfo RefDCMInfo;
  int nfiles;
  //int ndcmfiles, nthfile;

  dcmfile = argv[1];
  dcmdir = fio_dirname(dcmfile);
  printf("dcmfile = %s\n",dcmfile);
  printf("dcmdir = %s\n",dcmdir);
  if(! IsDICOM(dcmfile)){
    printf("ERROR: %s is not a dicom file\n",dcmfile);
    exit(1);
  }
  // SeriesNo = 0x20, 0x11
  GetDICOMInfo(dcmfile, &RefDCMInfo, FALSE, 1);


  // scan directory
  err=ScanDir(dcmdir, &FileNames, &nfiles);
  printf("Nfiles = %d\n",nfiles);

  


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

