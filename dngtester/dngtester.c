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

int SortDCMFileInfo(DICOMInfo **dcmfi_list, int nlist);
int CompareDCMFileInfo(const void *a, const void *b);
int DCMCountFrames(DICOMInfo **dcmfi_list, int nlist);

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
  DICOMInfo RefDCMInfo, **dcminfo;
  int nfiles, nframes,nslices, r, c, s, f;
  int ndcmfiles, nthfile, mritype=0,nvox;
  unsigned short *v16=NULL;
  //unsigned char *v8=NULL;
  DCM_ELEMENT *element;
  MATRIX *vox2ras,*CRScenter,*RAScenter;

  dcmfile = argv[1];
  dcmdir = fio_dirname(dcmfile);
  printf("dcmfile = %s\n",dcmfile);
  printf("dcmdir = %s\n",dcmdir);
  if(!IsDICOM(dcmfile)){
    printf("ERROR: %s is not a dicom file\n",dcmfile);
    exit(1);
  }
  // SeriesNo = 0x20, 0x11
  GetDICOMInfo(dcmfile, &RefDCMInfo, FALSE, 1);
  printf("Ref Series No = %d\n",RefDCMInfo.SeriesNumber);

  // scan directory
  err=ScanDir(dcmdir, &FileNames, &nfiles);
  if(err) exit(1);
  printf("Found %d files, checking for dicom\n",nfiles);

  ndcmfiles=0;
  for(nthfile = 0; nthfile < nfiles; nthfile ++){
    //printf("%d %s\n",nthfile,FileNames[nthfile]);
    if(!IsDICOM(FileNames[nthfile])) continue;
    ndcmfiles ++;
  }
  printf("Found %d dicom files\n",ndcmfiles);

  dcminfo = (DICOMInfo **)calloc(ndcmfiles,sizeof(DICOMInfo *));
  ndcmfiles=0;
  for(nthfile = 0; nthfile < nfiles; nthfile ++){
    //printf("%d %s\n",nthfile,FileNames[nthfile]);
    if(!IsDICOM(FileNames[nthfile])) continue;
    dcminfo[ndcmfiles] = (DICOMInfo *)calloc(1,sizeof(DICOMInfo));
    GetDICOMInfo(FileNames[nthfile], dcminfo[ndcmfiles], FALSE, 1);    
    ndcmfiles ++;
  }

  printf("Sorting\n");
  SortDCMFileInfo(dcminfo, ndcmfiles);
  for(nthfile = 0; nthfile < ndcmfiles; nthfile ++){
    printf("%d %s\n",nthfile,dcminfo[nthfile]->FileName);
  }
  printf("Counting Frames\n");
  nframes = DCMCountFrames(dcminfo, ndcmfiles);
  printf("nframes = %d\n",nframes);
  nslices = ndcmfiles/nframes;
  printf("nslices = %d\n",nslices);
  if(nslices*nframes != ndcmfiles){
    printf("ERROR: the number of frames * number of slices does\n"
	   "not equal the number of dicom files.\n");
    exit(1);
  }
  if(RefDCMInfo.BitsAllocated == 8)  mritype = MRI_UCHAR;
  if(RefDCMInfo.BitsAllocated == 16) mritype = MRI_SHORT;
  //mri = MRIallocHeader(ncols,nrows,nslices,MRI_SHORT);
  mri = MRIallocSequence(RefDCMInfo.Columns,RefDCMInfo.Rows,
			 nslices,mritype,nframes);
  if(mri==NULL){
    printf("ERROR: mri alloc failed\n");
    exit(1);
  }
  mri->tr      = RefDCMInfo.RepetitionTime;
  mri->te      = RefDCMInfo.EchoTime;
  mri->ti      = RefDCMInfo.InversionTime;
  mri->flip_angle = RefDCMInfo.FlipAngle;
  mri->width   = RefDCMInfo.Columns;
  mri->height  = RefDCMInfo.Rows;
  mri->depth   = nslices;
  mri->nframes = nframes;
  mri->fov     = RefDCMInfo.FieldOfView;
  mri->thick   = RefDCMInfo.SliceThickness;
  mri->xsize   = RefDCMInfo.xsize;
  mri->ysize   = RefDCMInfo.ysize;
  mri->zsize   = RefDCMInfo.SliceThickness;
  mri->x_r     = -RefDCMInfo.Vc[0];
  mri->x_a     = -RefDCMInfo.Vc[1];
  mri->x_s     = +RefDCMInfo.Vc[2];
  mri->y_r     = -RefDCMInfo.Vr[0];
  mri->y_a     = -RefDCMInfo.Vr[1];
  mri->y_s     = +RefDCMInfo.Vr[2];
  mri->z_r     = -RefDCMInfo.Vs[0];
  mri->z_a     = -RefDCMInfo.Vs[1];
  mri->z_s     = +RefDCMInfo.Vs[2];

  vox2ras = MRIxfmCRS2XYZ(mri, 0);
  vox2ras->rptr[1][4] = -dcminfo[0]->ImagePosition[0];
  vox2ras->rptr[2][4] = -dcminfo[0]->ImagePosition[1];
  vox2ras->rptr[3][4] = +dcminfo[0]->ImagePosition[2];
  CRScenter = MatrixZero(4,1,NULL);
  // These should not be -1, but consistent with old version
  CRScenter->rptr[1][1] = (RefDCMInfo.Columns-1)/2.0;
  CRScenter->rptr[2][1] = (RefDCMInfo.Rows-1)/2.0;
  CRScenter->rptr[3][1] = (nslices-1)/2.0;
  CRScenter->rptr[4][1] = 1;
  RAScenter = MatrixMultiply(vox2ras,CRScenter,NULL);
  mri->c_r = RAScenter->rptr[1][1];
  mri->c_a = RAScenter->rptr[2][1];
  mri->c_s = RAScenter->rptr[3][1];

  // Return here if only reading in the header

  nvox = RefDCMInfo.Columns * RefDCMInfo.Rows;
  if(RefDCMInfo.BitsAllocated == 16){
    nthfile = 0;
    for(s=0; s < nslices; s++){
      for(f=0; f < nframes; f++){
	printf("slice %2d, frame %2d\n",s,f);
	element = GetElementFromFile(dcminfo[nthfile]->FileName,0x7FE0,0x10);
	if(element == NULL){
	  printf("ERROR: reading pixel data from %s\n",dcminfo[nthfile]->FileName);
	  MRIfree(&mri);
	}
	v16 = (unsigned short *)(element->d.string);
	for(r=0; r < RefDCMInfo.Rows; r++) {
	  for(c=0; c < RefDCMInfo.Columns; c++) {
	    MRISseq_vox(mri,c,r,s,f) = *(v16++);
	  }
	}
	FreeElementData(element);
	free(element);
	nthfile++;
      } // frame
    } // slice
  } // 16 bit

  MRIwrite(mri,"tp1b.mgh");

  return(0);
}

/*----------------------------------------------------------*/
int CompareDCMFileInfo(const void *a, const void *b)
{

  DICOMInfo *dcmfi1;
  DICOMInfo *dcmfi2;
  int n;
  float dv[3], dvsum2, dot, actm1, actm2;

  dcmfi1 = *((DICOMInfo **) a);
  dcmfi2 = *((DICOMInfo **) b);

  /* ------ Sort by Slice Position -------- */
  /* Compute vector from the first to the second */
  dvsum2 = 0;
  for(n=0; n < 3; n++){
    dv[n] = dcmfi2->ImagePosition[n] - dcmfi1->ImagePosition[n];
    dvsum2 += (dv[n]*dv[n]);
  }
  for(n=0; n < 3; n++) dv[n] /= sqrt(dvsum2); /* normalize */
  /* Compute dot product with Slice Normal vector */
  dot = 0;
  for(n=0; n < 3; n++) dot += (dv[n] * dcmfi1->Vs[n]);
  // Now Sort by slice position.
  if(dot > +0.5) return(-1);
  if(dot < -0.5) return(+1);

  /* Sort by Series Number */
  if(dcmfi1->SeriesNumber < dcmfi2->SeriesNumber) return(-1);
  if(dcmfi1->SeriesNumber > dcmfi2->SeriesNumber) return(+1);

  /* Sort by Image Number (Temporal Sequence) */
  if(dcmfi1->ImageNumber < dcmfi2->ImageNumber) return(-1);
  if(dcmfi1->ImageNumber > dcmfi2->ImageNumber) return(+1);
  
  /* Sort by Echo Number  */
  if(dcmfi1->EchoNumber < dcmfi2->EchoNumber) return(-1);
  if(dcmfi1->EchoNumber > dcmfi2->EchoNumber) return(+1);
  
  /* Sort by Acquisition Time */
  sscanf(dcmfi1->AcquisitionTime,"%f",&actm1);
  sscanf(dcmfi2->AcquisitionTime,"%f",&actm2);
  if(actm1 < actm2) return(-1);
  if(actm1 > actm2) return(+1);

  printf("WARNING: files are not found to be different and "
         "cannot be sorted\n");
  printf("File1: %s\n",dcmfi1->FileName);
  printf("File2: %s\n",dcmfi2->FileName);

  return(0);
}

/*-----------------------------------------------------------*/
int SortDCMFileInfo(DICOMInfo **dcmfi_list, int nlist)
{
  qsort(dcmfi_list,nlist,sizeof(DICOMInfo **),CompareDCMFileInfo);
  return(0);
}

/*-----------------------------------------------------------*/
int DCMCountFrames(DICOMInfo **dcmfi_list, int nlist)
{
  // Assumes they have been sorted
  int nframes, nth, c, stop;
  double ImgPos0[3],d;
  DICOMInfo *dcmfi;

  nframes = 0;
  for(nth=0;nth < nlist; nth++){
    dcmfi = dcmfi_list[nth];
    if(nth==0){
      for(c=0; c < 3; c++) 
	ImgPos0[c] = dcmfi->ImagePosition[c];
    }
    stop=0;
    for(c=0; c < 3; c++){
      d = fabs(ImgPos0[c] - dcmfi->ImagePosition[c]);
      if(d > .001) {stop=1; break;}
    }
    if(stop) break;
    nframes++;
  }
  return(nframes);
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

