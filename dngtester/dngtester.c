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

MRI *ReverseMorph(char *subject);
int GCAMmorphSurfToAtlas(MRIS *mris, GCA_MORPH *gcam);

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

/*----------------------------------------*/
int main(int argc, char **argv)
{
  
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  subject     = argv[1];
  hemi        = argv[2];
  surfname    = argv[3];
  outsurfname = argv[4];
  subject = "fsr-tst";
  hemi = "lh";
  surfname = "white";
  outsurfname = "white.tal.nonlin";


  sprintf(tmpstr,"%s/%s/mri/transforms/talairach.m3z",SUBJECTS_DIR,subject);
  printf("Reading %s\n",tmpstr);
  gcam = GCAMreadAndInvert(tmpstr);
  if(gcam == NULL) exit(1);
  Gdiag_no = 1;


  if(!strcmp(surfname,outsurfname)){
    printf("ERROR: input and output surfaces are the same\n");
    exit(1);
  }

  printf("subj=%s, hemi=%s, surf=%s, out=%s\n",subject,hemi,surfname,outsurfname);
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
  printf("Reading %s\n",tmpstr);
  mris = MRISread(tmpstr);
  if(mris == NULL) exit(1);

  printf("Morphing\n");
  GCAMmorphSurfToAtlas(mris, gcam);

  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,outsurfname);
  printf("Writing to %s\n",tmpstr);
  MRISwrite(mris,tmpstr);

  exit(0);
  /*----------------------------------------------------------------*/

  fsh = getenv("FREESURFER_HOME");
  //set CPAtlas = $FREESURFER_HOME/average/$hemi.$GCS
  //set GCS     = curvature.buckner40.filled.desikan_killiany.gcs


  printf("Reading ico\n");
  mris = ReadIcoByOrder(7, 100);

  sprintf(tmpstr,"%s/average/%s.curvature.buckner40.filled.desikan_killiany.gcs",fsh,hemi);
  printf("Reading gcsa  %s\n",tmpstr);
  gcsa = GCSAread(tmpstr);
  mris->ct = gcsa->ct;

  printf("Building\n");
  GCSAbuildMostLikelyLabels(gcsa,mris);

  MRISmodeFilterAnnotations(mris, 2);
  sprintf(tmpstr,"./%s.aparc.annot",hemi);
  MRISwriteAnnotation(mris, tmpstr);



  exit(1);
  /*-------------------------------------------------*/

  //set GCA = RB40_talairach_2005_12_30.gca
  //$FREESURFER_HOME/average/$GCA
  printf("Reading %s\n",argv[2]) ;
  template = MRIreadHeader(argv[2],MRI_VOLUME_TYPE_UNKNOWN);
  if(!template) exit(1);

  printf("Reading %s\n",argv[1]) ;
  gca = GCAread(argv[1]) ;
  printf("Building\n");
  mri = GCAbuildMostLikelyVolume(gca, NULL) ;

  printf("Upsampling\n");
  V2Rsrc = MRIxfmCRS2XYZtkreg(mri);
  V2Rtemplate = MRIxfmCRS2XYZtkreg(template);
  V2V = MatrixMultiply(MatrixInverse(V2Rsrc,NULL),V2Rtemplate,NULL);
  MatrixPrint(stdout,V2V);

  mri2 = MRIallocSequence(template->width, template->height,template->depth,
                             MRI_UCHAR, 1);
  MRIcopyHeader(template, mri2) ;
  template->nframes = 1;
  MRIvol2Vol(mri, mri2, NULL, SAMPLE_NEAREST, 0);

  printf("Writings\n");
  MRIwrite(mri,"T1MLV.lowres.mgz");
  MRIwrite(mri2,"T1MLV.mgz");

  MRIfree(&mri);
  MRIfree(&mri2);

  mri = GCAbuildMostLikelyLabelVolume(gca);
  mri2 = MRIallocSequence(template->width, template->height,template->depth,
                             MRI_UCHAR, 1);
  MRIcopyHeader(template, mri2) ;
  template->nframes = 1;
  MRIvol2Vol(mri, mri2, NULL, SAMPLE_NEAREST, 0);

  MRIwrite(mri,"aseg.mlv.lowres.mgz");
  MRIwrite(mri2,"aseg.mlv.mgz");

  return(0);
}
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
int GCAMmorphSurfToAtlas(MRIS *mris, GCA_MORPH *gcam)
{
  char *sd, tmpstr[2000];
  MRI *mri;
  MATRIX *vox2ras, *ras2vox, *ras, *crs;
  int vtxno;
  VERTEX *v;
  extern char *subject;
  Real Acval,Arval,Asval,Mcval,Mrval,Msval;
  float Acval2,Arval2,Asval2;

  sd = getenv("SUBJECTS_DIR");
  //sprintf(tmpstr,"%s/%s/mri/orig.mgz",sd,mris->subject_name);
  sprintf(tmpstr,"%s/%s/mri/orig.mgz",sd,subject);
  mri = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
  if(mri==NULL){
    printf("ERROR: reading %s\n",tmpstr);
    return(1);
  }

  printf("Inverting \n");
  GCAMinvert(gcam, mri);

  vox2ras = MRIxfmCRS2XYZtkreg(mri);
  ras2vox = MatrixInverse(vox2ras,NULL);
  ras = MatrixAlloc(4,1,MATRIX_REAL);
  ras->rptr[4][1] = 1;
  crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[4][1] = 1;

  printf("Appling Inverse Morph \n");
  for(vtxno = 0; vtxno < mris->nvertices; vtxno++){

    // convert vertex RAS to CRS in anatomical volume
    v = &(mris->vertices[vtxno]);
    ras->rptr[1][1] = v->x;
    ras->rptr[2][1] = v->y;
    ras->rptr[3][1] = v->z;
    crs = MatrixMultiply(ras2vox,ras,crs);
    Acval = crs->rptr[1][1];
    Arval = crs->rptr[2][1];
    Asval = crs->rptr[3][1];

    // compute the CRS where this point will move to in the morph space
    MRIsampleVolume(gcam->mri_xind,Acval,Arval,Asval,&Mcval);
    MRIsampleVolume(gcam->mri_yind,Acval,Arval,Asval,&Mrval);
    MRIsampleVolume(gcam->mri_zind,Acval,Arval,Asval,&Msval);

    Mcval *= gcam->spacing;
    Mrval *= gcam->spacing;
    Msval *= gcam->spacing;

    // compute RAS of the morphed point
    crs->rptr[1][1] = Mcval;
    crs->rptr[2][1] = Mrval;
    crs->rptr[3][1] = Msval;
    ras = MatrixMultiply(vox2ras,crs,ras);

    if(Gdiag_no > 0){
      GCAMsampleMorph(gcam, Mcval, Mrval, Msval, &Acval2, &Arval2, &Asval2); // test
      printf("------------------------------------------------\n");
      printf("%5d   ras=(%6.1f,%6.1f,%6.1f), crsA =(%6.1f,%6.1f,%6.1f) \n"
	     "     crsM=(%6.1f,%6.1f,%6.1f), crsA2=(%6.1f,%6.1f,%6.1f) \n"
	     "     rasM=(%6.1f,%6.1f,%6.1f)\n",
	     vtxno, v->x,v->y,v->z, Acval,Arval,Asval, Mcval,Mrval,Msval,
	     Acval2,Arval2,Asval2, ras->rptr[1][1],ras->rptr[2][1],ras->rptr[3][1]);
      GCAMsampleMorphCheck(gcam, .01, Mcval, Mrval, Msval);
    }

    // pack it back into the vertex
    v->x = ras->rptr[1][1];
    v->y = ras->rptr[2][1];
    v->z = ras->rptr[3][1];
  }
  return(0);
}
/*----------------------------------------------------------------------------
  GCAMsampleMorphRAS() - given an RAS coord in the Morph space (xyzMorph),
  compute the RAS in the Anat space (xyzAnat). Has not been tested yet.
  ----------------------------------------------------------------------------*/
int GCAMsampleMorphRAS(GCA_MORPH *gcam, float xMorph, float yMorph, float zMorph,
		       float  *xAnat,  float  *yAnat,  float  *zAnat, 
		       MATRIX *vox2ras)
{
  static MATRIX *ras=NULL, *crs=NULL, *ras2vox=NULL;
  float  cMorph, rMorph, sMorph;
  float  cAnat, rAnat, sAnat;

  ras2vox = MatrixInverse(vox2ras,ras2vox);
  if(!ras){
    ras = MatrixAlloc(4,1,MATRIX_REAL);
    ras->rptr[4][1] = 1;
  }
  if(!crs){
    crs= MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
  }

  ras->rptr[1][1] = xMorph;
  ras->rptr[2][1] = yMorph;
  ras->rptr[3][1] = zMorph;
  crs = MatrixMultiply(ras2vox,ras,crs);
  cMorph = crs->rptr[1][1];
  rMorph = crs->rptr[2][1];
  sMorph = crs->rptr[3][1];
  
  err = GCAMsampleMorph(gcam, cMorph, rMorph, sMorph, &cAnat, &rAnat, &sAnat);
  if(err) return(1);

  crs->rptr[1][1] = cAnat;
  crs->rptr[2][1] = rAnat;
  crs->rptr[3][1] = sAnat;
  ras = MatrixMultiply(vox2ras,crs,ras);

  *xAnat = ras->rptr[1][1];
  *yAnat = ras->rptr[2][1];
  *zAnat = ras->rptr[3][1];

  return(0);
}



/*----------------------------------------------------------------------------*/


MRI *GCAMmorphToAtlas2(MRI *mri_src, GCA_MORPH *gcam, MRI *mri_morphed, int frame)
{
  int        width, height, depth, x, y, z, start_frame, end_frame,nframes ;
  float      xd, yd, zd ;
  Real       val, xoff, yoff, zoff ;

  if (frame >= 0 && frame < mri_src->nframes)
    start_frame = end_frame = frame ;
  else{
    start_frame = 0 ; end_frame = mri_src->nframes-1 ;
  }
  nframes = end_frame - start_frame;

  width  = gcam->width  * gcam->spacing ; 
  height = gcam->height * gcam->spacing ; 
  depth  = gcam->depth  * gcam->spacing ; 

  // GCAM is a non-linear voxel-to-voxel transform
  // it also assumes that the uniform voxel size
  if(mri_morphed) {
    if ( (mri_src->xsize != mri_src->ysize)
	 || (mri_src->xsize != mri_src->zsize)
	 || (mri_src->ysize != mri_src->zsize))		{
      ErrorExit(ERROR_BADPARM, "non-uniform volumes cannot be used for GCAMmorphToAtlas()\n");
    }
  }
  if(!mri_morphed) {
    mri_morphed = MRIallocSequence(width, height, depth, mri_src->type, nframes) ;
    MRIcopyHeader(mri_src, mri_morphed) ;
  }

  if(getenv("MGH_TAL"))	{
    xoff = -7.42 ;
    yoff = 24.88 ;
    zoff = -18.85 ;
    printf("INFO: adding MGH tal offset (%2.1f, %2.1f, %2.1f) to xform\n", xoff, yoff, zoff) ;
  }
  else  xoff = yoff = zoff = 0 ;

  // step through gcam volume 1mm at a time
  for (x = 0 ; x < width ; x++)	{
    for (y = 0 ; y < height ; y++)		{
      for (z = 0 ; z < depth ; z++)			{
	if (x == Gx && y == Gy && z == Gz) DiagBreak() ;

	if (!GCAMsampleMorph(gcam, (float)x*mri_src->thick, 
			     (float)y*mri_src->thick, (float)z*mri_src->thick, 
			     &xd, &yd, &zd)) {
	  xd /= mri_src->thick ; yd /= mri_src->thick ; zd /= mri_src->thick ; 
	  xd += xoff ; yd += yoff ; zd += zoff ;
	  for (frame = start_frame ; frame <= end_frame ; frame++) {
	    if (xd > -1 && yd > -1 && zd > 0 &&
		xd < mri_src->width && yd < mri_src->height && zd < mri_src->depth)
	      MRIsampleVolumeFrameType(mri_src, xd, yd, zd, frame, SAMPLE_TRILINEAR, &val) ;
	    else
	      val = 0.0 ;
	    MRIsetVoxVal(mri_morphed, x, y, z, frame-start_frame, val) ;
	  }
	}
      }
    }
  }

  // copy the gcam dst information to the morphed volume
  if (getenv("USE_AVERAGE305")) {
      fprintf(stderr, "INFO: Environmental variable USE_AVERAGE305 set\n");
      fprintf(stderr, "INFO: Modifying dst c_(r,a,s), using average_305 values\n");
      mri_morphed->c_r = -0.0950;
      mri_morphed->c_a = -16.5100;
      mri_morphed->c_s = 9.7500;
      mri_morphed->ras_good_flag = 1;
      // now we cache transform and thus we have to do the following whenever
      // we change direction cosines
      MRIreInitCache(mri_morphed);
    }
  else
    useVolGeomToMRI(&gcam->atlas, mri_morphed);

  return(mri_morphed) ;
}
double round(double x);
/*-------------------------------------------------------------------*/
MRI *ReverseMorph(char *subject){
  char *sd, tmpstr[2000];
  MRI *mri, *revmorph, *hit, *dist2;
  float cfov,rfov,sfov, c,r,s, cc,rr,ss;
  float d2,voxd2;
  int icc,irr,iss,voxhit,gcaerr;
  GCA_MORPH *gcam;
  MATRIX *vox2ras, *ras, *crs;

  sd = getenv("SUBJECTS_DIR");
  sprintf(tmpstr,"%s/%s/mri/orig.mgz",sd,subject);
  mri = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
  if(mri==NULL) return(NULL);

  sprintf(tmpstr,"%s/%s/mri/transforms/talairach.m3z",sd,subject);
  printf("Reading %s\n",tmpstr);
  gcam = GCAMread(tmpstr);
  if(gcam == NULL) {
    MRIfree(&mri);
    return(NULL);
  }

  revmorph = MRIallocSequence(mri->width,mri->height,mri->depth,MRI_FLOAT,3);
  hit = MRIallocSequence(mri->width,mri->height,mri->depth,MRI_FLOAT,1);
  dist2 = MRIallocSequence(mri->width,mri->height,mri->depth,MRI_FLOAT,1);

  vox2ras = MRIxfmCRS2XYZtkreg(mri);
  crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[4][1] = 1;
  ras = MatrixAlloc(4,1,MATRIX_REAL);

  cfov  = gcam->width  * gcam->spacing ; 
  rfov  = gcam->height * gcam->spacing ; 
  sfov  = gcam->depth  * gcam->spacing ; 

  // go through tal space
  for(c=0; c<cfov; c++){
    printf("c = %f\n",c);
    for(r=0; r<rfov; r++){
      for(s=0; s<sfov; s++){

	gcaerr = GCAMsampleMorph(gcam, c,r,s, &cc,&rr,&ss);
	if(!gcaerr){
	  if(cc < 0 || cc >= mri->width-1 ||
	     rr < 0 || rr >= mri->height-1 ||
	     ss < 0 || ss >= mri->depth-1 ) continue;
	  icc = round(cc);
	  irr = round(rr);
	  iss = round(ss);
	  d2 = pow(icc-cc,2) + pow(irr-rr,2) + pow(iss-ss,2);
	  voxhit = MRIgetVoxVal(hit,icc,irr,iss,0);
	  voxd2 = MRIgetVoxVal(dist2,icc,irr,iss,0); 
	  if(!voxhit || (voxhit && voxd2 > d2)){
	    MRIsetVoxVal(dist2,icc,irr,iss,0,d2);
	    crs->rptr[1][1] = c; 
	    crs->rptr[2][1] = r; 
	    crs->rptr[3][1] = s; 
	    ras = MatrixMultiply(vox2ras,crs,NULL);
	    
	    MRIsetVoxVal(revmorph,icc,irr,iss,0,ras->rptr[1][1]);
	    MRIsetVoxVal(revmorph,icc,irr,iss,1,ras->rptr[2][1]);
	    MRIsetVoxVal(revmorph,icc,irr,iss,2,ras->rptr[3][1]);
	  } // if vox hit
	  MRIsetVoxVal(hit,icc,irr,iss,0,1);
	} // gcaerr
      } // slice
    } // row 
  } // col

  MRIwrite(hit,"hit.mgh");
  MRIwrite(dist2,"dist2.mgh");
  MRIwrite(revmorph,"revmoph.mgh");

  return(revmorph);
}




