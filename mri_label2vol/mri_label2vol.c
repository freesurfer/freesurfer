/*
  Name:    mri_label2vol.c
  Author:  Douglas N. Greve 
  email:   analysis-bugs@nmr.mgh.harvard.edu
  Date:    2/27/02
  Purpose: Converts a label to a segmentation volume.
  $Id: mri_label2vol.c,v 1.3 2004/08/04 17:42:26 greve Exp $
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "error.h"
#include "diag.h"
#include "proto.h"

#include "mri.h"
#include "mrisutils.h"
#include "MRIio_old.h"
#include "mri_identify.h"
#include "mri2.h"
#include "matrix.h"
#include "version.h"
#include "registerio.h"
#include "resample.h"


#define PROJ_TYPE_NONE 0
#define PROJ_TYPE_ABS  1
#define PROJ_TYPE_FRAC 2

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
static int  checkhemi(char *hemi);
static int get_proj_type_id(char *projtype);
static int get_crs(MATRIX *Tras2vox, float x, float y, float z,
		   int *c, int *r, int *s, MRI *vol);
static int is_surface_label(LABEL *label);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_label2vol.c,v 1.3 2004/08/04 17:42:26 greve Exp $";
char *Progname = NULL;

char *LabelList[100];
int nlabels = 0;

char *TempVolId = NULL;
char *RegMatFile = NULL;
double FillThresh = 0.0;
double ProjDelta = .1;
double ProjStart = 0;
double ProjStop = 0;
char *ProjType = NULL;
char *OutVolId = NULL;
char *HitVolId = NULL;
char *subject = NULL;
char *hemi    = NULL;
char *SurfId  = "white";

MRI_SURFACE *Surf=NULL;
MRI *OutVol=NULL, *TempVol=NULL, *HitVol=NULL;

int debug;

LABEL *srclabel;
MATRIX *R, *Tvox2ras, *Tras2vox;
double ProjDepth;
int   ProjTypeId;
int   DoProj = 0;
char  *SUBJECTS_DIR = NULL;
char  *thicknessname = "thickness";
char  fname[1000];
double TempVoxVol;
double LabelVoxVol = 1;
double nHitsThresh;


/*---------------------------------------------------------------*/
int main(int argc, char **argv)
{
  int  nargs, nthlabel, float2int, err, nthpoint, vtxno;
  float ipr,bpr,intensity;
  char *regsubject;
  float x,y,z;
  int c,r,s, oob, nhits, nhitsmax, nhitsmax_label;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, 
      "$Id: mri_label2vol.c,v 1.3 2004/08/04 17:42:26 greve Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0];

  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);
  printf("%s\n",vcid);

  // Load the template volume
  TempVol = MRIreadHeader(TempVolId,MRI_VOLUME_TYPE_UNKNOWN);
  if(TempVol == NULL){
    printf("ERROR: reading %s header\n",TempVolId);
    exit(1);
  }
  Tvox2ras = MRIxfmCRS2XYZtkreg(TempVol); 
  Tras2vox = MatrixInverse(Tvox2ras,NULL);
  printf("Template RAS-to-Vox: --------\n");
  MatrixPrint(stdout,Tras2vox);

  TempVoxVol = TempVol->xsize * TempVol->ysize * TempVol->zsize;
  nHitsThresh = FillThresh*TempVoxVol/LabelVoxVol;
  printf("Template Voxel Volume: %g\n",TempVoxVol);
  printf("nHits Thresh: %g\n",nHitsThresh);

  if(RegMatFile != NULL){
    // Load registration matrix
    err = regio_read_register(RegMatFile, &regsubject, &ipr, &bpr, 
			      &intensity, &R, &float2int);
    if(err) exit(1);
    printf("RegMat: --------\n");
    MatrixPrint(stdout,R);
    MatrixMultiply(Tras2vox,R,Tras2vox);
  }
  printf("Label RAS-to-Vox: --------\n");
  MatrixPrint(stdout,Tras2vox);

  if(DoProj){
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if(SUBJECTS_DIR==NULL){
      fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
    // Load the surface used for projection
    Surf = MRISloadSurfSubject(subject,hemi,SurfId,SUBJECTS_DIR);
    if(Surf == NULL){
      printf("ERROR: could not load surface.\n");
      exit(1);
    }
    // Load the thickness for projection along the normal
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,
	    hemi,thicknessname);
    printf("Reading thickness %s\n",fname);
    MRISreadCurvatureFile(Surf, fname);
  }

  // Create hit volume based on template, one frame for each label
  HitVol = MRIallocSequence(TempVol->width, TempVol->height, 
          TempVol->depth, MRI_INT, nlabels );
  if(HitVol == NULL){
    printf("ERROR: could not alloc hit volume\n");
    exit(1);
  }
  MRIcopyHeader(TempVol,HitVol);
  HitVol->nframes = nlabels;

  // Go through each label
  for(nthlabel = 0; nthlabel < nlabels; nthlabel++){
    printf("nthlabel = %d\n",nthlabel);
    srclabel = LabelRead(NULL, LabelList[nthlabel]);
    if(srclabel == NULL){
      printf("ERROR reading %s\n",LabelList[nthlabel]);
      exit(1);
    }
    if(DoProj && !is_surface_label(srclabel)){
      printf("ERROR: label %s is not a surface label.\n",
	     LabelList[nthlabel]);
      exit(1);
    }
    // Go through each point in the label 
    for(nthpoint = 0; nthpoint < srclabel->n_points; nthpoint++){
      if(debug) printf("  nthpoint = %d\n",nthpoint);

      if(DoProj){ // Project along the surface normal
	vtxno = srclabel->lv[nthpoint].vno;
	ProjDepth = ProjStart;

	while(ProjDepth <= ProjStop){

	  if(ProjTypeId == PROJ_TYPE_ABS)
	    ProjNormDist(&x,&y,&z,Surf,vtxno,ProjDepth);
	  if(ProjTypeId == PROJ_TYPE_FRAC)
	    ProjNormFracThick(&x,&y,&z,Surf,vtxno,ProjDepth);
	  oob = get_crs(Tras2vox,x,y,z,&c,&r,&s,TempVol);
	  if(debug) printf("   ProjDepth %g   %g %g %g (%g)  %d %d %d   %d\n",
		 ProjDepth,x,y,z,Surf->vertices[vtxno].curv,c,r,s,oob);
	  if(oob) continue; // Out of the volume
	  MRIIseq_vox(HitVol,c,r,s,nthlabel) ++;
	  // Accumulate hit volume
	  ProjDepth += ProjDelta;
	  if(ProjDelta == 0) break; // only do once

	} // end loop through projection depths
      }// end Do Projection

      else{ // Simply compute crs from the label xyz
	x = srclabel->lv[nthpoint].x;
	y = srclabel->lv[nthpoint].y;
	z = srclabel->lv[nthpoint].z;
	oob = get_crs(Tras2vox,x,y,z,&c,&r,&s,TempVol);
	if(debug) printf("   %g %g %g   %d %d %d   %d\n",x,y,z,c,r,s,oob);
	if(oob) continue; // Out of the volume
	MRIIseq_vox(HitVol,c,r,s,nthlabel) ++;
      }
    } // end loop over label points

    LabelFree(&srclabel) ;
  } // End loop over labels


  if(HitVolId != NULL) MRIwrite(HitVol,HitVolId);

  // Create output volume based on template, but use 1 frame short
  OutVol = MRIallocSequence(TempVol->width, TempVol->height, 
          TempVol->depth, MRI_INT, 1);
  if(OutVol == NULL){
    printf("ERROR: could not alloc output volume\n");
    exit(1);
  }
  MRIcopyHeader(TempVol,OutVol);
  OutVol->nframes = 1;

  // Threshold hit volumes and set outvol to nthlabel+1
  for(c=0; c < OutVol->width; c++){
    for(r=0; r < OutVol->height; r++){
      for(s=0; s < OutVol->depth; s++){

	nhitsmax = 0;
	nhitsmax_label = -1;
	for(nthlabel = 0; nthlabel < nlabels; nthlabel++){
	  nhits = (int)MRIgetVoxVal(HitVol,c,r,s,nthlabel);
	  //nhits = MRIIseq_vox(HitVol,c,r,s,nthlabel);
	  if(nhits <= nHitsThresh) continue;
	  if(nhitsmax < nhits){
	    nhitsmax = nhits;
	    nhitsmax_label = nthlabel;
	  }
	}
	MRIIseq_vox(OutVol,c,r,s,0) = nhitsmax_label + 1;

      }
    }
  }

  // Save out volume
  MRIwrite(OutVol,OutVolId);

  printf("done \n");
  return(0);
}
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if(argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while(nargc > 0){

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;

    else if (!strcmp(option, "--label")){
      if(nargc < 1) argnerr(option,1);
      LabelList[nlabels] = pargv[0];
      nlabels++;
      nargsused = 1;
    }
    else if (!strcmp(option, "--temp")){
      if(nargc < 1) argnerr(option,1);
      TempVolId = pargv[0];
      nargsused = 1;
    }

    else if (!strcmp(option, "--reg")){
      if(nargc < 1) argnerr(option,1);
      RegMatFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--fillthresh")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&FillThresh);
      nargsused = 1;
    }
    else if (!strcmp(option, "--labvoxvol")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&LabelVoxVol);
      nargsused = 1;
    }
    else if (!strcmp(option, "--proj")){
      if(nargc < 4) argnerr(option,4);
      ProjType = pargv[0];
      sscanf(pargv[1],"%lf",&ProjStart);
      sscanf(pargv[2],"%lf",&ProjStop);
      sscanf(pargv[3],"%lf",&ProjDelta);
      ProjTypeId = get_proj_type_id(ProjType);
      if(ProjStart > ProjStop){
	printf("ERROR: projection start must be <= stop\n");
	exit(1);
      }
      if(ProjDelta == 0 && (ProjStart != ProjStop)){
	printf("ERROR: cannot spec projection delta = 0.\n");
	exit(1);
      }
      if(ProjDelta < 0){
	printf("ERROR: projection delta must be > 0\n");
	exit(1);
      }
      DoProj = 1;
      nargsused = 4;
    }
    else if (!strcmp(option, "--subject")){
      if(nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--hemi")){
      if(nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      checkhemi(hemi);
      nargsused = 1;
    }
    else if (!strcmp(option, "--hit")){
      if(nargc < 1) argnerr(option,1);
      HitVolId = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--o")){
      if(nargc < 1) argnerr(option,1);
      OutVolId = pargv[0];
      nargsused = 1;
    }

    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(singledash(option))
	fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void)
{
  printf("USAGE: mri_label2vol\n") ;
  printf("\n");
  printf("   --label labelid <--label labelid>  \n");
  printf("   --temp tempvolid : template volume\n");
  printf("   --reg regmat : VolXYZ = R*LabelXYZ\n");
  printf("   --fillthresh thresh : between 0 and 1 (def 0)\n");
  printf("   --labvoxvol voxvol : volume of each label point (def 1mm3)\n");
  printf("\n");
  printf("   --proj type start stop delta\n");
  printf("   --subject subjectid\n");
  printf("   --hemi hemi\n");
  printf("\n");
  printf("   --o volid : output volume\n");
  printf("   --hits hitvolid : each frame is nhits for a label\n");
  printf("\n");
  printf("   --version : print version and exit\n");
  printf("   --help \n");
  printf("\n");
  printf("%s\n",vcid);
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(

"
Help Outline:
  - SUMMARY
  - ARGUMENTS
  - RESOLVING MULTI-LABEL AMBIGUITIES
  - CHECKING YOUR RESULTS
  - EXAMPLES
  - KNOWN BUGS
  - SEE ALSO

SUMMARY

Converts a label or a set of labels into a volume. For a single lable,
the volume will be binary: 1 where the label is and 0 where it is not.
For multiple labels, the volume will be 0 where no labels were found
otherwise the value will the the label number. For a voxel to be
declared part of a label, it must have enough hits in the voxel and
it must have more hits than any other label.

ARGUMENTS

--label labelfile <--label labelfile>

Enter the name of the label file. For multiple labels, use multiple
--label flags. Labels can be created manually with tkmedit and
tksurfer or automatically from a subcortical segmentation or cortical
annotation. Labels are simple text files. The first line is a header.
Each following line contains data with 5 columns. The first column is
the vertex number of the label point. The next 3 columns are the X, Y,
and Z of the point. The last can be ignored. If the label is not a 
surface-based label, then the vertex number will be -1.


--temp tempvolid

Template volume. The output volume will have the same size and geometry
as the template. Template must have geometry information (ie, direction
cosines and voxel sizes). Required.

--reg regmatfile

tkregister-style registration matrix (see tkregister2 --help) which maps
the XYZ of the label to the XYZ of the template volume. If not specified,
then the identity is assumed.

--fillthresh thresh

Relative threshold which the number hits in a voxel must exceed for
the voxel to be considered a candidate for membership in the label. A
'hit' is when a label point falls into a voxel. thresh is a value
between 0 and 1 indicating the fraction of the voxel volume that must
be filled by label points. The voxel volume is determined from the
template. It is assumed that the each label point represents a voxel
1mm3 in size (which can be changed with --labvoxvol). So, the actual
number of hits needed to exceed threshold is
thresh*TempVoxVol/LabVoxVol. The default value is 0, meaning that any
label point that falls into a voxel makes that voxel a candidate for
membership in the label.  Note: a label must also have the most hits
in the voxel before that voxel will actually be assigned to the label
in the volume. Note: the label voxel volume becomes a little ambiguous
for surface labels, particularly when they are 'filled in' with
projection.

--labvoxvol voxvol

Volume covered by each label point. Default is 1mm3. This only affects
the fill threshold (--fillthresh). Note: the label voxel volume
becomes a little ambiguous for surface labels, particularly when they
are projected.

--proj type start stop delta

Project the label along the surface normal. type can be abs or frac. 
abs means that the start, stop, and delta are measured in mm. frac
means that start, stop, and delta are relative to the thickness at
each vertex. The label definition is changed to fill in label
points in increments of delta from start to stop. Requires subject
and hemi in order to load in a surface and thickness. Uses the
white surface. The label MUST have been defined on the surface.

--subject subjectid

FREESURFER subject identifier. Needed when using --proj.

--hemi hemi

Hemisphere to use loading the surface for --proj. Legal values are
lh and rh.

--o volid

Single frame output volume in which each voxel will have the number of
the label to which it is assigned (or 0 for no label). The label
number is the order in which it appears on the command-line.  Takes
any format accepted by mri_convert (eg, spm, analyze, bshort, mgh).

--hits hitvolid

Hit volume. This is a multi-frame volume, with one frame for each
label. The value at each voxel for a given frame is the number of hits
that voxel received for that label. This is mostly good as a debugging
tool, but you could use it to implement your own multi-label
arbitration routine. Or you could binarize to have each label
represented separately. Takes any format accepted by mri_convert (eg,
spm, analyze, bshort, mgh).

RESOLVING MULTI-LABEL AMBIGUITIES

When there are multiple lables, it is possible that more than one
label will map to a single voxel in the output volume. When this
happens, the voxel is assigned to the label with the most label
points in that voxel. Note that the voxel must still pass the 
fill threshold test in order to be considered part of the label.

CHECKING YOUR RESULTS

It is very important to check that the conversion of the label to the
volume was done correctly. It may be that it is way off or it could be
off around the edges. This is particularly true for surface-based
labels or when converting a label to a low-resolution space.
To check the result, load the orig volume into tkmedit. The orig
volume should be in the label space. Load the mri_label2vol output
volume as an overlay; this makes the labeled voxels appear as
'activity'.  Finally, load the label itself. You should see the label
(in green) sitting on top of the 'activity' of the labeled volume.
See EXAMPLE 1 for an example.


EXAMPLES

1. Convert a label into a binary mask in the functional space; require
that a functional voxel be filled at least 50%% by the label:

mri_label2vol 
  --label lh-avg_central_sulcus.label 
  --temp f_000.bshort 
  --reg register.dat 
  --fillthresh .5 
  --o cent-lh_000.bshort

To see how well the label is mapped into the functional volume, run

tkmedit bert orig 
  -overlay ./cent-lh_000.bshort 
  -overlay-reg ./register.dat -fthresh .5 -fmid 1

Then load the label with File->Label->LoadLabel. The label should
overlap with the overlay. The overlap will not be perfect but it
should be very close.

2. Convert a surface label into a binary mask in the functional space.
Fill in all the cortical gray matter. Require that a functional voxel
be filled at least 30%% by the label:

mri_label2vol 
  --label lh-avg_central_sulcus.label 
  --temp f_000.bshort 
  --reg register.dat 
  --fillthresh .3 
  --proj frac 0 1 .1 
  --subject bert --hemi lh
  --o cent-lh_000.bshort

3. Convert a surface label into a binary mask in the functional space.
Sample a 1mm ribbon 2mm below the gray/white surface. Do
not require a fill threshold:

mri_label2vol 
  --label lh-avg_central_sulcus.label 
  --temp f_000.bshort 
  --reg register.dat 
  --proj abs -3 -2 .1 
  --subject bert --hemi lh
  --o cent-lh_000.bshort

4. Convert two labels into a volume in the same space as the labels:

mri_label2vol 
  --label lh-avg_central_sulcus.label 
  --label lh-avg_calcarine_sulcus.label 
  --temp $SUBJECTS_DIR/bert/orig
  --o cent_calc.img

The voxels corresponding to lh-avg_central_sulcus.label will have a of 
value of 1 whereas those assigned to lh-avg_calcarine_sulcus.label will
have a value of 2.

KNOWN BUGS

1. When the output type is COR, all the voxels will be zero. The work-around
is to save it as some other type, then use mri_convert with --no_rescale 1
to convert it to COR.

SEE ALSO

mri_label2label, mri_cor2label, mri_annotation2label, mri_mergelabels,
tkregister2, mri_convert, tkmedit, tksurfer.

http://surfer.nmr.mgh.harvard.edu/docs/tkmedit_guide.html
http://surfer.nmr.mgh.harvard.edu/docs/tksurfer_doc.html

"


);


  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if(n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(nlabels == 0){
    printf("ERROR: you must spec at least one label\n");
    exit(1);
  }
  if(OutVolId == NULL){
    printf("ERROR: no output specified\n");
    exit(1);
  }
  if(DoProj){
    if(subject == NULL){
      printf("ERROR: subject needed in order to load surface for projection\n");
      exit(1);
    }
    if(hemi == NULL){
      printf("ERROR: hemi needed in order to load surface for projection\n");
      exit(1);
    }
  }
  if(subject != NULL && !DoProj)
    printf("INFO: subject not needed, igorning.\n");
  if(hemi != NULL && !DoProj)
    printf("INFO: hemi not needed, igorning.\n");

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  int n;
  fprintf(fp,"Number of labels: %d\n",nlabels);
  for(n=0;n<nlabels;n++)
    fprintf(fp,"%s\n",LabelList[n]);
  fprintf(fp,"Template Volume: %s\n",TempVolId);
  fprintf(fp,"Outut Volume: %s\n",OutVolId);
  fprintf(fp,"Registration File: %s\n",RegMatFile);
  fprintf(fp,"Fill Threshold: %g\n",FillThresh);
  fprintf(fp,"Label Vox Vol:  %g\n",LabelVoxVol);
  fprintf(fp,"ProjType:       %s\n",ProjType);
  fprintf(fp,"ProjTypeId:     %d\n",ProjTypeId);
  fprintf(fp,"ProjStart:      %g\n",ProjStart);
  fprintf(fp,"ProjStop:       %g\n",ProjStop);
  fprintf(fp,"ProjDelta:      %g\n",ProjDelta);
  fprintf(fp,"Subject:  %s\n",subject);
  fprintf(fp,"Hemi:     %s\n",hemi);

  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*------------------------------------------------------------*/
/* Exits if something is wrong with hemi. */
static int checkhemi(char *hemi)
{
  if(! strcmp(hemi,"lh")) return(0);
  if(! strcmp(hemi,"rh")) return(0);
  printf("ERROR: hemi = %s, must be lh or rh\n",hemi);
  exit(1);
}
/*------------------------------------------------------------*/
static int get_proj_type_id(char *projtype)
{
  if(projtype == NULL)          return(PROJ_TYPE_NONE);
  if(! strcmp(projtype,"none")) return(PROJ_TYPE_NONE);
  if(! strcmp(projtype,"abs"))  return(PROJ_TYPE_ABS);
  if(! strcmp(projtype,"frac")) return(PROJ_TYPE_FRAC);
  printf("ERROR: projtype = %s, must be abs or frac\n",projtype);
  exit(1);
}
/*------------------------------------------------------------*/
// Computes the nearest col, row, slice of the xyz (RAS) point.
// Returns 1 if the point is out of the volume. Should really
// do this in one of the util libraries.
static int get_crs(MATRIX *Tras2vox, float x, float y, float z,
		   int *c, int *r, int *s, MRI *vol)
{
  float cf, rf, sf;

  cf = Tras2vox->rptr[1][1]*x + 
       Tras2vox->rptr[1][2]*y +
       Tras2vox->rptr[1][3]*z +
       Tras2vox->rptr[1][4];

  rf = Tras2vox->rptr[2][1]*x + 
       Tras2vox->rptr[2][2]*y +
       Tras2vox->rptr[2][3]*z +
       Tras2vox->rptr[2][4];

  sf = Tras2vox->rptr[3][1]*x + 
       Tras2vox->rptr[3][2]*y +
       Tras2vox->rptr[3][3]*z +
       Tras2vox->rptr[3][4];

  *c = rint(cf);
  *r = rint(rf);
  *s = rint(sf);

  if(vol == NULL) return(0);

  if(*c < 0 || *c >= vol->width)  return(1);
  if(*r < 0 || *r >= vol->height) return(1);
  if(*s < 0 || *s >= vol->depth)  return(1);

  return(0);
}
/*------------------------------------------------------------*/
static int is_surface_label(LABEL *label)
{
  int nthpoint,vtxno;

  for(nthpoint = 0; nthpoint < label->n_points; nthpoint++){
    vtxno = label->lv[nthpoint].vno;
    if(vtxno == -1) return(0);
  }

  return(1);
}
