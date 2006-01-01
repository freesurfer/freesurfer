#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

#define NEW_VERSION_MAGIC_NUMBER  16777215 // defined in mrisurf.c

int IDisCurv(char *curvfile);
MRI *MRISreadCurvAsMRI(char *curvfile);


char *Progname = "dngtester";
char *subject=NULL, *hemi=NULL;
char *SUBJECTS_DIR = NULL;
char *curvfile = NULL;
MRI *mri;
MRIS *mris;
char tmpstr[2000];
int err;

/*----------------------------------------*/
int main(int argc, char **argv)
{
  subject = argv[1];
  hemi = argv[2];
  curvfile = argv[3];

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  sprintf(tmpstr,"%s/%s/surf/%s.orig",SUBJECTS_DIR,subject,hemi);
  printf("%s\n",tmpstr);
  mris = MRISread(tmpstr);

  sprintf(tmpstr,"%s.%s",hemi,curvfile);
  err = MRISreadCurvatureFile(mris, tmpstr);
  if(err) printf("curv err\n");
  mri = MRIcopyMRIS(NULL,mris,0,"curv");
  MRIwrite(mri,"./lh.curv0.mgh");

  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,curvfile);
  printf("%s\n",tmpstr);
  mri = MRISreadCurvAsMRI(tmpstr);
  MRIwrite(mri,"./lh.curv2.mgh");


  return(0);
}

/*----------------------------------------*/
int IDisCurv(char *curvfile)
{
  int magno;
  FILE *fp;

  fp = fopen(curvfile,"r");
  if (fp==NULL)return(0); // return quietly
  fread3(&magno,fp);
  fclose(fp) ;
  if(magno != NEW_VERSION_MAGIC_NUMBER) return(0);
  return(1);
}
/*----------------------------------------*/
MRI *MRISreadCurvAsMRI(char *curvfile)
{
  int    magno,k,vnum,fnum, vals_per_vertex ;
  float  curv;
  FILE   *fp;
  MRI *curvmri;

  if(!IDisCurv(curvfile)) return(NULL);
  
  fp = fopen(curvfile,"r");
  fread3(&magno,fp);
  
  vnum = freadInt(fp);
  fnum = freadInt(fp);
  vals_per_vertex = freadInt(fp) ;
  if (vals_per_vertex != 1){
    fclose(fp) ;
    printf("ERROR: MRISreadCurvAsMRI: %s, vals/vertex %d unsupported\n",
	   curvfile,vals_per_vertex);
    return(NULL);
  }

  curvmri = MRIalloc(vnum,1,1,MRI_FLOAT);
  for (k=0;k<vnum;k++){
    curv = freadFloat(fp) ;
    MRIsetVoxVal(curvmri,k,0,0,0,curv);
  }
  fclose(fp);

  return(curvmri) ;
}
