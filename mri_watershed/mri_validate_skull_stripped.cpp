
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include "mrisurf.h"
#include "utils.h"
#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "minc.h"
#include "region.h"
#include "machine.h"
#include "fio.h"
#include "mri_identify.h"
#include "mrisurf.h"
#include "fmriutils.h"
#include "gcsa.h"
#include "fsgdf.h"
#include "icosahedron.h"
#include "gca.h"
#include "gcamorph.h"
#include "DICOMRead.h"
#include "fsenv.h"
#include "qdecutils.h"
#include "dti.h"
#include "registerio.h"
#include "timer.h"
#include "evschutils.h"
#include "matrix.h"
#include "matfile.h"
#include "randomfields.h"

static void ErrorValidation(const char *string) {
  fprintf(stdout, "\nValidation Error: %s\n",string) ;
  exit(1) ;
}

const char *Progname = "Validation";
MRI  *mri_ref, *mri_tst;
int  width, height, depth, i, j, k, consistent=1;
int  nb_union, falsep, falsen;
float probap, proban, Ec, c, vox_ref, vox_tst;

int main(int argc, char **argv)
{
  if(argc != 4	){
    printf("USAGE: tntester mri_reference mri_test weight (>1) \n");
    exit(1);
  }
   
  c = atof(argv[3]);
  if (c<1)
    ErrorValidation("\n Problem in dimensions \n");
   
  //printf("Reading source \n");
  mri_ref=MRIread(argv[1]);  
  mri_tst=MRIread(argv[2]);  
  
  if (mri_ref->width == mri_tst->width && mri_ref->height == mri_tst->height && mri_ref->depth == mri_tst->depth){
    depth = mri_ref->width;
    height = mri_ref->height; 
    width= mri_ref->depth;
    }
  else
    ErrorValidation("\n Problem in dimensions \n");
   
  nb_union = 0;  falsep = 0;  falsen = 0;
  
  for ( k=0; k<depth; k++)
    for ( j=0; j<height; j++)
      for ( i=0; i<width; i++){
      vox_ref = MRIgetVoxVal(mri_ref, i, j, k ,0);      
      vox_tst = MRIgetVoxVal(mri_tst, i, j, k ,0);
      if (vox_ref>0 && vox_tst>0 && vox_ref!=vox_tst)
          consistent=0;
      if (vox_ref>0){
        nb_union ++;
        if (vox_tst==0)
	  falsen ++;  // miss
	}
      else
        if(vox_tst>0){
	  nb_union ++;
	  falsep ++;
	  }  
  }

  probap = (float)falsep/nb_union;
  proban = (float)falsen/nb_union;
  Ec = (probap+c*proban)/(1+c);
  
 // if (!consistent)
 //   printf(" The two images don't have the same intensity at some positive points \n\n");

  
//  printf("RESULT of Validation:\n");
//  printf("E(c)= %f  ", Ec);
 printf("   pf   =  %f\n", probap);
  printf("   pm   =  %f\n", proban);
//  printf("   J    =  %f\n", 1-probap-proban);

return(0);
}
