#include <stdio.h>
#include "mri.h"

char *Progname;

int main(int argc, char *argv[])
{
	MRI *mri;
	Progname=argv[0];
  printf("Generate a test volume.\n");
	
	mri=MRIalloc(20,20,20,MRI_UCHAR);

	MRIvox(mri,10,10,10)=255;
	MRIvox(mri,11,10,10)=255;
	MRIvox(mri,11,11,10)=255;
	MRIvox(mri,11,11,11)=255;
	MRIvox(mri,10,11,11)=255;
	MRIvox(mri,10,10,11)=255;


	MRIwrite(mri,argv[1]);

  return 0;
}
