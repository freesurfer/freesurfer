//
// regdat2xfm.cpp
//

#include <iostream>
#include <iomanip>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>

extern "C" {
#include "error.h"
#include "mri.h"
#include "transform.h"
#include "talairachex.h"

  char *Progname="regdat2xfm";
}

using namespace std;

void usage()
{
  cout << "regdat2xfm <srcvol> <targetvol> <regdat> <xfm>" << endl;
}

MATRIX *getRAS2RegRAS(MRI *mri)
{
  MATRIX *i2regRAS=MatrixAlloc(4, 4, MATRIX_REAL);
  // stupid one
  *MATRIX_RELT(i2regRAS, 1, 1) = -mri->xsize ; 
  *MATRIX_RELT(i2regRAS, 1, 4) = (mri->xsize)*(mri->width)/2.;
  *MATRIX_RELT(i2regRAS, 2, 3) = -mri->ysize ; 
  *MATRIX_RELT(i2regRAS, 2, 4) = -(mri->zsize)*(mri->depth)/2;
  *MATRIX_RELT(i2regRAS, 3, 2) = mri->zsize ;  
  *MATRIX_RELT(i2regRAS, 3, 4) = (mri->zsize)*(mri->height)/2;
  *MATRIX_RELT(i2regRAS, 4, 4) = 1 ;

  MATRIX *r_to_i = extract_r_to_i(mri);

  MATRIX *res = MatrixMultiply(i2regRAS, r_to_i, NULL);
  MatrixFree(&i2regRAS);
  MatrixFree(&r_to_i);
  return res;
}

int main(int argc, char *argv[])
{
  if (argc < 4)
  {
    usage();
    return -1;
  }
  MRI *src = MRIread(argv[1]);
  MRI *dst = MRIread(argv[2]);
  // get the mri2frmi
  MATRIX *regSrc2Dest = 0;
  int type = TransformFileNameType(argv[3]);
  if (type == REGISTER_DAT)
  {
    fMRI_REG *reg = StatReadRegistration(argv[3]);
    MatrixCopy(reg->mri2fmri, regSrc2Dest);
    StatFreeRegistration(&reg);
  }
  // get ras2RegRAS
  MATRIX *ras2RegRAS = getRAS2RegRAS(src);
  MATRIX *ras2RegRAS2 = getRAS2RegRAS(dst);
  MATRIX *regRAS2RAS = MatrixInverse(ras2RegRAS2, NULL);

  //       src  ->  RAS
  //        |1       |
  //        V        V
  //       src  ->  regRAS
  //        |        |
  //        |        | (mri2fmri)
  //        V        V
  //       dst  ->  regRAS
  //        | 1      |
  //        V        V
  //       dst  ->  RAS
  //////////////////////////////////////////
  MATRIX *tmp = MatrixMultiply(regSrc2Dest, ras2RegRAS, NULL);
  MATRIX *ras2RAS = MatrixMultiply(regRAS2RAS, tmp, NULL);

  MatrixFree(&regSrc2Dest);
  MatrixFree(&ras2RegRAS2);
  MatrixFree(&ras2RegRAS);
  MatrixFree(&tmp);

  // now write MNI xfm
  FILE *fp = fopen(argv[4], "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "could not open file %s",argv[4]));

  fprintf(fp, "MNI Transform File\n") ;
  // now saves src and dst in comment line
  fprintf(fp, "%sVolume: %s\n","%", argv[2]); // dst
  fprintf(fp, "%sVolume: %s\n","%", argv[1]); // src
  fprintf(fp, "\n") ;
  fprintf(fp, "Transform_Type = Linear;\n") ;
  fprintf(fp, "Linear_Transform =\n") ;

  for (int row = 1 ; row <= 3 ; row++)
  {
    fprintf(fp, "      %f       %f       %f       %f",
	    *MATRIX_RELT(ras2RAS,row,1), *MATRIX_RELT(ras2RAS,row,2), 
	    *MATRIX_RELT(ras2RAS,row,3), *MATRIX_RELT(ras2RAS,row,4)) ;
    if (row == 3)
      fprintf(fp, ";") ;
    fprintf(fp, "\n") ;
  }
  fclose(fp);
  MatrixFree(&ras2RAS);
  MRIfree(&src);
  MRIfree(&dst);
}


