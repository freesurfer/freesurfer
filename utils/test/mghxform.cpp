/**
 * @brief test xform name saving and reading
 *
 */
/*
 * Original Author: Y. Tosa
 *
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

extern "C"
{
#include "tags.h"
#include "matrix.h"
#include "mri.h"

const char *Progname="mghxform";
}

using namespace std;

void printLinearTransform(MRI *mri)
{
  MATRIX *mat = MatrixAlloc(4,4, MATRIX_REAL);
  for (int row = 1 ; row <= 3 ; row++)
  {
    *MATRIX_RELT(mat,row,1) = mri->linear_transform->m[0][row-1];
    *MATRIX_RELT(mat,row,2) = mri->linear_transform->m[1][row-1];
    *MATRIX_RELT(mat,row,3) = mri->linear_transform->m[2][row-1];
    *MATRIX_RELT(mat,row,4) = mri->linear_transform->m[3][row-1];
  }
  fprintf(stderr, "talairach transform\n");
  MatrixPrint(stdout, mat);
  MatrixFree(&mat);
}

int main(int argc, char *argv[])
{
  cout << "expanding the test data" << endl;
  system("gunzip -c orig.tar.gz | tar xvf - > /dev/null");
  // now we have talairach.xfm and orig cor files
  cout << "reading COR" << endl;
  MRI *mriCOR=MRIread((char*)"./orig");
  if (!mriCOR)
  {
    cerr << "could not read orig volume" << endl;
    exit(1);
    return -1;
  }
  printLinearTransform(mriCOR);
  // this should have read the xform
  cout << "writing mgh with xform info" << endl;
  MRIwrite(mriCOR, (char*)"./testxfm.mgh");
  MRIfree(&mriCOR);

  cout << "reading mgh with xform info" << endl;
  MRI *mriMGH = MRIread((char*)"./testxfm.mgh");
  cout << mriMGH->transform_fname << endl;
  printLinearTransform(mriMGH);
  if (strcmp(mriMGH->transform_fname, "./orig/../talairach.xfm"))
  {
    MRIfree(&mriMGH);

    cerr << "wrong filename for the transform" << endl;
    exit(1);
    return -1;
  }
  MRIfree(&mriMGH);
  system("rm -Rf ./orig/");
  system("rm ./testxfm.mgh");
  exit(0);
  return 0;
}


