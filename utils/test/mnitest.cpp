/**
 * @brief testing three kinds of mni xfm reading
 *
 */
/*
 * Original Author: Yasanari Tosa
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
#include <iomanip>

extern "C"
{

#include "macros.h"
#include "matrix.h"
#include "transform.h"

  const char *Progname="mnitest";
}

using namespace std;

bool isEqual(MATRIX *m1, MATRIX *m2)
{
  for (int j=1; j < 5; ++j)
    for (int i=1; i < 5; ++i)
    {
      if (!FZERO(*MATRIX_RELT(m1, i, j) - *MATRIX_RELT(m2, i, j)))
      {
        cout << "(" << i << ", " << j << ")" << endl;
        cout << *MATRIX_RELT(m1, i, j) <<
        "      " << *MATRIX_RELT(m2, i, j) << endl;
        return false;
      }
    }
  cout.flush();
  return true;
}

void initialize_matrices
(MATRIX *m1, MATRIX *m2, MATRIX *m3, MATRIX *m4, MATRIX *m5)
{
// m1
// 0.99999822116498 -0.00856134979192902 -0.00330484398336356 -9.16348127558962
// 0.00849241373202879 0.999797376141696 -0.0203387865527882 -40.2488790680195
// 0.00347816160954997 0.0203098637039503 0.999828007058502 -3.17331194010443;
  *MATRIX_RELT(m1, 1, 1) = 0.99999822116498;
  *MATRIX_RELT(m1, 1, 2) = -0.00856134979192902;
  *MATRIX_RELT(m1, 1, 3) = -0.00330484398336356;
  *MATRIX_RELT(m1, 1, 4) = -9.16348127558962;
  *MATRIX_RELT(m1, 2, 1) = 0.00849241373202879;
  *MATRIX_RELT(m1, 2, 2) = 0.999797376141696;
  *MATRIX_RELT(m1, 2, 3) = -0.0203387865527882;
  *MATRIX_RELT(m1, 2, 4) = -40.2488790680195;
  *MATRIX_RELT(m1, 3, 1) = 0.00347816160954997;
  *MATRIX_RELT(m1, 3, 2) = 0.0203098637039503;
  *MATRIX_RELT(m1, 3, 3) = 0.999828007058502;
  *MATRIX_RELT(m1, 3, 4) = -3.17331194010443;
  *MATRIX_RELT(m1, 4, 1) = 0;
  *MATRIX_RELT(m1, 4, 2) = 0;
  *MATRIX_RELT(m1, 4, 3) = 0;
  *MATRIX_RELT(m1, 4, 4) = 1;
  // m2
  *MATRIX_RELT(m2, 1, 1) = *MATRIX_RELT(m1, 1, 1);
  *MATRIX_RELT(m2, 1, 2) = *MATRIX_RELT(m1, 1, 2);
  *MATRIX_RELT(m2, 1, 3) = *MATRIX_RELT(m1, 1, 3);
  *MATRIX_RELT(m2, 1, 4) = *MATRIX_RELT(m1, 1, 4);
  *MATRIX_RELT(m2, 2, 1) = *MATRIX_RELT(m1, 2, 1);
  *MATRIX_RELT(m2, 2, 2) = *MATRIX_RELT(m1, 2, 2);
  *MATRIX_RELT(m2, 2, 3) = *MATRIX_RELT(m1, 2, 3);
  *MATRIX_RELT(m2, 2, 4) = *MATRIX_RELT(m1, 2, 4);
  *MATRIX_RELT(m2, 3, 1) = *MATRIX_RELT(m1, 3, 1);
  *MATRIX_RELT(m2, 3, 2) = *MATRIX_RELT(m1, 3, 2);
  *MATRIX_RELT(m2, 3, 3) = *MATRIX_RELT(m1, 3, 3);
  *MATRIX_RELT(m2, 3, 4) = *MATRIX_RELT(m1, 3, 4);
  *MATRIX_RELT(m2, 4, 1) = 0;
  *MATRIX_RELT(m2, 4, 2) = 0;
  *MATRIX_RELT(m2, 4, 3) = 0;
  *MATRIX_RELT(m2, 4, 4) = 1;
// m3
// -0.253240725075032 -0.942293197668388 -0.218980969016388 -2.65098829180222
// 0.966925445826678 -0.25365952341587 -0.0266838601651706 20.2288466190752
// -0.0304025883166679 -0.218495711189758 0.975364191886973 -18.0512767925497;
  *MATRIX_RELT(m3, 1, 1) = -0.253240725075032;
  *MATRIX_RELT(m3, 1, 2) = -0.942293197668388;
  *MATRIX_RELT(m3, 1, 3) = -0.218980969016388;
  *MATRIX_RELT(m3, 1, 4) = -2.65098829180222;
  *MATRIX_RELT(m3, 2, 1) = 0.966925445826678;
  *MATRIX_RELT(m3, 2, 2) = -0.25365952341587;
  *MATRIX_RELT(m3, 2, 3) = -0.0266838601651706;
  *MATRIX_RELT(m3, 2, 4) = 20.2288466190752;
  *MATRIX_RELT(m3, 3, 1) = -0.0304025883166679;
  *MATRIX_RELT(m3, 3, 2) = -0.218495711189758;
  *MATRIX_RELT(m3, 3, 3) =  0.975364191886973;
  *MATRIX_RELT(m3, 3, 4) = -18.0512767925497;
  *MATRIX_RELT(m3, 4, 1) = 0;
  *MATRIX_RELT(m3, 4, 2) = 0;
  *MATRIX_RELT(m3, 4, 3) = 0;
  *MATRIX_RELT(m3, 4, 4) = 1;
// m4
// 1.15407013893127 -0.0289636701345444 0.0497409962117672 -5.76248216629028
// 0.0179894287139177 0.94682639837265 0.133943930268288 -22.0936622619629
// -0.0524738132953644 -0.158202826976776 1.1253559589386 -0.760713577270508;
  *MATRIX_RELT(m4, 1, 1) = 1.15407013893127;
  *MATRIX_RELT(m4, 1, 2) = -0.0289636701345444;
  *MATRIX_RELT(m4, 1, 3) = 0.0497409962117672;
  *MATRIX_RELT(m4, 1, 4) = -5.76248216629028;
  *MATRIX_RELT(m4, 2, 1) = 0.0179894287139177;
  *MATRIX_RELT(m4, 2, 2) = 0.94682639837265;
  *MATRIX_RELT(m4, 2, 3) = 0.133943930268288;
  *MATRIX_RELT(m4, 2, 4) = -22.0936622619629;
  *MATRIX_RELT(m4, 3, 1) = -0.0524738132953644;
  *MATRIX_RELT(m4, 3, 2) = -0.158202826976776;
  *MATRIX_RELT(m4, 3, 3) = 1.1253559589386;
  *MATRIX_RELT(m4, 3, 4) = -0.760713577270508;
  *MATRIX_RELT(m4, 4, 1) = 0;
  *MATRIX_RELT(m4, 4, 2) = 0;
  *MATRIX_RELT(m4, 4, 3) = 0;
  *MATRIX_RELT(m4, 4, 4) = 1;
// m5
// 0.999989       0.000015       -0.000046       -0.369853
// 0.000008       1.000032       -0.000005       3.526909
// -0.000012       -0.000062       1.000003       -11.687472;
  *MATRIX_RELT(m5, 1, 1) = 0.999989;
  *MATRIX_RELT(m5, 1, 2) = 0.000015;
  *MATRIX_RELT(m5, 1, 3) = -0.000046;
  *MATRIX_RELT(m5, 1, 4) = -0.369853;
  ///////////////////////////////////////
  *MATRIX_RELT(m5, 2, 1) = 0.000008;
  *MATRIX_RELT(m5, 2, 2) = 1.000032;
  *MATRIX_RELT(m5, 2, 3) = -0.000005;
  *MATRIX_RELT(m5, 2, 4) = 3.526909;
  /////////////////////////////////
  *MATRIX_RELT(m5, 3, 1) = -0.000012;
  *MATRIX_RELT(m5, 3, 2) = -0.000062;
  *MATRIX_RELT(m5, 3, 3) = 1.000003;
  *MATRIX_RELT(m5, 3, 4) = -11.687472;
  /////////////////////////////////////
  *MATRIX_RELT(m5, 4, 1) = 0;
  *MATRIX_RELT(m5, 4, 2) = 0;
  *MATRIX_RELT(m5, 4, 3) = 0;
  *MATRIX_RELT(m5, 4, 4) = 1;
}

const char *file1="./bruce.xfm";
const char *file2="./bruce2.xfm";
const char *file3="./david.xfm";
const char *file4="./tosa.xfm";
const char *file5="./tosa2.xfm";

using namespace std;

int main(int argc, char *argv[])
{
  MATRIX *m1, *m2, *m3, *m4, *m5;

  m1 = MatrixAlloc(4, 4, MATRIX_REAL);
  m2 = MatrixAlloc(4, 4, MATRIX_REAL);
  m3 = MatrixAlloc(4, 4, MATRIX_REAL);
  m4 = MatrixAlloc(4, 4, MATRIX_REAL);
  m5 = MatrixAlloc(4, 4, MATRIX_REAL);

  initialize_matrices(m1, m2, m3, m4, m5);

  LTA *lta1, *lta2, *lta3, *lta4, *lta5;
  cout << "////////////////////////////////////" << endl;
  cout << "Old way of reading -----------------" << endl;
  cout << "////////////////////////////////////" << endl;
  cout << "bruce.xfm" << endl;
  lta1 = LTAread((char*)file1);
  LTAprint(stdout, lta1);
  // old way always convert it to vox-to-vox
#if 0
  if (!isEqual(lta1->xforms[0].m_L, m1))
  {
    cerr << "lta1 read wrong" << endl;
    return -1;
  }
#endif
  cout << "bruce2.xfm" << endl;
  lta2 = LTAread((char*)file2);
  LTAprint(stdout, lta2);
#if 0
  if (!isEqual(lta2->xforms[0].m_L, m2))
  {
    cerr << "lta1 read wrong" << endl;
    return -1;
  }
#endif
  cout << "david.xfm" << endl;
  lta3 = LTAread((char*)file3);
  LTAprint(stdout, lta3);
#if 0
  if (!isEqual(lta3->xforms[0].m_L, m3))
  {
    cerr << "lta1 read wrong" << endl;
    return -1;
  }
#endif
  cout << "tosa.xfm" << endl;
  lta4 = LTAread((char*)file4);
  LTAprint(stdout, lta4);
#if 0
  if (!isEqual(lta4->xforms[0].m_L, m4))
  {
    cerr << "lta1 read wrong" << endl;
    return -1;
  }
#endif
  LTAfree(&lta1);
  LTAfree(&lta2);
  LTAfree(&lta3);
  LTAfree(&lta4);

  cout << "////////////////////////////////////" << endl;
  cout << "New way of reading -----------------" << endl;
  cout << "////////////////////////////////////" << endl;
  cout << "bruce.xfm" << endl;
  cout << "////////////////////////////////////" << endl;
  lta1 = LTAreadEx(file1);
  LTAprint(stdout, lta1);
  if (!isEqual(lta1->xforms[0].m_L, m1))
  {
    cerr << "lta1 read wrong" << endl;
    return -1;
  }
  cout << "bruce2.xfm" << endl;
  cout << "////////////////////////////////////" << endl;
  lta2 = LTAreadEx(file2);
  LTAprint(stdout, lta2);
  if (!isEqual(lta2->xforms[0].m_L, m2))
  {
    cerr << "lta1 read wrong" << endl;
    return -1;
  }
  cout << "david.xfm" << endl;
  cout << "////////////////////////////////////" << endl;
  lta3 = LTAreadEx(file3);
  LTAprint(stdout, lta3);
  if (!isEqual(lta3->xforms[0].m_L, m3))
  {
    cerr << "lta1 read wrong" << endl;
    return -1;
  }
  cout << "tosa.xfm" << endl;
  cout << "////////////////////////////////////" << endl;
  lta4 = LTAreadEx(file4);
  LTAprint(stdout, lta4);
  if (!isEqual(lta4->xforms[0].m_L, m4))
  {
    cerr << "lta1 read wrong" << endl;
    return -1;
  }
  cout << "tosa2.xfm" << endl;
  cout << "////////////////////////////////////" << endl;
  lta5 = LTAreadEx(file5);
  LTAprint(stdout, lta5);
  MatrixPrint(stdout, m5);
  if (!isEqual(lta5->xforms[0].m_L, m5))
  {
    cerr << "lta1 read wrong" << endl;
    return -1;
  }

  LTAfree(&lta1);
  LTAfree(&lta2);
  LTAfree(&lta3);
  LTAfree(&lta4);
  LTAfree(&lta5);
  MatrixFree(&m1);
  MatrixFree(&m2);
  MatrixFree(&m3);
  MatrixFree(&m4);
  MatrixFree(&m5);
}
