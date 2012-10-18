/**
 * @file  LineProfTest.cpp
 * @brief Test for line profile library
 *
 */

/*
 * Original Author: Martin Reuter, Oct. 17th ,2012
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2012/10/18 14:33:42 $
 *    $Revision: 1.2 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "LineProf.h"

int main(int argc, char *argv[])
{

  LineProf::InitializePetsc();

  // setup polygon
  std::vector < std::vector < double > > points2d;
  std::vector < int > segment0;
  std::vector < int > segment1;
  std::vector < int > segmentL;
  std::vector < int > segmentR;
  
  int count = 0;
  std::vector < double > p(2);
  for (int i = -20; i<=20; i++)
  {
    p[0] = double(i);
    p[1] = 0.0;
    points2d.push_back(p);
    segment0.push_back(count);
    count++;
  }
  for (int i = -20; i<=20; i++)
  {
    p[0] = double(i);
    p[1] = 10.0;
    points2d.push_back(p);
    segment1.push_back(count);
    count++;
  }
  for (int i = 0; i<=10; i++)
  {
    p[0] = -20.0;
    p[1] = double(i);
    points2d.push_back(p);
    segmentL.push_back(count);
    count++;
  }
  for (int i = 0; i<=10; i++)
  {
    p[0] = 20.0;
    p[1] = double(i);
    points2d.push_back(p);
    segmentR.push_back(count);
    count++;
  }
  

  // init LineProf class
  LineProf LP (points2d,segment0,segment1,segmentL,segmentR);
  
  
  // solve Laplace
  int paddingL       = 10;
  int paddingR       = 10;
  double dresolution = 0.5;
  int convergence    = 8;
  LP.solveLaplace(paddingL,paddingR,dresolution,convergence);


  // compute line profiles
  int offset     = 10;
  double dspacing = 2;
  std::vector < std::vector < double > > referenceLine;
  std::vector < std::vector < std::vector < double > > > profiles;
  profiles = LP.ComputeProfiles(offset, dspacing, referenceLine);

  LineProf::FinalizePetsc();
  
}

