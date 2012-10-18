/**
 * @file  LineProfTest.cpp
 * @brief Test example for line profile library
 *
 */

/*
 * Original Author: Martin Reuter, Oct. 17th ,2012
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2012/10/18 17:41:05 $
 *    $Revision: 1.3 $
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

#include <cmath>

#include "LineProf.h"



void createTestSegment(double startx, double starty, double endx, double endy, double dist,
                       std::vector < std::vector < double > > & points2d,
                       std::vector < int >& segment)
{
  unsigned int count = points2d.size();
  double length = std::sqrt((endx-startx)*(endx-startx) + (endy-starty)*(endy-starty));
  int num = int(length / dist);
  
  std::vector < double > p(2);
  
  for (int i=0;i<=num;i++)
  {
    double d = double(i) / double(num);
    p[0] = (1-d) * startx + d * endx;
    p[1] = (1-d) * starty + d * endy;
    points2d.push_back(p);
    segment.push_back(count);
    count++;
  }
    
}

int main(int argc, char *argv[])
{

  // You need to initialize Petsc at the beginning of your program:
  LineProf::InitializePetsc();


  // setup parameters
  double laplaceResolution = 0.5;
  double referenceSize     = 0.1;
  double resolution        = laplaceResolution * referenceSize;
  double distance          = resolution / 3.0;
  double profileSpacing    = 2.0;
  double spacing           = profileSpacing * referenceSize;


  // Here we setup the polygon and 4 boundary segments
  // as a rectangle from 10..20, 6..8
  std::vector < std::vector < double > > points2d;
  std::vector < int > segment0;
  std::vector < int > segment1;
  std::vector < int > segmentL;
  std::vector < int > segmentR;  
  createTestSegment(10,6,20,6,distance,points2d,segment0);
  createTestSegment(10,8,20,8,distance,points2d,segment1);
  createTestSegment(10,6,10,8,distance,points2d,segmentL);
  createTestSegment(20,6,20,8,distance,points2d,segmentR);
  

  // We initialize LineProf class by passing the polygon
  LineProf LP (points2d,segment0,segment1,segmentL,segmentR);
  
  
  // Next we solve the Laplace on the domain
  int paddingL       = 10;
  int paddingR       = 10;
  int convergence    = 8;
  LP.solveLaplace(paddingL,paddingR,resolution,convergence);


  // And finally compute line profiles
  int offset     = 10;
  std::vector < std::vector < std::vector < double > > > profiles;
  profiles = LP.ComputeProfiles(offset, spacing);


  // Petsc needs to be finalized at the end of your program
  LineProf::FinalizePetsc();
  
}

