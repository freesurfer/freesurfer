/**
 * @brief Test example for line profile library
 *
 */

/*
 * Original Author: Martin Reuter, Oct. 17th ,2012
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

#include <math.h>

#include "LineProf.h"


/** Creates a line segment between start(xy) and end(xy)
  where the distance between points is approx dist.
  New points including start and end are appended to points2d
  and indices are appended to segment. */
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

/** Checks the profile lines for the rectangle.
  They should be vertical lines (fixed x coord) where
  the y coordinate is equally spaced between 6 and 8.*/
void checkProfiles(const std::vector < std::vector < std::vector < double > > >& profiles)
{

  if (profiles.size() != 40)
  {
    std::cerr <<"ERROR: expecting 40 profiles, but got " << profiles.size() << " !" << std::endl;
    exit(1);
  }

  double eps = 0.000001;
  double ystart = 6.0;
  double ydelta = 0.1;
  double xstart = 11.1;
  double xdelta = 0.2;
  double xp = 0.0;
  for (unsigned int p=0; p<profiles.size(); p++)
  {
    if (profiles[p].size() != 21)
    {
      std::cerr <<"ERROR: expecting 21 points in profile " << p << ", but got " << profiles[p].size() << " !" << std::endl;
      exit(1);
    }
    xp = profiles[p][0][0];
    //std::cout << " p " << p << " x " << xp << std::endl;
    if (fabs(xp - (p*xdelta+xstart)) > eps)
    {
        std::cerr << "ERROR: Profile " <<p <<"  x-coordinate difference!" << std::endl;
        std::cerr <<setprecision(9) << p*xdelta+xstart << " != "<< xp << std::endl;
        exit(1);      
    }
  
    for (unsigned int i=0; i<profiles[p].size(); i++)
    {
      if (profiles[p][i].size() != 2)
      {
        std::cerr <<"ERROR: expecting 2 coords in profile " << p << " point " << i << " , but got " << profiles[p][i].size() << " !" << std::endl;
        exit(1);
      }
      if (fabs(profiles[p][i][0] - xp) > eps)
      {
        std::cerr <<"ERROR: Profile " <<p <<"  point " << i << " x-coordinate difference!" << std::endl;
        std::cerr <<setprecision(9) << xp << " != "<< profiles[p][i][0] << std::endl;
        exit(1);
      }
      
      if (fabs(profiles[p][i][1] - (i*ydelta+ystart)) > eps)
      {
        std::cerr << "ERROR: Profile " <<p <<"  point " << i << " y-coordinate difference!" << std::endl;
        std::cerr <<setprecision(9) << i*ydelta+ystart << " != "<< profiles[p][i][1] << std::endl;
        exit(1);
      }
    }
  }
  std::cout << std::endl << "Profile check successfull ! " << std::endl<< std::endl;
  
}

void printProfiles(const std::vector < std::vector < std::vector < double > > >& profiles)
{
  for (unsigned int i = 0;i<profiles.size();i++)
  {
    cout << "Profile " << i << endl;
    for (unsigned int j = 0 ; j<profiles[i].size();j++)
    {
      cout << profiles[i][j][0] << " " << profiles[i][j][1] << endl;
    }
    cout << endl;
  }
}


int main(int argc, char *argv[])
{

  // You need to initialize Petsc at the beginning of your program:
  LineProf::InitializePetsc();


  // setup parameters:
  
  // voxel size (smallest voxel side length), needed to convert voxel to RAS
  // usually you get this from the image!
  double referenceSize     = 0.1;
  
  // we want to solve the laplace at twice the voxel resolution (half a voxel)
  // Unless you want to sample really densly, I think this can be keept at 0.5
  double laplaceResolution = 0.5;
   
  // spacing in voxels between line profiles (here every second voxel)
  double profileSpacing    = 2.0;

  // offsets (necessary as the line profiles are not straight close to the
  // side boundaries)
  int offset      = 10;
  int paddingL    = 10;
  int paddingR    = 10;
  
  // Parameter for the convergence of laplace solver (can probably be keept at 8)
  int convergence = 8;
  
  // compute RAS distances from the above using referenceSize:
  // RAS resolution for laplace solver
  double ras_resolution     = laplaceResolution * referenceSize; 
  // RAS sampling distance on input lines
  double ras_distance       = ras_resolution / 3.0;
  // RAS spacing between line profiles
  double ras_spacing        = profileSpacing * referenceSize;
  // RAS offset
  double ras_offset         = offset * referenceSize;


  // Here we setup the polygon and 4 boundary segments
  // as a rectangle from 10..20, 6..8 (ras coords)
  std::vector < std::vector < double > > points2d;
  std::vector < int > segment0;
  std::vector < int > segment1;
  std::vector < int > segmentL;
  std::vector < int > segmentR;  
  createTestSegment(10, 6, 20, 6, ras_distance, points2d, segment0);
  createTestSegment(10, 8, 20, 8, ras_distance, points2d, segment1);
  createTestSegment(10, 6, 10, 8, ras_distance, points2d, segmentL);
  createTestSegment(20, 6, 20, 8, ras_distance, points2d, segmentR);
  

  // We initialize LineProf class by passing the polygon
  LineProf LP (points2d, segment0, segment1, segmentL, segmentR);
  
  
  // Next we solve the Laplace on the domain
  LP.solveLaplace(paddingL, paddingR, ras_resolution, convergence);


  // And finally compute line profiles
  std::vector < std::vector < std::vector < double > > > profiles;
  profiles = LP.ComputeProfiles(ras_offset, ras_spacing);


  // We check the profiles (should be vertical lines in this rectangle)
  //printProfiles(profiles);
  checkProfiles(profiles);


  // Petsc should be finalized at the end of your program
  LineProf::FinalizePetsc();
  
}

