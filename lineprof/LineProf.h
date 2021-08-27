/**
 * @brief Class to interface with the laplace solver and line profiler
 *
 */

/*
 * Original Author: Martin Reuter
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
 
#ifndef _LineProf_h
#define _LineProf_h

#include <vector>
#include "Tracer.h"

class vtkPolyData;

/** \class LineProf
 * \brief Class to interface with the laplace solver and line profiler library
 * The library uses VTK ITK and Petsc to solve the laplace equation on a 2D
 * polygon (with 4 boundary segments).
 */
class LineProf
{

public:

  //! Petsc Initializer (call this at the beginning of your program once)
  static int InitializePetsc(bool bDoNotExitOnError = false);
  
  //! Petsc Finalizer (call this at the end of your program once)
  static int FinalizePetsc();

  //! Constructor from 2d points and 4 boundary segments
  LineProf(const std::vector < std::vector < double > >& points2d,
           const std::vector < int >& segment0,
           const std::vector < int >& segment1,
           const std::vector < int >& segmentL,
           const std::vector < int >& segmentR);
           
  //! Destructor
  ~LineProf() { if (_tracer) delete _tracer; };
    
  //! Laplace solver
  void solveLaplace(int paddingL, int paddingR,
                    double resolution,
		                int convergenceCriterion);

  //! Compute Profiles as array of lines (each line an array of 2d coordinates)
  std::vector < std::vector < std::vector < double > > >
  ComputeProfiles(double offset, double spacing);  

private:

  //! Samples points along midline
  std::vector < std::vector < double > >
  samplePointsMidline(double offset, double spacing);
    
  //! Computes isolines at levels in vec
  std::vector < std::vector < std::vector < double > > >
  ComputeIsolines(const std::vector < double >& vec, double x0, double y0);

  //! Converts 2D polygon to vtkPolyData
  vtkPolyData* GetPolyData();
  
  //! Converts 2D polygon to vtkPolyData
  std::vector < std::vector < double > > ConvertLine(const Tracer::LineType& line);
  
  double pointDistance(const std::vector < double >& pt1,
                       const std::vector < double >& pt2);
  
  std::vector < std::vector < double > > _points2d;
  std::vector < int > _segment0;
  std::vector < int > _segment1;
  std::vector < int > _segmentL;
  std::vector < int > _segmentR;
  Tracer* _tracer;

};


#endif
