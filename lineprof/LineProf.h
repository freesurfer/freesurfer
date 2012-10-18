#ifndef _LineProf_h
#define _LineProf_h

#include <vector>
#include "Tracer.h"

class vtkPolyData;

class LineProf
{

public:
  
  LineProf(const std::vector < std::vector < double > >& points2d,
           const std::vector < int >& segment0,
           const std::vector < int >& segment1,
           const std::vector < int >& segmentL,
           const std::vector < int >& segmentR);
  ~LineProf() { if (_tracer) delete _tracer; };
  
  void solveLaplace(int paddingL, int paddingR,
                    double dresolution,
		                int convergenceCriterion);

  std::vector < std::vector < std::vector < double > > >
  ComputeProfiles(int offset,	double dspacing, const std::vector < std::vector < double > >& referenceLine);
  
  std::vector < std::vector < std::vector < double > > >
  ComputeIsolines(const std::vector < double >& vec, double x0, double y0);

  static int InitializePetsc();
  static int FinalizePetsc();

private:
  
  vtkPolyData* GetPolyData();
  
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
