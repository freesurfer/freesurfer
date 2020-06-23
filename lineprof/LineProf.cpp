#include <math.h>

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

#include "LineProf.h"
#include "Tracer.h"
#include "Preprocessor.h"
#include "PetscSolver.h"

LineProf::LineProf(const std::vector < std::vector < double > >& points2d,
                   const std::vector < int >& segment0,
                   const std::vector < int >& segment1,
                   const std::vector < int >& segmentL,
                   const std::vector < int >& segmentR)
                  : _points2d(points2d), _segment0(segment0), _segment1(segment1),
                    _segmentL(segmentL), _segmentR(segmentR), _tracer(NULL)
{
// maybe some error check here? e.g. if points are really 2d?
}

void LineProf::solveLaplace(int paddingL, int paddingR, double dresolution, int convergenceCriterion)
{
  vtkPolyData* polyData = GetPolyData();

  //--------------------------
  // preprocessing
  // - create mask image
  // - initialize boundary conditions
  //
  Preprocessor pre;
  pre.SetInput( polyData );
  Preprocessor::SizeType padding;
  padding[0] = paddingL;
  padding[1] = paddingR;
  pre.SetPadding( padding );

  Preprocessor::ResolutionType resolution;  
  resolution.Fill(dresolution);
  pre.SetResolution(resolution);

  pre.Update();

  //-------------------------
  // solver step
  //
  PetscSolver solver;
  solver.SetInputData( pre.GetOutputData() );
  solver.SetInputMask( pre.GetOutputLabel(),
           pre.GetMaskInsideValue(),
           pre.GetMaskZeroValue()
           );
  solver.Update(convergenceCriterion);

  // post-processing
  
  // filter the mask
  typedef PetscSolver::MaskImageType MaskImageType;
  MaskImageType::Pointer mask = MaskImageType::New();
  mask->SetRegions( solver.GetOutputMask()->GetRequestedRegion() );
  mask->Allocate();
  mask->CopyInformation( solver.GetOutputMask() );
  typedef itk::ImageRegionConstIterator<MaskImageType> MaskConstIteratorType;
  typedef itk::ImageRegionIterator<MaskImageType> MaskIteratorType;
  MaskConstIteratorType inputIter(solver.GetOutputMask(), 
          solver.GetOutputMask()->GetRequestedRegion());
  MaskIteratorType outputIter(mask, mask->GetRequestedRegion() );
  for( inputIter.GoToBegin(), outputIter.GoToBegin();
       !outputIter.IsAtEnd(); 
       ++inputIter, ++outputIter )
    if ( inputIter.Get() == solver.GetInsideValue() )
      outputIter.Set(1);
    else outputIter.Set(0);

  // initialize tracer object
  //
  if (_tracer) delete _tracer;
  _tracer = new Tracer;
  _tracer->SetInputData(solver.GetOutputData());
  _tracer->SetInputMask(mask);
  _tracer->SetInputContours(polyData);
  
}

std::vector < std::vector < std::vector < double > > >
LineProf::ComputeProfiles(double offset,  double dspacing)
{
  if (_tracer == NULL)
  {
    std::cerr << "Line Profile Error: call solveLaplace before trying to compute profiles!" << std::endl;
    return std::vector < std::vector < std::vector < double > > > (0);
  }

  std::vector < std::vector < std::vector < double > > > lc;

//  // compute isoline at 0.5 level (in the middle)
//  std::vector < double > dv(1,0.5);
//  std::vector < std::vector < std::vector < double > > > vml;
//  vml = ComputeIsolines(dv,_points2d[_segment0[0]][0], _points2d[_segment0[0]][1]);
//  const std::vector < std::vector < double > > & referenceLine = vml[0];
//
//   // compute the points on the reference line
//   std::vector < std::vector < double > > points;
// 
//   
//   double ddist(0);
//   for( int count(offset), maxVal(referenceLine.size()-offset); count < maxVal; ++count )
//   {
//     ddist += pointDistance( referenceLine[count-1], referenceLine[count] );
//     if (ddist>dspacing)
//     {
//       points.push_back( referenceLine[count] );
//       ddist = 0;
//     }
//   } // next count

   // compute the points on the mid iso level (0.5)
   std::vector < std::vector < double > > points = samplePointsMidline(offset,dspacing);
  

  // compute the profiles for each of the points previously determined
  for(std::vector < std::vector < double > >::const_iterator cit = points.begin();
      cit != points.end(); ++cit)
      lc.push_back( ConvertLine(_tracer->ComputeProfile((*cit)[0], (*cit)[1]) ) );
    
  return lc;
}

//* Samples num points along the midline, leaving offset space on both sides */
std::vector < std::vector < double > >
LineProf::samplePointsMidline(double offset,  double dspacing)
{
  //std::cout << "LineProf::samplePointsMidline(offset " << offset << " ,  spacing " << dspacing << " )" <<std::endl;
  std::vector < std::vector < double > > points;

  // compute isoline at 0.5 level (in the middle)
  std::vector < double > dv(1,0.5);
  std::vector < std::vector < std::vector < double > > > vml;
  vml = ComputeIsolines(dv,_points2d[_segment0[0]][0], _points2d[_segment0[0]][1]);
  const std::vector < std::vector < double > > & rl = vml[0];

  // compute length of reference line
  unsigned int s = rl.size();
  double l = 0.0;
  for (unsigned int i = 1; i < s ; i++)
  {
    l += pointDistance(rl[i], rl[i-1]);
    //std::cout << rl[i][0]  << "  " << rl[i][1] << std::endl;
  }
  //std::cout << " length : " << l << std::endl;

  // length available after removing offsets (on both sides):
  double sl = l-2.0 * offset;
  if (sl <= 0) return points;

  // if sl > 0 we will place at least one point, 
  // this tells us how many more we will have:
  int num = (int)(sl/dspacing);
  
  // unused space:
  double remain = sl - num * dspacing;
  
  // where to place the first (at offset plus half the unused space)
  double lnext = offset + 0.5 * remain;
  double lstop = l - offset;
  double dist  = 0.0;
  double lastl = 0.0;
  std::vector < double > point(2);
  for (unsigned int i = 1; i < s ; i++)
  {
    //std::cout << "lnext: " << lnext << std::endl;
    if (lnext > lstop) break;
    dist = pointDistance(rl[i], rl[i-1]);
    //std::cout << "lastl: " << lastl << "  dist : " << dist << std::endl;
    if (lastl + dist > lnext) // subsample segment
    {
      double d = (lnext-lastl)/dist;
      //cout << " d :"<< d << endl;
      //cout << " x(i-1): "<<rl[i-1][0] << "  x(i): " << rl[i][0] <<endl;
      point[0] = rl[i-1][0] * (1-d) + rl[i][0] *d;
      point[1] = rl[i-1][1] * (1-d) + rl[i][1] *d;
      //cout << " add point " << point[0] << " " << point[1] << endl;
      points.push_back(point);
      lnext += dspacing; // try next point
      i--; // check this segment again
    }
    else
    {
      // switch to next
      lastl += dist;
    }
  }
  //std::cout << " computed number of points: " << points.size() << std::endl;
  return points;
}

/** Computes isolines at the levels specified in vec.
  x0 and y0 specify a point close to the starting point of the lines 
  to ensure consistent orientation */
std::vector < std::vector < std::vector < double > > >
LineProf::ComputeIsolines(const std::vector < double >& vec, double x0, double y0)
{
  //std::cout << "LineProf::computeIsolines (vec, x0 = " << x0 << " , y0 = " << y0 << " )" <<std::endl;
  if (_tracer == NULL)
  {
    std::cerr << "Line Profile Error: call solveLaplace before trying to compute isolines!" << std::endl;
    return std::vector < std::vector < std::vector < double > > > (0);
  }
  
  std::vector < std::vector < std::vector < double > > > lc;

  for( std::vector < double >::const_iterator cit = vec.begin();
       cit != vec.end();
       ++cit )
  {
    Tracer::LineType bufLine = _tracer->ComputeIsoline(*cit);
    double d1 = x0- bufLine.front()[0];
    double d2 = y0- bufLine.front()[1];
    double d3 = x0- bufLine.back()[0];
    double d4 = y0- bufLine.back()[1];
    if ( d1*d1 + d2 *d2 > d3 *d3 + d4 * d4 )
      std::reverse( bufLine.begin(), bufLine.end() );
    lc.push_back( ConvertLine(bufLine) );
  }

  return lc;
}


vtkPolyData* LineProf::GetPolyData() 
{
  vtkPolyData*  outData   = vtkPolyData::New();
  vtkPoints*    outPoints = vtkPoints::New();
  vtkCellArray* outLines  = vtkCellArray::New();

  for(unsigned int ui(0), noPts(_points2d.size()); ui < noPts; ++ui)
  {
    assert(_points2d[ui].size()==2);
    outPoints->InsertNextPoint(_points2d[ui][0],_points2d[ui][1],0.0);
  }

  unsigned int noPts;
  
  noPts = _segment0.size();
  outLines->InsertNextCell( noPts );
  for(unsigned int ui(0); ui < noPts; ++ui )
    outLines->InsertCellPoint( _segment0[ui] );
  
  noPts = _segment1.size();
  outLines->InsertNextCell( noPts );
  for(unsigned int ui(0); ui < noPts; ++ui )
    outLines->InsertCellPoint( _segment1[ui] );
  
  noPts = _segmentL.size();
  outLines->InsertNextCell( noPts );
  for(unsigned int ui(0); ui < noPts; ++ui )
    outLines->InsertCellPoint( _segmentL[ui] );
    
  noPts = _segmentR.size();
  outLines->InsertNextCell( noPts );
  for(unsigned int ui(0); ui < noPts; ++ui )
    outLines->InsertCellPoint( _segmentR[ui] );

  outData->SetPoints(outPoints);
  outData->SetLines(outLines);
  return outData;
}

std::vector < std::vector < double > >
LineProf::ConvertLine(const Tracer::LineType& line)
{
  std::vector < std::vector < double > > outLine;
  outLine.reserve( line.size() );
  
  for( Tracer::LineType::const_iterator cit = line.begin();
       cit != line.end(); ++cit)
  {
    std::vector < double >  dvec;
    for(unsigned int ui=0; ui < 2; ++ui)
      dvec.push_back( (*cit)[ui] );

    outLine.push_back(dvec);
  }

  return outLine;
}

double LineProf::pointDistance(const std::vector < double >& pt1,
                               const std::vector < double >& pt2)
{
  double dsum(0.0);
  double d;
  for(unsigned int ui =0; ui < 2; ++ui)
  {
    d = pt2[ui]-pt1[ui];
    dsum += d*d;
  }
  return std::sqrt(dsum);
}

int LineProf::InitializePetsc(bool DoNotExitOnError)
{
  PetscInitializeNoArguments();
  PetscMPIInt mpiSize;
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
  Tracer::DoNotExitOnError = DoNotExitOnError;
  return 1;
}

int LineProf::FinalizePetsc()
{
  PetscFinalize();
  return 1;
}

