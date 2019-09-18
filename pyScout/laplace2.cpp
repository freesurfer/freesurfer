
#include <math.h>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

#include "laplace2.h"

#include "preprocessor.h"
#include "solver.h"
#include "tracer.h"


//----------------------------------


double sqr(double x) { return x*x; }

// double
// distance(const Tracer::PointType& pt1,
// 	 const Tracer::PointType& pt2)
// {
//   double dsum(0);
//   for(unsigned int ui =0; 
//       ui < 2; ++ui)
//     { dsum += sqr( pt2[ui]-pt1[ui] ); }
//   return std::sqrt(dsum);
// }

double 
pointDistance(const DoubleVector& pt1,
	      const DoubleVector& pt2)
{
 double dsum(0);
  for(unsigned int ui =0; 
      ui < 2; ++ui)
    { dsum += sqr( pt2[ui]-pt1[ui] ); }
  return std::sqrt(dsum);
}

//-----------------------------------------------------------

int
InitializePackage()
{
  PetscInitialize(0, 0,(char*)0,NULL);
  PetscMPIInt mpiSize;
  MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);
}

int
Finalize()
{
  PetscFinalize();
}

Tracer
ComputeSolution(WrapperPolyData* wrapper,
		int pad1, int pad2,
		double dresolution,
		int convergenceCriterion
		)
{
  vtkPolyData* polyData = GetPolyData(wrapper);

  //--------------------------
  // preprocessing
  // - create mask image
  // - initialize boundary conditions
  //
  Preprocessor pre;
  pre.SetInput( polyData );
  Preprocessor::SizeType padding;
  padding[0] = pad1;
  padding[1] = pad2;
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
  Tracer tracer;
  tracer.SetInputData(solver.GetOutputData());
  tracer.SetInputMask(mask);
  tracer.SetInputContours(polyData);
  
  return tracer;
}

WrapperLineType
ConvertLine(const Tracer::LineType& line)
{
  WrapperLineType outLine;
  outLine.reserve( line.size() );
  for( Tracer::LineType::const_iterator cit = line.begin();
       cit != line.end(); ++cit)
    {
      DoubleVector dvec;
      for(unsigned int ui=0; 
	  ui < 2; ++ui)
	dvec.push_back( (*cit)[ui] );

      outLine.push_back(dvec);
    }

  return outLine;
}

WrapperLineSetType
ComputeProfiles(int offset,
		double dspacing,
		WrapperLineType& referenceLine,
		const Tracer& tracer)
{
  typedef Tracer::LineType LineType;
  typedef std::vector<LineType> LineContainerType;
  WrapperLineSetType lc;

  // compute the points on the reference line
  WrapperLineType points;
  
  double ddist(0);

  for( int count(offset), maxVal(referenceLine.size()-offset);
       count < maxVal;
       ++count )
    {
      ddist += pointDistance( referenceLine[count-1], referenceLine[count] );

      if (ddist>dspacing)
	{
	  points.push_back( referenceLine[count] );
	  ddist = 0;
	}
      
    } // next count

  // compute the profiles for each of the points previously determined
  for(WrapperLineType::const_iterator cit = points.begin();
      cit != points.end(); ++cit)
      lc.push_back( ConvertLine(tracer.ComputeProfile((*cit)[0], (*cit)[1]) ) );
    
  return lc;
}

WrapperLineSetType
ComputeIsolines(const DoubleVector& vec,
		const Tracer& tracer,
		double x0, double y0)
{
  WrapperLineSetType lc;

  for( DoubleVector::const_iterator cit = vec.begin();
       cit != vec.end();
       ++cit )
    {
      Tracer::LineType bufLine 
	= tracer.ComputeIsoline(*cit);
      if ( sqr(x0- bufLine.front()[0]) +
	   sqr(y0- bufLine.front()[1]) >
	   sqr(x0- bufLine.back()[0]) +
	   sqr(y0- bufLine.back()[1]) )
	std::reverse( bufLine.begin(), bufLine.end() );
      lc.push_back( ConvertLine(bufLine) );
    }

  return lc;
}

