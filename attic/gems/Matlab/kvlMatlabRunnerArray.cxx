#include "kvlMatlabRunnerArray.h"
#include "kvlReadCroppedImage.h"
#include "kvlGetImageBuffer.h"
#include "kvlSetImageBuffer.h"
#include "kvlSmoothImageBuffer.h"
#include "kvlReadMeshCollection.h"
#include "kvlGetMesh.h"
#include "kvlRasterizeAtlasMesh.h"
#include "kvlGetAlphasInMeshNodes.h"
#include "kvlSetAlphasInMeshNodes.h"
#include "kvlSmoothMeshCollection.h"
#include "kvlScaleMeshCollection.h"
#include "kvlWriteImage.h"
#include "kvlCreateImage.h"
#include "kvlClear.h"
#include "kvlSetMaximumNumberOfThreads.h"
#include "kvlScaleMesh.h"
#include "kvlSmoothMesh.h"
#include "kvlReadImage.h"
#include "kvlGetTransformMatrix.h"
#include "kvlGetMeshNodePositions.h"
#include "kvlSetMeshCollectionPositions.h"
#include "kvlWriteMeshCollection.h"
#include "kvlTransformMeshCollection.h"
#include "kvlSetKOfMeshCollection.h"
#include "kvlCreateTransform.h"
#include "kvlEvaluateMeshPositionWithEntropy.h"
#include "kvlSetMeshNodePositions.h"
#include "kvlChangeK.h"
#include "kvlCreateRGBImage.h"
#include "kvlWriteRGBImage.h"
#include "kvlGetCroppedRegion.h"
#include "kvlGetCostAndGradientCalculator.h"
#include "kvlGetAverageAtlasMeshPositionCostAndGradientCalculator.h"
#include "kvlEvaluateMeshPosition.h"
#include "kvlGetOptimizer.h"
#include "kvlStepOptimizer.h"
#include "kvlCreateMeshCollection.h"
#include "kvlDrawJacobianDeterminant.h"


namespace kvl
{

MatlabRunnerArray::Pointer MatlabRunnerArray::m_Instance = 0;

//
//
//
MatlabRunnerArray
::MatlabRunnerArray()
{
  m_Array.push_back( ReadCroppedImage::New().GetPointer() );
  m_Array.push_back( GetImageBuffer::New().GetPointer() );
  m_Array.push_back( SetImageBuffer::New().GetPointer() );
  m_Array.push_back( ReadMeshCollection::New().GetPointer() );
  m_Array.push_back( GetMesh::New().GetPointer() );
  m_Array.push_back( RasterizeAtlasMesh::New().GetPointer() );
  m_Array.push_back( SmoothImageBuffer::New().GetPointer() );
  m_Array.push_back( GetAlphasInMeshNodes::New().GetPointer() );
  m_Array.push_back( SetAlphasInMeshNodes::New().GetPointer() );
  m_Array.push_back( SmoothMeshCollection::New().GetPointer() );
  m_Array.push_back( ScaleMeshCollection::New().GetPointer() );
  m_Array.push_back( WriteImage::New().GetPointer() );
  m_Array.push_back( CreateImage::New().GetPointer() );
  m_Array.push_back( Clear::New().GetPointer() );
  m_Array.push_back( SetMaximumNumberOfThreads::New().GetPointer() );
  m_Array.push_back( ScaleMesh::New().GetPointer() );
  m_Array.push_back( SmoothMesh::New().GetPointer() );
  m_Array.push_back( ReadImage::New().GetPointer() );
  m_Array.push_back( GetTransformMatrix::New().GetPointer() );
  m_Array.push_back( GetMeshNodePositions::New().GetPointer() );
  m_Array.push_back( SetMeshCollectionPositions::New().GetPointer() );
  m_Array.push_back( WriteMeshCollection::New().GetPointer() );
  m_Array.push_back( TransformMeshCollection::New().GetPointer() );
  m_Array.push_back( SetKOfMeshCollection::New().GetPointer() );
  m_Array.push_back( CreateTransform::New().GetPointer() );
  m_Array.push_back( EvaluateMeshPositionWithEntropy::New().GetPointer() );
  m_Array.push_back( SetMeshNodePositions::New().GetPointer() );
  m_Array.push_back( ChangeK::New().GetPointer() );
  m_Array.push_back( CreateRGBImage::New().GetPointer() );
  m_Array.push_back( WriteRGBImage::New().GetPointer() );
  m_Array.push_back( GetCroppedRegion::New().GetPointer() );
  m_Array.push_back( GetCostAndGradientCalculator::New().GetPointer() );
  m_Array.push_back( GetAverageAtlasMeshPositionCostAndGradientCalculator::New().GetPointer() );
  m_Array.push_back( EvaluateMeshPosition::New().GetPointer() );
  m_Array.push_back( GetOptimizer::New().GetPointer() );
  m_Array.push_back( StepOptimizer::New().GetPointer() );
  m_Array.push_back( CreateMeshCollection::New().GetPointer() );
  m_Array.push_back( DrawJacobianDeterminant::New().GetPointer() );
}


//
// Return the single instance of the MatlabRunnerArray
//
MatlabRunnerArray::Pointer
MatlabRunnerArray
::New()
{
  
   if ( !MatlabRunnerArray::m_Instance )
    {
    MatlabRunnerArray::m_Instance = new MatlabRunnerArray;

   // Remove extra reference from construction.
    MatlabRunnerArray::m_Instance->UnRegister();
    }

  return MatlabRunnerArray::m_Instance;
 
}


//
//
//
bool 
MatlabRunnerArray
::Run( const std::string& runnerName, 
       int nlhs, mxArray* plhs[],
       int nrhs, const mxArray* prhs[] )
{
  
  bool executedCorrectRunner = false;
  
  for ( std::vector< MatlabRunner::Pointer >::const_iterator it = m_Array.begin();
        it != m_Array.end(); ++ it )
    {
      
    if ( !runnerName.compare( it->GetPointer()->GetNameOfClass() ) )
      {
      //std::cout << "Found a runner with name " << runnerName 
      //          << " : running it!" << std::endl;
      it->GetPointer()->Run( nlhs, plhs, nrhs, prhs );
      executedCorrectRunner = true;
      break;
      }
    }
  
  return executedCorrectRunner;
}



} // end namespace kvl


