#include "kvlAtlasMeshDeformationLevenbergMarquardt2.h"

#include "itkAffineTransform.h"


namespace kvl
{

//
//
//
AtlasMeshDeformationLevenbergMarquardt
::AtlasMeshDeformationLevenbergMarquardt()
{
  
}


//
//
//
AtlasMeshDeformationLevenbergMarquardt
::~AtlasMeshDeformationLevenbergMarquardt()
{
}



//
//
//
AtlasPositionGradientContainerType::Pointer
AtlasMeshDeformationLevenbergMarquardt
::GetStep( float lambda, bool verbose ) const
{

  itk::TimeProbe  timeProbe;
  timeProbe.Start();


  // Solve diagonal linear system
  AtlasPositionGradientContainerType::Pointer  step = AtlasPositionGradientContainerType::New();
  AtlasMesh::PointDataContainer::ConstIterator  paramIt 
            = this->GetFragmentProcessor().GetMesh()->GetPointData()->Begin();
  FragmentProcessorType::GradientContainerType::const_iterator  gradIt = m_Gradient.begin();
  FragmentProcessorType::DiagonalHessianContainerType::const_iterator  hessIt = m_Hessian.begin();
            
  for ( ; paramIt != this->GetFragmentProcessor().GetMesh()->GetPointData()->End(); 
        ++paramIt, ++gradIt, ++hessIt )
    {
    AtlasPositionGradientType  entry( 0.0f );
  
    if ( paramIt.Value().m_CanMoveX )
      {
      entry[ 0 ] = -( ( gradIt->second )[ 0 ] ) / ( ( 1 + lambda ) * ( ( hessIt->second )[ 0 ] ) );
      }

    if ( paramIt.Value().m_CanMoveY )
      {
      entry[ 1 ] = -( ( gradIt->second )[ 1 ] ) / ( ( 1 + lambda ) * ( ( hessIt->second )[ 1 ] ) );
      }

    if ( paramIt.Value().m_CanMoveZ )
      {
      entry[ 2 ] = -( ( gradIt->second )[ 2 ] ) / ( ( 1 + lambda ) * ( ( hessIt->second )[ 2 ] ) );
      }
      
    // Remember that we now have the step in the coordinate system parallell to the mesh. In order
    // to go back into image grid coordinate system, we need to multiply by T (from  x = T * u )
    if ( this->GetMeshToImageTransform() )
      {
      entry = ( this->GetMeshToImageTransform()->GetMatrix() ) * entry;
      }


    // std::cout << "Gradient in point with index " << pointIt.Index() << ": [ "
    //          << entry[ 0 ] << "   " << entry[ 1 ] << "   " << entry[ 2 ] << "]" << std::endl;

    step->InsertElement( paramIt.Index(), entry );

    } // End loop over all points
  

  timeProbe.Stop();
  std::cout << "Time taken to solve the Levenberg-Marquardt system of equations: " << timeProbe.GetMeanTime() << std::endl;

  return step;

}



} // end namespace kvl

