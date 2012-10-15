/**
 * @file  kvlAtlasMeshSmoother.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
 *    $Revision: 1.4 $
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
#include "kvlAtlasMeshSmoother.h"


namespace kvl
{


//
//
//
AtlasMeshSmoother
::AtlasMeshSmoother()
{
  m_MeshCollection = 0;
  m_Sigma = 1.0f;
}



//
//
//
AtlasMeshSmoother
::~AtlasMeshSmoother()
{
}




//
//
//
void
AtlasMeshSmoother
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
}




//
//
//
AtlasMeshCollection::Pointer
AtlasMeshSmoother
::GetSmoothedMeshCollection()
{

  // Sanity check on input
  if ( !m_MeshCollection )
  {
    itkExceptionMacro( << "No mesh collection set!" );
  }

  std::cout << "Smoothing mesh collection..." << std::flush;

  // Construct container to hold smoothed alphas
  AtlasMesh::PointDataContainer::Pointer  smoothedParameters = AtlasMesh::PointDataContainer::New();


  // Loop over all points
  for ( AtlasMesh::PointsContainer::ConstIterator centerRefPosIt = m_MeshCollection->GetReferencePosition()->Begin();
        centerRefPosIt != m_MeshCollection->GetReferencePosition()->End(); ++centerRefPosIt )
  {
    // Loop over all points in the mesh
    AtlasAlphasType  smoothedAlphas( m_MeshCollection->GetPointParameters()->Begin().Value().m_Alphas.Size() );
    smoothedAlphas.Fill( 0.0f );
    AtlasMesh::PointsContainer::ConstIterator  posIt = m_MeshCollection->GetReferencePosition()->Begin();
    AtlasMesh::PointDataContainer::ConstIterator  paramIt = m_MeshCollection->GetPointParameters()->Begin();
    for ( ; posIt != m_MeshCollection->GetReferencePosition()->End(); ++posIt, ++paramIt )
    {
      // Calculate the distance to the center point
      const  float  sqrDistance = pow( posIt.Value()[ 0 ] - centerRefPosIt.Value()[ 0 ], 2 ) +
                                  pow( posIt.Value()[ 1 ] - centerRefPosIt.Value()[ 1 ], 2 ) +
                                  pow( posIt.Value()[ 2 ] - centerRefPosIt.Value()[ 2 ], 2 );

      // If distance too big, ignore this point
      if ( sqrDistance > pow( 3 * m_Sigma, 2 ) )
      {
        continue;
      }


      // Add weighted contribution
      const float  weight = exp( -sqrDistance / pow( m_Sigma, 2 ) / 2 );
      smoothedAlphas += weight * paramIt.Value().m_Alphas;
    } // End loop over all points in the mesh

    // Normalize the smoothed alphas
    smoothedAlphas /= smoothedAlphas.sum();

    // Include the average in the smoothed alphas container
    PointParameters  smoothedEntry = m_MeshCollection->GetPointParameters()->ElementAt( centerRefPosIt.Index() );
    smoothedEntry.m_Alphas = smoothedAlphas;
    smoothedParameters->InsertElement( centerRefPosIt.Index(), smoothedEntry );
  } // End loop over all points

  // Construct a mesh collection to return
  AtlasMeshCollection::Pointer  result = AtlasMeshCollection::New();
  result->SetPointParameters( smoothedParameters );
  result->SetCells( m_MeshCollection->GetCells() );
  result->SetReferencePosition(  m_MeshCollection->GetReferencePosition() );
  result->SetK( m_MeshCollection->GetK() );
  result->SetPositions( m_MeshCollection->GetPositions() );

  std::cout << "done!" << std::endl;
  //m_MeshCollection->Write( "original.txt" );
  //result->Write( "smoothed.txt" );

  return result;
}



} // end namespace kvl
