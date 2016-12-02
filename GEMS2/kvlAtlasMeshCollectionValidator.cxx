#include "kvlAtlasMeshCollectionValidator.h"


namespace kvl
{


//
//
//
AtlasMeshCollectionValidator
::AtlasMeshCollectionValidator()
{
}



//
//
//
AtlasMeshCollectionValidator
::~AtlasMeshCollectionValidator()
{
}




//
//
//
void 
AtlasMeshCollectionValidator
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
}




//
//
//
bool
AtlasMeshCollectionValidator
::Validate( const AtlasMeshCollection* meshCollection )
{
  //std::cout << "VALIDATOR: starting to inspect the meshCollection " << meshCollection << std::endl;


  // Check correctness of meshCollection. Loop over all its cells, and check that all points to which they refer
  // actually exist
  for ( AtlasMeshCollection::CellsContainerType::ConstIterator cellIt = meshCollection->GetCells()->Begin();
        cellIt != meshCollection->GetCells()->End(); ++cellIt )
    {
    for ( AtlasMesh::CellType::PointIdIterator  pit = cellIt.Value()->PointIdsBegin();
          pit != cellIt.Value()->PointIdsEnd(); ++pit )
      {
      if ( !meshCollection->GetReferencePosition()->IndexExists( *pit ) )
        {
        std::cout << "Ouch! Cell with index " << cellIt.Index() << " refers to point with id " << *pit
                  << " which doesn't exist!" << std::endl;
        return false;
        }
      }

    }



  // Check that every tetrahedron has all its 4 triangles, 6 lines, and 4 vertices present
  bool  isValid = true;
  for ( AtlasMeshCollection::CellsContainerType::ConstIterator cellIt = meshCollection->GetCells()->Begin();
        cellIt != meshCollection->GetCells()->End(); ++cellIt )
    {
    if ( cellIt.Value()->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
      {
      continue;
      }

    AtlasMesh::CellType::PointIdConstIterator  pointIt = cellIt.Value()->PointIdsBegin();
    AtlasMesh::PointIdentifier  p0 = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier  p1 = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier  p2 = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier  p3 = *pointIt;

    int  numberOfVerticesFound = 0;
    int  numberOfLinesFound = 0;
    int  numberOfTrianglesFound = 0;

    // Loop over all cells, and check
    for ( AtlasMeshCollection::CellsContainerType::ConstIterator cellIt2 = meshCollection->GetCells()->Begin();
          cellIt2 != meshCollection->GetCells()->End(); ++cellIt2 )
      {
      // Check if all of the points are good
      bool  allThePointsAreGood = true;
      for ( AtlasMesh::CellType::PointIdConstIterator  pointIt = cellIt2.Value()->PointIdsBegin();
            pointIt != cellIt2.Value()->PointIdsEnd(); ++pointIt )
        {
        if ( ( *pointIt != p0 ) && ( *pointIt != p1 ) && ( *pointIt != p2 ) && ( *pointIt != p3 ) )
          {
          allThePointsAreGood = false;
          break;
          }

        }

      if ( allThePointsAreGood )
        {
        if ( cellIt2.Value()->GetType() == AtlasMesh::CellType::VERTEX_CELL )
          {
          numberOfVerticesFound++;
          }
        else if ( cellIt2.Value()->GetType() == AtlasMesh::CellType::LINE_CELL )
          {
          numberOfLinesFound++;
          }
        else if ( cellIt2.Value()->GetType() == AtlasMesh::CellType::TRIANGLE_CELL )
          {
          numberOfTrianglesFound++;
          }
        }

      } // End loop over all cells

      if ( numberOfVerticesFound != 4 )
        {
        std::cout << "VALIDATOR:  One or more of the vertices belong to tetrahedron with index " << cellIt.Index()
                  << " are missing" << std::endl;
        isValid = false;
        }
      if ( numberOfLinesFound != 6 )
        {
        std::cout << "VALIDATOR:  One or more of the lines belong to tetrahedron with index " << cellIt.Index()
                  << " are missing" << std::endl;
        isValid = false;
        }
      if ( numberOfTrianglesFound != 4 )
        {
        std::cout << "VALIDATOR:  One or more of the triangles belong to tetrahedron with index " << cellIt.Index()
                  << " are missing" << std::endl;
        isValid = false;
        }

    } // End loop over all tetrahedra  


  // Check that every triangle has all its 3 lines and 3 vertices present
  for ( AtlasMeshCollection::CellsContainerType::ConstIterator cellIt = meshCollection->GetCells()->Begin();
        cellIt != meshCollection->GetCells()->End(); ++cellIt )
    {
    if ( cellIt.Value()->GetType() != AtlasMesh::CellType::TRIANGLE_CELL )
      {
      continue;
      }

    AtlasMesh::CellType::PointIdConstIterator  pointIt = cellIt.Value()->PointIdsBegin();
    AtlasMesh::PointIdentifier  p0 = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier  p1 = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier  p2 = *pointIt;

    int  numberOfVerticesFound = 0;
    int  numberOfLinesFound = 0;

    // Loop over all cells, and check
    for ( AtlasMeshCollection::CellsContainerType::ConstIterator cellIt2 = meshCollection->GetCells()->Begin();
          cellIt2 != meshCollection->GetCells()->End(); ++cellIt2 )
      {
      // Check if all of the points are good
      bool  allThePointsAreGood = true;
      for ( AtlasMesh::CellType::PointIdConstIterator  pointIt = cellIt2.Value()->PointIdsBegin();
            pointIt != cellIt2.Value()->PointIdsEnd(); ++pointIt )
        {
        if ( ( *pointIt != p0 ) && ( *pointIt != p1 ) && ( *pointIt != p2 ) )
          {
          allThePointsAreGood = false;
          break;
          }

        }

      if ( allThePointsAreGood )
        {
        if ( cellIt2.Value()->GetType() == AtlasMesh::CellType::VERTEX_CELL )
          {
          numberOfVerticesFound++;
          }
        else if ( cellIt2.Value()->GetType() == AtlasMesh::CellType::LINE_CELL )
          {
          numberOfLinesFound++;
          }
        }

      } // End loop over all cells

      if ( numberOfVerticesFound != 3 )
        {
        std::cout << "VALIDATOR:  One or more of the vertices belong to triangle with index " << cellIt.Index()
                  << " are missing" << std::endl;
        isValid = false;
        }
      if ( numberOfLinesFound != 3 )
        {
        std::cout << "VALIDATOR:  One or more of the lines belong to triangle with index " << cellIt.Index()
                  << " are missing" << std::endl;
        isValid = false;
        }

    } // End loop over all triangles  

  // Check that every line has all its 2 vertices present
  for ( AtlasMeshCollection::CellsContainerType::ConstIterator cellIt = meshCollection->GetCells()->Begin();
        cellIt != meshCollection->GetCells()->End(); ++cellIt )
    {
    if ( cellIt.Value()->GetType() != AtlasMesh::CellType::LINE_CELL )
      {
      continue;
      }

    AtlasMesh::CellType::PointIdConstIterator  pointIt = cellIt.Value()->PointIdsBegin();
    AtlasMesh::PointIdentifier  p0 = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier  p1 = *pointIt;

    int  numberOfVerticesFound = 0;

    // Loop over all cells, and check
    for ( AtlasMeshCollection::CellsContainerType::ConstIterator cellIt2 = meshCollection->GetCells()->Begin();
          cellIt2 != meshCollection->GetCells()->End(); ++cellIt2 )
      {
      // Check if all of the points are good
      bool  allThePointsAreGood = true;
      for ( AtlasMesh::CellType::PointIdConstIterator  pointIt = cellIt2.Value()->PointIdsBegin();
            pointIt != cellIt2.Value()->PointIdsEnd(); ++pointIt )
        {
        if ( ( *pointIt != p0 ) && ( *pointIt != p1 ) )
          {
          allThePointsAreGood = false;
          break;
          }

        }

      if ( allThePointsAreGood )
        {
        if ( cellIt2.Value()->GetType() == AtlasMesh::CellType::VERTEX_CELL )
          {
          numberOfVerticesFound++;
          }
        }

      } // End loop over all cells

      if ( numberOfVerticesFound != 2 )
        {
        std::cout << "VALIDATOR:  One or more of the vertices belong to line with index " << cellIt.Index()
                  << " are missing" << std::endl;
        isValid = false;
        }

    } // End loop over all lines  


  // Finally, check that there no points that belong to no cell at all.
  for ( AtlasMeshCollection::PointsContainerType::ConstIterator  pointIt = meshCollection->GetReferencePosition()->Begin();
        pointIt != meshCollection->GetReferencePosition()->End();
        ++pointIt )
    {
    bool pointBelongsToSomeCell = false;
    for ( AtlasMeshCollection::CellsContainerType::ConstIterator cellIt = meshCollection->GetCells()->Begin();
          cellIt != meshCollection->GetCells()->End(); ++cellIt )
      {
      for ( AtlasMesh::CellType::PointIdConstIterator  pointIt2 = cellIt.Value()->PointIdsBegin();
            pointIt2 != cellIt.Value()->PointIdsEnd(); ++pointIt2 )
        {
        if ( *pointIt2 == pointIt.Index() )
          {
          pointBelongsToSomeCell = true;
          break;
          }
        }

      }

    if ( !pointBelongsToSomeCell )
      {
      std::cout << "VALIDATOR: point with id " << pointIt.Index() << " doesn't belong to any cell!" << std::endl;
      isValid = false;
      }

    }


  // Loop over all tetrahedra, and check if their volume is positive
  for ( AtlasMeshCollection::CellsContainerType::ConstIterator cellIt = meshCollection->GetCells()->Begin();
        cellIt != meshCollection->GetCells()->End(); ++cellIt )
    {
    if ( cellIt.Value()->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
      {
      continue;
      }


    // Retrieve the reference positions of the vertices
    AtlasMesh::CellType::PointIdIterator  pit = cellIt.Value()->PointIdsBegin();
    const float x0 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 0 ];
    const float y0 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 1 ];
    const float z0 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 2 ];
    ++pit;
    const float x1 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 0 ];
    const float y1 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 1 ];
    const float z1 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 2 ];
    ++pit;
    const float x2 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 0 ];
    const float y2 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 1 ];
    const float z2 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 2 ];
    ++pit;
    const float x3 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 0 ];
    const float y3 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 1 ];
    const float z3 = meshCollection->GetReferencePosition()->ElementAt( *pit )[ 2 ];


    // Calculate the volume of the tetrahedron in reference position. Lambda is the Jacobian of
    // the transform from a standarizedized tetrahedron, which has volume 1/6, to the tethrahedron
    // in reference position
    const float  lambda11 = -x0 + x1;
    const float  lambda21 = -y0 + y1;
    const float  lambda31 = -z0 + z1;
    const float  lambda12 = -x0 + x2;
    const float  lambda22 = -y0 + y2;
    const float  lambda32 = -z0 + z2;
    const float  lambda13 = -x0 + x3;
    const float  lambda23 = -y0 + y3;
    const float  lambda33 = -z0 + z3;
    const float  volume = ( lambda11 * ( lambda22*lambda33 - lambda32*lambda23 )
                            - lambda12 * ( lambda21*lambda33 - lambda31*lambda23 )
                            + lambda13 * ( lambda21*lambda32 - lambda31*lambda22 ) ) / 6;

    if ( volume <= 0  )
      {
      std::cout << "VALIDATOR: tetrahedron with id " << cellIt.Index() << " has negative reference volume!" << std::endl;
      isValid = false;
      }
    }  


  return isValid;
}



} // end namespace kvl
