#include "kvlAtlasMeshCollection.h"

#include <gzstream.h>
#include <fstream>

#include "vnl/vnl_inverse.h"
#include "vnl/vnl_matrix_fixed.h"
#include "kvlTetrahedronAspectRatio.h"



namespace kvl
{


//
//
//
AtlasMeshCollection
::AtlasMeshCollection()
{
  m_PointParameters = 0;
  m_Cells = 0;
  m_ReferenceTetrahedronInfos = 0;
  m_ReferencePosition = 0;
  m_K = 10;

  m_CellLinks = 0;
  
}



//
//
//
AtlasMeshCollection
::~AtlasMeshCollection()
{

  // Clean up cells if no-one is using them anymore. First clear the cached meshes,
  // so that we don't accidentally think someone else is still using the cells.
  m_Meshes.clear();
  if ( m_Cells )
    {
    if ( m_Cells->GetReferenceCount() == 1 ) 
      {
      
      AtlasMesh::CellsContainer::Iterator cellIt  = m_Cells->Begin();
      while( cellIt != m_Cells->End() )
        {
        delete cellIt.Value();
        ++cellIt; 
        }
      
      }
         
    }

}




//
//
//
void 
AtlasMeshCollection
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{



}



//
//
//
void
AtlasMeshCollection
::GenerateFromSingleMesh( AtlasMesh* mesh, unsigned int numberOfMeshes, double K )
{
  // Clean up cache of what we may have already
  m_ReferenceTetrahedronInfos = 0;
  m_Meshes.clear(); // delete all content from the mesh container
  m_CellLinks = 0;

  // Initialize topology
  m_Cells = mesh->GetCells();
    
  // Initialize parameters
  m_PointParameters = mesh->GetPointData();
  m_ReferencePosition = mesh->GetPoints();
  m_K = K;
 
  // Initialize positions by making copies of the mesh position
  PointsContainerType::ConstPointer  sourcePosition = mesh->GetPoints();
  m_Positions.clear();// delete all content from the positions container
  for ( unsigned int i = 0; i < numberOfMeshes; i++ )
    {
    // Create a new Points container for this mesh
    PointsContainerType::Pointer  target = PointsContainerType::New();
    // Go through all source points
    PointsContainerType::ConstIterator  sourceIt = sourcePosition->Begin();
    while ( sourceIt != sourcePosition->End() )
      {
      // insert source coords for this point into target
      target->InsertElement( sourceIt.Index(), sourceIt.Value() );
      ++sourceIt;
      }

    // Add the new mesh points container to the mesh collection
    m_Positions.push_back( target );
    }                                                
  
}
  
/*!
  \fn const AtlasMesh* AtlasMeshCollection::GetMesh()
  \brief Returns a pointer to the given mesh number. If the m_Meshes
  private variable has not been initialized, then it gets initialized.
  m_Positions, m_Cells, and m_PointParams must be initialized.
 */
const AtlasMesh*
AtlasMeshCollection
::GetMesh( unsigned int meshNumber ) const
{
  // Sanity check on requested mesh
  if ( meshNumber >= m_Positions.size() )
    {
    return 0;
    }

  // If cached meshes container is empty, create it
  if ( m_Meshes.size() == 0 )
    {
    // go thru each mesh in the collection
    for ( unsigned int i = 0; i < m_Positions.size(); i++ )
      {
      // Construct an atlas mesh
      AtlasMesh::Pointer  mesh = AtlasMesh::New();
      mesh->SetPoints( m_Positions[ i ] );
      mesh->SetCells( m_Cells );
      mesh->SetPointData( m_PointParameters );
      mesh->SetCellData( const_cast< CellDataContainerType* >( this->GetReferenceTetrahedronInfos() ) );

      // Add to the container
      m_Meshes.push_back( mesh );
      }
    }
    
  // Return the cached mesh
  return m_Meshes[ meshNumber ].GetPointer();

}

//
//
// 
AtlasMesh::ConstPointer
AtlasMeshCollection
::GetReferenceMesh() const
{
  // Construct an atlas mesh
  AtlasMesh::Pointer  mesh = AtlasMesh::New();
  mesh->SetPoints( m_ReferencePosition );
  mesh->SetCells( m_Cells );
  mesh->SetPointData( m_PointParameters );
  mesh->SetCellData( const_cast< CellDataContainerType* >( this->GetReferenceTetrahedronInfos() ) );

  return mesh.GetPointer();
}

//
//
//
/*!
  \fn const AtlasMeshCollection::CellDataContainerType* AtlasMeshCollection::GetReferenceTetrahedronInfos()
  \brief Returns a pointer to the m_ReferenceTetrahedronInfos private
  member. If the priv member is NULL, then computes the volume and matrix inverse.
 */
const AtlasMeshCollection::CellDataContainerType*
AtlasMeshCollection
::GetReferenceTetrahedronInfos() const
{
  if ( !m_ReferenceTetrahedronInfos )
    {
    // Calculate the tetrahedron volume parameters
    //std::cout << "Calculating the tetrahedron volume parameters...";
    m_ReferenceTetrahedronInfos = CellDataContainerType::New();
    for ( AtlasMesh::CellsContainer::Iterator cellIt = m_Cells->Begin();
          cellIt != m_Cells->End();
          ++cellIt )
      {
      AtlasMesh::CellType*  cell = cellIt.Value();
    
      if( cell->GetType() == AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
        // Create info structure and fill in all the fields one by one
        ReferenceTetrahedronInfo  info;

        // Retrieve the reference positions of the vertices
        AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
        const double x0 = m_ReferencePosition->ElementAt( *pit )[ 0 ];
        const double y0 = m_ReferencePosition->ElementAt( *pit )[ 1 ];
        const double z0 = m_ReferencePosition->ElementAt( *pit )[ 2 ];
        ++pit;
        const double x1 = m_ReferencePosition->ElementAt( *pit )[ 0 ];
        const double y1 = m_ReferencePosition->ElementAt( *pit )[ 1 ];
        const double z1 = m_ReferencePosition->ElementAt( *pit )[ 2 ];
        ++pit;
        const double x2 = m_ReferencePosition->ElementAt( *pit )[ 0 ];
        const double y2 = m_ReferencePosition->ElementAt( *pit )[ 1 ];
        const double z2 = m_ReferencePosition->ElementAt( *pit )[ 2 ];
        ++pit;
        const double x3 = m_ReferencePosition->ElementAt( *pit )[ 0 ];
        const double y3 = m_ReferencePosition->ElementAt( *pit )[ 1 ];
        const double z3 = m_ReferencePosition->ElementAt( *pit )[ 2 ];


        // Calculate the volume of the tetrahedron in reference position. Lambda is the Jacobian of
        // the transform from a standarizedized tetrahedron, which has volume 1/6, to the tethrahedron
        // in reference position
        const double  lambda11 = -x0 + x1;
        const double  lambda21 = -y0 + y1;
        const double  lambda31 = -z0 + z1;
        const double  lambda12 = -x0 + x2;
        const double  lambda22 = -y0 + y2;
        const double  lambda32 = -z0 + z2;
        const double  lambda13 = -x0 + x3;
        const double  lambda23 = -y0 + y3;
        const double  lambda33 = -z0 + z3;
        const double referenceVolume = ( lambda11 * ( lambda22*lambda33 - lambda32*lambda23 )
                                      - lambda12 * ( lambda21*lambda33 - lambda31*lambda23 )
                                      + lambda13 * ( lambda21*lambda32 - lambda31*lambda22 ) ) / 6;

        info.m_ReferenceVolumeTimesK = m_K * referenceVolume;

        // Let's precalculate inv( [ p0 p1 p2 p3; 1 1 1 1 ] ) and call it Z
        vnl_matrix_fixed< double, 4, 4 >  referenceMatrix;
        referenceMatrix.put( 0, 0, x0 );
        referenceMatrix.put( 1, 0, y0 );
        referenceMatrix.put( 2, 0, z0 );
        referenceMatrix.put( 3, 0, 1.0f );

        referenceMatrix.put( 0, 1, x1 );
        referenceMatrix.put( 1, 1, y1 );
        referenceMatrix.put( 2, 1, z1 );
        referenceMatrix.put( 3, 1, 1.0f );

        referenceMatrix.put( 0, 2, x2 );
        referenceMatrix.put( 1, 2, y2 );
        referenceMatrix.put( 2, 2, z2 );
        referenceMatrix.put( 3, 2, 1.0f );

        referenceMatrix.put( 0, 3, x3 );
        referenceMatrix.put( 1, 3, y3 );
        referenceMatrix.put( 2, 3, z3 );
        referenceMatrix.put( 3, 3, 1.0f );

	// DNG: might be faster to have a dedicated 4x4 inverse, if not there already
        vnl_matrix_fixed< double, 4, 4 >  inverseReferenceMatrix = vnl_inverse( referenceMatrix );
        info.m_Z11 = inverseReferenceMatrix.get( 0, 0 );
        info.m_Z21 = inverseReferenceMatrix.get( 1, 0 );
        info.m_Z31 = inverseReferenceMatrix.get( 2, 0 );
        info.m_Z41 = inverseReferenceMatrix.get( 3, 0 );

        info.m_Z12 = inverseReferenceMatrix.get( 0, 1 );
        info.m_Z22 = inverseReferenceMatrix.get( 1, 1 );
        info.m_Z32 = inverseReferenceMatrix.get( 2, 1 );
        info.m_Z42 = inverseReferenceMatrix.get( 3, 1 );

        info.m_Z13 = inverseReferenceMatrix.get( 0, 2 );
        info.m_Z23 = inverseReferenceMatrix.get( 1, 2 );
        info.m_Z33 = inverseReferenceMatrix.get( 2, 2 );
        info.m_Z43 = inverseReferenceMatrix.get( 3, 2 );
  
	if(0){//if(cellIt.Index() == 1652908){
         std::cout << "Calculated reference tetrahedron info as follows:" << std::endl;
         std::cout << "     p0: [" << x0 << ", " << y0 << ", " << z0 << "]" << std::endl;
         std::cout << "     p1: [" << x1 << ", " << y1 << ", " << z1 << "]" << std::endl;
         std::cout << "     p2: [" << x2 << ", " << y2 << ", " << z2 << "]" << std::endl;
         std::cout << "     p3: [" << x3 << ", " << y3 << ", " << z3 << "]" << std::endl;
         std::cout << "      first column of Z: [" << info.m_Z11 << ", "
                                                   << info.m_Z21 << ", "
                                                   << info.m_Z31 << ", "
                                                   << info.m_Z41 << "]'" << std::endl;
         std::cout << "     second column of Z: [" << info.m_Z12 << ", "
                                                   << info.m_Z22 << ", "
                                                   << info.m_Z32 << ", "
                                                   << info.m_Z42 << "]'" << std::endl;
         std::cout << "      third column of Z: [" << info.m_Z13 << ", "
                                                   << info.m_Z23 << ", "
                                                   << info.m_Z33 << ", "
                                                   << info.m_Z43 << "]'" << std::endl;
        std::cout << "      m_ReferenceVolumeTimesK: " << info.m_ReferenceVolumeTimesK << std::endl;
	}
        m_ReferenceTetrahedronInfos->InsertElement( cellIt.Index(), info );
        }
      
      } // End loop over all cells

    //std::cout << "...done!" << std::endl;
    }

  return m_ReferenceTetrahedronInfos;

}




//
//
// 
bool
AtlasMeshCollection
::Write( const char* fileName ) const
{

  // Only write if all fields are set
  if ( ( !m_PointParameters ) || ( !m_Cells ) || ( !m_ReferencePosition ) ||
       ( !m_Positions.size() ) )
    {
    std::cerr << "Not a complete mesh collection" << std::endl;
    return false;
    }
    
  // Open the file name
  std::string  zippedFileName = std::string( fileName ) + ".gz";
  ogzstream  out( zippedFileName.c_str() );
  if ( out.bad() )
    {
    std::cerr << "Can't open " << zippedFileName << " for writing." << std::endl;
    return false;
    }
    
  // Write some general mesh info
  out << "Number of points: " << m_PointParameters->Size() << std::endl;
  out << "Number of cells: " << m_Cells->Size() << std::endl;
  out << "Number of labels: " << m_PointParameters->Begin().Value().m_Alphas.size() << std::endl;
  out << "Number of meshes: " << this->GetNumberOfMeshes() << std::endl;

  // Write out reference position
  out << "Reference position " << std::endl;
  for ( PointsContainerType::ConstIterator  pointIt = m_ReferencePosition->Begin();
        pointIt != m_ReferencePosition->End();
        ++pointIt )
    {
    out << "   " << pointIt.Index() << "  "
        << pointIt.Value()[ 0 ] << "  " << pointIt.Value()[ 1 ] << "  " << pointIt.Value()[ 2 ] << std::endl;
    }
  
  // Write out K
  out << "K: " << m_K << std::endl;
  
    
  // Write out positions
  for ( unsigned int positionNumber = 0; positionNumber < m_Positions.size(); positionNumber++ )
    {
    out << "Position " << positionNumber << std::endl;
    
    PointsContainerType::Pointer  position = m_Positions[ positionNumber ];
    
    PointsContainerType::ConstIterator  pointIt = position->Begin();
    while ( pointIt != position->End() )
      {
      out << "   " << pointIt.Index() << "  " 
          << pointIt.Value()[ 0 ] << "  " << pointIt.Value()[ 1 ] << "  " << pointIt.Value()[ 2 ] << std::endl;
      
      ++pointIt;
      }
    
    }
  
  // Write out cells
  out << "Cells: " << std::endl;
  CellsContainerType::ConstIterator  cellIt = m_Cells->Begin();
  CellsContainerType::ConstIterator  cellEnd = m_Cells->End();
  while ( cellIt != cellEnd )
    {
    AtlasMesh::CellType*  cell = cellIt.Value();
    
    out << "   " << cellIt.Index() << "   ";
    
    if ( cell->GetType() == AtlasMesh::CellType::VERTEX_CELL )
      {
      out << "VERTEX   ";
      }
    else if ( cell->GetType() == AtlasMesh::CellType::LINE_CELL )
      {
      out << "LINE   ";
      }
    else if ( cell->GetType() == AtlasMesh::CellType::TRIANGLE_CELL )
      {
      out << "TRIANGLE   ";
      }
    else if ( cell->GetType() == AtlasMesh::CellType::TETRAHEDRON_CELL )
      {
      out << "TETRAHEDRON   ";
      }
    else
      {
      itkExceptionMacro( "Mesh collection may only contain vertices, lines, triangles, and tetrahedra." );
      }
        
    AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
    while ( pit != cell->PointIdsEnd() )
      {
      out << *pit << " ";
      ++pit;
      }
    out << std::endl;
    
    ++cellIt;
    }
  
    
  // Write out point parameters
  out << "Point parameters: " << std::endl;
  PointDataContainerType::ConstIterator  pointParamIt = m_PointParameters->Begin();
  while ( pointParamIt !=  m_PointParameters->End() )
    {
    out << "   " << pointParamIt.Index() << "   ";
    
    // Alphas
    for ( unsigned int i=0; i < pointParamIt.Value().m_Alphas.size(); i++ )
      {
      out << pointParamIt.Value().m_Alphas[ i ] << "  ";
      }
    
    
    // Can-change alphas
    if ( pointParamIt.Value().m_CanChangeAlphas )
      {
      out << "true  ";
      }
    else
      {
      out << "false  ";
      }
      
    // Can-move X
    if ( pointParamIt.Value().m_CanMoveX )
      {
      out << "true  ";
      }
    else
      {
      out << "false  ";
      }
      
    // Can-move Y
    if ( pointParamIt.Value().m_CanMoveY )
      {
      out << "true  ";
      }
    else
      {
      out << "false  ";
      }
    
    // Can-move Z
    if ( pointParamIt.Value().m_CanMoveZ )
      {
      out << "true  ";
      }
    else
      {
      out << "false  ";
      }
    
    out << std::endl;
    ++pointParamIt;
    }
        
  return true;
    
}


//
//
//
template< class T >
static bool GetValue( const std::string line, const std::string& searchString, T& value )
{
  std::size_t  foundPosition;
  std::string  valueString = line;
  
  if ( ( foundPosition = valueString.find( searchString ) ) == std::string::npos )
    return false;
    
  valueString.erase( foundPosition, searchString.size() ); 
  std::istringstream  valueStream( valueString );
  valueStream >> value;

  return true;
}


//  
//
//
bool
AtlasMeshCollection
::Read( const char* fileName )
{
  // Clean up cache of what we may have already
  m_ReferenceTetrahedronInfos = 0;
  m_Meshes.clear();
  m_CellLinks = 0;

#if 0
  std::string  zippedFileName = std::string( fileName ) + ".gz";
  igzstream  in( zippedFileName.c_str() );
  if ( in.bad() )
    {
    std::cerr << "Can't open " << zippedFileName << " for reading" << std::endl;
    return false;
    }
#else
  igzstream  in( fileName );
  if ( !in.rdbuf()->is_open() )
    {
    // Try the gz extenstion
    const std::string  zippedFileName = std::string( fileName ) + ".gz";
    in.open( zippedFileName.c_str() );
    if ( !in.rdbuf()->is_open() )
      {
      std::cerr << "Can't open " << fileName << " for reading" << std::endl;
      return false;
      }
    }
#endif
  
  const int size = 1023;
  char buffer[ size ];
  std::string  line;
  
  // Number of points
  if ( !in.getline( buffer, size ) )
    return false;
  line = buffer;
  unsigned int  numberOfPoints;
  if ( !GetValue( line, "Number of points:",  numberOfPoints ) )
    return false;

  // Number of cells
  if ( !in.getline( buffer, size ) )
    return false;
  line = buffer;
  unsigned int  numberOfCells;
  if ( !GetValue( line, "Number of cells:",  numberOfCells ) )
    return false;
  
  // Number of labels
  if ( !in.getline( buffer, size ) )
    return false;
  line = buffer;
  unsigned int  numberOfLabels;
  if ( !GetValue( line, "Number of labels:",  numberOfLabels ) )
    return false;

  // Number of meshes
  if ( !in.getline( buffer, size ) )
    return false;
  line = buffer;
  unsigned int  numberOfMeshes;
  if ( !GetValue( line, "Number of meshes:",  numberOfMeshes ) )
    return false;
  
  // Some output
  //std::cout << "numberOfPoints: " << numberOfPoints << std::endl;
  //std::cout << "numberOfCells: " << numberOfCells << std::endl;
  //std::cout << "numberOfLabels: " << numberOfLabels << std::endl;
  //std::cout << "numberOfMeshes: " << numberOfMeshes << std::endl;


  // Read reference position
  //std::cout << "Reading reference position " << std::endl;
  
  // Skip the line saying "Reference position"
  if ( !in.getline( buffer, size ) )
    return false;
  
#ifndef USE_DYNAMIC_MESH
  std::map< AtlasMesh::PointIdentifier, AtlasMesh::PointIdentifier >  pointIdCompressionLookupTable;
#endif  
  
  m_ReferencePosition = PointsContainerType::New();
  for ( unsigned int pointNumber = 0; pointNumber < numberOfPoints; pointNumber++ )
    {
    // Get the line
    if ( !in.getline( buffer, size ) )
      return false;
    line = buffer;
    std::istringstream  lineStream( line );
    AtlasMesh::PointType  point;
    AtlasMesh::PointIdentifier  id;
    lineStream >> id >> point[ 0 ] >> point[ 1 ] >> point[ 2 ];

#ifndef USE_DYNAMIC_MESH
    const AtlasMesh::PointIdentifier  originalId = id;
    id = pointIdCompressionLookupTable.size();
    pointIdCompressionLookupTable[ originalId ] = id;
#endif  

    // On Mac, istringstream cannot convert strings holding tiny float values into floats
    if ( lineStream.fail() )
      {
      itkExceptionMacro( "The following line could not be parsed into floats: " << line );
      }
    
    //std::cout << "   Adding point " << id << ": " << point << std::endl;
    m_ReferencePosition->InsertElement( id, point );
    }
  

  // Read K
  if ( !in.getline( buffer, size ) )
    return false;
  line = buffer;
  if ( !GetValue( line, "K:",  m_K ) )
    return false;
   
  
  // Read positions
  m_Positions.clear();
  for ( unsigned int i = 0; i < numberOfMeshes; i++ )
    {
    //std::cout << "Reading positions " << i << std::endl;
    
    // Skip the line saying "Position X"
    if ( !in.getline( buffer, size ) )
      return false;
    
    PointsContainerType::Pointer  position = PointsContainerType::New();
    for ( unsigned int pointNumber = 0; pointNumber < numberOfPoints; pointNumber++ )
      {
      // Get the line
      if ( !in.getline( buffer, size ) )
        return false;
      line = buffer;
      std::istringstream  lineStream( line );
      AtlasMesh::PointType  point;
      AtlasMesh::PointIdentifier  id;
      lineStream >> id >> point[ 0 ] >> point[ 1 ] >> point[ 2 ];
      
#ifndef USE_DYNAMIC_MESH
      id = pointIdCompressionLookupTable.find( id )->second;
#endif    
    
      // On Mac, istringstream cannot convert strings holding tiny float values into floats
      if ( lineStream.fail() )
        {
        itkExceptionMacro( "The following line could not be parsed into floats: " << line );
        }

      //std::cout << "   Adding point " << id << ": " << point << std::endl;
      position->InsertElement( id, point );
      }
      
      
    // Push back
    m_Positions.push_back( position );
    }                                                
   
  
  // Read the cells
  
  // Skip the line saying "Cells:"
  if ( !in.getline( buffer, size ) )
    return false;

  typedef itk::VertexCell< AtlasMesh::CellType >    VertexCell;
  typedef itk::LineCell< AtlasMesh::CellType >      LineCell;
  typedef itk::TriangleCell< AtlasMesh::CellType >  TriangleCell;
  typedef itk::TetrahedronCell< AtlasMesh::CellType >  TetrahedronCell;
  
#ifndef USE_DYNAMIC_MESH
  std::map< AtlasMesh::CellIdentifier, AtlasMesh::CellIdentifier >  cellIdCompressionLookupTable; 
#endif  
  
  m_Cells = CellsContainerType::New();
  for ( unsigned int cellNumber = 0; cellNumber < numberOfCells; cellNumber++ )
    {
    // Get the line
    if ( !in.getline( buffer, size ) )
      return false;
    line = buffer;
    std::istringstream  lineStream( line );
    AtlasMesh::CellIdentifier  cellId;
    std::string  cellType;
    lineStream >> cellId >> cellType;
    
#ifndef USE_DYNAMIC_MESH
    const AtlasMesh::CellIdentifier  originalCellId = cellId;
    cellId = cellIdCompressionLookupTable.size();
    cellIdCompressionLookupTable[ originalCellId ] = cellId;
#endif  
  
    //std::cout << "Detected cell with id: " << cellId << " of type: " << cellType << std::endl;     
    
    if ( cellType == "VERTEX" )
      {
      // Read the id of the point
      AtlasMesh::PointIdentifier  pointId0;
      lineStream >> pointId0;
    
#ifndef USE_DYNAMIC_MESH
      pointId0 = pointIdCompressionLookupTable.find( pointId0 )->second;
#endif    
      
      // Create a new vertex cell
      AtlasMesh::CellAutoPointer newCell;
      newCell.TakeOwnership( new VertexCell );
      newCell->SetPointId( 0, pointId0 );
      
      // Add the cell
      m_Cells->InsertElement( cellId, newCell.ReleaseOwnership() );
      }
    else if ( cellType == "LINE" )
      {
      // Read the id of the points
      AtlasMesh::PointIdentifier  pointId0;
      AtlasMesh::PointIdentifier  pointId1;
      lineStream >> pointId0 >> pointId1;

#ifndef USE_DYNAMIC_MESH
      pointId0 = pointIdCompressionLookupTable.find( pointId0 )->second;
      pointId1 = pointIdCompressionLookupTable.find( pointId1 )->second;
#endif    
    
      // Create a new line cell
      AtlasMesh::CellAutoPointer newCell;
      newCell.TakeOwnership( new LineCell );
      newCell->SetPointId( 0, pointId0 );
      newCell->SetPointId( 1, pointId1 );
      
      // Add the cell
      m_Cells->InsertElement( cellId, newCell.ReleaseOwnership() );
      
      }
    else if ( cellType == "TRIANGLE" )
      {
      // Read the id of the points
      AtlasMesh::PointIdentifier  pointId0;
      AtlasMesh::PointIdentifier  pointId1;
      AtlasMesh::PointIdentifier  pointId2;
      lineStream >> pointId0 >> pointId1 >> pointId2;
      
#ifndef USE_DYNAMIC_MESH
      pointId0 = pointIdCompressionLookupTable.find( pointId0 )->second;
      pointId1 = pointIdCompressionLookupTable.find( pointId1 )->second;
      pointId2 = pointIdCompressionLookupTable.find( pointId2 )->second;
#endif    
    
      // Create a new triangle cell
      AtlasMesh::CellAutoPointer newCell;
      newCell.TakeOwnership( new TriangleCell );
      newCell->SetPointId( 0, pointId0 );
      newCell->SetPointId( 1, pointId1 );
      newCell->SetPointId( 2, pointId2 );
      
      // Add the cell
      m_Cells->InsertElement( cellId, newCell.ReleaseOwnership() );

      }
    else
      {
      // Read the id of the points
      AtlasMesh::PointIdentifier  pointId0;
      AtlasMesh::PointIdentifier  pointId1;
      AtlasMesh::PointIdentifier  pointId2;
      AtlasMesh::PointIdentifier  pointId3;
      lineStream >> pointId0 >> pointId1 >> pointId2 >> pointId3;
      
#ifndef USE_DYNAMIC_MESH
      pointId0 = pointIdCompressionLookupTable.find( pointId0 )->second;
      pointId1 = pointIdCompressionLookupTable.find( pointId1 )->second;
      pointId2 = pointIdCompressionLookupTable.find( pointId2 )->second;
      pointId3 = pointIdCompressionLookupTable.find( pointId3 )->second;
#endif    
    
      // Create a new triangle cell
      AtlasMesh::CellAutoPointer newCell;
      newCell.TakeOwnership( new TetrahedronCell );
      newCell->SetPointId( 0, pointId0 );
      newCell->SetPointId( 1, pointId1 );
      newCell->SetPointId( 2, pointId2 );
      newCell->SetPointId( 3, pointId3 );
      
      // Add the cell
      m_Cells->InsertElement( cellId, newCell.ReleaseOwnership() );

      }
    
    }
    
    
      
  // Read point parameters
  
  // Skip the line saying "Point parameters:"
  if ( !in.getline( buffer, size ) )
    return false;
  
  m_PointParameters = PointDataContainerType::New();
  for ( unsigned int pointNumber = 0; pointNumber < numberOfPoints; pointNumber++ )
    {
     // Get the line
    if ( !in.getline( buffer, size ) )
      return false;
    line = buffer;
    std::istringstream  lineStream( line );
    AtlasMesh::PointIdentifier  pointId;
    lineStream >> pointId;
#ifndef USE_DYNAMIC_MESH
    pointId = pointIdCompressionLookupTable.find( pointId )->second;
#endif  
    AtlasAlphasType   alphas( numberOfLabels );
    AtlasMesh::PixelType  pointParameter;
    pointParameter.m_Alphas = alphas;
    for ( unsigned int labelNumber=0; labelNumber < numberOfLabels; labelNumber++ )
      {
      lineStream >> pointParameter.m_Alphas[ labelNumber ];

        // On Mac, istringstream cannot convert strings holding tiny float values into float
        if ( lineStream.fail() )
        {
          pointParameter.m_Alphas[ labelNumber ] = 0.0f;
          lineStream.clear();
        }

      }
    std::string  canChangeAlphasString;
    std::string  canMoveXString;
    std::string  canMoveYString;
    std::string  canMoveZString;
    lineStream >> canChangeAlphasString >> canMoveXString >> canMoveYString >> canMoveZString;
    if ( canChangeAlphasString == "true" )
      {
      pointParameter.m_CanChangeAlphas = true;
      }
    else
      {
      pointParameter.m_CanChangeAlphas = false;
      }
    if ( canMoveXString == "true" )
      {
      pointParameter.m_CanMoveX = true;
      }
    else
      {
      pointParameter.m_CanMoveX = false;
      }
    if ( canMoveYString == "true" )
      {
      pointParameter.m_CanMoveY = true;
      }
    else
      {
      pointParameter.m_CanMoveY = false;
      }
    if ( canMoveZString == "true" )
      {
      pointParameter.m_CanMoveZ = true;
      }
    else
      {
      pointParameter.m_CanMoveZ = false;
      }
   
    // Insert 
    m_PointParameters->InsertElement( pointId, pointParameter );
    }
        

  return true;
}





/*!
  \fn AtlasMeshCollection::Construct()
  \brief Construct a mesh collection out of whole cloth.  The mesh is
  constructed by defining cubes in a volume. Each cube is filled with
  five tetrahedra. meshSize is the number of cubes in each dim. The
  domainSize is the width of the volume in each dimension (units
  unclear, probably best to think of them as 1mm, or
  not). forceBorderVerticesToBackground will cause the vertices at the
  edge of the volume to the first class (which is usually the
  background).
*/
void
AtlasMeshCollection
::Construct( const unsigned int*  meshSize, const unsigned int*  domainSize,
             double K, 
             unsigned int numberOfClasses, unsigned int numberOfMeshes,
             bool  forceBorderVerticesToBackground )
{
  // Clean up cache of what we may have already
  m_ReferenceTetrahedronInfos = 0;
  m_Meshes.clear();
  m_CellLinks = 0;

  // Create a single mesh, then replicate it below
  // Use a mesh source to create the mesh
  MeshSourceType::Pointer  meshSource = MeshSourceType::New();
  for ( unsigned int  x = 0; x < meshSize[ 0 ]-1; x++ )
    {
    for ( unsigned int  y = 0; y < meshSize[ 1 ]-1; y++ )
      {
      for ( unsigned int  z = 0; z < meshSize[ 2 ]-1; z++ )
        {
	// Get the two extreme corners of the cube
        float  x1 = static_cast< float >( x ) * static_cast< float >( domainSize[ 0 ] - 1 )
                                              / static_cast< float >( meshSize[ 0 ] - 1 );
        float  y1 = static_cast< float >( y ) * static_cast< float >( domainSize[ 1 ] - 1 )
                                              / static_cast< float >( meshSize[ 1 ] - 1 );
        float  z1 = static_cast< float >( z ) * static_cast< float >( domainSize[ 2 ] - 1 )
                                              / static_cast< float >( meshSize[ 2 ] - 1 );

        float  x2 = static_cast< float >( x + 1 ) * static_cast< float >( domainSize[ 0 ] - 1 )
                                              / static_cast< float >( meshSize[ 0 ] - 1 );
        float  y2 = static_cast< float >( y + 1 ) * static_cast< float >( domainSize[ 1 ] - 1 )
                                              / static_cast< float >( meshSize[ 1 ] - 1 );
        float  z2 = static_cast< float >( z + 1 ) * static_cast< float >( domainSize[ 2 ] - 1 )
                                              / static_cast< float >( meshSize[ 2 ] - 1 );

#if 0
        x1 = static_cast< int >( x1 + 0.5 );
        y1 = static_cast< int >( y1 + 0.5 );
        z1 = static_cast< int >( z1 + 0.5 );

        x2 = static_cast< int >( x2 + 0.5 );
        y2 = static_cast< int >( y2 + 0.5 );
        z2 = static_cast< int >( z2 + 0.5 );
#endif

  
        const double  p0[] = { x1, y1, z1 };
        const double  p1[] = { x2, y1, z1 };
        const double  p2[] = { x1, y2, z1 };
        const double  p3[] = { x2, y2, z1 };
        const double  p4[] = { x1, y1, z2 };
        const double  p5[] = { x2, y1, z2 };
        const double  p6[] = { x1, y2, z2 };
        const double  p7[] = { x2, y2, z2 };
 
        // Each cube will be filled by 5 tetrahedra. There are two different configurations however that should alternate in a checker-board pattern.
        bool flippedConfiguration = true;
        if ( x % 2 )
          {
          flippedConfiguration = !flippedConfiguration;
          }
        if ( y % 2 )
          {
          flippedConfiguration = !flippedConfiguration;
          }
        if ( z % 2 )
          {
          flippedConfiguration = !flippedConfiguration;
          }

        this->FillCubeWithTetrahedra( meshSource, flippedConfiguration,
                                      p0, p1, p2, p3, p4, p5, p6, p7 );
        }
      }
    }

  // Alpha[Class] is the probability that the given point belongs to class Class.
  // Each point needs an alpha vector, but set up some defaults first
  // Assign flat alphas as a starting point. 
  AtlasAlphasType   flatAlphasEntry( numberOfClasses );
  flatAlphasEntry.Fill( 1.0f / static_cast< float >( numberOfClasses ) );
  
  // Vertices lying on the border can not move freely and belong to first (background) class
  AtlasAlphasType   borderAlphasEntry( numberOfClasses );
  if ( forceBorderVerticesToBackground )
    {
    borderAlphasEntry.Fill( 0.0f );
    borderAlphasEntry[ 0 ] = 1.0f;
    }
  else
    {
    borderAlphasEntry = flatAlphasEntry;
    }

  // Now go through all the points
  for ( AtlasMesh::PointsContainer::ConstIterator  pointIt = meshSource->GetOutput()->GetPoints()->Begin();
        pointIt != meshSource->GetOutput()->GetPoints()->End();
        ++pointIt )
    {

    // Goal here is to fill the pointParameters struct for this point
    AtlasMesh::PixelType  pointParameters;
    pointParameters.m_Alphas = flatAlphasEntry;
    pointParameters.m_CanChangeAlphas = true;
    
    // pointIt.Value() will return an array[3] of the coordinates of the point
    // Check whether the point is on the X extreme
    if ( ( pointIt.Value()[ 0 ] == 0 ) || ( pointIt.Value()[ 0 ] == ( domainSize[ 0 ] - 1 ) ) )
      {
      pointParameters.m_CanMoveX = false;
#if 0
      pointParameters.m_Alphas = borderAlphasEntry;
      pointParameters.m_CanChangeAlphas = false;
#endif
      }
    else
      {
      pointParameters.m_CanMoveX = true;
      }
    // Check whether the point is on the Y extreme    
    if ( ( pointIt.Value()[ 1 ] == 0 ) || ( pointIt.Value()[ 1 ] == ( domainSize[ 1 ] - 1 ) ) )
      {
      pointParameters.m_CanMoveY = false;
#if 0      
      pointParameters.m_Alphas = borderAlphasEntry;
      pointParameters.m_CanChangeAlphas = false;
#endif
      }
    else
      {
      pointParameters.m_CanMoveY = true;
      }
    // Check whether the point is on the Z extreme
    if ( ( pointIt.Value()[ 2 ] == 0 ) || ( pointIt.Value()[ 2 ] == ( domainSize[ 2 ] - 1 ) ) )
      {
      pointParameters.m_CanMoveZ = false;
#if 0      
      pointParameters.m_Alphas = borderAlphasEntry;
      pointParameters.m_CanChangeAlphas = false;
#endif
      }
    else
      {
      pointParameters.m_CanMoveZ = true;
      }
    // Insert the pointParameters struct into the mesh
    meshSource->GetOutput()->SetPointData( pointIt.Index(), pointParameters );
    }
    
  // Now simply create the mesh collection by repeating the mesh
  this->GenerateFromSingleMesh( meshSource->GetOutput(), numberOfMeshes, K );
}
  



//
//
//
bool
AtlasMeshCollection
::GetCollapsed( AtlasMesh::CellIdentifier edgeId,
                AtlasMeshCollection::Pointer& collapsed,
                std::set< AtlasMesh::CellIdentifier >& disappearingCells,
                AtlasMesh::CellIdentifier& unifiedVertexId, bool initializeAlphasToFlat ) const
{
  // Sanity check
  if ( !m_Cells->IndexExists( edgeId ) )
    return false;
    
  if ( m_Cells->ElementAt( edgeId )->GetType() !=  AtlasMesh::CellType::LINE_CELL )
    return false;
 
  // Retrieve the id's of the two vertices of the edge to be collapsed
  const AtlasMesh::CellType*  edge = m_Cells->ElementAt( edgeId );
  AtlasMesh::CellType::PointIdConstIterator  pointIt = edge->PointIdsBegin();
  AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
  ++pointIt;
  AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;
  
#if 0
  std::cout << "Collapsing edgeId: " << edgeId 
            << "\n     edgePoint0Id: " << edgePoint0Id 
            << "\n     edgePoint1Id: " << edgePoint1Id 
            << std::endl;  
#endif
            
  // Retrieve if the vertices can move in X, Y, and Z direction, and check if the
  // edge can be collapsed at all accordingly
  bool  canMoveX0 = m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveX;
  bool  canMoveY0 = m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveY;
  bool  canMoveZ0 = m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveZ;
  bool  canMoveX1 = m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveX;
  bool  canMoveY1 = m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveY;
  bool  canMoveZ1 = m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveZ;
  
  unsigned int numberOfRestrictionsOnEdgePoint0 = 0;
  if ( !canMoveX0 )
    numberOfRestrictionsOnEdgePoint0++;
  if ( !canMoveY0 )
    numberOfRestrictionsOnEdgePoint0++;
  if ( !canMoveZ0 )
    numberOfRestrictionsOnEdgePoint0++;
    
  unsigned int numberOfRestrictionsOnEdgePoint1 = 0;
  if ( !canMoveX1 )
    numberOfRestrictionsOnEdgePoint1++;
  if ( !canMoveY1 )
    numberOfRestrictionsOnEdgePoint1++;
  if ( !canMoveZ1 )
    numberOfRestrictionsOnEdgePoint1++;


  switch ( numberOfRestrictionsOnEdgePoint0 )
    {
    case 1:
      // The first point lies in a plane
      switch ( numberOfRestrictionsOnEdgePoint1 )
        {
        case 0:
          // The other point lies in the middle of the volume; no problem at all
          break;
        default:
          // The other point is somehow constraint. Make sure it lies in the same
          // plane (could be on edge or even corner point) as the first point
          if ( ( !canMoveX0 && canMoveX1 ) || ( !canMoveY0 && canMoveY1 ) || ( !canMoveZ0 && canMoveZ1 ) )
            {
            return false;
            }
        }
      break;
    case 2:
      // The first point lies on an edge
      switch ( numberOfRestrictionsOnEdgePoint1 )
        {
        case 0:
          // The other point lies in the middle of the volume; no problem at all
          break;
        case 1:
          // The other point lies in a plane. Make sure the first point lies in the
          // same plane
          if ( ( !canMoveX1 && canMoveX0 ) || ( !canMoveY1 && canMoveY0 ) || ( !canMoveZ1 && canMoveZ0 ) )
            {
            return false;
            }
          break;
        default:
          // The other point is lies on an edge or is even a corner point. Make sure they're both on
          // the same edge
          if ( ( !canMoveX0 && ( canMoveX1 || ( !canMoveZ0 && canMoveZ1 ) || ( !canMoveY0 && canMoveY1 ) ) ) ||
               ( !canMoveY0 && ( canMoveY1 || ( !canMoveZ0 && canMoveZ1 ) || ( !canMoveX0 && canMoveX1 ) ) ) ||
               ( !canMoveZ0 && ( canMoveZ1 || ( !canMoveX0 && canMoveX1 ) || ( !canMoveY0 && canMoveY1 ) ) ) )
            {
            return false;
            }
        }
      break;
    case 3:
      // Exclude collapsing two corner points
      if ( numberOfRestrictionsOnEdgePoint1 == 3 )
        {
        return false;
        }
    }

  //std::cout << "Decided that the edge can be collapsed" << std::endl;





  // Find all tetrahedra that have p0 or p1 as their vertex. Divide those into two categories:
  // the ones that will disappear (p0 *and* p1 as vertices), and the ones that have only one
  // (referred to as "affected tetrahedra" as those will actually change shape due to the edge
  //  collapse )

  // Create links back from points to cells.
  const AtlasMesh::CellLinksContainerPointer  cellLinks = this->GetCellLinks();

  // Intersection of all the cells that contain p0 or p1 will be cells that will disappear
  const std::set< AtlasMesh::CellIdentifier >&  cellsContainingEdgePoint0 = cellLinks->ElementAt( edgePoint0Id );
  const std::set< AtlasMesh::CellIdentifier >&  cellsContainingEdgePoint1 = cellLinks->ElementAt( edgePoint1Id );
  //std::cout << "   Found " << cellsContainingEdgePoint0.size()
  //          << " cells containing point with id " << edgePoint0Id << std::endl;
  //std::cout << "   Found " << cellsContainingEdgePoint0.size()
  //          << " cells containing point with id " << edgePoint1Id << std::endl;
  std::set< AtlasMesh::CellIdentifier >   initialDisappearingCells;
  std::set_intersection( cellsContainingEdgePoint0.begin(), cellsContainingEdgePoint0.end(),
                         cellsContainingEdgePoint1.begin(), cellsContainingEdgePoint1.end(),
                         std::inserter( initialDisappearingCells, initialDisappearingCells.begin() ) );

  // Loop over all the cells that are tagged for disappearing, and find the tetrahedra. For each such
  // tetrahedron, there are additional cells that will disappear: the one triangle containing p1 but
  // not p0, and its two edges that contain p1
  std::set< AtlasMesh::CellIdentifier >  extraDisappearingCells;
  for ( std::set< AtlasMesh::CellIdentifier >::const_iterator  disappearingIt = initialDisappearingCells.begin();
        disappearingIt != initialDisappearingCells.end(); ++disappearingIt )
    {
    // Get the cell
    const AtlasMesh::CellType*  cell = m_Cells->ElementAt( *disappearingIt );
  
    if ( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
      {
      continue;
      }

    // Retrieve ids of the points
    AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
    AtlasMesh::PointIdentifier  point0Id = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier  point1Id = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier  point2Id = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier  point3Id = *pointIt;

    // Look for the two points that are NOT p0 or p1
    AtlasMesh::PointIdentifier  firstOtherPointId;
    AtlasMesh::PointIdentifier  secondOtherPointId;
    std::vector< AtlasMesh::PointIdentifier >  pointIds;
    pointIds.push_back( point0Id );
    pointIds.push_back( point1Id );
    pointIds.push_back( point2Id );
    pointIds.push_back( point3Id );
    std::vector< AtlasMesh::PointIdentifier >::const_iterator it = pointIds.begin();
    for ( ; it != pointIds.end(); ++it )
      {
      if ( ( *it != edgePoint0Id ) && ( *it != edgePoint1Id ) )
        {
        firstOtherPointId = *it;
        break;
        }
      }
    for ( ++it; it != pointIds.end(); ++it )
      {
      if ( ( *it != edgePoint0Id ) && ( *it != edgePoint1Id ) )
        {
        secondOtherPointId = *it;
        break;
        }
      }

    //std::cout << "For disappearing tetrahdron with id " << *disappearingIt
    //          << ", we have for the two points that are NOT p0 or p1: "
    //          << firstOtherPointId << " and " << secondOtherPointId << std::endl;

    // For each of these other points, find the intersection between the cells containing
    // p1 and the cells containing the other point. There should only be one line in there,
    // and that line will disappear as well
    std::vector< AtlasMesh::PointIdentifier >  otherPointIds;
    otherPointIds.push_back( firstOtherPointId );
    otherPointIds.push_back( secondOtherPointId );
    for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  otherIt = otherPointIds.begin();
          otherIt != otherPointIds.end(); ++otherIt )
      {
      const std::set< AtlasMesh::CellIdentifier >&  cellsContainingOtherPoint
                      = cellLinks->ElementAt( *otherIt );
      //std::cout << "     Found " << cellsContainingOtherPoint.size()
      //          << " cells containing point with id " << *otherIt << std::endl;
      std::set< AtlasMesh::CellIdentifier >  cellsContainingOtherPointAndP1;
      std::set_intersection( cellsContainingOtherPoint.begin(), cellsContainingOtherPoint.end(),
                             cellsContainingEdgePoint1.begin(), cellsContainingEdgePoint1.end(),
                             std::inserter( cellsContainingOtherPointAndP1,
                                            cellsContainingOtherPointAndP1.begin() ) );
      
      for ( std::set< AtlasMesh::CellIdentifier >::const_iterator  it = cellsContainingOtherPointAndP1.begin();
            it != cellsContainingOtherPointAndP1.end(); ++it )
        {
        // Get the cell
        const AtlasMesh::CellType*  cell = m_Cells->ElementAt( *it );
      
        if ( cell->GetType() == AtlasMesh::CellType::LINE_CELL )
          {
          //std::cout << "        Found an extra disappearing line that contains both "
          //          << *otherIt << " and p1: " << *it << std::endl;
          extraDisappearingCells.insert( *it );
          break;
          }
        }  

      }  // End loop over each of the two "other" points


    // Find the intersection between the cells containing p1, the first other point, and the second
    // other point. There should only be one triangle in there, and that triangle will disappear as
    // well.
    const std::set< AtlasMesh::CellIdentifier >&  cellsContainingFirstOtherPoint
                      = cellLinks->ElementAt( firstOtherPointId );
    const std::set< AtlasMesh::CellIdentifier >&  cellsContainingSecondOtherPoint
                      = cellLinks->ElementAt( secondOtherPointId );
    std::set< AtlasMesh::CellIdentifier >  cellsContainingFirstOtherPointAndP1;
    std::set_intersection( cellsContainingFirstOtherPoint.begin(), cellsContainingFirstOtherPoint.end(),
                           cellsContainingEdgePoint1.begin(), cellsContainingEdgePoint1.end(),
                           std::inserter( cellsContainingFirstOtherPointAndP1,
                                          cellsContainingFirstOtherPointAndP1.begin() ) );
    std::set< AtlasMesh::CellIdentifier >  cellsContainingFirstAndSecondOtherPointAndP1;
    std::set_intersection( cellsContainingFirstOtherPointAndP1.begin(), cellsContainingFirstOtherPointAndP1.end(),
                           cellsContainingSecondOtherPoint.begin(), cellsContainingSecondOtherPoint.end(),
                           std::inserter( cellsContainingFirstAndSecondOtherPointAndP1,
                                          cellsContainingFirstAndSecondOtherPointAndP1.begin() ) );
      
    for ( std::set< AtlasMesh::CellIdentifier >::const_iterator  it = cellsContainingFirstAndSecondOtherPointAndP1.begin();
          it != cellsContainingFirstAndSecondOtherPointAndP1.end(); ++it )
      {
      // Get the cell
      const AtlasMesh::CellType*  cell = m_Cells->ElementAt( *it );

      if ( cell->GetType() == AtlasMesh::CellType::TRIANGLE_CELL )
        {
        //std::cout << "        Found an extra disappearing triangle that contains "
        //          << firstOtherPointId << ", " << secondOtherPointId << ", and p1: " << *it << std::endl;
        extraDisappearingCells.insert( *it );
        break;
        }
      }



    } // End loop over all disappearing tetrahedra


  // The vertex corresponding to edgePoint1Id is also disappearing
  for ( std::set< AtlasMesh::CellIdentifier >::const_iterator  it = cellsContainingEdgePoint1.begin();
        it != cellsContainingEdgePoint1.end(); ++it )
    {
    // Get the cell
    const AtlasMesh::CellType*  cell = m_Cells->ElementAt( *it );

    if ( cell->GetType() == AtlasMesh::CellType::VERTEX_CELL )
      {
      //std::cout << "        Found an extra disappearing vertex that contains p1" << std::endl;
      extraDisappearingCells.insert( *it );
      break;
      }
    }


  // Unify both the disappearingCells and the extraDisappearingCells
  disappearingCells.clear();
  std::set_union( initialDisappearingCells.begin(), initialDisappearingCells.end(),
                  extraDisappearingCells.begin(), extraDisappearingCells.end(),
                  std::inserter( disappearingCells, disappearingCells.begin() ) );
  
  

  // Look up the tetrahedra that will be affected by the edge collapse (i.e. those that have
  // either p0 or p1 as a vertex. You can do this by looping over all cells that contain either
  // p0 or p1, disregard the non-tetrahedral cells, and making sure the remaining tethrahedra
  // are not in the set of cells that will disappear.
  std::vector< AtlasMesh::CellIdentifier >   affectedTetrahedra;
  std::set< AtlasMesh::CellIdentifier >  cellsContainingP0OrP1;
  std::set_union( cellsContainingEdgePoint0.begin(), cellsContainingEdgePoint0.end(),
                  cellsContainingEdgePoint1.begin(), cellsContainingEdgePoint1.end(),
                  std::inserter( cellsContainingP0OrP1, cellsContainingP0OrP1.begin() ) );
  for ( std::set< AtlasMesh::CellIdentifier >::const_iterator  p0orP1It = cellsContainingP0OrP1.begin();
        p0orP1It != cellsContainingP0OrP1.end(); ++p0orP1It )
    {
    // Get the cell
    const AtlasMesh::CellType*  cell = m_Cells->ElementAt( *p0orP1It );
  
    if ( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
      {
      continue;
      }

    // Make sure this cell is not in the disappearingCells set
    if ( disappearingCells.find( *p0orP1It ) == disappearingCells.end() )
      {
      //std::cout << "Found a tetrahedron that is affected: " << *p0orP1It << std::endl;
      affectedTetrahedra.push_back( *p0orP1It );
      }


    } // End loop over all cells that containing either p0 or p1


 
 
    
      
    
  // For each mesh, propose three positions: either of the positions of the two edge points, and
  // their average. If none of the positions yields a valid (i.e. counter-clockwise triangles) mesh, the 
  // edge can not be collapsed.
  std::vector< AtlasMesh::PointType >  newPositionsOfUnifiedVertex;
  for ( unsigned int meshNumber=0; meshNumber < this->GetNumberOfMeshes()+1; meshNumber++ )
    {
    PointsContainerType::Pointer  thisPosition;
    if ( meshNumber < this->GetNumberOfMeshes() )
      {
      thisPosition = m_Positions[ meshNumber ];
      }
    else
      {
      thisPosition = m_ReferencePosition;
      }
    
#if 0
    std::cout << "Trying to find a valid position of the unified vertex for mesh number: " 
              << meshNumber << std::endl;
#endif
              
    // Retrieve the position of the two vertices of the edge to be collapsed
    AtlasMesh::PointType  edgePoint0 = thisPosition->ElementAt( edgePoint0Id );
    AtlasMesh::PointType  edgePoint1 = thisPosition->ElementAt( edgePoint1Id );
    
#if 0
    std::cout << "          edgePoint0: " << edgePoint0 << std::endl; 
    std::cout << "          edgePoint1: " << edgePoint1 << std::endl; 
#endif
    
    // Depending on the allowed movements of the two vertices, propose positions to be tried
    std::vector< AtlasMesh::PointType >  testPoints;
    if ( numberOfRestrictionsOnEdgePoint0 > numberOfRestrictionsOnEdgePoint1 )
      {
      testPoints.push_back( edgePoint0 );
      }
    else if ( numberOfRestrictionsOnEdgePoint0 < numberOfRestrictionsOnEdgePoint1 )
      {
      testPoints.push_back( edgePoint1 );
      }
    else
      {
      AtlasMesh::PointType  middlePoint;
      middlePoint[ 0 ] = ( edgePoint0[ 0 ] + edgePoint1[ 0 ] ) / 2;
      middlePoint[ 1 ] = ( edgePoint0[ 1 ] + edgePoint1[ 1 ] ) / 2;
      middlePoint[ 2 ] = ( edgePoint0[ 2 ] + edgePoint1[ 2 ] ) / 2;

      AtlasMesh::PointType  halfWayBetweenMiddlePointAndEdgePoint0;
      halfWayBetweenMiddlePointAndEdgePoint0[ 0 ] = edgePoint0[ 0 ] * 0.75 + edgePoint1[ 0 ] * 0.25;
      halfWayBetweenMiddlePointAndEdgePoint0[ 1 ] = edgePoint0[ 1 ] * 0.75 + edgePoint1[ 1 ] * 0.25;
      halfWayBetweenMiddlePointAndEdgePoint0[ 2 ] = edgePoint0[ 2 ] * 0.75 + edgePoint1[ 2 ] * 0.25;

      AtlasMesh::PointType  halfWayBetweenMiddlePointAndEdgePoint1;
      halfWayBetweenMiddlePointAndEdgePoint1[ 0 ] = edgePoint0[ 0 ] * 0.25 + edgePoint1[ 0 ] * 0.75;
      halfWayBetweenMiddlePointAndEdgePoint1[ 1 ] = edgePoint0[ 1 ] * 0.25 + edgePoint1[ 1 ] * 0.75;
      halfWayBetweenMiddlePointAndEdgePoint1[ 2 ] = edgePoint0[ 2 ] * 0.25 + edgePoint1[ 2 ] * 0.75;

      testPoints.push_back( middlePoint );
      testPoints.push_back( halfWayBetweenMiddlePointAndEdgePoint0 );
      testPoints.push_back( halfWayBetweenMiddlePointAndEdgePoint1 );
      testPoints.push_back( edgePoint0 );
      testPoints.push_back( edgePoint1 );
      }

    
    // Loop over all test points, and check if one gives us a valid mesh
    bool illegalPosition = false;
    unsigned int testPointNumber = 0;
    for ( ; testPointNumber < testPoints.size(); testPointNumber++ )
      {
      illegalPosition = false;
      AtlasMesh::PointType  newPoint = testPoints[ testPointNumber ];
      
      //std::cout << "AtlasMeshCollection:           --> testing newPoint: " << newPoint << std::endl;
      
      // Loop over all affected tetrahedra, and check validity under proposed position
      std::vector< AtlasMesh::CellIdentifier >::const_iterator  affectedIt = affectedTetrahedra.begin();
      while ( affectedIt != affectedTetrahedra.end() )
        {
        const AtlasMesh::CellType*  cell = m_Cells->ElementAt( *affectedIt );
        
        // Retrieve ids and original positions of the four vertices
        pointIt = cell->PointIdsBegin();
        AtlasMesh::PointIdentifier  point0Id = *pointIt;
        ++pointIt;
        AtlasMesh::PointIdentifier  point1Id = *pointIt;
        ++pointIt;
        AtlasMesh::PointIdentifier  point2Id = *pointIt;
        ++pointIt;
        AtlasMesh::PointIdentifier  point3Id = *pointIt;

        AtlasMesh::PointType  point0 = thisPosition->ElementAt( point0Id );
        AtlasMesh::PointType  point1 = thisPosition->ElementAt( point1Id );
        AtlasMesh::PointType  point2 = thisPosition->ElementAt( point2Id );
        AtlasMesh::PointType  point3 = thisPosition->ElementAt( point3Id );

        
        // We know that one of the vertices lies on the edge to be collapsed. Look it up
        // and change its position to the current test position
        if ( ( point0Id == edgePoint0Id ) || ( point0Id == edgePoint1Id ) )
          {
          point0 = newPoint;
          }
        else if ( ( point1Id == edgePoint0Id ) || ( point1Id == edgePoint1Id ) )
          {
          point1 = newPoint;
          }
        else if ( ( point2Id == edgePoint0Id ) || ( point2Id == edgePoint1Id ) )
          {
          point2 = newPoint;
          }
        else if ( ( point3Id == edgePoint0Id ) || ( point3Id == edgePoint1Id ) )
          {
          point3 = newPoint;
          }

        
        // Test if this tethrahedron is left-rotating. If it's not, stop looping over all affected tetrahedra.
        // Do this by calculating the volume of the tetrahedron; this should be positive.
        // In what follows, the matrix Lambda is the Jacobian of the transform from a standarized tetrahedron
        // ( ( 0 0 0 )^T, ( 1 0 0 )^T, ( 0 1 0 )^T, ( 0 0 1 )^T ), which has volume 1/6, to the actual tetrahedron
        const float  lambda11 = -point0[ 0 ] + point1[ 0 ];
        const float  lambda21 = -point0[ 1 ] + point1[ 1 ];
        const float  lambda31 = -point0[ 2 ] + point1[ 2 ];
        const float  lambda12 = -point0[ 0 ] + point2[ 0 ];
        const float  lambda22 = -point0[ 1 ] + point2[ 1 ];
        const float  lambda32 = -point0[ 2 ] + point2[ 2 ];
        const float  lambda13 = -point0[ 0 ] + point3[ 0 ];
        const float  lambda23 = -point0[ 1 ] + point3[ 1 ];
        const float  lambda33 = -point0[ 2 ] + point3[ 2 ];
        const float  volume = ( lambda11 * ( lambda22*lambda33 - lambda32*lambda23 ) 
                                - lambda12 * ( lambda21*lambda33 - lambda31*lambda23 )
                                + lambda13 * ( lambda21*lambda32 - lambda31*lambda22 ) ) / 6;
        if ( volume <= 0  )
          {
          //std::cout << "AtlasMeshCollection: volume of tetrahedron " << *affectedIt
          //          << " would be negative" << std::endl;
          illegalPosition = true;
          break;
          }
        
#if 1
        if ( meshNumber == this->GetNumberOfMeshes() )
          {
          const float badness = TetrahedronRadiusRatio( point0, point1, point2, point3 );
          //std::cout << "AtlasMeshCollection: badness of tetrahedron " << *affectedIt
          //          << " at reference position " << newPoint << " is " << badness << std::endl;
          // std::cout << "                      point0: " << point0 << std::endl;
          // std::cout << "                      point1: " << point1 << std::endl;
          // std::cout << "                      point2: " << point2 << std::endl;
          // std::cout << "                      point3: " << point3 << std::endl;
          if ( badness > 10 )
            {
            //std::cout << "referencePosition of tetrahedron " << *affectedIt
            //          << " would be too bad (" << badness << ")" << std::endl;
            illegalPosition = true;
            break;
            }

          }
#endif
          
        ++affectedIt;
        }
        

      // If this test point was successful, stop looking for other test points
      if ( !illegalPosition )
        {
        break;
        }
  
      } // End loop over test positions
  
    // If none of the test points was legal, return negative
    if ( illegalPosition )
      {
      //std::cout << "Couldn't find a suitable position. Collapsing unsuccessful" << std::endl;
      return false;
      }
      
    // Retrieve position
    newPositionsOfUnifiedVertex.push_back( testPoints[ testPointNumber ] );
    
    
    } // End loop over all meshes
     


  // Look up the vertex corresponding to edgePoint0Id
  for ( std::set< AtlasMesh::CellIdentifier >::const_iterator  it = cellsContainingEdgePoint0.begin();
        it != cellsContainingEdgePoint0.end(); ++it )
    {
    // Get the cell
    const AtlasMesh::CellType*  cell = m_Cells->ElementAt( *it );

    if ( cell->GetType() == AtlasMesh::CellType::VERTEX_CELL )
      {
      unifiedVertexId = *it;
      //std::cout << "        Found that unifiedVertexId is " << unifiedVertexId << std::endl;
      break;
      }
    }


  
  // Copy all positions of the points, except for the second point on the edge, which simply
  // disappears, and first point, which gets the previously determined new position
  //std::cout << "Creating positions" << std::endl;
  std::vector< PointsContainerType::Pointer >  collapsedPositions;
  PointsContainerType::Pointer   collapsedReferencePosition = nullptr;
  for ( unsigned int meshNumber = 0; meshNumber < this->GetNumberOfMeshes()+1; meshNumber++ )
    {
    PointsContainerType::Pointer  thisPosition;
    if ( meshNumber < this->GetNumberOfMeshes() )
      {
      thisPosition = m_Positions[ meshNumber ];
      }
    else
      {
      thisPosition = m_ReferencePosition;
      }
     
    PointsContainerType::Pointer  collapsedPosition = PointsContainerType::New();
    
    PointsContainerType::ConstIterator  positionIt = thisPosition->Begin();
    while ( positionIt !=  thisPosition->End() )
      {
      if ( positionIt.Index() ==  edgePoint0Id )
        {
        collapsedPosition->InsertElement( positionIt.Index(),  newPositionsOfUnifiedVertex[ meshNumber ] );
        }
      else if ( positionIt.Index() !=  edgePoint1Id )
        {
        collapsedPosition->InsertElement( positionIt.Index(),  positionIt.Value() );
        }      
    
      ++positionIt;  
      }

      
    if ( meshNumber < this->GetNumberOfMeshes() )
      {
      // Push back  
      collapsedPositions.push_back( collapsedPosition );
      }
    else
      {
      collapsedReferencePosition = collapsedPosition;
      }
      
    }
    
      
  // Copy all parameters of the points, except for the second point on the edge, which simply
  // disappears, and the first point on the edge, which will be the average alpha of the two
  // collapsing vertices. As for the moving options of the new point, this will be the most 
  // stringent of the two.
  //std::cout << "Creating point parameters" << std::endl;
  
  PointDataContainerType::Pointer  collapsedPointParameters = PointDataContainerType::New();
  PointDataContainerType::ConstIterator  pointParamIt = m_PointParameters->Begin();
  if ( !initializeAlphasToFlat )
    {
    while ( pointParamIt != m_PointParameters->End() )
      {
      if ( pointParamIt.Index() ==  edgePoint0Id )
        {
        PointParameters  pointParams;
        pointParams.m_Alphas = ( m_PointParameters->ElementAt( edgePoint0Id ).m_Alphas + 
                                m_PointParameters->ElementAt( edgePoint1Id ).m_Alphas ) / 2;
        pointParams.m_CanMoveX = ( m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveX && 
                                  m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveX );
        pointParams.m_CanMoveY = ( m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveY && 
                                  m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveY );
        pointParams.m_CanMoveZ = ( m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveZ &&
                                  m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveZ );

        //std::cout << "Point parameters of collapsed point: " << std::endl;
        //std::cout << "     Alphas: " << pointParams.m_Alphas << std::endl;
        //std::cout << "     CanMoveX: " << pointParams.m_CanMoveX << std::endl;
        //std::cout << "     CanMoveY: " << pointParams.m_CanMoveY << std::endl;
        //std::cout << "     CanMoveZ: " << pointParams.m_CanMoveZ << std::endl;

        collapsedPointParameters->InsertElement( pointParamIt.Index(),  pointParams );
        }
      else if ( pointParamIt.Index() !=  edgePoint1Id )
        {
        // Copy
        collapsedPointParameters->InsertElement( pointParamIt.Index(),  pointParamIt.Value() );
        }
      
      ++pointParamIt;
      }
    }
  else
    {
    const unsigned int  numberOfClasses = pointParamIt.Value().m_Alphas.Size();
    AtlasAlphasType   flatAlphasEntry( numberOfClasses );
    flatAlphasEntry.Fill( 1.0f / static_cast< float >( numberOfClasses ) );
      
    while ( pointParamIt != m_PointParameters->End() )
      {
      if ( pointParamIt.Index() ==  edgePoint0Id )
        {
        PointParameters  pointParams;
        
        pointParams.m_CanChangeAlphas = ( m_PointParameters->ElementAt( edgePoint0Id ).m_CanChangeAlphas && 
                                  m_PointParameters->ElementAt( edgePoint1Id ).m_CanChangeAlphas );
        
        if ( pointParams.m_CanChangeAlphas )
          {
          pointParams.m_Alphas = flatAlphasEntry;
          }
        else
          {
          pointParams.m_Alphas = pointParamIt.Value().m_Alphas;
          }
          

        pointParams.m_CanMoveX = ( m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveX && 
                                  m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveX );
        pointParams.m_CanMoveY = ( m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveY && 
                                  m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveY );
        pointParams.m_CanMoveZ = ( m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveZ &&
                                  m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveZ );

        //std::cout << "Point parameters of collapsed point: " << std::endl;
        //std::cout << "     Alphas: " << pointParams.m_Alphas << std::endl;
        //std::cout << "     CanChangeAlphas: " << pointParams.m_CanChangeAlphas << std::endl;
        //std::cout << "     CanMoveX: " << pointParams.m_CanMoveX << std::endl;
        //std::cout << "     CanMoveY: " << pointParams.m_CanMoveY << std::endl;
        //std::cout << "     CanMoveZ: " << pointParams.m_CanMoveZ << std::endl;

        collapsedPointParameters->InsertElement( pointParamIt.Index(),  pointParams );
        }
      else if ( pointParamIt.Index() !=  edgePoint1Id )
        {
        // Copy
        PointParameters  pointParams;
        
        pointParams.m_CanChangeAlphas = pointParamIt.Value().m_CanChangeAlphas;
        
        if ( pointParams.m_CanChangeAlphas )
          {
          pointParams.m_Alphas = flatAlphasEntry;
          }
        else
          {
          pointParams.m_Alphas = pointParamIt.Value().m_Alphas;
          }
          
        pointParams.m_CanMoveX = pointParamIt.Value().m_CanMoveX;
        pointParams.m_CanMoveY = pointParamIt.Value().m_CanMoveY;
        pointParams.m_CanMoveZ = pointParamIt.Value().m_CanMoveZ;
        
        collapsedPointParameters->InsertElement( pointParamIt.Index(), pointParams );
        }
      
      
      ++pointParamIt;
      }
    }    


  // Now loop over all cells, and simply copy unless the cell is disappearing. Also, correct all
  // references to the second edge point (which is disappearing) to a reference to the first edge 
  // point 
  //std::cout << "Creating cells" << std::endl;
  CellsContainerType::Pointer  collapsedCells = CellsContainerType::New();
  for ( AtlasMesh::CellsContainer::Iterator cellIt = m_Cells->Begin();
        cellIt != m_Cells->End(); ++cellIt )
    {
    bool cellIsDisappearing = false;
    std::set< AtlasMesh::CellIdentifier >::const_iterator  disappearIt = disappearingCells.begin();
    while (  disappearIt != disappearingCells.end() )
      {
      if ( cellIt.Index() == *disappearIt )
        {
        cellIsDisappearing = true;
        break;
        }
      
      ++disappearIt;
      }

    if ( !cellIsDisappearing )
      {
      const AtlasMesh::CellType*  cell = cellIt.Value();
      
      // Create a new cell of the correct type
      typedef itk::VertexCell< AtlasMesh::CellType >    VertexCell;
      typedef itk::LineCell< AtlasMesh::CellType >      LineCell;
      typedef itk::TriangleCell< AtlasMesh::CellType >  TriangleCell;
      typedef itk::TetrahedronCell< AtlasMesh::CellType >  TetrahedronCell;

      AtlasMesh::CellAutoPointer newCell;
      
      if ( cell->GetType() == AtlasMesh::CellType::VERTEX_CELL )
        {
        // Create a new vertex cell
        newCell.TakeOwnership( new VertexCell );
        }
      else if ( cell->GetType() == AtlasMesh::CellType::LINE_CELL  )
        {
        // Create a new line cell
        newCell.TakeOwnership( new LineCell );
        }
      else if ( cell->GetType() == AtlasMesh::CellType::TRIANGLE_CELL  )
        {
        // Create a new triangle cell
        newCell.TakeOwnership( new TriangleCell );
        }
      else
        {
        // Create a new tetrahedron cell
        newCell.TakeOwnership( new TetrahedronCell );
        }
      
      
      // Add points to the new cell: Copy in most cases, correct instances of second collapse-edge point
      int localId = 0;
      pointIt = cell->PointIdsBegin();
      while ( pointIt != cell->PointIdsEnd() )
        {
        AtlasMesh::PointIdentifier  pointId = *pointIt;
        
        // Correct if necessary
        if ( pointId == edgePoint1Id )
          {
          pointId = edgePoint0Id;
          }
        
        // 
        newCell->SetPointId( localId, pointId );
        
        ++localId;     
        ++pointIt;
        }
                            
                                    
      // Add the cell
      collapsedCells->InsertElement( cellIt.Index(), newCell.ReleaseOwnership() );
      }
    
    }
    
    
  // Create a new mesh collection to hold the result.
  collapsed = AtlasMeshCollection::New();
  collapsed->SetPointParameters( collapsedPointParameters );
  collapsed->SetCells( collapsedCells );
  collapsed->SetReferencePosition( collapsedReferencePosition );
  collapsed->SetPositions( collapsedPositions );
  collapsed->SetK( m_K );

  
  // Done.
  return true;
}





//
//
//
AtlasMeshCollection::Pointer
AtlasMeshCollection
::GetRegionGrown( AtlasMesh::CellIdentifier seedId, unsigned int radius, bool makeOuterPointsImmobile ) const
{
  
  // Sanity check
  if ( !( m_Cells->IndexExists( seedId ) ) )
    {
    // Cell simply doesn't exist.
    return nullptr;
    }
 
  // Get the cell
  const AtlasMesh::CellType*  cell = m_Cells->ElementAt( seedId );

  
  /**
   *
   * Part I. Create all the cells of the region grown 
   *
   */
  //std::cout << "Creating cells" << std::endl;
   
  
  // Create links back from points to cells. 
  const AtlasMesh::CellLinksContainerPointer  cellLinks = this->GetCellLinks();
  
  
  // Initialize the set of outer vertices with the those belonging to the seed cell
  std::set< AtlasMesh::PointIdentifier >  outerPoints;
  AtlasMesh::CellType::PointIdConstIterator  pit = cell->PointIdsBegin();
  while ( pit != cell->PointIdsEnd() )
    {
    outerPoints.insert( *pit );
    ++pit;
    }
    
    
  // Initialize the region grown cells with the seed cell and all its subcells
  CellsContainerType::Pointer  regionGrownCells = CellsContainerType::New();
  
  std::set< AtlasMesh::PointIdentifier >::const_iterator  outerPointIt = outerPoints.begin();
  while ( outerPointIt != outerPoints.end() )
    {
    
    // Get all neigboring cells of this outer vertex
    const std::set< AtlasMesh::CellIdentifier >&  cellNeighbors = cellLinks->ElementAt( *outerPointIt );

    //std::cout << "   Found " << cellNeighbors.size()
    //          << " neighboring cells of outer point " << *outerPointIt << std::endl;

    // Loop over all neighboring cells
    std::set< AtlasMesh::CellIdentifier >::const_iterator  neighborIt = cellNeighbors.begin();
    while ( neighborIt != cellNeighbors.end() )
      {
      // Don't do anything if the cell was already added
      if ( regionGrownCells->IndexExists( *neighborIt ) )
        {
        ++neighborIt;
        continue;
        }
      
      // Get the cell
      const AtlasMesh::CellType*  cell = m_Cells->ElementAt( *neighborIt );
       
      // Don't do anything if the cell contains other points than the ones in the original seed cell
      bool  containsWrongPoints = false;
      AtlasMesh::CellType::PointIdConstIterator  pit = cell->PointIdsBegin();
      while ( pit != cell->PointIdsEnd() )
        {
        if ( std::find( outerPoints.begin(), outerPoints.end(), *pit ) == outerPoints.end() )
          {
          containsWrongPoints = true;
          break;
          }
        
        ++pit;
        }
      if ( containsWrongPoints )
        {
        ++neighborIt;
        continue;
        }
        
    
 
      // Create a new cell of the correct type
      typedef itk::VertexCell< AtlasMesh::CellType >    VertexCell;
      typedef itk::LineCell< AtlasMesh::CellType >      LineCell;
      typedef itk::TriangleCell< AtlasMesh::CellType >  TriangleCell;
      typedef itk::TetrahedronCell< AtlasMesh::CellType >  TetrahedronCell;

      AtlasMesh::CellAutoPointer newCell;
      if ( cell->GetType() == AtlasMesh::CellType::VERTEX_CELL )
        {
        // Create a new vertex cell
        newCell.TakeOwnership( new VertexCell );
        }
      else if ( cell->GetType() == AtlasMesh::CellType::LINE_CELL )
        {
        // Create a new line cell
        newCell.TakeOwnership( new LineCell );
        }
      else if ( cell->GetType() == AtlasMesh::CellType::TRIANGLE_CELL )
        {
        // Create a new triangle cell
        newCell.TakeOwnership( new TriangleCell );
        }
      else 
        {
        // Create a new tetrahedron cell
        newCell.TakeOwnership( new TetrahedronCell );
        }

      
      // Add points to the new cell
      newCell->SetPointIds( cell->PointIdsBegin(),
                            cell->PointIdsEnd() );

      // Add the new cell
      regionGrownCells->InsertElement( *neighborIt, newCell.ReleaseOwnership() );

      ++neighborIt;
      }

    ++outerPointIt;
    }

      
      


  
  for ( unsigned int iterationNumber = 0; iterationNumber < radius; iterationNumber++ )
    {
    //std::cout << "iterationNumber: " << iterationNumber << std::endl;
    
    // Loop over all outer points, retrieve their neighboring cells, and add those to the 
    // region grown mesh collection, taking advantage of the fact that no duplicates will be
    // made anyway in the itkMapContainer. Whenever a vertex is added, check if this is really
    // a new vertex - if it is, it will form part of the new outer boundary
    std::set< AtlasMesh::PointIdentifier >  newOuterPoints;
    std::set< AtlasMesh::PointIdentifier >::const_iterator  outerPointIt = outerPoints.begin();
    while ( outerPointIt != outerPoints.end() )
      {
      // Get all neigboring cells of this outer vertex
      const std::set< AtlasMesh::CellIdentifier >&  cellNeighbors = cellLinks->ElementAt( *outerPointIt );

      //std::cout << "   Found " << cellNeighbors.size()
      //          << " neighboring cells of outer point " << *outerPointIt << std::endl;

      // Loop over all neighboring cells
      std::set< AtlasMesh::CellIdentifier >::const_iterator  neighborIt = cellNeighbors.begin();
      while ( neighborIt != cellNeighbors.end() )
        {
        //std::cout << "    Inspecting cell with id " << *neighborIt << std::endl;

        if ( regionGrownCells->IndexExists( *neighborIt ) )
          {
          ++neighborIt;
          continue;
          }
        
        const AtlasMesh::CellType*  cell = m_Cells->ElementAt( *neighborIt );

        // Create a new cell of the correct type
        typedef itk::VertexCell< AtlasMesh::CellType >    VertexCell;
        typedef itk::LineCell< AtlasMesh::CellType >      LineCell;
        typedef itk::TriangleCell< AtlasMesh::CellType >  TriangleCell;
        typedef itk::TetrahedronCell< AtlasMesh::CellType >  TetrahedronCell;

        AtlasMesh::CellAutoPointer newCell;
        if ( cell->GetType() == AtlasMesh::CellType::LINE_CELL  )
          {
          // Create a new line cell
          newCell.TakeOwnership( new LineCell );
          
          // We know that one the points lies on the outer boundary. The other one is new
          // and will become part of the new outer boundary
          AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
          AtlasMesh::PointIdentifier  pointId0 = *pointIt;
          ++pointIt;
          AtlasMesh::PointIdentifier  pointId1 = *pointIt;
          //std::cout << "      LineCell wiht id " << *neighborIt << " gave rise to adding point with id ";
          if ( *outerPointIt == pointId0 )
            {
            newOuterPoints.insert( pointId1 );
            //std::cout << pointId1;
            }
          else
            {
            newOuterPoints.insert( pointId0 );
            //std::cout << pointId0;
            }
          //std::cout << " to the outer points." << std::endl;
          
          }
        else if ( cell->GetType() == AtlasMesh::CellType::TRIANGLE_CELL )
          {
          // Create a new triangle cell
          newCell.TakeOwnership( new TriangleCell );
          }
        else if ( cell->GetType() == AtlasMesh::CellType::TETRAHEDRON_CELL )
          {
          // Create a new tetrahedron cell
          newCell.TakeOwnership( new TetrahedronCell );
          }

        
        // Add points to the new cell
        newCell->SetPointIds( cell->PointIdsBegin(),
                              cell->PointIdsEnd() );

        // Add the new cell
        regionGrownCells->InsertElement( *neighborIt, newCell.ReleaseOwnership() );

        ++neighborIt;
        }

      ++outerPointIt;
      }
      
    
    
    
    // Also add the vertex ids of the new outer points, as well as lines connecting two outer points
    outerPointIt = newOuterPoints.begin();
    while ( outerPointIt != newOuterPoints.end() )
      {
      // Get all neigboring cells of this outer vertex
      const std::set< AtlasMesh::CellIdentifier >&  cellNeighbors = cellLinks->ElementAt( *outerPointIt );
       
      // Get the one neighbor that is of type VERTEX_CELL, and add it to the region grown cells
      std::set< AtlasMesh::CellIdentifier >::const_iterator  neighborIt = cellNeighbors.begin();
      while ( neighborIt != cellNeighbors.end() )
        {
        const AtlasMesh::CellType*  cell = m_Cells->ElementAt( *neighborIt );

        if ( cell->GetType() == AtlasMesh::CellType::VERTEX_CELL )
          {
          // Create a new cell of the correct type
          typedef itk::VertexCell< AtlasMesh::CellType >    VertexCell;
          AtlasMesh::CellAutoPointer newCell;
          newCell.TakeOwnership( new VertexCell );
        
          // Add points to the new cell
          newCell->SetPointIds( cell->PointIdsBegin(),
                                cell->PointIdsEnd() );

          // Add the new cell
          regionGrownCells->InsertElement( *neighborIt, newCell.ReleaseOwnership() );
          
          }
        else if ( cell->GetType() == AtlasMesh::CellType::LINE_CELL )
          {
          // Retrieve the other point of this line
          AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
          AtlasMesh::PointIdentifier  pointId0 = *pointIt;
          ++pointIt;
          AtlasMesh::PointIdentifier  pointId1 = *pointIt;
          AtlasMesh::PointIdentifier  otherPointId = pointId0;
          if ( *outerPointIt == pointId0 )
            {
            otherPointId = pointId1;
            }

          // If this other point is also on the outer border, we have to add a line cell
          if ( std::find( newOuterPoints.begin(), newOuterPoints.end(), otherPointId ) != 
               newOuterPoints.end() )
            {
            if ( !( regionGrownCells->IndexExists( *neighborIt ) ) )
              {
               
              //std::cout << "Adding line " << *neighborIt
              //          << " because it lies between two outer points" << std::endl;
              
              // Create a new line cell
              typedef itk::LineCell< AtlasMesh::CellType >    LineCell;
              AtlasMesh::CellAutoPointer newCell;
              newCell.TakeOwnership( new LineCell );

              // Add points to the new cell
              newCell->SetPointIds( cell->PointIdsBegin(),
                                    cell->PointIdsEnd() );

              // Add the new cell
              regionGrownCells->InsertElement( *neighborIt, newCell.ReleaseOwnership() );
                       
              
              }
 
            }
            
                      
          }
        else if ( cell->GetType() == AtlasMesh::CellType::TRIANGLE_CELL )
          {
          // Retrieve the other two points of this triangle
          AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
          AtlasMesh::PointIdentifier  pointId0 = *pointIt;
          ++pointIt;
          AtlasMesh::PointIdentifier  pointId1 = *pointIt;
          ++pointIt;
          AtlasMesh::PointIdentifier  pointId2 = *pointIt;
          AtlasMesh::PointIdentifier  otherPointId1;
          AtlasMesh::PointIdentifier  otherPointId2;
          if ( *outerPointIt == pointId0 )
            {
            otherPointId1 = pointId1;
            otherPointId2 = pointId2;
            }
          else if  ( *outerPointIt == pointId1 )
            {
            otherPointId1 = pointId0;
            otherPointId2 = pointId2;
            }
          else
            {
            otherPointId1 = pointId0;
            otherPointId2 = pointId1;
            }

          // If the two other points are also on the outer border, we have to add a triangle cell
          if ( ( std::find( newOuterPoints.begin(), newOuterPoints.end(), otherPointId1 ) != newOuterPoints.end() ) &&
               ( std::find( newOuterPoints.begin(), newOuterPoints.end(), otherPointId2 ) != newOuterPoints.end() ) )
            {
            if ( !( regionGrownCells->IndexExists( *neighborIt ) ) )
              {
               
              //std::cout << "Adding triangle " << *neighborIt
              //          << " because it lies between three outer points" << std::endl;
              
              // Create a new line cell
              typedef itk::TriangleCell< AtlasMesh::CellType >    TriangleCell;
              AtlasMesh::CellAutoPointer newCell;
              newCell.TakeOwnership( new TriangleCell );

              // Add points to the new cell
              newCell->SetPointIds( cell->PointIdsBegin(),
                                    cell->PointIdsEnd() );

              // Add the new cell
              regionGrownCells->InsertElement( *neighborIt, newCell.ReleaseOwnership() );
                       
              
              }
 
            }
            
                      
          }
        
        ++neighborIt;
        }
      
      ++outerPointIt;
      }

      
    // Prepare for next iteration
    outerPoints = newOuterPoints;
    }
  
  //std::cout << "Created " << regionGrownCells->Size() << " cells." << std::endl;
    
#if 0  
  std::cout << "Outer points for mesh 0: " << std::endl;
  std::set< AtlasMesh::PointIdentifier >::const_iterator  outerPointIt = outerPoints.begin();
  while ( outerPointIt != outerPoints.end() )
    {
    std::cout << "    " << *outerPointIt << "   (" 
              << m_Positions[ 0 ]->ElementAt( *outerPointIt )[ 0 ]  << ", " 
              << m_Positions[ 0 ]->ElementAt( *outerPointIt )[ 1 ] << ")" << std::endl;
    
    ++outerPointIt;
    }
    
#endif
  
  
  /**
   *
   * Part II. Loop over all the cells of the region grown, and copy the appropriate elements
   *          of the mesh collection
   *
   */
  //std::cout << "Creating other mesh collection elements from the cells." << std::endl;
  
  // Create containers
  std::vector< PointsContainerType::Pointer >  regionGrownPositions;
  for ( unsigned int meshNumber = 0; meshNumber < this->GetNumberOfMeshes(); meshNumber++ )
    {
    regionGrownPositions.push_back( PointsContainerType::New() );
    }
  PointsContainerType::Pointer  regionGrownReferencePosition = PointsContainerType::New();
  PointDataContainerType::Pointer  regionGrownPointParameters = PointDataContainerType::New();

    
  // Loop over all cells
  CellsContainerType::ConstIterator  cellIt = regionGrownCells->Begin();
  CellsContainerType::ConstIterator  cellEnd = regionGrownCells->End();
  while ( cellIt != cellEnd )
    {
    AtlasMesh::CellType*  cell = cellIt.Value();
    
    if ( cell->GetType() == AtlasMesh::CellType::VERTEX_CELL )
      {
      // Retrieve point id
      AtlasMesh::PointIdentifier  pointId = *( cell->PointIdsBegin() );
      
      // Insert point for each mesh
      for ( unsigned int meshNumber = 0; meshNumber < this->GetNumberOfMeshes(); meshNumber++ )
        {
        regionGrownPositions[ meshNumber ]->InsertElement( pointId, 
                                                           m_Positions[ meshNumber ]->ElementAt( pointId ) );
        }

      // Insert point for reference position
      regionGrownReferencePosition->InsertElement( pointId, m_ReferencePosition->ElementAt( pointId ) );
      
      // Insert point parameters
      regionGrownPointParameters->InsertElement( pointId, m_PointParameters->ElementAt( pointId ) );
      
      }
    
    ++cellIt;
    }
  
  
  
  /**
   *
   * Part III. Make the outer points immobile if the user so desires
   *
   */ 
  if ( makeOuterPointsImmobile )
    {
    // Loop over all outer points
    std::set< AtlasMesh::PointIdentifier >::const_iterator  outerPointIt = outerPoints.begin();
    while ( outerPointIt != outerPoints.end() )
      {
      ( regionGrownPointParameters->ElementAt( *outerPointIt ) ).m_CanChangeAlphas = false;
      ( regionGrownPointParameters->ElementAt( *outerPointIt ) ).m_CanMoveX = false;
      ( regionGrownPointParameters->ElementAt( *outerPointIt ) ).m_CanMoveY = false;
      ( regionGrownPointParameters->ElementAt( *outerPointIt ) ).m_CanMoveZ = false;
      
      ++outerPointIt;
      }
    
    }
    
    
    
  // Collect everything and return
  AtlasMeshCollection::Pointer  regionGrown = AtlasMeshCollection::New();
  regionGrown->SetPointParameters( regionGrownPointParameters );
  regionGrown->SetCells( regionGrownCells );
  regionGrown->SetReferencePosition( regionGrownReferencePosition );
  regionGrown->SetPositions( regionGrownPositions );
  regionGrown->SetK( m_K );
  
  return regionGrown;

}

/*!
  \fn AtlasMeshCollection::PointerAtlasMeshCollection::GetUpsampled() const
  \brief Upsamples all meshes in the collection by dividing each
  tetrahedra into 5 smaller tetrahedra. Reconstruct all 8 corners of
  the original cube based on how the tessellation was originally done,
  split it into 2^3 = 8 subcubes, and full those subcubes with 5 small
  tetrahedra
*/
AtlasMeshCollection::Pointer
AtlasMeshCollection
::GetUpsampled() const
{
  // Look up what domain size is (DNG: why not just store it?)
  int  domainSize[] = { 0, 0, 0 };
  for ( AtlasMesh::PointsContainer::ConstIterator  it = this->GetReferencePosition()->Begin();
        it != this->GetReferencePosition()->End(); ++it )
    {
    for ( int i = 0; i < 3; i++ )
      {
      if ( it.Value()[ i ] > domainSize[ i ] )
        {
        domainSize[ i ] = static_cast< int >( it.Value()[ i ] + 0.5 );
        }
      }
    }
  for ( int i = 0; i < 3; i++ )
    {
    domainSize[ i ]++;
    }

  std::cout << "Domain size: " << std::endl;
  for ( int i = 0; i < 3; i++ )
    {
    std::cout << domainSize[ i ] << " ";
    }
  std::cout << std::endl;


  // Loop over all meshes, including reference
  std::vector< PointsContainerType::Pointer >  upsampledPositions;
  PointsContainerType::Pointer   upsampledReferencePosition = nullptr;
  CellsContainerType::Pointer  upsampledCells;
  for ( unsigned int meshNumber=0; meshNumber < this->GetNumberOfMeshes()+1; meshNumber++ )
    {
    //
    PointsContainerType::Pointer  thisPosition;
    if ( meshNumber < this->GetNumberOfMeshes() )
      {
      thisPosition = m_Positions[ meshNumber ];
      }
    else
      {
      thisPosition = m_ReferencePosition;
      }
    
    // Use a mesh source to create a mesh, simply splitting every tetrahedron in the original mesh into seven sub-tetrahedra
    typedef itk::AutomaticTopologyMeshSource< AtlasMesh >  MeshSourceType;
    MeshSourceType::Pointer  meshSource = MeshSourceType::New();

#if 0
    CellsContainerType::ConstIterator  cellIt = m_Cells->Begin();
    while ( cellIt != m_Cells->End() )
      {
      const AtlasMesh::CellType*  cell = cellIt.Value();
  
      if( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
        ++cellIt;
        continue;
        }
  
      // Retrieve ids of the points
      AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
      AtlasMesh::PointIdentifier  point0Id = *pointIt;
      ++pointIt;
      AtlasMesh::PointIdentifier  point1Id = *pointIt;
      ++pointIt;
      AtlasMesh::PointIdentifier  point2Id = *pointIt;
      ++pointIt;
      AtlasMesh::PointIdentifier  point3Id = *pointIt;
      
      // Retrieve positions of the points
      AtlasMesh::PointType  point0 = thisPosition->ElementAt( point0Id );
      AtlasMesh::PointType  point1 = thisPosition->ElementAt( point1Id );
      AtlasMesh::PointType  point2 = thisPosition->ElementAt( point2Id );
      AtlasMesh::PointType  point3 = thisPosition->ElementAt( point3Id );

      // Add seven subtetrahedra
      AtlasMesh::PointType  point01;
      point01[ 0 ] = ( point0[ 0 ] + point1[ 0 ] ) / 2;
      point01[ 1 ] = ( point0[ 1 ] + point1[ 1 ] ) / 2;
      point01[ 2 ] = ( point0[ 2 ] + point1[ 2 ] ) / 2;
      AtlasMesh::PointType  point12;
      point12[ 0 ] = ( point1[ 0 ] + point2[ 0 ] ) / 2;
      point12[ 1 ] = ( point1[ 1 ] + point2[ 1 ] ) / 2;
      point12[ 2 ] = ( point1[ 2 ] + point2[ 2 ] ) / 2;
      AtlasMesh::PointType  point02;
      point02[ 0 ] = ( point0[ 0 ] + point2[ 0 ] ) / 2;
      point02[ 1 ] = ( point0[ 1 ] + point2[ 1 ] ) / 2;
      point02[ 2 ] = ( point0[ 2 ] + point2[ 2 ] ) / 2;
      AtlasMesh::PointType  point03;
      point03[ 0 ] = ( point0[ 0 ] + point3[ 0 ] ) / 2;
      point03[ 1 ] = ( point0[ 1 ] + point3[ 1 ] ) / 2;
      point03[ 2 ] = ( point0[ 2 ] + point3[ 2 ] ) / 2;
      AtlasMesh::PointType  point13;
      point13[ 0 ] = ( point1[ 0 ] + point3[ 0 ] ) / 2;
      point13[ 1 ] = ( point1[ 1 ] + point3[ 1 ] ) / 2;
      point13[ 2 ] = ( point1[ 2 ] + point3[ 2 ] ) / 2;
      AtlasMesh::PointType  point23;
      point23[ 0 ] = ( point2[ 0 ] + point3[ 0 ] ) / 2;
      point23[ 1 ] = ( point2[ 1 ] + point3[ 1 ] ) / 2;
      point23[ 2 ] = ( point2[ 2 ] + point3[ 2 ] ) / 2;

      
      meshSource->AddTetrahedron( meshSource->AddPoint( point1 ),
                                  meshSource->AddPoint( point12 ),
                                  meshSource->AddPoint( point01 ),
                                  meshSource->AddPoint( point13 ) );
      meshSource->AddTetrahedron( meshSource->AddPoint( point01 ),
                                  meshSource->AddPoint( point12 ),
                                  meshSource->AddPoint( point0 ),
                                  meshSource->AddPoint( point13 ) );
      meshSource->AddTetrahedron( meshSource->AddPoint( point13 ),
                                  meshSource->AddPoint( point23 ),
                                  meshSource->AddPoint( point12 ),
                                  meshSource->AddPoint( point0 ) );
      meshSource->AddTetrahedron( meshSource->AddPoint( point0 ),
                                  meshSource->AddPoint( point12 ),
                                  meshSource->AddPoint( point02 ),
                                  meshSource->AddPoint( point23 ) );
      meshSource->AddTetrahedron( meshSource->AddPoint( point02 ),
                                  meshSource->AddPoint( point12 ),
                                  meshSource->AddPoint( point2 ),
                                  meshSource->AddPoint( point23 ) );
      meshSource->AddTetrahedron( meshSource->AddPoint( point13 ),
                                  meshSource->AddPoint( point23 ),
                                  meshSource->AddPoint( point0 ),
                                  meshSource->AddPoint( point03 ) );
      meshSource->AddTetrahedron( meshSource->AddPoint( point13 ),
                                  meshSource->AddPoint( point23 ),
                                  meshSource->AddPoint( point03 ),
                                  meshSource->AddPoint( point3 ) );
      
      
      // Go to next tetrahedron
      ++cellIt;
      } // End loop over all tetrahedra
#else


    float  precision[ 3 ];
    for ( int i = 0; i < 3; i++ )
      {
      precision[ i ] = 1e5 / static_cast< float >( domainSize[ i ] - 1 );
      }


    // Loop over all tetrahedra, making use of our prior knowledge about how these
    // tetrahedra were constructed in the first place. From
    // that information, we can reconstruct all 8 corners of the cube, split it into
    // 2^3 = 8 subcubes, and full those subcubes with 5 small tetrahedra
    int  counter = 0;
    AtlasMesh::PointType  p0;
    AtlasMesh::PointType  p1;
    AtlasMesh::PointType  p2;
    AtlasMesh::PointType  p3;
    AtlasMesh::PointType  p4;
    AtlasMesh::PointType  p5;
    AtlasMesh::PointType  p6;
    AtlasMesh::PointType  p7;
    for ( CellsContainerType::ConstIterator  cellIt = m_Cells->Begin();
          cellIt != m_Cells->End(); ++cellIt )
      {
      const AtlasMesh::CellType*  cell = cellIt.Value();
  
      if ( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
        continue;
        }

      AtlasMesh::CellType::PointIdConstIterator  pit = cell->PointIdsBegin();
      switch ( counter )
        {
        case 0:
          // Retrieve p0 as the first point in the tetrahedron
          p0 = thisPosition->ElementAt( *pit );
          
          // Retrieve p2 as the fourth point in the tetrahedron
          for ( int i = 0; i < 3; i++ )
            {
            ++pit;
            }
          p2 = thisPosition->ElementAt( *pit );
          
          break;
        case 1:
          // Retrieve p7 as the first point in the tetrahedron
          p7 = thisPosition->ElementAt( *pit );

          break;
        case 2:
          // Retrieve p1 as the third point in the tetrahedron
          for ( int i = 0; i < 2; i++ )
            {
            ++pit;
            }
          p1 = thisPosition->ElementAt( *pit );

          // Retrieve p5 as the fourth point in the tetrahedron
          ++pit;
          p5 = thisPosition->ElementAt( *pit );

          break; 
        case 3:
          // Retrieve p6 as the third point in the tetrahedron
          for ( int i = 0; i < 2; i++ )
            {
            ++pit;
            }
          p6 = thisPosition->ElementAt( *pit );

          // Retrieve p4 as the fourth point in the tetrahedron
          ++pit;
          p4 = thisPosition->ElementAt( *pit );

          break; 
        default:
          // Retrieve p3 as the third point in the tetrahedron
          for ( int i = 0; i < 2; i++ )
            {
            ++pit;
            }
          p3 = thisPosition->ElementAt( *pit );
          break; 
        }


      if ( counter != 4 )
        {
        counter++;
        continue;
        }

      // If we are here, it means we have encountered all 5 tetrahedra belonging to one cube,
      // and we have collected all corner coordinates
      // Reset the counter for the next cube, and carry on with the interesting stuff
      counter = 0;
      // std::cout << "Found cube with p0: " << p0 << std::endl;
      // std::cout << "                p1: " << p1 << std::endl;
      // std::cout << "                p2: " << p2 << std::endl;
      // std::cout << "                p3: " << p3 << std::endl;
      // std::cout << "                p4: " << p4 << std::endl;
      // std::cout << "                p5: " << p5 << std::endl;
      // std::cout << "                p6: " << p6 << std::endl;
      // std::cout << "                p7: " << p7 << std::endl;

#if 0
      // Construct points on the middle of each cube edge
      AtlasMesh::PointType  p01;
      for ( int i = 0; i < 3; i++ )
        {
        p01[ i ] = ( p0[ i ] + p1[ i ] ) / 2;
        }
      AtlasMesh::PointType  p02;
      for ( int i = 0; i < 3; i++ )
        {
        p02[ i ] = ( p0[ i ] + p2[ i ] ) / 2;
        }
      AtlasMesh::PointType  p13;
      for ( int i = 0; i < 3; i++ )
        {
        p13[ i ] = ( p1[ i ] + p3[ i ] ) / 2;
        }
      AtlasMesh::PointType  p23;
      for ( int i = 0; i < 3; i++ )
        {
        p23[ i ] = ( p2[ i ] + p3[ i ] ) / 2;
        }
      AtlasMesh::PointType  p15;
      for ( int i = 0; i < 3; i++ )
        {
        p15[ i ] = ( p1[ i ] + p5[ i ] ) / 2;
        }
      AtlasMesh::PointType  p37;
      for ( int i = 0; i < 3; i++ )
        {
        p37[ i ] = ( p3[ i ] + p7[ i ] ) / 2;
        }
      AtlasMesh::PointType  p26;
      for ( int i = 0; i < 3; i++ )
        {
        p26[ i ] = ( p2[ i ] + p6[ i ] ) / 2;
        }
      AtlasMesh::PointType  p04;
      for ( int i = 0; i < 3; i++ )
        {
        p04[ i ] = ( p0[ i ] + p4[ i ] ) / 2;
        }
      AtlasMesh::PointType  p57;
      for ( int i = 0; i < 3; i++ )
        {
        p57[ i ] = ( p5[ i ] + p7[ i ] ) / 2;
        }
      AtlasMesh::PointType  p67;
      for ( int i = 0; i < 3; i++ )
        {
        p67[ i ] = ( p6[ i ] + p7[ i ] ) / 2;
        }
      AtlasMesh::PointType  p46;
      for ( int i = 0; i < 3; i++ )
        {
        p46[ i ] = ( p4[ i ] + p6[ i ] ) / 2;
        }
      AtlasMesh::PointType  p45;
      for ( int i = 0; i < 3; i++ )
        {
        p45[ i ] = ( p4[ i ] + p5[ i ] ) / 2;
        }

      // Construct points in the middle of each cube face
      AtlasMesh::PointType  p0123;
      for ( int i = 0; i < 3; i++ )
        {
        p0123[ i ] = ( p0[ i ] + p1[ i ] + p2[ i ] + p3[ i ] ) / 4;
        }
      AtlasMesh::PointType  p1357;
      for ( int i = 0; i < 3; i++ )
        {
        p1357[ i ] = ( p1[ i ] + p3[ i ] + p5[ i ] + p7[ i ] ) / 4;
        }
      AtlasMesh::PointType  p2367;
      for ( int i = 0; i < 3; i++ )
        {
        p2367[ i ] = ( p2[ i ] + p3[ i ] + p6[ i ] + p7[ i ] ) / 4;
        }
      AtlasMesh::PointType  p0246;
      for ( int i = 0; i < 3; i++ )
        {
        p0246[ i ] = ( p0[ i ] + p2[ i ] + p4[ i ] + p6[ i ] ) / 4;
        }
      AtlasMesh::PointType  p0145;
      for ( int i = 0; i < 3; i++ )
        {
        p0145[ i ] = ( p0[ i ] + p1[ i ] + p4[ i ] + p5[ i ] ) / 4;
        }
      AtlasMesh::PointType  p4567;
      for ( int i = 0; i < 3; i++ )
        {
        p4567[ i ] = ( p4[ i ] + p5[ i ] + p6[ i ] + p7[ i ] ) / 4;
        }


      {
      for ( int i = 0; i < 3; i++ )
        {
        const float  precision = 1e0;

        p01[ i ] = static_cast< float >( static_cast< int >( p01[ i ] * precision + 0.5 ) ) / precision;
        p02[ i ] = static_cast< float >( static_cast< int >( p02[ i ] * precision + 0.5 ) ) / precision;
        p13[ i ] = static_cast< float >( static_cast< int >( p13[ i ] * precision + 0.5 ) ) / precision;
        p23[ i ] = static_cast< float >( static_cast< int >( p23[ i ] * precision + 0.5 ) ) / precision;
        p15[ i ] = static_cast< float >( static_cast< int >( p15[ i ] * precision + 0.5 ) ) / precision;
        p37[ i ] = static_cast< float >( static_cast< int >( p37[ i ] * precision + 0.5 ) ) / precision;
        p26[ i ] = static_cast< float >( static_cast< int >( p26[ i ] * precision + 0.5 ) ) / precision;
        p04[ i ] = static_cast< float >( static_cast< int >( p04[ i ] * precision + 0.5 ) ) / precision;
        p57[ i ] = static_cast< float >( static_cast< int >( p57[ i ] * precision + 0.5 ) ) / precision;
        p67[ i ] = static_cast< float >( static_cast< int >( p67[ i ] * precision + 0.5 ) ) / precision;
        p46[ i ] = static_cast< float >( static_cast< int >( p46[ i ] * precision + 0.5 ) ) / precision;
        p45[ i ] = static_cast< float >( static_cast< int >( p45[ i ] * precision + 0.5 ) ) / precision;

        p0123[ i ] = static_cast< float >( static_cast< int >( p0123[ i ] * precision + 0.5 ) ) / precision;
        p1357[ i ] = static_cast< float >( static_cast< int >( p1357[ i ] * precision + 0.5 ) ) / precision;
        p2367[ i ] = static_cast< float >( static_cast< int >( p2367[ i ] * precision + 0.5 ) ) / precision;
        p0246[ i ] = static_cast< float >( static_cast< int >( p0246[ i ] * precision + 0.5 ) ) / precision;
        p0145[ i ] = static_cast< float >( static_cast< int >( p0145[ i ] * precision + 0.5 ) ) / precision;
        p4567[ i ] = static_cast< float >( static_cast< int >( p4567[ i ] * precision + 0.5 ) ) / precision;
        }


      }




      // Construct point in the middle of the cube
      AtlasMesh::PointType  pMiddle;
      for ( int i = 0; i < 3; i++ )
        {
        pMiddle[ i ] = ( p0[ i ] + p1[ i ] + p2[ i ] + p3[ i ] + p4[ i ] + p5[ i ] + p6[ i ] + p7[ i ] ) / 8;
        }


      // Now fill each subcube with 5 tetrahedra
      this->FillCubeWithTetrahedra( meshSource, false,
                                    p0, p01, p02, p0123,
                                    p04, p0145, p0246, pMiddle );
      this->FillCubeWithTetrahedra( meshSource, true,
                                    p01, p1, p0123, p13,
                                    p0145, p15, pMiddle, p1357 );
      this->FillCubeWithTetrahedra( meshSource, true,
                                    p02, p0123, p2, p23,
                                    p0246, pMiddle, p26, p2367 );
      this->FillCubeWithTetrahedra( meshSource, false,
                                    p0123, p13, p23, p3,
                                    pMiddle, p1357, p2367, p37 );
      this->FillCubeWithTetrahedra( meshSource, true,
                                    p04, p0145, p0246, pMiddle,
                                    p4, p45, p46, p4567 );
      this->FillCubeWithTetrahedra( meshSource, false,
                                    p0145, p15, pMiddle, p1357,
                                    p45, p5, p4567, p57 );
      this->FillCubeWithTetrahedra( meshSource, false,
                                    p0246, pMiddle, p26, p2367,
                                    p46, p4567, p6, p67 );
      this->FillCubeWithTetrahedra( meshSource, true,
                                    pMiddle, p1357, p2367, p37,
                                    p4567, p57, p67, p7 );


#else

      int  integerP0[ 3 ];
      int  integerP1[ 3 ];
      int  integerP2[ 3 ];
      int  integerP3[ 3 ];
      int  integerP4[ 3 ];
      int  integerP5[ 3 ];
      int  integerP6[ 3 ];
      int  integerP7[ 3 ];
      for ( int i = 0; i < 3; i++ )
        {
        integerP0[ i ] = static_cast< int >( p0[ i ] * precision[ i ] + 0.5 );
        integerP1[ i ] = static_cast< int >( p1[ i ] * precision[ i ] + 0.5 );
        integerP2[ i ] = static_cast< int >( p2[ i ] * precision[ i ] + 0.5 );
        integerP3[ i ] = static_cast< int >( p3[ i ] * precision[ i ] + 0.5 );
        integerP4[ i ] = static_cast< int >( p4[ i ] * precision[ i ] + 0.5 );
        integerP5[ i ] = static_cast< int >( p5[ i ] * precision[ i ] + 0.5 );
        integerP6[ i ] = static_cast< int >( p6[ i ] * precision[ i ] + 0.5 );
        integerP7[ i ] = static_cast< int >( p7[ i ] * precision[ i ] + 0.5 );
        }


      // Construct points on the middle of each cube edge
      AtlasMesh::PointType  p01;
      for ( int i = 0; i < 3; i++ )
        {
        p01[ i ] = ( integerP0[ i ] + integerP1[ i ] ) / 2;
        }
      AtlasMesh::PointType  p02;
      for ( int i = 0; i < 3; i++ )
        {
        p02[ i ] = ( integerP0[ i ] + integerP2[ i ] ) / 2;
        }
      AtlasMesh::PointType  p13;
      for ( int i = 0; i < 3; i++ )
        {
        p13[ i ] = ( integerP1[ i ] + integerP3[ i ] ) / 2;
        }
      AtlasMesh::PointType  p23;
      for ( int i = 0; i < 3; i++ )
        {
        p23[ i ] = ( integerP2[ i ] + integerP3[ i ] ) / 2;
        }
      AtlasMesh::PointType  p15;
      for ( int i = 0; i < 3; i++ )
        {
        p15[ i ] = ( integerP1[ i ] + integerP5[ i ] ) / 2;
        }
      AtlasMesh::PointType  p37;
      for ( int i = 0; i < 3; i++ )
        {
        p37[ i ] = ( integerP3[ i ] + integerP7[ i ] ) / 2;
        }
      AtlasMesh::PointType  p26;
      for ( int i = 0; i < 3; i++ )
        {
        p26[ i ] = ( integerP2[ i ] + integerP6[ i ] ) / 2;
        }
      AtlasMesh::PointType  p04;
      for ( int i = 0; i < 3; i++ )
        {
        p04[ i ] = ( integerP0[ i ] + integerP4[ i ] ) / 2;
        }
      AtlasMesh::PointType  p57;
      for ( int i = 0; i < 3; i++ )
        {
        p57[ i ] = ( integerP5[ i ] + integerP7[ i ] ) / 2;
        }
      AtlasMesh::PointType  p67;
      for ( int i = 0; i < 3; i++ )
        {
        p67[ i ] = ( integerP6[ i ] + integerP7[ i ] ) / 2;
        }
      AtlasMesh::PointType  p46;
      for ( int i = 0; i < 3; i++ )
        {
        p46[ i ] = ( integerP4[ i ] + integerP6[ i ] ) / 2;
        }
      AtlasMesh::PointType  p45;
      for ( int i = 0; i < 3; i++ )
        {
        p45[ i ] = ( integerP4[ i ] + integerP5[ i ] ) / 2;
        }

      // Construct points in the middle of each cube face
      AtlasMesh::PointType  p0123;
      for ( int i = 0; i < 3; i++ )
        {
        p0123[ i ] = ( integerP0[ i ] + integerP1[ i ] + integerP2[ i ] + integerP3[ i ] ) / 4;
        }
      AtlasMesh::PointType  p1357;
      for ( int i = 0; i < 3; i++ )
        {
        p1357[ i ] = ( integerP1[ i ] + integerP3[ i ] + integerP5[ i ] + integerP7[ i ] ) / 4;
        }
      AtlasMesh::PointType  p2367;
      for ( int i = 0; i < 3; i++ )
        {
        p2367[ i ] = ( integerP2[ i ] + integerP3[ i ] + integerP6[ i ] + integerP7[ i ] ) / 4;
        }
      AtlasMesh::PointType  p0246;
      for ( int i = 0; i < 3; i++ )
        {
        p0246[ i ] = ( integerP0[ i ] + integerP2[ i ] + integerP4[ i ] + integerP6[ i ] ) / 4;
        }
      AtlasMesh::PointType  p0145;
      for ( int i = 0; i < 3; i++ )
        {
        p0145[ i ] = ( integerP0[ i ] + integerP1[ i ] + integerP4[ i ] + integerP5[ i ] ) / 4;
        }
      AtlasMesh::PointType  p4567;
      for ( int i = 0; i < 3; i++ )
        {
        p4567[ i ] = ( integerP4[ i ] + integerP5[ i ] + integerP6[ i ] + integerP7[ i ] ) / 4;
        }


      //
      for ( int i = 0; i < 3; i++ )
        {
        p0[ i ] = integerP0[ i ];
        p1[ i ] = integerP1[ i ];
        p2[ i ] = integerP2[ i ];
        p3[ i ] = integerP3[ i ];
        p4[ i ] = integerP4[ i ];
        p5[ i ] = integerP5[ i ];
        p6[ i ] = integerP6[ i ];
        p7[ i ] = integerP7[ i ];
        }

      // std::cout << "Found cube with p0: " << p0 << std::endl;
      // std::cout << "                p1: " << p1 << std::endl;
      // std::cout << "                p2: " << p2 << std::endl;
      // std::cout << "                p3: " << p3 << std::endl;
      // std::cout << "                p4: " << p4 << std::endl;
      // std::cout << "                p5: " << p5 << std::endl;
      // std::cout << "                p6: " << p6 << std::endl;
      // std::cout << "                p7: " << p7 << std::endl;

      // Construct point in the middle of the cube
      AtlasMesh::PointType  pMiddle;
      for ( int i = 0; i < 3; i++ )
        {
        pMiddle[ i ] = ( integerP0[ i ] + integerP1[ i ] + integerP2[ i ] + integerP3[ i ] + integerP4[ i ] + integerP5[ i ] + integerP6[ i ] + integerP7[ i ] ) / 8;
        }


      // Now fill each subcube with 5 tetrahedra
      this->FillCubeWithTetrahedra( meshSource, false,
                                    p0, p01, p02, p0123,
                                    p04, p0145, p0246, pMiddle );
      this->FillCubeWithTetrahedra( meshSource, true,
                                    p01, p1, p0123, p13,
                                    p0145, p15, pMiddle, p1357 );
      this->FillCubeWithTetrahedra( meshSource, true,
                                    p02, p0123, p2, p23,
                                    p0246, pMiddle, p26, p2367 );
      this->FillCubeWithTetrahedra( meshSource, false,
                                    p0123, p13, p23, p3,
                                    pMiddle, p1357, p2367, p37 );
      this->FillCubeWithTetrahedra( meshSource, true,
                                    p04, p0145, p0246, pMiddle,
                                    p4, p45, p46, p4567 );
      this->FillCubeWithTetrahedra( meshSource, false,
                                    p0145, p15, pMiddle, p1357,
                                    p45, p5, p4567, p57 );
      this->FillCubeWithTetrahedra( meshSource, false,
                                    p0246, pMiddle, p26, p2367,
                                    p46, p4567, p6, p67 );
      this->FillCubeWithTetrahedra( meshSource, true,
                                    pMiddle, p1357, p2367, p37,
                                    p4567, p57, p67, p7 );



#endif




#if 0
      {
      // Show all points in the current 
      for ( PointsContainerType::ConstIterator  it = meshSource->GetOutput()->GetPoints()->Begin();
            it != meshSource->GetOutput()->GetPoints()->End(); ++ it )
        {
        std::cout << "     Having now " << it.Index() << " -> " << it.Value() << std::endl;
        }

//       if ( ( meshNumber == 0 ) && ( meshSource->GetOutput()->GetPoints()->Size() > 53 ) )
//         {
//         std::cout << meshSource->GetOutput()->GetPoints()->ElementAt( 49 ) << std::endl;
//         std::cout << meshSource->GetOutput()->GetPoints()->ElementAt( 52 ) << std::endl;
//         std::cout << ( meshSource->GetOutput()->GetPoints()->ElementAt( 49 ) ==
//                        meshSource->GetOutput()->GetPoints()->ElementAt( 52 ) ) << std::endl;
//         std::cout << ( meshSource->GetOutput()->GetPoints()->ElementAt( 49 )[ 0 ] ==
//                        meshSource->GetOutput()->GetPoints()->ElementAt( 52 )[ 0 ] ) << std::endl;
//         std::cout << ( meshSource->GetOutput()->GetPoints()->ElementAt( 49 )[ 1 ] ==
//                        meshSource->GetOutput()->GetPoints()->ElementAt( 52 )[ 1 ] ) << std::endl;
//         std::cout << ( meshSource->GetOutput()->GetPoints()->ElementAt( 49 )[ 2 ] ==
//                        meshSource->GetOutput()->GetPoints()->ElementAt( 52 )[ 2 ] ) << std::endl;
//         std::cout << ( meshSource->GetOutput()->GetPoints()->ElementAt( 49 )[ 2 ] -
//                        meshSource->GetOutput()->GetPoints()->ElementAt( 52 )[ 2 ] ) << std::endl;
//        }

      }
#endif
      
      } // End loop over all tetrahedra

#if 1
    for ( PointsContainerType::Iterator  it = meshSource->GetOutput()->GetPoints()->Begin();
          it != meshSource->GetOutput()->GetPoints()->End(); ++ it )
      {
      it.Value()[ 0 ] /= precision [ 0 ];
      it.Value()[ 1 ] /= precision [ 1 ];
      it.Value()[ 2 ] /= precision [ 2 ];
      }
#endif

#endif
      
    // If this is the first mesh, remember the topology
    if ( meshNumber == 0 )
      {
      upsampledCells = meshSource->GetOutput()->GetCells();
      }
    
    // Remember the position
    if ( meshNumber < this->GetNumberOfMeshes() )
      {
      // Push back  
      upsampledPositions.push_back( meshSource->GetOutput()->GetPoints() );
      std::cout << "Having " << meshSource->GetOutput()->GetPoints()->Size() << " points for mesh number "
                << meshNumber << std::endl;
      }
    else
      {
      upsampledReferencePosition = meshSource->GetOutput()->GetPoints();
      std::cout << "Having " << meshSource->GetOutput()->GetPoints()->Size() << " points for reference mesh" << std::endl;
      }

    }  // End loop over all meshes
  

  // Double-check that all meshes have the same number of points. If not, than because of some weird numberical
  // errors, points are duplicated
  for ( unsigned int meshNumber=0; meshNumber < this->GetNumberOfMeshes(); meshNumber++ )
    {
    if ( upsampledPositions[ meshNumber ]->Size() != upsampledReferencePosition->Size() )
      {
      std::cerr << "Upsampling failed because of numerical inaccuracies!" << std::endl;
      return nullptr;
      }
    }
      
    
  // Assign flat alphas as a starting point. Vertices lying on the border can not move freely and belong to 
  // first class
  PointDataContainerType::Pointer  upsampledPointParameters = PointDataContainerType::New();
  const unsigned int  numberOfClasses = m_PointParameters->Begin().Value().m_Alphas.Size();
  
  AtlasAlphasType   flatAlphasEntry( numberOfClasses );
  flatAlphasEntry.Fill( 1.0f / static_cast< float >( numberOfClasses ) );
  
  AtlasAlphasType   borderAlphasEntry( numberOfClasses );
  borderAlphasEntry.Fill( 0.0f );
  borderAlphasEntry[ 0 ] = 1.0f;
  
  for ( AtlasMesh::PointsContainer::ConstIterator  pointIt = upsampledPositions[ 0 ]->Begin();
        pointIt != upsampledPositions[ 0 ]->End();
        ++pointIt )
    {
    AtlasMesh::PixelType  pointParameters;
    
    pointParameters.m_Alphas = flatAlphasEntry;
    pointParameters.m_CanChangeAlphas = true;
    
    if ( ( fabs( pointIt.Value()[ 0 ] ) < 1e-5 ) || ( fabs( pointIt.Value()[ 0 ] - ( domainSize[ 0 ] - 1 ) ) < 1e-5 ) )
      {
      pointParameters.m_CanMoveX = false;
#if 0
      pointParameters.m_Alphas = borderAlphasEntry;
      pointParameters.m_CanChangeAlphas = false;
#endif
      }
    else
      {
      pointParameters.m_CanMoveX = true;
      }
    
    if ( ( fabs( pointIt.Value()[ 1 ] ) < 1e-5 ) || ( fabs( pointIt.Value()[ 1 ] - ( domainSize[ 1 ] - 1 ) ) < 1e-5 ) )
      {
      pointParameters.m_CanMoveY = false;
#if 0
      pointParameters.m_Alphas = borderAlphasEntry;
      pointParameters.m_CanChangeAlphas = false;
#endif
      }
    else
      {
      pointParameters.m_CanMoveY = true;
      }

    if ( ( fabs( pointIt.Value()[ 2 ] ) < 1e-5 ) || ( fabs( pointIt.Value()[ 2 ] - ( domainSize[ 2 ] - 1 ) ) < 1e-5 ) )
      {
      pointParameters.m_CanMoveZ = false;
#if 0
      pointParameters.m_Alphas = borderAlphasEntry;
      pointParameters.m_CanChangeAlphas = false;
#endif
      }
    else
      {
      pointParameters.m_CanMoveZ = true;
      }
        
    upsampledPointParameters->InsertElement( pointIt.Index(), pointParameters );
    }

    
  // Collect everything and return
  AtlasMeshCollection::Pointer  upsampled = AtlasMeshCollection::New();
  upsampled->SetPointParameters( upsampledPointParameters );
  upsampled->SetCells( upsampledCells );
  upsampled->SetReferencePosition( upsampledReferencePosition );
  upsampled->SetPositions( upsampledPositions );
  upsampled->SetK( m_K );
  
  return upsampled;

}




//
//
//
AtlasMeshCollection::Pointer  
AtlasMeshCollection
::GetEdgeSplitted( AtlasMesh::CellIdentifier edgeId, 
                    AtlasMesh::CellIdentifier  newVertexId,
                    AtlasMesh::PointIdentifier  newPointId ) const
{
#if 0 
  AtlasMesh::CellIdentifier  oneOfTheTransverseEdgesIdDummy;
  AtlasMesh::CellIdentifier  transverseEdge0IdDummy;
  AtlasMesh::CellIdentifier  transverseEdge1IdDummy;

  return this->GetEdgeSplitted( edgeId, newVertexId, newPointId, transverseEdge0IdDummy, transverseEdge1IdDummy );
#else
  return nullptr;
#endif

}




//
//
//
AtlasMeshCollection::Pointer  
AtlasMeshCollection
::GetEdgeSwapped( AtlasMesh::CellIdentifier edgeId ) const
{
#if 0
  // Construct default new ids
  AtlasMesh::CellsContainer::ConstIterator  lastCellIt = m_Cells->End();
  lastCellIt--;
  AtlasMesh::PointsContainer::ConstIterator  lastPointIt = m_ReferencePosition->End();
  lastPointIt--;
  const AtlasMesh::CellIdentifier  newVertexId = lastCellIt.Index() + 1;
  const AtlasMesh::PointIdentifier  newPointId = lastPointIt.Index() + 1; 
    
  AtlasMesh::CellIdentifier  newEdgeIdDummy;
  return this->GetEdgeSwapped( edgeId, newVertexId, newPointId, newEdgeIdDummy );
#else
  return nullptr;
#endif
}




//
//
//
AtlasMeshCollection::Pointer  
AtlasMeshCollection
::GetEdgeSwapped( AtlasMesh::CellIdentifier edgeId, 
                   AtlasMesh::CellIdentifier  newVertexId,
                   AtlasMesh::PointIdentifier  newPointId,
                   AtlasMesh::CellIdentifier&  newEdgeId ) const
{

#if 0  
  // Sanity check
  if ( !m_Cells->IndexExists( edgeId ) )
    return 0;
    
  if ( m_Cells->ElementAt( edgeId )->GetType() !=  AtlasMesh::CellType::LINE_CELL )
    return 0;

  // Edges on borders are typically degenerate when swapped; simply forbid such swaps
  const AtlasMesh::CellType*  edge = m_Cells->ElementAt( edgeId );
  AtlasMesh::CellType::PointIdConstIterator  pointIt = edge->PointIdsBegin();
  AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
  ++pointIt;
  AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;
  bool  canMoveX0 = m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveX;
  bool  canMoveY0 = m_PointParameters->ElementAt( edgePoint0Id ).m_CanMoveY;
  bool  canMoveX1 = m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveX;
  bool  canMoveY1 = m_PointParameters->ElementAt( edgePoint1Id ).m_CanMoveY;

#if 0
  if ( !canMoveY0 && canMoveX1 )
    {
    return false;
    }
  if ( !canMoveX0 && canMoveY1 )
    {
    return false;
    }
  if ( !canMoveY1 && canMoveX0 )
    {
    return false;
    }
  if ( !canMoveX1 && canMoveY0 )
    {
    return false;
    }
#else
  if ( !canMoveX0 || !canMoveY0 || !canMoveX1 || !canMoveY1 )
    {
    return 0;
    }
#endif


  // Get edge splitted mesh, and remember the id's of the newly created transverse edges
  //std::cout << "Splitting edge " << edgeId << " (newVertexId: " << newVertexId 
  //          << ", newPointId: " << newPointId<< ")" << std::endl; 
  AtlasMesh::CellIdentifier  transverseEdge0Id;
  AtlasMeshCollection::Pointer  splitted = this->GetEdgeSplitted( edgeId, newVertexId, newPointId, 
                                                                    transverseEdge0Id, newEdgeId );
  if ( !splitted )
    {
    return 0;
    }

  //splitted->Write( "debugSplitted.txt" );


  // We are going to collapse one of the newly created transverse edges. But before we do that,
  // we have to make sure that the newly created vertex will move, not the other vertex. Therefore,
  // artificially set the other vertex to immoble before collapsing the edge. After collapsing, we
  // can restore its setting.
  pointIt = splitted->GetCells()->ElementAt( transverseEdge0Id )->PointIdsBegin();
  AtlasMesh::PointIdentifier  pointToTemporarilyImmbolizeId = *pointIt;
  if ( *pointIt == newPointId )
    {
    ++pointIt;
    pointToTemporarilyImmbolizeId = *pointIt;
    }
  const bool  canMoveXInPointToTemporarilyImmbolize = 
      splitted->GetPointParameters()->ElementAt( pointToTemporarilyImmbolizeId ).m_CanMoveX;
  const bool  canMoveYInPointToTemporarilyImmbolize = 
      splitted->GetPointParameters()->ElementAt( pointToTemporarilyImmbolizeId ).m_CanMoveY;
  splitted->GetPointParameters()->ElementAt( pointToTemporarilyImmbolizeId ).m_CanMoveX = false;
  splitted->GetPointParameters()->ElementAt( pointToTemporarilyImmbolizeId ).m_CanMoveY = false;
  

  // Now try to collapse one of the newly created transverse edges
  //std::cout << "Now trying to collapse edge with id: " << oneOfTheTransverseEdgesId << std::endl;
  AtlasMeshCollection::Pointer  collapsed;
  std::set< AtlasMesh::CellIdentifier >  disappearingCellsDummy;
  AtlasMesh::CellIdentifier  unifiedVertexId; 
  if ( !splitted->GetCollapsed( transverseEdge0Id, collapsed, disappearingCellsDummy, unifiedVertexId ) )
    {
    return 0;
    }

  //collapsed->Write( "debugCollapsed.txt" );
  
  // Restore the mobility settings of the unified vertex
  const AtlasMesh::PointIdentifier  mobilityToRestorePointId = 
      *( collapsed->GetCells()->ElementAt( unifiedVertexId )->PointIdsBegin() );
  collapsed->GetPointParameters()->ElementAt( mobilityToRestorePointId  ).m_CanMoveX = 
      canMoveXInPointToTemporarilyImmbolize;
  collapsed->GetPointParameters()->ElementAt( mobilityToRestorePointId ).m_CanMoveY = 
      canMoveYInPointToTemporarilyImmbolize;


  // If successful, return the result
  return collapsed;

#else
  return nullptr;

#endif
}                 




//
//
//
AtlasMeshCollection::Pointer
AtlasMeshCollection
::GetEdgeSplitted( AtlasMesh::CellIdentifier  edgeId, 
                    AtlasMesh::CellIdentifier  newVertexId,
                    AtlasMesh::PointIdentifier  newPointId,
                    AtlasMesh::CellIdentifier&  transverseEdge0Id,
                    AtlasMesh::CellIdentifier&  transverseEdge1Id ) const
{

#if 0

  /*
         p2                                   p2
        /|\                                  /|\
       / | \                                / ¦ \
  p0  /  |  \ p3            ======>    p0 /__|__\ p3
      \  |  /                              \  |p4/
       \ | /                                \ | /
        \|/                                  \|/
         p1                                   p1

  NOTE: in the remainder, p1 and p2 may be swapped!!!!

   */

  

  // Sanity check
  if ( !m_Cells->IndexExists( edgeId ) )
    return 0;
    
  if ( m_Cells->ElementAt( edgeId )->GetType() !=  AtlasMesh::CellType::LINE_CELL )
    return 0;


  // Initialize transverseEdgeId's
  transverseEdge0Id = 0;
  transverseEdge1Id = 0;

  // Get point ids of points 1 and 2, and create a point id for the new point 4
  const AtlasMesh::CellType*  edge = m_Cells->ElementAt( edgeId );
  AtlasMesh::CellType::PointIdConstIterator  pointIt = edge->PointIdsBegin();
  AtlasMesh::PointIdentifier  point1Id = *pointIt;
  ++pointIt;
  AtlasMesh::PointIdentifier  point2Id = *pointIt;
  AtlasMesh::PointIdentifier  point4Id = newPointId;
  

  // Create new reference position container, which is a simple copy of the original one
  // except for the newly added point 4
  PointsContainerType::Pointer  splittedReferencePosition = PointsContainerType::New();
  for ( PointsContainerType::ConstIterator refIt =  m_ReferencePosition->Begin(); 
        refIt != m_ReferencePosition->End(); ++refIt )
    {
    splittedReferencePosition->InsertElement( refIt.Index(), refIt.Value() );
    }
  AtlasMesh::PointType  referencePoint4;
  referencePoint4[ 0 ] = ( m_ReferencePosition->ElementAt( point1Id )[ 0 ] +
                            m_ReferencePosition->ElementAt( point2Id )[ 0 ] ) / 2;
  referencePoint4[ 1 ] = ( m_ReferencePosition->ElementAt( point1Id )[ 1 ] +
                            m_ReferencePosition->ElementAt( point2Id )[ 1 ] ) / 2;
  splittedReferencePosition->InsertElement( point4Id, referencePoint4 );


  // Do exactly the same for each of the mesh positions
  std::vector< PointsContainerType::Pointer >  splittedPositions;
  for ( unsigned int  meshNumber = 0; meshNumber < m_Positions.size(); meshNumber++ )
    {
    PointsContainerType::ConstPointer  thisPosition = m_Positions[ meshNumber ].GetPointer();
    PointsContainerType::Pointer  splittedPosition = PointsContainerType::New();
    for ( PointsContainerType::ConstIterator  posIt =  thisPosition->Begin(); 
          posIt != thisPosition->End(); ++posIt )
      {
      splittedPosition->InsertElement( posIt.Index(), posIt.Value() );
      }
    AtlasMesh::PointType  point4;
    point4[ 0 ] = ( thisPosition->ElementAt( point1Id )[ 0 ] +
                     thisPosition->ElementAt( point2Id )[ 0 ] ) / 2;
    point4[ 1 ] = ( thisPosition->ElementAt( point1Id )[ 1 ] +
                     thisPosition->ElementAt( point2Id )[ 1 ] ) / 2;
    splittedPosition->InsertElement( point4Id, point4 );
    
    splittedPositions.push_back( splittedPosition );
    }


  // Create a new point parameter container, which is a simple copy of the original one
  // except for the newly added point 4
  PointDataContainerType::Pointer  splittedPointParameters = PointDataContainerType::New();
  for ( PointDataContainerType::ConstIterator  paramIt = m_PointParameters->Begin();
        paramIt != m_PointParameters->End(); ++paramIt )
    {
    splittedPointParameters->InsertElement( paramIt.Index(), paramIt.Value() );
    }
  const bool  canMoveX1 = m_PointParameters->ElementAt( point1Id ).m_CanMoveX;
  const bool  canMoveY1 = m_PointParameters->ElementAt( point1Id ).m_CanMoveY;
  const bool  canMoveX2 = m_PointParameters->ElementAt( point2Id ).m_CanMoveX;
  const bool  canMoveY2 = m_PointParameters->ElementAt( point2Id ).m_CanMoveY;
  const bool  canMoveX4 = ( canMoveX1 || canMoveX2 );
  const bool  canMoveY4 = ( canMoveY1 || canMoveY2 );
  const bool  canChangeAlphas1 = m_PointParameters->ElementAt( point1Id ).m_CanChangeAlphas;
  const bool  canChangeAlphas2 = m_PointParameters->ElementAt( point2Id ).m_CanChangeAlphas;
  const bool  canChangeAlphas4 = ( canChangeAlphas1 || canChangeAlphas2 );
  AtlasAlphasType  alphas4 = m_PointParameters->ElementAt( point1Id ).m_Alphas;
  alphas4 += m_PointParameters->ElementAt( point2Id ).m_Alphas;
  alphas4 /= 2;
  PointParameters  pointParameters4;
  pointParameters4.m_Alphas = alphas4;
  pointParameters4.m_CanChangeAlphas = canChangeAlphas4;
  pointParameters4.m_CanMoveX = canMoveX4;
  pointParameters4.m_CanMoveY = canMoveY4;
  splittedPointParameters->InsertElement( point4Id, pointParameters4 );




  // Construct cells container for newly created mesh, and populate by looping over
  // cells from original cells container. At the same time, remember the ID's of the 
  // other points involved
  typedef itk::VertexCell< AtlasMesh::CellType >    VertexCell;
  typedef itk::LineCell< AtlasMesh::CellType >      LineCell;
  typedef itk::TriangleCell< AtlasMesh::CellType >  TriangleCell;
  AtlasMesh::CellIdentifier  vertex4Id = newVertexId;
  AtlasMesh::CellIdentifier  line14Id = newVertexId + 1;
  AtlasMesh::CellIdentifier  nextFreeCellId = newVertexId + 2;
  
  CellsContainerType::Pointer  splittedCells = CellsContainerType::New();

  AtlasMesh::CellAutoPointer  vertex4;
  vertex4.TakeOwnership( new VertexCell );
  vertex4->SetPointId( 0, point4Id );
  splittedCells->InsertElement( vertex4Id, vertex4.ReleaseOwnership() );

  AtlasMesh::CellAutoPointer  line14;
  line14.TakeOwnership( new LineCell );
  line14->SetPointId( 0, point1Id );
  line14->SetPointId( 1, point4Id );
  splittedCells->InsertElement( line14Id, line14.ReleaseOwnership() );


  for ( CellsContainerType::ConstIterator  cellIt = m_Cells->Begin();
         cellIt != m_Cells->End(); ++cellIt )
    {
    const AtlasMesh::CellType*  cell = cellIt.Value();
    
    if ( cell->GetType() == AtlasMesh::CellType::VERTEX_CELL )
      {
      // Copy vertex cell
      AtlasMesh::CellAutoPointer  newVertex;
      newVertex.TakeOwnership( new VertexCell );
      AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
      newVertex->SetPointId( 0, *pointIt );
      splittedCells->InsertElement( cellIt.Index(), newVertex.ReleaseOwnership() );
      }
    else if ( cell->GetType() == AtlasMesh::CellType::LINE_CELL  )
      {
      // Retrieve point ids of the vertices
      AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
      const AtlasMesh::PointIdentifier  idOfFirstPoint = *pointIt;
      ++pointIt;
      const AtlasMesh::PointIdentifier  idOfSecondPoint = *pointIt;

      if ( ( ( idOfFirstPoint == point1Id ) && ( idOfSecondPoint == point2Id ) ) ||
           ( ( idOfFirstPoint == point2Id ) && ( idOfSecondPoint == point1Id ) ) )
        {
        // 
        AtlasMesh::CellAutoPointer  newLine;
        newLine.TakeOwnership( new LineCell );
        newLine->SetPointId( 0, point2Id );
        newLine->SetPointId( 1, point4Id );
        splittedCells->InsertElement( cellIt.Index(), newLine.ReleaseOwnership() );
        }
      else
        {
        // Copy line cell
        AtlasMesh::CellAutoPointer  newLine;
        newLine.TakeOwnership( new LineCell );
        newLine->SetPointId( 0, idOfFirstPoint );
        newLine->SetPointId( 1, idOfSecondPoint );
        splittedCells->InsertElement( cellIt.Index(), newLine.ReleaseOwnership() );
        }

      }
    else
      {
      // Retrieve point ids of the vertices
      AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
      const AtlasMesh::PointIdentifier  idOfFirstPoint = *pointIt;
      ++pointIt;
      const AtlasMesh::PointIdentifier  idOfSecondPoint = *pointIt;
      ++pointIt;
      const AtlasMesh::PointIdentifier  idOfThirdPoint = *pointIt;

      const bool  configurationIsX12 = ( ( idOfSecondPoint == point1Id ) && ( idOfThirdPoint == point2Id ) );
      const bool  configurationIs12X = ( ( idOfFirstPoint == point1Id ) && ( idOfSecondPoint == point2Id ) );
      const bool  configurationIs2X1 = ( ( idOfFirstPoint == point2Id ) && ( idOfThirdPoint == point1Id ) );

      const bool  configurationIsX21 = ( ( idOfSecondPoint == point2Id ) && ( idOfThirdPoint == point1Id ) );
      const bool  configurationIs21X = ( ( idOfFirstPoint == point2Id ) && ( idOfSecondPoint == point1Id ) );
      const bool  configurationIs1X2 = ( ( idOfFirstPoint == point1Id ) && ( idOfThirdPoint == point2Id ) );

      if ( configurationIsX12 || configurationIs12X || configurationIs2X1 || 
           configurationIsX21 || configurationIs21X || configurationIs1X2 )
        {
        // Triangle 1
        AtlasMesh::CellAutoPointer  newTriangle1;
        newTriangle1.TakeOwnership( new TriangleCell );
        AtlasMesh::CellIdentifier  newTriangle1Id = nextFreeCellId;
        nextFreeCellId++;
        
        // Triangle 2
        AtlasMesh::CellAutoPointer  newTriangle2;
        newTriangle2.TakeOwnership( new TriangleCell );
        AtlasMesh::CellIdentifier  newTriangle2Id = nextFreeCellId;
        nextFreeCellId++;

        // Line
        AtlasMesh::CellAutoPointer  newLine;
        newLine.TakeOwnership( new LineCell );
        AtlasMesh::CellIdentifier  newLineId = nextFreeCellId;
        if ( transverseEdge0Id == 0 )
          {
          // First triangle that is split up
          transverseEdge0Id = newLineId;
          }
        else
          {
          // Second triangle that is split up - if it exists!!!
          transverseEdge1Id = newLineId;
          }
        nextFreeCellId++;
 
        // Fill in contents
        if ( configurationIsX12 )
          {
          newTriangle1->SetPointId( 0, idOfFirstPoint );
          newTriangle1->SetPointId( 1, point1Id );
          newTriangle1->SetPointId( 2, point4Id );

          newTriangle2->SetPointId( 0, idOfFirstPoint );
          newTriangle2->SetPointId( 1, point4Id );
          newTriangle2->SetPointId( 2, point2Id );
          
          newLine->SetPointId( 0, idOfFirstPoint );
          newLine->SetPointId( 1, point4Id );
          }
       else if ( configurationIs12X )
          {
          newTriangle1->SetPointId( 0, idOfThirdPoint );
          newTriangle1->SetPointId( 1, point1Id );
          newTriangle1->SetPointId( 2, point4Id );

          newTriangle2->SetPointId( 0, idOfThirdPoint );
          newTriangle2->SetPointId( 1, point4Id );
          newTriangle2->SetPointId( 2, point2Id );
          
          newLine->SetPointId( 0, idOfThirdPoint );
          newLine->SetPointId( 1, point4Id );
          }
        else if ( configurationIs2X1 )
          {
          newTriangle1->SetPointId( 0, idOfSecondPoint );
          newTriangle1->SetPointId( 1, point1Id );
          newTriangle1->SetPointId( 2, point4Id );

          newTriangle2->SetPointId( 0, idOfSecondPoint );
          newTriangle2->SetPointId( 1, point4Id );
          newTriangle2->SetPointId( 2, point2Id );
          
          newLine->SetPointId( 0, idOfSecondPoint );
          newLine->SetPointId( 1, point4Id );
          }
        else if ( configurationIsX21 )
          {
          newTriangle1->SetPointId( 0, idOfFirstPoint );
          newTriangle1->SetPointId( 1, point2Id );
          newTriangle1->SetPointId( 2, point4Id );

          newTriangle2->SetPointId( 0, idOfFirstPoint );
          newTriangle2->SetPointId( 1, point4Id );
          newTriangle2->SetPointId( 2, point1Id );
          
          newLine->SetPointId( 0, idOfFirstPoint );
          newLine->SetPointId( 1, point4Id );
          }
        else if ( configurationIs21X )
          {
          newTriangle1->SetPointId( 0, idOfThirdPoint );
          newTriangle1->SetPointId( 1, point2Id );
          newTriangle1->SetPointId( 2, point4Id );

          newTriangle2->SetPointId( 0, idOfThirdPoint );
          newTriangle2->SetPointId( 1, point4Id );
          newTriangle2->SetPointId( 2, point1Id );
          
          newLine->SetPointId( 0, idOfThirdPoint );
          newLine->SetPointId( 1, point4Id );
          }
        else if ( configurationIs1X2 )
          {
          newTriangle1->SetPointId( 0, idOfSecondPoint );
          newTriangle1->SetPointId( 1, point2Id );
          newTriangle1->SetPointId( 2, point4Id );

          newTriangle2->SetPointId( 0, idOfSecondPoint );
          newTriangle2->SetPointId( 1, point4Id );
          newTriangle2->SetPointId( 2, point1Id );
          
          newLine->SetPointId( 0, idOfSecondPoint );
          newLine->SetPointId( 1, point4Id );
          }
    

        // Add newly created cells
        splittedCells->InsertElement( newTriangle1Id, newTriangle1.ReleaseOwnership() );
        splittedCells->InsertElement( newTriangle2Id, newTriangle2.ReleaseOwnership() );
        splittedCells->InsertElement( newLineId, newLine.ReleaseOwnership() );
        }
      else
        {
        // Copy triangle cell
        AtlasMesh::CellAutoPointer  newTriangle;
        newTriangle.TakeOwnership( new TriangleCell );
        newTriangle->SetPointId( 0, idOfFirstPoint );
        newTriangle->SetPointId( 1, idOfSecondPoint );
        newTriangle->SetPointId( 2, idOfThirdPoint );
        splittedCells->InsertElement( cellIt.Index(), newTriangle.ReleaseOwnership() );
        }
      
      }
      
      

    } // End loop over all cells
    
    
    
  // Create a new mesh collection to hold the result.
  AtlasMeshCollection::Pointer  splitted = AtlasMeshCollection::New();
  splitted->SetPointParameters( splittedPointParameters );
  splitted->SetCells( splittedCells );
  splitted->SetReferencePosition( splittedReferencePosition );
  splitted->SetPositions( splittedPositions );
  splitted->SetK( m_K );

  return splitted;

#else
  return nullptr;
#endif

}

/*! 
  \fn void AtlasMeshCollection::FlattenAlphas()
  \brief Changes alphas at each point to have uniform alphas (ie, the
  prob of a given class at that point). Alphas at the edge of the FoV 
  have class0 (background) set to 1, all else 0.
*/

void
AtlasMeshCollection
::FlattenAlphas()
{
  // Create a flat alpha entry  
  const unsigned int  numberOfClasses = m_PointParameters->Begin().Value().m_Alphas.Size();
  AtlasAlphasType   flatAlphasEntry( numberOfClasses );
  flatAlphasEntry.Fill( 1.0f / static_cast< float >( numberOfClasses ) );

  // Create a border alphas entry (only first class is possible)
  AtlasAlphasType   borderAlphasEntry( numberOfClasses );
  borderAlphasEntry.Fill( 0.0f );
  borderAlphasEntry[ 0 ] = 1.0f;

  // Loop over all points, and replace the alphas with a flat alpha entry if the alphas can be changed,
  // or border alpha entry otherwise
  for ( PointDataContainerType::Iterator  pointParamIt = m_PointParameters->Begin();
        pointParamIt != m_PointParameters->End(); ++pointParamIt )
    {
    if ( pointParamIt.Value().m_CanChangeAlphas )
      {
      pointParamIt.Value().m_Alphas = flatAlphasEntry;
      }
    else
      {
      pointParamIt.Value().m_Alphas = borderAlphasEntry;
      }

    }
}

/*!
  \fn void AtlasMeshCollection::Transform( int meshNumber, const TransformType* transform )
  \brief Applies the given transform to the given mesh. The reference mesh is used when
  meshNumber = nMeshes+1. m_ReferenceTetrahedronInfos and m_Meshes are cleared forcing
  a regeneration of those data.
*/
void
AtlasMeshCollection
::Transform( int meshNumber, const TransformType* transform )
{
  // Sanity check on requested mesh
  if ( meshNumber >= static_cast< int >( m_Positions.size() ) )
    {
    itkExceptionMacro( "Can't transform mesh number " << meshNumber << " because there are only " << m_Positions.size() << " meshes!" );
    }

  if ( !transform )
    {
    std::cout << "No transform set to transform the mesh with!" << std::endl;
    return;
    }

  // Select the position container to work on
  PointsContainerType::Pointer  position;
  if ( meshNumber >= 0 )
    {
    position = m_Positions[ meshNumber ];
    }
  else
    {
    position = m_ReferencePosition;
    }

  // Now loop over all points, and alter their positions
  std::cout << "Transforming points" << std::endl;
  for ( PointsContainerType::Iterator  pointIt = position->Begin();
        pointIt != position->End();
        ++pointIt )
    {
    pointIt.Value() = transform->TransformPoint( pointIt.Value() );
    }


  // Cached data members are no longer valid
  m_ReferenceTetrahedronInfos = 0;
  m_Meshes.clear();

}



//
//
//
/*!
  \fn AtlasMesh::CellLinksContainerPointerAtlasMeshCollection::GetCellLinks() const
  \brief Returns m_CellLinks. Builds the cell links if they do not already exist.
  A "Cell Link" indicates which points belong to which cells (?).
*/
AtlasMesh::CellLinksContainerPointer
AtlasMeshCollection
::GetCellLinks() const
{
  if ( !m_CellLinks )
    {
    //std::cout << "Building cell links..." << std::endl;

    // Create links back from points to cells. 
    AtlasMesh::Pointer  tmpMesh = AtlasMesh::New();
    tmpMesh->SetCells( m_Cells );
  #if 1
    tmpMesh->SetPoints( m_Positions[ 0 ] );  // Bug in itk::Mesh::BuildCellLinks(): Won't do anything w/o points
  #endif
    tmpMesh->BuildCellLinks();
    m_CellLinks = tmpMesh->GetCellLinks();
    
    //std::cout << "... done!" << std::endl;
    }

  return m_CellLinks;
}


//
//
//
void
AtlasMeshCollection
::FillCubeWithTetrahedra( MeshSourceType* meshSource, bool flippedConfiguration,
                          const AtlasMesh::PointType&  p0,
                          const AtlasMesh::PointType&  p1,
                          const AtlasMesh::PointType&  p2,
                          const AtlasMesh::PointType&  p3,
                          const AtlasMesh::PointType&  p4,
                          const AtlasMesh::PointType&  p5,
                          const AtlasMesh::PointType&  p6,
                          const AtlasMesh::PointType&  p7 ) const
{
  // No matter the configuration, there are always five tethrahedra; the first point of
  // the first tehtrahedron is always p0, the fourth point of the first tethrahedron is p2, the
  // first point of the second tetrahedron is p7, the third point of the third tetrahedron
  // is p1, the fourth point of the third tetrahedron is p5, the third point of the fourth
  // tetrahedron is p6, the rourth point of the fourth tetrahedron is p4, and the third point
  // of the fifth tetrahedron is p3. These properties will prove useful when upsampling a
  // regular mesh collection.
  if ( flippedConfiguration )
    {
    meshSource->AddTetrahedron( meshSource->AddPoint( p0 ),
                                meshSource->AddPoint( p4 ),
                                meshSource->AddPoint( p1 ),
                                meshSource->AddPoint( p2 ) );
    meshSource->AddTetrahedron( meshSource->AddPoint( p7 ),
                                meshSource->AddPoint( p4 ),
                                meshSource->AddPoint( p2 ),
                                meshSource->AddPoint( p1 ) );
    meshSource->AddTetrahedron( meshSource->AddPoint( p7 ),
                                meshSource->AddPoint( p4 ),
                                meshSource->AddPoint( p1 ),
                                meshSource->AddPoint( p5 ) );
    meshSource->AddTetrahedron( meshSource->AddPoint( p2 ),
                                meshSource->AddPoint( p7 ),
                                meshSource->AddPoint( p6 ),
                                meshSource->AddPoint( p4 ) );
    meshSource->AddTetrahedron( meshSource->AddPoint( p2 ),
                                meshSource->AddPoint( p1 ),
                                meshSource->AddPoint( p3 ),
                                meshSource->AddPoint( p7 ) );

    }
  else
    {
    meshSource->AddTetrahedron( meshSource->AddPoint( p0 ),
                                meshSource->AddPoint( p6 ),
                                meshSource->AddPoint( p3 ),
                                meshSource->AddPoint( p2 ) );
    meshSource->AddTetrahedron( meshSource->AddPoint( p7 ),
                                meshSource->AddPoint( p5 ),
                                meshSource->AddPoint( p6 ),
                                meshSource->AddPoint( p3 ) );
    meshSource->AddTetrahedron( meshSource->AddPoint( p3 ),
                                meshSource->AddPoint( p0 ),
                                meshSource->AddPoint( p1 ),
                                meshSource->AddPoint( p5 ) );
    meshSource->AddTetrahedron( meshSource->AddPoint( p0 ),
                                meshSource->AddPoint( p5 ),
                                meshSource->AddPoint( p6 ),
                                meshSource->AddPoint( p4 ) );
    meshSource->AddTetrahedron( meshSource->AddPoint( p0 ),
                                meshSource->AddPoint( p5 ),
                                meshSource->AddPoint( p3 ),
                                meshSource->AddPoint( p6 ) );
    }



}




} // end namespace kvl
