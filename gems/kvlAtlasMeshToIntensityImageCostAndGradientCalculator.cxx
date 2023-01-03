#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"

#include <itkMath.h>
#include "vnl/vnl_matrix_fixed.h"
#include "kvlTetrahedronInteriorConstIterator.h"



namespace kvl
{

//
//
//
AtlasMeshToIntensityImageCostAndGradientCalculator
::AtlasMeshToIntensityImageCostAndGradientCalculator()
{
#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
  m_tetrahedronCnt = 0;
  m_zeroVoxel_tetrahedronCnt = 0;

  m_totalVoxel = 0;
  m_totalVoxelInTetrahedron = 0;
#endif
  /* 
   * 1. m_LikelihoodFilter is declared as LikelihoodImageFilterBase
   * 2. m_LikelihoodFilter is holding a GMMLikelihoodImageFilter object
   */ 

  m_LikelihoodFilter = gmmLikelihoodFilterType::New();

}


//
//
//
AtlasMeshToIntensityImageCostAndGradientCalculator
::~AtlasMeshToIntensityImageCostAndGradientCalculator()
{
}


//
//
//
void 
AtlasMeshToIntensityImageCostAndGradientCalculator
::SetParameters( const std::vector< vnl_vector< double > >& means,
                 const std::vector< vnl_matrix< double > >& variances,
                 const std::vector< double >&  mixtureWeights,
                 const std::vector< int >&  numberOfGaussiansPerClass )
{
   /* 
    * 1. m_LikelihoodFilter is declared as LikelihoodImageFilterBase
    * 2. m_LikelihoodFilter is holding a GMMLikelihoodImageFilter object
    * 3. dynamic_case m_LikelihoodFilter to GMMLikelihoodImageFilter, so it can call derived class method SetParameters()
    */ 
    dynamic_cast<gmmLikelihoodFilterType*>(m_LikelihoodFilter.GetPointer())
            ->SetParameters( means,
                             variances,
                             mixtureWeights,
                             numberOfGaussiansPerClass );
}


//
//
//
void 
AtlasMeshToIntensityImageCostAndGradientCalculator
::AddDataContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                    const AtlasMesh::PointType& p1,
                                    const AtlasMesh::PointType& p2,
                                    const AtlasMesh::PointType& p3,
                                    const AtlasAlphasType&  alphasInVertex0,
                                    const AtlasAlphasType&  alphasInVertex1,
                                    const AtlasAlphasType&  alphasInVertex2,
                                    const AtlasAlphasType&  alphasInVertex3,
                                    ThreadAccumDataType&  priorPlusDataCost,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex0,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex1,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex2,
                                    AtlasPositionGradientThreadAccumType&  gradientInVertex3 )
{
  
#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
  int voxelCnt = 0;
  m_tetrahedronCnt++;
#endif
  // Loop over all voxels within the tetrahedron and do The Right Thing  
  const int  numberOfClasses = alphasInVertex0.Size();
  TetrahedronInteriorConstIterator< LikelihoodFilterType::OutputPixelType >  it( m_LikelihoodFilter->GetOutput(), p0, p1, p2, p3 );
  for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
      if (alphasInVertex0[ classNumber ] != 0 || alphasInVertex1[ classNumber ] != 0 || alphasInVertex2[ classNumber ] != 0 || alphasInVertex3[ classNumber ] != 0)
        it.AddExtraLoading( alphasInVertex0[ classNumber ], 
                            alphasInVertex1[ classNumber ], 
                            alphasInVertex2[ classNumber ], 
                            alphasInVertex3[ classNumber ] );
    }  
    
  for ( ; !it.IsAtEnd(); ++it )
    {
#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
    voxelCnt++;
    it.m_totalVoxelInTetrahedron++;
#endif

    // Skip voxels for which nothing is known
    if ( it.Value().Size() == 0 )
      {
      //std::cout << "Skipping: " << it.Value().Size() << std::endl;
      continue;
      }
      
    //
    double likelihood = 0.0;
    double  xGradientBasis = 0.0;
    double  yGradientBasis = 0.0;
    double  zGradientBasis = 0.0;
    int classIdx = 0; 
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      if (alphasInVertex0[ classNumber ] != 0 || alphasInVertex1[ classNumber ] != 0 || alphasInVertex2[ classNumber ] != 0 || alphasInVertex3[ classNumber ] != 0)
	{
        // Get the Gaussian mixture model likelihood of this class at the intensity of this pixel
        const double mixture = it.Value()[ classNumber ];
        
        // Add contribution of the likelihood
        likelihood += mixture * it.GetExtraLoadingInterpolatedValue( classIdx );
      
        //
        xGradientBasis += mixture * it.GetExtraLoadingNextRowAddition( classIdx );
        yGradientBasis += mixture * it.GetExtraLoadingNextColumnAddition( classIdx );
        zGradientBasis += mixture * it.GetExtraLoadingNextSliceAddition( classIdx );

        classIdx++;
	}
      } // End loop over all classes
      
      
    //  Add contribution to log-likelihood
    likelihood = likelihood + 1e-15; //dont want to divide by zero
    priorPlusDataCost -= log( likelihood );


    //
    xGradientBasis /= likelihood;
    yGradientBasis /= likelihood;
    zGradientBasis /= likelihood;

    // Add contribution to gradient in vertex 0
    gradientInVertex0[ 0 ] += xGradientBasis * it.GetPi0();
    gradientInVertex0[ 1 ] += yGradientBasis * it.GetPi0();
    gradientInVertex0[ 2 ] += zGradientBasis * it.GetPi0();

    // Add contribution to gradient in vertex 1
    gradientInVertex1[ 0 ] += xGradientBasis * it.GetPi1();
    gradientInVertex1[ 1 ] += yGradientBasis * it.GetPi1();
    gradientInVertex1[ 2 ] += zGradientBasis * it.GetPi1();
    
    // Add contribution to gradient in vertex 2
    gradientInVertex2[ 0 ] += xGradientBasis * it.GetPi2();
    gradientInVertex2[ 1 ] += yGradientBasis * it.GetPi2();
    gradientInVertex2[ 2 ] += zGradientBasis * it.GetPi2();
    
    // Add contribution to gradient in vertex 3
    gradientInVertex3[ 0 ] += xGradientBasis * it.GetPi3();
    gradientInVertex3[ 1 ] += yGradientBasis * it.GetPi3();
    gradientInVertex3[ 2 ] += zGradientBasis * it.GetPi3();
    
    
    } // End loop over all voxels within the tetrahedron


#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
  if (voxelCnt == 0)
    m_zeroVoxel_tetrahedronCnt++;

  printf("[DEBUG] Iter #%2d, Rasterize #%2d, Tetrahedron #%8d: totalVoxel = %4d, InTetrahedron = %4d, %6.2f\%\n",
         m_Iterations, m_RasterizeCalls,
         m_tetrahedronCnt, it.m_totalVoxel, it.m_totalVoxelInTetrahedron, 
         (it.m_totalVoxel == 0) ? 0.0 : (double)it.m_totalVoxelInTetrahedron/it.m_totalVoxel*100);

  if (it.m_totalVoxel < it.m_totalVoxelInTetrahedron)
    printf("[DEBUG] ###### In Tetrahedron Voxel Count Error: Iter #%2d, Rasterize #%2d, Tetrahedron #%8d: totalVoxel = %4d, InTetrahedron = %4d, %6.2f\%\n",
           m_Iterations, m_RasterizeCalls,
           m_tetrahedronCnt, it.m_totalVoxel, it.m_totalVoxelInTetrahedron, 
           (it.m_totalVoxel == 0) ? 0.0 : (double)it.m_totalVoxelInTetrahedron/it.m_totalVoxel*100);

  if (voxelCnt == 0 && !(it.m_totalVoxel == 0 || it.m_totalVoxelInTetrahedron == 0))
  {
    printf("[DEBUG] ###### Voxel Zero Count Error: Iter #%2d, Rasterize #%2d, Tetrahedron #%8d: totalVoxel = %4d, InTetrahedron = %4d\n",
           m_Iterations, m_RasterizeCalls,
           m_tetrahedronCnt, it.m_totalVoxel, it.m_totalVoxelInTetrahedron); 
  }

  m_totalVoxel += it.m_totalVoxel;
  m_totalVoxelInTetrahedron += it.m_totalVoxelInTetrahedron;
  it.m_totalVoxel = 0;
  it.m_totalVoxelInTetrahedron = 0;
#endif  
}

} // end namespace kvl
