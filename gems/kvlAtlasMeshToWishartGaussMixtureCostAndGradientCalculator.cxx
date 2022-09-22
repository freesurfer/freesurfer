#include "kvlAtlasMeshToWishartGaussMixtureCostAndGradientCalculator.h"

#include <itkMath.h>
#include "vnl/vnl_matrix_fixed.h"
#include "kvlTetrahedronInteriorConstIterator.h"


namespace kvl
{

//
//
//
AtlasMeshToWishartGaussMixtureCostAndGradientCalculator
::AtlasMeshToWishartGaussMixtureCostAndGradientCalculator()
{

    m_LikelihoodFilter = gmmLikelihoodFilterType::New();
    m_WMMLikelihoodFilter = wmmLikelihoodFilterType::New();

}


//
//
//
AtlasMeshToWishartGaussMixtureCostAndGradientCalculator
::~AtlasMeshToWishartGaussMixtureCostAndGradientCalculator()
{
}


//
//
//
void
AtlasMeshToWishartGaussMixtureCostAndGradientCalculator
::SetGaussianImages( const std::vector< ImageType::ConstPointer >& images )
{
    //
    for ( unsigned int contrastNumber = 0; contrastNumber < images.size(); contrastNumber++ )
    {
        dynamic_cast<gmmLikelihoodFilterType*>(m_LikelihoodFilter.GetPointer())
                ->SetInput( contrastNumber, images[ contrastNumber ]);
    }

    m_WMMLikelihoodFilter -> SetGaussianImages(images);
}

//
//
//
void
AtlasMeshToWishartGaussMixtureCostAndGradientCalculator
::SetWishartImages( const std::vector< ImageType::ConstPointer >& images )
{
    //
    for ( unsigned int contrastNumber = 0; contrastNumber < images.size(); contrastNumber++ )
    {
        m_WMMLikelihoodFilter->SetInput( contrastNumber, images[ contrastNumber ]);
    }

}


//
//
//
void
AtlasMeshToWishartGaussMixtureCostAndGradientCalculator
::Rasterize( const AtlasMesh* mesh )
{

  // Make sure the likelihoods are up-to-date
  //m_LikelihoodFilter->SetNumberOfThreads( 1 );
  m_WMMLikelihoodFilter->Update();

  // Now rasterize
  Superclass::Rasterize( mesh );

}



//
//
//
void
AtlasMeshToWishartGaussMixtureCostAndGradientCalculator
::SetParameters( const std::vector< vnl_vector< double > >& means,
                 const std::vector< vnl_matrix< double > >& variances,
                 const std::vector< double >&  mixtureWeights,
                 const std::vector< int >&  numberOfGaussiansPerClass,
                 const std::vector< double >& degreesOfFreedom,
                 const std::vector< vnl_matrix< double > >& scaleMatrices,
                 const std::vector< double >&  wmmMixtureWeights,
                 const std::vector< int >&  numberOfWishartsPerClass,
                 const double& voxratio)
{
    dynamic_cast<gmmLikelihoodFilterType*>(m_LikelihoodFilter.GetPointer())
            ->SetParameters( means,
                             variances,
                             mixtureWeights,
                             numberOfGaussiansPerClass);

    m_WMMLikelihoodFilter ->SetParameters( means,
                             variances,
                             mixtureWeights,
                             numberOfGaussiansPerClass,
                             degreesOfFreedom,
                             scaleMatrices,
                             wmmMixtureWeights,
                             numberOfWishartsPerClass,
                             voxratio);
}



//
//
//
void
AtlasMeshToWishartGaussMixtureCostAndGradientCalculator
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

  // Loop over all voxels within the tetrahedron and do The Right Thing
  const unsigned long  numberOfClasses = alphasInVertex0.Size();
  TetrahedronInteriorConstIterator< LikelihoodFilterType::OutputPixelType >  gmm_it( m_LikelihoodFilter->GetOutput(), p0, p1, p2, p3 );
  TetrahedronInteriorConstIterator< LikelihoodFilterType::OutputPixelType >  wmm_it( m_WMMLikelihoodFilter->GetOutput(), p0, p1, p2, p3 );
  for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
      gmm_it.AddExtraLoading( alphasInVertex0[ classNumber ],
                          alphasInVertex1[ classNumber ],
                          alphasInVertex2[ classNumber ],
                          alphasInVertex3[ classNumber ] );
      wmm_it.AddExtraLoading( alphasInVertex0[ classNumber ],
                          alphasInVertex1[ classNumber ],
                          alphasInVertex2[ classNumber ],
                          alphasInVertex3[ classNumber ] );
    }

  for ( ; !gmm_it.IsAtEnd(); ++gmm_it)
    {
    // Skip voxels for which nothing is known
    if ( gmm_it.Value().Size() == 0 | wmm_it.Value().Size() == 0)
      {
      //std::cout << "Skipping: " << it.Value().Size() << std::endl;
        ++wmm_it;
        continue;
      }

    //
    double likelihoodAgregate = 0.0;
    double  maxExponent;
    double  xGradientBasis = 0.0;
    double  yGradientBasis = 0.0;
    double  zGradientBasis = 0.0;
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      // Get the Gaussian mixture model likelihood of this class at the intensity of this pixel
      const double gmmmixturelog = log(gmm_it.Value()[ classNumber ]);
      const double wmmmixturelog = wmm_it.Value()[ classNumber ];
      const double mixturelog = gmmmixturelog+wmmmixturelog;

      const double weightInterpolated = gmm_it.GetExtraLoadingInterpolatedValue( classNumber );
      const double weightNextRow = gmm_it.GetExtraLoadingNextRowAddition( classNumber );
      const double weightNextColumn = gmm_it.GetExtraLoadingNextColumnAddition( classNumber );
      const double weightNextSlice = gmm_it.GetExtraLoadingNextSliceAddition( classNumber );

      // Add contribution of the likelihood
      if (classNumber == 0)
      {
          maxExponent = mixturelog;
          likelihoodAgregate = weightInterpolated;

          xGradientBasis = weightNextRow;
          yGradientBasis = weightNextColumn;
          zGradientBasis = weightNextSlice;
      }
      else if (mixturelog>maxExponent)
      {
          likelihoodAgregate = likelihoodAgregate*exp(maxExponent-mixturelog)
                  + weightInterpolated;

          xGradientBasis = xGradientBasis*exp(maxExponent-mixturelog)
                  + weightNextRow;
          yGradientBasis = yGradientBasis*exp(maxExponent-mixturelog)
                  + weightNextColumn;
          zGradientBasis = zGradientBasis*exp(maxExponent-mixturelog)
                  + weightNextSlice;

          maxExponent = mixturelog;
      }
      else
      {
          likelihoodAgregate += weightInterpolated*exp(mixturelog-maxExponent);

          xGradientBasis += weightNextRow*exp(mixturelog-maxExponent);
          yGradientBasis += weightNextColumn*exp(mixturelog-maxExponent);
          zGradientBasis += weightNextSlice*exp(mixturelog-maxExponent);
      }

      if (weightInterpolated!=wmm_it.GetExtraLoadingInterpolatedValue( classNumber ) |
              weightNextRow!=wmm_it.GetExtraLoadingNextRowAddition( classNumber )|
              weightNextColumn!=wmm_it.GetExtraLoadingNextColumnAddition( classNumber )|
              weightNextSlice!=wmm_it.GetExtraLoadingNextSliceAddition( classNumber ) )
      {
          std::cout<<"classNumber= "<<classNumber<<std::endl;
          std::cout<<"likelihoodAgregate = "<<likelihoodAgregate<<std::endl;
          std::cout<<"gmmmixturelog = "<<gmmmixturelog<<std::endl;
          std::cout<<"wmmmixturelog = "<<wmmmixturelog<<std::endl;
          std::cout<<"gmm interp = "<<gmm_it.GetExtraLoadingInterpolatedValue( classNumber )<<std::endl;
          std::cout<<"wmm interp = "<<wmm_it.GetExtraLoadingInterpolatedValue( classNumber )<<std::endl;
          std::cout<<"gmm row = "<<gmm_it.GetExtraLoadingNextRowAddition( classNumber )<<std::endl;
          std::cout<<"wmm row = "<<wmm_it.GetExtraLoadingNextRowAddition( classNumber )<<std::endl;
          std::cout<<"gmm column = "<<gmm_it.GetExtraLoadingNextColumnAddition( classNumber )<<std::endl;
          std::cout<<"wmm column = "<<wmm_it.GetExtraLoadingNextColumnAddition( classNumber )<<std::endl;
          std::cout<<"gmm slice = "<<gmm_it.GetExtraLoadingNextSliceAddition( classNumber )<<std::endl;
          std::cout<<"wmm slice = "<<wmm_it.GetExtraLoadingNextSliceAddition( classNumber )<<std::endl;
          itkExceptionMacro(<<"Voxel location mismatch");
      }
      } // End loop over all classes


    //  Add contribution to log-likelihood
    if (likelihoodAgregate<1e-15)
    {
        likelihoodAgregate = likelihoodAgregate + 1e-15; //dont want to divide by zero
    }
    priorPlusDataCost -= log( likelihoodAgregate ) + maxExponent;


    //
    xGradientBasis /= likelihoodAgregate;
    yGradientBasis /= likelihoodAgregate;
    zGradientBasis /= likelihoodAgregate;

    // Add contribution to gradient in vertex 0
    gradientInVertex0[ 0 ] += xGradientBasis * gmm_it.GetPi0();
    gradientInVertex0[ 1 ] += yGradientBasis * gmm_it.GetPi0();
    gradientInVertex0[ 2 ] += zGradientBasis * gmm_it.GetPi0();

    // Add contribution to gradient in vertex 1
    gradientInVertex1[ 0 ] += xGradientBasis * gmm_it.GetPi1();
    gradientInVertex1[ 1 ] += yGradientBasis * gmm_it.GetPi1();
    gradientInVertex1[ 2 ] += zGradientBasis * gmm_it.GetPi1();

    // Add contribution to gradient in vertex 2
    gradientInVertex2[ 0 ] += xGradientBasis * gmm_it.GetPi2();
    gradientInVertex2[ 1 ] += yGradientBasis * gmm_it.GetPi2();
    gradientInVertex2[ 2 ] += zGradientBasis * gmm_it.GetPi2();

    // Add contribution to gradient in vertex 3
    gradientInVertex3[ 0 ] += xGradientBasis * gmm_it.GetPi3();
    gradientInVertex3[ 1 ] += yGradientBasis * gmm_it.GetPi3();
    gradientInVertex3[ 2 ] += zGradientBasis * gmm_it.GetPi3();



    ++wmm_it;
    } // End loop over all voxels within the tetrahedron

//  std::cout<<"got out of tet loop"<<std::endl;
}





} // end namespace kvl

