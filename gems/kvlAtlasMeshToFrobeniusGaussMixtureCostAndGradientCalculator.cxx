#include "kvlAtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator.h"

#include <itkMath.h>
#include "vnl/vnl_matrix_fixed.h"
#include "kvlTetrahedronInteriorConstIterator.h"


namespace kvl
{

//
//
//
AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator
::AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator()
{

    m_LikelihoodFilter = gmmLikelihoodFilterType::New();
    m_DiffusionLikelihoodFilter = frobmmLikelihoodFilterType::New();

}


//
//
//
AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator
::~AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator()
{
}


#if 0
// Same functionality can be achieved using SetImages() from its base class.
//
//
//
void
AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator
::SetGaussianImages( const std::vector< ImageType::ConstPointer >& images )
{
    //
    for ( unsigned int contrastNumber = 0; contrastNumber < images.size(); contrastNumber++ )
    {
        // call SetInput() method from itk::ImageToImageFilter
        m_LikelihoodFilter->SetInput( contrastNumber, images[ contrastNumber ]);
    }

    // this call won't update m_DSWbetaMMLikelihoodFilter images because it is operating on a different object
    //m_FrobMMLikelihoodFilter -> SetGaussianImages(images);
}
#endif

#if 0
// moved to base class SetDiffusionImages()
//
//
void
AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator
::SetFrobeniusImages( const std::vector< ImageType::ConstPointer >& images )
{
    //
    for ( unsigned int contrastNumber = 0; contrastNumber < images.size(); contrastNumber++ )
    {
        m_FrobMMLikelihoodFilter->SetInput( contrastNumber, images[ contrastNumber ]);
    }

}


//
//
//
void
AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator
::Rasterize( const AtlasMesh* mesh )
{

  // Make sure the likelihoods are up-to-date
  //m_LikelihoodFilter->SetNumberOfThreads( 1 );
  m_FrobMMLikelihoodFilter->Update();

  // Now rasterize
  Superclass::Rasterize( mesh );

}
#endif



//
//
//
void
AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator
::SetParameters( const std::vector< vnl_vector< double > >& means,
                 const std::vector< vnl_matrix< double > >& variances,
                 const std::vector< double >&  mixtureWeights,
                 const std::vector< int >&  numberOfGaussiansPerClass )
{
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
AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator
::SetDiffusionParameters( const int numberOfContrasts,
                           const std::vector< double >&  frobMixtureWeights,
                           const std::vector< int >&  numberOfFrobeniusPerClass,
                           const double& voxratio,
                           const std::vector< double >& frobVariance,
                           const std::vector< vnl_vector< double > >& frobMeans )

{
    dynamic_cast<frobmmLikelihoodFilterType*>(m_DiffusionLikelihoodFilter.GetPointer())
            ->SetParameters( numberOfContrasts,
                             frobVariance,
                             frobMeans,
                             frobMixtureWeights,
                             numberOfFrobeniusPerClass,
                             voxratio);
}



//
//
//
void
AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator
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
  TetrahedronInteriorConstIterator< LikelihoodFilterType::OutputPixelType >  fmm_it( m_DiffusionLikelihoodFilter->GetOutput(), p0, p1, p2, p3 );
  for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
      gmm_it.AddExtraLoading( alphasInVertex0[ classNumber ],
                          alphasInVertex1[ classNumber ],
                          alphasInVertex2[ classNumber ],
                          alphasInVertex3[ classNumber ] );
      fmm_it.AddExtraLoading( alphasInVertex0[ classNumber ],
                          alphasInVertex1[ classNumber ],
                          alphasInVertex2[ classNumber ],
                          alphasInVertex3[ classNumber ] );
    }

  for ( ; !gmm_it.IsAtEnd(); ++gmm_it)
    {
    // Skip voxels for which nothing is known
    if ( gmm_it.Value().Size() == 0 | fmm_it.Value().Size() == 0)
      {
      //std::cout << "Skipping: " << it.Value().Size() << std::endl;
        ++fmm_it;
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
      const double fmmmixturelog = fmm_it.Value()[ classNumber ];
      const double mixturelog = gmmmixturelog+fmmmixturelog;

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

      if (weightInterpolated!=fmm_it.GetExtraLoadingInterpolatedValue( classNumber ) |
              weightNextRow!=fmm_it.GetExtraLoadingNextRowAddition( classNumber )|
              weightNextColumn!=fmm_it.GetExtraLoadingNextColumnAddition( classNumber )|
              weightNextSlice!=fmm_it.GetExtraLoadingNextSliceAddition( classNumber ) )
      {
          std::cout<<"classNumber= "<<classNumber<<std::endl;
          std::cout<<"likelihoodAgregate = "<<likelihoodAgregate<<std::endl;
          std::cout<<"gmmmixturelog = "<<gmmmixturelog<<std::endl;
          std::cout<<"wmmmixturelog = "<<fmmmixturelog<<std::endl;
          std::cout<<"gmm interp = "<<gmm_it.GetExtraLoadingInterpolatedValue( classNumber )<<std::endl;
          std::cout<<"wmm interp = "<<fmm_it.GetExtraLoadingInterpolatedValue( classNumber )<<std::endl;
          std::cout<<"gmm row = "<<gmm_it.GetExtraLoadingNextRowAddition( classNumber )<<std::endl;
          std::cout<<"wmm row = "<<fmm_it.GetExtraLoadingNextRowAddition( classNumber )<<std::endl;
          std::cout<<"gmm column = "<<gmm_it.GetExtraLoadingNextColumnAddition( classNumber )<<std::endl;
          std::cout<<"wmm column = "<<fmm_it.GetExtraLoadingNextColumnAddition( classNumber )<<std::endl;
          std::cout<<"gmm slice = "<<gmm_it.GetExtraLoadingNextSliceAddition( classNumber )<<std::endl;
          std::cout<<"wmm slice = "<<fmm_it.GetExtraLoadingNextSliceAddition( classNumber )<<std::endl;
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



    ++fmm_it;
    } // End loop over all voxels within the tetrahedron

//  std::cout<<"got out of tet loop"<<std::endl;
}





} // end namespace kvl

