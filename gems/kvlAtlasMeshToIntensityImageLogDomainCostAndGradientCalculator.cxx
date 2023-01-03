#include "kvlAtlasMeshToIntensityImageLogDomainCostAndGradientCalculator.h"

#include <itkMath.h>
#include "vnl/vnl_matrix_fixed.h"
#include "kvlTetrahedronInteriorConstIterator.h"



namespace kvl
{

//
//
//
AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator
::AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator()
{
  /* 
   * 1. m_LikelihoodFilter is declared as LikelihoodImageFilterBase
   * 2. m_LikelihoodFilter is holding a GMMLikelihoodImageFilter object
   */ 

  m_LikelihoodFilter = gmmLikelihoodFilterType::New();

}


//
//
//
AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator
::~AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator()
{
}


#if 0
//
//
//
void 
AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator
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
#endif

//
//
//
void 
AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator
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
  const int numberOfClasses = alphasInVertex0.Size();
  TetrahedronInteriorConstIterator< LikelihoodFilterType::OutputPixelType >  it( m_LikelihoodFilter->GetOutput(), p0, p1, p2, p3 );
  for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
    it.AddExtraLoading( alphasInVertex0[ classNumber ],
                        alphasInVertex1[ classNumber ],
                        alphasInVertex2[ classNumber ],
                        alphasInVertex3[ classNumber ] );
    }

  for ( ; !it.IsAtEnd(); ++it )
    {
    // Skip voxels for which nothing is known
    if ( it.Value().Size() == 0 )
      {
      //std::cout << "Skipping: " << it.Value().Size() << std::endl;
      continue;
      }

    //
    double likelihoodAgregate = 0.0;
    double  maxExponent = -1000.0;
    double  xGradientBasis = 0.0;
    double  yGradientBasis = 0.0;
    double  zGradientBasis = 0.0;
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      // Get the Gaussian mixture model likelihood of this class at the intensity of this pixel
      const double mixture = it.Value()[ classNumber ];

      if (mixture>0)
      {
          const double mixturelog = log(mixture);


          const double weightInterpolated = it.GetExtraLoadingInterpolatedValue( classNumber );
          const double weightNextRow = it.GetExtraLoadingNextRowAddition( classNumber );
          const double weightNextColumn = it.GetExtraLoadingNextColumnAddition( classNumber );
          const double weightNextSlice = it.GetExtraLoadingNextSliceAddition( classNumber );

          // Add contribution of the likelihood
          if (maxExponent==-1000)
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
          // Add contribution of the likelihood
          //likelihood += mixture * it.GetExtraLoadingInterpolatedValue( classNumber );

          //
          //xGradientBasis += mixture * it.GetExtraLoadingNextRowAddition( classNumber );
          //yGradientBasis += mixture * it.GetExtraLoadingNextColumnAddition( classNumber );
          //zGradientBasis += mixture * it.GetExtraLoadingNextSliceAddition( classNumber );

          //if (abs( likelihood-(likelihoodAgregate*exp(maxExponent)))>(abs(likelihood)*1e-15))
          //{
          //std::cout<<"likelihood = "<<likelihood<<std::endl;
          //std::cout<<"likelihoodAgregate = "<<likelihood<<std::endl;
          //std::cout<<"mixturelog = "<<mixturelog<<std::endl;
          //std::cout<<"mixture = "<<likelihood<<std::endl;
          //itkExceptionMacro(<<"Voxel location mismatch");
          //}
      }
      } // End loop over all classes

    //  Add contribution to log-likelihood
    //likelihood = likelihood + 1e-15; //dont want to divide by zero
    //priorPlusDataCost -= log( likelihood );


    //
    //xGradientBasis /= likelihood;
    //yGradientBasis /= likelihood;
    //zGradientBasis /= likelihood;


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

}

} // end namespace kvl
