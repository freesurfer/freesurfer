#ifndef __kvlAtlasMeshToIntensityImageGradientCalculator_h
#define __kvlAtlasMeshToIntensityImageGradientCalculator_h

#include "kvlAtlasMeshRasterizor.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkAffineTransform.h"

#define PI 3.141592654

#define ENABLE_POISSON_RATHER_THAN_GAUSSIAN false

namespace kvl
{


namespace FragmentProcessor 
{

/**
 *
 */
class CalculateGradientToIntensityImage
{
public:
  
  typedef itk::Image< float, 3 >  ImageType;
  typedef itk::Image< AtlasAlphasType, 3 >  ProbabilityImageType;
  typedef itk::Image< unsigned char, 3 >  SegmentedImageType;

  
  CalculateGradientToIntensityImage() 
    {
    //m_Image = 0;
    m_ProbabilityImage = 0;
    m_SegmentedImage = 0;
      
    m_Mesh = 0;
    m_GradientInVertex0 = 0;
    m_GradientInVertex1 = 0;
    m_GradientInVertex2 = 0;
    m_GradientInVertex3 = 0;
    m_MinLogLikelihoodTimesPrior = 0;

    m_IgnoreDeformationPrior = false;
    m_OnlyDeformationPrior = false;

    m_mapCompToComp = 0; 

    if(ENABLE_POISSON_RATHER_THAN_GAUSSIAN)
       {
       cumsumLogX[0]=0;
       cumsumLogX[1]=0;
       for(int x=2; x<10000; x++)
         {
         cumsumLogX[x]=cumsumLogX[x-1]+log(x);
         }
       }
     }

  ~CalculateGradientToIntensityImage() {};


  void SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
    { m_mapCompToComp = mapCompToComp; }
  std::vector<unsigned char > * GetMapCompToComp()
    { return m_mapCompToComp; }

  void SetImages( std::vector<itk::Image< float, 3 >::Pointer> images ) 
    {
    m_Images = images;
    m_n_channels = images.size();
    }

  const itk::Image< float, 3 >::Pointer GetImage(int ind) const
    { 
    return m_Images[ind]; 
    }

  void SetProbabilityImage( const ProbabilityImageType* probabilityImage ) 
    {
    m_ProbabilityImage = probabilityImage;
    }

  const ProbabilityImageType* GetProbabilityImage() const
    { return m_ProbabilityImage; }
 
  void SetSegmentedImage( const SegmentedImageType* segmentedImage ) 
    {
    m_SegmentedImage = segmentedImage;
    }

  const SegmentedImageType* GetSegmentedImage() const
    { return m_SegmentedImage; }
  
 
  inline void operator()( const float& pi0, const float& pi1, const float& pi2, const float& pi3 )
    {
    // This will hold the precious gradient basis from which this voxel's contribution to each
    // of the tetrahedral nodes can be easily calculated
    double  tmpX = 0.0f;
    double  tmpY = 0.0f;
    double  tmpZ = 0.0f;

    if(m_OnlyDeformationPrior)
    {
      // Move on to the next pixel
      m_Index[ 0 ]++;
      return;
    }
    else if ( m_SegmentedImage )
      {

      if ( m_mapCompToComp == 0) // no collapsed labels
        {
        const SegmentedImageType::PixelType  classNumber = m_SegmentedImage->GetPixel( m_Index );
  
        // Collect the data terms for each vertex
        double  alpha0 = m_AlphasInVertex0[ classNumber ];
        double  alpha1 = m_AlphasInVertex1[ classNumber ];
        double  alpha2 = m_AlphasInVertex2[ classNumber ];
        double  alpha3 = m_AlphasInVertex3[ classNumber ];
  
        // Add contribution of the likelihood
        double  likelihood = alpha0 * pi0 + alpha1 * pi1 + alpha2 * pi2 + alpha3 * pi3 + 1e-15;
        double cost = -log( likelihood );
  
        //
        if (std::isnan(cost))  // Eugenio added this check
          {
             m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
            tmpX = 0; tmpY=0; tmpZ=0;
          }
        else
          {
            const double  tmp1 = 1.0 / likelihood;
            tmpX += tmp1 * ( m_XGradientBasis[ classNumber ] );
            tmpY += tmp1 * ( m_YGradientBasis[ classNumber ] );
            tmpZ += tmp1 * ( m_ZGradientBasis[ classNumber ] );
            m_MinLogLikelihoodTimesPrior += cost;
          }


        }
      else // version with collapsed labels
        {
        const SegmentedImageType::PixelType  collapsedClassNumber = m_SegmentedImage->GetPixel( m_Index );

        double  likelihood =  1e-15;
        for(int ind=0; ind<m_mapCompToComp[collapsedClassNumber].size(); ind++)
          {
          const SegmentedImageType::PixelType k = m_mapCompToComp[collapsedClassNumber][ind];
          // // std::cout  << " katigc " << ((int) collapsedClassNumber) << " " << ((int) k) << std::endl;
          double  alpha0 = m_AlphasInVertex0[ k ];
          double  alpha1 = m_AlphasInVertex1[ k ];
          double  alpha2 = m_AlphasInVertex2[ k ];
          double  alpha3 = m_AlphasInVertex3[ k ];

          likelihood += alpha0 * pi0 + alpha1 * pi1 + alpha2 * pi2 + alpha3 * pi3;

          tmpX += ( m_XGradientBasis[ k ] ); 
          tmpY += ( m_YGradientBasis[ k ] );
          tmpZ += ( m_ZGradientBasis[ k ] );

          }
        double cost = -log( likelihood );
  
        if (std::isnan(cost)) // Eugenio added this check
          {
            m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
            tmpX = 0; tmpY=0; tmpZ=0;
          }
        else
          {
            m_MinLogLikelihoodTimesPrior += cost;
            const double  tmp1 = 1.0 / likelihood;
            tmpX = tmpX * tmp1; 
            tmpY = tmpY * tmp1; 
            tmpZ = tmpZ * tmp1; 
          }
        }



      }
    else if ( m_ProbabilityImage )
      {
      // Loop over all classes
      const AtlasAlphasType&  weights = m_ProbabilityImage->GetPixel( m_Index );
      double  cost = 0.0f;
 
      if (isnan(weights[0]))
      {
        m_Index[ 0 ]++;
        return;
      }

      for ( unsigned int classNumber = 0; classNumber < m_AlphasInVertex0.Size(); classNumber++ )
        {
        // Get the weight of this class's contribution
        const double  weight = weights[ classNumber ];

        // Collect the data terms for each vertex
        double  alpha0 = m_AlphasInVertex0[ classNumber ];
        double  alpha1 = m_AlphasInVertex1[ classNumber ];
        double  alpha2 = m_AlphasInVertex2[ classNumber ];
        double  alpha3 = m_AlphasInVertex3[ classNumber ];

        // Add contribution of the likelihood
        double  likelihood = alpha0 * pi0 + alpha1 * pi1 + alpha2 * pi2 + alpha3 * pi3 + 1e-15;
        cost -= weight * log( likelihood );

        //
        const double  tmp1 = weight / likelihood;
        tmpX += tmp1 * ( m_XGradientBasis[ classNumber ] );
        tmpY += tmp1 * ( m_YGradientBasis[ classNumber ] );
        tmpZ += tmp1 * ( m_ZGradientBasis[ classNumber ] );

        }
      m_MinLogLikelihoodTimesPrior += cost;
      }
    else if (m_AlphasInVertex0[0]<-0.5) // special case: hack to handle partial voluming with a single Gaussian per vertex, with different means / variances
    {
      for(int c=0; c<m_n_channels; c++)
      {
     
        ImageType::PixelType intensity=m_Images[c]->GetPixel( m_Index );
        if (intensity > 0){

          double mu0 =  m_AlphasInVertex0[ 2*c+1 ];
          double mu1 =  m_AlphasInVertex1[ 2*c+1 ];
          double mu2 =  m_AlphasInVertex2[ 2*c+1 ];
          double mu3 =  m_AlphasInVertex3[ 2*c+1 ];

          double var0 =  m_AlphasInVertex0[ 2*c+2 ];
          double var1 =  m_AlphasInVertex1[ 2*c+2 ];
          double var2 =  m_AlphasInVertex2[ 2*c+2 ];
          double var3 =  m_AlphasInVertex3[ 2*c+2 ];

          double mu =  mu0 * pi0 + mu1 * pi1 + mu2 * pi2 + mu3 * pi3;
          double var =  var0 * pi0 + var1 * pi1 + var2 * pi2 + var3 * pi3;

          double likelihood = 1e-15 + exp( -0.5 / var * pow( intensity - mu , 2 )) / sqrt( 2 * PI * var );
          m_MinLogLikelihoodTimesPrior -= log( likelihood );

          double  intensityXGradientBasis = 0.0;
          double  intensityYGradientBasis = 0.0;
          double  intensityZGradientBasis = 0.0;

          //
          double tmp = 2*(intensity - mu) * var;
          intensityXGradientBasis += tmp * m_XGradientBasis[ 1 ];
          intensityYGradientBasis += tmp * m_YGradientBasis[ 1 ];
          intensityZGradientBasis += tmp * m_ZGradientBasis[ 1 ];

          tmp = (intensity - mu) * (intensity - mu) - var;
          intensityXGradientBasis += tmp * m_XGradientBasis[ 2 ];
          intensityYGradientBasis += tmp * m_YGradientBasis[ 2 ];
          intensityZGradientBasis += tmp * m_ZGradientBasis[ 2 ];
    

          tmpX += intensityXGradientBasis / var / var;
          tmpY += intensityYGradientBasis / var / var;
          tmpZ += intensityZGradientBasis / var / var;
        }
      }

    }
    else
      {
      //// std::cout << "!!!!!! Using original image to directly deform to !!!!!" << std::endl;
      vnl_vector<ImageType::PixelType> intensity_v; 
      std::vector<bool> isThere(m_n_channels);
      ImageType::PixelType intensity=0.0;
      // Eugenio: these 2 are only useful for multivariate Gaussian, but whatever...
      int nPresent=0;
      int index=0;

      if(m_n_channels > 1)
      {
       	//std::cout<<"More than one channel!"<<std::endl;
        int aux=1;
        vnl_vector<ImageType::PixelType> aux_v(m_n_channels,0.0);
	for(int c=0; c<m_n_channels; c++)
        {
          ImageType::PixelType p = m_Images[c]->GetPixel( m_Index );
          if(p>0){
            isThere[c]=true;
            aux_v[nPresent]=p;
            nPresent++;
            index+=aux;
          }
          else
          {
            isThere[c]=false;
          }
          aux = aux << 1;
        }
        if ( nPresent==0 )
        {
          // Move on to the next pixel
          m_Index[ 0 ]++;
          return;
        }
        else
        {
          intensity_v=aux_v.extract(nPresent);
        }
	
      }
      else
      {
        intensity=m_Images[0]->GetPixel( m_Index );
        if ( intensity == 0 )
        {
          // Move on to the next pixel
          m_Index[ 0 ]++;
          return;
        }
      }

      double likelihood = 0.0;
      double  intensityXGradientBasis = 0.0;
      double  intensityYGradientBasis = 0.0;
      double  intensityZGradientBasis = 0.0;
      for ( unsigned int classNumber = 0; classNumber < m_AlphasInVertex0.Size(); classNumber++ )
        {
        // Evaluate the Gaussian of this class at the intensity of this pixel
        double gauss = 0.0;
        if(m_n_channels == 1)
        {

        if(ENABLE_POISSON_RATHER_THAN_GAUSSIAN)
          {
          double csl=0;
          int x1=floor(intensity);
          if(intensity==x1)
            {
            csl=cumsumLogX[x1];
            }
          else
            {
            int x2=ceil(intensity);
            double a=intensity-x1;
            csl=(1-a)*cumsumLogX[x1]+a*cumsumLogX[x2];
            }
          gauss = exp (intensity * log( m_Means[ classNumber ][0] ) -  m_Means[ classNumber ][0] - csl );
          }
        else
          {
          gauss = exp( -0.5 * m_Precisions[classNumber][1][0][0] * pow( intensity - m_Means[ classNumber ][0] , 2 )) / 
                                      sqrt( 2 * PI / m_Precisions[ classNumber ][1][0][0] ); //NOTE: These are precisions now (precision = 1/variance) so we multiply
          }

	}
        else
        {
          
          vnl_vector<float>  dataV(intensity_v.size());
          int c=0;
          for(int i=0; i<m_n_channels; i++)
          {
            if(isThere[i])
            {
              dataV[c]=intensity_v[c]-m_Means[classNumber][i];
              c++;
            }
          }

          gauss= exp( -0.5 * dot_product(dataV,m_Precisions[classNumber][index]*dataV)) * OneOverSqrtDetCov[classNumber][index] * m_piTermMultiv[nPresent];
 
        }

        // Collect the data terms for each vertex
        double  alpha0 = m_AlphasInVertex0[ classNumber ];
        double  alpha1 = m_AlphasInVertex1[ classNumber ];
        double  alpha2 = m_AlphasInVertex2[ classNumber ];
        double  alpha3 = m_AlphasInVertex3[ classNumber ];

        // Add contribution of the likelihood
        likelihood += gauss * ( alpha0 * pi0 + alpha1 * pi1 + alpha2 * pi2 + alpha3 * pi3 );
	

        //
        intensityXGradientBasis += gauss * m_XGradientBasis[ classNumber ];
        intensityYGradientBasis += gauss * m_YGradientBasis[ classNumber ];
        intensityZGradientBasis += gauss * m_ZGradientBasis[ classNumber ];
        }
      likelihood = likelihood + 1e-15; //dont want to divide by zero
      m_MinLogLikelihoodTimesPrior -= log( likelihood );


      //
      tmpX = intensityXGradientBasis / likelihood;
      tmpY = intensityYGradientBasis / likelihood;
      tmpZ = intensityZGradientBasis / likelihood;

      if(0){//(m_Index[0]%50 == 0) && (m_Index[1]%50 == 0) && (m_Index[2]%50 == 0)){
	// std::cout<<"Gradients: ("<<tmpX<<","<<tmpY<<","<<tmpZ<<")"<<std::endl;
	// std::cout<<"Likelihood: "<<likelihood<<std::endl;
      }


      }

     
    // Add contribution to gradient in vertex 0
    AtlasPositionGradientType  gradientContributionToVertex0;
    gradientContributionToVertex0[ 0 ] = tmpX * pi0;
    gradientContributionToVertex0[ 1 ] = tmpY * pi0;
    gradientContributionToVertex0[ 2 ] = tmpZ * pi0;
    *m_GradientInVertex0 += gradientContributionToVertex0;

    // Add contribution to gradient in vertex 1
    AtlasPositionGradientType  gradientContributionToVertex1;
    gradientContributionToVertex1[ 0 ] = tmpX * pi1;
    gradientContributionToVertex1[ 1 ] = tmpY * pi1;
    gradientContributionToVertex1[ 2 ] = tmpZ * pi1;
    *m_GradientInVertex1 += gradientContributionToVertex1;

    // Add contribution to gradient in vertex 2
    AtlasPositionGradientType  gradientContributionToVertex2;
    gradientContributionToVertex2[ 0 ] = tmpX * pi2;
    gradientContributionToVertex2[ 1 ] = tmpY * pi2;
    gradientContributionToVertex2[ 2 ] = tmpZ * pi2;
    *m_GradientInVertex2 += gradientContributionToVertex2;

    // Add contribution to gradient in vertex 3
    AtlasPositionGradientType  gradientContributionToVertex3;
    gradientContributionToVertex3[ 0 ] = tmpX * pi3;
    gradientContributionToVertex3[ 1 ] = tmpY * pi3;
    gradientContributionToVertex3[ 2 ] = tmpZ * pi3;
    *m_GradientInVertex3 += gradientContributionToVertex3;


    // Move on to the next pixel
    m_Index[ 0 ]++;
    }
    
  inline void StartNewSpan( int x, int y, int z, const unsigned char* sourcePointer )
    {
    m_Index[ 0 ] = x;
    m_Index[ 1 ] = y;
    m_Index[ 2 ] = z;
    }
    
  inline bool StartNewTetrahedron( AtlasMesh::CellIdentifier cellId )
    {
    //
    // Cache relevant elements of the vertices of this triangle. The notation used is Y = [ p0 p1 p2 p3; 1 1 1 1 ]
    //
    AtlasMesh::CellAutoPointer  cell;
    m_Mesh->GetCell( cellId, cell );
          
    AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
    AtlasMesh::PointType  p;
    m_Mesh->GetPoint( *pit, &p );
    const double y11 = p[ 0 ];
    const double y21 = p[ 1 ];
    const double y31 = p[ 2 ];
    m_AlphasInVertex0 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_GradientInVertex0 = &( m_PositionGradient->ElementAt( *pit ) );
 
    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    const double y12 = p[ 0 ];
    const double y22 = p[ 1 ];
    const double y32 = p[ 2 ];
    m_AlphasInVertex1 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_GradientInVertex1 = &( m_PositionGradient->ElementAt( *pit ) );

    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    const double y13 = p[ 0 ];
    const double y23 = p[ 1 ];
    const double y33 = p[ 2 ];
    m_AlphasInVertex2 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_GradientInVertex2 = &( m_PositionGradient->ElementAt( *pit ) );
  
    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    const double y14 = p[ 0 ];
    const double y24 = p[ 1 ];
    const double y34 = p[ 2 ];
    m_AlphasInVertex3 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_GradientInVertex3 = &( m_PositionGradient->ElementAt( *pit ) );
  
    if ( !m_IgnoreDeformationPrior )
      {
      //
      // Retrieve reference triangle info components. Z is inv( [ p0 p1 p2 p3; 1 1 1 1 ] ) of the
      // tetrahedron in reference position
      //
      ReferenceTetrahedronInfo  info;
      m_Mesh->GetCellData( cellId, &info );

      const double  referenceVolumeTimesK = info.m_ReferenceVolumeTimesK;

      const double  z11 = info.m_Z11;
      const double  z21 = info.m_Z21;
      const double  z31 = info.m_Z31;
      const double  z41 = info.m_Z41;

      const double  z12 = info.m_Z12;
      const double  z22 = info.m_Z22;
      const double  z32 = info.m_Z32;
      const double  z42 = info.m_Z42;

      const double  z13 = info.m_Z13;
      const double  z23 = info.m_Z23;
      const double  z33 = info.m_Z33;
      const double  z43 = info.m_Z43;


      //
      // Now let's add Ashburner's prior cost for the tethrahedron deformation from its reference position
      //
      const double  m11 = z11*y11 + z21*y12 + z31*y13 + z41*y14;
      const double  m21 = z11*y21 + z21*y22 + z31*y23 + z41*y24;
      const double  m31 = z11*y31 + z21*y32 + z31*y33 + z41*y34;
      const double  m12 = z12*y11 + z22*y12 + z32*y13 + z42*y14;
      const double  m22 = z12*y21 + z22*y22 + z32*y23 + z42*y24;
      const double  m32 = z12*y31 + z22*y32 + z32*y33 + z42*y34;
      const double  m13 = z13*y11 + z23*y12 + z33*y13 + z43*y14;
      const double  m23 = z13*y21 + z23*y22 + z33*y23 + z43*y24;
      const double  m33 = z13*y31 + z23*y32 + z33*y33 + z43*y34;

      const double  detJ = m11 * ( m22*m33 - m32*m23 ) - m12 * ( m21*m33 - m31*m23 ) + m13 * ( m21*m32 - m31*m22 );
      if ( detJ <= 0 )
        {
        // std::cout << "Oooooops: tetrahedron " << cellId << " has managed to turn bad (detJ: " << detJ << ")" << std::endl;
        m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
        // std::cout << "PositionGradientCalculator: setting m_MinLogLikelihoodTimesPrior to: " << m_MinLogLikelihoodTimesPrior << std::endl;
	// std::cout << "Point data: "<<std::endl;
	// std::cout << "y11: "<<y11<<std::endl;
	// std::cout << "y21: "<<y21<<std::endl;
	// std::cout << "y31: "<<y31<<std::endl;
	// std::cout << "y12: "<<y12<<std::endl;
	// std::cout << "y22: "<<y22<<std::endl;
	// std::cout << "y32: "<<y32<<std::endl;
	// std::cout << "y13: "<<y13<<std::endl;
	// std::cout << "y23: "<<y23<<std::endl;
	// std::cout << "y33: "<<y33<<std::endl;
	// std::cout << "y14: "<<y14<<std::endl;
	// std::cout << "y24: "<<y24<<std::endl;
	// std::cout << "y34: "<<y34<<std::endl;
	// std::cout << "Reference triangle info components: "<<std::endl;
	// std::cout << "z11: "<<z11<<std::endl;
	// std::cout << "z21: "<<z21<<std::endl;
	// std::cout << "z31: "<<z31<<std::endl;
	// std::cout << "z41: "<<z41<<std::endl;
	// std::cout << "z12: "<<z12<<std::endl;
	// std::cout << "z22: "<<z22<<std::endl;
	// std::cout << "z32: "<<z32<<std::endl;
	// std::cout << "z42: "<<z42<<std::endl;
	// std::cout << "z13: "<<z13<<std::endl;
	// std::cout << "z23: "<<z23<<std::endl;
	// std::cout << "z33: "<<z33<<std::endl;
	// std::cout << "z43: "<<z43<<std::endl;
	AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
	// std::cout<<"Point id: "<<*(pit)<<std::endl;
    	++pit;
        // std::cout<<"Point id: "<<*(pit)<<std::endl;
    	++pit;
	// std::cout<<"Point id: "<<*(pit)<<std::endl;
        ++pit;
        // std::cout<<"Point id: "<<*(pit)<<std::endl;
        return false;
        }
      


      // Let's define K as inv( J ) * det( J )
      const double  k11 = ( m22*m33 - m23*m32 );
      const double  k12 = -( m12*m33 - m32*m13 );
      const double  k13 = ( m12*m23 - m22*m13 );
      const double  k21 = -( m21*m33 - m31*m23 );
      const double  k22 = ( m11*m33 - m13*m31 );
      const double  k23 = -( m11*m23 - m21*m13 );
      const double  k31 = ( m21*m32 - m31*m22 );
      const double  k32 = -( m11*m32 - m31*m12 );
      const double  k33 = ( m11*m22 - m12*m21 );

      // Trace of J' * J is actually the sum of the squares of the singular values of J: s1^2 + s2^2 + s3^2
      const double  sumOfSquaresOfSingularValuesOfJ = m11*m11 + m12*m12 + m13*m13 + m21*m21 + m22*m22 + m23*m23 + m31*m31 + m32*m32 + m33*m33;

      // Trace of ( inv(J) )' * inv( J ) is actually the sum of the squares of the reciprocals of the singular values
      // of J: 1/s1^2 + 1/s2^2 + 1/s3^2
      const double  traceOfKTransposeTimesK = ( k11*k11 + k12*k12 + k13*k13 + k21*k21 + k22*k22 + k23*k23 + k31*k31 + k32*k32 + k33*k33 );
      const double  sumOfSquaresOfReciprocalsOfSingularValuesOfJ = traceOfKTransposeTimesK / ( detJ * detJ );

      const double  priorCost = referenceVolumeTimesK * ( 1 + detJ ) *
                              ( sumOfSquaresOfSingularValuesOfJ + sumOfSquaresOfReciprocalsOfSingularValuesOfJ - 6 );
      m_MinLogLikelihoodTimesPrior += priorCost;
      //// std::cout << "PositionGradientCalculator: setting m_MinLogLikelihoodTimesPrior to: " << m_MinLogLikelihoodTimesPrior << std::endl;



      //
      // OK, now add contribution to derivatives of Ashburner's prior cost in each of the tetrahedron's vertices
      //
      const double  ddetJdm11 = m22*m33 - m32*m23;
      const double  ddetJdm21 = m13*m32 - m12*m33;
      const double  ddetJdm31 = m12*m23 - m13*m22;
      const double  ddetJdm12 = m31*m23 - m21*m33;
      const double  ddetJdm22 = m11*m33 - m13*m31;
      const double  ddetJdm32 = m13*m21 - m11*m23;
      const double  ddetJdm13 = m21*m32 - m31*m22;
      const double  ddetJdm23 = m12*m31 - m11*m32;
      const double  ddetJdm33 = m11*m22 - m12*m21;


      const double  tmp1 = referenceVolumeTimesK * ( sumOfSquaresOfSingularValuesOfJ + sumOfSquaresOfReciprocalsOfSingularValuesOfJ - 6 );
      const double  tmp2 = 2 * referenceVolumeTimesK * ( 1 + detJ );
      const double  tmp3 = pow( detJ,  3 );


      const double  dcostdm11 = tmp1 * ddetJdm11 + tmp2 * ( m11 + ( ( k22*m33 - k23*m23 - k32*m32 + k33*m22 ) * detJ - traceOfKTransposeTimesK * ddetJdm11 ) / tmp3 );
      const double  dcostdm21 = tmp1 * ddetJdm21 + tmp2 * ( m21 + ( ( -k21*m33 + k23*m13 + k31*m32 - k33*m12 ) * detJ - traceOfKTransposeTimesK * ddetJdm21 ) / tmp3 );
      const double  dcostdm31 = tmp1 * ddetJdm31 + tmp2 * ( m31 + ( ( k21*m23 - k22*m13 - k31*m22 + k32*m12 ) * detJ - traceOfKTransposeTimesK * ddetJdm31 ) / tmp3 );

      const double  dcostdm12 = tmp1 * ddetJdm12 + tmp2 * ( m12 + ( ( -k12*m33 + k13*m23 + k32*m31 - k33*m21 ) * detJ - traceOfKTransposeTimesK * ddetJdm12 ) / tmp3 );
      const double  dcostdm22 = tmp1 * ddetJdm22 + tmp2 * ( m22 + ( ( k11*m33 - k13*m13 - k31*m31 + k33*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm22 ) / tmp3 );
      const double  dcostdm32 = tmp1 * ddetJdm32 + tmp2 * ( m32 + ( ( -k11*m23 + k12*m13 + k31*m21 - k32*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm32 ) / tmp3 );

      const double  dcostdm13 = tmp1 * ddetJdm13 + tmp2 * ( m13 + ( ( k12*m32 - k13*m22 - k22*m31 + k23*m21 ) * detJ - traceOfKTransposeTimesK * ddetJdm13 ) / tmp3 );
      const double  dcostdm23 = tmp1 * ddetJdm23 + tmp2 * ( m23 + ( ( -k11*m32 + k13*m12 + k21*m31 -k23*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm23 ) / tmp3 );
      const double  dcostdm33 = tmp1 * ddetJdm33 + tmp2 * ( m33 + ( ( k11*m22 - k12*m12 - k21*m21 + k22*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm33 ) / tmp3 );


      const double  dcostdy11 = dcostdm11*z11 + dcostdm12*z12 + dcostdm13*z13;
      const double  dcostdy21 = dcostdm21*z11 + dcostdm22*z12 + dcostdm23*z13;
      const double  dcostdy31 = dcostdm31*z11 + dcostdm32*z12 + dcostdm33*z13;

      const double  dcostdy12 = dcostdm11*z21 + dcostdm12*z22 + dcostdm13*z23;
      const double  dcostdy22 = dcostdm21*z21 + dcostdm22*z22 + dcostdm23*z23;
      const double  dcostdy32 = dcostdm31*z21 + dcostdm32*z22 + dcostdm33*z23;

      const double  dcostdy13 = dcostdm11*z31 + dcostdm12*z32 + dcostdm13*z33;
      const double  dcostdy23 = dcostdm21*z31 + dcostdm22*z32 + dcostdm23*z33;
      const double  dcostdy33 = dcostdm31*z31 + dcostdm32*z32 + dcostdm33*z33;

      const double  dcostdy14 = dcostdm11*z41 + dcostdm12*z42 + dcostdm13*z43;
      const double  dcostdy24 = dcostdm21*z41 + dcostdm22*z42 + dcostdm23*z43;
      const double  dcostdy34 = dcostdm31*z41 + dcostdm32*z42 + dcostdm33*z43;


      // Add the stuff to the existing gradients in each of the tetrahedron's four vertices
      AtlasPositionGradientType  gradientContributionToVertex0;
      gradientContributionToVertex0[ 0 ] = dcostdy11;
      gradientContributionToVertex0[ 1 ] = dcostdy21;
      gradientContributionToVertex0[ 2 ] = dcostdy31;
      *m_GradientInVertex0 += gradientContributionToVertex0;
      
      AtlasPositionGradientType  gradientContributionToVertex1;
      gradientContributionToVertex1[ 0 ] = dcostdy12;
      gradientContributionToVertex1[ 1 ] = dcostdy22;
      gradientContributionToVertex1[ 2 ] = dcostdy32;
      *m_GradientInVertex1 += gradientContributionToVertex1;

      AtlasPositionGradientType  gradientContributionToVertex2;
      gradientContributionToVertex2[ 0 ] = dcostdy13;
      gradientContributionToVertex2[ 1 ] = dcostdy23;
      gradientContributionToVertex2[ 2 ] = dcostdy33;
      *m_GradientInVertex2 += gradientContributionToVertex2;

      AtlasPositionGradientType  gradientContributionToVertex3;
      gradientContributionToVertex3[ 0 ] = dcostdy14;
      gradientContributionToVertex3[ 1 ] = dcostdy24;
      gradientContributionToVertex3[ 2 ] = dcostdy34;
      *m_GradientInVertex3 += gradientContributionToVertex3;
      }

    //
    // Finally, precalculate some stuff that will be used over and over again, each time a voxel is visited
    //

    // Get Gamma, defined as Gamma = [ 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 1] * inv( Y )
    // where Y = [ p0 p1 p2 p3; 1 1 1 1 ]
    vnl_matrix_fixed< double, 4, 4 >  Y;
    Y.put( 0, 0, y11 );
    Y.put( 1, 0, y21 );
    Y.put( 2, 0, y31 );
    Y.put( 3, 0, 1.0f );

    Y.put( 0, 1, y12 );
    Y.put( 1, 1, y22 );
    Y.put( 2, 1, y32 );
    Y.put( 3, 1, 1.0f );

    Y.put( 0, 2, y13 );
    Y.put( 1, 2, y23 );
    Y.put( 2, 2, y33 );
    Y.put( 3, 2, 1.0f );

    Y.put( 0, 3, y14 );
    Y.put( 1, 3, y24 );
    Y.put( 2, 3, y34 );
    Y.put( 3, 3, 1.0f );

    vnl_matrix_fixed< double, 4, 4 >  invY = vnl_inverse( Y );
    const double  gamma11 = invY( 1, 0 );
    const double  gamma12 = invY( 1, 1 );
    const double  gamma13 = invY( 1, 2 );
    const double  gamma21 = invY( 2, 0 );
    const double  gamma22 = invY( 2, 1 );
    const double  gamma23 = invY( 2, 2 );
    const double  gamma31 = invY( 3, 0 );
    const double  gamma32 = invY( 3, 1 );
    const double  gamma33 = invY( 3, 2 );


    m_XGradientBasis = AtlasAlphasType( m_AlphasInVertex0.Size() );
    m_YGradientBasis = AtlasAlphasType( m_AlphasInVertex0.Size() );
    m_ZGradientBasis = AtlasAlphasType( m_AlphasInVertex0.Size() );
    for ( unsigned int labelNumber = 0; labelNumber < m_AlphasInVertex0.Size(); labelNumber++ )
      {
      double  alpha0 = m_AlphasInVertex0[ labelNumber ];
      double  alpha1 = m_AlphasInVertex1[ labelNumber ];
      double  alpha2 = m_AlphasInVertex2[ labelNumber ];
      double  alpha3 = m_AlphasInVertex3[ labelNumber ];
  

      m_XGradientBasis[ labelNumber ] = alpha1 * gamma11 + alpha2 * gamma21 + alpha3 * gamma31 - alpha0 * ( gamma11 + gamma21 + gamma31 );
      m_YGradientBasis[ labelNumber ] = alpha1 * gamma12 + alpha2 * gamma22 + alpha3 * gamma32 - alpha0 * ( gamma12 + gamma22 + gamma32 );
      m_ZGradientBasis[ labelNumber ] = alpha1 * gamma13 + alpha2 * gamma23 + alpha3 * gamma33 - alpha0 * ( gamma13 + gamma23 + gamma33 );
      }

    return true;
    }
    
  inline void SetMesh( const AtlasMesh* mesh )
    {
    m_Mesh = mesh;
    
    // Create a container to hold the position gradient
    m_PositionGradient = AtlasPositionGradientContainerType::New();
    
    // Initialize to zero
    AtlasPositionGradientType  zeroEntry( 0.0f );
    
    AtlasMesh::PointsContainer::ConstIterator pointIt = m_Mesh->GetPoints()->Begin();
    while ( pointIt != m_Mesh->GetPoints()->End() )
      {
      m_PositionGradient->InsertElement( pointIt.Index(), zeroEntry );
      ++pointIt;
      }
      
    m_MinLogLikelihoodTimesPrior = 0.0f;
    }
    
  const AtlasMesh* GetMesh() const
    {
    return m_Mesh;
    }
    
  void SetMeans( std::vector< vnl_vector<float> >& means )
    { m_Means = means; }

  //  
  void SetPrecisions( std::vector< vnl_matrix<float> >& precisions )
    {

      m_Precisions.resize(precisions.size()); 

      // We are going to compute 1/sqrt(det(COV)) for all possible covariances given all possible combinations of available channels
      // We use a binary representation for this. For instance, 6 = [1 1 0] means that we have channel 1 not available, but channels 2 and 3 available.
      OneOverSqrtDetCov.resize(precisions.size());
      int nCombos = (int)(pow(2,m_n_channels));
      std::vector<bool> presentChannels(m_n_channels);
      for(int classNumber=0; classNumber<precisions.size(); classNumber++)
      {
        vnl_matrix<float> FullCov=vnl_inverse<float>(precisions[classNumber]);
	OneOverSqrtDetCov[classNumber].resize(nCombos);
        m_Precisions[classNumber].resize(nCombos);
        OneOverSqrtDetCov[classNumber][0]=0;
        for(int n=1; n<nCombos; n++) 
        {
          // decode integer -> binary vector of present channels
          int k = n;
          int nPresent=0;
          for(int c=0; c<m_n_channels; c++)
          {
              if (k & 1) {presentChannels[c]=true; nPresent++; }
              else {presentChannels[c]=false;}
              k = (k >> 1);
          }

          // Extract sub-matrix
          vnl_matrix<float> PartialCov(nPresent,nPresent);
          int r=0; 
          for(int i=0; i<m_n_channels; i++)
          {
            if(presentChannels[i])
            {
              // copy from row i to row r              
              int c=0;
              for(int j=0; j<m_n_channels; j++)
              {
                if(presentChannels[j])
                {
                  // copy from i,j to r,c
                  PartialCov[r][c]=FullCov[i][j];
                  c++;
                }
              }

              r++;
            }
          }
          
          m_Precisions[classNumber][n]=vnl_inverse<float>(PartialCov);
          OneOverSqrtDetCov[classNumber][n]=1.0/sqrt(vnl_determinant(PartialCov));
        }
      }

      // We also compute the constant term for number of channels from 0 to m_n_channels
      m_piTermMultiv.resize(m_n_channels+1);
      for(int i=0; i<=m_n_channels; i++) m_piTermMultiv[i] = pow(2*PI,-0.5*i);
    }

  //
  void SetIgnoreDeformationPrior( bool ignoreDeformationPrior )
    { m_IgnoreDeformationPrior = ignoreDeformationPrior; }

  void SetOnlyDeformationPrior( bool onlyDeformationPrior )
    { m_OnlyDeformationPrior = onlyDeformationPrior; }

  const AtlasPositionGradientContainerType* GetPositionGradient() const
    {
    return m_PositionGradient;
    }

  AtlasPositionGradientContainerType* GetPositionGradient()
    {
    return m_PositionGradient;
    }
  
  double GetMinLogLikelihoodTimesPrior() const
    {
    //// std::cout << "PositionGradientCalculator: returning m_MinLogLikelihoodTimesPrior " << m_MinLogLikelihoodTimesPrior << std::endl;
    return m_MinLogLikelihoodTimesPrior;
    }
        
private:

    
  AtlasAlphasType  m_AlphasInVertex0;
  AtlasAlphasType  m_AlphasInVertex1;
  AtlasAlphasType  m_AlphasInVertex2;
  AtlasAlphasType  m_AlphasInVertex3;
  
  AtlasPositionGradientType*  m_GradientInVertex0;
  AtlasPositionGradientType*  m_GradientInVertex1;
  AtlasPositionGradientType*  m_GradientInVertex2;
  AtlasPositionGradientType*  m_GradientInVertex3;

  double cumsumLogX[10000];
    
  std::vector<itk::Image< float, 3 >::Pointer> m_Images;
  ProbabilityImageType::ConstPointer  m_ProbabilityImage;
  SegmentedImageType::ConstPointer  m_SegmentedImage;
  
  ImageType::IndexType  m_Index;
  AtlasMesh::ConstPointer  m_Mesh;
  AtlasPositionGradientContainerType::Pointer  m_PositionGradient;

  int m_n_channels;
  std::vector<float> m_piTermMultiv;

  std::vector< vnl_vector<float> > m_Means;
  std::vector<std::vector< vnl_matrix<float> > > m_Precisions;
  std::vector< std::vector<float> > OneOverSqrtDetCov;
  bool  m_IgnoreDeformationPrior;
  bool  m_OnlyDeformationPrior;

  AtlasAlphasType  m_XGradientBasis;
  AtlasAlphasType  m_YGradientBasis;
  AtlasAlphasType  m_ZGradientBasis;

  double  m_MinLogLikelihoodTimesPrior;

  std::vector<unsigned char > *m_mapCompToComp;

};

 
 

} // End namespace FragmentProcessor





/**
 *
 */
class AtlasMeshToIntensityImageGradientCalculator :
  public AtlasMeshRasterizor< FragmentProcessor::CalculateGradientToIntensityImage >
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshToIntensityImageGradientCalculator  Self;
  typedef AtlasMeshRasterizor< FragmentProcessor::CalculateGradientToIntensityImage >  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshToIntensityImageGradientCalculator, itk::Object );

  /** Some typedefs */
  typedef Superclass::FragmentProcessorType  FragmentProcessorType;
  typedef Superclass::LabelImageType  LabelImageType;
  typedef FragmentProcessorType::ImageType  ImageType;
  typedef FragmentProcessorType::ProbabilityImageType  ProbabilityImageType;
  typedef FragmentProcessorType::SegmentedImageType SegmentedImageType;
  typedef itk::AffineTransform< double, 3 >  TransformType;
      
  /** */
  void  SetImages( std::vector<itk::Image< float, 3 >::Pointer> images )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetImages( images );
      }

    }

  /** */
  const itk::Image< float, 3 >::Pointer GetImage(int ind) const
    { 
    return this->GetFragmentProcessor().GetImage(ind); 
    }


  /** */
  void  SetProbabilityImage( const ProbabilityImageType* probabilityImage )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetProbabilityImage( probabilityImage );  
      }
    }

  /** */
  const ProbabilityImageType*  GetProbabilityImage() const
    {
    return this->GetFragmentProcessor().GetProbabilityImage();
    }


  /** */
  void  SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetMapCompToComp( mapCompToComp );  
      }
    }

  /** */
  std::vector<unsigned char > * GetMapCompToComp()
    {
    return this->GetFragmentProcessor().GetMapCompToComp();
    }

  /** */
  void  SetSegmentedImage( const SegmentedImageType*  segmentedImage )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetSegmentedImage( segmentedImage );
      }
    }

  /** */
  const SegmentedImageType*  GetSegmentedImage() const
    {
    return this->GetFragmentProcessor().GetSegmentedImage();
    }


  /** */
  void SetMeans( std::vector< vnl_vector<float> >& means )
    { 
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetMeans( means ); 
      }
    }

  /** */ 
  void SetPrecisions( std::vector< vnl_matrix<float> >& precisions )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetPrecisions( precisions ); 
      }
    
    }

  /** */
  void  SetMeshToImageTransform( const TransformType* meshToImageTransform );

  /** */
  void  SetIgnoreDeformationPrior( bool ignoreDeformationPrior )
    { 
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetIgnoreDeformationPrior( ignoreDeformationPrior );
      }
    }  

/** */
  void  SetOnlyDeformationPrior( bool onlyDeformationPrior )
    { 
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetOnlyDeformationPrior( onlyDeformationPrior );
      }
    } 
      
  /** */
  double GetMinLogLikelihoodTimesPrior() const;
  
  /** */
  AtlasPositionGradientContainerType::ConstPointer GetPositionGradient() const
    {
    return this->GetPositionGradient();
    }
  
  /** */
  AtlasPositionGradientContainerType::Pointer GetPositionGradient();
  
  /**  */
  void Rasterize( const AtlasMesh* mesh )
    {
    // Rasterize using multithreading
    Superclass::Rasterize( mesh, true );
    }
  
protected:
  AtlasMeshToIntensityImageGradientCalculator();
  virtual ~AtlasMeshToIntensityImageGradientCalculator();

private:
  AtlasMeshToIntensityImageGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  typedef itk::Matrix< double >  SlidingBoundaryCorrectionMatrixType;
  SlidingBoundaryCorrectionMatrixType  m_SlidingBoundaryCorrectionMatrices[ 8 ]; 

};


} // end namespace kvl

#endif

