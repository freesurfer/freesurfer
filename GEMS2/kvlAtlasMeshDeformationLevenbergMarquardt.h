#ifndef __kvlAtlasMeshDeformationLevenbergMarquardt_h
#define __kvlAtlasMeshDeformationLevenbergMarquardt_h

#include "kvlAtlasMeshRasterizor.h"
#include "itkAffineTransform.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_matrix_fixed.h"
#include "gmm/gmm.h"
#include "itkTimeProbe.h"

#define PI 3.141592654

namespace kvl
{


namespace FragmentProcessor 
{

/**
 *
 */
class CalculateGradientAndHessian
{
public:

  typedef itk::Image< AtlasAlphasType, 3 >  ProbabilityImageType;
  typedef std::vector< double >  GradientType;
  typedef gmm::row_matrix< gmm::wsvector< double > > HessianType;
  typedef itk::AffineTransform< double, 3 >  TransformType;
  typedef itk::Image< float, 3 >  ImageType;
  typedef itk::Image< unsigned char, 3 >  SegmentedImageType;


  CalculateGradientAndHessian()
    {
    m_ProbabilityImage = 0;
    //m_Image = 0;
    m_SegmentedImage = 0;
    m_Mesh = 0;
    m_MinLogLikelihoodTimesPrior = 0;

    m_MeshToImageTransform = 0;

    m_mapCompToComp = 0;

    m_ContributionsFromPreviousTetrahdronFlushedOut = true;

    m_ComputeHessian = true;
    }

  ~CalculateGradientAndHessian() {};

  void SetImages( std::vector<itk::Image< float, 3 >::Pointer> images ) 
    {
    m_Images = images;
    m_n_channels = images.size();
    }

  const itk::Image< float, 3 >::Pointer GetImage(int ind) const
    { return m_Images[ind]; }

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

  const AtlasMesh*  GetMesh() const
    { return m_Mesh; }

  void SetMeans( std::vector< vnl_vector<float> >& means )
    { m_Means = means; }


  // Eugenio: this got a bit more complicated, since we now we want to precompute the 1/sqrt(determinant) for all possible submatrices of each covariance
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




  // Set/Get mapping of collapsed labels.
  void SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
    { m_mapCompToComp = mapCompToComp; }
  std::vector<unsigned char > * GetMapCompToComp()
    { return m_mapCompToComp; }

  inline void operator()( const float& pi0, const float& pi1, const float& pi2, const float& pi3 )
    {
    // This will hold the precious gradient basis from which this voxel's contribution to each
    // of the tetrahedral nodes can be easily calculated
    double  tmpX = 0.0f;
    double  tmpY = 0.0f;
    double  tmpZ = 0.0f;
    double  tmpXX = 0.0f;
    double  tmpXY = 0.0f;
    double  tmpXZ = 0.0f;
    double  tmpYY = 0.0f;
    double  tmpYZ = 0.0f;
    double  tmpZZ = 0.0f;

    if ( m_SegmentedImage && m_mapCompToComp)   // general case with collapsed labels
      {

      double  likelihood =  1e-15;
      // look around labels within this collapsed label
      const SegmentedImageType::PixelType  collapsedClassNumber = m_SegmentedImage->GetPixel( m_Index );
      for(int ind=0; ind<m_mapCompToComp[collapsedClassNumber].size(); ind++)
          {
           const SegmentedImageType::PixelType k = m_mapCompToComp[collapsedClassNumber][ind];

           // Collect the data terms for each vertex
           double  alpha0 = m_AlphasInVertex0[ k ];
           double  alpha1 = m_AlphasInVertex1[ k ];
           double  alpha2 = m_AlphasInVertex2[ k ];
           double  alpha3 = m_AlphasInVertex3[ k ];

           // Add contribution of the likelihood
           likelihood += ( alpha0 * pi0 + alpha1 * pi1 + alpha2 * pi2 + alpha3 * pi3 ) ;

           tmpX += ( m_XGradientBasis[ k ] ); 
           tmpY += ( m_YGradientBasis[ k ] );
           tmpZ += ( m_ZGradientBasis[ k ] );
           
           tmpXX += ( m_XGradientBasis[ k ] * m_XGradientBasis[ k ] );
           tmpXY += ( m_XGradientBasis[ k ] * m_YGradientBasis[ k ] );
           tmpXZ += ( m_XGradientBasis[ k ] * m_ZGradientBasis[ k ] );
           tmpYY += ( m_YGradientBasis[ k ] * m_YGradientBasis[ k ] );
           tmpYZ += ( m_YGradientBasis[ k ] * m_ZGradientBasis[ k ] );
           tmpZZ += ( m_ZGradientBasis[ k ] * m_ZGradientBasis[ k ] );
           }



      double cost = -log( likelihood );

      if (std::isnan(cost))  // Eugenio added this check
        {
          m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
          tmpX = 0; tmpY=0; tmpZ=0; tmpXX=1e-15; tmpYY=1e-15; tmpZZ=1e-15; tmpXY=0; tmpXZ=0; tmpYZ=0;
        }
      else
        {
          m_MinLogLikelihoodTimesPrior += cost;

          const double  tmp1 = 1.0 / likelihood;
          const double  tmp2 = 1.0 / ( likelihood * likelihood );
      
          tmpX = tmpX * tmp1; 
          tmpY = tmpY * tmp1; 
          tmpZ = tmpZ * tmp1;

          tmpXX = tmp2 * tmpXX;
          tmpXY = tmp2 * tmpXY;
          tmpXZ = tmp2 * tmpXZ;
          tmpYY = tmp2 * tmpYY;
          tmpYZ = tmp2 * tmpYZ;
          tmpZZ = tmp2 * tmpZZ; 
        }
        
      }
    else if ( m_SegmentedImage && m_mapCompToComp )  // no collapsed labels
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

      if (std::isnan(cost))  // Eugenio added this check
        {
          m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
          tmpX = 0; tmpY=0; tmpZ=0; tmpXX=1e-15; tmpYY=1e-15; tmpZZ=1e-15; tmpXY=0; tmpXZ=0; tmpYZ=0;
        }
      else
        {
          m_MinLogLikelihoodTimesPrior += cost;

          const double  tmp1 = 1.0 / likelihood;
          const double  tmp2 = 1.0 / ( likelihood * likelihood );

          tmpX += tmp1 * ( m_XGradientBasis[ classNumber ] );
          tmpY += tmp1 * ( m_YGradientBasis[ classNumber ] );
          tmpZ += tmp1 * ( m_ZGradientBasis[ classNumber ] );

          tmpXX += tmp2 * ( m_XGradientBasis[ classNumber ] * m_XGradientBasis[ classNumber ] );
          tmpXY += tmp2 * ( m_XGradientBasis[ classNumber ] * m_YGradientBasis[ classNumber ] );
          tmpXZ += tmp2 * ( m_XGradientBasis[ classNumber ] * m_ZGradientBasis[ classNumber ] );
          tmpYY += tmp2 * ( m_YGradientBasis[ classNumber ] * m_YGradientBasis[ classNumber ] );
          tmpYZ += tmp2 * ( m_YGradientBasis[ classNumber ] * m_ZGradientBasis[ classNumber ] );
          tmpZZ += tmp2 * ( m_ZGradientBasis[ classNumber ] * m_ZGradientBasis[ classNumber ] );
        }

      }
    else if ( m_ProbabilityImage )
      {
      // Loop over all classes
      const AtlasAlphasType&  weights = m_ProbabilityImage->GetPixel( m_Index );
      double  cost = 0.0f;
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
        const double  tmp2 = weight / ( likelihood * likelihood );
        tmpX += tmp1 * ( m_XGradientBasis[ classNumber ] );
        tmpY += tmp1 * ( m_YGradientBasis[ classNumber ] );
        tmpZ += tmp1 * ( m_ZGradientBasis[ classNumber ] );

        tmpXX += tmp2 * ( m_XGradientBasis[ classNumber ] * m_XGradientBasis[ classNumber ] );
        tmpXY += tmp2 * ( m_XGradientBasis[ classNumber ] * m_YGradientBasis[ classNumber ] );
        tmpXZ += tmp2 * ( m_XGradientBasis[ classNumber ] * m_ZGradientBasis[ classNumber ] );
        tmpYY += tmp2 * ( m_YGradientBasis[ classNumber ] * m_YGradientBasis[ classNumber ] );
        tmpYZ += tmp2 * ( m_YGradientBasis[ classNumber ] * m_ZGradientBasis[ classNumber ] );
        tmpZZ += tmp2 * ( m_ZGradientBasis[ classNumber ] * m_ZGradientBasis[ classNumber ] );
        }
      m_MinLogLikelihoodTimesPrior += cost;
      }
    else
      {
      //std::cout << "!!!!!! Using original image to directly deform to !!!!!" << std::endl;
      vnl_vector<ImageType::PixelType> intensity_v;
      ImageType::PixelType intensity=0.0f;
      std::vector<bool> isThere(m_n_channels);
      // Eugenio: these 2 are only useful for multivariate Gaussian, but whatever...
      int nPresent=0;
      int index=0;

      if(m_n_channels > 1)
      {
	//std::cout<<"More than one channel!"<<std::endl;
        int aux=1;
        vnl_vector<ImageType::PixelType> aux_v(m_n_channels,0.0f);
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

      double likelihood = 0.0f;
      double  intensityXGradientBasis = 0.0f;
      double  intensityYGradientBasis = 0.0f;
      double  intensityZGradientBasis = 0.0f;
      for ( unsigned int classNumber = 0; classNumber < m_AlphasInVertex0.Size(); classNumber++ )
        {
        // Evaluate the Gaussian of this class at the intensity of this pixel
        double gauss = 0.0f;
        if(m_n_channels == 1)
        {
          gauss = exp( -0.5 * m_Precisions[classNumber][1][0][0] * pow( intensity - m_Means[ classNumber ][0] , 2 )) / 
                                      sqrt( 2 * PI / m_Precisions[ classNumber ][1][0][0] ); //NOTE: These are precisions now (precision = 1/variance) so we multiply
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
        likelihood += gauss * ( alpha0 * pi0 + alpha1 * pi1 + alpha2 * pi2 + alpha3 * pi3 ) + 1e-15;

        //
        intensityXGradientBasis += gauss * m_XGradientBasis[ classNumber ];
        intensityYGradientBasis += gauss * m_YGradientBasis[ classNumber ];
        intensityZGradientBasis += gauss * m_ZGradientBasis[ classNumber ];
        }

      m_MinLogLikelihoodTimesPrior -= log( likelihood );


      //
      tmpX = intensityXGradientBasis / likelihood;
      tmpY = intensityYGradientBasis / likelihood;
      tmpZ = intensityZGradientBasis / likelihood;

      tmpXX = tmpX * tmpX;
      tmpXY = tmpX * tmpY;
      tmpXZ = tmpX * tmpZ;
      tmpYY = tmpY * tmpY;
      tmpYZ = tmpY * tmpZ;
      tmpZZ = tmpZ * tmpZ;
      }


    const double  pis[] = { pi0, pi1, pi2, pi3 };
    int  indexInHessianDataContribution = 0;
    for ( int vertexNumber = 0; vertexNumber < 4; vertexNumber++ )
      {
      const double  pi = pis[ vertexNumber ];

      // Add contribution to gradient
      m_GradientContribution[ vertexNumber * 3 ]     += tmpX * pi;
      m_GradientContribution[ vertexNumber * 3 + 1 ] += tmpY * pi;
      m_GradientContribution[ vertexNumber * 3 + 2 ] += tmpZ * pi;

      if ( m_ComputeHessian )
        {
        for ( int vertexNumber2 = vertexNumber; vertexNumber2 < 4; vertexNumber2++ )
          {
          // Add contribution to Hessian in which second-order terms are ignored
          const double  pi2 = pis[ vertexNumber2 ];
          const double  piTimesPi2 = pi * pi2;

          m_HessianDataContribution[ indexInHessianDataContribution ]   += tmpXX * piTimesPi2;
          m_HessianDataContribution[ indexInHessianDataContribution+1 ] += tmpXY * piTimesPi2;
          m_HessianDataContribution[ indexInHessianDataContribution+2 ] += tmpXZ * piTimesPi2;
          m_HessianDataContribution[ indexInHessianDataContribution+3 ] += tmpYY * piTimesPi2;
          m_HessianDataContribution[ indexInHessianDataContribution+4 ] += tmpYZ * piTimesPi2;
          m_HessianDataContribution[ indexInHessianDataContribution+5 ] += tmpZZ * piTimesPi2;

          indexInHessianDataContribution += 6;
          }
        } // End test if we need to compute the Hessian

      }

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
    //std::cout << "Starting tetrahedron with id " << cellId << std::endl;

    this->FlushContributionsFromPreviousTetrahedron();

    for ( int i = 0; i < 60; i++ )
      {
      m_HessianDataContribution[ i ] = 0.0f;
      }
    for ( int i = 0; i < 12; i++ )
      {
      m_GradientContribution[ i ] = 0.0f;
      }


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
    m_IndicesForVertex0[ 0 ] = m_ComponentXLookupTable[ *pit ];
    m_IndicesForVertex0[ 1 ] = m_ComponentYLookupTable[ *pit ];
    m_IndicesForVertex0[ 2 ] = m_ComponentZLookupTable[ *pit ];

    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    const double y12 = p[ 0 ];
    const double y22 = p[ 1 ];
    const double y32 = p[ 2 ];
    m_AlphasInVertex1 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_IndicesForVertex1[ 0 ] = m_ComponentXLookupTable[ *pit ];
    m_IndicesForVertex1[ 1 ] = m_ComponentYLookupTable[ *pit ];
    m_IndicesForVertex1[ 2 ] = m_ComponentZLookupTable[ *pit ];

    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    const double y13 = p[ 0 ];
    const double y23 = p[ 1 ];
    const double y33 = p[ 2 ];
    m_AlphasInVertex2 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_IndicesForVertex2[ 0 ] = m_ComponentXLookupTable[ *pit ];
    m_IndicesForVertex2[ 1 ] = m_ComponentYLookupTable[ *pit ];
    m_IndicesForVertex2[ 2 ] = m_ComponentZLookupTable[ *pit ];

    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    const double y14 = p[ 0 ];
    const double y24 = p[ 1 ];
    const double y34 = p[ 2 ];
    m_AlphasInVertex3 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_IndicesForVertex3[ 0 ] = m_ComponentXLookupTable[ *pit ];
    m_IndicesForVertex3[ 1 ] = m_ComponentYLookupTable[ *pit ];
    m_IndicesForVertex3[ 2 ] = m_ComponentZLookupTable[ *pit ];




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
//      std::cout << "Oooooops: tetrahedron " << cellId << " has managed to turn bad (detJ: " << detJ << ")" << std::endl;

      // std::cout << std::endl;
      // std::cout << "YfirstThreeRows: [ " << y11 << " " << y12 << " " << y13 << " " << y14 << ";" << std::endl;
      // std::cout << "                   " << y21 << " " << y22 << " " << y23 << " " << y24 << ";" << std::endl;
      // std::cout << "                   " << y31 << " " << y32 << " " << y33 << " " << y34 << "]" << std::endl;
      // 
      // std::cout << std::endl;
      // std::cout << "ZfirstThreeColumns: [ " << z11 << " " << z12 << " " << z13 << ";" << std::endl;
      // std::cout << "                      " << z21 << " " << z22 << " " << z23 << ";" << std::endl;
      // std::cout << "                      " << z31 << " " << z32 << " " << z33 << ";" << std::endl;
      // std::cout << "                      " << z41 << " " << z42 << " " << z43 << "]" << std::endl;

      m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
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
    //std::cout << "PositionGradientCalculator: setting m_MinLogLikelihoodTimesPrior to: " << m_MinLogLikelihoodTimesPrior << std::endl;



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


    double  dcostdy11 = dcostdm11*z11 + dcostdm12*z12 + dcostdm13*z13;
    double  dcostdy21 = dcostdm21*z11 + dcostdm22*z12 + dcostdm23*z13;
    double  dcostdy31 = dcostdm31*z11 + dcostdm32*z12 + dcostdm33*z13;

    double  dcostdy12 = dcostdm11*z21 + dcostdm12*z22 + dcostdm13*z23;
    double  dcostdy22 = dcostdm21*z21 + dcostdm22*z22 + dcostdm23*z23;
    double  dcostdy32 = dcostdm31*z21 + dcostdm32*z22 + dcostdm33*z23;

    double  dcostdy13 = dcostdm11*z31 + dcostdm12*z32 + dcostdm13*z33;
    double  dcostdy23 = dcostdm21*z31 + dcostdm22*z32 + dcostdm23*z33;
    double  dcostdy33 = dcostdm31*z31 + dcostdm32*z32 + dcostdm33*z33;

    double  dcostdy14 = dcostdm11*z41 + dcostdm12*z42 + dcostdm13*z43;
    double  dcostdy24 = dcostdm21*z41 + dcostdm22*z42 + dcostdm23*z43;
    double  dcostdy34 = dcostdm31*z41 + dcostdm32*z42 + dcostdm33*z43;


    if ( m_MeshToImageTransform )
      {
      // If we just use the dcostdy's defined above, we calculate gradients of the prior term w.r.t. to
      // the image grid coordinate system. In order to implement the sliding boundary conditions of
      // the mesh deformation, we'd like the gradients w.r.t. the mesh coordinate system. To accomplish
      // this, simply transform things around
      const double  dcostdu11 = ( m_MeshToImageTransform->GetMatrix() )( 0, 0 ) * dcostdy11 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 0 ) * dcostdy21 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 0 ) * dcostdy31;
      const double  dcostdu21 = ( m_MeshToImageTransform->GetMatrix() )( 0, 1 ) * dcostdy11 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 1 ) * dcostdy21 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 1 ) * dcostdy31;
      const double  dcostdu31 = ( m_MeshToImageTransform->GetMatrix() )( 0, 2 ) * dcostdy11 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 2 ) * dcostdy21 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 2 ) * dcostdy31;

      const double  dcostdu12 = ( m_MeshToImageTransform->GetMatrix() )( 0, 0 ) * dcostdy12 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 0 ) * dcostdy22 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 0 ) * dcostdy32;
      const double  dcostdu22 = ( m_MeshToImageTransform->GetMatrix() )( 0, 1 ) * dcostdy12 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 1 ) * dcostdy22 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 1 ) * dcostdy32;
      const double  dcostdu32 = ( m_MeshToImageTransform->GetMatrix() )( 0, 2 ) * dcostdy12 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 2 ) * dcostdy22 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 2 ) * dcostdy32;

      const double  dcostdu13 = ( m_MeshToImageTransform->GetMatrix() )( 0, 0 ) * dcostdy13 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 0 ) * dcostdy23 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 0 ) * dcostdy33;
      const double  dcostdu23 = ( m_MeshToImageTransform->GetMatrix() )( 0, 1 ) * dcostdy13 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 1 ) * dcostdy23 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 1 ) * dcostdy33;
      const double  dcostdu33 = ( m_MeshToImageTransform->GetMatrix() )( 0, 2 ) * dcostdy13 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 2 ) * dcostdy23 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 2 ) * dcostdy33;

      const double  dcostdu14 = ( m_MeshToImageTransform->GetMatrix() )( 0, 0 ) * dcostdy14 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 0 ) * dcostdy24 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 0 ) * dcostdy34;
      const double  dcostdu24 = ( m_MeshToImageTransform->GetMatrix() )( 0, 1 ) * dcostdy14 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 1 ) * dcostdy24 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 1 ) * dcostdy34;
      const double  dcostdu34 = ( m_MeshToImageTransform->GetMatrix() )( 0, 2 ) * dcostdy14 +
                               ( m_MeshToImageTransform->GetMatrix() )( 1, 2 ) * dcostdy24 +
                               ( m_MeshToImageTransform->GetMatrix() )( 2, 2 ) * dcostdy34;


      dcostdy11 = dcostdu11;
      dcostdy21 = dcostdu21;
      dcostdy31 = dcostdu31;

      dcostdy12 = dcostdu12;
      dcostdy22 = dcostdu22;
      dcostdy32 = dcostdu32;

      dcostdy13 = dcostdu13;
      dcostdy23 = dcostdu23;
      dcostdy33 = dcostdu33;

      dcostdy14 = dcostdu14;
      dcostdy24 = dcostdu24;
      dcostdy34 = dcostdu34;
      }

    // Add contributions of prior to gradient and hessian
    double  priorGradient[ 12 ];
    priorGradient[ 0 ]  = dcostdy11;
    priorGradient[ 1 ]  = dcostdy21;
    priorGradient[ 2 ]  = dcostdy31;
    priorGradient[ 3 ]  = dcostdy12;
    priorGradient[ 4 ]  = dcostdy22;
    priorGradient[ 5 ]  = dcostdy32;
    priorGradient[ 6 ]  = dcostdy13;
    priorGradient[ 7 ]  = dcostdy23;
    priorGradient[ 8 ]  = dcostdy33;
    priorGradient[ 9 ]  = dcostdy14;
    priorGradient[ 10 ] = dcostdy24;
    priorGradient[ 11 ] = dcostdy34;

    const double  tmp = 1.0f / ( 2 * priorCost + 1e-15 );

    int  indexInHessianPriorContribution = 0;
    for ( int i = 0; i < 12; i++ )
      {
      const double  gradientComponent = priorGradient[ i ];

      // Gradient
      m_GradientContribution[ i ] += gradientComponent;

      if ( m_ComputeHessian )
        {
        for ( int j = i; j < 12; j++ )
          {
          const double  gradientComponent2 = priorGradient[ j ];

          // Hessian
          m_HessianPriorContribution[ indexInHessianPriorContribution ] = tmp * gradientComponent * gradientComponent2;
          indexInHessianPriorContribution++;
          }
        } // End test if we need to compute the Hessian

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
    double  gamma11 = invY( 1, 0 );
    double  gamma12 = invY( 1, 1 );
    double  gamma13 = invY( 1, 2 );
    double  gamma21 = invY( 2, 0 );
    double  gamma22 = invY( 2, 1 );
    double  gamma23 = invY( 2, 2 );
    double  gamma31 = invY( 3, 0 );
    double  gamma32 = invY( 3, 1 );
    double  gamma33 = invY( 3, 2 );


    if ( m_MeshToImageTransform )
      {
      // If we just use the gamma's defined above, we calculate gradients of the image term w.r.t. to
      // the image grid coordinate system. In order to implement the sliding boundary conditions of
      // the mesh deformation, we'd like the gradients w.r.t. the mesh coordinate system. To accomplish
      // this, simply transform things around
      const double  transformedGamma11 = ( m_MeshToImageTransform->GetMatrix() )( 0, 0 ) * gamma11 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 1, 0 ) * gamma12 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 2, 0 ) * gamma13;
      const double  transformedGamma21 = ( m_MeshToImageTransform->GetMatrix() )( 0, 0 ) * gamma21 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 1, 0 ) * gamma22 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 2, 0 ) * gamma23;
      const double  transformedGamma31 = ( m_MeshToImageTransform->GetMatrix() )( 0, 0 ) * gamma31 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 1, 0 ) * gamma32 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 2, 0 ) * gamma33;

      const double  transformedGamma12 = ( m_MeshToImageTransform->GetMatrix() )( 0, 1 ) * gamma11 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 1, 1 ) * gamma12 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 2, 1 ) * gamma13;
      const double  transformedGamma22 = ( m_MeshToImageTransform->GetMatrix() )( 0, 1 ) * gamma21 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 1, 1 ) * gamma22 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 2, 1 ) * gamma23;
      const double  transformedGamma32 = ( m_MeshToImageTransform->GetMatrix() )( 0, 1 ) * gamma31 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 1, 1 ) * gamma32 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 2, 1 ) * gamma33;

      const double  transformedGamma13 = ( m_MeshToImageTransform->GetMatrix() )( 0, 2 ) * gamma11 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 1, 2 ) * gamma12 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 2, 2 ) * gamma13;
      const double  transformedGamma23 = ( m_MeshToImageTransform->GetMatrix() )( 0, 2 ) * gamma21 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 1, 2 ) * gamma22 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 2, 2 ) * gamma23;
      const double  transformedGamma33 = ( m_MeshToImageTransform->GetMatrix() )( 0, 2 ) * gamma31 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 1, 2 ) * gamma32 +
                                        ( m_MeshToImageTransform->GetMatrix() )( 2, 2 ) * gamma33;


      gamma11 = transformedGamma11;
      gamma12 = transformedGamma12;
      gamma13 = transformedGamma13;

      gamma21 = transformedGamma21;
      gamma22 = transformedGamma22;
      gamma23 = transformedGamma23;

      gamma31 = transformedGamma31;
      gamma32 = transformedGamma32;
      gamma33 = transformedGamma33;
      }



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



    //
    m_ContributionsFromPreviousTetrahdronFlushedOut = false; // We've started a new tetrahedron, so it's can't have been flushe out yet

    //std::cout << "Done setting up tetrahedron" << std::endl;
    return true;
    }

  inline void SetMesh( const AtlasMesh* mesh )
    {
    m_ContributionsFromPreviousTetrahdronFlushedOut = true;
    m_MinLogLikelihoodTimesPrior = 0.0f;

    //std::cout << "mesh: " << mesh << std::endl;
    //std::cout << "m_Mesh.GetPointer(): " << m_Mesh.GetPointer() << std::endl;
    if ( mesh == m_Mesh.GetPointer() )
      {
      //std::cout << "Not initializing any more because that's costly" << std::endl;
#if 1    
    
      gmm::fill( m_Gradient, 0.0 );

      for ( std::vector< gmm::wsvector< double > >::iterator  hesIt = m_Hessian.begin();
            hesIt != m_Hessian.end(); ++hesIt )
        {
        for ( gmm::linalg_traits< gmm::wsvector< double > >::iterator  hesRowIt = gmm::vect_begin( *hesIt );
              hesRowIt != gmm::vect_end( *hesIt ); ++hesRowIt )
          {
          *hesRowIt = 1e-12;
          }
        }
#else
      gmm::clear( m_Gradient );
      gmm::clear( m_Hessian );
#endif

      return;
      }



    m_Mesh = mesh;

    // Precalculate lookup table for the mesh nodes, mapping each component (x, y, and z)
    // to a continuous index to be used in the Levenberg-Marquardt calculations. In order to
    // respect the sliding boundary conditions of mesh deformation, certain components
    // are "knocked out", i.e., their step sizes are kept to zero and therefore not part
    // of the Levenberg-Marquardt step size calculation. In order to implement this efficiently,
    // we'll assign them a sentinel value: we make the Hessian and the gradient one index
    // bigger than they're supposed to be, map all "knocked out" componenents to the last
    // entry of these bigger Hessian and gradients, and then strip off the extra parts when
    // we actually solve the Levenberg-Marquardt linear system of equations.
    m_ComponentXLookupTable.clear();
    m_ComponentYLookupTable.clear();
    m_ComponentZLookupTable.clear();
    int  currentIndex = 0;
    for ( AtlasMesh::PointDataContainer::ConstIterator  paramIt = m_Mesh->GetPointData()->Begin();
          paramIt != m_Mesh->GetPointData()->End(); ++paramIt )
      {
      if ( paramIt.Value().m_CanMoveX )
        {
        m_ComponentXLookupTable[ paramIt.Index() ] = currentIndex;
        currentIndex++;
        }
      else
        {
        m_ComponentXLookupTable[ paramIt.Index() ] = -1; // Sentinel value; to be replace later
        }

      if ( paramIt.Value().m_CanMoveY )
        {
        m_ComponentYLookupTable[ paramIt.Index() ] = currentIndex;
        currentIndex++;
        }
      else
        {
        m_ComponentYLookupTable[ paramIt.Index() ] = -1; // Sentinel value; to be replace later
        }

      if ( paramIt.Value().m_CanMoveZ )
        {
        m_ComponentZLookupTable[ paramIt.Index() ] = currentIndex;
        currentIndex++;
        }
      else
        {
        m_ComponentZLookupTable[ paramIt.Index() ] = -1; // Sentinel value; to be replace later
        }

      }


    const int  numberOfEntries = currentIndex + 1; // Extra entry for "knocked out" components

    // Now that we now how many entries we'll have in the expanded Hessian and gradient,
    // replace the sentinel value with the last index
    for ( std::map< AtlasMesh::PointIdentifier, int >::iterator  it = m_ComponentXLookupTable.begin();
          it != m_ComponentXLookupTable.end(); ++it )
      {
      if ( it->second == -1 )
        {
        it->second = ( numberOfEntries - 1 );
        }
      }
    for ( std::map< AtlasMesh::PointIdentifier, int >::iterator  it = m_ComponentYLookupTable.begin();
          it != m_ComponentYLookupTable.end(); ++it )
      {
      if ( it->second == -1 )
        {
        it->second = ( numberOfEntries - 1 );
        }
      }
    for ( std::map< AtlasMesh::PointIdentifier, int >::iterator  it = m_ComponentZLookupTable.begin();
          it != m_ComponentZLookupTable.end(); ++it )
      {
      if ( it->second == -1 )
        {
        it->second = ( numberOfEntries - 1 );
        }
      }



    //std::cout << "numberOfEntries: " << numberOfEntries << std::endl;


    // Now that we have the lookup tables, we know the size of the Hessian matrix
    // and gradient vector. Allocate them.
    m_Gradient = GradientType( numberOfEntries );
    m_Hessian = HessianType( numberOfEntries, numberOfEntries );

#if 1    
    // Prefill the Hessian with tiny values in order to generate the correct sparseness structure
    // across threads - since each thread only visits a subset of tetrahedra and therefore has
    // only a subset of non-zero entries in the Hessian, it would otherwise be a hassle to combine
    // them later on
    for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = mesh->GetCells()->Begin();
          cellIt != mesh->GetCells()->End(); ++cellIt )
      {
      AtlasMesh::CellType*  cell = cellIt.Value();

      if( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
        continue;
        }
          
      AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
      int  indicesForVertex0[ 3 ];
      indicesForVertex0[ 0 ] = m_ComponentXLookupTable[ *pit ];
      indicesForVertex0[ 1 ] = m_ComponentYLookupTable[ *pit ];
      indicesForVertex0[ 2 ] = m_ComponentZLookupTable[ *pit ];


      ++pit;
      int  indicesForVertex1[ 3 ];
      indicesForVertex1[ 0 ] = m_ComponentXLookupTable[ *pit ];
      indicesForVertex1[ 1 ] = m_ComponentYLookupTable[ *pit ];
      indicesForVertex1[ 2 ] = m_ComponentZLookupTable[ *pit ];


      ++pit;
      int  indicesForVertex2[ 3 ];
      indicesForVertex2[ 0 ] = m_ComponentXLookupTable[ *pit ];
      indicesForVertex2[ 1 ] = m_ComponentYLookupTable[ *pit ];
      indicesForVertex2[ 2 ] = m_ComponentZLookupTable[ *pit ];


      ++pit;
      int  indicesForVertex3[ 3 ];
      indicesForVertex3[ 0 ] = m_ComponentXLookupTable[ *pit ];
      indicesForVertex3[ 1 ] = m_ComponentYLookupTable[ *pit ];
      indicesForVertex3[ 2 ] = m_ComponentZLookupTable[ *pit ];


      const int  componentXIndices[] = { indicesForVertex0[ 0 ], indicesForVertex1[ 0 ],
                                        indicesForVertex2[ 0 ], indicesForVertex3[ 0 ] };
      const int  componentYIndices[] = { indicesForVertex0[ 1 ], indicesForVertex1[ 1 ],
                                        indicesForVertex2[ 1 ], indicesForVertex3[ 1 ] };
      const int  componentZIndices[] = { indicesForVertex0[ 2 ], indicesForVertex1[ 2 ],
                                        indicesForVertex2[ 2 ], indicesForVertex3[ 2 ] };


      for ( int vertexNumber = 0; vertexNumber < 4; vertexNumber++ )
        {
        const int  indexX = componentXIndices[ vertexNumber ];
        const int  indexY = componentYIndices[ vertexNumber ];
        const int  indexZ = componentZIndices[ vertexNumber ];


        for ( int vertexNumber2 = 0; vertexNumber2 < 4; vertexNumber2++ )
          {
          const int  index2X = componentXIndices[ vertexNumber2 ];
          const int  index2Y = componentYIndices[ vertexNumber2 ];
          const int  index2Z = componentZIndices[ vertexNumber2 ];

          // Fill in first row of the 3x3 matrix
          const double eps = 1e-12;
          m_Hessian( indexX, index2X ) = eps;
          m_Hessian( indexX, index2Y ) = eps;
          m_Hessian( indexX, index2Z ) = eps;

          // Fill in second row of the 3x3 matrix
          m_Hessian( indexY, index2X ) = eps;
          m_Hessian( indexY, index2Y ) = eps;
          m_Hessian( indexY, index2Z ) = eps;

          // Fill in third row of the 3x3 matrix
          m_Hessian( indexZ, index2X ) = eps;
          m_Hessian( indexZ, index2Y ) = eps;
          m_Hessian( indexZ, index2Z ) = eps;

          }

        }
    
    
      } // End loop over all tetrahedra
#endif      
      
    }

  void  SetMeshToImageTransform( const TransformType* meshToImageTransform )
    { m_MeshToImageTransform = meshToImageTransform; }

  const TransformType* GetMeshToImageTransform() const
    { return m_MeshToImageTransform; }

  const GradientType& GetGradient() const
    {
    ( const_cast< CalculateGradientAndHessian* >( this ) )->FlushContributionsFromPreviousTetrahedron();

    return m_Gradient;
    }

  const HessianType& GetHessian() const
    {
    ( const_cast< CalculateGradientAndHessian* >( this ) )->FlushContributionsFromPreviousTetrahedron();

    return m_Hessian;
    }

  double GetMinLogLikelihoodTimesPrior() const
    {
    //std::cout << "PositionGradientCalculator: returning m_MinLogLikelihoodTimesPrior " << m_MinLogLikelihoodTimesPrior << std::endl;
    return m_MinLogLikelihoodTimesPrior;
    }

  const std::map< AtlasMesh::PointIdentifier, int >&  GetComponentXLookupTable() const
    {
    return m_ComponentXLookupTable;
    }

  const std::map< AtlasMesh::PointIdentifier, int >&  GetComponentYLookupTable() const
    {
    return m_ComponentYLookupTable;
    }

  const std::map< AtlasMesh::PointIdentifier, int >&  GetComponentZLookupTable() const
    {
    return m_ComponentZLookupTable;
    }

  void  SetComputeHessian( bool computeHessian )
    {
    m_ComputeHessian = computeHessian;
    }  


private:
  // If a symmetric matrix A = [ a11 a12 a13 ...; a12 a22 a23 ...; ] is stored as a vector b = [ a11 a12 a13 a22 a23 ... ],
  // what is the index in b corresponding to matrix element ( i ,j )?
  static int  GetIndexInUpperDiagonalMatrixRepresentation( int rowNumber, int colNumber, int matrixSize )
    {
    int  index = 0;
    if ( rowNumber <= colNumber )
      {
      const int  offset = static_cast< int >( rowNumber * ( rowNumber + 1 ) / 2.0f + ( matrixSize - rowNumber ) * rowNumber + 0.5 );
      index = offset + ( colNumber - rowNumber );
      }
    else
      {
      const int  offset = static_cast< int >( colNumber * ( colNumber + 1) / 2.0f + ( matrixSize - colNumber ) * colNumber + 0.5 );
      index = offset + ( rowNumber - colNumber );
      }

    return index;
    }

  //
  void  FlushContributionsFromPreviousTetrahedron()
    {
    if ( m_ContributionsFromPreviousTetrahdronFlushedOut )
      {
      // Already flushed out
      return;
      }  

    // Flush out 
    const int  componentXIndices[] = { m_IndicesForVertex0[ 0 ], m_IndicesForVertex1[ 0 ],
                                       m_IndicesForVertex2[ 0 ], m_IndicesForVertex3[ 0 ] };
    const int  componentYIndices[] = { m_IndicesForVertex0[ 1 ], m_IndicesForVertex1[ 1 ],
                                       m_IndicesForVertex2[ 1 ], m_IndicesForVertex3[ 1 ] };
    const int  componentZIndices[] = { m_IndicesForVertex0[ 2 ], m_IndicesForVertex1[ 2 ],
                                       m_IndicesForVertex2[ 2 ], m_IndicesForVertex3[ 2 ] };


    for ( int vertexNumber = 0; vertexNumber < 4; vertexNumber++ )
      {
      const int  indexX = componentXIndices[ vertexNumber ];
      const int  indexY = componentYIndices[ vertexNumber ];
      const int  indexZ = componentZIndices[ vertexNumber ];

      // Add contribution to gradient 
      m_Gradient[ indexX ] += m_GradientContribution[ vertexNumber * 3 ];
      m_Gradient[ indexY ] += m_GradientContribution[ vertexNumber * 3 + 1 ];
      m_Gradient[ indexZ ] += m_GradientContribution[ vertexNumber * 3 + 2 ];

      if ( m_ComputeHessian )
        {
        for ( int vertexNumber2 = 0; vertexNumber2 < 4; vertexNumber2++ )
          {
          const int  index2X = componentXIndices[ vertexNumber2 ];
          const int  index2Y = componentYIndices[ vertexNumber2 ];
          const int  index2Z = componentZIndices[ vertexNumber2 ];

          // Map ( vertexNumber, vertexNumber2 ) into an index of upper diagonal representation of a 4x4 matrix
          const int  fourTimesFourIndex = GetIndexInUpperDiagonalMatrixRepresentation( vertexNumber, vertexNumber2, 4 );

          // Index of first element of the 6 unique values corresponding to this pair of vertices in m_HessianDataContribution
          // is simply the above index times 6
          const int  startIndexInHessianDataContribution = fourTimesFourIndex * 6;


          // Fill in first row of the 3x3 matrix
          m_Hessian( indexX, index2X ) += m_HessianDataContribution[ startIndexInHessianDataContribution ] +
                                              m_HessianPriorContribution[ GetIndexInUpperDiagonalMatrixRepresentation( vertexNumber*3,
                                                                                                                        vertexNumber2*3, 12 ) ];
          m_Hessian( indexX, index2Y ) += m_HessianDataContribution[ startIndexInHessianDataContribution+1 ] +
                                              m_HessianPriorContribution[ GetIndexInUpperDiagonalMatrixRepresentation( vertexNumber*3,
                                                                                                                        vertexNumber2*3+1, 12 ) ];
          m_Hessian( indexX, index2Z ) += m_HessianDataContribution[ startIndexInHessianDataContribution+2 ] +
                                              m_HessianPriorContribution[ GetIndexInUpperDiagonalMatrixRepresentation( vertexNumber*3,
                                                                                                                        vertexNumber2*3+2, 12 ) ];

          // Fill in second row of the 3x3 matrix
          m_Hessian( indexY, index2X ) += m_HessianDataContribution[ startIndexInHessianDataContribution+1 ] +
                                              m_HessianPriorContribution[ GetIndexInUpperDiagonalMatrixRepresentation( vertexNumber*3+1,
                                                                                                                        vertexNumber2*3, 12 ) ];
          m_Hessian( indexY, index2Y ) += m_HessianDataContribution[ startIndexInHessianDataContribution+3 ] +
                                              m_HessianPriorContribution[ GetIndexInUpperDiagonalMatrixRepresentation( vertexNumber*3+1,
                                                                                                                        vertexNumber2*3+1, 12 ) ];
          m_Hessian( indexY, index2Z ) += m_HessianDataContribution[ startIndexInHessianDataContribution+4 ] +
                                              m_HessianPriorContribution[ GetIndexInUpperDiagonalMatrixRepresentation( vertexNumber*3+1,
                                                                                                                        vertexNumber2*3+2, 12 ) ];

          // Fill in third row of the 3x3 matrix
          m_Hessian( indexZ, index2X ) += m_HessianDataContribution[ startIndexInHessianDataContribution+2 ] +
                                              m_HessianPriorContribution[ GetIndexInUpperDiagonalMatrixRepresentation( vertexNumber*3+2,
                                                                                                                        vertexNumber2*3, 12 ) ];
          m_Hessian( indexZ, index2Y ) += m_HessianDataContribution[ startIndexInHessianDataContribution+4 ] +
                                              m_HessianPriorContribution[ GetIndexInUpperDiagonalMatrixRepresentation( vertexNumber*3+2,
                                                                                                                        vertexNumber2*3+1, 12 ) ];
          m_Hessian( indexZ, index2Z ) += m_HessianDataContribution[ startIndexInHessianDataContribution+5 ] +
                                              m_HessianPriorContribution[ GetIndexInUpperDiagonalMatrixRepresentation( vertexNumber*3+2,
                                                                                                                        vertexNumber2*3+2, 12 ) ];

          }

        } // End test if we will bother about computing the Hessian

      }
      
    m_ContributionsFromPreviousTetrahdronFlushedOut = true;
    }

 

  ProbabilityImageType::ConstPointer  m_ProbabilityImage;
  std::vector<ImageType::Pointer> m_Images;
  SegmentedImageType::ConstPointer  m_SegmentedImage;
  ImageType::IndexType  m_Index;

  int m_n_channels;
  std::vector<float> m_piTermMultiv;
  std::vector<unsigned char> *m_mapCompToComp;

  std::vector< vnl_vector<float> > m_Means;
  std::vector< std::vector<vnl_matrix<float> > > m_Precisions;
  std::vector< std::vector<float> > OneOverSqrtDetCov;

  AtlasAlphasType  m_AlphasInVertex0;
  AtlasAlphasType  m_AlphasInVertex1;
  AtlasAlphasType  m_AlphasInVertex2;
  AtlasAlphasType  m_AlphasInVertex3;

  AtlasMesh::ConstPointer  m_Mesh;

  AtlasAlphasType  m_XGradientBasis;
  AtlasAlphasType  m_YGradientBasis;
  AtlasAlphasType  m_ZGradientBasis;

  double  m_MinLogLikelihoodTimesPrior;

  //
  GradientType  m_Gradient;
  HessianType  m_Hessian;

  TransformType::ConstPointer  m_MeshToImageTransform;

  // Lookup table mapping point indices into consecutive indices
  std::map< AtlasMesh::PointIdentifier, int >  m_ComponentXLookupTable;
  std::map< AtlasMesh::PointIdentifier, int >  m_ComponentYLookupTable;
  std::map< AtlasMesh::PointIdentifier, int >  m_ComponentZLookupTable;

  int  m_IndicesForVertex0[ 3 ];
  int  m_IndicesForVertex1[ 3 ];
  int  m_IndicesForVertex2[ 3 ];
  int  m_IndicesForVertex3[ 3 ];

  // Temporary storage for holding the contribution of one tetrahedron to the
  // Hessian and gradient:
  //   * since the prior term contributions are symmetric 12 x 12 matrices
  //     (12 because of 4 nodes per tetrahedron, and 3 coordinates per node), we
  //     only calculate the upper triangle of this matrix. This results in
  //     ( 12 * 13 ) / 2 = 78 unique entries.
  //   * since the data term contributions are not only symmetric 12 x 12 matrices,
  //     but every 3x3 submatrix (each corresponding to one pair of nodes) is also
  //     symmetric within itself, we only calculate the ( 3 * 4 ) / 2 = 6 unique
  //     entries for every ( 4 * 5 ) / 2 = 10 unique pairs, resulting in 60 unique
  //     entries.
  // Upon completion of each tetrahedron, these contributions are flushed out into
  // the correct entries in the global Hessian and gradient. In our implementation,
  // we don't actually have a function that's called upon completion of a tetrahedron.
  // We therefore cheat a little and flush out stuff when a new tetrahedron starts.
  // Since the last tetrahedron doesn't get flushed out this way, it is done when
  // people ask for the results
  double  m_HessianPriorContribution[ 78 ]; // Upper diagonal of a 12 x 12 matrix (78 = 12 * 13 / 2 )
  double  m_HessianDataContribution[ 60 ];
  double  m_GradientContribution[ 12 ];
  
  bool  m_ContributionsFromPreviousTetrahdronFlushedOut;
  bool  m_ComputeHessian;

};

 
 

} // End namespace FragmentProcessor





/**
 *
 */
class AtlasMeshDeformationLevenbergMarquardt :
  public AtlasMeshRasterizor< FragmentProcessor::CalculateGradientAndHessian >
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationLevenbergMarquardt  Self;
  typedef AtlasMeshRasterizor< FragmentProcessor::CalculateGradientAndHessian >  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationLevenbergMarquardt, itk::Object );

  /** Some typedefs */
  typedef Superclass::FragmentProcessorType  FragmentProcessorType;
  typedef Superclass::LabelImageType  LabelImageType;
  typedef FragmentProcessorType::ProbabilityImageType  ProbabilityImageType;
  typedef FragmentProcessorType::SegmentedImageType  SegmentedImageType;
  typedef FragmentProcessorType::ImageType  ImageType;
  typedef FragmentProcessorType::TransformType  TransformType;
  typedef FragmentProcessorType::HessianType  HessianType;
  typedef FragmentProcessorType::GradientType  GradientType;

  
  /** */
  void  SetMaximumNumberOfCGIterations( int maximumNumberOfCGIterations )
    { m_MaximumNumberOfCGIterations = maximumNumberOfCGIterations; }
    
  int GetMaximumNumberOfCGIterations() const
    { return m_MaximumNumberOfCGIterations; }

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
  std::vector<unsigned char > *  GetMapCompToComp() 
    {
    return this->GetFragmentProcessor().GetMapCompToComp();
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
  void  SetImages( std::vector<itk::Image< float, 3 >::Pointer> images ) 
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetImages(images );
      }
    }

  /** */
  const ImageType::Pointer  GetImage(int ind) const
    {
    return this->GetFragmentProcessor().GetImage(ind);
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
  void  SetMeshToImageTransform( const TransformType* meshToImageTransform )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetMeshToImageTransform( meshToImageTransform );
      }
    }

  const TransformType* GetMeshToImageTransform() const
    {
    return this->GetFragmentProcessor().GetMeshToImageTransform();
    }

  /** */
  double GetMinLogLikelihoodTimesPrior() const
    {
    double  minLogLikelihoodTimesPrior = 0;
    for ( std::vector< FragmentProcessorType >::const_iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      if ( it->GetMinLogLikelihoodTimesPrior() == itk::NumericTraits< double >::max() )
        {
        return itk::NumericTraits< double >::max();  
        }
  
      minLogLikelihoodTimesPrior += it->GetMinLogLikelihoodTimesPrior();
      }
      
      return minLogLikelihoodTimesPrior;
    }
  
  /** */
  AtlasPositionGradientContainerType::Pointer  GetStep( float lambda, bool verbose=false ) const;
  
  /** */
  void  SetPeerToCopyHessianFrom( const Self*  peerToCopyHessianFrom )
    {
    m_PeerToCopyHessianFrom = peerToCopyHessianFrom;
    
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      if ( !m_PeerToCopyHessianFrom )
        {
        it->SetComputeHessian( true );  
        }
      else
        {
        it->SetComputeHessian( false );
        }
      }
    
    }
    
  /** */
  const Self*  GetPeerToCopyHessianFrom() const
    {
    return m_PeerToCopyHessianFrom;
    }
    
  
  /**  */
  void Rasterize( const AtlasMesh* mesh );

protected:
  AtlasMeshDeformationLevenbergMarquardt();
  virtual ~AtlasMeshDeformationLevenbergMarquardt();

private:
  AtlasMeshDeformationLevenbergMarquardt(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  
  /** Internal structure used for passing information to the threading library */
  struct ThreadStruct
    {
    std::vector< FragmentProcessorType >*  m_FragmentProcessors;
    std::set< gmm::size_type >  m_RowNumbersToAdd;
    };

  /** Static function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to ThreadedGenerateData(). */
  static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );
  
  //
  static bool AddHessianRows( std::vector< FragmentProcessorType >*  fragmentProcessors, 
                              std::set< gmm::size_type >&  rowNumbersToAdd, 
                              int numberOfRows, int threadId );
  
  mutable double  m_OldLambda;
  int  m_MaximumNumberOfCGIterations;
  
  Self::ConstPointer  m_PeerToCopyHessianFrom;
  
};


} // end namespace kvl

#endif

