#ifndef __kvlAtlasMeshDeformationLevenbergMarquardt_h
#define __kvlAtlasMeshDeformationLevenbergMarquardt_h

#include "kvlAtlasMeshRasterizor.h"
#include "itkAffineTransform.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkTimeProbe.h"


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
  typedef std::map< AtlasMesh::PointIdentifier, itk::FixedArray< double, 3 > >   GradientContainerType;
  typedef std::map< AtlasMesh::PointIdentifier, itk::FixedArray< double, 3 > >   DiagonalHessianContainerType;
  typedef itk::AffineTransform< float, 3 >  TransformType;
  typedef itk::Image< float, 3 >  ImageType;
  

  CalculateGradientAndHessian()
    {
    m_ProbabilityImage = 0;
    m_Image = 0;
    m_Mesh = 0;
    m_MinLogLikelihoodTimesPrior = 0;

    m_MeshToImageTransform = 0;

    m_GradientInVertex0 = 0;
    m_GradientInVertex1 = 0;
    m_GradientInVertex2 = 0;
    m_GradientInVertex3 = 0;

    m_HessianInVertex0 = 0;
    m_HessianInVertex1 = 0;
    m_HessianInVertex2 = 0;
    m_HessianInVertex3 = 0;

    m_UseProbabilityImage = false;
    }

  ~CalculateGradientAndHessian() {}

  void SetImage( const ImageType* image ) 
    {
    m_Image = image;
    }

  const ImageType* GetImage() const
    { return m_Image; }

  void SetProbabilityImage( const ProbabilityImageType* probabilityImage )
    {
    m_ProbabilityImage = probabilityImage;
    }

  const ProbabilityImageType* GetProbabilityImage() const
    { return m_ProbabilityImage; }

  void  SetUseProbabilityImage( bool  useProbabilityImage )
    {
    m_UseProbabilityImage = useProbabilityImage;
    }

  bool  GetUseProbabilityImage() const
    {
    return m_UseProbabilityImage;
    }

  const AtlasMesh*  GetMesh() const
    { return m_Mesh; }

  void SetMeans( const itk::Array< float >& means )
    { m_Means = means; }

  //
  void SetVariances( const itk::Array< float >& variances )
    { m_Variances = variances; }

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

    if ( m_UseProbabilityImage )
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

      const ImageType::PixelType  intensity = m_Image->GetPixel( m_Index );

      // Skip zero entries
      if ( intensity == 0 )
        {
        // Move on to the next pixel
        m_Index[ 0 ]++;

        return;
        }

      double likelihood = 0.0f;
      double  intensityXGradientBasis = 0.0f;
      double  intensityYGradientBasis = 0.0f;
      double  intensityZGradientBasis = 0.0f;
      for ( unsigned int classNumber = 0; classNumber < m_AlphasInVertex0.Size(); classNumber++ )
        {
        // Evaluate the Gaussian of this class at the intensity of this pixel
        const double  gauss = exp( -pow( intensity - m_Means[ classNumber ] , 2 ) /
                                    2 / m_Variances[ classNumber ] ) / 
                                    sqrt( 2 * 3.14 * m_Variances[ classNumber ] );

        // Collect the data terms for each vertex
        double  alpha0 = m_AlphasInVertex0[ classNumber ];
        double  alpha1 = m_AlphasInVertex1[ classNumber ];
        double  alpha2 = m_AlphasInVertex2[ classNumber ];
        double  alpha3 = m_AlphasInVertex3[ classNumber ];

        // Add contribution of the likelihood
        likelihood += gauss * ( alpha0 * pi0 + alpha1 * pi1 + alpha2 * pi2 + alpha3 * pi3 ) + 1e-5;

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



    // Add data contribution to gradient
    ( *m_GradientInVertex0 )[0] += tmpX * pi0;
    ( *m_GradientInVertex0 )[1] += tmpY * pi0;
    ( *m_GradientInVertex0 )[2] += tmpZ * pi0;

    ( *m_GradientInVertex1 )[0] += tmpX * pi1;
    ( *m_GradientInVertex1 )[1] += tmpY * pi1;
    ( *m_GradientInVertex1 )[2] += tmpZ * pi1;
    
    ( *m_GradientInVertex2 )[0] += tmpX * pi2;
    ( *m_GradientInVertex2 )[1] += tmpY * pi2;
    ( *m_GradientInVertex2 )[2] += tmpZ * pi2;

    ( *m_GradientInVertex3 )[0] += tmpX * pi3;
    ( *m_GradientInVertex3 )[1] += tmpY * pi3;
    ( *m_GradientInVertex3 )[2] += tmpZ * pi3;

    // Add data contribution to Hessian
    ( *m_HessianInVertex0 )[0] += tmpXX * pow( pi0, 2 );
    ( *m_HessianInVertex0 )[1] += tmpYY * pow( pi0, 2 );
    ( *m_HessianInVertex0 )[2] += tmpZZ * pow( pi0, 2 );

    ( *m_HessianInVertex1 )[0] += tmpXX * pow( pi1, 2 );
    ( *m_HessianInVertex1 )[1] += tmpYY * pow( pi1, 2 );
    ( *m_HessianInVertex1 )[2] += tmpZZ * pow( pi1, 2 );
    
    ( *m_HessianInVertex2 )[0] += tmpXX * pow( pi2, 2 );
    ( *m_HessianInVertex2 )[1] += tmpYY * pow( pi2, 2 );
    ( *m_HessianInVertex2 )[2] += tmpZZ * pow( pi2, 2 );

    ( *m_HessianInVertex3 )[0] += tmpXX * pow( pi3, 2 );
    ( *m_HessianInVertex3 )[1] += tmpYY * pow( pi3, 2 );
    ( *m_HessianInVertex3 )[2] += tmpZZ * pow( pi3, 2 );


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

    //
    // Cache relevant elements of the vertices of this tetrahedron. The notation used is Y = [ p0 p1 p2 p3; 1 1 1 1 ]
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
    m_GradientInVertex0 = &( m_Gradient[ *pit ] );
    m_HessianInVertex0 = &( m_Hessian[ *pit ] );
    
    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    const double y12 = p[ 0 ];
    const double y22 = p[ 1 ];
    const double y32 = p[ 2 ];
    m_AlphasInVertex1 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_GradientInVertex1 = &( m_Gradient[ *pit ] );
    m_HessianInVertex1 = &( m_Hessian[ *pit ] );

    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    const double y13 = p[ 0 ];
    const double y23 = p[ 1 ];
    const double y33 = p[ 2 ];
    m_AlphasInVertex2 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_GradientInVertex2 = &( m_Gradient[ *pit ] );
    m_HessianInVertex2 = &( m_Hessian[ *pit ] );

    ++pit;
    m_Mesh->GetPoint( *pit, &p );
    const double y14 = p[ 0 ];
    const double y24 = p[ 1 ];
    const double y34 = p[ 2 ];
    m_AlphasInVertex3 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_GradientInVertex3 = &( m_Gradient[ *pit ] );
    m_HessianInVertex3 = &( m_Hessian[ *pit ] );



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
      std::cout << "Oooooops: tetrahedron " << cellId << " has managed to turn bad (detJ: " << detJ << ")" << std::endl;

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
    ( *m_GradientInVertex0 )[0] += dcostdy11;
    ( *m_GradientInVertex0 )[1] += dcostdy21;
    ( *m_GradientInVertex0 )[2] += dcostdy31;
    
    ( *m_GradientInVertex1 )[0] += dcostdy12;
    ( *m_GradientInVertex1 )[1] += dcostdy22;
    ( *m_GradientInVertex1 )[2] += dcostdy32;
    
    ( *m_GradientInVertex2 )[0] += dcostdy13;
    ( *m_GradientInVertex2 )[1] += dcostdy23;
    ( *m_GradientInVertex2 )[2] += dcostdy33;
    
    ( *m_GradientInVertex3 )[0] += dcostdy14;
    ( *m_GradientInVertex3 )[1] += dcostdy24;
    ( *m_GradientInVertex3 )[2] += dcostdy34;

    const double  tmp = 1.0f / ( 2 * priorCost + 1e-15 );

    ( *m_HessianInVertex0 )[0] += tmp * pow( dcostdy11, 2 );
    ( *m_HessianInVertex0 )[1] += tmp * pow( dcostdy21, 2 );
    ( *m_HessianInVertex0 )[2] += tmp * pow( dcostdy31, 2 );
    
    ( *m_HessianInVertex1 )[0] += tmp * pow( dcostdy12, 2 );
    ( *m_HessianInVertex1 )[1] += tmp * pow( dcostdy22, 2 );
    ( *m_HessianInVertex1 )[2] += tmp * pow( dcostdy32, 2 );
    
    ( *m_HessianInVertex2 )[0] += tmp * pow( dcostdy13, 2 );
    ( *m_HessianInVertex2 )[1] += tmp * pow( dcostdy23, 2 );
    ( *m_HessianInVertex2 )[2] += tmp * pow( dcostdy33, 2 );
    
    ( *m_HessianInVertex3 )[0] += tmp * pow( dcostdy14, 2 );
    ( *m_HessianInVertex3 )[1] += tmp * pow( dcostdy24, 2 );
    ( *m_HessianInVertex3 )[2] += tmp * pow( dcostdy34, 2 );
    
    

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

    //std::cout << "Done setting up tetrahedron" << std::endl;

    return true;
    }

  inline void SetMesh( const AtlasMesh* mesh )
    {
    m_MinLogLikelihoodTimesPrior = 0.0f;

    m_Mesh = mesh;

    // Create containers for the gradient and (diagonal) hessian in each vertex, zero-filled each
    m_Gradient.clear();
    m_Hessian.clear();
    itk::FixedArray< double, 3 >  zeroEntry;
    zeroEntry.Fill( 0.0 );
    for ( AtlasMesh::PointsContainer::ConstIterator  pointIt = m_Mesh->GetPoints()->Begin();
          pointIt != m_Mesh->GetPoints()->End(); ++pointIt )
      {
      m_Gradient[ pointIt.Index() ] = zeroEntry;
      m_Hessian[ pointIt.Index() ] = zeroEntry;
      }

    //std::cout << "Done setting mesh" << std::endl;
    }

  void  SetMeshToImageTransform( const TransformType* meshToImageTransform )
    { m_MeshToImageTransform = meshToImageTransform; }

  const TransformType* GetMeshToImageTransform() const
    { return m_MeshToImageTransform; }

  const GradientContainerType& GetGradient() const
    {
    return m_Gradient;
    }

  const DiagonalHessianContainerType& GetHessian() const
    {
    return m_Hessian;
    }

  double GetMinLogLikelihoodTimesPrior() const
    {
    //std::cout << "PositionGradientCalculator: returning m_MinLogLikelihoodTimesPrior " << m_MinLogLikelihoodTimesPrior << std::endl;
    return m_MinLogLikelihoodTimesPrior;
    }


private:


  ProbabilityImageType::ConstPointer  m_ProbabilityImage;
  ImageType::ConstPointer  m_Image;
  ImageType::IndexType  m_Index;

  bool m_UseProbabilityImage;

  itk::Array< float >  m_Means;
  itk::Array< float >  m_Variances;

  AtlasAlphasType  m_AlphasInVertex0;
  AtlasAlphasType  m_AlphasInVertex1;
  AtlasAlphasType  m_AlphasInVertex2;
  AtlasAlphasType  m_AlphasInVertex3;
  
  itk::FixedArray< double, 3 >*  m_GradientInVertex0;
  itk::FixedArray< double, 3 >*  m_GradientInVertex1;
  itk::FixedArray< double, 3 >*  m_GradientInVertex2;
  itk::FixedArray< double, 3 >*  m_GradientInVertex3;

  itk::FixedArray< double, 3 >*  m_HessianInVertex0;
  itk::FixedArray< double, 3 >*  m_HessianInVertex1;
  itk::FixedArray< double, 3 >*  m_HessianInVertex2;
  itk::FixedArray< double, 3 >*  m_HessianInVertex3;

  
  AtlasMesh::ConstPointer  m_Mesh;
  GradientContainerType  m_Gradient;
  DiagonalHessianContainerType  m_Hessian;
  
  AtlasAlphasType  m_XGradientBasis;
  AtlasAlphasType  m_YGradientBasis;
  AtlasAlphasType  m_ZGradientBasis;

  double  m_MinLogLikelihoodTimesPrior;

  TransformType::ConstPointer  m_MeshToImageTransform;


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
  typedef FragmentProcessorType::ImageType  ImageType;
  typedef FragmentProcessorType::TransformType  TransformType;

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
  void  SetImage( const ImageType* image )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetImage( image );
      }
    }

  /** */
  const ImageType*  GetImage() const
    {
    return this->GetFragmentProcessor().GetImage();
    }

  /** */
  void SetMeans( const itk::Array< float >& means )
    { 
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetMeans( means ); 
      }
    }

  /** */ 
  void SetVariances( const itk::Array< float >& variances )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetVariances( variances ); 
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
  void  SetUseProbabilityImage( bool  useProbabilityImage )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetUseProbabilityImage( useProbabilityImage );
      }
    }

  /** */
  bool  GetUseProbabilityImage() const
    {
    return this->GetFragmentProcessor().GetUseProbabilityImage();
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
  

  /**  */
  void Rasterize( const AtlasMesh* mesh )
    {
    // Rasterize
    itk::TimeProbe  timeProbe;
    timeProbe.Start();
    Superclass::Rasterize( mesh, true );
    timeProbe.Stop();
    std::cout << "Time taken to rasterize Levenberg-Marquardt: " << timeProbe.GetMeanTime() << std::endl;
    
    // Collect results of all threads, by adding the contributions from the second-and-up threads to the
    // default one (i.e., the first thread)
    timeProbe = itk::TimeProbe();
    timeProbe.Start();
    m_Gradient = this->GetFragmentProcessor().GetGradient();
    m_Hessian = this->GetFragmentProcessor().GetHessian();

    std::vector< FragmentProcessorType >::const_iterator it = this->GetFragmentProcessors().begin();
    ++it; // Go to the second thread
    for ( ; it != this->GetFragmentProcessors().end(); ++it )
      {
      // Add gradient of this fragment processor
      FragmentProcessorType::GradientContainerType::const_iterator  sourceGradIt = it->GetGradient().begin();
      FragmentProcessorType::GradientContainerType::iterator  targetGradIt = m_Gradient.begin();
      for ( ; targetGradIt != m_Gradient.end(); ++targetGradIt, ++sourceGradIt )
        {
        for ( int i = 0; i < 3; i++ )
          {
          ( targetGradIt->second )[i] += ( sourceGradIt->second )[i];
          }
        }  

      // Add hessian of this fragment processor
      FragmentProcessorType::DiagonalHessianContainerType::const_iterator  sourceHessIt = it->GetHessian().begin();
      FragmentProcessorType::DiagonalHessianContainerType::iterator  targetHessIt = m_Hessian.begin();
      for ( ; targetHessIt != m_Hessian.end(); ++targetHessIt, ++sourceHessIt )
        {
        for ( int i = 0; i < 3; i++ )
          {
          ( targetHessIt->second )[i] += ( sourceHessIt->second )[i];
          }
        }  
      

      } // End loop over all fragment processors
    
    timeProbe.Stop();
    std::cout << "Time taken to add contributions: " << timeProbe.GetMeanTime() << std::endl;

    
    }

protected:
  AtlasMeshDeformationLevenbergMarquardt();
  virtual ~AtlasMeshDeformationLevenbergMarquardt();

private:
  AtlasMeshDeformationLevenbergMarquardt(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  FragmentProcessorType::GradientContainerType  m_Gradient;
  FragmentProcessorType::DiagonalHessianContainerType  m_Hessian;
  
  
};


} // end namespace kvl

#endif

