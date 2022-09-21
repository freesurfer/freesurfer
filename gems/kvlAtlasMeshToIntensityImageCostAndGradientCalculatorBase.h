#ifndef KVLATLASMESHTOINTENSITYIMAGECOSTANDGRADIENTCALCULATORBASE_H
#define KVLATLASMESHTOINTENSITYIMAGECOSTANDGRADIENTCALCULATORBASE_H

#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "itkImage.h"
#include "kvlLikelihoodImageFilterBase.h"


namespace kvl
{


class AtlasMeshToIntensityImageCostAndGradientCalculatorBase:
        public AtlasMeshPositionCostAndGradientCalculator
{
public:

    /** Standard class typedefs */
    typedef AtlasMeshToIntensityImageCostAndGradientCalculatorBase  Self;
    typedef AtlasMeshPositionCostAndGradientCalculator Superclass;
    typedef itk::SmartPointer< Self >  Pointer;
    typedef itk::SmartPointer< const Self >  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( AtlasMeshToIntensityImageCostAndGradientCalculator, AtlasMeshPositionCostAndGradientCalculator );

    /** Some typedefs */
    typedef itk::Image< float, 3 >  ImageType;

    /** */
    void SetImages( const std::vector< ImageType::ConstPointer >& images );

    /** */
    void Rasterize( const AtlasMesh* mesh );

protected:
    AtlasMeshToIntensityImageCostAndGradientCalculatorBase();
    virtual ~AtlasMeshToIntensityImageCostAndGradientCalculatorBase();

    void AddDataContributionOfTetrahedron( const AtlasMesh::PointType& p0,
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
                                           AtlasPositionGradientThreadAccumType&  gradientInVertex3 ) override;

    //
    typedef LikelihoodImageFilterBase< ImageType >  LikelihoodFilterType;
    LikelihoodFilterType::Pointer  m_LikelihoodFilter;

private:
    AtlasMeshToIntensityImageCostAndGradientCalculatorBase(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented


};


} // end namespace kvl

#endif // KVLATLASMESHTOINTENSITYIMAGECOSTANDGRADIENTCALCULATORBASE_H
