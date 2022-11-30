#ifndef GEMS_PYKVLCALCULATOR_H
#define GEMS_PYKVLCALCULATOR_H

#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"
#include "kvlAtlasMeshToIntensityImageLogDomainCostAndGradientCalculator.h"
#include "kvlConditionalGaussianEntropyCostAndGradientCalculator.h"
#include "kvlMutualInformationCostAndGradientCalculator.h"
#include "kvlAtlasMeshToPointSetCostAndGradientCalculator.h"
#include "kvlAverageAtlasMeshPositionCostAndGradientCalculator.h"

#include "kvlAtlasMeshCollection.h"
#include "pyKvlImage.h"
#include "pyKvlNumpy.h"
#include "pybind11/pybind11.h"
#include "pyKvlMesh.h"

class KvlCostAndGradientCalculator {

public:
    kvl::AtlasMeshPositionCostAndGradientCalculator::Pointer calculator;

    KvlCostAndGradientCalculator(std::string typeName,
                                 std::vector<KvlImage> images,
                                 std::string boundaryCondition,
                                 KvlTransform transform=KvlTransform(nullptr),
                                 py::array_t<double> means=py::array_t<double>(),
                                 py::array_t<double> variances=py::array_t<double>(),
                                 py::array_t<float> mixtureWeights=py::array_t<float>(),
                                 py::array_t<int> numberOfGaussiansPerClass=py::array_t<int>(),
                                 py::array_t<double> targetPoints=py::array_t<double>()
    ){
        if (typeName == "AtlasMeshToIntensityImage" || typeName == "AtlasMeshToIntensityImageLogDomain") {

            py::buffer_info means_info  = means.request();

            // Retrieve means if they are provided
            const int  numberOfGaussians = means_info.shape[0];
            const int  numberOfContrasts  = means_info.shape[1];

            std::vector< vnl_vector< double > >  means_converted;

            for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ ) {
                vnl_vector< double >  mean_converted( numberOfContrasts, 0.0f );

                for ( int contrastNumber = 0; contrastNumber < numberOfContrasts; contrastNumber++ ) {
                    mean_converted[ contrastNumber ] = means.at(gaussianNumber, contrastNumber);
                }
                means_converted.push_back( mean_converted );
            }
            // Retrieve variances if they are provided
            std::vector< vnl_matrix< double > >  variances_converted;
            for ( unsigned int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ ) {
                vnl_matrix< double >  variance( numberOfContrasts, numberOfContrasts, 0.0f );
                for ( unsigned int row = 0; row < numberOfContrasts; row++ ) {
                  for ( unsigned int col = 0; col < numberOfContrasts; col++ ) {
                    variance[ row ][ col ] = variances.at(gaussianNumber, row, col);
                    }
                  }
                  variances_converted.push_back( variance );
            }
            // Retrieve mixtureWeights if they are provided
            std::vector< double >  mixtureWeights_converted  = std::vector< double >( numberOfGaussians, 0.0f );
            for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ ) {
                mixtureWeights_converted[gaussianNumber] = ( mixtureWeights.at(gaussianNumber));
            }
            // Retrieve numberOfGaussiansPerClass if they are provided
            const int  numberOfClasses = numberOfGaussiansPerClass.request().shape[0];
            std::vector< int >  numberOfGaussiansPerClass_converted = std::vector< int >( numberOfClasses, 0 );
            for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ ) {
                numberOfGaussiansPerClass_converted[ classNumber ] = numberOfGaussiansPerClass.at(classNumber);
            }

            kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::Pointer myCalculator;
	    if (typeName == "AtlasMeshToIntensityImage")
            {
              myCalculator = kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::New();
            }
            else
	    {
              myCalculator = kvl::AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator::New();
            }

            std::vector< ImageType::ConstPointer> images_converted;
            for(auto image: images){
                ImageType::ConstPointer constImage = static_cast< const ImageType* >( image.m_image.GetPointer() );
                images_converted.push_back( constImage );
            }
            myCalculator->SetImages( images_converted );
            myCalculator->SetParameters( means_converted, variances_converted, mixtureWeights_converted, numberOfGaussiansPerClass_converted );
            calculator = myCalculator;

        } else if (typeName == "MutualInformation") {

            kvl::MutualInformationCostAndGradientCalculator::Pointer myCalculator = kvl::MutualInformationCostAndGradientCalculator::New();
            myCalculator->SetImage( images[ 0 ].m_image );
            calculator = myCalculator;

        } else if (typeName == "PointSet") {

            kvl::AtlasMeshToPointSetCostAndGradientCalculator::Pointer myCalculator = kvl::AtlasMeshToPointSetCostAndGradientCalculator::New();
            kvl::AtlasMesh::PointsContainer::Pointer alphaTargetPoints = kvl::AtlasMesh::PointsContainer::New();
            CreatePointSetFromNumpy(alphaTargetPoints, targetPoints);
            myCalculator->SetTargetPoints( alphaTargetPoints );
            calculator = myCalculator;

        } else {
            throw std::invalid_argument("Calculator type not supported.");
        }

        // Specify the correct type of boundary condition
        switch( boundaryCondition[ 0 ] )
        {
            case 'S':
            {
                std::cout << "SLIDING" << std::endl;
                calculator->SetBoundaryCondition( kvl::AtlasMeshPositionCostAndGradientCalculator::SLIDING );

                // Retrieve transform if one is provided
                TransformType::ConstPointer  constTransform = 0;
                constTransform = static_cast< const TransformType* >( transform.m_transform.GetPointer() );

                if ( constTransform.GetPointer() )
                {
                    calculator->SetMeshToImageTransform( constTransform );
                }
                break;
            }
            case 'A':
            {
                std::cout << "AFFINE" << std::endl;
                calculator->SetBoundaryCondition( kvl::AtlasMeshPositionCostAndGradientCalculator::AFFINE );
                break;
            }
            case 'T':
            {
                std::cout << "TRANSLATION" << std::endl;
                calculator->SetBoundaryCondition( kvl::AtlasMeshPositionCostAndGradientCalculator::TRANSLATION );
                break;
            }
            case 'N':
            {
                std::cout << "NONE" << std::endl;
                calculator->SetBoundaryCondition( kvl::AtlasMeshPositionCostAndGradientCalculator::NONE );
                break;
            }
            default:
            {
                throw std::invalid_argument( "Boundary condition type not supported." );
            }
        }
    }

    KvlCostAndGradientCalculator(KvlMeshCollection meshCollection,
                                 double K0,
                                 double K1,
                                 KvlTransform transform
    ){
        kvl::AverageAtlasMeshPositionCostAndGradientCalculator::Pointer myCalculator = kvl::AverageAtlasMeshPositionCostAndGradientCalculator::New();
        myCalculator->SetBoundaryCondition(kvl::AtlasMeshPositionCostAndGradientCalculator::SLIDING);
        // Retrieve transform if provided
        TransformType::ConstPointer constTransform = static_cast<const TransformType*>(transform.m_transform.GetPointer());
        if (constTransform.GetPointer()) myCalculator->SetMeshToImageTransform( constTransform );
        // Apply positions and Ks
        kvl::AtlasMeshCollection::Pointer meshCollectionPtr = meshCollection.GetMeshCollection();
        std::vector<double> Ks = {K0};
        std::vector<kvl::AtlasMesh::PointsContainer::ConstPointer> positions = {meshCollectionPtr->GetReferencePosition()};
        for (int  meshNumber = 0; meshNumber < meshCollectionPtr->GetPositions().size(); meshNumber++) {
            positions.push_back(meshCollectionPtr->GetPositions()[meshNumber].GetPointer()); 
            Ks.push_back(K1);
        }
        myCalculator->SetPositionsAndKs(positions, Ks);
        calculator = myCalculator;
    }

    std::pair<double, py::array_t<double>> EvaluateMeshPosition(const KvlMesh &mesh) {
        calculator->Rasterize( mesh.mesh );
        const double cost = calculator->GetMinLogLikelihoodTimesPrior();
        kvl::AtlasPositionGradientContainerType::ConstPointer gradient = calculator->GetPositionGradient();

        const size_t numberOfNodes = gradient->Size();
        auto const data = new double[numberOfNodes*3];
        auto data_it = data;
        for ( kvl::AtlasPositionGradientContainerType::ConstIterator  it = gradient->Begin();
              it != gradient->End(); ++it )
        {
            *data_it++ = it.Value()[0];
            *data_it++ = it.Value()[1];
            *data_it++ = it.Value()[2];
        }
        py::array_t<double> gradient_np = createNumpyArrayCStyle({numberOfNodes, 3}, data);
        return {cost, gradient_np};
    };
};

#endif //GEMS_PYKVLCALCULATOR_H
