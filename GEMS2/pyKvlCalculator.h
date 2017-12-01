#ifndef GEMS_PYKVLCALCULATOR_H
#define GEMS_PYKVLCALCULATOR_H

#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"
#include "kvlConditionalGaussianEntropyCostAndGradientCalculator.h"
#include "kvlMutualInformationCostAndGradientCalculator.h"

#include "kvlAtlasMeshCollection.h"
#include "pyKvlImage.h"
#include "pyKvlNumpy.h"
#include "pybind11/pybind11.h"
#include "pyKvlMesh.h"

class KvlCostAndGradientCalculator {

public:
    kvl::AtlasMeshPositionCostAndGradientCalculator::Pointer calculator;

    KvlCostAndGradientCalculator(std::string typeName, std::vector<KvlImage> images, std::string boundaryCondition){
        switch( typeName[ 0 ] )
        {
            // TODO: implement AtlasMeshToIntensityImage as needed for samseg
//            case 'A':
//            {
//                std::cout << "AtlasMeshToIntensityImage" << std::endl;
//                kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::Pointer  myCalculator
//                        = kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::New();
//                myCalculator->SetImages( images );
//                myCalculator->SetParameters( means, variances, mixtureWeights, numberOfGaussiansPerClass );
//                calculator = myCalculator;
//                break;
//            }
//            case 'C':
//            {
//                std::cout << "ConditionalGaussianEntropy" << std::endl;
//                kvl::ConditionalGaussianEntropyCostAndGradientCalculator::Pointer  myCalculator
//                        = kvl::ConditionalGaussianEntropyCostAndGradientCalculator::New();
//                myCalculator->SetImage( images[ 0 ].imageHandle  );
//                calculator = myCalculator;
//                break;
//            }
            case 'M':
            {
                std::cout << "MutualInformation" << std::endl;
                kvl::MutualInformationCostAndGradientCalculator::Pointer  myCalculator
                        = kvl::MutualInformationCostAndGradientCalculator::New();
                myCalculator->SetImage( images[ 0 ].imageHandle );
                calculator = myCalculator;
                break;
            }
            default:
            {
                throw std::invalid_argument( "Calculator type not supported." );
            }
        }


        // Specify the correct type of boundary condition
        switch( boundaryCondition[ 0 ] )
        {
//            case 'S':
//            {
//                std::cout << "SLIDING" << std::endl;
//                calculator->SetBoundaryCondition( kvl::AtlasMeshPositionCostAndGradientCalculator::SLIDING );
//                if ( constTransform.GetPointer() )
//                {
//                    calculator->SetMeshToImageTransform( constTransform );
//                }
//                break;
//            }
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
    std::pair<double, py::array_t<double>> EvaluateMeshPosition(const KvlMesh &mesh){
        calculator->Rasterize( mesh.mesh );
        const double cost = calculator->GetMinLogLikelihoodTimesPrior();
        kvl::AtlasPositionGradientContainerType::ConstPointer gradient = calculator->GetPositionGradient();

        const size_t numberOfNodes = gradient->Size();
        auto const data = new double[numberOfNodes*3];
        auto data_it = data;
        for ( kvl::AtlasPositionGradientContainerType::ConstIterator  it = gradient->Begin();
              it != gradient->End(); ++it, ++data_it )
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
