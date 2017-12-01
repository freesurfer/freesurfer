#ifndef GEMS_PYKVLOPTIMIZER_H
#define GEMS_PYKVLOPTIMIZER_H

#include "kvlAtlasMeshCollection.h"
#include "pyKvlImage.h"
#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"
#include "kvlConditionalGaussianEntropyCostAndGradientCalculator.h"
#include "kvlMutualInformationCostAndGradientCalculator.h"

#include "kvlAtlasMeshDeformationFixedStepGradientDescentOptimizer.h"
#include "kvlAtlasMeshDeformationGradientDescentOptimizer.h"
#include "kvlAtlasMeshDeformationConjugateGradientOptimizer.h"
#include "kvlAtlasMeshDeformationLBFGSOptimizer.h"

#include "pyKvlNumpy.h"
#include "pybind11/pybind11.h"
#include "pyKvlMesh.h"

namespace py = pybind11;

typedef itk::Image< float, 3 >  ImageType;
typedef ImageType::Pointer ImagePointer;

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

class KvlOptimizer {
    kvl::AtlasMeshDeformationOptimizer::Pointer optimizer;
public:
    KvlOptimizer(std::string typeName,
                               const KvlMesh &mesh,
                               const KvlCostAndGradientCalculator &calculator,
                               std::map<std::string, double> arguments){
        switch( typeName[ 0 ] )
        {
//            case 'F':
//            {
//                std::cout << "FixedStepGradientDescent" << std::endl;
//                AtlasMeshDeformationFixedStepGradientDescentOptimizer::Pointer  myOptimizer
//                        = AtlasMeshDeformationFixedStepGradientDescentOptimizer::New();
//                myOptimizer->SetStepSize( 1.0 );
//                optimizer = myOptimizer;
//                break;
//            }
//            case 'G':
//            {
//                std::cout << "GradientDescent" << std::endl;
//                AtlasMeshDeformationGradientDescentOptimizer::Pointer  myOptimizer
//                        = AtlasMeshDeformationGradientDescentOptimizer::New();
//                optimizer = myOptimizer;
//                break;
//            }
//            case 'C':
//            {
//                std::cout << "ConjugateGradient" << std::endl;
//                AtlasMeshDeformationConjugateGradientOptimizer::Pointer  myOptimizer
//                        = AtlasMeshDeformationConjugateGradientOptimizer::New();
//                optimizer = myOptimizer;
//                break;
//            }
            case 'L':
            {
                std::cout << "L-BFGS" << std::endl;
                kvl::AtlasMeshDeformationLBFGSOptimizer::Pointer  myOptimizer
                        = kvl::AtlasMeshDeformationLBFGSOptimizer::New();
                optimizer = myOptimizer;
                break;
            }
            default:
            {
                throw std::invalid_argument( "Optimizer type not supported." );
            }
        }


        // Parse additional options. Format is always [ 'someString', double ]
        for(auto const& arg : arguments){

            const std::string optionName = arg.first;
            const double optionValue = arg.second;

            switch( optionName[ 0 ] )
            {
                case 'V':
                {
                    std::cout << "Verbose: " << optionValue << std::endl;
                    if ( optionValue )
                    {
                        optimizer->SetVerbose( true );
                    }

                    break;
                }
                case 'M':
                {
                    if ( optionName.substr( 0, 8 ) == "MaximalD" )
                    {
                        std::cout << "MaximalDeformationStopCriterion: " << optionValue << std::endl;
                        optimizer->SetMaximalDeformationStopCriterion( optionValue );
                    }
                    else if ( optionName.substr( 0, 8 ) == "MaximumN" )
                    {
                        const int  maximumNumberOfIterations = static_cast< int >( optionValue );
                        std::cout << "MaximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
                        optimizer->SetMaximumNumberOfIterations( maximumNumberOfIterations );
                    }
                    else
                    {
                        std::ostringstream  errorStream;
                        errorStream << "optionName: " << optionName << " not understood";
                        throw std::invalid_argument( errorStream.str().c_str() );
                        break;
                    } // End figuring out which "M" you mean

                    break;
                }
                case 'L':
                {
                    std::cout << "LineSearchMaximalDeformationIntervalStopCriterion: " << optionValue << std::endl;
                    optimizer->SetLineSearchMaximalDeformationIntervalStopCriterion( optionValue );
                    break;
                }
                case 'B':
                {
                    kvl::AtlasMeshDeformationLBFGSOptimizer::Pointer  myOptimizer
                            = dynamic_cast< kvl::AtlasMeshDeformationLBFGSOptimizer* >( optimizer.GetPointer() );
                    if ( myOptimizer )
                    {
                        const int  maximumMemoryLength = static_cast< int >( optionValue );
                        std::cout << "BFGS-MaximumMemoryLength: " << maximumMemoryLength << std::endl;
                        myOptimizer->SetMaximumMemoryLength( maximumMemoryLength );
                    }
                    else
                    {
                        std::cout << "BFGS-MaximumMemoryLength only applies to BFGS optimizer" << std::endl;
                    }

                    break;
                }
                default:
                {
                    std::ostringstream  errorStream;
                    errorStream << "optionName: " << optionName << " not understood";
                    throw std::invalid_argument( errorStream.str().c_str() );
                }
            }

        }
        kvl::AtlasMesh::ConstPointer constMesh = mesh.mesh;
        kvl::AtlasMesh::Pointer mutableMesh = const_cast< kvl::AtlasMesh* >( constMesh.GetPointer() );
        // Pass the mesh and calculator to it
        optimizer->SetMesh( mutableMesh );
        optimizer->SetCostAndGradientCalculator( calculator.calculator );
    }
    std::pair<double, double> StepOptimizer(){
        const double  maximalDeformation = optimizer->Step();
        const double  minLogLikelihoodTimesPrior = optimizer->GetMinLogLikelihoodTimesPrior();
        return {minLogLikelihoodTimesPrior, maximalDeformation};
    }
};
#endif //GEMS_PYKVLOPTIMIZER_H
