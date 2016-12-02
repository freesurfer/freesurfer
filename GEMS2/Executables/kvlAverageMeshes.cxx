#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMesh.h"
#include "kvlAtlasMeshAveragingConjugateGradientOptimizer.h"  // I just copied kvlAtlasMeshDeformationConjugateGradientOptimizer, and will modify it as I need to loop over meshes, etc..

int main( int argc, char** argv )   
{
  // Sanity check on input
  if ( argc < 6 )
    {
    std::cerr << "Usage: "<< argv[ 0 ] << " inputMeshCollection initialMeshCollection Katlas KtimePoints outputMeshCollection [numberOfMeshesToAverage]" << std::endl;
    std::cerr << "   initialMeshCollection is expected to have initial position as position 0 " << std::endl;
    return -1;
    }


  // Read mesh collections from file
  kvl::AtlasMeshCollection::Pointer  inputCollection =  kvl::AtlasMeshCollection::New();
  inputCollection->Read( argv[ 1 ] );


  kvl::AtlasMeshCollection::Pointer  initialCollection =  kvl::AtlasMeshCollection::New();
  initialCollection->Read( argv[ 2 ] );

  std::istringstream  KAstream( argv[ 3 ] );
  float  Katlas;
  KAstream >> Katlas;

  std::istringstream  KTPstream( argv[ 4 ] );
  float  KtimePoints;
  KTPstream >> KtimePoints;

  if ( argc > 6 )
  {
    int  NmeshesToAverage;
    std::istringstream  MTAstream( argv[ 6 ] );
    MTAstream >> NmeshesToAverage;
    inputCollection->ResizePositions(NmeshesToAverage);
  }

  kvl::AtlasMeshAveragingConjugateGradientOptimizer::Pointer optimizer =   kvl::AtlasMeshAveragingConjugateGradientOptimizer::New();
  optimizer->SetInitialMeshCollection(initialCollection); 
  optimizer->SetInputMeshCollection(inputCollection);
//    optimizer->SetMaximalDeformationStopCriterion( 0.00001 );
//  optimizer->SetMeshToImageTransform( transform );
  optimizer->SetKatlas(Katlas);
  optimizer->SetKtimePoints(KtimePoints);
  optimizer->SetVerbose(true);

//  printf("Ks: %f, %f \n",Katlas,KtimePoints); 

for (int j=0; j<10; j++)
{
  bool ready=false;
  int i=0;
  while (!ready)
  {
      float maximalDeformation = optimizer->PerformOneIteration();
      float minLogLikelihoodTimesPrior = optimizer->GetMinLogLikelihoodTimesPrior();

      printf("Iteration %d-%d: took a maximal step of %.3f and the new cost is now %.4f \n",j+1,i+1,maximalDeformation, minLogLikelihoodTimesPrior);

      i=i+1;
  
      if (maximalDeformation < 5e-4 || i>=20)
      {
        ready=true;
      }
   }

   optimizer -> Initialize();
   optimizer->SetInitialMeshCollection(initialCollection);
/*
   optimizer =   kvl::AtlasMeshAveragingConjugateGradientOptimizer::New();
   optimizer->SetInitialMeshCollection(initialCollection); 
   optimizer->SetInputMeshCollection(inputCollection);
//    optimizer->SetMaximalDeformationStopCriterion( 0.00001 );
//  optimizer->SetMeshToImageTransform( transform );
   optimizer->SetKatlas(Katlas);
   optimizer->SetKtimePoints(KtimePoints);
   optimizer->SetVerbose(true);
*/
}


  // Write the result out
  if ( !initialCollection->Write( argv[5] ) )
    {
    std::cerr << "Could not write mesh collection to " << argv[5] << ".gz" << std::endl;
    exit( -1 );
    }
  std::cout << "Just wrote mesh collection to " << argv[5] << ".gz" << std::endl;


  return 0;
};
