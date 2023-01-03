#include "kvlAtlasMeshBuilder.h"
#include "itkCommand.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "kvlCompressionLookupTable.h"
#include "itkMGHImageIOFactory.h"
#include "mri.h"
#include "utils.h"
#include "argparse.h"
#include "fio.h"
#include "samseg-atlas.help.xml.h"

namespace kvl
{

class  BuilderCommand : public itk::Command
{
public:
  /** Standard class typedefs. */
  typedef BuilderCommand         Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self ) ;

  /** Run-time type information (and related methods). */
  itkTypeMacro( BuilderCommand, itk::Command );

  /** Abstract method that defines the action to be taken by the command. */
  virtual void Execute( itk::Object *caller, const itk::EventObject & event )
    {
    this->Execute(  const_cast< const Object* >( caller ), event );
    }

  /** Abstract method that defines the action to be taken by the command.
   * This variant is expected to be used when requests comes from a
   * const Object */
  virtual void Execute(const itk::Object *caller, const itk::EventObject & event )
    {

    if ( ( typeid( event ) == typeid( itk::IterationEvent ) ) ||
         ( typeid( event ) == typeid( itk::EndEvent ) ) )
      {
      kvl::AtlasMeshBuilder::ConstPointer  builder = static_cast< const kvl::AtlasMeshBuilder* >(  caller );

      printf("Writing out the current mesh collection\n");
      PrintMemUsage(stdout);

      // Write out the current mesh collection
      std::ostringstream  fileNameStream;
      fileNameStream << m_LogDirectory << "/CurrentMeshCollection" << builder->GetIterationNumber();
      std::string  fileName = fileNameStream.str();
      if ( !builder->GetCurrentMeshCollection()->Write( fileName.c_str() ) )
        {
        std::cerr <<  "Couldn't write mesh collection to file " << fileName << std::endl;
        exit( -1 );
        }
      std::cout << "Just wrote to " << fileName << std::endl;
    
      // Also append the current minLogLikelihood
      std::string  logFileName = m_LogDirectory;
      logFileName += "/logFile.txt";
      std::ofstream  out( logFileName.c_str(), std::ios::app );
      if ( out.bad() )
        {
        std::cerr <<  "Couldn't write to log file"  << std::endl;
        exit( -1 );
        }

      double  currentDataCost = 0;
      double  currentAlphasCost = 0;
      std::cout << "Computing current data and alphas cost..." << std::endl;
      itk::TimeProbe  probe;
      probe.Start();
      builder->GetCurrentDataAndAlphasCost( currentDataCost, currentAlphasCost );
      probe.Stop();
      std::cout << "...done" << std::endl;
      std::cout << "Took " << probe.GetMean() << " seconds to compute current data and alphas cost" << std::endl;
      double  currentPositionCost = builder->GetCurrentPositionCost();
      const double  currentCost = currentDataCost + currentAlphasCost + currentPositionCost;
      out << builder->GetIterationNumber() << "   "
           << currentDataCost << "   "
           << currentAlphasCost << "   "
           << currentPositionCost << "   "
           << currentCost << std::endl;

      }
    else if ( typeid( event ) == typeid( EdgeAnalysisProgressEvent ) ) 
      {
      // Show some progress
      kvl::AtlasMeshBuilder::ConstPointer  builder = static_cast< const kvl::AtlasMeshBuilder* >(  caller );
      std::cout << "Progress: " << builder->GetProgress() * 100 << "%" << std::endl;
      }

    }

  // 
  void  SetLogDirectory( const std::string&  logDirectory )
    { m_LogDirectory = logDirectory; }

protected:
  BuilderCommand() {}
  virtual ~BuilderCommand() {}

private:
  BuilderCommand(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::string  m_LogDirectory;
};

};



int main( int argc, char** argv )
{
  int  numberOfUpsamplingSteps;
  unsigned int  meshSizeX;
  unsigned int  meshSizeY;
  unsigned int  meshSizeZ;
  double  stiffness;
  unsigned int  numberOfIterations;
  double  edgeCollapseEncouragmentFactor;
  std::string  OutDir; // outdir
  std::string  initMeshFile; // outdir
  int nthreads = 1;

  ArgumentParser parser;
  parser.addHelp(samseg_atlas_help_xml, samseg_atlas_help_xml_len);
  parser.addArgument("--nup", 1, Int);
  parser.addArgument("--iters", 1, Int);
  parser.addArgument("--threads", 1, Int);
  parser.addArgument("--meshsize", 3, Int);
  parser.addArgument("--edgecollapse", 1, Float);
  parser.addArgument("--stiffness", 1, Float);
  parser.addArgument("-o","--outputdir", 1, String, true);
  parser.addArgument("--init", 1, String, false);
  parser.addArgument("-i","--inputs", '+', String, true);
  parser.parse(argc, argv);

  OutDir = parser.retrieve<std::string>("outputdir");
  initMeshFile = parser.retrieve<std::string>("init");
  numberOfUpsamplingSteps = parser.retrieve<int>("nup");
  numberOfIterations = parser.retrieve<int>("iters");
  edgeCollapseEncouragmentFactor = parser.retrieve<float>("edgecollapse");
  stiffness = parser.retrieve<float>("stiffness");
  std::vector<int> meshSize = parser.retrieve<std::vector<int>>("meshsize");
  meshSizeX = meshSize[0];
  meshSizeY = meshSize[1];
  meshSizeZ = meshSize[2];
  std::vector<std::string> inputList = parser.retrieve<std::vector<std::string>>("inputs");

  // nthreads starts as 1. If ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS exists, then it takes that. If 
  // the user explicitly specifies --threads, then it takes that
  if(getenv("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS") != NULL) sscanf(getenv("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"),"%d",&nthreads);
  if(parser.exists("threads")) nthreads = parser.retrieve<int>("threads");
  char tmpstr[2000];
  sprintf(tmpstr,"%d",nthreads);
  setenv("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS",tmpstr,1);

  std::cout << "samseg-atlas Command line params:" << std::endl;
  std::cout << "  numberOfUpsamplingSteps:        " << numberOfUpsamplingSteps << std::endl;
  std::cout << "  meshSizeX:                      " << meshSizeX << std::endl;
  std::cout << "  meshSizeY:                      " << meshSizeY << std::endl;
  std::cout << "  meshSizeZ:                      " << meshSizeZ << std::endl;
  std::cout << "  stiffness:                      " << stiffness << std::endl;
  std::cout << "  numberOfIterations:             " << numberOfIterations << std::endl;
  std::cout << "  edgeCollapseEncouragmentFactor: " << edgeCollapseEncouragmentFactor << std::endl;
  std::cout << "  OutDir:                   " << OutDir << std::endl;
  printf("threads %s\n",getenv("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"));

  PrintMemUsage(stdout);

  // Check whether the output dir exists
  if(fio_FileExistsReadable(OutDir.c_str())){
    printf("ERROR: output folder %s exists, delete to rerun\n",OutDir.c_str());
    exit(1);
  }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  // Read the input images
  typedef kvl::CompressionLookupTable::ImageType  LabelImageType;
  std::vector< LabelImageType::ConstPointer >  labelImages;
  for(int nthinput; nthinput < inputList.size(); nthinput++ ){
    std::cout << "Reading input image: " << inputList[nthinput] << std::endl;
    // Read the input image
    typedef itk::ImageFileReader< LabelImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName(inputList[nthinput]);
    reader->Update();
    LabelImageType::ConstPointer  labelImage = reader->GetOutput();
    // Override the spacing and origin since at this point we can't deal with that
    // Or throw an error?
    const double spacing[] = { 1, 1, 1 };
    const double origin[] = { 0, 0, 0 };
    const_cast< LabelImageType* >( labelImage.GetPointer() )->SetSpacing( spacing );
    const_cast< LabelImageType* >( labelImage.GetPointer() )->SetOrigin( origin );
    // Remember this image
    labelImages.push_back( labelImage );
  }

  // Create outdir
  if(!(itksys::SystemTools::MakeDirectory( OutDir.c_str()))){
    std::cerr << "Couldn't create log directory" << std::endl;
    exit(1);
  }
  // Start the logging file
  std::string  logFileName = OutDir + "/logFile.txt";
  std::ofstream  out( logFileName.c_str(), std::ios::out );
  if (out.bad()){
    std::cerr << "Couldn't write to log file" << std::endl;
    exit(1);
  }

  // Build a lookup table that maps the original intensities onto class numbers starting
  // at 0 and densely packed
  kvl::CompressionLookupTable::Pointer  lookupTable = kvl::CompressionLookupTable::New();
  lookupTable->Construct( labelImages );

  // Write the label images
  typedef itk::ImageFileWriter< LabelImageType >  WriterType;
  for ( int labelImageNumber = 0; labelImageNumber < labelImages.size();labelImageNumber++ ) {
    std::ostringstream  fileNameStream;
    fileNameStream << OutDir << "/labelImage";
    if(labelImageNumber < 10) fileNameStream << "0";
    fileNameStream << labelImageNumber << ".mhd";
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( labelImages[ labelImageNumber ] );
    writer->SetFileName( fileNameStream.str().c_str() ); 
    writer->Write();
  }
  // Write the lookup table
  std::ostringstream  fileNameStream;
  fileNameStream << OutDir << "/compressionLookupTable.txt";
  lookupTable->Write( fileNameStream.str().c_str() );

  out << "iterationNumber     dataCost     alphasCost     positionCost    totalCost"
      << std::endl;
  out << "------------------------------------------------------------------------------" << std::endl;

  // Check for explicit start mesh
  kvl::AtlasMeshCollection::Pointer explicitStartCollection = 0;
  if(!initMeshFile.empty()){
    explicitStartCollection = kvl::AtlasMeshCollection::New();
    if(!explicitStartCollection->Read( initMeshFile.c_str())){
      std::cerr << "Couldn't read mesh from file " << initMeshFile << std::endl;
      exit(1);
    } 
    else{
      std::cout << "explicitStartCollection found; reading from: " << initMeshFile << std::endl;
    }
  }

  PrintMemUsage(stdout);
    
  // Set up the builder
  kvl::AtlasMeshBuilder::Pointer  builder = kvl::AtlasMeshBuilder::New();
  const itk::Size< 3 >  initialSize = { meshSizeX,  meshSizeY, meshSizeZ };
  std::vector< double >  initialStiffnesses( numberOfUpsamplingSteps+1, stiffness );
  builder->SetUp( labelImages, lookupTable, initialSize, initialStiffnesses, numberOfIterations);
  builder->SetVerbose( false );

  // Add some observers/callbacks
  kvl::BuilderCommand::Pointer  builderCommand = kvl::BuilderCommand::New();
  builderCommand->SetLogDirectory( OutDir );
  builder->AddObserver( itk::StartEvent(), builderCommand );
  builder->AddObserver( itk::IterationEvent(), builderCommand );
  builder->AddObserver( itk::EndEvent(), builderCommand );
  builder->AddObserver( kvl::EdgeAnalysisProgressEvent(), builderCommand );
    
  // Let the beast go
  builder->Build(explicitStartCollection, edgeCollapseEncouragmentFactor);

  PrintMemUsage(stdout);
  
  return 0;
};

