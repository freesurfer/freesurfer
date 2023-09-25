#include "kvlAtlasMeshBuilder.h"
#include "itkCommand.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "kvlCompressionLookupTable.h"
#include "itkMGHImageIOFactory.h"



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

  // Sanity check on input
  if ( argc < 8 )
    {
    std::cerr << "Usage: " << argv[ 0 ] << " numberOfUpsamplingSteps meshSizeX meshSizeY meshSizeZ stiffness numberOfIterations edgeCollapseFactor logDirectory fileName1 [ fileName2 ... ]" << std::endl;

    return -1;
    }
  printf("\n\n\n\n\n\n\n\nConsider using samseg-atlas instead ... continuing after 2 sec\n\n\n\n\n\n\n\n\n");
  sleep(2);
    
  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  // Retrieve the input parameters
  std::ostringstream  inputParserStream;
  for ( int argumentNumber = 1; argumentNumber < 9; argumentNumber++ ) 
    {
    inputParserStream << argv[ argumentNumber ] << " ";
    }
  std::istringstream  inputStream( inputParserStream.str().c_str() );
  int  numberOfUpsamplingSteps;
  unsigned int  meshSizeX;
  unsigned int  meshSizeY;
  unsigned int  meshSizeZ;
  double  stiffness;
  unsigned int  numberOfIterations;
  double  edgeCollapseEncouragmentFactor;
  std::string  logDirectory;
  inputStream >> \
    numberOfUpsamplingSteps >> \
    meshSizeX >> meshSizeY >> meshSizeZ >> \
    stiffness >> \
    numberOfIterations >> \
    edgeCollapseEncouragmentFactor >> \
    logDirectory;

  std::cout << "kvlBuildAtlasMesh Command line params:" << std::endl;
  std::cout << "  numberOfUpsamplingSteps:        " << numberOfUpsamplingSteps << std::endl;
  std::cout << "  meshSizeX:                      " << meshSizeX << std::endl;
  std::cout << "  meshSizeY:                      " << meshSizeY << std::endl;
  std::cout << "  meshSizeZ:                      " << meshSizeZ << std::endl;
  std::cout << "  stiffness:                      " << stiffness << std::endl;
  std::cout << "  numberOfIterations:             " << numberOfIterations << std::endl;
  std::cout << "  edgeCollapseEncouragmentFactor: " << edgeCollapseEncouragmentFactor << std::endl;
  std::cout << "  logDirectory:                   " << logDirectory << std::endl;
  
  // Read the input images
  typedef kvl::CompressionLookupTable::ImageType  LabelImageType;
  std::vector< LabelImageType::ConstPointer >  labelImages;
  for ( int argumentNumber = 9; argumentNumber < argc; argumentNumber++ )
    {
    std::cout << "Reading input image: " << argv[ argumentNumber ] << std::endl;
    // Read the input image
    typedef itk::ImageFileReader< LabelImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( argv[ argumentNumber ] );
    reader->Update();
    LabelImageType::ConstPointer  labelImage = reader->GetOutput();

    // Over-ride the spacing and origin since at this point we can't deal with that
    const double spacing[] = { 1, 1, 1 };
    const double origin[] = { 0, 0, 0 };
    const_cast< LabelImageType* >( labelImage.GetPointer() )->SetSpacing( spacing );
    const_cast< LabelImageType* >( labelImage.GetPointer() )->SetOrigin( origin );

    // Remember this image
    labelImages.push_back( labelImage );
    }

  // Build a lookup table that maps the original intensities onto class numbers starting
  // at 0 and densely packed
  kvl::CompressionLookupTable::Pointer  lookupTable = kvl::CompressionLookupTable::New();
  lookupTable->Construct( labelImages );
  //lookupTable->Write( "compressionLookupTable.txt" );


  // Set up the builder
  kvl::AtlasMeshBuilder::Pointer  builder = kvl::AtlasMeshBuilder::New();
  const itk::Size< 3 >  initialSize = { meshSizeX,  meshSizeY, meshSizeZ };
  std::vector< double >  initialStiffnesses( numberOfUpsamplingSteps+1, stiffness );
  builder->SetUp( labelImages, lookupTable, initialSize, initialStiffnesses, numberOfIterations);
  builder->SetVerbose( false );

  // Add some observers/callbacks
  kvl::BuilderCommand::Pointer  builderCommand = kvl::BuilderCommand::New();
  builderCommand->SetLogDirectory( logDirectory );
  builder->AddObserver( itk::StartEvent(), builderCommand );
  builder->AddObserver( itk::IterationEvent(), builderCommand );
  builder->AddObserver( itk::EndEvent(), builderCommand );
  builder->AddObserver( kvl::EdgeAnalysisProgressEvent(), builderCommand );

  // Make sure the log directory exists and is empty
  if ( itksys::SystemTools::FileIsDirectory( logDirectory.c_str() ) )
    {
    // Logging directory exists already. Remove
    if ( !( itksys::SystemTools::RemoveADirectory( logDirectory.c_str() ) ) )
      {
      std::cerr << "Couldn't remove existing log directory" << std::endl;
      return -1;
      }
    }

  if ( !( itksys::SystemTools::MakeDirectory( logDirectory.c_str() ) ) )
    {
    std::cerr << "Couldn't create log directory" << std::endl;
    return -1;
    }

  // Write the lookup table
  std::ostringstream  fileNameStream;
  fileNameStream << logDirectory << "/compressionLookupTable.txt";
  lookupTable->Write( fileNameStream.str().c_str() );

    
  // Write the label images
  typedef itk::ImageFileWriter< LabelImageType >  WriterType;
  for ( int labelImageNumber = 0; 
        labelImageNumber < labelImages.size();
        labelImageNumber++ )
    {
    std::ostringstream  fileNameStream;
    fileNameStream << logDirectory << "/labelImage";
    if ( labelImageNumber < 10 )
      {
      fileNameStream << "0";
      }
    fileNameStream << labelImageNumber << ".mhd";
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( labelImages[ labelImageNumber ] );
    writer->SetFileName( fileNameStream.str().c_str() ); 
    writer->Write();
    }


  // Also start the loggin file
  std::string  logFileName = logDirectory;
  logFileName += "/logFile.txt";
  std::ofstream  out( logFileName.c_str(), std::ios::out );
  if ( out.bad() )
    {
    std::cerr << "Couldn't write to log file" << std::endl;
    return -1;
    }
  out << "iterationNumber     dataCost     alphasCost     positionCost    totalCost"
      << std::endl;
  out << "------------------------------------------------------------------------------" << std::endl;


  // If explicitStartCollection exists in the current directory, use it
  kvl::AtlasMeshCollection::Pointer  explicitStartCollection = nullptr;
  const std::string  explicitStartCollectionFileName = "explicitStartCollection.gz";
  if ( itksys::SystemTools::FileExists( explicitStartCollectionFileName.c_str(), true ) )
    {
    explicitStartCollection = kvl::AtlasMeshCollection::New();
    if ( !explicitStartCollection->Read( explicitStartCollectionFileName.c_str() ) )
      {
      std::cerr << "Couldn't read mesh from file " << explicitStartCollectionFileName << std::endl;
      return -1;
      } 
    else
      {
      std::cout << "explicitStartCollection found; reading from: " << explicitStartCollectionFileName << std::endl;
      }
    }

  // If edgeCollapseEncouragmentFactor.txt exists in the current directory, read it's content
  //double  edgeCollapseEncouragmentFactor = 1.0;
  const std::string  edgeCollapseEncouragmentFactorFileName = "edgeCollapseEncouragmentFactor.txt";
  //if ( itksys::SystemTools::FileExists( edgeCollapseEncouragmentFactorFileName.c_str(), true ) )
  //  {
    std::ifstream  fs( edgeCollapseEncouragmentFactorFileName.c_str() );
    if ( !( fs.fail() ) )
      {
      std::cout << "Reading " << edgeCollapseEncouragmentFactorFileName << std::endl;

      std::string line;
      if ( std::getline( fs, line ) )
        {
        //std::ostringstream  inputParserStream;
        //inputParserStream << line;
        //std::istringstream  inputStream( inputParserStream.str().c_str() );
        std::istringstream  inputStream( line.c_str() );
        inputStream >> edgeCollapseEncouragmentFactor;
        std::cout << "Using edgeCollapseEncouragmentFactor: " << edgeCollapseEncouragmentFactor << std::endl;
        }
      else
        {
        std::cerr << "Couldn't read edgeCollapseEncouragmentFactor from file: " 
                  << edgeCollapseEncouragmentFactorFileName << std::endl;
        return -1;
        }
        
      }  
   
  //  }
    
    
  // Let the beast go
  builder->Build( explicitStartCollection, edgeCollapseEncouragmentFactor );

  
  return 0;
};

