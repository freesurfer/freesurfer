#include "kvlAtlasMeshBuilder.h"
#include "itkCommand.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "kvlCompressionLookupTable.h"
#include "kvlAtlasMeshCollectionPositionCostCalculator2.h"
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

    if ( ( typeid( event ) == typeid( itk::StartEvent ) ) ||
          ( typeid( event ) == typeid( itk::IterationEvent ) ) ||
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
#if 0        
      const float  currentDataCost = builder->GetCurrentDataCost();
      const float  currentAlphasCost = builder->GetCurrentAlphasCost();
#else
      float  currentDataCost = 0;
      float  currentAlphasCost = 0;
      std::cout << "Computing current data and alphas cost..." << std::endl;
      itk::TimeProbe  probe;
      probe.Start();
      builder->GetCurrentDataAndAlphasCost( currentDataCost, currentAlphasCost );
      probe.Stop();
      std::cout << "...done" << std::endl;
      std::cout << "Took " << probe.GetMean() << " seconds to compute current data and alphas cost" << std::endl;
#endif
      
#if 0
      const float  currentPositionCost = builder->GetCurrentPositionCost();
      const float  currentCost = builder->GetCurrentCost();
#else
      const float  currentPositionCost = 0;
      const float  currentCost = currentDataCost + currentAlphasCost;
#endif
      out << builder->GetIterationNumber() << "   "
           << currentDataCost << "   "
           << currentAlphasCost << "   "
           << currentPositionCost << "   "
           << currentCost << std::endl;

      }
    else if ( ( typeid( event ) == typeid( StartEdgeAnalysisEvent ) ) || 
                  ( typeid( event ) == typeid( EndEdgeAnalysisEvent ) ) )
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
  if ( argc < 9 )
    {
    std::cerr << "Usage: " << argv[ 0 ] << " collapsedLabelFile('null' if no compressed labels) numberOfUpsamplingSteps meshSizeX meshSizeY meshSizeZ stiffness logDirectory fileName1 [ fileName2 ... ]" << std::endl;

    return -1;
    }
    
  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  // Retrieve the input parameters
  std::ostringstream  inputParserStream;
  for ( int argumentNumber = 1; argumentNumber < 8; argumentNumber++ ) 
    {
    inputParserStream << argv[ argumentNumber ] << " ";
    }
  std::istringstream  inputStream( inputParserStream.str().c_str() );
  std::string  collapsedLabelFile;
  unsigned int  numberOfUpsamplingSteps;
  unsigned int  meshSizeX;
  unsigned int  meshSizeY;
  unsigned int  meshSizeZ;
  float  stiffness;
  std::string  logDirectory;
  inputStream >>  collapsedLabelFile  >> numberOfUpsamplingSteps >>  meshSizeX >> meshSizeY >> meshSizeZ >> stiffness >> logDirectory;

  // read in list of collapsed labels
  std::ifstream clfs(collapsedLabelFile.c_str());
  std::vector<unsigned short> collapsedLabs;
  std::vector<unsigned short> mappedCollapsedLabs;
  std::vector<std::vector<unsigned short> > collapsedLabList;

  if(!(clfs.fail()))
  {
    std::string line;
    unsigned short comp, lab;
    int count=256; // collapsed labels go from 255 down
    while (std::getline(clfs, line))
    {
      count--;
      unsigned char mcl = (unsigned char) count;
      mappedCollapsedLabs.push_back(mcl);
      std::vector<unsigned short> ll;
      std::ostringstream  inputParserStream;
      inputParserStream << line;
      std::istringstream  inputStream( inputParserStream.str().c_str() );
      if (!(inputStream >>  comp  >> lab))
      {
        std::cerr << "Error reading collapsed label file" << std::endl;
        return -1;
      }
      collapsedLabs.push_back(comp);
      ll.push_back(lab);
      while(inputStream >>  lab) { ll.push_back(lab); }
      collapsedLabList.push_back(ll);
    }
  }

 
  // Read the images (in original ushort format)
  try
    {
    typedef kvl::CompressionLookupTable::ImageType  InputImageType;
    std::vector< InputImageType::ConstPointer >  originalImages;
    for ( int argumentNumber = 8; argumentNumber < argc; argumentNumber++ )
      {
      // Read the input image
      typedef itk::ImageFileReader< InputImageType >  ReaderType;
      ReaderType::Pointer  reader = ReaderType::New();
      reader->SetFileName( argv[ argumentNumber ] );
      reader->Update();
      InputImageType::ConstPointer  originalImage = reader->GetOutput();

      // Over-ride the spacing and origin since at this point we can't deal with that
      const double spacing[] = { 1, 1, 1 };
      const double origin[] = { 0, 0, 0 };
      const_cast< InputImageType* >( originalImage.GetPointer() )->SetSpacing( spacing );
      const_cast< InputImageType* >( originalImage.GetPointer() )->SetOrigin( origin );

      // Remember this image
      originalImages.push_back( originalImage );
      }

    // Build a lookup table that maps the original intensities onto uchar starting
    // at 0 and densely packed
    kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
    std::vector<std::vector<unsigned short> > modifiedCollapsedLabList = compressor->Construct( originalImages,collapsedLabs,collapsedLabList ); // 

    compressor->Write( "compressionLookupTable.txt" );

    // Collect the label images resulting from pushing the original images through the
    // lookup table
    typedef kvl::CompressionLookupTable::CompressedImageType  OutputImageType;
    std::vector< OutputImageType::ConstPointer >  labelImages;
    for ( std::vector< InputImageType::ConstPointer >::const_iterator it = originalImages.begin();
          it != originalImages.end(); ++it )
      {
        labelImages.push_back( ( compressor->CompressImage( ( *it ).GetPointer() , collapsedLabs, mappedCollapsedLabs) ).GetPointer() );
    }


    // Now build the lookup table from compressed / collapsed to compressed  + compressed to compress (this is just the identity)
    std::vector<unsigned char > mapCompToComp[256];
    for(int i=0; i<256; i++)
    {
      std::vector<unsigned char > v;
      if ( i<compressor->GetCompressionLookupTable().size()){
        v.push_back(i);
      }
      else if( i > 255-modifiedCollapsedLabList.size())
      {
        int ind=255-i;
        for (int j=0; j<modifiedCollapsedLabList[ind].size(); j++)
          v.push_back(modifiedCollapsedLabList[ind][j]);
      }
      mapCompToComp[i]=v; 
//      if(mapCompToComp[i].size()>0) { std::cout << i << ": "; for(int k=0; k<mapCompToComp[i].size(); k++) std::cout << ((int)mapCompToComp[i][k]) << " "; std::cout << std::endl; }
    }

    // Set up the builder
    kvl::AtlasMeshBuilder::Pointer  builder = kvl::AtlasMeshBuilder::New();
    builder->SetMapCompToComp(mapCompToComp);
    builder->SetLabelImages( labelImages );
    builder->SetNumberOfUpsamplingSteps( numberOfUpsamplingSteps );
    unsigned int  size[ 3 ] = { meshSizeX,  meshSizeY, meshSizeZ };
    float  stiffnesses[ 5 ] = { stiffness, stiffness, stiffness, stiffness, stiffness };
    builder->SetInitialMesh( size, stiffnesses );
    kvl::BuilderCommand::Pointer  builderCommand = kvl::BuilderCommand::New();
    builderCommand->SetLogDirectory( logDirectory );
    builder->AddObserver( itk::StartEvent(), builderCommand );
    builder->AddObserver( itk::IterationEvent(), builderCommand );
    builder->AddObserver( itk::EndEvent(), builderCommand );
    builder->AddObserver( kvl::StartEdgeAnalysisEvent(), builderCommand );
    builder->AddObserver( kvl::EndEdgeAnalysisEvent(), builderCommand );

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

    // Write the label images
    typedef itk::ImageFileWriter< OutputImageType >  WriterType;
    for ( unsigned int labelImageNumber = 0; 
          labelImageNumber < builder->GetLabelImages().size();
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
      writer->SetInput( builder->GetLabelImages()[ labelImageNumber ] );
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

    // Don't use position cost
    kvl::AtlasMeshCollectionPositionCostCalculator::SetReturnZero( true );


#if 1
    // If explicitStartCollection exists in the current directory, use it
    kvl::AtlasMeshCollection::Pointer  explicitStartCollection = 0;
    //const std::string  explicitStartCollectionFileName = "explicitStartCollection.txt.gz";
    const std::string  explicitStartCollectionFileName = "explicitStartCollection.gz";
    if ( itksys::SystemTools::FileExists( explicitStartCollectionFileName.c_str(), true ) )
      {
      explicitStartCollection = kvl::AtlasMeshCollection::New();
      if ( !explicitStartCollection->Read( explicitStartCollectionFileName.c_str() ) )
        {
        std::cerr << "Couldn't read mesh from file " << explicitStartCollectionFileName << std::endl;
        return -1;
        }
      }

    // Let the beast go
    builder->Build( explicitStartCollection );

#else
    // Let the beast go
    builder->Build();
#endif

    //builder->BuildWithPriorityQueue();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << e << std::endl;
    }

  return 0;
};

