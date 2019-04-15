#include <sstream>
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdlib.h> // getenv
#include <unistd.h>
#include <sys/utsname.h>

#include "FreeSurferExecutable.h"


FreeSurferExecutable::FreeSurferExecutable( int inArgs, char ** iaArgs ) {
  m_Parser = new CommandParser( inArgs, iaArgs );
  
  m_CommandLineArguments = GetCommandArguments( inArgs, iaArgs );
}

FreeSurferExecutable::~FreeSurferExecutable() {
  delete m_Parser;
  m_Parser = NULL;
  
// TODO: delete arguments  
}

void
FreeSurferExecutable::SetArgument( std::vector< Argument* > *iArguments, 
  std::string iFlag,
  std::string iShortDescription,
  std::string iLongDescription,
  std::string iExample,
  std::string iErrorMessage ) {

  Argument *argument = new Argument();
  argument->flag = iFlag;
  argument->shortDescription = iShortDescription;
  argument->longDescription = iLongDescription;
  argument->example = iExample;
  argument->errorMessage = iErrorMessage;

  iArguments->insert( iArguments->end(), argument );
} 

void 
FreeSurferExecutable::SetNextRequiredArgument( std::string iFlag, 
  std::string iShortDescription, 
  std::string iLongDescription,
  std::string iExample,
  std::string iErrorMessage ) {

  SetArgument( &m_RequiredArguments, iFlag, iShortDescription, iLongDescription,
    iExample, iErrorMessage );

}

void
FreeSurferExecutable::SetNextOptionalArgument( std::string iFlag,
  std::string iShortDescription,
  std::string iLongDescription ) {

  SetArgument( &m_OptionalArguments, iFlag, iShortDescription, iLongDescription,
    "", "" );

}

void
FreeSurferExecutable::SetName( std::string iName, 
                                  std::string iDescription ) {
  m_Name = iName;
  m_Description = iDescription;
}

std::string
FreeSurferExecutable::GetName() {
  return m_Name;
}

std::string
FreeSurferExecutable::GetSynposis() {
  std::string synopsis = m_Name;
  
  for( std::vector< Argument* >::const_iterator arguments = m_RequiredArguments.begin(); 
    arguments != m_RequiredArguments.end(); arguments++ ) {

    Argument *argument = *arguments;
    synopsis = synopsis + " " + argument->flag + 
               " <" + argument->shortDescription + ">";
  }  
  
  synopsis = synopsis + " [<options>]";
  
  return synopsis;
}

void
FreeSurferExecutable::SetOutput( std::string iOutput ) {
  m_Output = iOutput;
}

std::string
FreeSurferExecutable::GetExampleUsage() {
  std::string example = m_Name;
  
  for( std::vector< Argument* >::const_iterator arguments = m_RequiredArguments.begin(); 
    arguments != m_RequiredArguments.end(); arguments++ ) {

    Argument *argument = *arguments;
    example = example + " " + argument->flag + " " + argument->example;
  }  
  
  return example;
}

void
FreeSurferExecutable::SetVersion( std::string iVersion ) {
  m_Version = iVersion;
}

std::string
FreeSurferExecutable::GetVersion() {
  return m_Version;
}

void
FreeSurferExecutable::SetBugEmail( std::string iEmail ) {
  m_BugEmail = iEmail;
}

void
FreeSurferExecutable::PrintArgumentsHelp( 
  std::vector< Argument* > iArguments ) {
  
  for( std::vector< Argument* >::const_iterator arguments = iArguments.begin(); 
    arguments != iArguments.end(); arguments++ ) {

    Argument *argument = *arguments;
    std::cout << "   " << argument->flag << " <" << argument->shortDescription << ">" << std::endl;

    if( argument->longDescription != "" ) {
      std::cout << FormatString( argument->longDescription, "     " ) << std::endl;
    }

    std::cout << std::endl;
  }  
}

void
FreeSurferExecutable::PrintHelp() {
  std::string spaces = "   ";
  
  std::cout << "Name\n" << std::endl;
  std::cout << spaces << m_Name << " - " << m_Description << std::endl;

  std::cout << std::endl;

  std::cout << "Synopsis" << std::endl;
  std::cout << FormatString( GetSynposis(), spaces ) << std::endl;
 
  std::cout << std::endl;
 
  std::cout << "Required Arguments" << std::endl;
  PrintArgumentsHelp( m_RequiredArguments );

  std::cout << "Optional Arguments" << std::endl;
  PrintArgumentsHelp( m_OptionalArguments );
  
  std::cout << "Output" << std::endl;
  std::cout << spaces << m_Output << std::endl;  
  
  std::cout << std::endl;
  
  std::cout << "Example" << std::endl;
  std::cout << FormatString( GetExampleUsage(), spaces ) << std::endl;

  std::cout << std::endl;
  
  std::cout << "Version" << std::endl;
  std::cout << spaces << m_Version << std::endl;

  std::cout << std::endl;

  std::cout << "Reporting Bugs" << std::endl;
  std::cout << spaces << "Report bugs to <" << m_BugEmail << ">" << std::endl;

  std::cout << std::endl;
}

void 
FreeSurferExecutable::PrintError( std::string iError ) {
  std::cout << "ERROR: " << iError << std::endl;
}

void 
FreeSurferExecutable::PrintErrorArgument( std::string iArgument, 
    std::string iDetailedMessage ) {
    
  std::string errorMessage = iArgument + " argument missing";

  PrintError( errorMessage );
  PrintError( iDetailedMessage );
}

std::string*
FreeSurferExecutable::GetRequiredArguments() {
  std::string *requiredArgumentValues = 
    new std::string[ m_RequiredArguments.size() ];
  
  int i=0;
  for( std::vector< Argument* >::const_iterator arguments = m_RequiredArguments.begin(); 
       arguments != m_RequiredArguments.end(); i++, arguments++ ) {

    Argument *argument = *arguments;
    char *argumentStr = m_Parser->GetArgument( argument->flag.c_str() );

    if( argumentStr == NULL ) {
      PrintErrorArgument( argument->flag, argument->errorMessage );
      // TODO: throw an error
      throw "invalid";
    }
    
    requiredArgumentValues[i] = argumentStr;
  }
  
  return requiredArgumentValues;
}

void 
FreeSurferExecutable::WriteData( const std::string fileName, 
  const double *data, const int nData) {

  std::ofstream output( fileName.c_str() );
  
  for( int cData=0; cData<nData; cData++ ) {
    output << data[ cData ] << std::endl;
  }
  
  output.close();
    
}

void 
FreeSurferExecutable::WriteData( const std::string fileName, 
  double **dataArray, const int nRows, const int nCols ) {
    
  std::ofstream output( fileName.c_str() );
  
  for( int cRow=0; cRow<nRows; cRow++ ) {
    for( int cCol=0; cCol<nCols; cCol++ ) {
      output << dataArray[ cRow ][ cCol ] << "   ";
    }    
    output << std::endl;
  }
  
  output.close();
    
}

void 
FreeSurferExecutable::WriteData( const std::string fileName, 
  const float *data, const int nRows, const int nCols) {

  std::ofstream output( fileName.c_str() );
  
  for( int cRow=0; cRow<nRows; cRow++ ) {
    for( int cCol=0; cCol<nCols; cCol++ ) {
      const int index = cCol + cRow * nCols;      
      output << data[ index ] << "   ";
    }
    output << std::endl;
  }
  
  output.close();    
}

void 
FreeSurferExecutable::WriteDataAppend( const std::string fileName, 
  double **dataArray, const int nRows, const int nCols ) {
    
  std::ofstream output( fileName.c_str(), std::ios::out | std::ios::app );
  
  for( int cRow=0; cRow<nRows; cRow++ ) {
    for( int cCol=0; cCol<nCols; cCol++ ) {
      output << dataArray[ cRow ][ cCol ] << "   ";
    }    
    output << std::endl;
  }
  
  output.close();
    
}

void 
FreeSurferExecutable::WriteDataAppend( const std::string fileName, 
  const float *data, const int nRows, const int nCols) {

  std::ofstream output( fileName.c_str(), std::ios::out | std::ios::app );
  
  for( int cRow=0; cRow<nRows; cRow++ ) {
    for( int cCol=0; cCol<nCols; cCol++ ) {
      const int index = cCol + cRow * nCols;      
      output << data[ index ] << "   ";
    }
    output << std::endl;
  }
  
  output.close();    
}

void 
FreeSurferExecutable::ReadData( const char* fileName, 
  const int cCols, std::vector< double* > *oData ) {
  
  std::string s;
  std::ifstream dataFile( fileName );

  double value = 0.0;
  
  // row of data
  double *dataRow = NULL;

  // current index of the data
  int nData = 0;
  
  // read all the samples
  while( dataFile >> s ) {
    
    // position in columns
    const int nCol = nData % cCols;
    
    // new row of data
    if( nCol == 0 ) {
      dataRow = new double[ cCols ];
      oData->push_back( dataRow );
    }
    
    // convert the value to double format
    std::stringstream stream;
    stream << s;
    stream >> value;
    
    dataRow[ nCol ] = value;
    
    nData++;
  }
  
}

MRI* 
FreeSurferExecutable::ReadMRI( const char* fileName ) {

  // read the eigenvectors
  char* nonConstFileName = strdup( fileName );
  MRI* volume = ::MRIread( nonConstFileName );
  free( nonConstFileName );

  if ( volume == NULL ) {
    // TODO: throw an error
    std::cerr << "mri read failed..." << std::endl;
  }
  
  return volume;
}

std::string 
FreeSurferExecutable::GetCommandArguments( int argc, char * argv[] ) {  

  std::ostringstream output;

  for( int cArg=0; cArg<argc; cArg++ ) {
    output << argv[ cArg ] << " ";
  }
  output << std::endl;
  
  return output.str();
}

std::string 
FreeSurferExecutable::GetFieldAndParameter( 
  std::string field, std::string parameter ) {
    
  std::ostringstream output;
  output << field << ": \n" << "  " << parameter << std::endl << std::endl;
  return output.str();
}

std::string 
FreeSurferExecutable::GetCurrentDirectory() {
  
  const int nChars = 2000;
  char cwd[ nChars ];
  
  getcwd( cwd, nChars );
  
  return std::string( cwd );
}

std::string 
FreeSurferExecutable::GetOperatingSystem() {  
  utsname uts;
  uname( &uts );    
  std::string name( uts.sysname );
  return name;
}

std::string
FreeSurferExecutable::GetHostName() {
  utsname uts;
  uname( &uts );
  std::string name( uts.nodename );
  return name;
}

std::string
FreeSurferExecutable::GetMachine() {
  utsname uts;
  uname( &uts );
  std::string machine( uts.machine );
  return machine;
}

std::string
FreeSurferExecutable::GetEnvironmentVariable( std::string variable ) {
  char *result = getenv( variable.c_str() );
  std::string stringResult;
  if( result == NULL ) {
    stringResult = "not set";
  } else {
    stringResult = std::string( result );
  }
  return stringResult;
}

std::string 
FreeSurferExecutable::GetFreeSurferHome() {
  return GetEnvironmentVariable( "FREESURFER_HOME" );
}

std::string 
FreeSurferExecutable::GetSubjectsDirectory() {
  return GetEnvironmentVariable( "SUBJECTS_DIR" );
}

std::string
FreeSurferExecutable::GetFileExtension( std::string fileName ) {

  // take in part from itk  
  std::string::size_type ext = fileName.rfind( '.' );
  std::string exts = fileName.substr( ext );
  if( exts == ".gz" ) {
    std::string::size_type dotpos = fileName.rfind( '.', ext-1 );
    if( dotpos != std::string::npos ) {
      exts = fileName.substr(dotpos);
    }
  }  
  
  return exts;
}

void 
FreeSurferExecutable::PrintUsageError( std::string error ) {
  std::cerr << "\n*** usage error: " << error << std::endl << std::endl;
}

bool 
FreeSurferExecutable::IsWhitespace( const char x ) {
    return x == ' ' || x == '\n' || x == '\r' || x == '\t';
}

std::string 
FreeSurferExecutable::FormatString( 
  std::string toBeFormatted, const std::string spacing ) {

  const int endOfScreen = 80 - spacing.size();
  std::string formatted = spacing;
  
  int countDown = endOfScreen;
  for( std::string::iterator stringIt = toBeFormatted.begin();
    stringIt != toBeFormatted.end(); stringIt++ ) {

    if( countDown == 0 ){
      while( !IsWhitespace( *stringIt ) && stringIt != toBeFormatted.begin() ) {
        stringIt--;
        formatted.resize( formatted.size() - 1 );
      }
      stringIt++;
      formatted += "\n" + spacing;
      countDown = endOfScreen;
    }
        
    formatted += *stringIt;

    countDown--;    
  }
    
  return formatted;
}
