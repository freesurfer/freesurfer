#ifndef FREE_SURFER_EXECUTABLE_H
#define FREE_SURFER_EXECUTABLE_H

#include <string>
#include <vector>

 
#include "mri.h"


#include "CommandParser.h"

class FreeSurferExecutable {

private:
  struct Argument {
    std::string flag; 
    std::string shortDescription;
    std::string longDescription;
    std::string example;
    std::string errorMessage;
  };

protected:
  std::string m_Name;
  std::string m_Description;
  std::string m_Output;
  std::string m_Version;
  std::string m_BugEmail;
  
  std::vector< Argument* > m_RequiredArguments;
  std::vector< Argument* > m_OptionalArguments;
  
  CommandParser *m_Parser;
  
  std::string m_CommandLineArguments;
  
  void SetArgument( std::vector< Argument* > *iArguments,
    std::string iFlag,
    std::string iShortDescription,
    std::string iLongDescription,
    std::string iExample,
    std::string iErrorMessage );
    
  void PrintArgumentsHelp( std::vector< Argument* > iArguments );

  std::string GetFieldAndParameter( std::string field, std::string parameter );
  std::string GetCurrentDirectory();
  std::string GetOperatingSystem();
  std::string GetHostName();
  std::string GetMachine();
  std::string GetEnvironmentVariable( std::string variable );
  std::string GetFreeSurferHome();
  std::string GetSubjectsDirectory();
  std::string GetFileExtension( std::string fileName );

  void PrintUsageError( std::string error );
  
  static bool IsWhitespace( const char x );

  static std::string FormatString( std::string toBeFormated, 
    const std::string spacing );
  
public:

  FreeSurferExecutable( int inArgs, char ** iaArgs );  
  virtual ~FreeSurferExecutable();
  
  void PrintHelp();
  
  void SetNextRequiredArgument( std::string iFlag, std::string iShortDescription, 
    std::string iLongDescription, std::string iExample, 
    std::string iErrorMessage );

  void SetNextOptionalArgument( std::string iFlag, std::string iShortDescription,
    std::string iLongDescription );

  void SetName( std::string iName, std::string iDescription );
  std::string GetName();

  std::string GetSynposis();

  void SetOutput( std::string iOutput );
    
  std::string GetExampleUsage();

  void SetVersion( std::string iVersion );
  std::string GetVersion();
  
  void SetBugEmail( std::string iEmail );
  
  static void PrintError( std::string iError );

  static void PrintErrorArgument( std::string iArgument, 
    std::string iDetailedMessage );
    
  static void WriteData( const std::string fileName, const double *data, 
    const int nData);

  static void WriteData( const std::string fileName, 
    double **dataArray, const int nRows, const int cCols );

  static void WriteData( const std::string fileName, 
    const float *data, const int nRows, const int nCols );

  static void WriteDataAppend( const std::string fileName, 
    double **dataArray, const int nRows, const int cCols );

  static void WriteDataAppend( const std::string fileName, 
    const float *data, const int nRows, const int nCols );
    
  static void ReadData( const char* fileName, const int cCols, 
    std::vector< double* > *oData );
    
  static MRI* ReadMRI( const char* fileName );  

  std::string GetCommandArguments( int argc, char * argv[] );

  /**
   * Returns the input arguments for the flags in the same order that the
   * flags were originally input input.
   */
  std::string* GetRequiredArguments();
  
  virtual bool FillArguments() = 0;
  virtual void Run() = 0;
    
};

#endif
