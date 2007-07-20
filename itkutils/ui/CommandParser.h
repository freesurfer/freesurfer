#ifndef COMMAND_PARSER_H
#define COMMAND_PARSER_H

#include <vector>

/**
 * 
 */
class CommandParser {

private:

  int m_nArgs;
  char **m_aArgs;
  
public:

  static const int FLAG_NOT_FOUND = -1;

  CommandParser( int inArgs, char ** iaArgs);
  ~CommandParser();
  
  int GetNumberOfArguments();
  char **GetArguments();

  int GetPosition( const char *iFlag );  
  char *GetArgument( const char *iFlag );
  int GetArgumentInt( const char *iFlag );
  std::vector< int > GetArgumentIntVector( const char *iFlag );
  double GetArgumentDouble( const char *iFlag );
  bool GetArgumentBoolean( const char *iFlag );

};

#endif
