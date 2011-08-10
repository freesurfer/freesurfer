#ifndef _SBL_FILE_SYSTEM_H_
#define _SBL_FILE_SYSTEM_H_
#include <sbl/core/Array.h>
#include <sbl/core/String.h>
namespace sbl {


/*! \file FileSystem.h
    \brief The FileSystem module is used for manipulating files and directories.
    These functions represent paths with forward-slash (/) delimiters, 
    not back-slash (\) delimiters, even in Windows.
*/


//-------------------------------------------
// FILE SYSTEM UTILS
//-------------------------------------------


/// obtain a list of files in the given path that match the given prefix and suffix;
/// (both prefix and suffix are allowed to be blank (matches any name))
Array<String> dirFileList( const String &path, const String &prefix, const String &suffix );


/// obtain the file extension from a file name
String fileExtension( const String &fileName );


/// obtain the part before the file extension
String removeFileExtension( const String &fileName );


/// returns the path portion of a file name
String pathPart( const String &fileName );


/// true if file can be open for reading
bool fileExists( const String &fileName );


/// delete the file
void deleteFile( const String &fileName );


/// copy a file to a new name or location
void copyFile( const String &srcFileName, const String &destFileName );


/// move a file to a new name or location
void moveFile( const String &srcFileName, const String &destFileName );


/// append a line to a text file (e.g. for logging)
void appendText( const String &fileName, const String &str );


/// create a directory (assumes all directories exist up to the one being created)
void createDir( const String &path );


/// maintain a fixed number of files with the given extension and path by removing the oldest matching files 
void removeOldestFiles( const String &path, const String &extension, int keepCount );


/// copy the file to a backup file 
void backUpFile( const String &fileName );


/// returns true if the filename has an absolute (non-relative) path
bool isAbsoluteFileName( const String &fileName );


/// the timestamp of the last modification to the given file
int fileModificationTimestamp( const String &fileName );


//-------------------------------------------
// OPERATING SYSTEM UTILITIES
//-------------------------------------------


/// run a shell command
void runSystemCommand( const String &command );


} // end namespace sbl
#endif // _SBL_FILE_SYSTEM_H_

