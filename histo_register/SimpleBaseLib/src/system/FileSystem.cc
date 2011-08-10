// Licensed under MIT license; see license.txt.

#include <sbl/system/FileSystem.h>
#include <sbl/core/File.h>
#include <sbl/core/StringUtil.h>
#include <sbl/math/VectorUtil.h> // for RemoveOldestFiles
#include <stdlib.h>
#include <sys/stat.h>
#ifdef WIN32
    #include <external/win_dirent.h>
#else
    #include <dirent.h>
    #include <unistd.h>
#endif
namespace sbl {


//-------------------------------------------
// FILE SYSTEM UTILITIES
//-------------------------------------------


/// obtain a list of files in the given path that match the given prefix and suffix;
/// (both prefix and suffix are allowed to be blank (matches any name))
Array<String> dirFileList( const String &path, const String &prefix, const String &suffix ) {
    Array<String> dirList;
    DIR *dir = opendir( path.c_str() );
    if (dir) {
        dirent *ent = readdir( dir );
        while (ent) {
            String *dirEnt = new String( ent->d_name );
            if (dirEnt->startsWith( prefix ) && dirEnt->endsWith( suffix ))
                dirList.append( dirEnt );
            else
                delete dirEnt;
            ent = readdir( dir );
        }
        closedir( dir );
    }
    return sort( dirList );
}


/// obtain the file extension from a file name
String fileExtension( const String &fileName ) {
    int dotPos = fileName.lastCharPos( '.' );
    if (dotPos > 0)
        return fileName.rightOf( dotPos );
    else
        return String();
}


/// obtain the part before the file extension
String removeFileExtension( const String &fileName ) {
    int dotPos = fileName.lastCharPos( '.' );
    if (dotPos > 0)
        return fileName.leftOf( dotPos );
    else
        return fileName;
}


/// returns the path portion of a file name
String pathPart( const String &fileName ) {
    String path;
    int pos1 = fileName.lastCharPos( '/' );
    int pos2 = fileName.lastCharPos( '\\' );
    if (pos1 >= 0 || pos2 >= 0) {
        if (pos1 > pos2) 
            path = fileName.leftOf( pos1 );
        else
            path = fileName.leftOf( pos2 );
    }
    return path;
}


/// true if file can be open for reading
// fix(faster): don't open file
bool fileExists( const String &fileName ) {
    File file( fileName, FILE_READ, FILE_BINARY );
    return file.openSuccess();
}


/// delete the file
void deleteFile( const String &fileName ) {
#if WIN32
    DeleteFileA( fileName.c_str() );        
#else
    unlink( fileName.c_str() );
#endif
}


/// copy a file to a new name or location
void copyFile( const String &srcFileName, const String &destFileName ) {
#ifdef WIN32
    String srcFileNameClean = srcFileName.replace( '/', '\\' );
    String destFileNameClean = destFileName.replace( '/', '\\' );
    runSystemCommand( sprintF( "copy %s %s", srcFileNameClean.c_str(), destFileNameClean.c_str() ) );
#else
    runSystemCommand( sprintF( "cp -f %s %s", srcFileName.c_str(), destFileName.c_str() ) );
#endif
}


/// move a file to a new name or location
void moveFile( const String &srcFileName, const String &destFileName ) {
#ifdef WIN32
    String srcFileNameClean = srcFileName.replace( '/', '\\' );
    String destFileNameClean = destFileName.replace( '/', '\\' );
    runSystemCommand( sprintF( "move %s %s", srcFileNameClean.c_str(), destFileNameClean.c_str() ) );
#else
    runSystemCommand( sprintF( "mv -f %s %s", srcFileName.c_str(), destFileName.c_str() ) );
#endif
}


/// append a line to a text file (e.g. for logging)
void appendText( const String &fileName, const String &str ) {
    File file( fileName, FILE_APPEND, FILE_TEXT ); 
    if (file.openSuccess())
        file.writeF( str.c_str() );
}


/// create a directory (assumes all directories exist up to the one being created)
void createDir( const String &path ) {
    DIR *dir = opendir( path.c_str() );
    if (dir == NULL) {
#ifdef WIN32
        String pathClean = path.replace( '/', '\\' );
        runSystemCommand( String( "mkdir " ) + pathClean );
#else
        runSystemCommand( String( "mkdir " ) + path );
#endif
    } else {
        closedir( dir );
    }
}


/// maintain a fixed number of files with the given extension and path by removing the oldest matching files 
void removeOldestFiles( const String &path, const String &extension, int keepCount ) {
    Array<String> dirList = dirFileList( path, "", extension );
    if (dirList.count() > keepCount) {

        // get file modification timestamp for each file
        VectorI fileMod( dirList.count() );
        for (int i = 0; i < dirList.count(); i++) 
            fileMod[ i ] = fileModificationTimestamp( path + "/" + dirList[ i ] );

        // sort by file modification; oldest first
        VectorI sortInd = sortIndex( fileMod );
        int removeCount = dirList.count() - keepCount;
        for (int i = 0; i < removeCount; i++) {
            int index = sortInd[ i ];
            //disp( 1, "remove: %s, mod: %s", dirList[ index ].c_str(), TsToStrLocal( fileMod[ index ] ).c_str() );
            deleteFile( path + "/" + dirList[ index ] );
        }
        disp( 3, "removed %d files from %s", removeCount, path.c_str() );
    }
}


/// copy the file to a backup file 
void backUpFile( const String &fileName ) {
    String noExt = fileName.leftOfLast( '.' );
    String ext = fileName.rightOfLast( '.' );
    String newFileName;
    int index = 1;
    do {
        newFileName = noExt + sprintF( "04d", index ) + ext;
        index++;
    } while (fileExists( newFileName ));
    copyFile( fileName, newFileName );
}


/// returns true if the filename has an absolute (non-relative) path
// fix(clean): use better mechanism
bool isAbsoluteFileName( const String &fileName ) {
#ifdef WIN32
    return fileName.contains( ':' );
#else
    return fileName.startsWith( "/" );
#endif
}


/// the timestamp of the last modification to the given file
int fileModificationTimestamp( const String &fileName ) {
    struct stat fileInfo;
    if (stat( fileName.c_str(), &fileInfo) == -1) {
        warning( "fileModificationTimestamp: stat error" );
        return 0;
    }
    return (int) fileInfo.st_mtime;
}


//-------------------------------------------
// OPERATING SYSTEM UTILITIES
//-------------------------------------------


/// run a shell command
void runSystemCommand( const String &command ) {
    String commandCopy( command );
#ifdef WIN32
    commandCopy.replaceInPlace( '/', '\\' );
#endif
    if (system( commandCopy.c_str() ) != 0) {
        warning( "system call returned non-zero" );
    }
}


} // end namespace sbl

