// Licensed under MIT license; see license.txt.

#include <sbl/core/Display.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
namespace sbl {


// number of characters per indentation level
const int g_indentSize = 4;


// don't show or log disp() messages with indent greater than this threshold
int g_maxIndent = 100;


// various possible destinations for output
FILE *g_dispFile1 = stdout;
FILE *g_dispFile2 = NULL;
FILE *g_errorFile1 = stderr;
FILE *g_errorFile2 = NULL;
void (*g_dispCallback)( const char *str ) = NULL;
void (*g_statusCallback)( const char *str ) = NULL;
void (*g_errorCallback)( const char *str ) = NULL;
void (*g_progressCallback)( int index, int count ) = NULL;


// display and log a string 
void dispStr( const char *str ) {
    if (g_dispFile1)
        fprintf( g_dispFile1, "%s\n", str );
    if (g_dispFile2)
        fprintf( g_dispFile2, "%s\n", str );
    if (g_dispCallback)
        g_dispCallback( str );
}


// display and log an error/warning message
void dispErrorStr( const char *str ) {
    if (g_errorFile1)
        fprintf( g_errorFile1, "%s\n", str );
    if (g_errorFile2)
        fprintf( g_errorFile2, "%s\n", str );
    if (g_errorCallback)
        g_errorCallback( str );
}


/// display text, with printf-style formatting
void disp( int indent, const char *str, ... ) {

    // check indent level
    if (indent > g_maxIndent)
        return;

    // do printf
    char displayBuf[ 10000 ];
    va_list argList;
    va_start( argList, str );
    vsprintf( displayBuf, str, argList ); 

    // if indented
    if (indent) {
        int indentation = indent * g_indentSize;
        char *indentStr = new char [strlen( displayBuf ) + indentation + 1];
        for (int i = 0; i < indentation; i++)
            indentStr[i] = ' ';
        strcpy( indentStr + indentation, displayBuf );
        dispStr( indentStr );
        delete [] indentStr;

    // if not indented
    } else {
        dispStr( displayBuf );
    }
}



/// display a short-lived progress message (not logged)
void status( const char *str, ... ) {

    // do printf
    char displayBuf[ 10000 ];
    va_list argList;
    va_start( argList, str );
    vsprintf( displayBuf, str, argList ); 

    // display the status message
    if (g_statusCallback) {
        g_statusCallback( displayBuf );
    } else {
        printf( "\r%s", displayBuf );
    }
}


/// display a warning, with printf-style formatting
void warning( const char *str, ... ) {
    char displayBuf[ 10000 ];

    // do printf
    va_list argList;
    va_start( argList, str );
    vsprintf( displayBuf, str, argList ); 

    // add warning text
    char fullDisplayBuf[ 10000 ];
    sprintf( fullDisplayBuf, "warning: %s", displayBuf );
    dispErrorStr( fullDisplayBuf );
}


/// display a fatal error, with printf-style formatting (and kill the program)
void fatalError( const char *str, ... ) {
    char displayBuf[ 10000 ];

    // do printf
    va_list argList;
    va_start( argList, str );
    vsprintf( displayBuf, str, argList ); 

    // add warning text
    char fullDisplayBuf[ 10000 ];
    sprintf( fullDisplayBuf, "fatal error: %s", displayBuf );
    dispErrorStr( fullDisplayBuf );

    // kill program
    exit( 1 );
}


/// displays progress bar (assumes index in [0, count - 1]);
/// closes progress bar when index == count - 1; 
/// if count == -1, assumes unknown number if items
void progress( int index, int count ) {
    if (g_progressCallback)
        g_progressCallback( index, count );
}


/// don't show or log disp() messages with indent greater than this threshold
void setMaxDisplayedIndent( int maxIndent ) {
    g_maxIndent = maxIndent;
}


/// set a custom handler for disp() messages
void setDispCallback( void (*dispCallback)( const char *str ) ) {
    g_dispCallback = dispCallback;
}


/// set a custom handler for disp() messages
void setStatusCallback( void (*statusCallback)( const char *str ) ) {
    g_statusCallback = statusCallback;
}


/// set a custom handler for warning() and fatalError() messages
void setErrorCallback( void (*errorCallback)( const char *str ) ) {
    g_errorCallback = errorCallback;
}


/// set a custom handler for progress() calls
void setProgressCallback( void (*progressCallback)( int index, int count ) ) {
    g_progressCallback = progressCallback;
}


// get a file handle for output (stdout, stderr, or actual file)
FILE *outputFile( const char *fileName ) {
    FILE *file = NULL;
    if (fileName) {
        if (strcmp( fileName, "stdout" ) == 0) {
            file = stdout;
        } else if (strcmp( fileName, "stderr" ) == 0) {
            file = stderr;
        } else if (fileName[ 0 ]) {
            file = fopen( fileName, "w" );
        }
    }
    return file;
}


/// set output files for disp() messages; can use "stdout" / "stderr"
void setDispFileNames( const char *fileName1, const char *fileName2 ) {

    // close old files
    if (g_dispFile1 && g_dispFile1 != stdout && g_dispFile1 != stderr)
        fclose( g_dispFile1 );
    if (g_dispFile2 && g_dispFile2 != stdout && g_dispFile2 != stderr)
        fclose( g_dispFile2 );

    // open new files
    g_dispFile1 = outputFile( fileName1 );
    g_dispFile2 = outputFile( fileName2 );
}


/// set output files for warning() and fatalError() messages
void setErrorFileNames( const char *fileName1, const char *fileName2 ) {

    // close old files
    if (g_errorFile1 && g_errorFile1 != stdout && g_errorFile1 != stderr)
        fclose( g_errorFile1 );
    if (g_errorFile2 && g_errorFile2 != stdout && g_errorFile2 != stderr)
        fclose( g_errorFile2 );

    // open new files
    g_errorFile1 = outputFile( fileName1 );
    g_errorFile2 = outputFile( fileName2 );
}


} // end namespace sbl

