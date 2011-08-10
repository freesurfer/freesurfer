#ifndef _SBL_DISPLAY_H_
#define _SBL_DISPLAY_H_
namespace sbl {


/*! \file Display.h
    \brief The Display module provides functions for outputting text from commands
    (including logging and error messages).  The functions can be used with the GUI 
    modules to display diagnostic information in the user interface.
*/


/// display text, with printf-style formatting
void disp( int indent, const char *str, ... );


/// display a short-lived progress message (not logged)
void status( const char *str, ... );


/// display a warning, with printf-style formatting
void warning( const char *str, ... );


/// display a fatal error, with printf-style formatting (and kill the program)
void fatalError( const char *str, ... );


/// displays progress bar (assumes index in [0, count - 1]);
/// closes progress bar when index == count - 1; 
/// if count == -1, assumes unknown number if items
void progress( int index, int count = -1 );
inline void progressDone() { progress( 0, 1 ); }


/// don't show or log disp() messages with indent greater than this threshold
void setMaxDisplayedIndent( int maxIndent );


/// a debug-only assertion
#ifdef _DEBUG
#define assertDebug( condition ) if (!(condition)) fatalError( "assert failed at %s:%d: %s", __FILE__, __LINE__, #condition );
#else
#define assertDebug( condition )
#endif 


/// an assertion included in both debug and release mode
#define assertAlways( condition ) if (!(condition)) fatalError( "assert failed at %s:%d", __FILE__, __LINE__, #condition );


/// set a custom handler for disp() messages
void setDispCallback( void (*dispCallback)( const char *str ) );


/// set a custom handler for status() messages
void setStatusCallback( void (*statusCallback)( const char *str ) );


/// set a custom handler for warning() and fatalError() messages
void setErrorCallback( void (*errorCallback)( const char *str ) );


/// set a custom handler for progress() calls
void setProgressCallback( void (*progressCallback)( int index, int count ) );


/// set output files for disp() messages; can use "stdout" / "stderr"
void setDispFileNames( const char *fileName1, const char *fileName2 );


/// set output files for warning() and fatalError() messages; can use "stdout" / "stderr"
void setErrorFileNames( const char *fileName1, const char *fileName2 );


/// a generic interface to be implemented by a class that displays type T
template<typename T> class Display {
public:

    // avoid non-virtual distructor warnings
    virtual ~Display() {}

    /// implement this to display and object
    virtual void display( T obj ) = 0;
};


} // end namespace sbl
#endif // _SBL_DISPLAY_H_

