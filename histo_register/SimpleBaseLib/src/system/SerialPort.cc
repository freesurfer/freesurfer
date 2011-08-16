// Licensed under MIT license; see license.txt.

#include <sbl/system/SerialPort.h>
#include <sbl/core/Display.h>
#include <sbl/system/Timer.h>
#include <stdio.h>
#ifdef WIN32
    #include <Windows.h>
#endif 
namespace sbl {


//-------------------------------------------
// SERIAL PORT CLASS
//-------------------------------------------


/// open serial port; for windows, serial port name is of form "COM2"
SerialPort::SerialPort( const String &comPortName, int baudRate ) {
    m_comOpened = false;
    m_verbose = false;

    // allocate a buffer for incoming data
    m_bufferLen = 1000;
    m_buffer = new unsigned char[ m_bufferLen ];
    m_bufferPos = 0;

// windows-specific serial port opening code
#ifdef WIN32

    // open comm port file handle
    m_comFile = CreateFileA( comPortName.c_str(), 
        GENERIC_READ | GENERIC_WRITE,
        0, // exclusive access
        NULL, // no security
        OPEN_EXISTING,
        0, // no overlapped I/O
        NULL); // null template 
    if (m_comFile == INVALID_HANDLE_VALUE) {
        warning( "error: failed to open comm file" );
        return;
    }

    // create and clear comm buffers
    if (SetupComm( m_comFile, 128, 128) == false) 
        warning( "error: failed to set up comm buffers" );
    if (PurgeComm( m_comFile, PURGE_TXABORT | PURGE_TXCLEAR ) == false)
        warning( "error: failed to purge comm buffers" );

    // set timeout values (msecs)
    // timeout of 0 -> non-blocking
    COMMTIMEOUTS timeouts;
    timeouts.ReadIntervalTimeout = MAXDWORD;
    timeouts.ReadTotalTimeoutConstant = 0;
    timeouts.ReadTotalTimeoutMultiplier = 0;
    timeouts.WriteTotalTimeoutConstant = 0;
    timeouts.WriteTotalTimeoutMultiplier = 0;
    if (SetCommTimeouts( m_comFile, &timeouts ) == false) 
        warning( "error: error setting comm timeouts" );

    // configure comm port
    DCB dcb;
    if (GetCommState( m_comFile, &dcb )) {
        dcb.BaudRate = baudRate;
        dcb.ByteSize = 8;
        dcb.Parity = NOPARITY;
        dcb.StopBits = ONESTOPBIT;
        if (SetCommState( m_comFile, &dcb) == false)
            warning( "error: error setting comm settings" );
    } else {
        warning( "error: error getting comm settings" );
    }
    m_comOpened = true;
#endif // WIN32
}


// basic destructor
SerialPort::~SerialPort() {
    delete [] m_buffer;
    if (m_comOpened) {
#ifdef WIN32
        if (PurgeComm( m_comFile, PURGE_TXCLEAR | PURGE_RXCLEAR ) && m_verbose)
            warning( "error: failed to purge comm buffers" );
        CloseHandle( m_comFile );
#endif // WIN32
    }
}


/// read a byte; non-blocking; returns -1 if no byte to read 
int SerialPort::readByte() {
    unsigned char byte = 0;
    unsigned long readCount = (unsigned long)-1;
#ifdef WIN32
    ReadFile( m_comFile, &byte, 1, &readCount, NULL );
#endif // WIN32
    int val = byte;
    if (readCount == 0)
        val = -1;
    return val;
}


/// send a byte to the serial port
void SerialPort::writeByte( unsigned char byte, bool debug ) {
#ifdef WIN32
    //unsigned long error = 0;
    //ClearCommError( m_comFile, &error, NULL); // fix(clean): remove this?
    unsigned long writeCount = -1;
    WriteFile( m_comFile, &byte, 1, &writeCount, NULL );
#endif // WIN32
    if (debug) {
        if (byte >= 'A' && byte <= 'Z')
            disp( 1, "[%c]", byte );
        else
            disp( 1, "(%d)", byte );
    }
}


/// read a string up to the specified stop symbol
String SerialPort::readString( int stopByte ) {
    m_bufferPos = 0;
    readToBuffer( stopByte );
    m_buffer[ m_bufferPos ] = 0;
    return String( (char *) m_buffer );
}


/// write a string (assumes 8-bit char values)
void SerialPort::writeString( const String &str ) {
    for (int i = 0; i < str.length(); i++) 
        writeByte( str.c_str()[ i ] );
}


/// append read data onto internal buffer, until reaching the specified stop byte
void SerialPort::readToBuffer( int stopByte ) {
    int byte = readByte();
    while (byte >= 0) {
        m_buffer[ m_bufferPos++ ] = byte;
        if (m_bufferPos >= m_bufferLen) // fix(later): be smarter about this
            m_bufferPos = 0;
        if (byte == stopByte)
            break;
        byte = readByte(); // returns -1 if no byte to read; breaks loop
    }
    assertAlways( m_bufferPos < m_bufferLen );
}


} // end namespace sbl

