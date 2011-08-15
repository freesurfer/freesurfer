#ifndef _SBL_SERIAL_PORT_H_
#define _SBL_SERIAL_PORT_H_
#include <sbl/core/String.h>
#include <sbl/math/Vector.h>
namespace sbl {


/// The SerialPort class provides a simple platform-independent interface for serial devices.
class SerialPort {
public:

    /// open serial port; for windows, serial port name is of form "COM2"
    explicit SerialPort( const String &comPortName, int baudRate );

    // basic destructor
    ~SerialPort();

    // true if serial port opened successfully
    inline bool openSuccess() const { return m_comOpened; }

    //-------------------------------------------
    // READ / WRITE BYTES
    //-------------------------------------------

    /// read a byte; non-blocking; returns -1 if no byte to read
    int readByte();

    /// send a byte to the serial port
    void writeByte( unsigned char byte, bool debug = false );

    //-------------------------------------------
    // READ / WRITE STRINGS
    //-------------------------------------------

    /// read a string up to the specified stop symbol
    String readString( int stopByte );

    /// write a string (assumes 8-bit char values)
    void writeString( const String &str );

    //-------------------------------------------
    // READ TO BUFFER
    //-------------------------------------------

    /// append read data onto internal buffer, until reaching the specified stop byte
    void readToBuffer( int stopByte );

    inline const unsigned char *buffer() const { return m_buffer; }

    inline int bufferLen() const { return m_bufferPos; }

    /// reset/clear the buffer
    inline void resetBuffer() { m_bufferPos = 0; }

private:

    // the serial port connection
    bool m_comOpened;
    void *m_comFile;

    // configuration
    bool m_verbose;

    // a buffer of received data
    unsigned char *m_buffer;
    int m_bufferPos;
    int m_bufferLen;

    // disable copy constructor and assignment operator
    SerialPort( const SerialPort &x );
    SerialPort &operator=( const SerialPort &x );
};


} // end namespace sbl
#endif // _SBL_SERIAL_PORT_H_

