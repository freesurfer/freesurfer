// Licensed under MIT license; see license.txt.

#include <sbl/system/Microcontroller.h>
#include <sbl/core/Command.h>
#include <sbl/other/Plot.h>
#include <sbl/system/SerialPort.h>
#include <sbl/system/Timer.h>
#ifdef USE_GUI
    #include <sbl/gui/MainWindow.h>
#endif
namespace sbl {


//-------------------------------------------
// SIMPLE MESSAGE PROTOCOL
//-------------------------------------------


// constants for simple message protocol
#define PROTOCOL_VERSION '0'
#define PROTOCOL_TERMINATOR 255


/* 

standard message names:

"er" echo request
"dr" device info request
"ep" echo reply
"di" device info
"sr" request sensor values
"ss" stream sensor values
"sv" sensor values

*/


// compute a checksum of the buffer; output in [0, 254]
unsigned char computeCheckSum( const unsigned char *buf, int len ) {
    int checkSum = 0;
    for (int i = 0; i < len; i++)
        checkSum ^= buf[ i ];
    if (checkSum == 255) 
        checkSum = 254;
    return checkSum;
}


/// read a message from the serial port using a simple message protocol
bool readSerialMessage( SerialPort &serialPort, String &messageName, VectorI &values ) {

    // check for end-of-message
    bool messageGood = false;
    serialPort.readToBuffer( PROTOCOL_TERMINATOR );
    const unsigned char *message = serialPort.buffer();
    int length = serialPort.bufferLen();
    if (length && message[ length - 1 ] == PROTOCOL_TERMINATOR) {

        // remove data from serial port buffer
        serialPort.resetBuffer();

        // check message header/footer
        unsigned char protocolVersion = message[ 0 ];
        unsigned char payloadLength = message[ 3 ];
        unsigned char localCheckSum = computeCheckSum( message, length - 2 );
        unsigned char remoteCheckSum = message[ length - 2 ];
        unsigned char terminator = message[ length - 1 ];
        if (protocolVersion != PROTOCOL_VERSION) {
            warning( "invalid protocol version: %d (expected: %d)", protocolVersion, PROTOCOL_VERSION );
        } else if (terminator != 255) {
            warning( "invalid terminatore: %d", terminator );
        } else if (localCheckSum != remoteCheckSum) {
            warning( "invalid check sum: %d (remote: %d)", localCheckSum, remoteCheckSum );
        } else if (payloadLength != length - 6) {
            warning( "invalid payload length: %d (message length: %d)", payloadLength, length );
        } else {

            // extract the message payload
            messageName = String( 0, 3 );
            messageName.set( 0, message[ 1 ] );
            messageName.set( 1, message[ 2 ] );
            int valueCount = payloadLength / 2;
            int pos = 4;
            for (int i = 0; i < valueCount; i++) {
                int vLower = message[ pos++ ];
                int vUpper = message[ pos++ ];
                values.append( vLower + (vUpper << 7) );
            }
            messageGood = true;
        }
    }
    return messageGood;
}


/// write a message to the serial port using a simple message protocol
void writeSerialMessage( SerialPort &serialPort, const String &messageName, int valueCount, int value0, int value1, int value2 ) {
    assertAlways( messageName.length() == 2 );
    assertAlways( valueCount >= 0 && valueCount < 3 );
    assertAlways( serialPort.openSuccess() );

    // will hold message as constructed
    unsigned char buf[ 20 ];
    int pos = 0;

    // add header
    buf[ pos++ ] = PROTOCOL_VERSION;
    buf[ pos++ ] = (unsigned char) messageName.get( 0 );
    buf[ pos++ ] = (unsigned char) messageName.get( 1 );
    buf[ pos++ ] = valueCount * 2;

    // add values
    if (valueCount >= 1) {
        buf[ pos++ ] = value0 & 127;
        buf[ pos++ ] = value0 >> 7;
    }
    if (valueCount >= 2) {
        buf[ pos++ ] = value1 & 127;
        buf[ pos++ ] = value1 >> 7;
    }
    if (valueCount >= 3) {
        buf[ pos++ ] = value2 & 127;
        buf[ pos++ ] = value2 >> 7;
    }

    // add checksum and terminator
    unsigned checkSum = computeCheckSum( buf, pos );
    buf[ pos++ ] = checkSum;
    buf[ pos++ ] = 255;

    // send the message
    for (int i = 0; i < pos; i++) 
        serialPort.writeByte( buf[ i ] );
}


//-------------------------------------------
// READ SENSOR VALUES
//-------------------------------------------


/// request sensor values from microcontroller;
/// returns empty vector on error
VectorI readSensorValues( SerialPort &serialPort ) {

    // send sensor request message
    writeSerialMessage( serialPort, "sr" );

    // wait a moment for reply
    delaySeconds( 0.1f );

    // read message
    String messageName;
    VectorI values;
    readSerialMessage( serialPort, messageName, values );
    if (values.length())
        return subset( values, 1, values.length() - 1 ); // exclude first value, since it is the time diff
    return values;
}


// read values from sensors
void readSerialSensors( Config &conf ) {

    // get command parameters
    const String &portName = conf.readString( "portName" );
    int baudRate = conf.readInt( "baudRate", 9600 );
    if (conf.initialPass())
        return;

    // open serial port
    SerialPort serialPort( portName, baudRate );
    if (serialPort.openSuccess() == false)
        return;

    // read and display sensor values
    VectorI values = readSensorValues( serialPort );
    if (values.length()) {
        disp( 1, "values: %s", toString( values ).c_str() );
    } else {
        disp( 1, "no values" );
    }
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initMicrocontroller() {
    registerCommand( "ssread", readSerialSensors );
}


} // end namespace sbl

