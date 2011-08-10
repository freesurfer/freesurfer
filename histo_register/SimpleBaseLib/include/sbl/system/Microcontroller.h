#ifndef _SBL_MICROCONTROLLER_H_
#define _SBL_MICROCONTROLLER_H_
#include <sbl/math/Vector.h>
#include <sbl/system/SerialPort.h>
namespace sbl {


/*! \file Microcontroller.h
    \brief The Microcontroller module gives an interface for sending/receiving messages
    from a microcontroller.  The messages are transmitted with a simple binary protocol
    that includes basic error checking.
*/


// register commands, etc. defined in this module
void initMicrocontroller();


/// read a message from the serial port using a simple message protocol
bool readSerialMessage( SerialPort &serialPort, String &messageName, VectorI &values );


/// write a message to the serial port using a simple message protocol
void writeSerialMessage( SerialPort &serialPort, const String &messageName, int valueCount = 0, int value0 = 0, int value1 = 0, int value2 = 0 );


/// request sensor values from microcontroller;
/// returns empty vector on error
VectorI readSensorValues( SerialPort &serialPort );


} // end namespace sbl
#endif // _SBL_MICROCONTROLLER_H_

