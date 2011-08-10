#ifndef _SBL_SIGNAL_H_
#define _SBL_SIGNAL_H_
namespace sbl {


/*! \file Signal.h
    \brief The Signal module is used to capture CTRL-C keypresses.  
    After initSignal() is called, a single CTRL-C will attempt to cancel the 
    currently running command and a triple CTRL-C will terminate the program.
*/


// register commands, etc. defined in this module
void initSignal();


} // end namespace sbl
#endif // _SBL_SIGNAL_H_

