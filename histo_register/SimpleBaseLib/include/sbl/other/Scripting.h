#ifndef _SBL_SCRIPTING_H_
#define _SBL_SCRIPTING_H_
namespace sbl {


/*! \file Scripting.h
    \brief The Scripting module provides a simple Python scripting interface.
    Python scripts can use sbl.py to execute commands within an SBL-based C++ program.
*/


// register commands, etc. defined in this module
void initScripting();


} // end namespace sbl
#endif // _SBL_SCRIPTING_H_

