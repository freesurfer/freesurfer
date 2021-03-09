/**
 * @brief simple "messaging" class
 *
 *  `c_SMessage' is a simple "messaging" class that contains a string
 *  payload (the message body), an optional stream specifier and an
 *  optional formatting enumeration.
 *
 *  Typically, c_SMessage is embedded within other objects, and allows
 *  a convenient way to encapsulate string-type data that might be
 *  ultimately displayed in a terminal or parsed by some GUI display
 *  method.
 */
/*
 * Original Author: Rudolph Pienaar
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef __C_SMESSAGE_H__
#define __C_SMESSAGE_H__

#include <iostream>
#include <fstream>
#include <string>
//#include <cstdio>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "c_SSocket.h"

typedef enum {
    eSM_devnull, eSM_raw, eSM_column
} e_SMessageFormat;

typedef enum {
    eSM_cpp,      // C++ style I/O
    eSM_c,        // C style I/O
    eSS           // Simple Socket style I/O
} e_SMessageIO;

typedef enum {
    eAppend, eOverwrite
} e_FileMode;

typedef enum {
  eSM_start, eSM_stop
} e_SMTimerAction;

const int SMessage_STACKDEPTH = 64;

class C_SMessage {

    // data structures

    protected:
    //
    // generic object structures - used for internal bookkeeping
    // and debugging / automated tracing methods. The stackDepth
    // and str_proc[] variables are maintained by the debug_push|pop
    // methods
    //
    string  str_obj;                    // name of object class
    string  str_name;                   // name of object variable
    int     id;                         // id of agent
    int     iter;                       // current iteration in an
                                        //+ arbitrary processing scheme
    int     verbosity;                  // debug related value for object
    int     warnings;                   // show warnings (and warnings level)
    int     stackDepth;                 // current pseudo stack depth
    string  str_proc[SMessage_STACKDEPTH];
    // execution procedure stack

    //
    // actual data pertinent to this class
    //
    bool                b_canPrint;     // overall control flag: if FALSE
                                        //+ no output generated on called
                                        //+ functions. Essentially a 'kill 
                                        //+ switch'.
    bool                b_syslogPrepend;// default syslog prepend control --
                                        //+ can be overriden
    e_SMessageFormat    e_format;       // formatting enum of the body
    string              str_payload;    // the body of the object
    string              str_syslogID;   // a string used for syslogIDing
    e_SMessageIO        e_IO;           // IO type - either C style or
                                        // C++ style

    // column width specifiers: left, center, and right
    //+ in most cases, lw and rw are used.
    int                 lw;
    int                 cw;
    int                 rw;

    // output file concerns
    string              str_filename;   // output filename
    FILE*               pFILE_out;      // C-style (console) stream
    std::ofstream       ofs_out;        // File stream for C++ output
    bool                b_fileSpecified;// bool flag that is used to
                                        //+ determine how to stream

    // output UDP socket (if spec'd)
    bool                b_socketCreated;
    c_SSocket_UDP_transmit*  pcSS;

    // methods

    public:
    //
    // constructor / destructor block
    //
    void    core_construct( string  astr_name       = "unnamed",
                            int     a_id            = -1,
                            int     a_iter          = 0,
                            int     a_verbosity     = 0,
                            int     a_warnings      = 0,
                            int     a_stackDepth    = 0,
                            string  astr_proc       = "noproc");

    C_SMessage(
        string                  astr_body       = "",
        e_SMessageFormat        ae_format       = eSM_raw,
        FILE*                   apFILE_out      = stdout,
        e_SMessageIO            ae_IO           = eSM_c
        );

    C_SMessage(
        string                  astr_body       = "",
        e_SMessageFormat        ae_format       = eSM_raw,
        string                  astr_filename   = "stdout",
        e_SMessageIO            ae_IO           = eSM_cpp,
        e_FileMode              ae_fileMode     = eAppend
        );
    ~C_SMessage();

    //
    // error / warn / print block
    //
    void        debug_push(     string          astr_currentProc);
    void        debug_pop();

    void        error(
        string                  astr_msg        = "Some error has occurred",
        int                     code            = -1
        );
    void        warn(
        string                  astr_msg        = "Some warning has occurred",
        int                     code            = -1
        );
    void        function_trace( 
        string                  astr_msg,
        string                  astr_separator
        );

    //
    // access block
    //
    void    print();                        // print object
    int     stackDepth_get()            const   {return stackDepth;};
    void    stackDepth_set(int anum)            {stackDepth = anum;};
    int     iter_get()                  const   {return iter;};
    void    iter_set(int anum)                  {iter = anum;};
    int     id_get()                    const   {return id;};
    void    id_set(int anum)                    {id = anum;};
    int     verbosity_get()             const   {return verbosity;};
    void    verbosity_set(int anum)             {verbosity = anum;};
    int     warnings_get()              const   {return warnings;};
    void    warnings_set(int anum)              {warnings = anum;};
    string  str_obj_get()               const   {return str_obj;};
    void    str_obj_set(string astr)            {str_obj = astr;};
    string  str_name_get()              const   {return str_name;};
    void    str_name_set(string astr)           {str_name = astr;};
    string  str_proc_get()              const   {return str_proc[stackDepth_get()];};
    void    str_proc_set(int depth, string astr) {str_proc[depth] = astr;};
    string str_filename_get()           const   {return str_filename;};
    FILE*   pFILE_out_get()             const   {return pFILE_out;};
    void    pFILE_out_set(FILE* apFILE)         {pFILE_out = apFILE;};
    string str_syslogID_get()           const   {return str_syslogID;};
    void str_syslogID_set(string astr)          {str_syslogID = astr;};
    string  str_payload_get()           const   {return str_payload;};
    void    str_payload_set(string astr)        {str_payload = astr;};
    void    str_payload_clear()                 {str_payload = "";};
    void    str_payload_append(string astr)     {str_payload.append(astr);};
    void    str_payload_prepend(string astr)    {str_payload.insert(0, astr);};
    e_SMessageFormat e_format_get()     const   {return e_format;};
    void    e_format_set(e_SMessageFormat ae_format) {e_format  = ae_format;};
    e_SMessageIO e_IO_get()             const   {return e_IO;};
    void e_IO_set(e_SMessageIO ae_IO)           {e_IO  = ae_IO;};

    void lw_set(int aval) {lw = aval;};
    void rw_set(int aval) {rw = aval;};
    void cw_set(int aval) {cw = aval;};
    void b_canPrint_set(        bool    abval)  {b_canPrint = abval;};
    void b_syslogPrepend_set(   bool    abval)  {b_syslogPrepend = abval;};
    
    //
    // miscellaneous block - the main functionality provided by
    // this class
    //
    int         printf(         const char*             format, ...);
    int         lprintf(        const char*             format, ...);
    int         colprintf(
                                const char*             pch_lstr,
                                const char*             format, ...
        );
    int         pprintf(        const char*             format, ...);
    int         plprintf(       const char*             format, ...);
    int         pcolprintf(
                                const char*             pch_lstr,
                                const char*             format, ...
        );
    int         lprintf(        string&                 str_bufferOut,
                                const char*             format, ...);
    int         colprintf(
                                string&                 str_bufferOut,
                                const char*             pch_lstr,
                                const char*             format, ...
        );
    bool        timer(          e_SMTimerAction         e_action);
    string      syslog_prepend();
    void        dump(           bool                    ab_syslogPrepend = false,
                                string                  astr_outOfBand   = ""
                );
    bool        file_changeTo(  string                  astr_filename,
                                e_SMessageIO            ae_IO           = eSM_cpp);
};

#endif //__C_SMESSAGE_H__
