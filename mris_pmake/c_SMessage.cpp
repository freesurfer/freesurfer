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


#include <iostream>
#include <sstream>
using namespace std;

#include <sys/times.h>
#include <time.h>

// #include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <c_SMessage.h>

//
//\\\***
// C_SMessage definitions ****>>>>
/////***
//

void C_SMessage::debug_push( string astr_currentProc) {
  //
  // ARGS
  // astr_currentProc         in      method name to
  //                                  "push" on the "stack"
  //
  // DESC
  // This attempts to keep a simple record of methods that
  // are called. Note that this "stack" is severely crippled in
  // that it has no "memory" - names pushed on overwrite those
  // currently there.
  //

  if (stackDepth_get() >= SMessage_STACKDEPTH-1)
    error("Out of str_proc stack depth");
  stackDepth_set(stackDepth_get()+1);
  str_proc_set(stackDepth_get(), astr_currentProc);
}

void C_SMessage::debug_pop() {
  //
  // DESC
  // "pop" the stack. Since the previous name has been
  // overwritten, there is no restoration, per se. The
  // only important parameter really is the stackDepth.
  //

  stackDepth_set(stackDepth_get()-1);
}

void C_SMessage::error( string astr_msg,
                        int  code) {
  //
  // ARGS
  // atr_msg          in              message to dump to stderr
  // code             in              error code
  //
  // DESC
  // Print error related information. This routine throws an exception
  // to the class itself, allowing for coarse grained, but simple
  // error flagging.
  //

  cerr << "\nFatal error encountered.\n";
  cerr << "\tC_SMessage object `"       << str_name;
  cerr << "' (id: " << id << ")\n";
  cerr << "\tCurrent function: "        << str_obj      << "::";
  cerr << str_proc_get() << "\n";
  cerr << "\t" << astr_msg << "\n";
  cerr << "Throwing an exception to (this) with code " << code << "\n\n";
  throw(this);
}

void C_SMessage::warn( string astr_msg,
                       int  code) {
  //
  // ARGS
  // atr_msg          in              message to dump to stderr
  // code             in              error code
  //
  // DESC
  // Print error related information. Conceptually identical to
  // the `error' method, but no expection is thrown.
  //

  cerr << "\nWarning.\n";
  cerr << "\tC_SMessage object `"       << str_name;
  cerr << "' (id: " << id << ")\n";
  cerr << "\tCurrent function: "        << str_obj      << "::";
  cerr << str_proc_get() << "\n";
  cerr << "\t" << astr_msg << "(code: " << code << ")\n";
}

void C_SMessage::function_trace(string astr_msg,
                                string astr_separator) {
  //
  // ARGS
  // astr_msg         in      message to print
  // astr_separator   in      separator char between successive calls
  //
  // DESC
  // This method allows for run-time debugging-related information
  // in a particular class and method to be displayed to stderr.
  //

  string str_tab                      = "";
  static string str_objectName        = "";
  static string str_funcName          = "";

  if (verbosity_get() >= stackDepth_get()) {
    cerr << astr_separator;
    for (int i = 0; i < stackDepth_get(); i++)
      str_tab += "\t";
    if (str_objectName != str_name_get() )
      cerr << "\nC_SMessage`" << str_name_get();
    cerr << "' (id: " << id_get() << ")\n";
    if (str_funcName != str_proc_get()) {
      cerr << "\n" << str_tab << "Current function: " << str_obj;
      cerr << "::" << str_proc_get() << endl;
      cerr << "\tverbosity = "    << verbosity_get();
      cerr << ", stackDepth = "   << stackDepth_get() << endl;
    }
    cerr << "\n" << str_tab << astr_msg;
  }
  str_objectName      = str_name_get();
  str_funcName        = str_proc_get();
}

void    C_SMessage::core_construct(     string  astr_name,
                                        int     a_id,
                                        int     a_iter,
                                        int     a_verbosity,
                                        int     a_warnings,
                                        int     a_stackDepth,
                                        string  astr_proc) 
{
    //
    // ARGS
    // astr_name        in              name of object
    // a_id             in              id of object
    // a_iter           in              current iteration in arbitrary scheme
    // a_verbosity      in              verbosity of object
    // a_stackDepth     in              stackDepth
    // astr_proc        in              current proc that has been "debug_push"ed
    //
    // DESC
    // Simply fill in the core values of the object with some defaults
    //
    // HISTORY
    //  o Initial design and coding
    //
    // 24 September 2001
    //  o Added syslogID
    //

    str_obj                     = "C_SMessage";
    str_name                    = astr_name;
    id                          = a_id;
    iter                        = a_iter;
    verbosity                   = a_verbosity;
    warnings                    = a_warnings;
    stackDepth                  = a_stackDepth;
    str_proc[stackDepth]        = astr_proc;

    b_canPrint                  = true;
    b_syslogPrepend             = false;
    
    str_syslogID                = "";
    b_fileSpecified             = false;
    pFILE_out                   = NULL;

    b_socketCreated             = false;
    pcSS                        = NULL;

    lw                          = 40;
    rw                          = 40;
    cw                          = 40;
}

C_SMessage::C_SMessage(
    string                      astr_body,
    e_SMessageFormat            ae_format,
    FILE*                       apFILE_out,
    e_SMessageIO                ae_IO) {
  //
  // ARGS
  // astr_body          in              initial value of the payload
  // ae_format          in              format descriptor for the payload
  // pFILE_out          in              default stream (stdout) to which
  //                                    + payload is sent
  // ae_IO              in              IO style (C or C++)
  //
  // DESC
  // Constructor
  // Creates and initialises default object to stdout, C-style,
  // raw payload dumping.
  //
  // PRECONDITIONS
  // o This constructor is only for C-style C_SMessage objects!
  //
  // POSTCONDITIONS
  // o C-style C_SMessage object is created.
  //
  // WARNING... WARNING
  // o This constructor is due to become depreciated. It is much safer
  //   to use the main constructor which actually does error checking
  //   on filenames that are passed to it.
  //
  // HISTORY
  // 13 September 2000
  // o Initial design and coding.
  //   C-Style file handling is used since I couldn't find a simple
  //   way to use C++ streams and default parameters and using a unified
  //   interface over cout and text files
  //
  // 27 June 2001
  // o Expanded handling of C++ style IO.
  //
  // 05 April 2005
  // o Added c_SSocket.
  //

  core_construct();

  if (ae_IO != eSM_c)
    error("This constructor can only be used for C-style C_SMessage objects", -1);

  str_payload.assign(astr_body);
  e_format              = ae_format;
  pFILE_out             = apFILE_out;
  e_IO                  = ae_IO;
}


C_SMessage::C_SMessage(
    string                      astr_body,
    e_SMessageFormat            ae_format,
    string                      astr_filename,
    e_SMessageIO                ae_IO,
    e_FileMode                  ae_fileMode
) {
    //
    // ARGS
    //  astr_body               in              initial value of the payload
    //  ae_format               in              payload format descriptor
    //  astr_filename           in              filename in which payload must be
    //                                          + appended. If <ae_IO> is 'eSS'
    //                                          + then filename is assumed to
    //                                          + be in "hostname:port"  format.
    //  ae_IO                   in              IO style (C, C++, or SSocket)
    //  ae_fileMode             in              overwrite/append
    //
    // DESC
    //  Constructor
    //  Creates and initialises object.
    //
    // HISTORY
    // 13 September 2000
    //  o Initial design and coding.
    //    C-Style file handling is used since I couldn't find a simple
    //    way to use C++ streams and default parameters and using a unified
    //    interface over cout and text files
    //
    // 27 June 2001
    //  o Expanded to handle C++ style IO
    //
    // 13 July 2001
    //  o Adopted new C++ approach
    //
    // 15 March 2002
    // o Changed C++ output concerns for g++-3.x. Final referencing for dumping
    //   to stream is now handled exclusively by the dump() method. The constructor
    //   merely opens an output file if so specified and sets a boolean flag
    //   correspondingly.
    //
    // 05 April 2005
    //  o Incorporated c_SSocket UDP transmit class
    //
    // 09 December 2009
    //  o file append/overwrite
    //  

    core_construct();
    debug_push("C_SMessage");

    FILE*                       pFILE_stream;
    str_filename                = astr_filename;
    string                      str_fileMode    = "";
    b_fileSpecified             = false;

    if (str_filename  == "/dev/null")
        e_format        = eSM_devnull;
    else {
        str_payload.assign(astr_body);
        e_format        = ae_format;
        e_IO            = ae_IO;

        switch (e_IO) {
        case eSM_cpp:
            pFILE_out = NULL;
            if ( (str_filename != "stdout") && str_filename != "stderr" ) {
                if(ae_fileMode == eAppend)
                    ofs_out.open(astr_filename.c_str(), ios::app);
                else
                    ofs_out.open(astr_filename.c_str(), ios::out);
                if (!ofs_out)
                    error("Could not create output file:" + astr_filename, -1);
                b_fileSpecified = true;
            }
        break;

        case eSM_c:
            if (astr_filename == "stdout")
                pFILE_out       = stdout;
            else if (astr_filename == "stderr")
                pFILE_out        = stderr;
            else {
                if(ae_fileMode == eAppend)     str_fileMode    = "a";
                else                            str_fileMode    = "w";
                if ( (pFILE_stream = fopen(astr_filename.c_str(),
                                            str_fileMode.c_str())) == NULL) {
                    string  str_error = "Cannot open file " + astr_filename;
                    error(str_error, -1);
            }
            pFILE_out  = pFILE_stream;
            }
        break;

        case eSS:
            string  str_remotehost;         // host to receive UDP transmission
            string  str_remoteport;         // port on remote hostname
            int     pos;                    // position of separating ":"
            int     port;
            pos = astr_filename.find_first_of(":");
            string str_msg =  "For SSocket-type target, destination string ";
            str_msg += "must be '<hostname>:<port>'.\n";
            str_msg += "\tThe ':' char has to be present.";
            if ((unsigned)pos == (unsigned) string::npos)
                error(str_msg);
            str_remotehost = astr_filename.substr(0, pos);
            str_remoteport = astr_filename.substr(pos+1);
            port = atoi(str_remoteport.c_str());
            pcSS = new c_SSocket_UDP_transmit(str_remotehost, port);
            b_socketCreated = true;
            break;
        }
    }
    debug_pop();
}

C_SMessage::~C_SMessage() {
  switch (e_IO) {
  case eSM_c:
    fprintf(pFILE_out, "\n\n");
    fflush(pFILE_out);
    fclose(pFILE_out);
    break;
  case eSM_cpp:
    if ( b_fileSpecified ) {
      ofs_out.flush();
      ofs_out.close();
    }
    break;
  case eSS:
    delete pcSS;
    b_socketCreated = false;
    break;
  }
}

void C_SMessage::dump(
    bool                ab_syslogPrepend,
    string              astr_outOfBand
) {
  //
  // ARGS
  // ab_syslogPrepend   in              bool flag. If true
  //                                    do syslog_prepend
  //                                    with write operation
  // astr_outOfBound    in              if length is non-zero, will
  //                                    dump this string
  //                                    instead of str_payload
  //
  // DESC
  //  "dump" the payload (or the outOfBand string) to the specified file
  // descriptor
  //
  // PRECONDITIONS
  // o str_payload should be defined (unless need to dump outOfBand string)
  //
  // POSTCONDITIONS
  // o if ab_syslogPrepend is true, will call syslogPrepend()
  //
  // HISTORY
  // 13 September 2000
  //  o Initial design and coding.
  //
  // 18 June 2001
  //  o Added return values on fprintf
  //   Seems that in some cases output buffers (if file buffers) are
  //    not updated correctly.
  //
  // 27 June 2001
  //  o Expanded to accommodate both C and C++ style output
  //
  // 13 July 2001
  //  o Adopted new C++ approach
  //
  // 23 September 2001
  // o Added ab_syslogPrepend flag
  // o Added astr_outOfBound parameter
  //
  // 15 March 2002
  // o Changed C++ output design to conform with g++-3.x
  //
  // 25 March 2002
  //  o After almost a week of experimenting with g++-3.0.4 decided
  //    to revert back to gcc-2.95.3. Fortunately, this only
  //    entailed adding a (std::ostream&) typecast to std::cout
  //
  // 05 April 2005
  // o Added c_SSocket.
  //

  string str_syslogPrepend  = "";
  std::ostream& sout  = (b_fileSpecified) ?
                        ofs_out : (str_filename == "stdout") ?
                        (std::ostream&)std::cout :
                        (std::ostream&)std::cerr;

  if (e_format_get() != eSM_devnull) {
    if (ab_syslogPrepend)
      str_syslogPrepend = syslog_prepend();
    switch (e_IO) {
    case eSS:
      if (astr_outOfBand.length())
        pcSS->sendStr(str_syslogPrepend + astr_outOfBand);
      else
        pcSS->sendStr(str_syslogPrepend + str_payload);
      break;
    case eSM_cpp:
      if (astr_outOfBand.length())
        sout << str_syslogPrepend << astr_outOfBand;
      else
        sout << str_syslogPrepend << str_payload;
      sout.flush();
      break;
    case eSM_c:
      if (astr_outOfBand.length())
        fprintf(pFILE_out, "%s%s",
                str_syslogPrepend.c_str(),
                astr_outOfBand.c_str());
      else
        fprintf(pFILE_out, "%s%s",
                str_syslogPrepend.c_str(),
                str_payload.c_str());
      fflush(pFILE_out);
      break;
    }
  }
}

string
C_SMessage::syslog_prepend() {
    //
    // DESC
    // Returns a string that conforms (roughly) to the first few
    // fields of syslog type output.
    //
    // This method is usually called by some object that contains
    // a C_SMessage object. The data field, str_syslogID, is used
    // to convey optional "external" information relevant to the
    // syslog stamp.
    //
    // PRECONDITIONS
    //  o Attempts to emulate POSIX (Linux) type of syslogging.
    //  o Assumes the existence of HOSTNAME env variable.
    //
    // POSTCONDITIONS
    //  o creates a string of the following format:
    //   [$date] [$HOSTNAME] astr_ID
    //  o note that each time this method is called, it will prepend
    //   the str_payload variable with a syslog entry... so don't call
    //   it multiple times on the same payload!
    //
    // HISTORY
    // 24 September 2001
    // o Initial design and coding.
    //
    // 08 December 2005
    // o Note that in some architectures / machines the call to
    //  to getenv("HOSTNAME") fails in strange and bizarre ways!
    //
    // 01 December 2009
    // o getenv("HOSTNAME") borks on ubuntu despite the 'try' block.
    //   Rewrote using gethostname().
    //

    time_t              time_now                = time(NULL);
    stringstream        sstream("");
    string              str_hostname("");
    char                pch_buffer[65536];
    //int                 ret                     = 0;

    try {
        //ret = gethostname(pch_buffer, 65536);
        gethostname(pch_buffer, 65536);
        str_hostname = pch_buffer;
    } catch (...) {
        str_hostname = "";
    }

    // Time stamp issues
    string  str_time = ctime(&time_now);
    // strip out trailing `YYYY\n'
    str_time.erase(str_time.length()-5);
    sstream << str_time << " ";

    sstream <<  str_hostname << " " << str_syslogID << "\t";

    return(sstream.str());
}

bool
C_SMessage::timer(
  e_SMTimerAction  e_timerAction
) {
  //
  // ARGS
  //  e_timerAction     in      enumerated timer
  //                            action to take.
  //
  // DESC
  //  This method does some simple timer-type processing.
  //
  // PRECONDITIONS
  // o Call with `eSM_start' before calling with `eSM_stop'.
  // o This method might be Linux-specific. Check time() calls
  //  for other OSes.
  //
  // POSTCONDITIONS
  // o If called with `e_stop' first, will do nothing and
  //  return `false'; else will return `true'.
  // o This method will do a syslog_prepend()
  //
  // HISTORY
  // 24 September 2001
  // o Initial design and coding.
  //
  // 26 September 2001
  // o Fixed up bug with time issues (used wrong underlying structures).
  //

  static time_t tt_begin, tt_end;
  float        f_processtime;
  stringstream sstream("");

  switch (e_timerAction) {
  case eSM_start:
    time(&tt_begin);
    sstream << "Timer START." << endl;
    dump(true, sstream.str());
    break;
  case eSM_stop:
    time(&tt_end);
    sstream << "Timer STOP. ";
    f_processtime = difftime(tt_end, tt_begin);
    sstream << "Total number of seconds: " << f_processtime << endl;
    dump(true, sstream.str());
    break;
  }
  return true;
}

bool
C_SMessage::file_changeTo(
    string                      astr_filename,
    e_SMessageIO                ae_IO) {
    //
    // ARGS
    //  astr_filename           in              new filename to channel output
    //  ae_IO                   in              type of channel (C or C++)
    //
    // DESC
    //  This method associates the internal file pointer
    //  with the passed filename. It's primary purpose is
    //  to allow a convenient mechanism to re-direct the
    //  stream of an already created C_SMessage object.
    //
    //  This method understands two "special" filenames -
    //  "stdout" and "stderr" which, although passed as
    //  strings, are interpreted as file pointers.
    //
    //  If the <ae_IO> is 'eSS', a socket is assumed. The
    // <astr_filename> is interpreted as '<host>:<port>'.
    //
    // HISTORY
    // 14 September 2000
    //  o Initial design and coding.
    //
    // 18 June 2001
    //  o Added "/dev/null" as special file type (changes e_format and
    //    blocks all write attempts on channel during dump()).
    //
    // 27 June 2001
    //  o Expanded to accommodate C and C++ style output
    //
    // 05 July 2001
    //  o Closing file streams only if != (stdout || stderr)
    //  o C++ style open mode set to ios:app
    //
    // 13 July 2001
    //  o Adopted new C++ approach
    //
    // 15 March 2002
    // o Changed C++ approach for g++-3.x.
    //
    // 05 April 2005
    //  o c_SSocket incorporation.
    //

    debug_push("file_changeTo");

    // Close legacy files
    if ( b_fileSpecified ) {
        if (ofs_out.is_open())
        ofs_out.close();
        if (pFILE_out) {
        fflush(pFILE_out);
        fclose(pFILE_out);
        }
    }
    b_fileSpecified  = false;

    e_IO  = ae_IO;
    str_filename = astr_filename;

    switch (e_IO) {
    case eSM_c:
        FILE*       pFILE_stream;
        if (astr_filename == "stdout")
        pFILE_out       = stdout;
        else if (astr_filename == "stderr")
        pFILE_out       = stderr;
        else if (astr_filename == "/dev/null")
        e_format = eSM_devnull;
        else {
        if ( (pFILE_stream = fopen(astr_filename.c_str(), "a")) == NULL) {
            string  str_error = "Cannot open file " + astr_filename;
            error(str_error);
        }
        pFILE_out           = pFILE_stream;
        }
        break;
    case eSM_cpp:
        if ( str_filename == "/dev/null")
        e_format = eSM_devnull;

        if ( (str_filename != "stdout") && str_filename != "stderr" ) {
        ofs_out.open(astr_filename.c_str(), ios::app);
        if (!ofs_out)
            error("Could not create output file:" + astr_filename, -1);
        b_fileSpecified = true;
        }
        break;
    case eSS:
        string str_remotehost;  // host to receive UDP transmission
        string str_remoteport;  // port on remote hostname
        int  pos;               // position of separating ":"
        int  port;
        pos = astr_filename.find_first_of(":");
        string str_msg =  "For SSocket-type target, destination string ";
        str_msg += "must be '<hostname>:<port>'.\n";
        str_msg += "\tThe ':' char has to be present.";
        if ((unsigned)pos == (unsigned) string::npos) error(str_msg);
        str_remotehost = astr_filename.substr(0, pos);
        str_remoteport = astr_filename.substr(pos+1);
        port = atoi(str_remoteport.c_str());
        if (b_socketCreated) delete pcSS;
        pcSS = new c_SSocket_UDP_transmit(str_remotehost, port);
        break;
    }
    debug_pop();
    return true;
}

int
C_SMessage::printf(const char* format, ...) {
    char        pch_buffer[65536];
    int         ret             = 0;
    va_list     vp_arg;
    string      str_syslog      = "";
    string      str_buffer      = "";
    
    va_start(vp_arg, format);
    vsnprintf(pch_buffer, 65536, format, vp_arg);
    va_end(vp_arg);
    if(b_canPrint) {
        str_buffer      = pch_buffer;
        ret             = fprintf(pFILE_out, "%s", str_buffer.c_str());
    }
    fflush(pFILE_out);
    return ret;
}

int
C_SMessage::lprintf(const char* format, ...) {
    char        pch_buffer[65536];
    int         ret             = 0;
    va_list     vp_arg;
    string      str_syslog      = "";
    string      str_buffer      = "";
    
    va_start(vp_arg, format);
    vsnprintf(pch_buffer, 65536, format, vp_arg);
    va_end(vp_arg);
    if(b_canPrint) {
        if(b_syslogPrepend) {
            str_syslog      = syslog_prepend();
            str_buffer      = str_syslog + " " + pch_buffer;
        } else
            str_buffer      = pch_buffer;
        ret = fprintf(pFILE_out, "%-*s", lw, str_buffer.c_str());
    }
    fflush(pFILE_out);
    return ret;
}

int
C_SMessage::pprintf(const char* format, ...) {
    char        pch_buffer[65536];
    char        pch_bufferOut[131072];
    int         ret;
    va_list     vp_arg;
    string      str_syslog      = "";
    string      str_buffer      = "";

    va_start(vp_arg, format);
    vsnprintf(pch_buffer, 65536, format, vp_arg);
    va_end(vp_arg);
    if(b_syslogPrepend) {
        str_syslog      = syslog_prepend();
        str_buffer      = str_syslog + " " + pch_buffer;
    } else
        str_buffer      = pch_buffer;
    ret = sprintf(pch_bufferOut, "%s", str_buffer.c_str());
    str_payload         += pch_bufferOut;
    return ret;
}

int
C_SMessage::plprintf(const char* format, ...) {
    char        pch_buffer[65536];
    char        pch_bufferOut[131072];
    int         ret;
    va_list     vp_arg;
    string      str_syslog      = "";
    string      str_buffer      = "";

    va_start(vp_arg, format);
    vsnprintf(pch_buffer, 65536, format, vp_arg);
    va_end(vp_arg);
    if(b_syslogPrepend) {
        str_syslog      = syslog_prepend();
        str_buffer      = str_syslog + " " + pch_buffer;
    } else 
        str_buffer      = pch_buffer;
    ret = sprintf(pch_bufferOut, "%-*s", lw, str_buffer.c_str());
    str_payload         += pch_bufferOut;
    return ret;
}

int
C_SMessage::lprintf(
    string&     str_bufferOut,
    const char* format, ...)
{
    char        pch_buffer[65536];
    char        pch_bufferOut[131072];
    int         ret;
    va_list     vp_arg;
    string      str_syslog      = "";
    string      str_buffer      = "";

    va_start(vp_arg, format);
    vsnprintf(pch_buffer, 65536, format, vp_arg);
    va_end(vp_arg);
    if(b_syslogPrepend) {
        str_syslog      = syslog_prepend();
        str_buffer      = str_syslog + " " + pch_buffer;
    } else
        str_buffer      = pch_buffer;
    ret = sprintf(pch_bufferOut, "%-*s", lw, str_buffer.c_str());
    str_bufferOut  = pch_bufferOut;
    return ret;
}

int
C_SMessage::pcolprintf(
        const char*             pch_lstr,
        const char*             format, ...
        )
{
    char        pch_buffer[65536];
    char        pch_bufferOut[131072];
    int         retlw, retrw;
    int         len;
    va_list     vp_arg;
    string      str_syslog      = "";
    string      str_buffer      = "";

    va_start(vp_arg, format);
    vsnprintf(pch_buffer, 65536, format, vp_arg);
    va_end(vp_arg);
    if(b_syslogPrepend) {
        str_syslog      = syslog_prepend();
        str_buffer      = str_syslog + " " + pch_lstr;
    } else
        str_buffer      = pch_lstr;
    retlw = sprintf(pch_bufferOut, "%-*s", lw, str_buffer.c_str());
    str_payload.append(pch_bufferOut);
    len = strlen(pch_buffer);
    if(pch_buffer[len-1] == '\n') {
        pch_buffer[len-1] = '\0';
        retrw = sprintf(pch_bufferOut, "%-*s\n", rw, pch_buffer);
    } else
        retrw = sprintf(pch_bufferOut, "%-*s", rw, pch_buffer);
    str_payload.append(pch_bufferOut);
    return retlw + retrw;
}

int
C_SMessage::colprintf(
        string&                 str_bufferOut,
        const char*             pch_lstr,
        const char*             format, ...
        )
{
    char        pch_buffer[65536];
    char        pch_bufferOut[131072];
    int         retlw, retrw;
    int         len             = 0;
    va_list     vp_arg;
    string      str_syslog      = "";
    string      str_buffer      = "";

    va_start(vp_arg, format);
    vsnprintf(pch_buffer, 65536, format, vp_arg);
    va_end(vp_arg);
    if(b_syslogPrepend) {
        str_syslog      = syslog_prepend();
        str_buffer      = str_syslog + " " + pch_lstr;
    } else
        str_buffer      = pch_lstr;
    retlw = sprintf(pch_bufferOut, "%-*s", lw, str_buffer.c_str());
    str_bufferOut       = pch_bufferOut;
    len = strlen(pch_buffer);
    if(pch_buffer[len-1] == '\n') {
        pch_buffer[len-1] = '\0';
        retrw = sprintf(pch_bufferOut, "%-*s\n", rw, pch_buffer);
    } else
        retrw = sprintf(pch_bufferOut, "%-*s", rw, pch_buffer);
    str_bufferOut       += pch_bufferOut;
    return retlw + retrw;
}

int
C_SMessage::colprintf(
        const char*             pch_lstr,
        const char*             format, ...
        )
{
    char        pch_buffer[65536];
    int         retlw           = 0;
    int         retrw           = 0;
    int         len             = 0;
    va_list     vp_arg;
    string      str_syslog      = "";
    string      str_buffer      = "";

    va_start(vp_arg, format);
    vsnprintf(pch_buffer, 65536, format, vp_arg);
    va_end(vp_arg);
    if(b_canPrint) {
        if(b_syslogPrepend) {
            str_syslog      = syslog_prepend();
            str_buffer      = str_syslog + " " + pch_lstr;
        } else
            str_buffer      = pch_lstr;
        retlw   = fprintf(pFILE_out, "%-*s", lw, str_buffer.c_str());
        len     = strlen(pch_buffer);
        if(pch_buffer[len-1] == '\n') {
            pch_buffer[len-1] = '\0';
            retrw =  fprintf(pFILE_out, "%-*s\n", rw, pch_buffer);
        } else
            retrw =  fprintf(pFILE_out, "%-*s", rw, pch_buffer);
    }
    fflush(pFILE_out);
    return retlw + retrw;
}
