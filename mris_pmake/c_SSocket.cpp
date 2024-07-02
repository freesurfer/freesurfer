/**
 * @brief defines a simple wrapper class around standard Berkeley sockets.
 *
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
using namespace std;

#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/time.h>
#include <netdb.h>
#include <fcntl.h>

#if defined __sun__
#include <sys/file.h>
#endif

#include <c_SSocket.h>

#define FILE_BUF 32768

//
//\\\---
// Base class issues --->>>
/////---
//

//-----------------------------------------------------------
// Base socket constructors/destructor, output dump routines
//-----------------------------------------------------------

void c_SSocket::debug_push(string astr_currentProc) {
  if (stackDepth_get() >= SSTACKDEPTH-1) {
    cout << "Current stackDepth:\t" << stackDepth_get()             << endl;
    cout << "stackdepth limit:\t"   << SSTACKDEPTH                  << endl;
    for (int i=0; i<SSTACKDEPTH; i++)
      cout << "Stack depth: " << i << "\t" << str_proc_get(i)      << endl;
    error("Out of str_proc stack depth");
  }
  stackDepth_set(stackDepth_get()+1);
  str_proc_set(stackDepth_get(), astr_currentProc);
}

void c_SSocket::debug_pop() {
  stackDepth--;
}

void
c_SSocket::error(
  string  astr_msg  /* = ""  Error message */,
  int  code  /* = 1  Error ID  */) {

  // comment this line out - it is causing compiler errors for rocky8 gcc10 
  //extern int errno;

  cerr << "\nFatal error encountered.\n";
  cerr << "\tSSocket `" << str_name << "' (id: " << id << ")\n";
  cerr << "\tCurrent function: " << str_obj << "::" << str_proc_get() << "\n";
  cerr << "\t" << astr_msg << "\n";
  //cerr << "\t" << "System errno:\t" << errno << endl;
  perror("Error returned from system");
  print();
  cerr << "Throwing exception to (this) with code " << code << "\n\n";
  throw(this);
}

void
c_SSocket::warn(
  string astr_msg /* = ""  Warning message */,
  int  code  /* = 1  Warning code  */) {
  if (warnings) {
    cerr << "\nWarning.\n";
    cerr << "\tSSocket `" << str_name << "' (id: " << id << ")\n";
    cerr << "\tCurrent function: " << str_obj << "::" << str_proc_get() << "\n";
    cerr << "\t" << astr_msg << " (warning code: " << code << ")\n";
  }
}

void
c_SSocket::function_trace(
  string  astr_msg /* = ""  Trace message */,
  string  astr_separator /* = ""  Separator message */) {
  int i;
  string              str_tab                 = "";
  static string       str_objectName          = "";
  static string       str_funcName            = "";

  if (verbosity>=stackDepth) {
    cerr << astr_separator;
    for (i=0; i<stackDepth; i++)
      str_tab += "    ";
    if (str_objectName != str_name)
      cerr << "\nSSocket `" << str_name_get() << "' (id: " << id << ")" << endl;
    if (str_funcName != str_proc_get()) {
      cerr << "\n" << str_tab << "Current function: " << str_obj << "::";
      cerr << str_proc_get();
    }
    cerr << "\n" << str_tab << astr_msg << endl;
  }
  str_objectName      = str_name_get();
  str_funcName        = str_proc_get();
}

void
c_SSocket::core_construct(
  string          astr_name       /* = "unnamed" */,
  int             a_id            /* = -1  */,
  int             a_verbosity     /* = 0  */,
  int             a_warnings      /* = 0  */,
  int             a_stackDepth    /* = 0  */,
  string          astr_proc       /* = "noproc" */) {
  //
  // ARGS
  // astr_name  in              name of object
  //  a_id   in              internal id number of object
  //  a_verbosity  in  verbosity level
  // a_warnings  in  warnings level
  // a_stackDepth  in  stack depth (debugging)
  // astr_proc  in  currently executing proc (debugging)
  //
  // DESC
  //  Common core statements for all constructors
  //

  str_name    = astr_name;
  id  = a_id;
  verbosity = a_verbosity;
  warnings = a_warnings;
  stackDepth = a_stackDepth;

  str_proc_set(stackDepth_get(), "no name");
  str_obj = "C_SSocket";
}

c_SSocket::c_SSocket(
  int a_protocol,
  string astr_name  /* = "Unnamed SSocket" */,
  int  a_id   /* = 0   */) {
  //
  // Constructor
  //
  // ARGS
  // astr_name  in  Name of this SSocket
  // a_id   in  id of this SSocket
  //
  // DESC
  //  Base SSocket constructor.

  core_construct(astr_name, a_id);
  debug_push("c_SSocket");

  int sock;
  sock  = socket(AF_INET, a_protocol, 0);
  protocol = a_protocol;
  if (sock < 0)
    error("Problem creating socket");
  sockID = sock;
  debug_pop();
}


c_SSocket::c_SSocket(const c_SSocket &c_SSocket) {
  //
  // Copy constructor
  //

  debug_push("c_SSocket (copy constructor)");

  error("Copy constructor not yet defined");

  debug_pop();

}

c_SSocket & c_SSocket::operator=(const c_SSocket & c_SSocket) {
  //
  // Overloaded (=) operator
  //

  debug_push("operator=");

  error("Overloaded operator= not yet defined");

  return *this;

  debug_pop();
}

c_SSocket::~c_SSocket(void) {
  shutdown(sockID, 0);
  close(sockID);
}

void
c_SSocket::print() {
  cout << "Socket object name:\t"    << str_name  << endl;
  cout << "Socket object id:\t" << id   << endl;
  cout << "Socket object type:\t" << str_obj << endl;
  cout << "Socket protocol:\t" << ((protocol == SOCK_DGRAM) ? "Datagram" : "Stream") << endl;
  cout << "Socket ID:\t\t"  << sockID  << endl;
}

//
//\\\***
// c_SSocket_UDP definitions ****>>>>
/////***
//

c_SSocket_UDP::c_SSocket_UDP(
  string astr_name   /* = "Unnamed SSocket" */,
  int  aid    /* = 0   */) :

    c_SSocket(SOCK_DGRAM, astr_name, aid) {
  //
  // ARGS
  // astr_name  in  name of object variable
  // aid   in  id of object variable
  //
  // DESC
  // Basically a thin "fall-through" constructor to the base
  // SSocket class.
  //
  // HISTORY
  // 28 August 2001
  // o Initial design and coding.
  //

  debug_push("c_SSocket_UDP");
  readBufSize = 65536;
  debug_pop();
}

c_SSocket_UDP::~c_SSocket_UDP() {
  //
  // Destructor
  //
}

void
c_SSocket_UDP::print() {
  //
  // DESC
  // Simple info print method
  //

  c_SSocket::print();
}

int
c_SSocket_UDP::sendStr(
  string astr_text,
  bool ab_localecho /* = false */
) {
  //
  // ARGS
  //  astr_text in text to be sent
  //  ab_localecho in local echo of text
  //
  // DESC
  // "Transmits" the passed text over the socket.
  //
  // PRECONDITIONS
  // SSocket must be defined
  //
  // POSTCONDITIONS
  // The result of the socket system call `write' is returned.
  //

  debug_push("sendStr");

  int ret;

  ret = sendto(sockID, astr_text.c_str(), astr_text.length()+1, 0,
               (struct sockaddr*) &STsin_name, sizeof(STsin_name));

  if (ret < 0)
    error("Some error occurred while writing on stream socket");

  if (ab_localecho)
    cout << astr_text << endl;

  debug_pop();
  return ret;
}

int
c_SSocket_UDP::sendFile(
  string astr_filename,
  bool ab_localecho /* = false */
) {
  //
  // ARGS
  //  astr_filename  in name of file whose contents
  //      are to be sent over the socket
  // ab_localecho  in if true, echo text locally as well
  //
  // DESC
  // Reads in the contents of a filename and transmits
  // them to the remote host. The file's contents are
  // read into a single string and the whole string is
  // sent at once.
  //
  // PRECONDITIONS
  // SSocket must be defined.
  //
  // POSTCONDITIONS
  // return:
  //   1  success
  //   0  failure
  //
  //
  // TODO
  // o Dynamically allocating memory for file read buffer?
  // o Send packets of information instead of the entire
  //   file?
  //

  debug_push("sendStr");

  FILE*  pFILE_stream;
  char  pch_filebuf[1024];
  char ch;
  int  i = 0;
  int  ret;

  strcpy(pch_filebuf, "");
  if ( (pFILE_stream=fopen(astr_filename.c_str(), "r")) == NULL) {
    warn("Error opening file");
    return 0;
  }
  while (!feof(pFILE_stream)) {
    ret = fscanf(pFILE_stream, "%c", &ch);
    pch_filebuf[i++] = ch;
  }
  pch_filebuf[i] = '\0';
  ret = sendStr(pch_filebuf, 0);

  debug_pop();

  return ret;
}

int
c_SSocket_UDP::recv_timeout(
  string&   astr_payload,
  struct sockaddr*  apST_saFrom /* = NULL */,
  int*    ap_fromLen      /* = NULL */) {
  //
  // ARGS
  //  astr_payload out payload received
  //  apST_saFrom     out     If non-NULL, the sockaddr information from
  //                           the connecting host is returned with the
  //                           payload.
  //  ap_fromLen out     If pST_saFrom is non-NULL, this parameter
  //                           contains the length of the pST_saFrom
  //                           structure.
  //
  // DESC
  //  Listens on a socket (using timeout concerns).
  //
  //  If the last two parameters are non-NULL, the transmitting host's
  //  sockaddr address structure is returned.
  //
  // HISTORY
  // 23 November 2004
  // o Fixed minor bug regarding return value of rval for zero FD_ISSET().
  //
  // PRECONDITIONS
  // SSocket must be defined.
  //
  // POSTCONDITIONS
  // returns:
  //     val from read system call
  //

  debug_push("recv_timeout");

  int   rval = 0;
  static char* pch_buf = new char[readBufSize];

  astr_payload = "";

  FD_ZERO(&fd_ready);
  FD_SET(sockID, &fd_ready);
  STtval_timeout.tv_sec = timeoutSec;
  STtval_timeout.tv_usec = timeoutUsec;

  if (select(sockID+1, &fd_ready, NULL, NULL, &STtval_timeout) < 0)
    error("Problem with `select' system call");

  if (FD_ISSET(sockID, &fd_ready)) {
    rval = read(sockID, pch_buf, readBufSize);
    if (rval < 0)
      error("Problem reading Datagram packet");
    else if (rval>0) {
      astr_payload += pch_buf;
      bzero(pch_buf, sizeof(*pch_buf));
      debug_pop();
      return rval;
    }
  }

  debug_pop();
  return rval;
}

int
c_SSocket_UDP::recv_nonblocking(
  string&   astr_payload,
  struct sockaddr*  apST_saFrom /* = NULL */,
  int*    ap_fromLen      /* = NULL */) {
  //
  // ARGS
  //  astr_payload out payload received
  //  apST_saFrom     out     If non-NULL, the sockaddr information from
  //                           the connecting host is returned with the
  //                           payload.
  //  ap_fromLen out     If pST_saFrom is non-NULL, this parameter
  //                           contains the length of the pST_saFrom
  //                           structure.
  //
  // DESC
  //  Listens on a socket (non blocking).
  //
  //  If the last two parameters are non-NULL, the transmitting host's
  //  sockaddr address structure is returned.
  //
  // HISTORY
  // 17 November 2004
  // o pch_buf needs to be explicitly declared as a static array!
  //
  // PRECONDITIONS
  // SSocket must be defined and be configured to be non-blocking.
  //
  // POSTCONDITIONS
  // returns:
  //     val from read system call
  //

  debug_push("recv_nonblocking");

  static int   rval;
  static char* pch_buf = new char[readBufSize];

  astr_payload = "";

  rval = read(sockID, pch_buf, readBufSize);
  if (rval == 0) {
    cout << "Closing socket and returning\n";
    close(sockID);
    debug_pop();
    return 0;
  } else if (rval>0) {
    astr_payload += pch_buf;
    debug_pop();
    return rval;
  }

  debug_pop();
  return rval;
}

int
c_SSocket_UDP::recv(
  string&   astr_msg,
  struct sockaddr*  apST_saFrom /* = NULL */,
  int*   ap_fromLen      /* = NULL */
) {
  //
  // ARGS
  //  astr_payload out payload received
  //  apST_saFrom     out     If non-NULL, the sockaddr information from
  //                           the connecting host is returned with the
  //                           payload.
  //  ap_fromLen out     If pST_saFrom is non-NULL, this parameter
  //                           contains the length of the pST_saFrom
  //                           structure.
  //
  // DESC
  //  Dispatch method. Depending on value of timeout, calls either
  // recv_nonblocking or rec_timeout.
  //
  //  If the last two parameters are non-NULL, the transmitting host's
  //  sockaddr address structure is returned.
  //
  // PRECONDITIONS
  // SSocket must be defined and be configured to be non-blocking.
  //
  // POSTCONDITIONS
  // returns:
  //     val from read system call
  //
  // HISTORY
  // 21 November 2001
  //  o Initial design and coding.
  //

  debug_push("recv");
  int rval;

  if (timeoutSec || timeoutUsec)
    rval = recv_timeout( astr_msg,
                         apST_saFrom,
                         ap_fromLen);
  else
    rval = recv_nonblocking( astr_msg,
                             apST_saFrom,
                             ap_fromLen);

  debug_pop();
  return rval;
}


//
//\\\***
// c_SSocket_UDP_transmit definitions ****>>>>
/////***
//

c_SSocket_UDP_transmit::c_SSocket_UDP_transmit(
  string astr_remoteHostName,
  int a_port,
  string astr_name   /* = "Unnamed SSocket" */,
  int  aid    /* = 0   */
) : c_SSocket_UDP(astr_name, aid) {
  //
  // ARGS
  // astr_remoteHostName in  name of host to transmit to
  // a_port   in  port on remote host
  // astr_name  in  name of object variable
  // a_id   in  id of object variable
  //
  // DESC
  // Constructor of transmit UDP sockets
  //

  str_obj = "c_SSocket_UDP_transmit";
  debug_push(str_obj);

  port  = a_port;
  str_remoteHostName = astr_remoteHostName;

  SThost_hp       = new (struct hostent);
  SThost_hp           = gethostbyname(astr_remoteHostName.c_str());
  if (SThost_hp == NULL)
    error(astr_remoteHostName + ": unknown host");
  bzero(&STsin_name, sizeof(STsin_name));
  bcopy(SThost_hp->h_addr, &STsin_name.sin_addr, SThost_hp->h_length);
  STsin_name.sin_family = SThost_hp->h_addrtype;
  STsin_name.sin_port   = htons(port);

  debug_pop();
}

c_SSocket_UDP_transmit::~c_SSocket_UDP_transmit() {
  //
  // Destructor
  //
}

void
c_SSocket_UDP_transmit::print() {
  c_SSocket_UDP::print();
  cout << "Remote host:\t\t"  << str_remoteHostName << endl;
  cout << "Using port:\t\t"   << port   << endl;
}

//
//\\\***
// c_SSocket_UDP_receive definitions ****>>>>
/////***
//

c_SSocket_UDP_receive::c_SSocket_UDP_receive(
  int a_port   /* = 1701    */,
  int a_timeoutSec /* = 0     */,
  int a_timeoutUsec /* = 0     */,
  string astr_name  /* = "Unnamed SSocket_UDP_receive" */,
  int  aid   /* = 0     */) :

    c_SSocket_UDP(astr_name, aid) {
  //
  // ARGS
  //  a_port  in port to listen on
  // a_timeoutSec in timeout value (sec) on port
  // a_timeoutUsec in timeout value (usec) on port
  // astr_name in name of object variable
  // aid  in id of object variable
  //
  // DESC
  // `server'-type sockets can only be created on localhost
  //
  // HISTORY
  // 17 November 2004
  // o Added usec timeout
  //

  str_obj = "c_SSocket_UDP_receive";
  debug_push(str_obj);

  port = a_port;
  timeoutSec = a_timeoutSec;
  timeoutUsec = a_timeoutUsec;

  int length;

  // Name socket using wildcards
  STsin_name.sin_family         = AF_INET;
  STsin_name.sin_addr.s_addr    = INADDR_ANY;
  STsin_name.sin_port           = htons(port);

#if 0
  if (bind(sockID, (struct sockaddr *)&STsin_name, sizeof(STsin_name))) {
    char pch_error[1024];
    sprintf(pch_error, "Problem binding to socket %d", port);
    perror("System error");
    error(pch_error);
  }
#endif

  // Find assigned port
  length = sizeof(STsin_name);
#if defined __linux__ || defined __sun__ || defined Darwin
  if (getsockname(sockID,  (struct sockaddr *) &STsin_name,
                  (socklen_t *)  &length))
    error("Problem getting socket name");
#else
  if (getsockname(sockID,  (struct sockaddr *) &STsin_name,
                  (int *)   &length))
    error("Problem getting socket name");
#endif // __linux__ || __sun__

  port = ntohs(STsin_name.sin_port);

  // If timeout is zero, set socket to be non-blocking
  if (!(timeoutSec||timeoutUsec))
    if (fcntl(sockID, F_SETFL, O_NONBLOCK)<0)
      error("Problem with fcntl on socket");

  debug_pop();
}

c_SSocket_UDP_receive::~c_SSocket_UDP_receive() {
  //
  // Destructor
  //
}

void c_SSocket_UDP_receive::print() {
  c_SSocket_UDP::print();
  cout << "timeoutSec:\t\t"   << timeoutSec  << endl;
  cout << "timeoutUsec:\t\t"   << timeoutUsec  << endl;
  cout << "Using port:\t\t"   << port  << endl;
}

