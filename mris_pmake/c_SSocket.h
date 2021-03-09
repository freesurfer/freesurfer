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

#ifndef __C_SOCK_H__
#define __C_SOCK_H__

#include <iostream>
#include <string>
using namespace std;

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <unistd.h>
#include <netinet/in.h>
#include <netdb.h>

const int       SSTACKDEPTH     = 64;

class c_SSocket {

  //
  // Data members
  //

public:

  // type and var info
  string      str_obj;                                // name of object class
  string str_name;                               // name of object variable
  int  id;                                 // id of socket
  int  verbosity;                         // Debug related value for object
  int  warnings;                         // Show warnings
  int  stackDepth;                         // Current procedure stackDepth
  string str_proc[SSTACKDEPTH];                  // Used to track the current
  //      procedure being
  //      executed
protected:

  // base class info
  int          protocol;  // datagram or stream
  int                 sockID;                 // the socket descriptor
  string  str_buf;                // payload buffer
  struct sockaddr_in  STsin_name;            // sockaddr_in structure
  // (shared by Tx/Rx)

  //
  // Method members
  //

public:

  //
  // Constructor / destructor block
  //
  void core_construct(
    string          astr_name       = "unnamed",
    int             a_id            = -1,
    int             a_verbosity     = 0,
    int             a_warnings      = 0,
    int             a_stackDepth    = 0,
    string          astr_proc       = "noproc"
  );
  c_SSocket(
    int a_protocol,
    string astr_name  = "Unnamed SSocket",
    int  aid   = 0
  );
  ~c_SSocket(void);
  c_SSocket(
    const c_SSocket & c_SSocket
  );
  c_SSocket & operator=(const c_SSocket & c_SSocket);

  //
  // Error / warn /  print block
  //
  void         debug_push(     string  astr_currentProc);
  void         debug_pop();
  void         error( string  astr_msg         = "",
                      int     code             = 1);
  void         warn(   string  astr_msg         = "",
                       int     code             = 1);
  void         function_trace( string  astr_msg = "",
                               string  astr_separator = "");
  void   print();

  //
  // Access "housekeeping" state info
  //

  const string        str_obj_get()           const {
    return str_obj;
  };
  const string        str_name_get()  const {
    return str_name;
  };
  int         id_get()  const {
    return id;
  };
  int         verbosity_get()  const {
    return verbosity;
  };
  int           warnings_get()          const {
    return warnings;
  };
  int           stackDepth_get() const {
    return stackDepth;
  };

  void         str_obj_set(string astr_val) {
    str_obj = astr_val;
  };
  void         str_name_set(string astr_val) {
    str_name = astr_val;
  };
  void  str_proc_set(int depth, string astr_proc) {
    str_proc[depth] = astr_proc;
  } ;
  void  id_set(int value) {
    id = value ;
  } ;
  void  verbosity_set(int value) {
    verbosity = value ;
  } ;
  void  warnings_set(int value) {
    warnings = value ;
  } ;
  void  stackDepth_set(int value) {
    stackDepth = value ;
  } ;
  const string        str_proc_get()
  const {
    return str_proc[stackDepth_get()];
  };
  const string        str_proc_get(int i) {
    return str_proc[i];
  };

  //
  // get/set member routines
  //
  int         sockID_get()          const {
    return sockID;
  };
  int  protocol_get()  const {
    return protocol;
  };

};

class c_SSocket_UDP : public c_SSocket {

protected:
  fd_set  fd_ready;  // file descriptor (for receive)
  struct timeval STtval_timeout;  // timeval structure (for receive)
  int                 timeoutSec;  // timeout (in sec) of socket on listen
  int                 timeoutUsec;  // timeout (in usec) of socket on listen
  int   port;   // port socket is associated with
  int   readBufSize;  // size of buffer

public:
  c_SSocket_UDP(
    string astr_name  = "Unnamed SSocket",
    int  aid   = 0
  );
  ~c_SSocket_UDP(void);

  //
  // Access block
  //
  int  port_get()
  const {
    return port;
  };
  int  timeoutSec_get()
  const {
    return timeoutSec;
  };
  int  timeoutUsec_get()
  const {
    return timeoutUsec;
  };
  int  readBufSize_get()
  const {
    return readBufSize;
  };
  void readBufSize_set( int aval) {
    readBufSize = aval;
  };

  //
  // Functional block
  //
  void print();
  int  sendStr(
    string astr_msg,
    bool ab_localecho  = false
  );

  int  sendFile(
    string astr_fileName,
    bool ab_localecho = false
  );
  int  recv_nonblocking(
    string&  astr_msg,
    struct sockaddr*  apST_saFrom = NULL,
    int*  ap_fromLen      = NULL
  );
  int  recv_timeout(
    string&  astr_msg,
    struct sockaddr*  apST_saFrom = NULL,
    int*  ap_fromLen      = NULL
  );
  int  recv(
    string&  astr_msg,
    struct sockaddr*  apST_saFrom = NULL,
    int*  ap_fromLen      = NULL
  );

};

class c_SSocket_UDP_receive : public c_SSocket_UDP {

protected:

public:
  //
  // Constructor / destructor block
  //
  c_SSocket_UDP_receive(
    int a_port,
    int a_timeoutSec = 0,
    int a_timeoutUsec = 0,
    string astr_name  = "Unnamed SSocket_UDP_receive",
    int  aid   = 0
  );
  ~c_SSocket_UDP_receive(void);

  //
  // Functional block
  //
  void print();

};

class c_SSocket_UDP_transmit : public c_SSocket_UDP {

protected:
  string  str_remoteHostName;     // remotehost name for "write" sockets
  struct hostent*     SThost_hp;  // hostent structure containing data
  // about remote host
public:
  //
  // Constructor / destructor block
  //
  c_SSocket_UDP_transmit(
    string astr_remoteHostName,
    int a_port,
    string astr_name  = "Unnamed SSocket",
    int  aid   = 0
  );
  ~c_SSocket_UDP_transmit(void);

  //
  // Functional block
  //
  void print();

};


#endif //  __C_SOCK_H__
