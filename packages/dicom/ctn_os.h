/*
          Copyright (C) 1993, 1994, RSNA and Washington University

          The software and supporting documentation for the Radiological
          Society of North America (RSNA) 1993, 1994 Digital Imaging and
          Communications in Medicine (DICOM) Demonstration were developed
          at the
                  Electronic Radiology Laboratory
                  Mallinckrodt Institute of Radiology
                  Washington University School of Medicine
                  510 S. Kingshighway Blvd.
                  St. Louis, MO 63110
          as part of the 1993, 1994 DICOM Central Test Node project for, and
          under contract with, the Radiological Society of North America.

          THIS SOFTWARE IS MADE AVAILABLE, AS IS, AND NEITHER RSNA NOR
          WASHINGTON UNIVERSITY MAKE ANY WARRANTY ABOUT THE SOFTWARE, ITS
          PERFORMANCE, ITS MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR
          USE, FREEDOM FROM ANY COMPUTER DISEASES OR ITS CONFORMITY TO ANY
          SPECIFICATION. THE ENTIRE RISK AS TO QUALITY AND PERFORMANCE OF
          THE SOFTWARE IS WITH THE USER.

          Copyright of the software and supporting documentation is
          jointly owned by RSNA and Washington University, and free access
          is hereby granted as a license to use this software, copy this
          software and prepare derivative works based upon this software.
          However, any distribution of this software source code or
          supporting documentation or derivative works (source code and
          supporting documentation) must include the three paragraphs of
          the copyright notice.
*/
/* Copyright marker.  Copyright will be inserted above.  Do not remove */
/*
** @$=@$=@$=
*/
/*
**    DICOM 93
**       Electronic Radiology Laboratory
**     Mallinckrodt Institute of Radiology
**  Washington University School of Medicine
**
** Module Name(s):
** Author, Date: Stephen M. Moore, 15-Apr-93
** Intent:  This header defines public typedefs for the DICOM
**   software produced at the Mallinckrodt Institute of
**   Radiology.  It also defines unique identifiers
**   for standard classes and objects defined by the
**   standard.
** Last Update:  $Author: nicks $, $Date: 2006/12/29 02:09:01 $
** Source File:  $RCSfile: ctn_os.h,v $
** Revision:  $Revision: 1.3 $
** Status:  $State: Exp $
*/

#ifndef CTNOS_IS_IN
#define CTNOS_IS_IN 1

#ifdef  __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <signal.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#include <winsock.h>
#include <process.h>
#include <sys/timeb.h>
#include <direct.h>

  typedef SOCKET CTN_SOCKET;
#define CTN_BAD_SOCKET INVALID_SOCKET

#else
#include <unistd.h>
#include <sys/file.h>
#include <sys/socket.h>
  /*#include <sys/param.h>*/
#include <sys/time.h>
#include <sys/wait.h>
#include <netinet/in.h>
#include <netdb.h>
#include <sys/param.h>
#include <sys/utsname.h>
#include <dirent.h>

  typedef int CTN_SOCKET;
#define CTN_BAD_SOCKET -1
#endif

#ifdef SOLARIS
#include <fcntl.h>
#endif

#ifdef USEREGCOMP
#include <regex.h>
#endif


#ifdef  __cplusplus
}
#endif

#endif
