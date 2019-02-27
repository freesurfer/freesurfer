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
**    DICOM 93
**       Electronic Radiology Laboratory
**     Mallinckrodt Institute of Radiology
**  Washington University School of Medicine
**
** Module Name(s):
** Author, Date: Stephen M. Moore, 26-Apr-93
** Intent:
** This file defines private structures for the DICOM information
** object package.
** Last Update:  $Author: nicks $, $Date: 2006/12/29 02:09:01 $
** Source File:  $RCSfile: dcmprivate.h,v $
** Revision:  $Revision: 1.4 $
** Status:  $State: Exp $
*/

#ifdef  __cplusplus
extern "C"
{
#endif

  typedef struct
  {
    void *reserved[2];
    unsigned short group;
    /*    unsigned long groupLength; */
    unsigned long baseLength;
    int longVRAttributes;
    LST_HEAD *elementList;
  }
  PRV_GROUP_ITEM;

  typedef struct
  {
    void *reserved[2];
    DCM_ELEMENT element;
    int byteOrder;
    off_t dataOffset;
    off_t currentOffset;
    size_t allocatedDataLength;
    size_t originalDataLength;
    size_t paddedDataLength;
    int fragmentFlag;
  }
  PRV_ELEMENT_ITEM;

#define DCM_OBJUNDEFINED 0x01
#define DCM_OBJCOMMAND 0x02
#define DCM_OBJDATASET 0x03

  typedef struct
  {
    void *reserved[2];
    char keyType[32];
    DCM_OBJECTTYPE objectType;
    int accessMethod;
    CTNBOOLEAN deleteFlag;
    CTNBOOLEAN groupLengthFlag;
    unsigned long objectSize;
    unsigned long offset;
    unsigned long pixelSize;
    unsigned long pixelOffset;
    unsigned short pixelBitsAllocated;
    unsigned short pixelRepresentation;
    PRV_GROUP_ITEM *groupCtx;
    PRV_ELEMENT_ITEM *elementCtx;
    int fd;
    char fileName[1024];
    void *userCtx;
    CONDITION(*rd) (void *ctx, void *buf, int toRead, int *bytesRead);
    CONDITION(*sk) (void *ctx, int offset, int flag);
    LST_HEAD *groupList;
    CTNBOOLEAN preambleFlag;
    unsigned char preamble[DCM_PREAMBLELENGTH];
    unsigned long dataOptions;
    unsigned long metaHeaderLength;
    int longVRAttributes;
    char waveformDataVR[DICOM_CS_LENGTH+1];
  }
  PRIVATE_OBJECT;

#define KEY_DCM_OBJECT "KEY ACR NEMA V3 OBJECT"

#define DCM_FILE_ACCESS  1
#define DCM_MEMORY_ACCESS 2

  typedef union {
    unsigned short s;
    unsigned char u[2];
  }   SHORT_WORD;

  typedef union {
#ifdef __alpha
    unsigned int l;
#else
    unsigned long l;
#endif
    unsigned char u[4];
  }   LONG_WORD;

#define GET_SHORT_SAME_ORDER(A,B) {  \
 SHORT_WORD sss;    \
 sss.u[0] = (A)[0];   \
 sss.u[1] = (A)[1];   \
 (B) = sss.s;    \
}

#define GET_SHORT_REVERSE_ORDER(A,B) {  \
 SHORT_WORD sss;    \
 sss.u[0] = (A)[1];   \
 sss.u[1] = (A)[0];   \
 (B) = sss.s;    \
}

#define GET_LONG_SAME_ORDER(A,B) {  \
 LONG_WORD lll;    \
 lll.u[0] = (A)[0];   \
 lll.u[1] = (A)[1];   \
 lll.u[2] = (A)[2];   \
 lll.u[3] = (A)[3];   \
 (B) = lll.l;    \
}

#define GET_LONG_REVERSE_ORDER(A,B) {  \
 LONG_WORD lll;    \
 lll.u[0] = (A)[3];   \
 lll.u[1] = (A)[2];   \
 lll.u[2] = (A)[1];   \
 lll.u[3] = (A)[0];   \
 (B) = lll.l;    \
}

#ifdef  __cplusplus
}
#endif
