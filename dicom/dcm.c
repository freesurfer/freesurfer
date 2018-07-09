/* Modified Dec 15 2014 by DNG to handle double and float. The changes
were not extensive. Not sure how robust it is but it seems to work. */

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
**        DICOM 93
**         Electronic Radiology Laboratory
**       Mallinckrodt Institute of Radiology
**    Washington University School of Medicine
**
** Module Name(s):  DCM_OpenFile
**      DCM_CreateObject
**      DCM_CloseObject
**      DCM_AddElement
**      DCM_AddSequenceElement
**      DCM_RemoveElement
**      DCM_GetElementValue
**      DCM_GetElementSize
**      DCM_ScanParseObject
**      DCM_ImportStream
**      DCM_ExportStream
**      DCM_GetObjectSize
**      DCM_DumpElements
**      DCM_Debug
**      DCM_WriteFile
**      DCM_ModifyElements
**      DCM_ParseObject
**      DCM_RemoveGroup
**      DCM_GetSequenceList
**      DCM_ComputeExportSize
**  private functions
**      newElementItem
**      findCreateGroup
**      insertNewElement
**      updateObjectType
**      updateSpecialElements
**      exportFixedFields
**      exportData
**      fileSize
**      swapInPlace
**      checkObject
**      writeFile
**      countBytes
**      exportStream
**      verifyFormat
**      readFile
** Author, Date:  Stephen M. Moore, 26-Apr-93
** Intent:
**  This file contains the routines which implement a facility for
**  manipulating DICOM V3 information objects.  These routines parse
**  and construct NEMA objects (or streams).  Users are able to
**  read/write/examine individual attributes.  The package uses those
**  individual elements to construct an internal memory representation of
**  an object.  This representation is a linked list of the individual
**  attributes.  The user of the package can add attributes to an object,
**  remove an element from an object, query the object about an attribute,
**  and convert the object to and from its "stream" representation.
**  In addition, the package can parse a file which contains a stream
**  and create its internal object.
** Last Update:   $Author: greve $, $Date: 2014/12/15 22:39:45 $
** Source File:   $RCSfile: dcm.c,v $
** Revision:    $Revision: 1.33 $
** Status:    $State: Exp $
*/

#include <sys/fcntl.h>
#include <ctype.h>

#ifdef SunOS
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#define BYTEORDER_SAME      1
#define BYTEORDER_REVERSE   2
#define NATIVE_ORDER      BYTEORDER_SAME

#include "dicom_platform.h"
#ifdef BIG_ENDIAN_ARCHITECTURE
#define LITTLE_ORDER  BYTEORDER_REVERSE
#define BIG_ORDER BYTEORDER_SAME
#endif
#ifdef LITTLE_ENDIAN_ARCHITECTURE
#define LITTLE_ORDER  BYTEORDER_SAME
#define BIG_ORDER BYTEORDER_REVERSE
#endif
#ifndef LITTLE_ORDER
#error LITTLE_ORDER not defined!
#endif
#ifndef BIG_ORDER
#error BIG_ORDER not defined!
#endif

#include "ctn_os.h"
#include "dicom.h"
#include "condition.h"
#include "lst.h"
#include "dicom_uids.h"
#include "dicom_objects.h"
#include "dcmprivate.h"

static CTNBOOLEAN debug = FALSE;/* Flag for debugging messages to stdout */

// undeclared ones
#ifndef Darwin
#ifndef SunOS
#ifndef Windows_NT
extern void swab(const void *from, void *to, size_t n);
#endif
#endif
#endif

/* Prototypes for internal functions
 */
static CONDITION
newElementItem(DCM_ELEMENT * src, CTNBOOLEAN allocateData,
               PRV_ELEMENT_ITEM ** dst);
static CONDITION
findCreateGroup(PRIVATE_OBJECT ** object, unsigned short group,
                PRV_GROUP_ITEM ** groupPtr);
static CONDITION
insertNewElement(PRIVATE_OBJECT ** object,
                 DCM_ELEMENT * element);
static CONDITION
updateObjectType(PRIVATE_OBJECT ** object,
                 DCM_ELEMENT * element);
static CONDITION
updateSpecialElements(PRIVATE_OBJECT ** object,
                      PRV_ELEMENT_ITEM * item);
static void
exportFixedFields(DCM_ELEMENT * element,
                  unsigned char *b, U32 length, int byteOrder,
                  CTNBOOLEAN explicitVR, U32 * rtnLength);
static CONDITION
exportData(PRIVATE_OBJECT ** object, PRV_ELEMENT_ITEM * item,
           unsigned char *src,
           unsigned char *dst, U32 length, int byteOrder,
           U32 * rtnLength);
#ifdef MACOS
static long fileSize(int fd);
#else
static int fileSize(int fd);
#endif
static void swapInPlace(PRIVATE_OBJECT ** object, DCM_ELEMENT * e);
static CONDITION checkObject(PRIVATE_OBJECT ** object, const char *caller);
static CONDITION
writeFile(void *buffer, U32 length, int flag, void /* int */ *fd);
static CONDITION
countBytes(void *buffer, U32 length, int flag,
           void /* unsigned long */ *sizePtr);
static CONDITION
exportStream(DCM_OBJECT ** callerObject, unsigned long opt,
             void *buffer, U32 bufferlength, 
	     DCM_EXPORT_STREAM_CALLBACK* callback,
             void *ctx, int sequenceLevel);

static CONDITION
verifyFormat(PRV_ELEMENT_ITEM * item);
#if 0
static CONDITION
readFile(char *name, unsigned char *callerBuf, int fd, long size,
         off_t fileOffset, int recursionLevel,
         unsigned long opt, DCM_OBJECT ** callerObject,
         U32 * scannedLength, CTNBOOLEAN * remainOpenFlag,
         void *ctx,
         CONDITION(*rd) (void *ctx, void *buf, int toRead, int *bytesRead),
         CONDITION(*sk) (void *ctx, int offset, int flag));
#endif
static CONDITION
readFile1(const char *name, unsigned char *callerBuf, int fd, U32 size,
          off_t * fileOffset, int recursionLevel,
          unsigned long opt, PRIVATE_OBJECT ** parentObject,
          DCM_OBJECT ** callerObject,
          U32 * scannedLength, CTNBOOLEAN * remainOpenFlag,
          void *ctx,
          CONDITION(*rd) (void *ctx, void *buf, int toRead, int *bytesRead),
          CONDITION(*sk) (void *ctx, int offset, int flag));

static PRV_ELEMENT_ITEM *locateElement(PRIVATE_OBJECT ** obj, DCM_TAG tag);
static void computeVM(PRIVATE_OBJECT ** object, DCM_ELEMENT * element);
static void ctxSensitiveLookup(PRIVATE_OBJECT ** object,
                               DCM_ELEMENT * element);
static CONDITION
copyData(PRIVATE_OBJECT ** object, PRV_ELEMENT_ITEM * item,
         DCM_ELEMENT * to, U32 * rtnLength);
static CONDITION
readLengthToEnd(int fd, const char *fileName,
                unsigned long opt, U32 * lengthToEnd);
#ifdef LITTLE_ENDIAN_ARCHITECTURE
static void swapATGroupElement(DCM_ELEMENT * e);
#endif
static void
dumpBinaryData(void *d, DCM_VALUEREPRESENTATION vr,
               long vm, long vmLimit);
static void
compareGroup(PRV_GROUP_ITEM * g1, PRV_GROUP_ITEM * g2,
             void (*callback) (const DCM_ELEMENT * e1,
                               const DCM_ELEMENT * e2,
                               void *ctx),
             void *ctx);
static void remapFileName(const char *name, char *mapName);


/* DCM_OpenFile
**
** Purpose:
**  This function opens a file that conforms to the DICOM V3 standard.
**  The file is parsed and an internal representation of the data is created.
**  A handle is returned to the caller to allow access to the data.
**
** Parameter Dictionary:
**   name ASCIZ string giving path name to file to be opened.
**   opt  BITMASK giving options used when opening file and
**    interpreting data.  Current options have to do with the
**    byte order of the data in the file.  Legal masks for
**    this field are:
**      DCM_ORDERNATIVE
**      DCM_ORDERLITTLEENDIAN
**      DCM_ORDERBIGENDIAN
**   object Pointer to handle to return to caller  giving caller
**    future access to the data in the object.
**
** Return Values:
**
**  DCM_ELEMENTCREATEFAILED
**  DCM_ELEMENTLENGTHERROR
**  DCM_ELEMENTOUTOFORDER
**  DCM_FILEACCESSERROR
**  DCM_FILEOPENFAILED
**  DCM_ILLEGALOPTION
**  DCM_ILLEGALSTREAMLENGTH
**  DCM_LISTFAILURE
**  DCM_NORMAL
**  DCM_OBJECTCREATEFAILED
**  DCM_UNEVENELEMENTLENGTH
**
** Algorithm:
**  Determine file byte order (per caller's options)
**  Create new ACR object
**  Open file read only.
**  Determine size of file
**  While you have not reached end of file
**      Read next (group, element, length) triple
**      Lookup data element in dictionary
**      Read data value according to byte order and (group,element)
**      Check to see that data element is in numerical order
**      Add data element to linked list
**  End while
*/

CONDITION
DCM_OpenFile(const char *name, unsigned long opt, DCM_OBJECT ** callerObject) {
  CONDITION cond;
  int fd;
  off_t fileOffset = 0;
  U32 lengthToEnd;
  U32 size;
  CTNBOOLEAN
  remainFileOpen = FALSE; /* Leave file open after parse ? */

  if ((opt & (DCM_ORDERMASK | DCM_FILEFORMATMASK)) == 0)
    return COND_PushCondition(DCM_ILLEGALOPTION,
                              DCM_Message(DCM_ILLEGALOPTION), "Byte order",
                              "DCM_OpenFile");

#ifdef _MSC_VER
  fd = open(name, O_RDONLY | O_BINARY);
#else
  fd = open(name, O_RDONLY);
#endif
  if ((fd < 0) && ((opt & DCM_FILENAMEMASK) == DCM_TRYFILENAMECHANGE)) {
    char mapName[1024];
    perror(name);
    remapFileName(name, mapName);
#ifdef _MSC_VER
    fd = open(mapName, O_RDONLY | O_BINARY);
#else
    fd = open(mapName, O_RDONLY);
    if (fd < 0) {
      strcat(mapName, ".");
      fd = open(mapName, O_RDONLY);
    }
#endif
  }
  if (fd < 0) {
    return COND_PushCondition(DCM_FILEOPENFAILED,
                              DCM_Message(DCM_FILEOPENFAILED), name,
                              "DCM_OpenFile");
  }
  size = fileSize(fd);
  if (size <= 0)
    return DCM_FILEACCESSERROR;

  if ((opt & DCM_LENGTHTOENDMASK) == DCM_USELENGTHTOEND) {
    cond = readLengthToEnd(fd, name,
                           opt & (~DCM_LENGTHTOENDMASK), &lengthToEnd);
    if (cond != DCM_NORMAL) {
      (void) close(fd);
      return COND_PushCondition(DCM_FILEOPENFAILED,
                                DCM_Message(DCM_FILEOPENFAILED), name,
                                "DCM_OpenFile");
    }
    size = lengthToEnd;
    fileOffset = 24;
    (void) lseek(fd, 24, SEEK_SET);
  }
#ifdef OLDSMM
  cond = readFile(name, NULL, fd, size, 0, 0, opt, callerObject, NULL,
                  &remainFileOpen, NULL, NULL, NULL);
#endif
  cond = readFile1(name, NULL, fd, size, &fileOffset, 0, opt, NULL,
                   callerObject, NULL, &remainFileOpen, NULL, NULL, NULL);
  if ((cond != DCM_NORMAL) || !remainFileOpen)
    (void) close(fd);
  if (cond != DCM_NORMAL) {
    if (debug)
      DCM_DumpElements(callerObject, 1);
    return COND_PushCondition(DCM_FILEOPENFAILED,
                              DCM_Message(DCM_FILEOPENFAILED), name,
                              "DCM_OpenFile");
  } else
    return DCM_NORMAL;
}

CONDITION
DCM_ReadStream(DCM_OBJECT ** callerObject, unsigned long opt, long size,
               void *ctx,
               CONDITION(*rd) (void *ctx, void *buf, int toRead,
                               int *bytesRead),
               CONDITION(*sk) (void *ctx, int offset, int flag)) {
  CONDITION cond;
  int fd = -1;
  CTNBOOLEAN
  remainFileOpen = FALSE; /* Leave file open after parse ? */
  off_t fileOffset = 0;

  if ((opt & (DCM_ORDERMASK | DCM_FILEFORMATMASK)) == 0)
    return COND_PushCondition(DCM_ILLEGALOPTION,
                              DCM_Message(DCM_ILLEGALOPTION),
                              "Byte order",
                              "DCM_ReadStream");

  cond = readFile1("", NULL, fd, size, &fileOffset, 0, opt, NULL,
                   callerObject, NULL, &remainFileOpen, ctx, rd, sk);
  if (cond != DCM_NORMAL)
    return COND_PushCondition(DCM_READSTREAMFAILED,
                              DCM_Message(DCM_READSTREAMFAILED),
                              "DCM_ReadStream");
  else
    return DCM_NORMAL;
}

/* DCM_CreateObject
**
** Purpose:
**  This function creates a new object and initializes some
**  of the fields in the object
**
** Parameter Dictionary:
**  object    Pointer to caller's handle for object to be created.
**  opt   Flag with options used when creating object.
**      The only option that we use now is DCM_NOGROUPLENGTH.
**
** Return Values:
**  DCM_NORMAL
**  DCM_OBJECTCREATEFAILED
**  DCM_LISTFAILURE
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

CONDITION
DCM_CreateObject(DCM_OBJECT ** object, unsigned long opt) {
  PRIVATE_OBJECT
  * obj;

  if (object == NULL) {
    (void) COND_PushCondition(DCM_NULLADDRESS,
                              DCM_Message(DCM_NULLADDRESS),
                              "DCM_CreateObject");
    return COND_PushCondition(DCM_OBJECTCREATEFAILED,
                              DCM_Message(DCM_OBJECTCREATEFAILED),
                              "DCM_CreateObject");
  }
  obj = (PRIVATE_OBJECT *) CTN_MALLOC(sizeof(PRIVATE_OBJECT));
  if (obj == NULL) {
    (void) COND_PushCondition(DCM_MALLOCFAILURE,
                              DCM_Message(DCM_MALLOCFAILURE),
                              sizeof(PRIVATE_OBJECT),
                              "DCM_CreateObject");
    *object = NULL;
    return COND_PushCondition(DCM_OBJECTCREATEFAILED,
                              DCM_Message(DCM_OBJECTCREATEFAILED),
                              "DCM_CreateObject");
  }
  (void) memset(obj, 0, sizeof(PRIVATE_OBJECT));
  (void) strcpy(obj->keyType, KEY_DCM_OBJECT);


  obj->objectType = DCM_OBJECTUNKNOWN;
  obj->accessMethod = DCM_MEMORY_ACCESS;
  obj->deleteFlag = FALSE;
  if ((opt & DCM_GROUPLENGTHMASK) == DCM_NOGROUPLENGTH)
    obj->groupLengthFlag = FALSE;
  else
    obj->groupLengthFlag = TRUE;
  obj->objectSize = 0;
  obj->offset = 0;
  obj->pixelSize = 0;
  obj->pixelOffset = 0;
  obj->pixelBitsAllocated = 0;
  obj->pixelRepresentation = 0xffff;
  obj->groupCtx = NULL;
  obj->elementCtx = NULL;
  obj->fd = -1;
  obj->fileName[0] = '\0';
  obj->preambleFlag = FALSE;
  obj->preamble[0] = '\0';
  obj->dataOptions = 0;
  obj->metaHeaderLength = 0xffffffff;
  obj->longVRAttributes = 0;
  obj->waveformDataVR[0] = '\0';

  obj->groupList = LST_Create();
  if (obj->groupList == NULL) {
    CTN_FREE(obj);
    *object = NULL;
    return COND_PushCondition(DCM_LISTFAILURE,
                              DCM_Message(DCM_LISTFAILURE),
                              "DCM_CreateObject");
  }
  *object = (DCM_OBJECT *) obj;
  return DCM_NORMAL;
}


/* DCM_CloseObject
**
** Purpose:
**  Close an information object by freeing memory allocated to it and
**  destroying caller's reference to it.
**
** Parameter Dictionary:
**  callerObject  Address of caller's pointer to a DCM_OBJECT.
**
** Return Values:
**
**  DCM_FILEDELETEFAILED
**  DCM_ILLEGALOBJECT
**  DCM_LISTFAILURE
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Algorithm:
*/

CONDITION
DCM_CloseObject(DCM_OBJECT ** callerObject) {
  CONDITION
  cond;
  PRV_GROUP_ITEM
  * group;
  PRV_ELEMENT_ITEM
  * element;
  PRIVATE_OBJECT
  ** object;
  DCM_SEQUENCE_ITEM
  * sequenceItem;
  DCM_FRAGMENT_ITEM* fragmentItem;

  if (debug)
    fprintf(stderr, "Starting DCM_CloseObject\n");

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_CloseObject");
  if (cond != DCM_NORMAL)
    return cond;

  if ((*object)->fd > 0)
    (void) close((*object)->fd);

  if (debug)
    fprintf(stderr, "DCM_CloseObject: Legal object and file closed\n");

  while ((group = (PRV_GROUP_ITEM*)LST_Pop(&(*object)->groupList)) != NULL) {
    if (debug)
      fprintf(stderr, "DCM_CloseObject: group %04x\n", group->group);

    while ((element = (PRV_ELEMENT_ITEM*)LST_Pop(&group->elementList)) != NULL) {
      if (debug)
        fprintf(stderr, "DCM_CloseObject: Element %08x\n",
                (unsigned int) element->element.tag);
      if (element->element.representation == DCM_SQ) {
        if (debug)
          fprintf(stderr, "Sequence List Address: %lx\n",
                  (unsigned long) element->element.d.sq);
        if (element->element.d.sq != NULL) {
          while ((sequenceItem = (DCM_SEQUENCE_ITEM*)LST_Pop(&element->element.d.sq)) != NULL) {
            (void) DCM_CloseObject(&sequenceItem->object);
            CTN_FREE(sequenceItem);
          }
          (void) LST_Destroy(&element->element.d.sq);
        }
      } else if (element->fragmentFlag) {
        if (debug)
          fprintf(stderr, "Fragment List Address: %lx\n",
                  (unsigned long) element->element.d.fragments);
        if (element->element.d.fragments != NULL) {
          while ((fragmentItem = (DCM_FRAGMENT_ITEM*)LST_Pop(&element->element.d.fragments))
                 != NULL) {
            CTN_FREE(fragmentItem);
          }
          (void) LST_Destroy(&element->element.d.fragments);
        }
      }
      if (debug)
        fprintf(stderr, "DCM_CloseObject: free %8lx\n",
                (unsigned long) element);

      CTN_FREE(element);
    }
    cond = LST_Destroy(&group->elementList);
    if (cond != LST_NORMAL)
      return COND_PushCondition(DCM_LISTFAILURE,
                                DCM_Message(DCM_LISTFAILURE),
                                "DCM_CloseObject");
    CTN_FREE(group);
  }
  cond = LST_Destroy(&(*object)->groupList);
  if (cond != LST_NORMAL)
    return COND_PushCondition(DCM_LISTFAILURE,
                              DCM_Message(DCM_LISTFAILURE),
                              "DCM_CloseObject");

  cond = DCM_NORMAL;
  if ((*object)->deleteFlag) {
    if (unlink((*object)->fileName) != 0) {
      /**
         (void) COND_PushCondition(DCM_FILEDELETEFAILED, strerror(errno));
      **/
      cond = COND_PushCondition(DCM_FILEDELETEFAILED,
                                DCM_Message(DCM_FILEDELETEFAILED),
                                (*object)->fileName, strerror(errno),
                                "DCM_CloseObject");
    }
  }
  CTN_FREE(*object);
  *object = NULL;
  return cond;
}

/* DCM_AddElement
**
** Purpose:
**  Add an element to an existing DCM object
**
** Parameter Dictionary:
**  object    Pointer to caller's existing DCM object.
**  element   Pointer to DCM element to be added to object
**
** Return Values:
**
**  DCM_ILLEGALOBJECT
**  DCM_ILLEGALREPRESENTATION
**  DCM_INSERTFAILED
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Algorithm:
**  Check caller's object to make certain it is a legal object
**  Check element to see that caller is not trying to add
**      group length or length to end data elements.
**  Lookup element in the dictionary
**  If element is not in the dictionary, use caller's definition
**  If element is in the dictionary, make certain caller used
**      proper definitions or left things undefined.
**  Call findCreateGroup to
**      - create new group if this group does not exist
**      - create length to end element if necessary
**      - update object size for this object
**      - set CURRENT pointer in linked list to head of this group
**  Call insertNewElement to
**      - create a copy of the caller's element
**      - insert copy into linked list
**  Call updateObjectType to
**      - update this object as type COMMAND, DATASET, MESSAGE
*/

CONDITION
DCM_AddElement(DCM_OBJECT ** callerObject, DCM_ELEMENT * element) {
  CONDITION
  cond;
  DCM_ELEMENT
  localElement;
  PRIVATE_OBJECT
  ** object;
  PRV_GROUP_ITEM
  * groupItem;

  object = (PRIVATE_OBJECT **) callerObject;

  cond = checkObject(object, "DCM_AddElement");
  if (cond != DCM_NORMAL)
    return cond;

  if ((DCM_TAG_ELEMENT(element->tag) == 0x0000))
    return COND_PushCondition(DCM_ILLEGALADD,
                              DCM_Message(DCM_ILLEGALADD),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_AddElement");


  localElement = *element;

  cond = DCM_LookupElement(&localElement);
  if (cond != DCM_NORMAL) {
    (void) COND_PopCondition(0);
    localElement = *element;
  } else {
    if (localElement.representation == DCM_OT ||
        localElement.representation == DCM_CTX)
      localElement.representation = element->representation;
    if (element->representation != DCM_UN &&
        element->representation != localElement.representation) {
      return COND_PushCondition(DCM_ILLEGALREPRESENTATION,
                                DCM_Message(DCM_ILLEGALREPRESENTATION),
                                DCM_TAG_GROUP(element->tag),
                                DCM_TAG_ELEMENT(element->tag),
                                "DCM_AddElement");
    }
  }

  cond = findCreateGroup(object, DCM_TAG_GROUP(localElement.tag), &groupItem);
  if (cond != DCM_NORMAL)
    return COND_PushCondition(DCM_INSERTFAILED,
                              DCM_Message(DCM_INSERTFAILED),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_AddElement");

  cond = insertNewElement(object, &localElement);
  if (cond != DCM_NORMAL)
    return COND_PushCondition(DCM_INSERTFAILED,
                              DCM_Message(DCM_INSERTFAILED),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_AddElement");

  cond = updateObjectType(object, &localElement);
  if (cond != DCM_NORMAL)
    return COND_PushCondition(DCM_INSERTFAILED,
                              DCM_Message(DCM_INSERTFAILED),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_AddElement");

  return DCM_NORMAL;
}

/* DCM_AddSequenceElement
**
** Purpose:
**  Add a sequence element to an existing DCM object.  This
**  function takes ownership of the caller's sequence list
**  when it adds the element to the object.  The caller's
**  copy of the sequence list is removed.
**
** Parameter Dictionary:
**  object    Pointer to caller's existing DCM object.
**  element   Pointer to DCM element to be added to object
**
** Return Values:
**
**  DCM_ILLEGALOBJECT
**  DCM_ILLEGALREPRESENTATION
**  DCM_INSERTFAILED
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Algorithm:
*/

CONDITION
DCM_AddSequenceElement(DCM_OBJECT ** callerObject, DCM_ELEMENT * element) {
  CONDITION cond;
  DCM_ELEMENT localElement;
  PRIVATE_OBJECT **object;
  PRV_GROUP_ITEM *groupItem;

  object = (PRIVATE_OBJECT **) callerObject;

  cond = checkObject(object, "DCM_AddSequenceElement");
  if (cond != DCM_NORMAL)
    return cond;

  if ((DCM_TAG_ELEMENT(element->tag) == 0x0000))
    return COND_PushCondition(DCM_ILLEGALADD,
                              DCM_Message(DCM_ILLEGALADD),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_AddElement");


  localElement = *element;

  cond = DCM_LookupElement(&localElement);
  if (cond != DCM_NORMAL) {
    (void) COND_PopCondition(0);
    localElement = *element;
  } else {
    localElement.representation = element->representation;
  }
  if (localElement.representation != DCM_SQ) {
    return COND_PushCondition(DCM_NOTASEQUENCE,
                              DCM_Message(DCM_NOTASEQUENCE),
                              DCM_TAG_GROUP(localElement.tag),
                              DCM_TAG_ELEMENT(localElement.tag),
                              "DCM_AddSequenceElement");
  }
  cond = findCreateGroup(object, DCM_TAG_GROUP(localElement.tag), &groupItem);
  if (cond != DCM_NORMAL)
    return COND_PushCondition(DCM_INSERTFAILED,
                              DCM_Message(DCM_INSERTFAILED),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_AddSequenceElement");

  cond = insertNewElement(object, &localElement);
  if (cond != DCM_NORMAL)
    return COND_PushCondition(DCM_INSERTFAILED,
                              DCM_Message(DCM_INSERTFAILED),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_AddElement");

  cond = updateObjectType(object, &localElement);
  if (cond != DCM_NORMAL)
    return COND_PushCondition(DCM_INSERTFAILED,
                              DCM_Message(DCM_INSERTFAILED),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_AddSequenceElement");

  /*
   * We have taken ownership of the sequence list, so zero out caller's
   * copy
   */
  element->d.sq = NULL;

  return DCM_NORMAL;
}

/* DCM_RemoveElement
**
** Purpose:
**  This function removes a single element from an information object.
**
** Parameter Dictionary:
**  callerObject    Handle to the object
**  tag     The tag of the element to be removed
**
** Return Values:
**
**  DCM_ELEMENTNOTFOUND
**  DCM_ILLEGALOBJECT
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

CONDITION
DCM_RemoveElement(DCM_OBJECT ** callerObject, DCM_TAG tag) {
  PRIVATE_OBJECT
  ** object;
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem,
  *groupLengthItem;
  CONDITION
  cond;
  CTNBOOLEAN
  flag;
  unsigned short
  group,
  element;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_RemoveElement");
  if (cond != DCM_NORMAL)
    return cond;

  group = DCM_TAG_GROUP(tag);
  element = DCM_TAG_ELEMENT(tag);

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&((*object)->groupList));
  if (groupItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              group, element,
                              "DCM_RemoveElement");

  (void) LST_Position(&((*object)->groupList), groupItem);

  flag = FALSE;
  while ((groupItem != NULL) && (flag == FALSE)) {
    if (groupItem->group == group)
      flag = TRUE;
    else
      groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
  }
  if (flag == FALSE)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              group, element,
                              "DCM_RemoveElement");

  elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
  if (elementItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              group, element,
                              "DCM_RemoveElement");

  (void) LST_Position(&groupItem->elementList, elementItem);

  groupLengthItem = elementItem;
  if (DCM_TAG_ELEMENT(groupLengthItem->element.tag) != 0x0000)
    groupLengthItem = NULL;


  flag = FALSE;
  while ((elementItem != NULL) && (flag == FALSE)) {
    if (DCM_TAG_ELEMENT(elementItem->element.tag) == element)
      flag = TRUE;
    else
      elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
  }

  if (flag == FALSE)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              group, element,
                              "DCM_RemoveElement");

  if (groupItem->baseLength != DCM_UNSPECIFIEDLENGTH) {
    groupItem->baseLength -= elementItem->paddedDataLength + 2 + 2 + 4;
    if (groupLengthItem != NULL) {
      *groupLengthItem->element.d.ul = groupItem->baseLength;
    }
  }
  if ((*object)->objectSize != DCM_UNSPECIFIEDLENGTH)
    (*object)->objectSize -= elementItem->paddedDataLength + 2 + 2 + 4;
  if (elementItem->element.representation == DCM_OW ||
      elementItem->element.representation == DCM_OB ||
      elementItem->element.representation == DCM_SQ) {
    groupItem->longVRAttributes--;
    (*object)->longVRAttributes--;
  }
  (void) LST_Remove(&(groupItem->elementList), LST_K_AFTER);
  CTN_FREE(elementItem);
  return DCM_NORMAL;
}

/* DCM_GetElementValue
**
** Purpose:
**  This function retrieves the data from one data element and
**  returns it in a buffer allocated by the caller.  In the event
**  the data is larger than the caller's buffer, multiple calls
**  are used to retrieve the data.
**
** Parameter Dictionary:
**  object    Pointer to user's object containing desired element
**  element   DCM_ELEMENT structure containing (group,element)
**      specification of desired data element
**  rtnLength Pointer to caller variable to hold length of
**      data returned by this call.
**  ctx   Pointer to context variable used for multiple
**      calls to this function.  Caller should set the
**      pointer to NULL before the first call and not
**      touch the pointer again.
**
** Return Values:
**
**  DCM_CANNOTGETSEQUENCEVALUE
**  DCM_ELEMENTNOTFOUND
**  DCM_GETINCOMPLETE
**  DCM_ILLEGALCONTEXT
**  DCM_ILLEGALOBJECT
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Algorithm:
**  Check caller's object to make certain it is a legal DCM object.
**  Find head of the object linked list.
**  Search linked list sequentially until object is found or end
**  of list reached.
**  If end of list
**    return DCM_ELEMENTNOTFOUND
**  If CTX pointer containts NULL
**      Begin copy from beginning of element
**  else
**      Begin copy from address in CTX pointer
**  Copy data from element data area to user buffer
**  If copy is incomplete (remaining data longer than caller's buffer)
**      Update CTX pointer to point to next uncopied part of data
**      Return DCM_GETINCOMPLETE
**  else
**      Update CTX pointer to point past data area.
**      Return DCM_NORMAL
*/

CONDITION
DCM_GetElementValue(DCM_OBJECT ** callerObject, DCM_ELEMENT * element,
                    U32 * rtnLength, void **ctx) {
  PRIVATE_OBJECT
  ** object;
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem;
  ssize_t
  nBytes;
  CONDITION
  cond;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_GetElementValue");
  if (cond != DCM_NORMAL)
    return cond;

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&(*object)->groupList);
  if (groupItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_GetElementValue");

  (void) LST_Position(&(*object)->groupList, groupItem);
  while (groupItem != NULL) {
    if (groupItem->group == DCM_TAG_GROUP(element->tag))
      break;

    groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
  }
  if (groupItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_GetElementValue");

  elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
  if (elementItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_GROUP(element->tag),
                              "DCM_GetElementValue");

  (void) LST_Position(&groupItem->elementList, elementItem);
  while (elementItem != NULL) {
    if (elementItem->element.tag == element->tag) {
      unsigned char *p;
      U32 l;

      if (element->representation == DCM_SQ)
        return COND_PushCondition(DCM_CANNOTGETSEQUENCEVALUE,
                                  DCM_Message(DCM_CANNOTGETSEQUENCEVALUE),
                                  element->tag, "DCM_GetElementValue");

      p = (unsigned char*)*ctx;
      if ((long) p > elementItem->element.length)
        return COND_PushCondition(DCM_ILLEGALCONTEXT,
                                  DCM_Message(DCM_ILLEGALCONTEXT),
                                  "DCM_GetElementValue");

      l = MIN(element->length, (elementItem->element.length - (long) p));

      *rtnLength = l;
      {
        if (elementItem->element.d.ot == NULL) {
          if ((*object)->fd != -1) {
            (void) lseek((*object)->fd,
                         elementItem->dataOffset + (long) p, SEEK_SET);
            nBytes = read((*object)->fd, element->d.ot, (int) l);
          } else {
            (*object)->sk((*object)->userCtx,
                          (long) (elementItem->dataOffset + (long) p),
                          SEEK_SET);
	    int nBytes_int;
            cond = (*object)->rd((*object)->userCtx, element->d.ot, l,
                                 &nBytes_int);
	    nBytes = (ssize_t)nBytes_int;
          }
          if ((unsigned) nBytes != l) {
            return COND_PushCondition(DCM_FILEACCESSERROR,
                                      DCM_Message(DCM_FILEACCESSERROR),
                                      (*object)->fileName,
                                      "DCM_GetElementValue");
          }
#ifdef LITTLE_ENDIAN_ARCHITECTURE
          if (elementItem->element.representation == DCM_AT) {
            DCM_ELEMENT e;
            e = elementItem->element;
            e.length = l;
            e.d.ot = element->d.ot;
            swapATGroupElement(&e);
          }
#endif
          if (elementItem->byteOrder == BYTEORDER_REVERSE) {
            DCM_ELEMENT e;
            e = elementItem->element;
            e.length = l;
            e.d.ot = element->d.ot;
            swapInPlace(object, &e);
          }
        } else {
          unsigned char *q;
          q = (unsigned char *) elementItem->element.d.ot +
              (long) p;
          (void) memcpy(element->d.ot, q, l);
          if (elementItem->byteOrder == BYTEORDER_REVERSE) {
            DCM_ELEMENT e;
            e = elementItem->element;
            e.length = l;
            e.d.ot = element->d.ot;
            swapInPlace(object, &e);
          }
        }
        p += l;
        *ctx = (void *) p;
        if ((long) p == elementItem->element.length)
          return DCM_NORMAL;
        else
          return DCM_GETINCOMPLETE;
      }

    }
    elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
  }
  return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                            DCM_Message(DCM_ELEMENTNOTFOUND),
                            DCM_TAG_GROUP(element->tag),
                            DCM_TAG_ELEMENT(element->tag),
                            "DCM_GetElementValue");
}

char*
DCM_GetString(DCM_OBJECT** callerObject, DCM_TAG tag) {
  DCM_ELEMENT e;
  CONDITION cond;
  char* s;
  char tmp[64] = "";
  char b[64] = "";

  e.tag = tag;
  cond = DCM_GetElement(callerObject, tag, &e);
  if (cond != DCM_NORMAL) {
    COND_PopCondition(TRUE);
    return 0;
  }

  if (DCM_IsString(e.representation)) {
    s = (char*)malloc(e.length + 1);
    e.d.string = s;
    cond = DCM_ParseObject(callerObject, &e, 1, 0, 0, 0);
    if (cond != DCM_NORMAL) {
      free(s);
      s = 0;
    }
    return s;
  }

  if (e.representation == DCM_SQ) {
    return 0;
  }

  if (e.length > sizeof(b))
    return 0;

  e.d.ot = b;
  cond = DCM_ParseObject(callerObject, &e, 1, 0, 0, 0);
  if (cond != DCM_NORMAL) {
    COND_PopCondition(TRUE);
    return 0;
  }

  switch (e.representation) {
  case DCM_AT:
    strcpy(tmp, "<Unimplemented>");
    break;
  case DCM_FD:
    sprintf(tmp, "%f", *e.d.fd);
    break;
  case DCM_FL:
    sprintf(tmp, "%f", *e.d.fl);
    break;
  case DCM_SL:
    sprintf(tmp, "%d", *e.d.sl);
    break;
  case DCM_SQ:
    strcpy(tmp, "<Unimplemented>");
    break;
  case DCM_SS:
    sprintf(tmp, "%d", *e.d.ss);
    break;
  case DCM_UL:
    sprintf(tmp, "%u", *e.d.ul);
    break;
  case DCM_UN:
    strcpy(tmp, "<Unimplemented>");
    break;
  case DCM_US:
    sprintf(tmp, "%d", *e.d.us);
    break;
    /*case DCM_UNKNOWN:*/
  case DCM_RET:
  case DCM_CTX:
  case DCM_OB:
  case DCM_OW:
  case DCM_DLM:
  default:
    strcpy(tmp, "<Unimplemented>");
    break;
  }

  s = (char*) malloc(strlen(tmp) + 1);
  strcpy(s, tmp);

  return s;
}



CONDITION
DCM_GetElementValueOffset(DCM_OBJECT ** callerObject, DCM_ELEMENT * element,
                          unsigned long offset) {
  PRIVATE_OBJECT **object;
  PRV_ELEMENT_ITEM *elementItem;
  int nBytes;
  CONDITION cond;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_GetElementValue");
  if (cond != DCM_NORMAL)
    return cond;

  elementItem = locateElement(object, element->tag);
  if (elementItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(element->tag),
                              DCM_TAG_ELEMENT(element->tag),
                              "DCM_GetElementValueOffset");


  {
    unsigned char *p;
    U32 l;

    if (element->representation == DCM_SQ)
      return COND_PushCondition(DCM_CANNOTGETSEQUENCEVALUE,
                                DCM_Message(DCM_CANNOTGETSEQUENCEVALUE),
                                element->tag, "DCM_GetElementValueOffset");

    p = (unsigned char *) offset;
    ;
    if ((long) p > elementItem->element.length)
      return COND_PushCondition(DCM_BADOFFSET,
                                DCM_Message(DCM_BADOFFSET),
                                (int) offset,
                                (int) elementItem->element.length,
                                "DCM_GetElementValueLength");

    l = element->length;
    if (l + offset > elementItem->element.length) {
      return COND_PushCondition(DCM_BADLENGTH,
                                DCM_Message(DCM_BADLENGTH),
                                (int) offset, (int) l,
                                (int) elementItem->element.length,
                                "DCM_GetElementValueLength");
    } {
      if (elementItem->element.d.ot == NULL) {
        if ((*object)->fd != -1) {
          (void) lseek((*object)->fd,
                       elementItem->dataOffset + (long) p, SEEK_SET);
          nBytes = read((*object)->fd, element->d.ot, (int) l);
        } else {
          (*object)->sk((*object)->userCtx,
                        (long) (elementItem->dataOffset + (long) p),
                        SEEK_SET);
          cond = (*object)->rd((*object)->userCtx, element->d.ot, l,
                               &nBytes);
        }
        if ((unsigned) nBytes != l) {
          return COND_PushCondition(DCM_FILEACCESSERROR,
                                    DCM_Message(DCM_FILEACCESSERROR),
                                    (*object)->fileName,
                                    "DCM_GetElementValueValue");
        }
#ifdef LITTLE_ENDIAN_ARCHITECTURE
        if (elementItem->element.representation == DCM_AT) {
          DCM_ELEMENT e;
          e = elementItem->element;
          e.length = l;
          e.d.ot = element->d.ot;
          swapATGroupElement(&e);
        }
#endif
        if (elementItem->byteOrder == BYTEORDER_REVERSE) {
          DCM_ELEMENT e;
          e = elementItem->element;
          e.length = l;
          e.d.ot = element->d.ot;
          swapInPlace(object, &e);
        }
      } else {
        unsigned char *q;
        q = (unsigned char *) elementItem->element.d.ot +
            (long) p;
        (void) memcpy(element->d.ot, q, l);
        if (elementItem->byteOrder == BYTEORDER_REVERSE) {
          DCM_ELEMENT e;
          e = elementItem->element;
          e.length = l;
          e.d.ot = element->d.ot;
          swapInPlace(object, &e);
        }
      }
      return DCM_NORMAL;
    }

  }
}



/* DCM_GetElementSize
**
** Purpose:
**  Return the size of one data element in an ACR object.
**
** Parameter Dictionary:
**  object    Pointer to caller's ACR object
**  element   Pointer to ACR element that defines data element
**      of interest by specifying (group,element) pair
**  rtnLength Pointer to caller variable to hold returned
**      length of data element
**
** Return Values:
**
**  DCM_NORMAL
**  DCM_NULLOBJECT
**  DCM_ILLEGALOBJECT
**  DCM_ELEMENTNOTFOUND
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

CONDITION
DCM_GetElementSize(DCM_OBJECT ** callerObject, DCM_TAG tag,
                   U32 * rtnLength) {
  PRIVATE_OBJECT
  ** object;
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem;
  CONDITION
  cond;
  CTNBOOLEAN
  flag;
  unsigned short
  group,
  element;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_GetElementSize");
  if (cond != DCM_NORMAL)
    return cond;

  group = DCM_TAG_GROUP(tag);
  element = DCM_TAG_ELEMENT(tag);

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&((*object)->groupList));
  if (groupItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              group, element,
                              "DCM_GetElementSize");

  (void) LST_Position(&((*object)->groupList), groupItem);

  flag = FALSE;
  while ((groupItem != NULL) && (flag == FALSE)) {
    if (groupItem->group == group)
      flag = TRUE;
    else
      groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
  }
  if (flag == FALSE)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              group, element,
                              "DCM_GetElementSize");

  elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
  if (elementItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              group, element,
                              "DCM_GetElementSize");

  (void) LST_Position(&groupItem->elementList, elementItem);

  flag = FALSE;
  while ((elementItem != NULL) && (flag == FALSE)) {
    if (elementItem->element.tag == tag)
      flag = TRUE;
    else
      elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
  }

  if (flag == FALSE)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              group, element,
                              "DCM_GetElementSize");


  *rtnLength = elementItem->element.length;
  return DCM_NORMAL;
}


/* DCM_ScanParseObject
**
** Purpose:
**  DCM_ScanParseObject is used to allow a caller to examine every
**  element in a DICOM object and to parse the elements in the object.
**  The caller passes a list of elements to be parsed from the object.
**  This function examines each element in the object in order
**  (ascending group/element).  If the element in the object is found
**  in the caller's parse list, the element is parsed (and a value
**  placed in storage allocated by the caller).  If the element is
**  not found in the caller's list, a callback function is invoked
**  to notify the caller of the element.  When the callback function
**  is invoked, the arguments are:
**    DCM_ELEMENT *e  Pointer to the individual element
**    void *ctx Caller's context information
**
**  This function is very useful for determining exactly which
**  elements are present in an object without having to ask for
**  each one individually.
**
** Parameter Dictionary:
**  callerObject  Pointer to caller's DICOM object
**  buf   Unused
**  bufferSizd  Unused
**  vector    A vector of elements which are to be parsed.  An entry
**      in the vector gives the tag and describes where the
**      parsed data is to be stored.
**  vectorLength  Number of entries in the vector.
**  callback  Caller function invoked for an element that is in
**      the object but is not found in caller's list.
**  ctx   Context information that is passed to callback function.
**
** Return Values:
**
**  DCM_NORMAL
**  DCM_NULLOBJECT
**  DCM_ILLEGALOBJECT
**
** Algorithm:
*/

CONDITION
DCM_ScanParseObject(DCM_OBJECT ** callerObject, void *buf, size_t bufferSize,
                    DCM_FLAGGED_ELEMENT * elementVector, int vectorLength,
                    CONDITION(*callback) (const DCM_ELEMENT* e, void* ctx),
                    void *ctx) {
  PRIVATE_OBJECT
  ** object;
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem;
  CONDITION
  cond;
  CTNBOOLEAN
  done = FALSE;
  DCM_ELEMENT
  e;
  int
  i;
  CTNBOOLEAN
  found;
  U32
  l=0;
  char
  *p;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_ScanParseObject");
  if (cond != DCM_NORMAL)
    return cond;

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&((*object)->groupList));
  (void) LST_Position(&((*object)->groupList), groupItem);
  while (groupItem != NULL && !done) {
    elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
    (void) LST_Position(&groupItem->elementList, elementItem);
    while (elementItem != NULL && !done) {
      for (found = FALSE, i = 0; !found && i < vectorLength; i++) {
        if (elementItem->element.tag == elementVector[i].e.tag) {
          found = TRUE;
          (void)copyData(object,elementItem,&elementVector[i].e, &l);
          *elementVector[i].flagAddress |= elementVector[i].flag;

          if (DCM_IsString(elementVector[i].e.representation)) {
            elementVector[i].e.d.string[l] = '\0';
            p = elementVector[i].e.d.string + l - 1;
            while (p >= elementVector[i].e.d.string && (*p == ' '))
              *p-- = '\0';
          }
        }
      }
      if (!found) {
        e = elementItem->element;
        cond = callback(&e, ctx);
        if (cond != DCM_NORMAL)
          done = TRUE;
      }
      elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
    }
    groupItem = (PRV_GROUP_ITEM*)LST_Next(&((*object)->groupList));
  }
  return DCM_NORMAL;
}

/* DCM_ImportStream
**
** Purpose:
**  Import data from memory in DCM stream format and create
**  an internal memory representation of the object.
**
** Parameter Dictionary:
**  buf   Pointer to caller's buffer containing ACR NEMA data
**  length    Length of input data in bytes
**  opt   Bitmask giving options for interpreting data.
**      Legal values specify the order of the bytes in the data
**        ACR_ORDERNATIVE
**        ACR_ORDERLITTLEENDIAN
**        ACR_ORDERBIGENDIAN
**  object    Pointer to object created and returned by this function
**
** Return Values:
**
**  DCM_ELEMENTCREATEFAILED
**  DCM_ELEMENTLENGTHERROR
**  DCM_ELEMENTOUTOFORDER
**  DCM_FILEACCESSERROR
**  DCM_ILLEGALOPTION
**  DCM_ILLEGALSTREAMLENGTH
**  DCM_LISTFAILURE
**  DCM_NORMAL
**  DCM_OBJECTCREATEFAILED
**  DCM_UNEVENELEMENTLENGTH
**
** Algorithm:
**  call private import stream function which handles recursion
**
*/

CONDITION
DCM_ImportStream(unsigned char *buf, unsigned long length,
                 unsigned long opt, DCM_OBJECT ** callerObject) {
#ifdef DEBUG
  if (debug)
    (void) fprintf(stderr, "DCM_ImportStream, %ld bytes\n", length);
#endif

  if ((opt & DCM_ORDERMASK) == 0)
    return COND_PushCondition(DCM_ILLEGALOPTION,
                              DCM_Message(DCM_ILLEGALOPTION), "Byte order",
                              "DCM_ImportStream");

  return readFile1("", buf, -1, length, 0, 0, opt, NULL,
                   callerObject, NULL, NULL,
                   NULL, NULL, NULL);
}

/* DCM_ExportStream
**
** Purpose:
**  Export a DICOM object into the stream format suitable
**  for network transmission or disk storage.
**
** Parameter Dictionary:
**  object    Pointer to caller's DICOM object
**  opt   Bitmask giving options for exporting data.  Legal
**      options give the byte order of exported data:
**        DCM_ORDERNATIVE
**        DCM_ORDERLITTLEENDIAN
**        DCM_ORDERBIGENDIAN
**  buffer    Pointer to caller's buffer to hold next slug
**      of DCM stream data.
**  bufferlength  Length of caller's buffer to hold stream data.
**  returnlength  Pointer to caller's variable into which we write
**      the amount of data exported.
**  ctx   Pointer to context variable we maintain to keep
**      track of our location in export process.
**
** Return Values:
**
**  DCM_FILEACCESSERROR
**  DCM_ILLEGALOBJECT
**  DCM_LISTFAILURE
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

CONDITION
DCM_ExportStream(DCM_OBJECT ** callerObject, unsigned long opt,
                 void *buffer, unsigned long bufferlength,
                 DCM_EXPORT_STREAM_CALLBACK* callback,
                 void *ctx) {


  PRIVATE_OBJECT
  ** object;
  CONDITION
  cond;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_ExportStream");
  if (cond != DCM_NORMAL)
    return cond;

  return exportStream(callerObject, opt, buffer, bufferlength, callback,
                      ctx, 0);
}

/* DCM_GetObjectSize
**
** Purpose:
**  Return the size of a DICOM object when it is represented in
**  stream format.
**
** Parameter Dictionary:
**  object    Pointer to caller's DICOM object
**  returnlength  Pointer to unsigned long variable to hold length of
**      object
**
** Return Values:
**
**  DCM_ILLEGALOBJECT
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

CONDITION
DCM_GetObjectSize(DCM_OBJECT ** callerObject, unsigned long *returnlength) {
  PRIVATE_OBJECT
  ** object;
  CONDITION
  cond;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_GetObjectSize");
  if (cond != DCM_NORMAL)
    return cond;

  *returnlength = (*object)->objectSize;
  return DCM_NORMAL;
}

/* DCM_DumpElements
**
** Purpose:
**  Dump a short description of each data element in an object to
**  stdout (for use as a debugging tool).
**
** Parameter Dictionary:
**  object    Pointer to caller's handle for DCM object to be dumped
**  vm    Limit on the value multiplicity for printing
**      binary data.  Print no more than vm values.
**
** Return Values:
**
**  DCM_ILLEGALOBJECT
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Algorithm:
**  Check caller's handle to make certain it is a legal DCM object
**  Print object type (COMMAND, DATASET, MESSAGE) and size in bytes
**  For each GROUP in the object linked list
**      For each ELEMENT ITEM in the group linked list
**    print group, element, size, description
**    print some or all of data based on data element representation
**      (ASCII number, ASCII text, binary)
**      End for
**  End for
*/
CONDITION
DCM_DumpElements(DCM_OBJECT ** callerObject, long vm) {
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem;
  PRIVATE_OBJECT
  ** object;
  CONDITION
  cond;
  DCM_SEQUENCE_ITEM
  * sq;
  char
  scratch[128];
  int
  stringLength;

  object = (PRIVATE_OBJECT **) callerObject;

  cond = checkObject(object, "DCM_DumpElements");
  if (cond != DCM_NORMAL)
    return cond;

  printf("\nDCM Dump Elements\n");
  switch ((*object)->objectType) {
  case DCM_OBJECTUNKNOWN:
    printf("Object type: UNKNOWN\n");
    break;
  case DCM_OBJECTCOMMAND:
    printf("Object type: COMMAND\n");
    break;
  case DCM_OBJECTIMAGE:
    printf("Object type: IMAGE\n");
    break;
  case DCM_OBJECTELEMENTLIST:
    printf("Object type: ELEMENT LIST\n");
    break;
  default:
    printf("Object type: Unknown (error)\n");
    break;
  }
  printf("Object size: %ld\n", (*object)->objectSize);

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&(*object)->groupList);
  if (groupItem != NULL)
    (void) LST_Position(&(*object)->groupList, groupItem);

  while (groupItem != NULL) {
    printf("Group: %04x, Length: %8ld\n", groupItem->group,
           groupItem->baseLength);
    elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
    if (elementItem != NULL)
      (void) LST_Position(&groupItem->elementList, elementItem);
    while (elementItem != NULL) {
      (void) printf("%04x %04x %8d ",
                    DCM_TAG_GROUP(elementItem->element.tag),
                    DCM_TAG_ELEMENT(elementItem->element.tag),
                    elementItem->element.length);
      (void) printf("//%31s//", elementItem->element.description);
      if (elementItem->element.d.ot == NULL)
        (void) printf("Data on disk\n");
      else {
        switch (elementItem->element.representation) {
        case DCM_AE:
        case DCM_AS:
        case DCM_CS:
        case DCM_DA:
        case DCM_DT:
          stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
          strncpy(scratch, elementItem->element.d.string, stringLength);
          scratch[stringLength] = '\0';
          (void) printf("%s\n", scratch);
          break;
        case DCM_DD:
          (void) printf("Unimplemented\n");
          break;
        case DCM_FD:
          (void) printf("%lf\n",*elementItem->element.d.fd);
          break;
        case DCM_FL:
          (void) printf("%lf\n",*elementItem->element.d.fd);
          break;
        case DCM_DS:
        case DCM_IS:
        case DCM_LO:
        case DCM_LT:
        case DCM_PN:
        case DCM_SH:
        case DCM_UT:
          stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
          strncpy(scratch, elementItem->element.d.string, stringLength);
          scratch[stringLength] = '\0';
          (void) printf("%s\n", scratch);
          break;
        case DCM_SL:
          (void) printf("%8x %d\n", *elementItem->element.d.sl,
                        *elementItem->element.d.sl);
          if (vm > 1)
            dumpBinaryData(elementItem->element.d.ot,
                           elementItem->element.representation,
                           elementItem->element.length / sizeof(U32), vm);
          break;
        case DCM_SS:
          (void) printf("%4x %d\n", *elementItem->element.d.ss,
                        *elementItem->element.d.ss);
          if (vm > 1)
            dumpBinaryData(elementItem->element.d.ot,
                           elementItem->element.representation,
                           elementItem->element.length / sizeof(short), vm);
          break;
        case DCM_SQ:
          (void) printf("SEQUENCE\n");
          sq = (DCM_SEQUENCE_ITEM*)LST_Head(&elementItem->element.d.sq);
          if (sq != NULL)
            (void) LST_Position(&elementItem->element.d.sq, sq);
          printf("DCM Dump Sequence\n");
          while (sq != NULL) {
            (void) DCM_DumpElements(&sq->object, vm);
            sq = (DCM_SEQUENCE_ITEM*)LST_Next(&elementItem->element.d.sq);
          }
          printf("DCM Dump Sequence Complete\n");
          break;
        case DCM_ST:
          stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
          strncpy(scratch, elementItem->element.d.string, stringLength);
          scratch[stringLength] = '\0';
          (void) printf("%s\n", scratch);
          break;
        case DCM_TM:
        case DCM_UI:
          stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
          strncpy(scratch, elementItem->element.d.string, stringLength);
          scratch[stringLength] = '\0';
          (void) printf("%s\n", scratch);
          break;
        case DCM_AT:
        case DCM_UL:
          (void) printf("%8x %d\n", *elementItem->element.d.ul,
                        *elementItem->element.d.ul);
          if (vm > 1)
            dumpBinaryData(elementItem->element.d.ot,
                           elementItem->element.representation,
                           elementItem->element.length / sizeof(U32), vm);
          break;
        case DCM_US:
          (void) printf("%4x %d\n", *elementItem->element.d.us,
                        *elementItem->element.d.us);
          if (vm > 1)
            dumpBinaryData(elementItem->element.d.ot,
                           elementItem->element.representation,
                           elementItem->element.length / \
                           sizeof(unsigned short), vm);
          break;
        case DCM_OB:
        case DCM_UN:
          dumpBinaryData(elementItem->element.d.ot,
                         elementItem->element.representation,
                         elementItem->element.length , 8);
          break;

        case DCM_OT:
        case DCM_OW:
          /*case DCM_UNKNOWN:*/
        case DCM_RET:
          (void) printf("Unimplemented\n");
          break;
        default:
          (void) printf("Some unimplemented logic if here\n");
          break;
        }
      }
      elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
    }
    groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
  }

  printf("DCM Dump Elements Complete\n\n");
  return DCM_NORMAL;
}

CONDITION
DCM_FormatElements(DCM_OBJECT ** callerObject, long vm, const char* prefix) {
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem;
  PRIVATE_OBJECT
  ** object;
  CONDITION
  cond;
  DCM_SEQUENCE_ITEM
  * sq;
  char
  scratch[128];
  int
  stringLength;
  char localPrefix[128];

  object = (PRIVATE_OBJECT **) callerObject;

  cond = checkObject(object, "DCM_DumpElements");
  if (cond != DCM_NORMAL)
    return cond;

  printf("\n%sDCM Dump Elements\n", prefix);
  switch ((*object)->objectType) {
  case DCM_OBJECTUNKNOWN:
    printf("%sObject type: UNKNOWN\n", prefix);
    break;
  case DCM_OBJECTCOMMAND:
    printf("%sObject type: COMMAND\n", prefix);
    break;
  case DCM_OBJECTIMAGE:
    printf("%sObject type: IMAGE\n", prefix);
    break;
  case DCM_OBJECTELEMENTLIST:
    printf("%sObject type: ELEMENT LIST\n", prefix);
    break;
  default:
    printf("%sObject type: Unknown (error)\n", prefix);
    break;
  }
  printf("%sObject size: %ld\n", prefix, (*object)->objectSize);

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&(*object)->groupList);
  if (groupItem != NULL)
    (void) LST_Position(&(*object)->groupList, groupItem);

  while (groupItem != NULL) {
    printf("%sGroup: %04x, Length: %8lu\n", prefix, groupItem->group,
           groupItem->baseLength);
    elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
    if (elementItem != NULL)
      (void) LST_Position(&groupItem->elementList, elementItem);
    while (elementItem != NULL) {
      (void) printf("%s%04x %04x %8d ",
                    prefix,
                    DCM_TAG_GROUP(elementItem->element.tag),
                    DCM_TAG_ELEMENT(elementItem->element.tag),
                    elementItem->element.length);
      (void) printf("//%31s//", elementItem->element.description);
      if (elementItem->element.d.ot == NULL)
        (void) printf("Data on disk\n");
      else {
        switch (elementItem->element.representation) {
        case DCM_AE:
        case DCM_AS:
        case DCM_CS:
        case DCM_DA:
        case DCM_DT:
          stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
          strncpy(scratch, elementItem->element.d.string, stringLength);
          scratch[stringLength] = '\0';
          (void) printf("%s\n", scratch);
          break;
        case DCM_DD:
        case DCM_FD:
        case DCM_FL:
          (void) printf("Unimplemented\n");
          break;
        case DCM_DS:
        case DCM_IS:
        case DCM_LO:
        case DCM_LT:
        case DCM_PN:
        case DCM_SH:
        case DCM_UT:
          stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
          strncpy(scratch, elementItem->element.d.string, stringLength);
          scratch[stringLength] = '\0';
          (void) printf("%s\n", scratch);
          break;
        case DCM_SL:
          (void) printf("%8x %d\n", (unsigned int) *elementItem->element.d.sl,
                        *elementItem->element.d.sl);
          if (vm > 1)
            dumpBinaryData(elementItem->element.d.ot,
                           elementItem->element.representation,
                           elementItem->element.length / sizeof(U32), vm);
          break;
        case DCM_SS:
          (void) printf("%4x %d\n", *elementItem->element.d.ss,
                        *elementItem->element.d.ss);
          if (vm > 1)
            dumpBinaryData(elementItem->element.d.ot,
                           elementItem->element.representation,
                           elementItem->element.length / sizeof(short), vm);
          break;
        case DCM_SQ:
          (void) printf("SEQUENCE\n");
          sq = (DCM_SEQUENCE_ITEM*)LST_Head(&elementItem->element.d.sq);
          if (sq != NULL)
            (void) LST_Position(&elementItem->element.d.sq, sq);
          printf("%sDCM Dump Sequence\n", prefix);
          strcpy(localPrefix, prefix);
          strcat(localPrefix, " ");
          while (sq != NULL) {
            (void) DCM_FormatElements(&sq->object, vm, localPrefix);
            sq = (DCM_SEQUENCE_ITEM*)LST_Next(&elementItem->element.d.sq);
          }
          printf("%sDCM Dump Sequence Complete\n", prefix);
          break;
        case DCM_ST:
          stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
          strncpy(scratch, elementItem->element.d.string, stringLength);
          scratch[stringLength] = '\0';
          (void) printf("%s\n", scratch);
          break;
        case DCM_TM:
        case DCM_UI:
          stringLength = MIN(sizeof(scratch) - 1, elementItem->element.length);
          strncpy(scratch, elementItem->element.d.string, stringLength);
          scratch[stringLength] = '\0';
          (void) printf("%s\n", scratch);
          break;
        case DCM_AT:
        case DCM_UL:
          (void) printf("%8x %u\n", (unsigned int) *elementItem->element.d.ul,
                        *elementItem->element.d.ul);
          if (vm > 1)
            dumpBinaryData(elementItem->element.d.ot,
                           elementItem->element.representation,
                           elementItem->element.length / sizeof(U32), vm);
          break;
        case DCM_US:
          (void) printf("%4x %d\n", *elementItem->element.d.us,
                        *elementItem->element.d.us);
          if (vm > 1)
            dumpBinaryData(elementItem->element.d.ot,
                           elementItem->element.representation,
                           elementItem->element.length / \
                           sizeof(unsigned short), vm);
          break;
        case DCM_OT:
        case DCM_OW:
        case DCM_OB:
          /*case DCM_UNKNOWN:*/
        case DCM_RET:
          (void) printf("Unimplemented\n");
          break;
        default:
          (void) printf("Some unimplemented logic if here\n");
          break;
        }
      }
      elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
    }
    groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
  }

  printf("%sDCM Dump Elements Complete\n\n", prefix);
  return DCM_NORMAL;
}

/* DCM_Debug
**
** Purpose:
**  To enable the debugging facility
**
** Parameter Dictionary:
**  flag  CTNBOOLEAN variable TRUE if caller wants to turn on debugging
**    info; FALSE otherwise
**
** Return Values:
**  None
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

void
DCM_Debug(CTNBOOLEAN flag) {
  debug = flag;
}

/* DCM_WriteFile
**
** Purpose:
**  Export an object from the internal representation and
**  write the stream representation to a file.
**
** Parameter Dictionary:
**  object    DCM_OBJECT which is to be written to the file
**  opt   Bitmask giving options for exporting data.  Legal
**      options give the byte order of exported data:
**        DCM_ORDERNATIVE
**        DCM_ORDERLITTLEENDIAN
**        DCM_ORDERBIGENDIAN
**  file    ASCIIZ name of the file to be created.
**
** Return Values:
**
**  DCM_FILEACCESSERROR
**  DCM_FILECREATEFAILED
**  DCM_ILLEGALOBJECT
**  DCM_LISTFAILURE
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

CONDITION
DCM_WriteFile(DCM_OBJECT ** callerObject, unsigned long opt, const char *file) {
  PRIVATE_OBJECT
  ** object;
  int
  fd;
  unsigned char
  buf[2048];
  CONDITION
  cond;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_WriteFile");
  if (cond != DCM_NORMAL)
    return cond;
#ifdef MACOS
  fd = open(file, O_WRONLY | O_CREAT | O_TRUNC);
#elif _MSC_VER
  fd = _open(file, O_WRONLY | O_CREAT | O_TRUNC | O_BINARY,
             _S_IREAD | _S_IWRITE);
#else
  fd = open(file, O_WRONLY | O_CREAT | O_TRUNC, 0666);
#endif
  if (fd < 0) {
    return COND_PushCondition(DCM_FILECREATEFAILED,
                              DCM_Message(DCM_FILECREATEFAILED),
                              file, strerror(errno),
                              "DCM_WriteFile");
  }
  cond = DCM_ExportStream(callerObject, opt, buf,
                          (unsigned long) sizeof(buf), writeFile, &fd);
  if (cond != DCM_NORMAL)
    return cond;

  (void) close(fd);
  return DCM_NORMAL;
}

/* DCM_ModifyElements
**
** Purpose:
**
** Parameter Dictionary:
**  callerObject    Handle to user's DICOM object to be modified
**  vector      Mandatory elements that need to be stored
**        in the object
**  count     Number of such mandatory elements
**  flaggedVector   Optional elements
**  flaggedCount    Number of such optional elements
**  updateCount   Total number of elements updated (returned to
**        caller)
**
** Return Values:
**
**  DCM_ILLEGALOBJECT
**  DCM_ILLEGALREPRESENTATION
**  DCM_INSERTFAILED
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Algorithm:
**  Check caller's object to make certain it is a legal DCM object.
**  Find head of the object linked list.
**  Search linked list sequentially until object is found or end
**  of list reached.
**  If end of list
**    return DCM_ELEMENTNOTFOUND
**  If CTX pointer containts NULL
**      Begin copy from beginning of element
**  else
**      Begin copy from address in CTX pointer
**  Copy data from element data area to user buffer
**  If copy is incomplete (remaining data longer than caller's buffer)
**      Update CTX pointer to point to next uncopied part of data
**      Return DCM_GETINCOMPLETE
**  else
**      Update CTX pointer to point past data area.
**      Return DCM_NORMAL
*/

CONDITION
DCM_ModifyElements(DCM_OBJECT ** callerObject,
                   DCM_ELEMENT * vector, int count,
                   DCM_FLAGGED_ELEMENT * flaggedVector, int flaggedCount,
                   int *updateCount) {
  PRIVATE_OBJECT
  ** object;
  CONDITION
  cond;
  DCM_ELEMENT
  e;
  int
  c = 0;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_ModifyElement");
  if (cond != DCM_NORMAL)
    return cond;

  while (count-- > 0) {
    cond = DCM_RemoveElement(callerObject, vector->tag);
    if (cond != DCM_NORMAL)
      (void) COND_PopCondition(FALSE);

    e = *vector;
    if (DCM_IsString(e.representation))
      e.length = strlen(e.d.string);

    cond = DCM_AddElement(callerObject, &e);
    if (cond != DCM_NORMAL)
      return cond;

    c++;
    vector++;
  }

  while (flaggedCount-- > 0) {
    if ((*(flaggedVector->flagAddress) & flaggedVector->flag) != 0) {
      cond = DCM_RemoveElement(callerObject, flaggedVector->e.tag);
      if (cond != DCM_NORMAL)
        (void) COND_PopCondition(FALSE);

      e = flaggedVector->e;
      if (DCM_IsString(e.representation))
        e.length = strlen(e.d.string);
      cond = DCM_AddElement(callerObject, &e);
      if (cond != DCM_NORMAL)
        return cond;
      c++;
    }
    flaggedVector++;
  }

  if (updateCount != NULL)
    *updateCount = c;
  return DCM_NORMAL;
}


/* DCM_ParseObject
**
** Purpose:
**  Parse the object and store the mandatory and optional elements in
**  different vectors.
**
** Parameter Dictionary:
**      callerObject            Handle to user's DICOM object to be modified
**      vector                  Mandatory elements that need to be stored
**                              in the object
**      count                   Number of such mandatory elements
**      flaggedVector           Optional elements
**      flaggedCount            Number of such optional elements
**      parseCount              Total number of elements parsed (returned to
**                              caller)
**
** Return Values:
**
**  DCM_CANNOTGETSEQUENCEVALUE
**  DCM_ELEMENTNOTFOUND
**  DCM_GETINCOMPLETE
**  DCM_ILLEGALCONTEXT
**  DCM_ILLEGALOBJECT
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/
CONDITION
DCM_ParseObject(DCM_OBJECT ** callerObject, DCM_ELEMENT * vector,
                int count,
                DCM_FLAGGED_ELEMENT * flaggedVector, int flagCount,
                int *parseCount) {
  PRIVATE_OBJECT
  ** object;
  CONDITION
  cond;
  void
  *ctx;
  U32
  l;
  int
  c = 0;
  char
  *p;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_ParseObject");
  if (cond != DCM_NORMAL)
    return cond;

  while (count-- > 0) {
    ctx = NULL;
    cond = DCM_GetElementValue(callerObject, vector, &l, &ctx);
    if (cond != DCM_NORMAL)
      return cond;
    if (DCM_IsString(vector->representation)) {
      vector->d.string[l] = '\0';
      p = vector->d.string + l - 1;
      while (p >= vector->d.string && (*p == ' '))
        *p-- = '\0';
    }
    c++;
    vector++;
  }

  while (flagCount-- > 0) {
    ctx = NULL;
    cond = DCM_GetElementValue(callerObject, &flaggedVector->e, &l, &ctx);
    if (cond != DCM_NORMAL) {
      (void) COND_PopCondition(FALSE);
    } else {
      c++;
      if (DCM_IsString(flaggedVector->e.representation)) {
        flaggedVector->e.d.string[l] = '\0';
        p = flaggedVector->e.d.string + l - 1;
        while (p >= flaggedVector->e.d.string && (*p == ' '))
          *p-- = '\0';
      }
      *(flaggedVector->flagAddress) |= flaggedVector->flag;
    }
    flaggedVector++;
  }

  if (parseCount != NULL)
    *parseCount = c;
  return DCM_NORMAL;
}


/* DCM_RemoveGroup
**
** Purpose:
**  Remove an element with the given group number from the object
**
** Parameter Dictionary:
**  callerObject    Handle to caller's object
**  group     Group number of the element to be removed.
**
** Return Values:
**
**  DCM_GROUPNOTFOUND
**  DCM_ILLEGALOBJECT
**  DCM_LISTFAILURE
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

CONDITION
DCM_RemoveGroup(DCM_OBJECT ** callerObject, unsigned short group) {
  PRIVATE_OBJECT
  ** object;
  CONDITION
  cond;
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem;
  CTNBOOLEAN
  found = FALSE;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_RemoveGroup");
  if (cond != DCM_NORMAL)
    return cond;

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&(*object)->groupList);
  if (groupItem == NULL)
    return COND_PushCondition(DCM_GROUPNOTFOUND,
                              DCM_Message(DCM_GROUPNOTFOUND),
                              (int) group, "DCM_RemoveGroup");

  (void) LST_Position(&(*object)->groupList, groupItem);

  while (!found && (groupItem != NULL)) {
    if (groupItem->group == group)
      found = TRUE;
    else
      groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
  }
  if (groupItem == NULL)
    return COND_PushCondition(DCM_GROUPNOTFOUND,
                              DCM_Message(DCM_GROUPNOTFOUND),
                              (int) group, "DCM_RemoveGroup");


  while ((elementItem = (PRV_ELEMENT_ITEM*)LST_Pop(&groupItem->elementList)) != NULL)
    CTN_FREE(elementItem);

  groupItem = (PRV_GROUP_ITEM*)LST_Remove(&(*object)->groupList, LST_K_AFTER);
  cond = LST_Destroy(&groupItem->elementList);
  if (cond != LST_NORMAL)
    return COND_PushCondition(DCM_LISTFAILURE,
                              DCM_Message(DCM_LISTFAILURE),
                              "DCM_RemoveGroup");
  CTN_FREE(groupItem);
  return DCM_NORMAL;
}

/* DCM_GetSequenceList
**
** Purpose:
**  Obtain the sequence list from the DICOM object corresponding to the
**  tag value.
**
** Parameter Dictionary:
**  object    Handle to the DICOM object
**  tag   Tag number of the sequence list element to be obtained
**      from the DICOM object
**  list    Holds the sequence list. Returned to the caller.
**
** Return Values:
**
**  DCM_ELEMENTNOTFOUND
**  DCM_ILLEGALOBJECT
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

CONDITION
DCM_GetSequenceList(DCM_OBJECT ** object, DCM_TAG tag, LST_HEAD ** list) {
  PRIVATE_OBJECT
  ** obj;
  CONDITION
  cond;
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem;
  CTNBOOLEAN
  found = FALSE;

  obj = (PRIVATE_OBJECT **) object;
  cond = checkObject(obj, "DCM_GetSequenceList");
  if (cond != DCM_NORMAL)
    return cond;

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&(*obj)->groupList);
  if (groupItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_ELEMENT(tag),
                              "DCM_GetSequenceList");

  (void) LST_Position(&(*obj)->groupList, groupItem);
  while (groupItem != NULL) {
    if (groupItem->group == DCM_TAG_GROUP(tag))
      break;

    groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*obj)->groupList);
  }
  if (groupItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_ELEMENT(tag),
                              "DCM_GetSequenceList");

  elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
  if (elementItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_GROUP(tag),
                              "DCM_GetSequenceTag");

  (void) LST_Position(&groupItem->elementList, elementItem);
  while (!found && (elementItem != NULL)) {
    if (elementItem->element.tag == tag) {
      *list = elementItem->element.d.sq;
      found = TRUE;
    }
    elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
  }
  if (found)
    return DCM_NORMAL;
  else
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_ELEMENT(tag),
                              "DCM_GetSequenceList");
}

CONDITION
DCM_GetSequenceElement(DCM_OBJECT ** object, DCM_TAG top, DCM_ELEMENT * e) {
  PRIVATE_OBJECT **obj;
  CONDITION cond;
  PRV_ELEMENT_ITEM *elementItem;
  DCM_SEQUENCE_ITEM *seqItem;

  // CTNBOOLEAN found = FALSE;

  obj = (PRIVATE_OBJECT **) object;
  cond = checkObject(obj, "DCM_GetSequenceElement");
  if (cond != DCM_NORMAL)
    return cond;

  elementItem = locateElement(obj, top);
  if (elementItem == NULL) {
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(top),
                              DCM_TAG_ELEMENT(top),
                              "DCM_GetElementSequence");
  }
  if (elementItem->element.representation != DCM_SQ) {
    return COND_PushCondition(DCM_UNEXPECTEDREPRESENTATION,
                              DCM_Message(DCM_UNEXPECTEDREPRESENTATION),
                              "DCM_GetSequenceElement", "sequence");
  }
  seqItem = (DCM_SEQUENCE_ITEM*)LST_Head(&elementItem->element.d.sq);
  cond = DCM_ParseObject(&seqItem->object, e, 1, NULL, 0, NULL);
  return cond;

#if 0
  return DCM_NORMAL;
#endif
}

/* DCM_GetElementValueList
**
** Purpose:
**
** Parameter Dictionary:
**  Define the parameters to the function
**
** Return Values:
**
**  DCM_ELEMENTNOTFOUND
**  DCM_ILLEGALOBJECT
**  DCM_LISTFAILURE
**  DCM_MALLOCFAILURE
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/
CONDITION
DCM_GetElementValueList(DCM_OBJECT ** object, DCM_TAG tag,
                        size_t structureSize, long stringOffset,
                        LST_HEAD ** list) {
  PRIVATE_OBJECT
  ** obj;
  CONDITION
  cond;
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem;
  CTNBOOLEAN
  found = FALSE;
  char
  *src,
  *dst,
  *p;
  U32
  l;

  obj = (PRIVATE_OBJECT **) object;
  cond = checkObject(obj, "DCM_GetSequenceList");
  if (cond != DCM_NORMAL)
    return cond;

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&(*obj)->groupList);
  if (groupItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_ELEMENT(tag),
                              "DCM_GetSequenceList");

  (void) LST_Position(&(*obj)->groupList, groupItem);
  while (groupItem != NULL) {
    if (groupItem->group == DCM_TAG_GROUP(tag))
      break;

    groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*obj)->groupList);
  }
  if (groupItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_ELEMENT(tag),
                              "DCM_GetSequenceList");

  elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
  if (elementItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_GROUP(tag),
                              "DCM_GetSequenceTag");

  (void) LST_Position(&groupItem->elementList, elementItem);
  while (!found && (elementItem != NULL)) {
    if (elementItem->element.tag == tag) {
      found = TRUE;
    } else
      elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
  }
  if (!found)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_ELEMENT(tag),
                              "DCM_GetElementValueList");

  if (!DCM_IsString(elementItem->element.representation)) {
    return COND_PushCondition(DCM_UNEXPECTEDREPRESENTATION,
                              DCM_Message(DCM_UNEXPECTEDREPRESENTATION),
                              "DCM_GetElementValueList",
                              "string");
  }
  src = elementItem->element.d.string;
  l = elementItem->element.length;
  while (l > 0) {
    while (l > 1 && (*src == ' ' || *src == DCM_DELIMITOR)) {
      l--;
      src++;
    }
    if ((l == 1) && (*src == ' ' || *src == DCM_DELIMITOR))
      l--;

    if (l != 0) {
      p = (char*)CTN_MALLOC(structureSize);
      if (p == NULL)
        return COND_PushCondition(DCM_MALLOCFAILURE,
                                  DCM_Message(DCM_MALLOCFAILURE),
                                  structureSize,
                                  "DCM_GetElementValueList");
      dst = p + stringOffset;
      while ((l > 1) && (*src != DCM_DELIMITOR)) {
        *dst++ = *src++;
        l--;
      }
      if ((l == 1) && (*src != ' ')) {
        *dst++ = *src++;
        l--;
      }
      *dst = '\0';
      ;
      cond = LST_Enqueue(list, p);
      if (cond != LST_NORMAL)
        return COND_PushCondition(DCM_LISTFAILURE,
                                  DCM_Message(DCM_LISTFAILURE),
                                  "DCM_GetElementValueList");
    }
  }
  return DCM_NORMAL;
}

/* DCM_AddElementList
**
** Purpose:
**  Add an element list to the DICOM object
**
** Parameter Dictionary:
**  callerObject    Handle to object to which the element is to be
**        added
**  element     The element in which the string obtained from
**        the list is to be stored. Finally the element
**        is added to the DICOM object.
**  list      List of structures , each containing a string
**        starting at some offset specified by the
**        parameter "offset"
**  offset      Offset in each individual structure (see
**        explanation for parameter list)
**
** Return Values:
**
**  DCM_ILLEGALOBJECT
**  DCM_ILLEGALREPRESENTATION
**  DCM_INSERTFAILED
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/
CONDITION
DCM_AddElementList(DCM_OBJECT ** callerObject, DCM_ELEMENT * element,
                   LST_HEAD * list, long offset) {
  DCM_ELEMENT
  e;      /* Local copy of caller's element */
  CONDITION
  cond;
  char
  *s;

  e = *element;
  cond = DCM_ListToString(list, offset, &s);
  if (cond != DCM_NORMAL)
    return cond;

  e.d.string = s;
  e.length = strlen(s);
  cond = DCM_AddElement(callerObject, &e);
  CTN_FREE(s);
  return cond;
}

/* DCM_GetElement
**
** Purpose:
**  Get the element with the specified tag number from the given DICOM
**  object
**
** Parameter Dictionary:
**  callerObject    Handle to the DICOM object
**  tag     Tag number of the element to be obtained
**        from the object
**  element     The element to be returned
**
** Return Values:
**
**  DCM_ELEMENTNOTFOUND
**  DCM_ILLEGALOBJECT
**  DCM_NORMAL
**  DCM_NULLOBJECT
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/
CONDITION
DCM_GetElement(DCM_OBJECT ** callerObject, DCM_TAG tag, DCM_ELEMENT * element) {
  PRIVATE_OBJECT
  ** obj;
  CONDITION
  cond;
  PRV_ELEMENT_ITEM
  * elementItem;

  obj = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(obj, "DCM_GetElementVM");
  if (cond != DCM_NORMAL)
    return cond;

  elementItem = locateElement(obj, tag);
  if (elementItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_ELEMENT(tag),
                              "DCM_GetElementVM");
  *element = elementItem->element;
  element->d.ot = NULL;
  return DCM_NORMAL;
}

CONDITION
DCM_ComputeExportLength(DCM_OBJECT ** callerObject, unsigned long opt,
                        unsigned long *length) {
  PRIVATE_OBJECT
  ** object;
  unsigned char
  buf[2048];
  CONDITION
  cond;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_ComputeExportSize");
  if (cond != DCM_NORMAL)
    return cond;

  *length = 0;
  cond = DCM_ExportStream(callerObject, opt, buf,
                          (unsigned long) sizeof(buf), countBytes, length);
  if (cond != DCM_NORMAL)
    return cond;

  return DCM_NORMAL;
}

CONDITION
DCM_CompareAttributes(DCM_OBJECT ** o1, DCM_OBJECT ** o2,
                      void (*callback) (const DCM_ELEMENT * e1,
                                        const DCM_ELEMENT * e2,
                                        void *ctx),
                      void *ctx) {
  PRIVATE_OBJECT **object1,
  **object2;
  PRV_GROUP_ITEM *groupItem1,
  *groupItem2;
  CONDITION cond;

  object1 = (PRIVATE_OBJECT **) o1;
  cond = checkObject(object1, "DCM_CompareAttributes");
  if (cond != DCM_NORMAL)
    return cond;

  object2 = (PRIVATE_OBJECT **) o2;
  cond = checkObject(object1, "DCM_CompareAttributes");
  if (cond != DCM_NORMAL)
    return cond;

  groupItem1 = (PRV_GROUP_ITEM*)LST_Head(&(*object1)->groupList);
  if (groupItem1 != NULL)
    (void) LST_Position(&(*object1)->groupList, groupItem1);

  groupItem2 = (PRV_GROUP_ITEM*)LST_Head(&(*object2)->groupList);
  if (groupItem2 != NULL)
    (void) LST_Position(&(*object2)->groupList, groupItem2);


  while (groupItem1 != NULL) {
    if (groupItem2 == NULL) {
      compareGroup(groupItem1, NULL, callback, ctx);
      groupItem1 = (PRV_GROUP_ITEM*)LST_Next(&(*object1)->groupList);
    } else if (groupItem1->group == groupItem2->group) {
      compareGroup(groupItem1, groupItem2, callback, ctx);
      groupItem1 = (PRV_GROUP_ITEM*)LST_Next(&(*object1)->groupList);
      groupItem2 = (PRV_GROUP_ITEM*)LST_Next(&(*object2)->groupList);
    } else if (groupItem1->group > groupItem2->group) {
      compareGroup(NULL, groupItem2, callback, ctx);
      groupItem2 = (PRV_GROUP_ITEM*)LST_Next(&(*object2)->groupList);
    } else {
      compareGroup(groupItem1, NULL, callback, ctx);
      groupItem1 = (PRV_GROUP_ITEM*)LST_Next(&(*object1)->groupList);
    }
  }

  while (groupItem2 != NULL) {
    compareGroup(NULL, groupItem2, callback, ctx);
    groupItem2 = (PRV_GROUP_ITEM*)LST_Next(&(*object2)->groupList);
  }
  return DCM_NORMAL;
}

CTNBOOLEAN
DCM_GroupPresent(DCM_OBJECT ** o1, U16 group) {
  PRIVATE_OBJECT **object;
  PRV_GROUP_ITEM * item;
  CONDITION cond;
  CTNBOOLEAN tooFar = FALSE;

  object = (PRIVATE_OBJECT **) o1;
  cond = checkObject(object, "DCM_CompareAttributes");
  if (cond != DCM_NORMAL)
    return FALSE;


  item = (PRV_GROUP_ITEM*)LST_Head(&(*object)->groupList);
  if (item != NULL)
    (void) LST_Position(&(*object)->groupList, item);

  while (item != NULL && !tooFar) {
    if (item->group == group) {
      return TRUE;
    } else if (item->group > group) {
      tooFar = TRUE;
    } else {
      item = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
    }
  }
  return FALSE;
}

/*     ------------------------------------------------------------
**  Private functions below here
*/

/* newElementItem
**
** Purpose:
**  Create a new element item suitable for placing in the linked
**  list representation of an ACR object.  Copy data from an
**  existing element, but skip the actual data field.
**  Describe the purpose of the function
**
** Parameter Dictionary:
**  src Pointer to source element that is to be copied
**  dst Pointer to pointer to destination element which is allocated
**    by this routine and filled in appropriately.
**
** Return Values:
**  DCM_NORMAL
**  DCM_ELEMENTCREATEFAILED
**
** Algorithm:
**  Allocate new element item of size:
**    Size PRV_ELEMENT_ITEM + length of data value
**  Copy data from caller's DCM_ELEMENT into newly created
**    PRV_ELEMENT_ITEM.
**  Point data value of newly created PRV_ELEMENT_ITEM to part of the
**  allocated space (just past the end of the PRV_ELEMENT_ITEM).
*/
static CONDITION
newElementItem(DCM_ELEMENT * src, CTNBOOLEAN allocateData,
               PRV_ELEMENT_ITEM ** dst) {
  long l;

  if (allocateData && (src->representation != DCM_SQ)) {
    l = src->length;
    if (l & 1)
      l++;
  } else
    l = 0;

  if (debug) {
    fprintf(stderr, "newElementItem: CTN_MALLOC %8d %8d ", (int)l,
            (int)(sizeof(PRV_ELEMENT_ITEM) + l));
  }
  *dst = (PRV_ELEMENT_ITEM *) CTN_MALLOC(sizeof(PRV_ELEMENT_ITEM) + l);
  if (debug)
    fprintf(stderr, "%8lx\n", (unsigned long) *dst);

  if (*dst == NULL) {
    return COND_PushCondition(DCM_ELEMENTCREATEFAILED,
                              DCM_Message(DCM_ELEMENTCREATEFAILED),
                              "newElementItem",
                              DCM_TAG_GROUP(src->tag),
                              DCM_TAG_ELEMENT(src->tag),
                              l);
  }
  memset(*dst, 0, sizeof(PRV_ELEMENT_ITEM));
  (*dst)->element = *src;
  (*dst)->byteOrder = NATIVE_ORDER;
  (*dst)->allocatedDataLength = (size_t) l;
  (*dst)->originalDataLength = src->length;
  (*dst)->paddedDataLength = src->length;
  if (allocateData)
    (*dst)->element.d.ot = ((char *) (*dst)) + sizeof(PRV_ELEMENT_ITEM);
  else
    (*dst)->element.d.ot = NULL;

  (*dst)->fragmentFlag = 0;
  return DCM_NORMAL;
}

/* findCreateGroup
**
** Purpose:
**  Find the group in the DCM object corresponding to the group
**  passed by the caller.  If the group does not yet exist, create
**  a new group.  Set the CURRENT pointer in the linked list
**  to point at that group.
**
** Parameter Dictionary:
**  object    Pointer to caller's DCM object
**  group   Group number to locate/create
**  groupPtr  Mechanism for returning pointer to located group
**
** Return Values:
**
**  DCM_ELEMENTCREATEFAILED
**  DCM_LISTFAILURE
**  DCM_NORMAL
**
** Algorithm:
**  Set ITEM to head of linked list of ACR object
**  Set CURRENT item in linked list to ITEM
**  Search sequentially through linked list until:
**      - Reach exisiting group that matches caller's group
**      - Reach a group with larger group number than caller's group
**      - Reach end of linked list
**  Each time you move to a new item, update CURRENT to point to that item
**  If reached existing group
**      return
**  If reached a group with larger group number than caller's group,
**      Insert new group with Group Length Element (0000) just before
**      the group with the larger group number.
**      Set CURRENT pointer in linked list to point at new group
**      If group is COMMAND or IDENTIFYING,
**    Insert Length to End Element
**      Return
**  If reached end of the linked list
**      Append new group with Group Length Element (0000) to the end
**      of the linked list.
**      Set CURRENT pointer in linked list to point at new group
**      If group is COMMAND or IDENTIFYING,
**    Insert Length to End Element
**      Return
**
*/

static CONDITION
findCreateGroup(PRIVATE_OBJECT ** object, unsigned short group,
                PRV_GROUP_ITEM ** groupItem) {
  PRV_GROUP_ITEM
  * item;
  CONDITION
  cond;
  CTNBOOLEAN
  tooFar = FALSE;

  item = (PRV_GROUP_ITEM*)LST_Head(&(*object)->groupList);
  if (item != NULL)
    (void) LST_Position(&(*object)->groupList, item);

  while (item != NULL && !tooFar) {
    if (item->group == group) {
      *groupItem = item;
      return DCM_NORMAL;
    } else if (item->group > group) {
      tooFar = TRUE;
    } else {
      item = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
    }
  }

  {
    U32 l;
    PRV_GROUP_ITEM *newGroupItem;
    DCM_ELEMENT groupLength = {0, DCM_UL, "", 1, sizeof(l),{NULL}};
    PRV_ELEMENT_ITEM *groupLengthItem;

    newGroupItem = (PRV_GROUP_ITEM*)CTN_MALLOC(sizeof(*newGroupItem));
    if (newGroupItem == NULL)
      return COND_PushCondition(DCM_ELEMENTCREATEFAILED,
                                DCM_Message(DCM_ELEMENTCREATEFAILED),
                                "findCreateGroup",
                                group, 0xffff, sizeof(*newGroupItem));


    *groupItem = newGroupItem;
    newGroupItem->group = group;
    newGroupItem->baseLength = 0;
    newGroupItem->longVRAttributes = 0;
    newGroupItem->elementList = LST_Create();
    if (newGroupItem->elementList == NULL)
      return COND_PushCondition(DCM_LISTFAILURE,
                                DCM_Message(DCM_LISTFAILURE),
                                "findCreateGroup");

    if (tooFar)
      cond = LST_Insert(&(*object)->groupList, newGroupItem, LST_K_BEFORE);
    else
      cond = LST_Enqueue(&(*object)->groupList, newGroupItem);
    if (cond != LST_NORMAL)
      return COND_PushCondition(DCM_LISTFAILURE,
                                DCM_Message(DCM_LISTFAILURE),
                                "findCreateGroup");
    (void) LST_Position(&(*object)->groupList, newGroupItem);
    if (cond != LST_NORMAL)
      return COND_PushCondition(DCM_LISTFAILURE,
                                DCM_Message(DCM_LISTFAILURE),
                                "findCreateGroup");

    groupLength.d.ul = &l;
    l = 0;
    if ((*object)->groupLengthFlag) {
      groupLength.tag = DCM_MAKETAG(group, 0);
      cond = newElementItem(&groupLength, TRUE, &groupLengthItem);
      (void) memcpy(groupLengthItem->element.d.ot, &l, sizeof(l));

      if (LST_Insert(&newGroupItem->elementList,
                     groupLengthItem,
                     LST_K_AFTER) !=
          LST_NORMAL)
        return COND_PushCondition(DCM_LISTFAILURE,
                                  DCM_Message(DCM_LISTFAILURE),
                                  "findCreateGroup");

      (*object)->objectSize += 8 + groupLengthItem->element.length;
    }
  }
  return DCM_NORMAL;
}

/* insertNewElement
**
** Purpose:
**  Create a new DCM_ELEMENT item using a copy of the caller's
**  DCM_ELEMENT and insert it into the ACR object's linked list.
**
** Parameter Dictionary:
**  object    Pointer to caller's ACR_OBJECT
**  element   Pointer to caller's ACR_ELEMENT to be copied
**      and inserted into linked list.
**
** Return Values:
**
**  DCM_ELEMENTCREATEFAILED
**  DCM_LISTFAILURE
**  DCM_NORMAL
**  DCM_UNEVENELEMENTLENGTH
**  DCM_BADELEMENTINGROUP
**
** Algorithm:
**  Call newElementItem to create a copy of the DCM_ELEMENT
**  Copy caller's data into data area allocated by newElementItem
**  Increment object's OBJECTSIZE field by size of new element
**  Use CURRENT pointer in DCM object linked list to get pointer
**  to the group where we insert this element
**  Update Group Length by adding size of new element
**  Search sequentially through linked list until we reach:
**      - End of linked list
**      - A different group
**      - An element in the same group with a larger element number
**  If reached end of linked list
**      Append new ACR_ELEMENTITEM to end of linked list
**  If reached a different group
**      Insert new ACR_ELEMENTITEM just before new group
**  If reached an element in the same group with a larger element number
**      Insert new ACR_ELEMENTITEM just before the "larger" element
*/
static CONDITION
insertNewElement(PRIVATE_OBJECT ** object, DCM_ELEMENT * element) {
  PRV_ELEMENT_ITEM
  * nextItem,
  *newItem;
  PRV_GROUP_ITEM
  * groupItem;
  CONDITION
  cond;
  char
  *p;

  cond = newElementItem(element, TRUE, &newItem);
  if (cond != DCM_NORMAL) {
    return cond;
  }
  newItem->byteOrder = DCM_ORDERNATIVE;
  if ((newItem->element.length & 1) != 0) {
    if (newItem->element.representation == DCM_AE) {
      p = newItem->element.d.string;  /* repair, check for 16 */
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_AS) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_CS) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_DA) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_DS) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_IS) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_LT) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_LO) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_PN) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_SH) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_ST) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_TM) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_UI) {
      p = newItem->element.d.string;
      p[newItem->element.length] = '\0';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_UT) {
      p = newItem->element.d.string;
      p[newItem->element.length] = ' ';
      newItem->paddedDataLength = element->length + 1;
      (void) memcpy(p, element->d.string, element->length);
    } else if (newItem->element.representation == DCM_SQ) {
      /*      newItem->element.length = 0xffffffff; */
      newItem->element.d.sq = element->d.sq;
    } else {
      CTN_FREE(newItem);
      return COND_PushCondition(DCM_UNEVENELEMENTLENGTH,
                                DCM_Message(DCM_UNEVENELEMENTLENGTH),
                                DCM_TAG_GROUP(element->tag),
                                DCM_TAG_ELEMENT(element->tag),
                                element->length,
                                "insertNewElement");
    }
  } else if (newItem->element.representation != DCM_SQ) {
    (void) memcpy(newItem->element.d.ot, element->d.ot, element->length);
  } else {
    /*  newItem->element.length = 0xffffffff; */
    newItem->element.d.sq = element->d.sq;
  }
  if ((*object)->objectSize != DCM_UNSPECIFIEDLENGTH)
    (*object)->objectSize += 8 + newItem->paddedDataLength;

  /* repair */
  cond = updateSpecialElements(object, newItem);
  if (cond != DCM_NORMAL)
    return cond;

  groupItem = (PRV_GROUP_ITEM*)LST_Current(&(*object)->groupList);
  if (groupItem == NULL)
    return COND_PushCondition(DCM_LISTFAILURE,
                              DCM_Message(DCM_LISTFAILURE),
                              "insertNewElement");

  if (groupItem->baseLength != DCM_UNSPECIFIEDLENGTH)
    groupItem->baseLength += 2 + 2 + 4 + newItem->paddedDataLength;

  if (newItem->element.representation == DCM_OW ||
      newItem->element.representation == DCM_OB ||
      newItem->element.representation == DCM_SQ) {
    groupItem->longVRAttributes++;
    (*object)->longVRAttributes++;
  }
  if ((nextItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList)) == NULL) {
    cond = LST_Enqueue(&groupItem->elementList, newItem);
    if (cond != LST_NORMAL)
      return COND_PushCondition(DCM_LISTFAILURE,
                                DCM_Message(DCM_LISTFAILURE),
                                "insertNewElement");
    else
      return DCM_NORMAL;
  }
  (void) LST_Position(&groupItem->elementList, nextItem);
  if (DCM_TAG_ELEMENT(nextItem->element.tag) == 0x0000)
    (void) memcpy(nextItem->element.d.ot, &groupItem->baseLength,
                  sizeof(groupItem->baseLength));

  /*  Now, search through the linked list for a place to insert/append
  **  this new item.
  */

  while (nextItem != NULL) {
    if (DCM_TAG_GROUP(element->tag) !=
        DCM_TAG_GROUP(nextItem->element.tag)) {
      return COND_PushCondition(DCM_BADELEMENTINGROUP,
                                DCM_Message(DCM_BADELEMENTINGROUP),
                                DCM_TAG_GROUP(nextItem->element.tag),
                                DCM_TAG_ELEMENT(nextItem->element.tag),
                                groupItem->group, "insertNewElement");
    } else if (DCM_TAG_ELEMENT(element->tag) <
               DCM_TAG_ELEMENT(nextItem->element.tag)) {
      cond = LST_Insert(&groupItem->elementList, newItem, LST_K_BEFORE);
      if (cond != LST_NORMAL)
        return COND_PushCondition(DCM_LISTFAILURE,
                                  DCM_Message(DCM_LISTFAILURE),
                                  "insertNewElement");
      else
        return DCM_NORMAL;
    }
    nextItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
  }

  /*  If we fall out of the loop, we must have reached the end of
  **  the group.  Add the element to the end of the list of elements
  **  in this group.
  */

  cond = LST_Enqueue(&groupItem->elementList, newItem);
  if (cond != LST_NORMAL)
    return COND_PushCondition(DCM_LISTFAILURE,
                              DCM_Message(DCM_LISTFAILURE),
                              "insertNewElement");
  else
    return DCM_NORMAL;
}

static CONDITION
insertThisElementItem(PRIVATE_OBJECT ** object, PRV_ELEMENT_ITEM* newItem) {
  PRV_ELEMENT_ITEM * nextItem;
  PRV_GROUP_ITEM * groupItem = 0;
  CONDITION cond;

  /* repair */
  cond = updateSpecialElements(object, newItem);
  if (cond != DCM_NORMAL)
    return cond;

  cond = findCreateGroup(object, DCM_TAG_GROUP(newItem->element.tag),
                         &groupItem);

  if (groupItem == NULL)
    return COND_PushCondition(DCM_LISTFAILURE,
                              DCM_Message(DCM_LISTFAILURE),
                              "insertThisElementItem");

  if (groupItem->baseLength != DCM_UNSPECIFIEDLENGTH)
    groupItem->baseLength += 2 + 2 + 4 + newItem->paddedDataLength;

  if (newItem->element.representation == DCM_OW ||
      newItem->element.representation == DCM_OB ||
      newItem->element.representation == DCM_SQ) {
    groupItem->longVRAttributes++;
    (*object)->longVRAttributes++;
  }

  if ((nextItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList)) == NULL) {
    cond = LST_Enqueue(&groupItem->elementList, newItem);
    if (cond != LST_NORMAL)
      return COND_PushCondition(DCM_LISTFAILURE,
                                DCM_Message(DCM_LISTFAILURE),
                                "insertThisElementItem");
    else
      return DCM_NORMAL;
  }

  (void) LST_Position(&groupItem->elementList, nextItem);
  if (DCM_TAG_ELEMENT(nextItem->element.tag) == 0x0000)
    (void) memcpy(nextItem->element.d.ot, &groupItem->baseLength,
                  sizeof(groupItem->baseLength));

  /*  Now, search through the linked list for a place to insert/append
  **  this new item.
  */

  while (nextItem != NULL) {
    if (DCM_TAG_GROUP(newItem->element.tag) !=
        DCM_TAG_GROUP(nextItem->element.tag)) {
      return COND_PushCondition(DCM_BADELEMENTINGROUP,
                                DCM_Message(DCM_BADELEMENTINGROUP),
                                DCM_TAG_GROUP(nextItem->element.tag),
                                DCM_TAG_ELEMENT(nextItem->element.tag),
                                groupItem->group, "insertThisElementItem");
    } else if (DCM_TAG_ELEMENT(newItem->element.tag) <
               DCM_TAG_ELEMENT(nextItem->element.tag)) {
      cond = LST_Insert(&groupItem->elementList, newItem, LST_K_BEFORE);
      if (cond != LST_NORMAL)
        return COND_PushCondition(DCM_LISTFAILURE,
                                  DCM_Message(DCM_LISTFAILURE),
                                  "insertThisElementItem");
      else
        return DCM_NORMAL;
    }
    nextItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
  }

  /*  If we fall out of the loop, we must have reached the end of
  **  the group.  Add the element to the end of the list of elements
  **  in this group.
  */

  cond = LST_Enqueue(&groupItem->elementList, newItem);
  if (cond != LST_NORMAL)
    return COND_PushCondition(DCM_LISTFAILURE,
                              DCM_Message(DCM_LISTFAILURE),
                              "insertThisElementItem");
  else
    return DCM_NORMAL;
}

/* updateObjectType
**
** Purpose:
**  Possibly modify the objectType field of an DCM object to identify
**  the object as COMMAND, DATASET or MESSAGE.
**
** Parameter Dictionary:
**  object    Pointer to caller's PRIVATE object to be updated
**  element   Pointer to DCM_ELEMENT which will be added to
**      the object and possibly cause a change in the
**      type of the object.
**
** Return Values:
**  DCM_NORMAL
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

static CONDITION
updateObjectType(PRIVATE_OBJECT ** object, DCM_ELEMENT * element) {
  switch ((*object)->objectType) {
  case DCM_OBJECTUNKNOWN:
    if (DCM_TAG_GROUP(element->tag) == DCM_GROUPCOMMAND)
      (*object)->objectType = DCM_OBJECTCOMMAND;
    else
      (*object)->objectType = DCM_OBJECTELEMENTLIST;
    break;
  case DCM_OBJECTCOMMAND:
    if (DCM_TAG_GROUP(element->tag) != DCM_GROUPCOMMAND)
      (*object)->objectType = DCM_OBJECTELEMENTLIST;
    break;
  case DCM_OBJECTELEMENTLIST:
  case DCM_OBJECTIMAGE:
    break;
  default:
    break;
  }
  return DCM_NORMAL;
}

/* updateSpecialElements
**
** Purpose:
**  Update special elements in a DICOM object when a new data element
**  is added to the object.  These special fields are used by other
**  parts of the package which have to refer to those fields and wish
**  to do so without searching through the entire list.  This could
**  get messy and is a candidate for redesign.
**
** Parameter Dictionary:
**  object    Pointer to caller's PRIVATE DICOM object
**  element   Pointer to DCM_ELEMENT that is being added to
**      the DICOM object
**
** Return Values:
**  DCM_NORMAL
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/
static CONDITION
updateSpecialElements(PRIVATE_OBJECT ** object,
                      PRV_ELEMENT_ITEM * item) {
  int idx;
#ifdef BIG_ENDIAN_ARCHITECTURE
  unsigned char* b1;
  unsigned char* b2;
#endif

  switch (item->element.tag) {
  case DCM_IMGBITSALLOCATED:
    (*object)->pixelBitsAllocated = *item->element.d.us;
    break;
  case DCM_IMGPIXELREPRESENTATION:
    (*object)->pixelRepresentation = *item->element.d.us;
    break;
  case DCM_METAGROUPLENGTH:
    (*object)->metaHeaderLength = *item->element.d.ul;
    //printf("(*object)->metaHeaderLength: %d",(*object)->metaHeaderLength);
#ifdef BIG_ENDIAN_ARCHITECTURE
    b1 = (unsigned char *)&(*object)->metaHeaderLength;
    b2 = (unsigned char *)item->element.d.ul;
    *b1++ = b2[2];
    *b1++ = b2[3];
    *b1++ = b2[0];
    *b1++ = b2[1];
#endif
    break;
  case DCM_METATRANSFERSYNTAX:
    if (strcmp(item->element.d.string, DICOM_TRANSFERLITTLEENDIAN) == 0) {
      (*object)->dataOptions = DCM_ORDERLITTLEENDIAN;
    } else if (strcmp(item->element.d.string,
                      DICOM_TRANSFERLITTLEENDIANEXPLICIT) == 0) {
      (*object)->dataOptions = DCM_EXPLICITLITTLEENDIAN;
    } else if (strcmp(item->element.d.string,
                      DICOM_TRANSFERBIGENDIANEXPLICIT) == 0) {
      (*object)->dataOptions = DCM_EXPLICITBIGENDIAN;
    } else {  /* Must be an encapsulated transfer syntax */
      (*object)->dataOptions = DCM_EXPLICITLITTLEENDIAN;
    }
    break;
  case DCM_MAKETAG(0x003a, 0x0103):
          strncpy((*object)->waveformDataVR, item->element.d.string,
                  item->element.length);
    (*object)->waveformDataVR[item->element.length] = '\0';
    idx = item->element.length - 1;
    while (idx >= 0 && (*object)->waveformDataVR[idx] == ' ') {
      (*object)->waveformDataVR[idx] = '\0';
      idx--;
    }
    break;
  default:
    break;
  }
  return DCM_NORMAL;
}

typedef struct {
  DCM_VALUEREPRESENTATION representation;
  char code[3];
}
VRMAP;

static VRMAP vrMap[] = {
                         {
                           DCM_AE, "AE"
                         },
                         {DCM_AS, "AS"},
                         {DCM_AT, "AT"},
                         {DCM_CS, "CS"},
                         {DCM_DA, "DA"},
                         {DCM_DD, "DD"},
                         {DCM_DS, "DS"},
                         {DCM_FD, "FD"},
                         {DCM_FL, "FL"},
                         {DCM_IS, "IS"},
                         {DCM_LO, "LO"},
                         {DCM_LT, "LT"},
                         {DCM_OT, "OT"},
                         {DCM_SH, "SH"},
                         {DCM_SL, "SL"},
                         {DCM_SQ, "SQ"},
                         {DCM_SS, "SS"},
                         {DCM_ST, "ST"},
                         {DCM_TM, "TM"},
                         {DCM_UI, "UI"},
                         {DCM_UL, "UL"},
                         {DCM_UN, "UN"},
                         {DCM_US, "US"},
                         {DCM_UT, "UT"},
                         /*{DCM_UNKNOWN, "UK"},*/
                         {DCM_RET, "RT"},
                         {DCM_CTX, "  "},
                         {DCM_PN, "PN"},
                         {DCM_OB, "OB"},
                         {DCM_OW, "OW"},
                         {DCM_DT, "DT"},
                         {DCM_DLM, ""}
                       };

static VRMAP *
lookupVRCode(const char *code) {
  int i;

  for (i = 0; i < (int) DIM_OF(vrMap); i++) {
    if (strcmp(code, vrMap[i].code) == 0)
      return &vrMap[i];
  }

  return NULL;
}

static void
mapVRtoASCII(DCM_VALUEREPRESENTATION vr, char *s) {
  int i;

  for (i = 0; i < (int) DIM_OF(vrMap); i++) {
    if (vr == vrMap[i].representation) {
      strcpy(s, vrMap[i].code);
      return;
    }
  }

  strcpy(s, "");
  return;
}

static void
exportVRLength(DCM_ELEMENT * e, unsigned char *b, int byteOrder,
               U32 * rtnLength) {
  const char *c = "xx";
  unsigned char *p;
  U16 shortLength;
  DCM_VALUEREPRESENTATION vr;

  vr = e->representation;
  if (e->tag == DCM_MAKETAG(0x003a, 0x1000))
    vr = DCM_OB;

  { unsigned int i;
    for (i = 0; i < DIM_OF(vrMap); i++) {
      if (vr == vrMap[i].representation) {
        c = vrMap[i].code;
        break;
      }
    }
  }

  *b++ = *c++;
  *b++ = *c++;
  *rtnLength += 2;

  if (vr == DCM_OB || vr == DCM_OW || vr == DCM_SQ || vr == DCM_UN) {
    *b++ = 0x00;
    *b++ = 0x00;
    if (byteOrder == BYTEORDER_SAME) {
      p = (unsigned char *) &e->length;
      *b++ = *p++;
      *b++ = *p++;
      *b++ = *p++;
      *b++ = *p++;
    } else {
      p = (unsigned char *) &e->length;
      *b++ = p[3];
      *b++ = p[2];
      *b++ = p[1];
      *b++ = p[0];
    }
    *rtnLength += 6;
  } else {
    shortLength = (U16) e->length;
    if (byteOrder == BYTEORDER_SAME) {
      p = (unsigned char *) &shortLength;
      *b++ = *p++;
      *b++ = *p++;
    } else {
      p = (unsigned char *) &shortLength;
      *b++ = p[1];
      *b++ = p[0];
    }
    *rtnLength += 2;
  }
}

static CONDITION
exportPreamble(PRIVATE_OBJECT ** obj, unsigned char *dst,
               U32 bufferLength, U32 * rtnLength) {
  *rtnLength = 0;
  if (bufferLength < (DCM_PREAMBLELENGTH + 4))
    return COND_PushCondition(DCM_EXPORTBUFFERTOOSMALL,
                              DCM_Message(DCM_EXPORTBUFFERTOOSMALL),
                              (int) bufferLength,
                              "exportPreamble");

  (void) memcpy(dst, (*obj)->preamble, DCM_PREAMBLELENGTH);
  dst += DCM_PREAMBLELENGTH;
  (void) memcpy(dst, "DICM", 4);
  *rtnLength += DCM_PREAMBLELENGTH + 4;

  return DCM_NORMAL;
}

/* exportFixedFields
**
** Purpose:
**  This function exports the fixed length fields of an DCM_ELEMENT
**  to the caller's buffer if there is sufficient space in the
**  caller's buffer.
**
** Parameter Dictionary:
**  element   Pointer to the actual data element to be exported
**  b   Pointer to the caller's buffer to hold exported data
**  length    Length of the remaining space in the caller's
**      buffer
**  byteOrder flag giving the order of the bytes as they are
**      exported.  Should be one of:
**        BYTEORDER_SAME
**        BYTEORDER_REVERSE
**  rtnLength Pointer to caller variable to hold the length
**      of the data exported.  The length of the data
**      exported will be 0 if the caller's buffer is
**      too small to hold the fixed length fields.
**
** Return Values:
**  None
**
** Algorithm:
**  If caller buffer is too small to hold all fixed length fields
**      Place 0 in caller's rtnLength variable
**      return
**  Else
**      If byteOrder is the same
**    Copy fixed length fields in same byte order
**      Else
**    Copy fixed length fields in reverse byte order
**      Set caller's rtnLength variable to 8 (short, short, long)
*/

static void
exportFixedFields(DCM_ELEMENT * e,
                  unsigned char *b, U32 length, int byteOrder,
                  CTNBOOLEAN explicitVR, U32 * rtnLength) {
  unsigned char
  *p;
  unsigned short
  group,
  element;
  U32
  minimumLength;

  group = DCM_TAG_GROUP(e->tag);
  element = DCM_TAG_ELEMENT(e->tag);
  if (e->representation == DCM_DLM)
    explicitVR = FALSE;

  minimumLength = sizeof(group) + sizeof(element) + sizeof(e->length);
  if (explicitVR)
    minimumLength += 4;

  *rtnLength = 0;
  if (length >= minimumLength) {
    if (byteOrder == BYTEORDER_SAME) {
      p = (unsigned char *) &group;
      *b++ = *p++;
      *b++ = *p++;
      p = (unsigned char *) &element;
      *b++ = *p++;
      *b++ = *p++;
      *rtnLength += 4;
      if (explicitVR) {
        exportVRLength(e, b, byteOrder, rtnLength);
      } else {
        p = (unsigned char *) &e->length;
        *b++ = *p++;
        *b++ = *p++;
        *b++ = *p++;
        *b++ = *p++;
        *rtnLength += 4;
      }
    } else {
      p = (unsigned char *) &group;
      *b++ = p[1];
      *b++ = p[0];
      p = (unsigned char *) &element;
      *b++ = p[1];
      *b++ = p[0];
      *rtnLength += 4;
      if (explicitVR) {
        exportVRLength(e, b, byteOrder, rtnLength);
      } else {
        p = (unsigned char *) &e->length;
        *b++ = p[3];
        *b++ = p[2];
        *b++ = p[1];
        *b++ = p[0];
        *rtnLength += 4;
      }
    }
  }
}

/* exportData
**
** Purpose:
**  Export the data part of a DCM_ELEMENT.  This function exports
**  all or part of the data portion of an DCM_ELEMENT in the byte order
**  requested by the caller.  The caller specifies the byte order
**  in a formal argument.  The function uses context information to
**  know where to start the export in one data element.  The function
**  does not update the context information but does return the
**  number of bytes exported.
**
** Parameter Dictionary:
**  object    Pointer to the caller's ACR object which is
**      being exported.
**  element   Pointer to the ACR_ELEMENT that is being exported
**  b   Pointer to the caller's buffer to hold the
**      exported data.
**  length    Length of the caller's buffer to hold the data.
**  byteOrder Flag giving the order of the bytes in the exported
**      stream.  Flag should be one of:
**          BYTEORDER_SAME
**          BYTEORDER_REVERSE
**  rtnLength Pointer to caller variable to hold number of bytes
**      that are actually exported.
**
** Return Values:
**
**  DCM_FILEACCESSERROR
**  DCM_NORMAL
**
** Algorithm
**
**  Set caller's rtnLength variable to 0
**  Export data based on representation of data element
**  CASE 16 bit binary:
**      While (length >= 2)
**    If (byte order is same OR element is 8 bit pixel data)
**        Copy 2 bytes to output area
**        Increment input/output pointers by 2
**    Else
**        Copy and swap 2 bytes to output area
**        Increment input/output pointers by 2
**    Endif
**    Decrement length by 2
**    Increment caller's rtnLength by 2
**      End while
**
**  CASE 32 bit binary:
**      While (length >= 4)
**    If (byte order is same)
**        Copy 4 bytes to output area
**        Increment input/output pointers by 4
**    Else
**        Copy and swap 4 bytes to output area
**        Increment input/output pointers by 4
**    Endif
**    Decrement length by 4
**    Increment caller's rtnLength by 4
**
**  CASE ascii text, ascii numeric, or unknown:
**      Use memcpy to copy as of the remaining data as will fit
**    in the caller's buffer.
**      Set caller's rtnLength to the amount of data copied.
**
*/
union {
  unsigned short sh[2];
  unsigned char ch[4];
}
static groupElement = {};

static CONDITION
exportData(PRIVATE_OBJECT ** object, PRV_ELEMENT_ITEM * item,
           unsigned char *src,
           unsigned char *b, U32 length, int byteOrder,
           U32 * rtnLength) {
  /* repair OT for pixel data*/
  unsigned char
  *p;
  DCM_TAG
  * tag;
  DCM_ELEMENT
  * element;
  int nBytes;
  CONDITION cond;

  element = &item->element;

  *rtnLength = 0;
  if (element->d.ot == NULL) {
    if ((*object)->fd != -1) {
      (void) lseek((*object)->fd, item->currentOffset, SEEK_SET);
      nBytes = read((*object)->fd, b, (int) length);
    } else {
      (*object)->sk((*object)->userCtx, item->currentOffset, SEEK_SET);
      cond = (*object)->rd((*object)->userCtx, b, (long) length, &nBytes);
    }
    if ((U32) nBytes != length) {
      char b[512];
      sprintf(b, "byte count: %d %d, errno: %d", nBytes, length, errno);
      (void) COND_PushCondition(DCM_GENERALWARNING,
                                DCM_Message(DCM_GENERALWARNING),
                                "exportData", b);
      return COND_PushCondition(DCM_FILEACCESSERROR,
                                DCM_Message(DCM_FILEACCESSERROR),
                                (*object)->fileName,
                                "exportData");
    }
#ifdef LITTLE_ENDIAN_ARCHITECTURE
    if (item->element.representation == DCM_AT) {
      DCM_ELEMENT e;
      e = *element;
      e.length = length;
      e.d.ot = b;
      swapATGroupElement(&e);
    }
#endif
    if (byteOrder != item->byteOrder) {
      DCM_ELEMENT e;
      e = *element;
      e.length = length;
      e.d.ot = b;
      swapInPlace(object, &e);
    }
    *rtnLength = (U32) nBytes;
    item->currentOffset += nBytes;
  } else {
    p = src;
    switch (element->representation) {
    case DCM_AE:
    case DCM_AS:
    case DCM_CS:
    case DCM_DA:
    case DCM_DT:
    case DCM_DD:
    case DCM_DS:
    case DCM_FD:
    case DCM_IS:
    case DCM_LO:
    case DCM_LT:
    case DCM_OB:
    case DCM_OT:
    case DCM_PN:
    case DCM_SH:
    case DCM_SQ:
    case DCM_ST:
    case DCM_TM:
    case DCM_UI:
    case DCM_UT:
      (void) memcpy(b, p, length);
      *rtnLength = length;
      break;
    case DCM_AT:
      tag = (DCM_TAG *) p;
      while (length >= 4) {
        groupElement.sh[0] = DCM_TAG_GROUP(*tag);
        groupElement.sh[1] = DCM_TAG_ELEMENT(*tag);
        if (byteOrder == BYTEORDER_SAME) {
          *b++ = groupElement.ch[0];  /* Copy the group */
          *b++ = groupElement.ch[1];
          *b++ = groupElement.ch[2];  /* Now, the element */
          *b++ = groupElement.ch[3];
        } else {
          *b++ = groupElement.ch[1];  /* Copy the group */
          *b++ = groupElement.ch[0];
          *b++ = groupElement.ch[3];  /* Now, the element */
          *b++ = groupElement.ch[2];
        }
        tag++;
        length -= 4;
        *rtnLength += 4;
      }
      break;
    case DCM_SL:
    case DCM_UL:
    case DCM_FL:
      while (length >= 4) {
        if (byteOrder == BYTEORDER_SAME) {
          *b++ = *p++;
          *b++ = *p++;
          *b++ = *p++;
          *b++ = *p++;
        } else {
          *b++ = p[3];
          *b++ = p[2];
          *b++ = p[1];
          *b++ = p[0];
          p += 4;
        }
        length -= 4;
        *rtnLength += 4;
      }
      break;
    case DCM_SS:
    case DCM_US:
    case DCM_OW:
      /*
       * Temorary hack by Nilesh to support memory mapping for testing
       * purposes.
       */
      length &= ~1;
      *rtnLength += length;
      if (element->tag == DCM_PXLPIXELDATA) {
        if (byteOrder == item->byteOrder)
          (void) memcpy(b, p, length);
        else
#ifdef SOLARIS
          swab((char *) p, (char *) b, length);
#elif defined AIXV3
          swab((short *) p, (short *) b, length);
#elif defined MACOS
          /* Not Yet Defined */
#else
          swab(p, b, (size_t)length);
#endif
      } else {
        if (byteOrder == BYTEORDER_SAME)
          (void) memcpy(b, p, length);
        else
#ifdef SOLARIS
          swab((char *) p, (char *) b, length);
#elif defined AIXV3
          swab((short *) p, (short *) b, length);
#elif defined MACOS
          /* Not Yet Defined */
#else
          swab(p, b, (size_t)length);
#endif
      }
      break;
      /*case DCM_UNKNOWN:*/
    case DCM_UN:
    default:
#if 0
      fprintf(stderr, "Should not get to default in exportData: %08x\n",
              element->tag);
#endif
      (void) memcpy(b, p, length);
      *rtnLength = length;
      break;
    }
  }
  return DCM_NORMAL;
}

static CONDITION
exportEncapsulatedPixels(PRIVATE_OBJECT ** object, PRV_ELEMENT_ITEM * item,
                         unsigned char* buffer, U32 bufferlength,
                         DCM_EXPORT_STREAM_CALLBACK* callback,
                         void* ctx) {
  DCM_ELEMENT * element;
  ssize_t nBytes;
  CONDITION cond;
  U32 toExport;
  U32 length;
  DCM_FRAGMENT_ITEM* fragmentItem = 0;
  DCM_ELEMENT e;
  U32 rtnLength = 0;

  element = &item->element;
  if (element->d.ot == NULL) {
    if ((*object)->fd != -1) {
      /* Seek to the beginning of the data. Have to back up 12 bytes to
      ** get the pixel tag, VR, etc
      */
      (void) lseek((*object)->fd, item->dataOffset-12, SEEK_SET);
    } else {
      (*object)->sk((*object)->userCtx, item->dataOffset-12, SEEK_SET);
    }

    toExport = item->originalDataLength + 12;
    while (toExport > 0) {
      length = (toExport < bufferlength) ? toExport : bufferlength;

      if ((*object)->fd != -1) {
        nBytes = read((*object)->fd, buffer, length);
      } else {
        int nBytes_int;
        cond = (*object)->rd((*object)->userCtx, buffer,
                             (long) length, &nBytes_int);
	nBytes = (ssize_t)nBytes_int;
      }
      if ((U32) nBytes != length) {
        char b[512];
        sprintf(b, "byte count: %lu %u, errno: %d", nBytes, length, errno);
        (void) COND_PushCondition(DCM_GENERALWARNING,
                                  DCM_Message(DCM_GENERALWARNING),
                                  "exportEncapsualtedPixels", b);
        return COND_PushCondition(DCM_FILEACCESSERROR,
                                  DCM_Message(DCM_FILEACCESSERROR),
                                  (*object)->fileName,
                                  "exportEncapsualtedPixels");
      }
      cond = callback(buffer, length, 0, ctx);
      if (cond != DCM_NORMAL) {
        return COND_PushCondition(DCM_CALLBACKABORTED,
                                  DCM_Message(DCM_CALLBACKABORTED),
                                  "exportStream");
      }
      toExport -= length;
    }
  } else {
    if (item->fragmentFlag != 1) {
      return COND_PushCondition(DCM_NOFRAGMENTSINOBJECT,
                                "DCM Exporting pixels but "
                                "did not find expected fragments in object");
    }
    e.tag = DCM_PXLPIXELDATA;
    e.d.ot = 0;
    e.representation = DCM_OB;
    e.length = 0xffffffff;
    exportFixedFields(&e, buffer, bufferlength,
                      LITTLE_ORDER /*byteOrder*/,
                      1 /* explicitV*/,
                      &rtnLength);
    toExport = rtnLength;
    e.tag = 0xfffee000;
    e.length = 0;
    e.representation = DCM_DLM;
    e.d.ot = 0;
    exportFixedFields(&e, buffer+toExport, bufferlength,
                      LITTLE_ORDER /*byteOrder*/,
                      1 /* explicitV*/,
                      &rtnLength);
    toExport += rtnLength;

    cond = callback(buffer, toExport, 0, ctx);
    if (cond != DCM_NORMAL) {
      return COND_PushCondition(DCM_CALLBACKABORTED,
                                DCM_Message(DCM_CALLBACKABORTED),
                                "exportEncapsulatedPixels");
    }

    fragmentItem = (DCM_FRAGMENT_ITEM*)LST_Head(&item->element.d.fragments);
    (void)LST_Position(&item->element.d.fragments, fragmentItem);
    while (fragmentItem != NULL) {
      printf("Fragment size: %6d\n", fragmentItem->length);
      e.tag = 0xfffee000;
      e.length = fragmentItem->length;
      e.representation = DCM_DLM;
      exportFixedFields(&e, buffer, bufferlength,
                        LITTLE_ORDER /*byteOrder*/,
                        1 /* explicitV*/,
                        &rtnLength);
      cond = callback(buffer, rtnLength, 0, ctx);
      if (cond != DCM_NORMAL) {
        return COND_PushCondition(DCM_CALLBACKABORTED,
                                  DCM_Message(DCM_CALLBACKABORTED),
                                  "exportEncapsulatedPixels");
      }
      cond = callback(fragmentItem->fragment, fragmentItem->length, 0, ctx);
      if (cond != DCM_NORMAL) {
        return COND_PushCondition(DCM_CALLBACKABORTED,
                                  DCM_Message(DCM_CALLBACKABORTED),
                                  "exportEncapsulatedPixels");
      }

      fragmentItem = (DCM_FRAGMENT_ITEM*)LST_Next(&item->element.d.fragments);
    }
    e.tag = 0xfffee0dd;
    e.length = 0;
    e.representation = DCM_DLM;
    e.d.ot = 0;
    exportFixedFields(&e, buffer, bufferlength,
                      LITTLE_ORDER /*byteOrder*/,
                      1 /* explicitV*/,
                      &rtnLength);
    cond = callback(buffer, rtnLength, 0, ctx);
    if (cond != DCM_NORMAL) {
      return COND_PushCondition(DCM_CALLBACKABORTED,
                                DCM_Message(DCM_CALLBACKABORTED),
                                "exportEncapsulatedPixels");
    }
  }
  return DCM_NORMAL;
}

static CONDITION
exportPixels(PRIVATE_OBJECT ** object, PRV_ELEMENT_ITEM * item,
             int encapsulatedPixels,
             unsigned char* buffer, U32 bufferlength,
             DCM_EXPORT_STREAM_CALLBACK* callback,
             void* ctx,
             int byteOrder, int explicitVR) {
  DCM_ELEMENT * element;
  CONDITION cond;
  U32 bytesExported = 0;
  U32 exportLength = 0;
  U32 rtnLength;
  U32 remainingData;
  unsigned char* dst;
  unsigned char* src;
  U32 c;

  if (encapsulatedPixels) {
    return exportEncapsulatedPixels(object, item, buffer,
                                    bufferlength, callback, ctx);
  }

  element = &item->element;
  rtnLength = 0;
  dst = buffer;
  c = bufferlength;
  exportFixedFields(element, dst, bufferlength, byteOrder,
                    explicitVR, &rtnLength);
  dst += rtnLength;
  c -= rtnLength;
  bytesExported = rtnLength;

  remainingData = element->length;
  src = (unsigned char*)element->d.ot;
  item->currentOffset = item->dataOffset;

  while (remainingData > 0) {
    if (debug) {
      fprintf(stderr, "Export: (%08x) %d\n",
              (unsigned int) element->tag, element->length);
    }

    if (element->d.ot != NULL) {
      remainingData = element->length -
                      (src - ((unsigned char *) element->d.ot));
    } else {
      remainingData = element->length -
                      (item->currentOffset - item->dataOffset);
    }

    exportLength = (remainingData < c) ? remainingData : c;
    cond = exportData(object, item, src, dst,
                      exportLength, byteOrder, &rtnLength);
    if (cond != DCM_NORMAL)
      return cond;

    src += rtnLength;
    dst += rtnLength;
    bytesExported += rtnLength;
    c -= rtnLength;

    if (c <= 20) {
      cond = callback(buffer, bytesExported, 0, ctx);
      if (cond != DCM_NORMAL) {
        return COND_PushCondition(DCM_CALLBACKABORTED,
                                  DCM_Message(DCM_CALLBACKABORTED),
                                  "exportPixels");
      }
      bytesExported = 0;
      c = bufferlength;
      dst = (unsigned char *) buffer;
    }
  }
  if (bytesExported > 0) {
    cond = callback(buffer, bytesExported, 0, ctx);
    if (cond != DCM_NORMAL) {
      return COND_PushCondition(DCM_CALLBACKABORTED,
                                DCM_Message(DCM_CALLBACKABORTED),
                                "exportPixels");
    }
  }

  return DCM_NORMAL;

#if 0
  if (element->d.ot == NULL) {
    if ((*object)->fd != -1) {
      /* Seek to the beginning of the data. Have to back up 12 bytes to
      ** get the pixel tag, VR, etc
      */
      (void) lseek((*object)->fd, item->dataOffset-12, SEEK_SET);
    } else {
      (*object)->sk((*object)->userCtx, item->dataOffset-12, SEEK_SET);
    }

    toExport = item->originalDataLength + 12;
    while (toExport > 0) {
      length = (toExport < bufferlength) ? toExport : bufferlength;

      if ((*object)->fd != -1) {
        nBytes = read((*object)->fd, buffer, length);
      } else {
        cond = (*object)->rd((*object)->userCtx, buffer,
                             (long) length, &nBytes);
      }
      if ((U32) nBytes != length) {
        char b[512];
        sprintf(b, "byte count: %d %d, errno: %d", nBytes, length, errno);
        (void) COND_PushCondition(DCM_GENERALWARNING,
                                  DCM_Message(DCM_GENERALWARNING),
                                  "exportPixels", b);
        return COND_PushCondition(DCM_FILEACCESSERROR,
                                  DCM_Message(DCM_FILEACCESSERROR),
                                  (*object)->fileName,
                                  "exportPixels");
      }
      cond = callback(buffer, length, 0, ctx);
      if (cond != DCM_NORMAL) {
        return COND_PushCondition(DCM_CALLBACKABORTED,
                                  DCM_Message(DCM_CALLBACKABORTED),
                                  "exportStream");
      }
      toExport -= length;
    }
  } else {}
  return DCM_NORMAL;
#endif

}

/* fileSize
**
** Purpose:
**  Determine the file size of a file on an open descriptor.
**
** Parameter Dictionary:
**  fd  File descriptor for an open file.
**
** Return Values:
**  the size of the open file in bytes (nonnegative)
**  a negative status value returned by fstat
**
** Algorithm:
**  call unix fstat system call to get file size.
**  if successful call
**      return file size
**  else
**      return status returned by fstat call (-1)
*/
#ifdef MACOS
static long
#else
static int
#endif
fileSize(int fd) {
  int
  status;
  struct stat
        im_stat;

  status = fstat(fd, &im_stat);
  if (status < 0) {
    return status;
  } else
    return im_stat.st_size;
}

/* swapInPlace
**
** Purpose:
**  Swap data in place for byte order adjustment.  Bytes are swapped
**  for data with representations of DCM_US and DCM_UL (binary values).
**  Special care is taken with pixel data which may be 8 bits.
**
** Parameter Dictionary:
**  object    Pointer to caller's DCM object containing the
**      element with the data to be swapped
**  element   Pointer to DCM_ELEMENT that contains the data to be
**      swapped.
**
** Return Values:
**  None
**
** Algorithm:
**  If (element->representation is 16 bit binary)
**      If (element is pixel data and pixel data is not 16 bits)
**    return
**      Swap in place short integers for this element.
**  Else if (element->representation is 32 bit binary)
**      Swap in place long integers for this element
*/

static void
swapInPlace(PRIVATE_OBJECT ** object, DCM_ELEMENT * e) {
  U32
  length;
  unsigned char
  tmp,
  *p1;

  length = e->length;
  p1 = (unsigned char*)e->d.ot;
  if (e->representation == DCM_US || e->representation == DCM_SS ||
      e->representation == DCM_OW || e->representation == DCM_AT) {
    if (e->tag == DCM_PXLPIXELDATA &&
        (*object)->pixelBitsAllocated != 16)
      return;

    while (length > 0) {
      tmp = p1[0];
      p1[0] = p1[1];
      p1[1] = tmp;
      p1 += 2;
      length -= 2;
    }
  } else if (e->representation == DCM_UL || e->representation == DCM_SL) {

    // there are cases where the representaiton is UL but length were
    // not divisible by 4!   In this way, at least it reads the length
    // specified

    while (length > 0) {
      tmp = p1[0];
      p1[0] = p1[1];
      p1[1] = tmp;
      p1 += 2;
      length -= 2;
    }

#if 0
    while (length > 0) {
      tmp = p1[0];
      p1[0] = p1[3];
      p1[3] = tmp;
      tmp = p1[1];
      p1[1] = p1[2];
      p1[2] = tmp;
      length -= 4;
      p1 += 4;
    }
#endif
  }
}


/* checkObject
**
** Purpose:
**  Examine a PRIVATE OBJECT to see if it looks like is has the proper
**  fields defined.  This function is used to try to make certain that
**  users call the DCM routines with the proper objects.  If the object
**  is legal, the function returns DCM_NORMAL.  If the object is not
**  legal, the function will return an error.
**
** Parameter Dictionary:
**  object    PRIVATE_OBJECT to be examined by this function
**  caller    Name of the function (ASCIZ) that called this
**      function.  In case of failure, this becomes part of
**      the error message that is pushed on the stack.
**
** Return Values:
**  DCM_NORMAL
**  DCM_NULLOBJECT
**  DCM_ILLEGALOBJECT
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

static CONDITION
checkObject(PRIVATE_OBJECT ** object, const char *caller) {
  if (object == NULL)
    return COND_PushCondition(DCM_NULLOBJECT, DCM_Message(DCM_NULLOBJECT),
                              caller);
  if (*object == NULL)
    return COND_PushCondition(DCM_NULLOBJECT, DCM_Message(DCM_NULLOBJECT),
                              caller);

  if (strcmp((*object)->keyType, KEY_DCM_OBJECT) != 0)
    return COND_PushCondition(DCM_ILLEGALOBJECT,
                              DCM_Message(DCM_ILLEGALOBJECT), caller);
  return DCM_NORMAL;
}


/* writeFile
**
** Purpose:
**  Write the data in the buffer into the file specified by the file
**  descriptor
**
** Parameter Dictionary:
**  buffer    Buffer holding the information to be written
**  length    Length of the buffer
**  flag    Unused
**  fd    File descriptor
**
** Return Values:
**
**  DCM_FILEIOERROR
**  DCM_NORMAL
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

static CONDITION
writeFile(void *buffer, U32 length, int flag,
          void /* int */ *fdPtr) {
  int
  bytesWritten;
  int *fd;

  fd = (int *) fdPtr;

  bytesWritten = write(*fd, buffer, (int) length);
  if (bytesWritten != (int) length)
    return COND_PushCondition(DCM_FILEIOERROR,
                              DCM_Message(DCM_FILEIOERROR), "",
                              strerror(errno),
                              "writeFile");
  else
    return DCM_NORMAL;
}

static CONDITION
countBytes(void *buffer, U32 length, int flag,
           void /* unsigned long */ *sizePtr) {
  unsigned long *size;

  size = (unsigned long *) sizePtr;

  *size += length;

  return DCM_NORMAL;
}

static CONDITION
setFileOptions(DCM_OBJECT ** obj, unsigned long *opt) {
  CONDITION cond;
  char xferSyntax[DICOM_UI_LENGTH + 1];
  DCM_ELEMENT e = {DCM_METATRANSFERSYNTAX, DCM_UI, "", 1, sizeof(xferSyntax),
                   {NULL}};

  e.d.string = xferSyntax;
  cond = DCM_ParseObject(obj, &e, 1, NULL, 0, NULL);
  if (cond != DCM_NORMAL)
    return cond;

  *opt = 0;
  if (strcmp(xferSyntax, DICOM_TRANSFERLITTLEENDIAN) == 0) {
    *opt = DCM_ORDERLITTLEENDIAN;
  } else if (strcmp(xferSyntax, DICOM_TRANSFERLITTLEENDIANEXPLICIT) == 0) {
    *opt = DCM_EXPLICITLITTLEENDIAN;
  } else if (strcmp(xferSyntax, DICOM_TRANSFERBIGENDIANEXPLICIT) == 0) {
    *opt = DCM_EXPLICITBIGENDIAN;
  } else {  /* Must be an encapsulated xfer syntax */
    *opt = DCM_ENCAPSULATEDPIXELS;
  }

  return DCM_NORMAL;
}

static CONDITION
extractFileOptions(unsigned long opt, CTNBOOLEAN * part10File,
                   CTNBOOLEAN * explicitVR, int *byteOrder,
                   CTNBOOLEAN* encapsulatedPixels) {
  *part10File = *explicitVR = FALSE;

  if ((opt & DCM_FILEFORMATMASK) == DCM_PART10FILE) {
    *part10File = TRUE;
    opt &= ~DCM_ORDERMASK;
    opt |= DCM_EXPLICITLITTLEENDIAN;
  }
  if ((opt & DCM_ORDERMASK) == 0)
    return COND_PushCondition(DCM_ILLEGALOPTION,
                              DCM_Message(DCM_ILLEGALOPTION), "Byte order",
                              "extractFileOptions");

  switch (opt & DCM_ORDERMASK) {
  case DCM_ORDERNATIVE:
    *byteOrder = NATIVE_ORDER;
    *encapsulatedPixels = FALSE;
    break;
  case DCM_ORDERLITTLEENDIAN:
    *byteOrder = LITTLE_ORDER;
    *encapsulatedPixels = FALSE;
    break;
  case DCM_EXPLICITLITTLEENDIAN:
    *byteOrder = LITTLE_ORDER;
    *explicitVR = TRUE;
    *encapsulatedPixels = FALSE;
    break;
  case DCM_ORDERBIGENDIAN:
    *byteOrder = BIG_ORDER;
    *encapsulatedPixels = FALSE;
    break;
  case DCM_EXPLICITBIGENDIAN:
    *byteOrder = BIG_ORDER;
    *explicitVR = TRUE;
    *encapsulatedPixels = FALSE;
    break;
  case DCM_ENCAPSULATEDPIXELS:
    *byteOrder = LITTLE_ORDER;
    *explicitVR = TRUE;
    *encapsulatedPixels = TRUE;
    break;
  default:
    *byteOrder = LITTLE_ORDER;
    *encapsulatedPixels = FALSE;
    break;
  }

  return DCM_NORMAL;
}

static U32
computeGroupLength(PRV_GROUP_ITEM * groupItem,
                   CTNBOOLEAN explicitVR) {
  return (explicitVR) ?
         groupItem->baseLength + 4 * groupItem->longVRAttributes :
         groupItem->baseLength;

}

/* exportStream
**
** Purpose:
**      Export a DICOM object into the stream format suitable
**      for network transmission or disk storage.
**
** Parameter Dictionary:
**  callerObject    Handle to caller's DICOM object
**  opt     Bit mask giving options for exporting data
**  buffer      Pointer to caller's buffer to hold next chunk
**        of DCM stream data
**  bufferlength    Length of caller's buffer to hold stream data
**  callback    Callback routine to be called.
**  ctx     Pointer to context variable we maintain to keep
**        track of our location in export process.
**  sequenceLevel   Current level in the sequence hierarchy
**
** Return Values:
**
**  DCM_FILEACCESSERROR
**  DCM_ILLEGALOBJECT
**  DCM_LISTFAILURE
**  DCM_NORMAL
**  DCM_NULLOBJECT
**  DCM_CALLBACKABORTED
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/

static CONDITION
exportStream(DCM_OBJECT ** callerObject, unsigned long opt,
             void *buffer, U32 bufferlength,
             DCM_EXPORT_STREAM_CALLBACK* callback,
             void *ctx, int sequenceLevel) {
  PRIVATE_OBJECT
  ** object;
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem;
  DCM_ELEMENT
  element;
  int
  byteOrder;
  int
  lastFlag = 0;
  unsigned char
  *src,
  *dst;
  U32
  c,
  bytesExported = 0,
                  rtnLength,
                  remainingData,
                  exportLength;
  CONDITION
  cond;
  DCM_SEQUENCE_ITEM
  * sequenceItem;
  DCM_ELEMENT
  itemMarker = {
                 DCM_DLMITEM, DCM_DLM, "", 1, DCM_UNSPECIFIEDLENGTH, {NULL}
               },
               itemDelimiter = {
                                 DCM_DLMITEMDELIMITATIONITEM, DCM_DLM, "", 1, 0, {NULL}
                               },
                               sequenceDelimiter = {
                                                     DCM_DLMSEQUENCEDELIMITATIONITEM, DCM_DLM, "", 1, 0, {NULL}
                                                   };
  CTNBOOLEAN
  unspecifiedSQLength = FALSE,
                        explicitVR = FALSE,
                                     part10File = FALSE,
                                                  encapsulatedPixels = FALSE;
  unsigned long fileOptions = 0;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "exportStream");
  if (cond != DCM_NORMAL)
    return cond;

  if ((opt & DCM_FILEFORMATMASK) == DCM_PART10FILE) {
    part10File = TRUE;
    opt &= ~DCM_ORDERMASK;
    opt |= DCM_EXPLICITLITTLEENDIAN;
    cond = setFileOptions(callerObject, &fileOptions);
    if (cond != DCM_NORMAL)
      return cond;
  }
  if ((opt & DCM_ORDERMASK) == 0)
    return COND_PushCondition(DCM_ILLEGALOPTION,
                              DCM_Message(DCM_ILLEGALOPTION),
                              "Byte order",
                              "exportStream");

  switch (opt & DCM_ORDERMASK) {
  case DCM_ORDERNATIVE:
    byteOrder = NATIVE_ORDER;
    break;
  case DCM_ORDERLITTLEENDIAN:
    byteOrder = LITTLE_ORDER;
    break;
  case DCM_EXPLICITLITTLEENDIAN:
    byteOrder = LITTLE_ORDER;
    explicitVR = TRUE;
    break;
  case DCM_ORDERBIGENDIAN:
    byteOrder = BIG_ORDER;
    break;
  case DCM_EXPLICITBIGENDIAN:
    byteOrder = BIG_ORDER;
    explicitVR = TRUE;
    break;
  case DCM_ENCAPSULATEDPIXELS:
    byteOrder = LITTLE_ORDER;
    explicitVR = TRUE;
    encapsulatedPixels = TRUE;
    break;
  default:
    byteOrder = LITTLE_ORDER;
    break;
  }

  /*  We are making this step mandatory for now (smm)*/

  opt &= ~DCM_SEQUENCELENGTHMASK;
  opt |= DCM_UNSPECIFIEDLENGTHFLAG;

  /*  End of something that is out of place */

  if ((opt & DCM_SEQUENCELENGTHMASK) == DCM_UNSPECIFIEDLENGTHFLAG)
    unspecifiedSQLength = TRUE;

  dst = (unsigned char *) buffer;
  c = bufferlength;

  if (part10File) {
    cond = exportPreamble(object, dst, c, &rtnLength);
    if (cond != DCM_NORMAL)
      return cond;

    dst += rtnLength;
    c -= rtnLength;
    bytesExported += rtnLength;
  }
  if (sequenceLevel != 0) {
    if (!unspecifiedSQLength)
      itemMarker.length = (*object)->objectSize;
    exportFixedFields(&itemMarker, dst, bufferlength, byteOrder,
                      explicitVR, &rtnLength);
    dst += rtnLength;
    c -= rtnLength;
    bytesExported += rtnLength;
  }
  groupItem = (PRV_GROUP_ITEM*)LST_Head(&(*object)->groupList);

  /*  Take this code out to allow empty groups. */
#if 0
  if (groupItem == NULL)
    return COND_PushCondition(DCM_LISTFAILURE,
                              DCM_Message(DCM_LISTFAILURE),
                              "exportStream");
#endif
  if (groupItem != NULL)
    (void) LST_Position(&(*object)->groupList, groupItem);

  while (groupItem != NULL) {
    if (part10File && groupItem->group != DCM_GROUPFILEMETA) {
      if (opt != fileOptions) {
        opt = fileOptions;
        cond = extractFileOptions(opt, &part10File,
                                  &explicitVR, &byteOrder,
                                  &encapsulatedPixels);
        if (cond != DCM_NORMAL)
          return cond;
      }
    }
    elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
    if (elementItem != NULL) {
      (void) LST_Position(&groupItem->elementList, elementItem);
      if (DCM_TAG_ELEMENT(elementItem->element.tag) == 0x0000) {
        U32 l;
        l = computeGroupLength(groupItem, explicitVR);
        *elementItem->element.d.ul = l;

        /* We have some problems computing group length
           for groups with sequences.
           ** For now, just remove this attribute, except
           for group 0000 and 0002.
        */
        if (groupItem->group != 0x0000 && groupItem->group != 0x0002)
          elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
      }
    }
    while (elementItem != NULL) {
      if (c <= 20) {
        cond = callback(buffer, bytesExported, 0, ctx);
        if (cond != DCM_NORMAL)
          return COND_PushCondition(DCM_CALLBACKABORTED,
                                    DCM_Message(DCM_CALLBACKABORTED),
                                    "exportStream");

        bytesExported = 0;
        c = bufferlength;
        dst = (unsigned char *) buffer;
      }
      element = elementItem->element;

      if (element.tag == DCM_PXLPIXELDATA) {
        /* Lots of special rules for pixel data. Handle separately */
        /* First, dump the current buffer */
        cond = callback(buffer, bytesExported, 0, ctx);
        if (cond != DCM_NORMAL)
          return COND_PushCondition(DCM_CALLBACKABORTED,
                                    DCM_Message(DCM_CALLBACKABORTED),
                                    "exportStream");

        cond = exportPixels(object, elementItem, encapsulatedPixels,
                            (unsigned char*)buffer, bufferlength, callback,
                            ctx, byteOrder, explicitVR);
        if (cond != DCM_NORMAL)
          return cond;

        bytesExported = 0;
        c = bufferlength;
        dst = (unsigned char *) buffer;
        rtnLength = 0;
      } else if (element.representation == DCM_SQ) {
        if (unspecifiedSQLength)
          element.length = DCM_UNSPECIFIEDLENGTH;

        exportFixedFields(&element, dst, bufferlength, byteOrder,
                          explicitVR, &rtnLength);
      } else {
        element.length = elementItem->paddedDataLength;
        exportFixedFields(&element, dst, bufferlength, byteOrder,
                          explicitVR, &rtnLength);
      }
      dst += rtnLength;
      c -= rtnLength;
      bytesExported += rtnLength;

      remainingData = element.length;
      src = (unsigned char*)element.d.ot;
      elementItem->currentOffset = elementItem->dataOffset;

      if (element.tag == DCM_PXLPIXELDATA) {
        /* Then, we did that above */
        ;
      } else if (element.representation == DCM_SQ) {

        cond = callback(buffer, bytesExported, 0, ctx);
        if (cond != DCM_NORMAL)
          return COND_PushCondition(DCM_CALLBACKABORTED,
                                    DCM_Message(DCM_CALLBACKABORTED),
                                    "exportStream");

        bytesExported = 0;
        c = bufferlength;
        dst = (unsigned char *) buffer;

        if (element.d.sq != NULL) {
          sequenceItem = (DCM_SEQUENCE_ITEM*)LST_Head(&element.d.sq);
          if (sequenceItem != NULL)
            (void) LST_Position(&element.d.sq, sequenceItem);
          while (sequenceItem != NULL) {
            cond = exportStream(&sequenceItem->object, opt,
                                buffer, bufferlength, callback, ctx,
                                sequenceLevel + 1);
            if (cond != DCM_NORMAL)
              return cond;
            sequenceItem = (DCM_SEQUENCE_ITEM*)LST_Next(&element.d.sq);
          }
        }
        if (element.length == DCM_UNSPECIFIEDLENGTH) {
          sequenceDelimiter.length = 0;
          exportFixedFields(&sequenceDelimiter, dst, bufferlength,
                            byteOrder, explicitVR, &rtnLength);
          dst += rtnLength;
          c -= rtnLength;
          bytesExported += rtnLength;
        }
      } else {
        while (remainingData > 0) {
          if (debug)
            fprintf(stderr, "Export: (%08x) %d\n",
                    (unsigned int) element.tag, element.length);
          if (element.d.ot != NULL)
            remainingData = element.length -
                            (src - ((unsigned char *) element.d.ot));
          else
            remainingData = element.length -
                            (elementItem->currentOffset - elementItem->dataOffset);

          exportLength = (remainingData < c) ? remainingData : c;
          cond = exportData(object, elementItem, src, dst,
                            exportLength, byteOrder, &rtnLength);
          if (cond != DCM_NORMAL)
            return cond;

          src += rtnLength;
          dst += rtnLength;
          bytesExported += rtnLength;
          c -= rtnLength;

          if (c <= 20) {
            cond = callback(buffer, bytesExported, 0, ctx);
            if (cond != DCM_NORMAL)
              return COND_PushCondition(DCM_CALLBACKABORTED,
                                        DCM_Message(DCM_CALLBACKABORTED),
                                        "exportStream");

            bytesExported = 0;
            c = bufferlength;
            dst = (unsigned char *) buffer;
          }
        }
      }
      elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
    }
    groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
  }
  if ((sequenceLevel != 0) && unspecifiedSQLength) {
    if (c <= 20) {
      cond = callback(buffer, bytesExported, 0, ctx);
      if (cond != DCM_NORMAL)
        return COND_PushCondition(DCM_CALLBACKABORTED,
                                  DCM_Message(DCM_CALLBACKABORTED),
                                  "exportStream");

      bytesExported = 0;
      c = bufferlength;
      dst = (unsigned char *) buffer;
    }
    exportFixedFields(&itemDelimiter, dst, bufferlength, byteOrder,
                      explicitVR, &rtnLength);
    dst += rtnLength;
    c -= rtnLength;
    bytesExported += rtnLength;
  }
  lastFlag = (sequenceLevel == 0) ? 1 : 0;
  cond = callback(buffer, bytesExported, lastFlag, ctx);
  if (cond != DCM_NORMAL)
    return COND_PushCondition(DCM_CALLBACKABORTED,
                              DCM_Message(DCM_CALLBACKABORTED),
                              "exportStream");

  return DCM_NORMAL;
}

/* verifyFormat
**
** Purpose:
**  This routine verifies the format of the
**  data value of attributes according to
**  the DICOM v3 standard.
**
** Parameter Dictionary:
**   element    Pointer to the DCM_ELEMENT containing the
**              element to be examined.
**
** Return Values:
**  DCM_NORMAL
**
** Algorithm:
**  switch(representation) {
**  case DCM_DA:
**      Retain all characters that are digits, '-', or '\'
**      If the resulting string is of odd length
**    Pad the string with ' '
**      break;
**  case DCM_TM:
**      Retain all characters that are digits, '.', '-', or '\'
**      If the resulting string is of odd length
**    Pad the string with ' '
**      break;
**  case DCM_CS, DCM_AS, DCM_DS, DCM_IS, DCM_LO, DCM_SH, DCM_UT:
**      Delete all the leading and trailing spaces.
**      If the resulting string is of odd length
**    Pad the string with ' '
**      break;
**  case DCM_LT, DCM_ST, DCM_PN:
**      Delete all the trailing spaces.
**      If the resulting string is of odd length
**    Pad the string with ' '
**      break;
**  }
*/

static CONDITION
verifyFormat(PRV_ELEMENT_ITEM * item) {
  int
  i,
  l;
  char
  *src,
  *dst,
  *p;
  DCM_ELEMENT
  * element;
  CTNBOOLEAN
  stopFlag = FALSE;

  element = &item->element;
  if (element->length > 0) {
    switch (element->representation) {
    case DCM_DA:
      src = dst = element->d.string;
      l = (int) element->length;
      for (i = 0; i < l; i++) {
        if (isdigit(*src) || (*src == '-') || (*src == '\\')) {
          *dst++ = *src++;
        } else {
          src++;
          element->length--;
        }
      }
      item->paddedDataLength = element->length;
      if (element->length & 1) {
        *dst = ' ';
        item->paddedDataLength++;
      }
      break;
    case DCM_TM:
      l = (int) element->length;
      src = dst = element->d.string;
      for (i = 0; i < l; i++) {
        if (isdigit(*src) || \
            (*src == '.') || \
            (*src == '-') || \
            (*src == '\\')) {
          *dst++ = *src++;
        } else {
          src++;
          element->length--;
        }
      }
      item->paddedDataLength = element->length;
      if (element->length & 1) {
        *dst = ' ';
        item->paddedDataLength++;
      }
      break;
      /*
       * Both the leading and trailing spaces are non-significant.
       */
    case DCM_CS:
    case DCM_AS:
    case DCM_DS:
    case DCM_IS:
    case DCM_LO:
    case DCM_SH:
    case DCM_UT:
      l = (int) element->length;
      src = dst = element->d.string;
      for (i = 0; i < l; i++) {
        if ((*src == ' ') && !stopFlag) {
          src++;
          element->length--;
        } else {
          stopFlag = TRUE;
          *dst++ = *src++;
        }
      }
      /*
       * Right now, dst points to the char follows the last char in the
       * string.
       */
      stopFlag = FALSE;
      l = (int) element->length;
      p = dst - 1;
      for (i = l; (i > 0) && !stopFlag; i--) {
        if ((*p == ' ') && !stopFlag) {
          p--;
          dst--;
          element->length--;
        } else
          stopFlag = TRUE;
      }
      item->paddedDataLength = element->length;
      if (element->length & 1) {
        *dst = ' ';
        item->paddedDataLength++;
      }
      break;
      /*
       * The trailing spaces are non-significant.
       */
    case DCM_LT:
    case DCM_ST:
      l = (int) element->length;
      src = element->d.string + l - 1;
      for (i = l; (i > 0) && !stopFlag; i--) {
        if ((*src == ' ') && !stopFlag) {
          src--;
          element->length--;
        } else
          stopFlag = TRUE;
      }
      item->paddedDataLength = element->length;
      if (element->length & 1) {
        *++src = ' ';
        item->paddedDataLength++;
      }
      break;
    case DCM_PN:
      /*
       * Strip off the trailing spaces.
       */
      l = (int) element->length;
      src = element->d.string + l - 1;
      for (i = l; (i > 0) && !stopFlag; i--) {
        if ((*src == ' ') && !stopFlag) {
          src--;
          element->length--;
        } else
          stopFlag = TRUE;
      }
      /*
       * Convert the name to the standard V3 format.
       */
      src = dst = element->d.string;
      l = element->length;
      for (i = 0; i < l;) {
        if ((src[i] == ',') || (src[i] == '^')) {
          *dst++ = '^';
          i++;
          while ((i < l) && (src[i] == ' ')) {
            element->length--;
            i++;
          }
        } else {
          *dst++ = src[i++];
        }
      }

      item->paddedDataLength = element->length;
      if (element->length & 1) {
        *dst = ' ';
        item->paddedDataLength++;
      }
      break;
    case DCM_UI:
      if (element->d.string[element->length - 1] == '\0')
        element->length--;
      if (element->d.string[element->length - 1] == ' ') {
        element->d.string[element->length - 1] = '\0';
        element->length--;
      }
      break;
    default:
      break;
    }
  }
  return DCM_NORMAL;
}

/* readFile
**
** Purpose:
**  Read DICOM object from a file
**
** Parameter Dictionary:
**  name      Name of the file
**  callerBuf   Buffer from which to read the object
**  fd      File descriptor
**  size      Size of the file
**  fileOffset    Offset in the file from which point read starts
**  recursionLevel    Level of recursion
**  opt     Indicates in what byte order to read
**  callerObject    The object into which the contents are stored
**  scannedLength   Length of data scanned
**  remainOpenFlag    Indicates whether the file remains open
**
** Return Values:
**
**  DCM_ELEMENTCREATEFAILED
**  DCM_ELEMENTLENGTHERROR
**  DCM_ELEMENTOUTOFORDER
**  DCM_FILEACCESSERROR
**  DCM_ILLEGALSTREAMLENGTH
**  DCM_LISTFAILURE
**  DCM_NORMAL
**  DCM_OBJECTCREATEFAILED
**  DCM_UNEVENELEMENTLENGTH
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/
#if 0
static CONDITION
readFile(char *name, unsigned char *callerBuf, int fd, long size,
         off_t fileOffset, int recursionLevel,
         unsigned long opt, DCM_OBJECT ** callerObject,
         U32 * scannedLength, CTNBOOLEAN * remainOpenFlag,
         void *ctx,
         CONDITION(*rd) (void *ctx, void *buf, int toRead, int *bytesRead),
         CONDITION(*sk) (void *ctx, int offset, int flag)) {
  CONDITION
  cond;
  int
  byteOrder;
  long
  lastGroup = -1,
              lastElement = -1;
  U32
  sequenceLength,
  localLength;
  PRIVATE_OBJECT
  ** object;
  PRV_GROUP_ITEM
  * groupItem = NULL;
  unsigned short
  group,
  element,
  tagGroup,
  tagElement;
  DCM_ELEMENT
  e,
  tagE;
  CTNBOOLEAN
  pixelFlag,
  convertFlag = FALSE,
                done = FALSE,
                       knownLength = TRUE,
                                     sequenceDone = FALSE,
                                                    createGroupFlag,
                                                    explicitVR = FALSE;
  unsigned char
  buf[8],
  *ptr;
  int
  nBytes;
  PRV_ELEMENT_ITEM
  * elementItem = NULL;
  DCM_OBJECT
  * sequenceObject;
  DCM_SEQUENCE_ITEM
  * sequenceItem;
  CTNBOOLEAN
  fileFlag = TRUE;

  if (callerBuf != NULL) {
    ptr = callerBuf;
    fileFlag = FALSE;
  } else
    ptr = buf;

  switch (opt & DCM_ORDERMASK) {
  case DCM_ORDERNATIVE:
    byteOrder = NATIVE_ORDER;
    break;
  case DCM_ORDERLITTLEENDIAN:
    byteOrder = LITTLE_ORDER;
    break;
  case DCM_EXPLICITLITTLEENDIAN:
    byteOrder = LITTLE_ORDER;
    explicitVR = TRUE;
    break;
  case DCM_ORDERBIGENDIAN:
    byteOrder = BIG_ORDER;
    break;
  case DCM_EXPLICITBIGENDIAN:
    byteOrder = BIG_ORDER;
    explicitVR = TRUE;
    break;
  default:
    byteOrder = NATIVE_ORDER;
    break;
  }
  if ((opt & DCM_CONVERTMASK) == DCM_FORMATCONVERSION)
    convertFlag = TRUE;

  if (scannedLength != NULL)
    *scannedLength = 0;

  cond = DCM_CreateObject(callerObject, opt);
  if (cond != DCM_NORMAL)
    return cond;

  object = (PRIVATE_OBJECT **) callerObject;
  if (fileFlag)
    strcpy((*object)->fileName, name);

  (*object)->fd = -1;
  (*object)->rd = rd;
  (*object)->sk = sk;
  (*object)->userCtx = ctx;
  if (size == (long) DCM_UNSPECIFIEDLENGTH)
    knownLength = FALSE;

  if ((fileFlag) && \
      ((opt & DCM_DELETEMASK) == DCM_DELETEONCLOSE) && \
      (recursionLevel == 0))
    (*object)->deleteFlag = TRUE;

  if (knownLength && (size == 0))
    done = TRUE;

  while (!done) {

    if ((size < 8) && knownLength) {
      if (debug)
        (void) DCM_DumpElements(callerObject, 0);
      (void) DCM_CloseObject(callerObject);
      return COND_PushCondition(DCM_ILLEGALSTREAMLENGTH,
                                DCM_Message(DCM_ILLEGALSTREAMLENGTH), size,
                                "readFile");
    }
    if (fileFlag) {
      if (fd != -1) {
        nBytes = read(fd, buf, 4);
      } else {
        cond = (*object)->rd((*object)->userCtx, buf, 4, &nBytes);
      }

      if (nBytes != 4)
        return COND_PushCondition(DCM_FILEACCESSERROR,
                                  DCM_Message(DCM_FILEACCESSERROR), name,
                                  "readFile");
      ptr = buf;
    }
    if (knownLength)
      size -= 4;
    fileOffset += (off_t) 4;
    if (scannedLength != NULL)
      (*scannedLength) += 4;
    (*object)->objectSize += 4;

    if (byteOrder == BYTEORDER_SAME) {
      GET_SHORT_SAME_ORDER(ptr, group);
      GET_SHORT_SAME_ORDER(ptr + 2, element);
      e.tag = DCM_MAKETAG(group, element);
    } else {
      GET_SHORT_REVERSE_ORDER(ptr, group);
      GET_SHORT_REVERSE_ORDER(ptr + 2, element);
      e.tag = DCM_MAKETAG(group, element);
    }
    ptr += 4;

    if (explicitVR) {
      if (fileFlag) {
        if (fd != -1) {
          nBytes = read(fd, buf, 4);
        } else {
          cond = (*object)->rd((*object)->userCtx, buf, 4, &nBytes);
        }

        if (nBytes != 4)
          return COND_PushCondition(DCM_FILEACCESSERROR,
                                    DCM_Message(DCM_FILEACCESSERROR), name,
                                    "readFile");
        ptr = buf;
      }
      if (knownLength)
        size -= 4;
      fileOffset += (off_t) 4;
      if (scannedLength != NULL)
        (*scannedLength) += 4;
      (*object)->objectSize += 4;
      if ((strncmp((char *) ptr, "OB", 2) == 0) ||
          (strncmp((char *) ptr, "OW", 2) == 0) ||
        (strncmp((char *) ptr, "SQ", 2) == 0)) {}
      else {}
    }
    else {

      if (fileFlag) {
        if (fd != -1) {
          nBytes = read(fd, buf, 4);
        } else {
          cond = (*object)->rd((*object)->userCtx, buf, 4, &nBytes);
        }

        if (nBytes != 4)
          return COND_PushCondition(DCM_FILEACCESSERROR,
                                    DCM_Message(DCM_FILEACCESSERROR), name,
                                    "readFile");
        ptr = buf;
      }
      if (knownLength)
        size -= 4;
      fileOffset += (off_t) 4;
      if (scannedLength != NULL)
        (*scannedLength) += 4;
      (*object)->objectSize += 4;


      if (byteOrder == BYTEORDER_SAME) {
        GET_LONG_SAME_ORDER(ptr, e.length);
      } else {
        GET_LONG_REVERSE_ORDER(ptr, e.length);
      }
      ptr += 4;
    }

    if (((e.length & 1) != 0) && (e.length != DCM_UNSPECIFIEDLENGTH)) {
      if (debug)
        (void) DCM_DumpElements(callerObject, 0);
      (void) DCM_CloseObject(callerObject);
      return COND_PushCondition(DCM_UNEVENELEMENTLENGTH,
                                DCM_Message(DCM_UNEVENELEMENTLENGTH),
                                group, element, e.length,
                                "readFile");
    }
    if ((e.length != (U32) DCM_UNSPECIFIEDLENGTH) && (e.length > (U32) size)) {
      if (debug)
        (void) DCM_DumpElements(callerObject, 0);
      (void) DCM_CloseObject(callerObject);
      return COND_PushCondition(DCM_ELEMENTLENGTHERROR,
                                DCM_Message(DCM_ELEMENTLENGTHERROR),
                                group, element, e.length, size, "readFile");
    }
    if ((e.tag == DCM_DLMITEMDELIMITATIONITEM) ||
        (e.tag == DCM_DLMSEQUENCEDELIMITATIONITEM)) {
      return DCM_NORMAL;
    }
    cond = DCM_LookupElement(&e);
    if (cond != DCM_NORMAL)
      (void) COND_PopCondition(0);
    if (e.representation == DCM_CTX)
      ctxSensitiveLookup(object, &e);

    if (e.representation == DCM_SQ) {
      cond = newElementItem(&e, FALSE, &elementItem);
      if (cond != DCM_NORMAL)
        return cond;
      elementItem->element.d.sq = LST_Create();
      if (elementItem->element.d.sq == NULL)
        return COND_PushCondition(DCM_LISTFAILURE,
                                  DCM_Message(DCM_LISTFAILURE), "readFile");

      localLength = elementItem->element.length;
      sequenceDone = (localLength == 0);

      while (!sequenceDone) {
        if (debug)
          fprintf(stderr, "Sequence Length: %d %x\n", localLength,
                  localLength);
        if (fileFlag) {
          if (fd != -1) {
            nBytes = read(fd, buf, 8);
          } else {
            cond = (*object)->rd((*object)->userCtx, buf, 8, &nBytes);
          }
          if (nBytes != 8)
            return COND_PushCondition(DCM_FILEACCESSERROR,
                                      DCM_Message(DCM_FILEACCESSERROR), name,
                                      "readFile");
          ptr = buf;
        }
        if (size != (long) DCM_UNSPECIFIEDLENGTH)
          size -= 8;
        fileOffset += (off_t) 8;
        if (scannedLength != NULL)
          (*scannedLength) += 8;
        (*object)->objectSize += 8;
        if (localLength != DCM_UNSPECIFIEDLENGTH)
          localLength -= 8;

        if (byteOrder == BYTEORDER_SAME) {
          GET_SHORT_SAME_ORDER(ptr, tagGroup);
          GET_SHORT_SAME_ORDER(ptr + 2, tagElement);
          tagE.tag = DCM_MAKETAG(tagGroup, tagElement);
          GET_LONG_SAME_ORDER(ptr + 4, tagE.length);
        } else {
          GET_SHORT_REVERSE_ORDER(ptr, tagGroup);
          GET_SHORT_REVERSE_ORDER(ptr + 2, tagElement);
          tagE.tag = DCM_MAKETAG(tagGroup, tagElement);
          GET_LONG_REVERSE_ORDER(ptr + 4, tagE.length);
        }
        ptr += 8;
        if (debug)
          fprintf(stderr, "Sequence item: %4x %4x %d (%x)\n",
                  tagGroup, tagElement, tagE.length, tagE.length);
        if (tagE.tag == DCM_DLMITEM) {
          /*        if (size != DCM_UNSPECIFIEDLENGTH)
                    size -= 8;
          */
          /*        fileOffset += 8;
           */
          cond = readFile(name,
                          (fileFlag) ? NULL : ptr,
                          fd, tagE.length,
                          fileOffset, recursionLevel + 1, opt,
                          &sequenceObject, &sequenceLength,
                          remainOpenFlag, ctx, rd, sk);
          if (cond == DCM_NORMAL) {
            sequenceItem = CTN_MALLOC(sizeof(*sequenceItem));
            if (sequenceItem == NULL)
              return COND_PushCondition(DCM_MALLOCFAILURE,
                                        DCM_Message(DCM_MALLOCFAILURE),
                                        sizeof(*sequenceItem), "readFile");

            sequenceItem->object = sequenceObject;
            cond = LST_Enqueue(&elementItem->element.d.sq,
                               sequenceItem);
            if (cond != LST_NORMAL)
              return COND_PushCondition(DCM_LISTFAILURE,
                                        DCM_Message(DCM_LISTFAILURE),
                                        "readFile");
            if (size != (long) DCM_UNSPECIFIEDLENGTH)
              size -= sequenceLength;
            fileOffset += (off_t) sequenceLength;
            if (scannedLength != NULL)
              *scannedLength += sequenceLength;
            (*object)->objectSize += sequenceLength;
            if (localLength != DCM_UNSPECIFIEDLENGTH)
              localLength -= sequenceLength;
            ptr += sequenceLength;
          } else
            return cond;
        } else {
          sequenceDone = TRUE;
        }
        if (localLength == 0)
          sequenceDone = TRUE;
      }
    } else {
      pixelFlag = (e.tag == DCM_PXLPIXELDATA);
      cond = newElementItem(&e, (pixelFlag == FALSE), &elementItem);
      if (cond != DCM_NORMAL) {
        (void) DCM_CloseObject(callerObject);
        return cond;
      }
      if (pixelFlag) {
        if (fileFlag)
          *remainOpenFlag = TRUE;
        elementItem->byteOrder = byteOrder;
        elementItem->dataOffset = fileOffset;
        elementItem->currentOffset = 0;
        if (fileFlag)
          elementItem->element.d.ot = NULL;
        else
          elementItem->element.d.ot = (void *) ptr;
        if ((*object)->pixelBitsAllocated == 8)
          elementItem->element.representation = DCM_OB;
        else
          elementItem->element.representation = DCM_OW;
        if (fileFlag) {
          if (fd != -1) {
            (void) lseek(fd, (off_t) elementItem->element.length, SEEK_CUR);
          } else {
            (*object)->sk((*object)->userCtx,
                          elementItem->element.length, SEEK_CUR);
          }
          (*object)->fd = fd;
        }
      } else {
        if (fileFlag) {
          if (fd != -1) {
            nBytes = read(fd, elementItem->element.d.ot,
                          (int) elementItem->element.length);
          } else {
            cond = (*object)->rd((*object)->userCtx,
                                 elementItem->element.d.ot,
                                 (long) elementItem->element.length, &nBytes);
          }
          if (nBytes != (int) elementItem->element.length) {
            (void) DCM_CloseObject(callerObject);
            return COND_PushCondition(DCM_FILEACCESSERROR,
                                      DCM_Message(DCM_FILEACCESSERROR),
                                      name, "readFile");
          }
        } else {
          (void) memcpy(elementItem->element.d.ot, ptr,
                        elementItem->element.length);
          ptr += elementItem->originalDataLength;
        }

#ifdef LITTLE_ENDIAN_ARCHITECTURE
        if (elementItem->element.representation == DCM_AT)
          swapATGroupElement(&elementItem->element);
#endif
        if (byteOrder != BYTEORDER_SAME)
          swapInPlace(object, &elementItem->element);
        if (convertFlag) {
          cond = verifyFormat(elementItem);
          if (cond != DCM_NORMAL)
            return cond;
        }
      }
      if (size != (long) DCM_UNSPECIFIEDLENGTH)
        size -= elementItem->originalDataLength;
      fileOffset += (off_t) elementItem->originalDataLength;
      if (scannedLength != NULL)
        (*scannedLength) += elementItem->originalDataLength;

      elementItem->paddedDataLength = elementItem->element.length;
      if (elementItem->paddedDataLength & 1)
        elementItem->paddedDataLength += 1;
      (*object)->objectSize += elementItem->paddedDataLength;
    }

    computeVM(object, &elementItem->element);

    if ((long) DCM_TAG_GROUP(e.tag) == lastGroup) {
      if ((long) DCM_TAG_ELEMENT(e.tag) <= lastElement)
        return COND_PushCondition(DCM_ELEMENTOUTOFORDER,
                                  DCM_Message(DCM_ELEMENTOUTOFORDER),
                                  group, element, "readFile");
    } else if ((long) DCM_TAG_GROUP(e.tag) > lastGroup) {}
    else {
      return COND_PushCondition(DCM_ELEMENTOUTOFORDER,
                                DCM_Message(DCM_ELEMENTOUTOFORDER),
                                group, element,
                                "readFile");
    }
    lastGroup = (long) group;
    lastElement = (long) element;

    if (groupItem == NULL)
      createGroupFlag = TRUE;
    else if (groupItem->group != group)
      createGroupFlag = TRUE;
    else
      createGroupFlag = FALSE;

    if (createGroupFlag == TRUE) {
      groupItem = CTN_MALLOC(sizeof(*groupItem));
      if (groupItem == NULL) {
        (void) DCM_CloseObject(callerObject);
        return COND_PushCondition(DCM_ELEMENTCREATEFAILED,
                                  DCM_Message(DCM_ELEMENTCREATEFAILED),
                                  "readFile",
                                  group, 0xffff, sizeof(*groupItem));
      }
      groupItem->group = group;
      groupItem->baseLength = 0;
      groupItem->longVRAttributes = 0;
      groupItem->elementList = LST_Create();
      if (groupItem->elementList == NULL) {
        (void) DCM_CloseObject(callerObject);
        return COND_PushCondition(DCM_LISTFAILURE,
                                  DCM_Message(DCM_LISTFAILURE),
                                  "readFile");
      }
      if (LST_Enqueue(&(*object)->groupList, groupItem) != LST_NORMAL) {
        (void) DCM_CloseObject(callerObject);
        return COND_PushCondition(DCM_LISTFAILURE,
                                  DCM_Message(DCM_LISTFAILURE),
                                  "readFile");
      }
    }
    if (element != 0x0000)
      groupItem->baseLength += 8 + elementItem->paddedDataLength;
    if ((element == 0x0000) && ((*object)->groupLengthFlag == FALSE)) {
      CTN_FREE(elementItem);
    } else {
      cond = LST_Enqueue(&groupItem->elementList, elementItem);
      if (cond != LST_NORMAL) {
        (void) DCM_CloseObject(callerObject);
        return COND_PushCondition(DCM_LISTFAILURE,
                                  DCM_Message(DCM_LISTFAILURE),
                                  "readFile");
      }
      cond = updateObjectType(object, &elementItem->element); /* repair */

      cond = updateSpecialElements(object, elementItem);  /* repair */
    }

    if (size == 0)
      done = TRUE;

#ifdef DEBUG
    if (debug) {
      /*lint -e644 */
      (void) fprintf(stderr,
                     "Address: %px Group %2x, element %2x, length %ld ",
                     elementItem,
                     DCM_TAG_GROUP(elementItem->element.tag),
                     DCM_TAG_ELEMENT(elementItem->element.tag),
                     elementItem->element.length);
      /*lint +e644 */
      (void) fprintf(stderr, "Object size: %d\n", (*object)->objectSize);
    }
#endif
  }

  groupItem = LST_Head(&(*object)->groupList);
  if (groupItem != NULL) {
    (void) LST_Position(&(*object)->groupList, groupItem);
    while (groupItem != NULL) {
      elementItem = LST_Head(&groupItem->elementList);
      if (elementItem != NULL) {
        if (DCM_TAG_ELEMENT(elementItem->element.tag) == 0x0000) {
          *elementItem->element.d.ul = groupItem->baseLength;
        }
      }
      groupItem = LST_Next(&(*object)->groupList);
    }
  }
  return DCM_NORMAL;
}
#endif

static CONDITION
readPreamble(const char *name, unsigned char **ptr, int fd, U32 * size,
             off_t * fileOffset, CTNBOOLEAN knownLength,
             PRIVATE_OBJECT ** object, U32 * scannedLength) {
  int nBytes,
  tmp;
  CONDITION cond;
  char label[4];

  if (*size == 0)
    return DCM_STREAMCOMPLETE;

  if ((*size < DCM_PREAMBLELENGTH + 4) && knownLength) {
    if (debug)
      (void) DCM_DumpElements((DCM_OBJECT **) object, 0);
    (void) DCM_CloseObject((DCM_OBJECT **) object);
    return COND_PushCondition(DCM_ILLEGALSTREAMLENGTH,
                              DCM_Message(DCM_ILLEGALSTREAMLENGTH), *size,
                              "readPreamble");
  }
  if (*ptr == NULL) {
    if (fd != -1) {
      nBytes = read(fd, (*object)->preamble, DCM_PREAMBLELENGTH);
      nBytes += read(fd, label, sizeof(label));
    } else {
      cond = (*object)->rd((*object)->userCtx, (*object)->preamble,
                           DCM_PREAMBLELENGTH, &nBytes);
      cond = (*object)->rd((*object)->userCtx, label,
                           sizeof(label), &tmp);
      nBytes += tmp;
    }

    if (nBytes != DCM_PREAMBLELENGTH + sizeof(label))
      return COND_PushCondition(DCM_FILEACCESSERROR,
                                DCM_Message(DCM_FILEACCESSERROR), name,
                                "readPreamble");
  } else {
    (void) memcpy((*object)->preamble, *ptr, DCM_PREAMBLELENGTH);
    (void) memcpy(label, (*ptr) + DCM_PREAMBLELENGTH, sizeof(label));
  }

  if (knownLength)
    *size -= DCM_PREAMBLELENGTH + sizeof(label);
  if (fileOffset != NULL) {
    *fileOffset += (off_t) DCM_PREAMBLELENGTH + sizeof(label);
  }
  if (*ptr != NULL)
    *ptr += DCM_PREAMBLELENGTH + sizeof(label);
  (*object)->objectSize += DCM_PREAMBLELENGTH + sizeof(label);

  if (strncmp(label, "DICM", 4) != 0)
    return 0;

  (*object)->preambleFlag = TRUE;
  return DCM_NORMAL;
}


static CONDITION
readGroupElement(const char *name, unsigned char **ptr, int fd, U32 * size,
                 off_t * fileOffset, CTNBOOLEAN knownLength, int byteOrder,
                 CTNBOOLEAN explicitVR, CTNBOOLEAN acceptVRMismatch,
                 PRIVATE_OBJECT ** object, U32 * scannedLength,
                 DCM_ELEMENT * e) {
  unsigned char *localPtr;
  unsigned char buf[4];
  int nBytes;
  CONDITION cond;
  unsigned short group,
  element;

  if (*size == 0)
    return DCM_STREAMCOMPLETE;

  if ((*size < 4) && knownLength) {
    if (debug)
      (void) DCM_DumpElements((DCM_OBJECT **) object, 0);
    (void) DCM_CloseObject((DCM_OBJECT **) object);
    return COND_PushCondition(DCM_ILLEGALSTREAMLENGTH,
                              DCM_Message(DCM_ILLEGALSTREAMLENGTH), *size,
                              "readFile");
  }
  if (*ptr == NULL) {
    if (fd != -1) {
#ifdef _WIN32
      nBytes = _read(fd, buf, 4);
#else
      nBytes = read(fd, buf, 4);
#endif
    } else {
      cond = (*object)->rd((*object)->userCtx, buf, 4, &nBytes);
    }

    if (nBytes != 4) {
      perror("");
      printf("Bytes read: %d\n", nBytes);

      return COND_PushCondition(DCM_FILEACCESSERROR,
                                DCM_Message(DCM_FILEACCESSERROR), name,
                                "readGroupElement");
    }
    localPtr = buf;
  } else {
    localPtr = *ptr;
  }

  if (knownLength)
    *size -= 4;
  if (fileOffset != NULL) {
    *fileOffset += (off_t) 4;
  }
  if (scannedLength != NULL)
    (*scannedLength) += 4;
  (*object)->objectSize += 4;

  if (byteOrder == BYTEORDER_SAME) {
    GET_SHORT_SAME_ORDER(localPtr, group);
    GET_SHORT_SAME_ORDER(localPtr + 2, element);
    e->tag = DCM_MAKETAG(group, element);
  } else {
    GET_SHORT_REVERSE_ORDER(localPtr, group);
    GET_SHORT_REVERSE_ORDER(localPtr + 2, element);
    e->tag = DCM_MAKETAG(group, element);
  }
  if (*ptr != NULL)
    *ptr += 4;

  if (debug)
    fprintf(stderr, "%04x %04x ", group, element);

  cond = DCM_LookupElement(e);
  if (cond != DCM_NORMAL)
    (void) COND_PopCondition(0);
  if (e->representation == DCM_CTX)
    ctxSensitiveLookup(object, e);

  return DCM_NORMAL;
}

static CONDITION
readVRLength(const char *name, unsigned char **ptr, int fd, U32 * size,
             off_t * fileOffset,
             CTNBOOLEAN knownLength, int byteOrder, CTNBOOLEAN explicitVR,
             CTNBOOLEAN acceptVRMismatch,
             PRIVATE_OBJECT ** object, U32 * scannedLength, DCM_ELEMENT * e) {
  unsigned char *localPtr;
  unsigned char buf[4];
  char vrCode[3];
  VRMAP *vrPtr;
  int nBytes;
  CONDITION cond;
  CTNBOOLEAN calculatedLength = FALSE;

  if (*size == 0)
    return DCM_STREAMCOMPLETE;

  if ((*size < 4) && knownLength) {
    if (debug)
      (void) DCM_DumpElements((DCM_OBJECT **) object, 0);
    (void) DCM_CloseObject((DCM_OBJECT **) object);
    return COND_PushCondition(DCM_ILLEGALSTREAMLENGTH,
                              DCM_Message(DCM_ILLEGALSTREAMLENGTH), *size,
                              "readVRLength");
  }
  if (*ptr == NULL) {
    if (fd != -1) {
#ifdef _WIN32
      nBytes = _read(fd, buf, 4);
#else
      nBytes = read(fd, buf, 4);
#endif
    } else {
      cond = (*object)->rd((*object)->userCtx, buf, 4, &nBytes);
    }

    if (nBytes != 4)
      return COND_PushCondition(DCM_FILEACCESSERROR,
                                DCM_Message(DCM_FILEACCESSERROR), name,
                                "readVRLength");
    localPtr = buf;
  } else {
    localPtr = *ptr;
    *ptr += 4;
  }

  if (knownLength)
    *size -= 4;
  if (fileOffset != NULL) {
    *fileOffset += (off_t) 4;
  }
  if (scannedLength != NULL)
    (*scannedLength) += 4;
  (*object)->objectSize += 4;

  e->length = 0;
  if (e->representation == DCM_DLM) {
    explicitVR = FALSE; /* Special rule for delimitors */
  }
  if (explicitVR) {
    vrCode[0] = localPtr[0];
    vrCode[1] = localPtr[1];
    vrCode[2] = '\0';
    vrPtr = lookupVRCode(vrCode);
    if (vrPtr == NULL)
      return COND_PushCondition(DCM_UNRECOGNIZEDVRCODE,
                                DCM_Message(DCM_UNRECOGNIZEDVRCODE), vrCode,
                                "readVRLength");

    if (vrPtr->representation != e->representation) {
      if (vrPtr->representation == DCM_OB) {
        /* This is probably from the waveform supplement where they */
        /* transmit as OB and expect us to pull it out later */
        /* We will just keep our VR which was based on context in */
        /* the object */
        e->representation = vrPtr->representation;
      } else if (e->representation == DCM_UN ||
                 e->representation == DCM_CTX ||
                 e->representation == DCM_RET ||
                 vrPtr->representation == DCM_OW ||
                 acceptVRMismatch) {  /* Believe input */
        e->representation = vrPtr->representation;
      } else {
        if (e->tag != DCM_PXLPIXELDATA)
          return COND_PushCondition(DCM_VRMISMATCH,
                                    DCM_Message(DCM_VRMISMATCH),
                                    vrCode, e->tag);
      }
    }
    if (vrPtr->representation != DCM_OW &&
        vrPtr->representation != DCM_OB &&
        vrPtr->representation != DCM_UN &&
        vrPtr->representation != DCM_UT &&
        vrPtr->representation != DCM_SQ) {
      unsigned short shortLength;
      if (byteOrder == BYTEORDER_SAME) {
        GET_SHORT_SAME_ORDER(localPtr + 2, shortLength);
      } else {
        GET_SHORT_REVERSE_ORDER(localPtr + 2, shortLength);
      }
      e->length = shortLength;
      /*      if (*ptr != NULL)
       *ptr += 4;
       */
      calculatedLength = TRUE;
    } else {
      if ((*size < 4) && knownLength) {
        if (debug)
          (void) DCM_DumpElements((DCM_OBJECT **) object, 0);
        (void) DCM_CloseObject((DCM_OBJECT **) object);
        return COND_PushCondition(DCM_ILLEGALSTREAMLENGTH,
                                  DCM_Message(DCM_ILLEGALSTREAMLENGTH),
                                  *size,
                                  "readVRLength");
      }
      if (*ptr == NULL) {
        if (fd != -1) {
#ifdef _WIN32
          nBytes = _read(fd, buf, 4);
#else
          nBytes = read(fd, buf, 4);
#endif
        } else {
          cond = (*object)->rd((*object)->userCtx, buf, 4, &nBytes);
        }

        if (nBytes != 4)
          return COND_PushCondition(DCM_FILEACCESSERROR,
                                    DCM_Message(DCM_FILEACCESSERROR),
                                    name,
                                    "readVRLength");
        localPtr = buf;
      } else {
        localPtr = *ptr;
        *ptr += 4;
      }

      if (knownLength)
        *size -= 4;
      if (fileOffset != 0) {
        *fileOffset += (off_t) 4;
      }
      if (scannedLength != NULL)
        (*scannedLength) += 4;
      (*object)->objectSize += 4;
    }
  }
  if (!calculatedLength) {
    if (byteOrder == BYTEORDER_SAME) {
      GET_LONG_SAME_ORDER(localPtr, e->length);
    } else {
      GET_LONG_REVERSE_ORDER(localPtr, e->length);
    }
    /*  if (*ptr != NULL)
     *ptr += 4;
     */
  }
  if (debug) {
    char localVR[10];
    mapVRtoASCII(e->representation, localVR);
    fprintf(stderr, "%2s %6d %06x %s\n",
            localVR, e->length, (unsigned int) *fileOffset,
            e->description);
  }
  if (((e->length & 1) != 0) && (e->length != DCM_UNSPECIFIEDLENGTH)) {
    if (debug)
      (void) DCM_DumpElements((DCM_OBJECT **) object, 0);
    (void) DCM_CloseObject((DCM_OBJECT **) object);
    return COND_PushCondition(DCM_UNEVENELEMENTLENGTH,
                              DCM_Message(DCM_UNEVENELEMENTLENGTH),
                              DCM_TAG_GROUP(e->tag), DCM_TAG_ELEMENT(e->tag),
                              e->length, "readFile");
  }
  if ((e->length != (U32) DCM_UNSPECIFIEDLENGTH) \
      && (e->length > (U32) (*size))) {
    if (debug)
      (void) DCM_DumpElements((DCM_OBJECT **) object, 0);
    (void) DCM_CloseObject((DCM_OBJECT **) object);
    return COND_PushCondition(DCM_ELEMENTLENGTHERROR,
                              DCM_Message(DCM_ELEMENTLENGTHERROR),
                              DCM_TAG_GROUP(e->tag), DCM_TAG_ELEMENT(e->tag),
                              e->length, *size, "readFile");
  }
  return DCM_NORMAL;
}

static CONDITION
readSequence(const char *name, unsigned char **ptr, int fd, U32 * size,
             off_t * fileOffset, int recursionLevel, unsigned long opt,
             int byteOrder, CTNBOOLEAN explicitVR, CTNBOOLEAN acceptVRMismatch,
             CTNBOOLEAN fileFlag, CTNBOOLEAN * remainOpenFlag,
             CTNBOOLEAN convertFlag, PRIVATE_OBJECT ** object,
             U32 * scannedLength, DCM_ELEMENT * e,
             PRV_ELEMENT_ITEM ** elementItem) {
  CTNBOOLEAN knownLength = TRUE;
  CONDITION cond;
  U32 sequenceLength;

  U32 localLength;
  CTNBOOLEAN sequenceDone;
  DCM_ELEMENT tagE;
  DCM_OBJECT
  * sequenceObject;
  DCM_SEQUENCE_ITEM *sequenceItem;
  CONDITION flag;
  unsigned char *localPtr;
  off_t itemTagOffset=0;

  if (*size == (long) DCM_UNSPECIFIEDLENGTH)
    knownLength = FALSE;

  cond = newElementItem(e, FALSE, elementItem);
  if (cond != DCM_NORMAL)
    return cond;
  (*elementItem)->element.d.sq = LST_Create();
  if ((*elementItem)->element.d.sq == NULL)
    return COND_PushCondition(DCM_LISTFAILURE,
                              DCM_Message(DCM_LISTFAILURE), "readSequence");

  localLength = (*elementItem)->element.length;
  sequenceDone = (localLength == 0);

  while (!sequenceDone) {
    if (debug)
      fprintf(stderr, "Sequence Length: %d %x\n", localLength,
              (unsigned int) localLength);

    sequenceLength = 0;
    if (fileOffset != NULL) {
      itemTagOffset = *fileOffset;
    }

    flag = readGroupElement(name, ptr, fd, &localLength,
                            fileOffset, knownLength,
                            byteOrder, explicitVR, acceptVRMismatch,
                            object, &sequenceLength, &tagE);
    if (flag == DCM_STREAMCOMPLETE)
      break;
    else if (flag != DCM_NORMAL)
      return flag;

    flag = readVRLength(name, ptr, fd, &localLength, fileOffset, knownLength,
                        byteOrder, explicitVR, acceptVRMismatch, object,
                        &sequenceLength, &tagE);
    if (flag != DCM_NORMAL)
      return flag;

    if (*size != (long) DCM_UNSPECIFIEDLENGTH)
      *size -= sequenceLength;
    if (scannedLength != NULL)
      *scannedLength += sequenceLength;

    sequenceLength = 0;


    if (debug)
      fprintf(stderr, "Sequence item: %4x %4x %d (%x)\n",
              DCM_TAG_GROUP(tagE.tag),
              DCM_TAG_ELEMENT(tagE.tag), tagE.length,
              (unsigned int) tagE.length);
    if (tagE.tag == DCM_DLMITEM) {
      localPtr = *ptr;
      cond = readFile1(name,
                       localPtr,
                       fd, tagE.length,
                       fileOffset, recursionLevel + 1, opt,
                       object, &sequenceObject, &sequenceLength,
                       remainOpenFlag, (*object)->userCtx, (*object)->rd,
                       (*object)->sk);
      if (!fileFlag)
        *ptr = localPtr + sequenceLength;
      if (cond == DCM_NORMAL) {
        sequenceItem = (DCM_SEQUENCE_ITEM*)CTN_MALLOC(sizeof(*sequenceItem));
        if (sequenceItem == NULL)
          return COND_PushCondition(DCM_MALLOCFAILURE,
                                    DCM_Message(DCM_MALLOCFAILURE),
                                    sizeof(*sequenceItem), "readFile");

        ((PRIVATE_OBJECT *) sequenceObject)->offset = itemTagOffset;
        sequenceItem->object = sequenceObject;
        cond = LST_Enqueue(&(*elementItem)->element.d.sq,
                           sequenceItem);
        if (cond != LST_NORMAL)
          return COND_PushCondition(DCM_LISTFAILURE,
                                    DCM_Message(DCM_LISTFAILURE), "readFile");
        if (*size != (long) DCM_UNSPECIFIEDLENGTH)
          *size -= sequenceLength;
        if (scannedLength != NULL)
          *scannedLength += sequenceLength;
        (*object)->objectSize += sequenceLength;
        if (localLength != DCM_UNSPECIFIEDLENGTH)
          localLength -= sequenceLength;
      } else
        return cond;
    } else {
      sequenceDone = TRUE;
    }
    if (localLength == 0)
      sequenceDone = TRUE;
  }
  return DCM_NORMAL;
}

static CONDITION
scanCompressedPixels(const char *name, unsigned char **ptr, int fd, U32 * size,
                     off_t * fileOffset, int recursionLevel, unsigned long opt,
                     int byteOrder, CTNBOOLEAN explicitVR,
                     CTNBOOLEAN acceptVRMismatch,
                     CTNBOOLEAN fileFlag, CTNBOOLEAN * remainOpenFlag,
                     CTNBOOLEAN convertFlag, PRIVATE_OBJECT ** object,
                     U32 * scannedLength, DCM_ELEMENT * e,
                     PRV_ELEMENT_ITEM ** elementItem) {
  CTNBOOLEAN knownLength = TRUE;
  U32 sequenceLength;
  U32 scannedBytes = 0;

  U32 localLength;
  CTNBOOLEAN sequenceDone;
  DCM_ELEMENT tagE;
  CONDITION flag;
  unsigned char *localPtr;

  if (*size == (long) DCM_UNSPECIFIEDLENGTH)
    knownLength = FALSE;

  localLength = (*elementItem)->element.length;
  sequenceDone = (localLength == 0);

  while (!sequenceDone) {
    sequenceLength = 0;
    flag = readGroupElement(name, ptr, fd, &localLength, fileOffset,
                            FALSE, byteOrder, explicitVR, acceptVRMismatch,
                            object, &sequenceLength, &tagE);
    if (flag == DCM_STREAMCOMPLETE)
      break;
    else if (flag != DCM_NORMAL)
      return flag;

    flag = readVRLength(name, ptr, fd, &localLength, fileOffset, knownLength,
                        byteOrder, explicitVR, acceptVRMismatch, object,
                        &sequenceLength, &tagE);
    if (flag != DCM_NORMAL)
      return flag;

    if (*size != (long) DCM_UNSPECIFIEDLENGTH)
      *size -= sequenceLength;
    if (scannedLength != NULL)
      *scannedLength += sequenceLength;
    scannedBytes += sequenceLength;

    if (debug)
      fprintf(stderr, "Sequence item: %4x %4x %d (%x)\n",
              DCM_TAG_GROUP(tagE.tag),
              DCM_TAG_ELEMENT(tagE.tag), tagE.length,
              (unsigned int) tagE.length);
    if (tagE.tag == DCM_DLMITEM) {
      localPtr = *ptr;
      if (tagE.length != 0) {
        lseek(fd, tagE.length, SEEK_CUR);
        *fileOffset += tagE.length;
        if (*size != (long) DCM_UNSPECIFIEDLENGTH)
          *size -= tagE.length;
        if (scannedLength != NULL)
          *scannedLength += tagE.length;
      }
    } else {
      sequenceDone = TRUE;
    }
    if (localLength == 0)
      sequenceDone = TRUE;

    if (debug)
      fprintf(stderr, "Scanned Bytes: %d\n", scannedBytes);
  }
  if ((scannedBytes & 1) != 0) {
    lseek(fd, 1, SEEK_CUR);
    *fileOffset += 1;
    if (*size != (long) DCM_UNSPECIFIEDLENGTH)
      *size -= 1;
  }
  return DCM_NORMAL;
}

static CONDITION
readData(const char *name, unsigned char **ptr, int fd, U32 * size,
         off_t * fileOffset,
         CTNBOOLEAN knownLength, int byteOrder, CTNBOOLEAN explicitVR,
         CTNBOOLEAN acceptVRMismatch,
         CTNBOOLEAN fileFlag, CTNBOOLEAN * remainOpenFlag,
         CTNBOOLEAN convertFlag, PRIVATE_OBJECT ** object,
         U32 * scannedLength, DCM_ELEMENT * e,
         PRV_ELEMENT_ITEM ** elementItem) {
  CTNBOOLEAN pixelFlag;
  CONDITION cond;
  int nBytes;

  pixelFlag = (e->tag == DCM_PXLPIXELDATA);
  cond = newElementItem(e, ((pixelFlag == FALSE) || (fileFlag == FALSE)),
                        elementItem);
  if (cond != DCM_NORMAL) {
    (void) DCM_CloseObject((DCM_OBJECT **) object);
    return cond;
  }
  if (pixelFlag && fileFlag) {
    /*  if (fileFlag) */
    *remainOpenFlag = TRUE;
    (*elementItem)->byteOrder = byteOrder;
    (*elementItem)->dataOffset = *fileOffset;
    (*elementItem)->currentOffset = 0;
    (*elementItem)->element.d.ot = NULL;
    if ((*object)->pixelBitsAllocated == 8)
      (*elementItem)->element.representation = DCM_OB;
    else
      (*elementItem)->element.representation = DCM_OW;
    if (fileFlag) {
      if (fd != -1) {
        if ((*elementItem)->element.length != DCM_UNSPECIFIEDLENGTH)
          (void) lseek(fd,
                       (off_t) (*elementItem)->element.length,
                       SEEK_CUR);
        else {
          U32 l1 = 0;
          U32 s1;
          off_t f1 = 0;

          s1 = *size;
          scanCompressedPixels("", ptr, fd,
                               &s1, /* size */
                               &f1, /* fileOffset */
                               0, 0,
                               byteOrder, explicitVR,
                               acceptVRMismatch,
                               fileFlag, remainOpenFlag,
                               convertFlag, object,
                               &l1, /* scannedLength */
                               e, elementItem);
          (*elementItem)->originalDataLength = l1;
          (*elementItem)->paddedDataLength = l1;
        }
      } else {
        (*object)->sk((*object)->userCtx,
                      (*elementItem)->element.length, SEEK_CUR);
      }
      (*object)->fd = fd;
    }
  } else {
    if (fileFlag) {
      if (fd != -1) {
        nBytes = read(fd, (*elementItem)->element.d.ot,
                      (int) (*elementItem)->element.length);
      } else {
        cond = (*object)->rd((*object)->userCtx,
                             (*elementItem)->element.d.ot,
                             (long) (*elementItem)->element.length, &nBytes);
      }
      if (nBytes != (int) (*elementItem)->element.length) {
        (void) DCM_CloseObject((DCM_OBJECT **) object);
        return COND_PushCondition(DCM_FILEACCESSERROR,
                                  DCM_Message(DCM_FILEACCESSERROR),
                                  name, "readFile");
      }
    } else {
      (void) memcpy((*elementItem)->element.d.ot, *ptr,
                    (*elementItem)->element.length);
      *ptr += (*elementItem)->originalDataLength;
    }
#ifdef LITTLE_ENDIAN_ARCHITECTURE
    if ((*elementItem)->element.representation == DCM_AT)
      swapATGroupElement(&(*elementItem)->element);
#endif
    if (byteOrder != BYTEORDER_SAME)
      swapInPlace(object, &(*elementItem)->element);
    if (convertFlag) {
      cond = verifyFormat(*elementItem);
      if (cond != DCM_NORMAL)
        return cond;
    }
  }
  if (*size != (long) DCM_UNSPECIFIEDLENGTH)
    *size -= (*elementItem)->originalDataLength;
  if (fileOffset != NULL) {
    *fileOffset += (off_t) (*elementItem)->originalDataLength;
  }
  if (scannedLength != NULL)
    (*scannedLength) += (*elementItem)->originalDataLength;

  if ((*elementItem)->element.length != DCM_UNSPECIFIEDLENGTH) {
    (*elementItem)->paddedDataLength = (*elementItem)->element.length;
  }
  if (((*elementItem)->paddedDataLength != DCM_UNSPECIFIEDLENGTH) &&
      ((*elementItem)->paddedDataLength & 1) )
    (*elementItem)->paddedDataLength += 1;
  (*object)->objectSize += (*elementItem)->paddedDataLength;

  return DCM_NORMAL;

}

static CONDITION
checkAttributeOrder(DCM_ELEMENT * e, long *lastGroup, long *lastElement,
                    CTNBOOLEAN allowRepeatElements) {
  unsigned short group;
  unsigned short element;

  group = DCM_TAG_GROUP(e->tag);
  element = DCM_TAG_ELEMENT(e->tag);

  if ((long) group == *lastGroup) {
    if (((long) element == *lastElement) && allowRepeatElements) {
      return DCM_REPEATEDELEMENT;
    }
    if ((long) element <= *lastElement)
      return COND_PushCondition(DCM_ELEMENTOUTOFORDER,
                                DCM_Message(DCM_ELEMENTOUTOFORDER),
                                group, element, "checkAttributeOrder");
  } else if ((long) group > *lastGroup) {}
  else {
    return COND_PushCondition(DCM_ELEMENTOUTOFORDER,
                              DCM_Message(DCM_ELEMENTOUTOFORDER),
                              group, element, "checkAttributeOrder");
  }
  *lastGroup = (long) group;
  *lastElement = (long) element;

  return DCM_NORMAL;
}

static CONDITION
handleGroupItem(PRIVATE_OBJECT ** obj, PRV_GROUP_ITEM ** groupItem,
                unsigned short group) {
  CTNBOOLEAN createGroupFlag;

  if (*groupItem == NULL)
    createGroupFlag = TRUE;
  else if ((*groupItem)->group != group)
    createGroupFlag = TRUE;
  else
    createGroupFlag = FALSE;

  if (createGroupFlag == TRUE) {
    *groupItem = (PRV_GROUP_ITEM*)CTN_MALLOC(sizeof(**groupItem));
    if (*groupItem == NULL) {
      (void) DCM_CloseObject((DCM_OBJECT **) obj);
      return COND_PushCondition(DCM_ELEMENTCREATEFAILED,
                                DCM_Message(DCM_ELEMENTCREATEFAILED),
                                "handleGroupItem",
                                group, 0xffff, sizeof(**groupItem));
    }
    (*groupItem)->group = group;
    (*groupItem)->baseLength = 0;
    (*groupItem)->longVRAttributes = 0;
    (*groupItem)->elementList = LST_Create();
    if ((*groupItem)->elementList == NULL) {
      (void) DCM_CloseObject((DCM_OBJECT **) obj);
      return COND_PushCondition(DCM_LISTFAILURE,
                                DCM_Message(DCM_LISTFAILURE),
                                "handleGroupItem");
    }
    if (LST_Enqueue(&(*obj)->groupList, *groupItem) != LST_NORMAL) {
      (void) DCM_CloseObject((DCM_OBJECT **) obj);
      return COND_PushCondition(DCM_LISTFAILURE,
                                DCM_Message(DCM_LISTFAILURE),
                                "handleGroupItem");
    }
  }
  return DCM_NORMAL;
}

static CONDITION
readFile1(const char *name, unsigned char *callerBuf, int fd, U32 size,
          off_t * fileOffset, int recursionLevel,
          unsigned long opt, PRIVATE_OBJECT ** parentObject,
          DCM_OBJECT ** callerObject,
          U32 * scannedLength, CTNBOOLEAN * remainOpenFlag,
          void *ctx,
          CONDITION(*rd) (void *ctx, void *buf, int toRead, int *bytesRead),
          CONDITION(*sk) (void *ctx, int offset, int flag)) {
  CONDITION cond = DCM_NORMAL;
  int byteOrder = NATIVE_ORDER;
  long lastGroup = -1, lastElement = -1;
  U32 sequenceLength = 0, scannedSequenceLength = 0;
  PRIVATE_OBJECT** object = NULL;
  PRV_GROUP_ITEM* groupItem = NULL;
  DCM_ELEMENT e;
  CTNBOOLEAN convertFlag = FALSE,
    done = FALSE,
    knownLength = TRUE,
    explicitVR = FALSE,
    acceptVRMismatch = FALSE,
    part10Flag = FALSE;
  unsigned char*ptr = NULL;
  PRV_ELEMENT_ITEM* elementItem = NULL;
  CTNBOOLEAN fileFlag = TRUE;
  CONDITION flag = DCM_NORMAL;
  CTNBOOLEAN allowRepeatElements = FALSE;

  ptr = callerBuf;
  if (ptr != NULL)
    fileFlag = FALSE;

  if ((opt & DCM_FILEFORMATMASK) == DCM_PART10FILE) {
    part10Flag = TRUE;
    opt &= ~DCM_ORDERMASK;
    opt &= ~DCM_FILEFORMATMASK;
    opt |= DCM_EXPLICITLITTLEENDIAN;
  }
  if ((opt & DCM_SPECIALFORMATMASK) == DCM_EFILM) {
    part10Flag = TRUE;
    opt &= ~DCM_ORDERMASK;
    opt &= ~DCM_FILEFORMATMASK;
    opt |= DCM_ORDERLITTLEENDIAN;
  }
  if ((opt & DCM_REPEATELEMENTSMASK) == DCM_ALLOWREPEATELEMENTS) {
    allowRepeatElements = TRUE;
  }

  switch (opt & DCM_ORDERMASK) {
  case DCM_ORDERNATIVE:
    byteOrder = NATIVE_ORDER;
    break;
  case DCM_ORDERLITTLEENDIAN:
    byteOrder = LITTLE_ORDER;
    break;
  case DCM_EXPLICITLITTLEENDIAN:
    byteOrder = LITTLE_ORDER;
    explicitVR = TRUE;
    break;
  case DCM_ORDERBIGENDIAN:
    byteOrder = BIG_ORDER;
    break;
  case DCM_EXPLICITBIGENDIAN:
    byteOrder = BIG_ORDER;
    explicitVR = TRUE;
    break;
  default:
    byteOrder = NATIVE_ORDER;
    break;
  }
  if ((opt & DCM_CONVERTMASK) == DCM_FORMATCONVERSION)
    convertFlag = TRUE;
  if ((opt & DCM_VRMASK) == DCM_ACCEPTVRMISMATCH)
    acceptVRMismatch = TRUE;

  if (scannedLength != NULL)
    *scannedLength = 0;

  cond = DCM_CreateObject(callerObject, opt);
  if (cond != DCM_NORMAL)
    return cond;

  object = (PRIVATE_OBJECT **) callerObject;
  if (fileFlag)
    strcpy((*object)->fileName, name);

  (*object)->fd = -1;
  (*object)->rd = rd;
  (*object)->sk = sk;
  (*object)->userCtx = ctx;
  (*object)->dataOptions = 0;
  if (size == (long) DCM_UNSPECIFIEDLENGTH)
    knownLength = FALSE;

  if ((fileFlag) \
      && ((opt & DCM_DELETEMASK) == DCM_DELETEONCLOSE)\
      && (recursionLevel == 0))
    (*object)->deleteFlag = TRUE;

  if (parentObject != NULL)
    (*object)->pixelRepresentation = (*parentObject)->pixelRepresentation;

  if (recursionLevel == 0 && part10Flag) {
    flag = readPreamble(name, &ptr, fd, &size, fileOffset, knownLength,
                        object, scannedLength);
    if (flag != DCM_NORMAL)
      goto abort;
  }
  while (!done) {
    flag = readGroupElement(name, &ptr, fd, &size, fileOffset, knownLength,
                            byteOrder, explicitVR, acceptVRMismatch, object,
                            scannedLength, &e);
    if (flag == DCM_STREAMCOMPLETE)
      break;
    else if (flag != DCM_NORMAL)
      goto abort;
#if 0
    if (e.tag == DCM_MAKETAG(0x7fe0, 0x0010)) {
      fprintf(stderr, "Found pixels\n");
    }
#endif
    flag = readVRLength(name, &ptr, fd, &size, fileOffset, knownLength,
                        byteOrder, explicitVR, acceptVRMismatch, object,
                        scannedLength, &e);
    if (flag != DCM_NORMAL)
      goto abort;

    if ((e.representation == DCM_UN) &&
        (e.length == DCM_UNSPECIFIEDLENGTH)) {
      e.representation = DCM_SQ;
    }
#ifndef SMM
    if ((e.tag == DCM_DLMITEMDELIMITATIONITEM) ||
        (e.tag == DCM_DLMSEQUENCEDELIMITATIONITEM)) {
      return DCM_NORMAL;
    }
#else
    if (e.tag == DCM_DLMITEMDELIMITATIONITEM) {
      (*object)->objectSize -= 8;
      return DCM_NORMAL;
    }
    if (e.tag == DCM_DLMSEQUENCEDELIMITATIONITEM)
      return DCM_NORMAL;
#endif

    if (e.representation == DCM_SQ) {
      sequenceLength = e.length;
      scannedSequenceLength = 0;
      flag = readSequence(name, &ptr, fd, &sequenceLength,
                          fileOffset, recursionLevel, opt,
                          byteOrder, explicitVR, acceptVRMismatch,
                          fileFlag, remainOpenFlag,
                          convertFlag, object, &scannedSequenceLength,
                          &e, &elementItem);
      if (flag != DCM_NORMAL)
        goto abort;
      if (size != (long) DCM_UNSPECIFIEDLENGTH)
        size -= scannedSequenceLength;
      if (scannedLength != NULL)
        *scannedLength += scannedSequenceLength;

    } else {

      flag = readData(name, &ptr, fd, &size, fileOffset, knownLength,
                      byteOrder, explicitVR, acceptVRMismatch, fileFlag,
                      remainOpenFlag, convertFlag,
                      object, scannedLength, &e, &elementItem);
      if (flag != DCM_NORMAL)
        goto abort;
    }
    computeVM(object, &elementItem->element);

    cond = checkAttributeOrder(&e,
                               &lastGroup,
                               &lastElement,
                               allowRepeatElements);
    if (cond != DCM_NORMAL) {
      if (cond == DCM_REPEATEDELEMENT) {
        CTN_FREE(elementItem);
        continue;
      } else {
        CTN_FREE(elementItem);
        return cond;
      }
    }

    cond = handleGroupItem(object, &groupItem, DCM_TAG_GROUP(e.tag));
    if (cond != DCM_NORMAL)
      /* goto abort; ASG */ return cond;

    if (DCM_TAG_ELEMENT(e.tag) != 0x0000) {
      groupItem->baseLength += 8 + elementItem->paddedDataLength;
      if (elementItem->element.representation == DCM_OB ||
          elementItem->element.representation == DCM_OW ||
          elementItem->element.representation == DCM_SQ) {
        groupItem->longVRAttributes++;
        (*object)->longVRAttributes++;
      }
    }
    if ((DCM_TAG_ELEMENT(e.tag) == 0x0000) \
        && ((*object)->groupLengthFlag == FALSE)) {
      CTN_FREE(elementItem);
    } else {
      cond = LST_Enqueue(&groupItem->elementList, elementItem);
      if (cond != LST_NORMAL) {
        (void) DCM_CloseObject(callerObject);
        return COND_PushCondition(DCM_LISTFAILURE,
                                  DCM_Message(DCM_LISTFAILURE),
                                  "readFile");
      }
      cond = updateObjectType(object, &elementItem->element); /* repair */

      cond = updateSpecialElements(object, elementItem);  /* repair */
    }

    if (size == 0)
      done = TRUE;

    if (part10Flag) {
      if ((*object)->objectSize == \
          (DCM_PREAMBLELENGTH + 4 + 12 + (*object)->metaHeaderLength)) {
        opt &= ~DCM_ORDERMASK;
        opt |= (*object)->dataOptions & DCM_ORDERMASK;
        explicitVR = FALSE;
        switch (opt & DCM_ORDERMASK) {
        case DCM_ORDERNATIVE:
          byteOrder = NATIVE_ORDER;
          break;
        case DCM_ORDERLITTLEENDIAN:
          byteOrder = LITTLE_ORDER;
          break;
        case DCM_EXPLICITLITTLEENDIAN:
          byteOrder = LITTLE_ORDER;
          explicitVR = TRUE;
          break;
        case DCM_ORDERBIGENDIAN:
          byteOrder = BIG_ORDER;
          break;
        case DCM_EXPLICITBIGENDIAN:
          byteOrder = BIG_ORDER;
          explicitVR = TRUE;
          break;
        default:
          byteOrder = LITTLE_ORDER;
          explicitVR = TRUE;
          break;
        }
      }
    }
  }

#ifdef SMM
#endif

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&(*object)->groupList);
  if (groupItem != NULL) {
    (void) LST_Position(&(*object)->groupList, groupItem);
    while (groupItem != NULL) {
      elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
      if (elementItem != NULL) {
        if (DCM_TAG_ELEMENT(elementItem->element.tag) == 0x0000) {
          *elementItem->element.d.ul = groupItem->baseLength;
        }
      }
      groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
    }
  }
  return DCM_NORMAL;

abort:
  return flag;
}

/* locateElement
**
** Purpose:
**  Locate the DICOM element with the specified tag number in the given
**  DICOM object
**
** Parameter Dictionary:
**  obj   Handle to the DICOM object to be searched
**  tag   Tag number of the element to be searched
**
** Return Values:
**  Pointer to the element if found else NULL
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/
static PRV_ELEMENT_ITEM *
locateElement(PRIVATE_OBJECT ** obj, DCM_TAG tag) {
  PRV_GROUP_ITEM
  * groupItem;
  PRV_ELEMENT_ITEM
  * elementItem;
  CTNBOOLEAN
  found = FALSE;

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&(*obj)->groupList);
  if (groupItem == NULL)
    return NULL;

  (void) LST_Position(&(*obj)->groupList, groupItem);
  while (groupItem != NULL) {
    if (groupItem->group == DCM_TAG_GROUP(tag))
      break;

    groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*obj)->groupList);
  }
  if (groupItem == NULL)
    return NULL;

  elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
  if (elementItem == NULL)
    return NULL;

  (void) LST_Position(&groupItem->elementList, elementItem);
  while (!found && (elementItem != NULL)) {
    if (elementItem->element.tag == tag) {
      found = TRUE;
    } else
      elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
  }
  if (found)
    return elementItem;
  else
    return NULL;
}

/* computeVM
**
** Purpose:
**  Compute the multiplicity of the specified element in the DICOM
**  object
**
** Parameter Dictionary:
**  object    Handle to the DICOM object
**  element   Element whose value multiplicity is to be found out.
**
** Return Values:
**  None
**
** Notes:
**
** Algorithm:
**  Description of the algorithm (optional) and any other notes.
*/
static void
computeVM(PRIVATE_OBJECT ** object, DCM_ELEMENT * element) {
  char
  *c;
  int
  i;

  switch (element->representation) {
  case DCM_AE:    /* Application Entity */
  case DCM_AS:    /* Age string */
  case DCM_CS:    /* Control string */
  case DCM_DA:    /* Date */
  case DCM_DS:    /* Decimal string */
  case DCM_DT:    /* Date/Time */
  case DCM_IS:    /* Integer string */
  case DCM_LO:    /* Long string */
  case DCM_PN:    /* Person Name */
  case DCM_SH:    /* Short string */
  case DCM_TM:    /* Time */
  case DCM_UI:    /* Unique identifier (UID) */
  case DCM_UT:    /* Unlimited text */
    element->multiplicity = 1;
    c = element->d.string;
    for (i = 0; i < (int) element->length; i++)
      if (*c++ == '\\')
        element->multiplicity++;
    break;

  case DCM_FD:    /* Floating double */
    element->multiplicity = element->length / 8;
    break;
  case DCM_AT:    /* Attribute tag */
  case DCM_FL:    /* Float */
  case DCM_SL:    /* Signed long */
  case DCM_UL:    /* Unsigned long */
    element->multiplicity = element->length / 4;
    break;
  case DCM_SS:    /* Signed short */
  case DCM_US:    /* Unsigned short */
    element->multiplicity = element->length / 2;
    break;
  case DCM_LT:    /* Long text */
  case DCM_OT:    /* Other binary value */
  case DCM_SQ:    /* Sequence of items */
  case DCM_ST:    /* Short text */
    /*case DCM_UNKNOWN:*/
  case DCM_UN:
  case DCM_RET:
  case DCM_CTX:
  case DCM_DD:    /* Data set */
  default:
    element->multiplicity = 1;
    break;
  }
}
/* ctxSensitiveLookup
**
** Purpose:
**  Lookup representation of elements that are context sensitive
**
** Parameter Dictionary:
**  object    Handle to the DICOM object containing this element.
**  element   Element who representation is to be determined.
**
** Return Values:
**  None
**
** Notes:
**
*/
static void
ctxSensitiveLookup(PRIVATE_OBJECT ** object, DCM_ELEMENT * element) {
  switch (element->tag) {
  case DCM_IMGSMALLESTIMAGEPIXELVALUE:
  case DCM_IMGLARGESTIMAGEPIXELVALUE:
  case DCM_IMGSMALLESTPIXELVALUESERIES:
  case DCM_IMGLARGESTPIXELVALUESERIES:
  case DCM_IMGSMALLESTIMAGEPIXELVALUEPLANE:
  case DCM_IMGLARGESTIMAGEPIXELVALUEPLANE:
  case DCM_IMGLUTDESCRIPTOR:
  case DCM_IMGLUTDATA:
  case DCM_IMGLOOKUPDATARED:
  case DCM_IMGLOOKUPDATAGREEN:
  case DCM_IMGLOOKUPDATABLUE:
    if ((*object)->pixelRepresentation == 0x0000)
      element->representation = DCM_US;
    else if ((*object)->pixelRepresentation == 0x0001)
      element->representation = DCM_SS;
    else
      element->representation = DCM_US;
    break;
  case DCM_MAKETAG(0x003a, 0x1000):
          if (strcmp((*object)->waveformDataVR, "SS") == 0)
            element->representation = DCM_SS;
    break;

  default:
    break;
  }
}

static CONDITION
copyData(PRIVATE_OBJECT ** object, PRV_ELEMENT_ITEM * from,
         DCM_ELEMENT * to, U32 * rtnLength) {
  unsigned char *p = NULL;
  U32 l;
  int nBytes;
  CONDITION cond;

  if (from->element.representation == DCM_SQ)
    return COND_PushCondition(DCM_CANNOTGETSEQUENCEVALUE,
                              DCM_Message(DCM_CANNOTGETSEQUENCEVALUE),
                              from->element.tag, "copyData (DCM internal)");

  l = MIN(from->element.length, to->length);
  if (rtnLength != NULL)
    *rtnLength = l;

  if (from->element.d.ot == NULL) {
    if ((*object)->fd != -1) {
      (void) lseek((*object)->fd,
                   from->dataOffset + (long) p, SEEK_SET);
      nBytes = read((*object)->fd, to->d.ot, (int) l);
    } else {
      (*object)->sk((*object)->userCtx,
                    (long) (from->dataOffset + (long) p), SEEK_SET);
      cond = (*object)->rd((*object)->userCtx, to->d.ot, (long) l, &nBytes);
    }
    if (nBytes != (int) l) {
      return COND_PushCondition(DCM_FILEACCESSERROR,
                                DCM_Message(DCM_FILEACCESSERROR),
                                (*object)->fileName,
                                "copyData (DCM internal)");
    }
#ifdef LITTLE_ENDIAN_ARCHITECTURE
    if (from->element.representation == DCM_AT) {
      DCM_ELEMENT e;
      e = from->element;
      e.length = l;
      e.d.ot = to->d.ot;
      swapATGroupElement(&e);
    }
#endif
    if (from->byteOrder == BYTEORDER_REVERSE) {
      DCM_ELEMENT e;
      e = from->element;
      e.length = l;
      e.d.ot = to->d.ot;
      swapInPlace(object, &e);
    }
  } else {
    unsigned char *q;
    q = (unsigned char *) from->element.d.ot +
        (long) p;
    (void) memcpy(to->d.ot, q, l);
  }
  p += l;
  if ((long) p == from->element.length)
    return DCM_NORMAL;
  else
    return DCM_GETINCOMPLETE;
}

static CONDITION
readLengthToEnd(int fd, const char *fileName,
                unsigned long opt, U32 * lengthToEnd) {
  unsigned char buf[24];
  DCM_OBJECT *obj;
  CONDITION cond;
  DCM_ELEMENT e = {DCM_MAKETAG(0x0008, 0x0001), DCM_UL, "", 1, 4, {NULL}};
  void *ctx = NULL;
  U32 rtnLength = 0;

  if (read(fd, buf, 24) != 24)
    return COND_PushCondition(DCM_FILEACCESSERROR,
                              DCM_Message(DCM_FILEACCESSERROR), fileName,
                              "(DCM)readLengthToEnd");

  cond = DCM_ImportStream(buf, 24, opt, &obj);
  if (cond != DCM_NORMAL)
    return cond;

  e.d.ul = lengthToEnd;
  cond = DCM_GetElementValue(&obj, &e, &rtnLength, &ctx);

  (void) DCM_CloseObject(&obj);

  return cond;
}

#ifdef LITTLE_ENDIAN_ARCHITECTURE
static void
swapATGroupElement(DCM_ELEMENT * e) {
  U32
  length;
  unsigned short
  tmp,
  *us;

  length = e->length;
  us = e->d.us;
  while (length >= 4) {
    tmp = us[0];
    us[0] = us[1];
    us[1] = tmp;
    us += 2;
    length -= 4;
  }
}
#endif

static void
dumpSS(short *ss, long vm) {
  long index = 0;
  while (index < vm) {
    printf("%7d ", *(ss++));
    if ((++index) % 8 == 0)
      printf("\n");
  }
  printf("\n");
}

static void
dumpSL(S32 * sl, long vm) {
  long index = 0;
  while (index < vm) {
    printf("%7d ", *(sl++));
    if ((++index) % 8 == 0)
      printf("\n");
  }
  printf("\n");
}

static void
dumpUS(unsigned short *us, long vm) {
  long index = 0;
  while (index < vm) {
    printf("%7d ", *(us++));
    if ((++index) % 8 == 0)
      printf("\n");
  }
  printf("\n");
}
static void
dumpUL(U32 * ul, long vm) {
  long index = 0;
  while (index < vm) {
    printf("%7u ", *(ul++));
    if ((++index) % 8 == 0)
      printf("\n");
  }
  printf("\n");
}
static void
dumpOB(unsigned char* c, long vm) {
  long index = 0;
  while (index < vm) {
    printf("%02x ", *(c++));
    if ((++index) % 8 == 0)
      printf("\n");
  }
  printf("\n");
}

static void
dumpBinaryData(void *d, DCM_VALUEREPRESENTATION vr, long vm,
               long vmLimit) {
  vm = (vm < vmLimit) ? vm : vmLimit;

  if (vm <= 1)
    return;


  switch (vr) {
  case DCM_SL:
    dumpSL((S32 *) d, vm);
    break;
  case DCM_UL:
    dumpUL((U32 *) d, vm);
    break;
  case DCM_SS:
    dumpSS((short *) d, vm);
    break;
  case DCM_US:
    dumpUS((unsigned short *) d, vm);
    break;
  case DCM_OB:
  case DCM_UN:
    dumpOB((unsigned char*) d, vm);
    break;
  default:
    break;
  }
}

static void
compareGroup(PRV_GROUP_ITEM * g1, PRV_GROUP_ITEM * g2,
             void (*callback) (const DCM_ELEMENT * e1,
                               const DCM_ELEMENT * e2,
                               void *ctx),
             void *ctx) {
  PRV_ELEMENT_ITEM *e1 = NULL,
                         *e2 = NULL;

  if (g1 != NULL) {
    e1 = (PRV_ELEMENT_ITEM*)LST_Head(&g1->elementList);
    if (e1 != NULL)
      LST_Position(&g1->elementList, e1);
  }
  if (g2 != NULL) {
    e2 = (PRV_ELEMENT_ITEM*)LST_Head(&g2->elementList);
    if (e2 != NULL)
      LST_Position(&g2->elementList, e2);
  }
  while (e1 != NULL) {
    if (e2 == NULL) {
      callback(&e1->element, NULL, ctx);
      e1 = (PRV_ELEMENT_ITEM*)LST_Next(&g1->elementList);
    } else if (e1->element.tag == e2->element.tag) {
      callback(&e1->element, &e2->element, ctx);
      e1 = (PRV_ELEMENT_ITEM*)LST_Next(&g1->elementList);
      e2 = (PRV_ELEMENT_ITEM*)LST_Next(&g2->elementList);
    } else if (e1->element.tag < e2->element.tag) {
      callback(&e1->element, NULL, ctx);
      e1 = (PRV_ELEMENT_ITEM*)LST_Next(&g1->elementList);
    } else {
      callback(NULL, &e2->element, ctx);
      e2 = (PRV_ELEMENT_ITEM*)LST_Next(&g2->elementList);
    }
  }

  while (e2 != NULL) {
    callback(NULL, &e2->element, ctx);
    e2 = (PRV_ELEMENT_ITEM*)LST_Next(&g2->elementList);
  }
}

static void
remapFileName(const char *name, char *mapName) {
  char c;

  while ((c = *name++) != '\0') {
    if (c == '\\')
      *mapName++ = '/';
    else if (isupper(c))
      *mapName++ = tolower(c);
    else
      *mapName++ = c;
  }
  *mapName = '\0';
}

static void
copySequence(PRIVATE_OBJECT ** dstObj, DCM_ELEMENT * e) {
  LST_HEAD *lst;
  DCM_SEQUENCE_ITEM *sqItem=NULL;
  DCM_ELEMENT newElement;

  lst = LST_Create();
  if (e->d.sq != NULL) {
    sqItem = (DCM_SEQUENCE_ITEM*)LST_Head(&e->d.sq);
    (void) LST_Position(&e->d.sq, sqItem);
  }
  while (sqItem != NULL) {
    DCM_OBJECT *copy;
    DCM_SEQUENCE_ITEM *copyItem;

    DCM_CopyObject(&sqItem->object, &copy);
    copyItem = (DCM_SEQUENCE_ITEM *)malloc(sizeof(*copyItem));
    copyItem->object = copy;
    (void) LST_Enqueue(&lst, copyItem);

    sqItem = (DCM_SEQUENCE_ITEM*)LST_Next(&e->d.sq);
  }

  memset(&newElement, 0, sizeof(newElement));
  newElement.tag = e->tag;
  newElement.representation = e->representation;
  newElement.d.sq = lst;
  DCM_AddSequenceElement((DCM_OBJECT **) dstObj, &newElement);
}

/* Restart public functions */

CONDITION
DCM_GetCompressedValue(DCM_OBJECT ** callerObject, DCM_TAG tag, void *buf,
                       size_t bufSize, DCM_GET_COMPRESSED_CALLBACK* callback,
                       void *ctx) {
  PRIVATE_OBJECT
  ** object;
  PRV_ELEMENT_ITEM
  * elementItem;
  S32 nBytes;
  S32 toRead;
  CONDITION cond;
  // int doneFlag = 0;
  size_t elementLength;
  unsigned char *ptr;
  U32 size = 0;
  off_t fileOffset = 0;
  unsigned long opt = 0;
  int byteOrder;
  int explicitVR;
  CTNBOOLEAN acceptVRMismatch = FALSE;
  DCM_ELEMENT e;
  U32 sequenceLength = 0;
  CONDITION flag;
  int index = 0;
  CTNBOOLEAN firstBuffer = TRUE;
  U32 *offsetBuffer = NULL;
  U32 offsetBufferCount = 0;
  U32 streamOffset = 0;
  int startOfFragment = 1;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_GetCompressedValue");
  if (cond != DCM_NORMAL)
    return cond;

  elementItem = locateElement(object, tag);

  if (elementItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_ELEMENT(tag),
                              "DCM_GetEncodedValue");

  elementLength = elementItem->originalDataLength;
  ptr = NULL;     /* Means reading from a file */
  size = DCM_UNSPECIFIEDLENGTH;
  fileOffset = elementItem->dataOffset;

  opt |= (*object)->dataOptions & DCM_ORDERMASK;
  explicitVR = FALSE;
  switch (opt & DCM_ORDERMASK) {
  case DCM_ORDERNATIVE:
    byteOrder = NATIVE_ORDER;
    break;
  case DCM_ORDERLITTLEENDIAN:
    byteOrder = LITTLE_ORDER;
    break;
  case DCM_EXPLICITLITTLEENDIAN:
    byteOrder = LITTLE_ORDER;
    explicitVR = TRUE;
    break;
  case DCM_ORDERBIGENDIAN:
    byteOrder = BIG_ORDER;
    break;
  case DCM_EXPLICITBIGENDIAN:
    byteOrder = BIG_ORDER;
    explicitVR = TRUE;
    break;
  default:
    byteOrder = LITTLE_ORDER;
    explicitVR = TRUE;
    break;
  }
  if ((opt & DCM_VRMASK) == DCM_ACCEPTVRMISMATCH)
    acceptVRMismatch = TRUE;

  (void) lseek((*object)->fd, elementItem->dataOffset, SEEK_SET);
  while (elementLength != 0) {
    sequenceLength = 0;
    memset(&e, 0, sizeof(e));
    flag = readGroupElement("", &ptr, (*object)->fd, &size, &fileOffset,
                            FALSE, byteOrder, explicitVR, acceptVRMismatch,
                            object, &sequenceLength, &e);
    if (flag == DCM_STREAMCOMPLETE)
      break;
    else if (flag != DCM_NORMAL)
      return flag;

    flag = readVRLength("", &ptr, (*object)->fd, &size, &fileOffset,
                        FALSE,  /* Known length */
                        byteOrder, explicitVR, acceptVRMismatch, object,
                        &sequenceLength, &e);
    if (flag != DCM_NORMAL)
      return flag;

    elementLength -= sequenceLength + e.length;

    if (firstBuffer) {
      firstBuffer = FALSE;
      if (e.length != 0) {
        offsetBuffer = (U32 *)CTN_MALLOC(e.length);
        offsetBufferCount = e.length / sizeof(U32);
        if (offsetBuffer == NULL)
          exit(1);  /* repair */
        nBytes = read((*object)->fd, offsetBuffer, e.length);
        if (nBytes != (ssize_t)e.length) {
          exit(1);  /* repair */
        }
        if (byteOrder == BYTEORDER_REVERSE) {
          DCM_ELEMENT offsetBufferElement;
          memset(&offsetBufferElement, 0, sizeof(DCM_ELEMENT));
          offsetBufferElement.length = e.length;
          offsetBufferElement.d.ul = offsetBuffer;
          offsetBufferElement.representation = DCM_UL;
          swapInPlace(object, &offsetBufferElement);
        }
        callback(offsetBuffer, e.length, index, 1, 0, 1, ctx);
        streamOffset = 0;
      } else {
        streamOffset = 0xffffffff;
      }
    } else {
      U32 l = e.length;
      U32 j;
      int lastIndex;

      lastIndex = index;
      for (j = 0; j < offsetBufferCount; j++) {
        if (streamOffset == offsetBuffer[j])
          index = j + 1;
      }
      startOfFragment = 1;
      while (l != 0) {
        toRead = MIN(bufSize, l);
        nBytes = read((*object)->fd, buf, toRead);
        if (nBytes != toRead) {
          exit(1);  /* repair */
        }
        callback(buf, toRead, index,
                 (index != lastIndex) ? 1 : 0,
                 0, startOfFragment, ctx);
        l -= toRead;
        lastIndex = index;  /* Guarantee first flag is off */
        startOfFragment = 0;
      }
      streamOffset += sequenceLength + e.length;
    }
    fileOffset += e.length;
    index++;
  }
  callback(buf, 0, index, 0, 1, 1, ctx);
  return DCM_NORMAL;
}

CONDITION
DCM_PrintSequenceList(DCM_OBJECT ** object, DCM_TAG tag) {
  PRIVATE_OBJECT **obj,
  *sqObject;
  CONDITION cond = DCM_NORMAL;
  PRV_ELEMENT_ITEM *elementItem;
  LST_HEAD *lst;
  DCM_SEQUENCE_ITEM *sqItem;

  obj = (PRIVATE_OBJECT **) object;
  cond = checkObject(obj, "DCM_PrintSequenceList");
  if (cond != DCM_NORMAL)
    return cond;

  elementItem = locateElement(obj, tag);

  if (elementItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_ELEMENT(tag),
                              "DCM_PrintSequenceList");

  lst = elementItem->element.d.sq;
  sqItem = (DCM_SEQUENCE_ITEM *)LST_Head(&lst);
  (void) LST_Position(&lst, sqItem);
  while (sqItem != NULL) {
    sqObject = (PRIVATE_OBJECT *) sqItem->object;
    printf("size: %6ld offset: %6ld, pixel offset: %6ld\n",
           sqObject->objectSize,
           sqObject->offset,
           sqObject->pixelOffset);
    sqItem = (DCM_SEQUENCE_ITEM *)LST_Next(&lst);
  }
  return cond;
}

CONDITION
DCM_GetSequenceByOffset(DCM_OBJECT ** object,
                        DCM_TAG tag,
                        unsigned long offset,
                        DCM_OBJECT ** rtnObject) {
  PRIVATE_OBJECT **obj,
  *sqObject;
  CONDITION cond;
  PRV_ELEMENT_ITEM *elementItem;
  LST_HEAD *lst;
  DCM_SEQUENCE_ITEM *sqItem;

  obj = (PRIVATE_OBJECT **) object;
  cond = checkObject(obj, "DCM_PrintSequenceList");
  if (cond != DCM_NORMAL)
    return cond;

  elementItem = locateElement(obj, tag);

  if (elementItem == NULL)
    return COND_PushCondition(DCM_ELEMENTNOTFOUND,
                              DCM_Message(DCM_ELEMENTNOTFOUND),
                              DCM_TAG_GROUP(tag),
                              DCM_TAG_ELEMENT(tag),
                              "DCM_PrintSequenceList");

  lst = elementItem->element.d.sq;
  sqItem = (DCM_SEQUENCE_ITEM *)LST_Head(&lst);
  (void) LST_Position(&lst, sqItem);
  while (sqItem != NULL) {
    sqObject = (PRIVATE_OBJECT *) sqItem->object;
    if (sqObject->offset == offset) {
      *rtnObject = sqItem->object;
      return DCM_NORMAL;
    }
    sqItem = (DCM_SEQUENCE_ITEM *)LST_Next(&lst);
  }
  return 0;
}

CONDITION
DCM_CopyObject(DCM_OBJECT ** src, DCM_OBJECT ** dst) {
  PRIVATE_OBJECT **srcObj;
  PRIVATE_OBJECT *dstObj;
  PRV_GROUP_ITEM *groupItem;
  PRV_ELEMENT_ITEM *elementItem;

  if (src == NULL) {
    (void) COND_PushCondition(DCM_NULLADDRESS,
                              DCM_Message(DCM_NULLADDRESS),
                              "DCM_CopyObject");
    return COND_PushCondition(DCM_OBJECTCREATEFAILED,
                              DCM_Message(DCM_OBJECTCREATEFAILED),
                              "DCM_CopyObject");
  }
  dstObj = (PRIVATE_OBJECT *) CTN_MALLOC(sizeof(PRIVATE_OBJECT));
  if (dstObj == NULL) {
    (void) COND_PushCondition(DCM_MALLOCFAILURE,
                              DCM_Message(DCM_MALLOCFAILURE),
                              sizeof(PRIVATE_OBJECT),
                              "DCM_CopyObject");
    *dst = NULL;
    return COND_PushCondition(DCM_OBJECTCREATEFAILED,
                              DCM_Message(DCM_OBJECTCREATEFAILED),
                              "DCM_CopyObject");
  }
  (void) memset(dstObj, 0, sizeof(PRIVATE_OBJECT));
  (void) strcpy(dstObj->keyType, KEY_DCM_OBJECT);

  dstObj->objectType = DCM_OBJECTUNKNOWN;
  dstObj->accessMethod = DCM_MEMORY_ACCESS;
  dstObj->deleteFlag = FALSE;
  dstObj->groupLengthFlag = FALSE;
  dstObj->objectSize = 0;
  dstObj->offset = 0;
  dstObj->pixelSize = 0;
  dstObj->pixelOffset = 0;
  dstObj->pixelBitsAllocated = 0;
  dstObj->pixelRepresentation = 0xffff;
  dstObj->groupCtx = NULL;
  dstObj->elementCtx = NULL;
  dstObj->fd = -1;
  dstObj->fileName[0] = '\0';
  dstObj->preambleFlag = FALSE;
  dstObj->preamble[0] = '\0';
  dstObj->dataOptions = 0;
  dstObj->metaHeaderLength = 0xffffffff;
  dstObj->longVRAttributes = 0;
  dstObj->waveformDataVR[0] = '\0';

  dstObj->groupList = LST_Create();
  if (dstObj->groupList == NULL) {
    CTN_FREE(dstObj);
    *dst = NULL;
    return COND_PushCondition(DCM_LISTFAILURE,
                              DCM_Message(DCM_LISTFAILURE),
                              "DCM_CreateObject");
  }
  srcObj = (PRIVATE_OBJECT **) src;

  groupItem = (PRV_GROUP_ITEM *)LST_Head(&(*srcObj)->groupList);
  if (groupItem != NULL)
    (void) LST_Position(&(*srcObj)->groupList, groupItem);

  while (groupItem != NULL) {
    elementItem = (PRV_ELEMENT_ITEM *)LST_Head(&groupItem->elementList);
    if (elementItem != NULL)
      (void) LST_Position(&groupItem->elementList, elementItem);
    while (elementItem != NULL) {
      if (elementItem->element.representation == DCM_SQ) {
        copySequence(&dstObj, &elementItem->element);
      } else {
        void* pvoid = (void*)& dstObj;
        DCM_AddElement((DCM_OBJECT **) pvoid, &elementItem->element);
      }
      elementItem = (PRV_ELEMENT_ITEM *)LST_Next(&groupItem->elementList);
    }
    groupItem = (PRV_GROUP_ITEM *)LST_Next(&(*srcObj)->groupList);
  }

  *dst = (DCM_OBJECT *) dstObj;
  return DCM_NORMAL;
}

CONDITION
DCM_MergeObject(DCM_OBJECT ** src, DCM_OBJECT ** dst) {
  PRIVATE_OBJECT **srcObj;
  PRIVATE_OBJECT *dstObj;
  PRV_GROUP_ITEM *groupItem;
  PRV_ELEMENT_ITEM *elementItem;

  if (src == NULL) {
    (void) COND_PushCondition(DCM_NULLADDRESS,
                              DCM_Message(DCM_NULLADDRESS),
                              "DCM_MergeObject");
    return COND_PushCondition(DCM_OBJECTCREATEFAILED,
                              DCM_Message(DCM_OBJECTCREATEFAILED),
                              "DCM_MergeObject");
  }
  dstObj = *((PRIVATE_OBJECT **)dst);
  if (dstObj == NULL) {
    (void) COND_PushCondition(DCM_MALLOCFAILURE,
                              DCM_Message(DCM_MALLOCFAILURE),
                              sizeof(PRIVATE_OBJECT),
                              "DCM_MergeObject");
    *dst = NULL;
    return COND_PushCondition(DCM_OBJECTCREATEFAILED,
                              DCM_Message(DCM_OBJECTCREATEFAILED),
                              "DCM_MergeObject");
  }
  srcObj = (PRIVATE_OBJECT **) src;

  groupItem = (PRV_GROUP_ITEM *)LST_Head(&(*srcObj)->groupList);
  if (groupItem != NULL)
    (void) LST_Position(&(*srcObj)->groupList, groupItem);

  while (groupItem != NULL) {
    elementItem = (PRV_ELEMENT_ITEM *)LST_Head(&groupItem->elementList);
    if (elementItem != NULL)
      (void) LST_Position(&groupItem->elementList, elementItem);
    while (elementItem != NULL) {
      if (elementItem->element.representation == DCM_SQ) {
        copySequence(&dstObj, &elementItem->element);
      } else {
        void* pvoid = (void*) & dstObj;
        DCM_AddElement((DCM_OBJECT **) pvoid, &elementItem->element);
      }
      elementItem = (PRV_ELEMENT_ITEM *)LST_Next(&groupItem->elementList);
    }
    groupItem = (PRV_GROUP_ITEM *)LST_Next(&(*srcObj)->groupList);
  }

  /**dst = (DCM_OBJECT *) dstObj;*/
  return DCM_NORMAL;
}


CONDITION
DCM_GetFirstElement(DCM_OBJECT ** callerObject, DCM_ELEMENT** e) {
  PRIVATE_OBJECT** object;
  PRV_GROUP_ITEM* groupItem;
  PRV_ELEMENT_ITEM* elementItem;
  CONDITION cond;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_GetFirstElement");
  if (cond != DCM_NORMAL)
    return cond;

  groupItem = (PRV_GROUP_ITEM*)LST_Head(&(*object)->groupList);

  if (groupItem == NULL) {
    *e = 0;
    return DCM_EMPTYOBJECT;
  }
  (void) LST_Position(&(*object)->groupList, groupItem);
  (*object)->groupCtx = groupItem;

  elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
  (*object)->elementCtx = elementItem;
  if (elementItem == NULL) {
    return DCM_GetNextElement(callerObject, e);
  }

  *e = &elementItem->element;
  return DCM_NORMAL;
}

CONDITION
DCM_GetNextElement(DCM_OBJECT ** callerObject, DCM_ELEMENT** e) {
  PRIVATE_OBJECT** object;
  PRV_GROUP_ITEM* groupItem;
  PRV_ELEMENT_ITEM* elementItem;
  CONDITION cond;

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_GetNextElement");
  if (cond != DCM_NORMAL)
    return cond;

  groupItem = (*object)->groupCtx;
  elementItem = (*object)->elementCtx;

  if (elementItem != 0) {
    (void)LST_Position(&groupItem->elementList, elementItem);
    elementItem = (PRV_ELEMENT_ITEM*)LST_Next(&groupItem->elementList);
  }

  if (elementItem == 0) {
    (void)LST_Position(&(*object)->groupList, groupItem);
    groupItem = (PRV_GROUP_ITEM*)LST_Next(&(*object)->groupList);
    if (groupItem != 0) {
      elementItem = (PRV_ELEMENT_ITEM*)LST_Head(&groupItem->elementList);
    }
  }

  if (groupItem == 0) {
    *e = 0;
    return DCM_GETNEXTELEMENTCOMPLETE;
  }

  (*object)->groupCtx = groupItem;
  (*object)->elementCtx = elementItem;

  if (elementItem == 0)
    return DCM_GetNextElement(callerObject, e);

  *e = &elementItem->element;
  return DCM_NORMAL;
}

CONDITION
DCM_AddFragment(DCM_OBJECT** callerObject, void* fragment, U32 fragmentLength) {
  PRIVATE_OBJECT** object;
  PRV_ELEMENT_ITEM* elementItem;
  PRV_ELEMENT_ITEM* newItem;
  CONDITION cond;
  PRV_GROUP_ITEM *groupItem = 0;
  DCM_FRAGMENT_ITEM* fragmentItem;
  U32 mallocLength;

  if ((fragmentLength & 1) != 0) {
    return COND_PushCondition(DCM_UNEVENFRAGMENTLENGTH,
                              DCM_Message(DCM_UNEVENFRAGMENTLENGTH),
                              fragmentLength, "DCM_AddFragment");
  }

  object = (PRIVATE_OBJECT **) callerObject;
  cond = checkObject(object, "DCM_AddFragment");
  if (cond != DCM_NORMAL)
    return cond;

  cond = findCreateGroup(object, 0x7fe0, &groupItem);
  if (cond != DCM_NORMAL)
    return COND_PushCondition(DCM_INSERTFAILED,
                              DCM_Message(DCM_INSERTFAILED),
                              0x7fe0, 0x0010, "DCM_AddFragment");

  elementItem = locateElement(object, 0x7fe00010);
  if (elementItem == NULL) {
    DCM_ELEMENT e;
    memset(&e, 0, sizeof(e));
    e.tag = DCM_PXLPIXELDATA;
    e.representation = DCM_OB;
    e.multiplicity = 1;
    e.length = 0;
    e.d.fragments = 0;
    cond = newElementItem(&e, FALSE, &newItem);
    if (cond != DCM_NORMAL)
      return cond;
    newItem->element.d.fragments = LST_Create();
    if (newItem->element.d.fragments == NULL) {
      return COND_PushCondition(DCM_LISTFAILURE,
                                DCM_Message(DCM_LISTFAILURE),
                                "DCM_AddFragment");
    }
    cond = insertThisElementItem(object, newItem);
    if (cond != DCM_NORMAL)
      return cond;
  }

  elementItem = locateElement(object, 0x7fe00010);

  mallocLength = sizeof(DCM_FRAGMENT_ITEM) + fragmentLength;
  fragmentItem = (DCM_FRAGMENT_ITEM*)CTN_MALLOC(mallocLength);
  if (fragmentItem == NULL) {
    return COND_PushCondition(DCM_MALLOCFAILURE,
                              DCM_Message(DCM_MALLOCFAILURE), mallocLength,
                              "DCM_AddFragment");
  }

  fragmentItem->fragment = (unsigned char*)fragmentItem +
                           (long)sizeof(DCM_FRAGMENT_ITEM);
  fragmentItem->length = fragmentLength;
  memcpy(fragmentItem->fragment, fragment, fragmentLength);
  elementItem->fragmentFlag = 1;
  LST_Enqueue(&elementItem->element.d.fragments, fragmentItem);

  return DCM_NORMAL;
}
