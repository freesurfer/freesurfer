/**
 * @brief define BYTE_ORDER, LITTLE_ENDIAN, BIG_ENDIAN and BYTE_ORDER
 *
 * SGI  /usr/include/sys/endian.h (decides on _MIPSEB or _MIPSEL)
 * Mac /usr/include/machine/endian.h (decides on __ppc__ or __i386__)
 * Linux  /usr/include/endian.h -> /usr/include/bits/endian.h
 * Solaris /usr/include/sys/isa_defs.h
 *
 * BSD defines non-underscore variables
 */
/*
 * Original Author: Yasunari Tosa
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

/* See if it is set in the imake config first */
#ifndef mghendian_h
#define mghendian_h

/* allow override */
#ifndef BYTE_ORDER

/////////////Linux////////////////////////////
#ifdef __linux__
#include <endian.h>

#ifndef BYTE_ORDER
#define BYTE_ORDER __BYTE_ORDER
#endif

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN __LITTLE_ENDIAN
#endif

#ifndef BIG_ENDIAN
#define BIG_ENDIAN __BIG_ENDIAN
#endif

#endif

/////////////Windows Cygwin////////////////////////////
#ifdef Windows_NT
#include <endian.h>

#ifndef BYTE_ORDER
#define BYTE_ORDER __BYTE_ORDER
#endif

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN __LITTLE_ENDIAN
#endif

#ifndef BIG_ENDIAN
#define BIG_ENDIAN __BIG_ENDIAN
#endif

#endif

////////////MacOS X and BSD ////////////////////////////
#if defined(__APPLE__) || defined(__NetBSD__) || defined(__OpenBSD__)
#include <machine/endian.h>
#endif

////////////Solaris 2.5.1//////////////////////
#ifdef sun

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN 1234
#endif

#ifndef BIG_ENDIAN
#define BIG_ENDIAN    4321
#endif

#include <sys/isa_defs.h>
/* only defines one of _LITTLE_ENDIAN or _BIG_ENDIAN */
#ifdef _LITTLE_ENDIAN
#define BYTE_ORDER LITTLE_ENDIAN
#endif

#ifdef _BIG_ENDIAN
#define BYTE_ORDER BIG_ENDIAN
#endif

#endif

/////////////IRIX  ////////////////////////////
#if defined(__sgi) || defined(Mips)
#include <sys/endian.h>
#endif

///////////////////////////////////////////////
#ifndef BYTE_ORDER
#error "Unknown OS to mghendian.h.   Please report"
#endif

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN 1234
#endif

#ifndef BIG_ENDIAN
#define BIG_ENDIAN    4321
#endif

///////////////////////////////////////////////////
#endif /* BYTE_ORDER */
///////////////////////////////////////////////////

/* final check   bomb if not */
#if !defined(BYTE_ORDER) || \
    (BYTE_ORDER != BIG_ENDIAN && BYTE_ORDER != LITTLE_ENDIAN)
/*  && BYTE_ORDER != PDP_ENDIAN) */
#error "Undefined or invalid BYTE_ORDER";
#endif

#endif /* mghendian.h */
