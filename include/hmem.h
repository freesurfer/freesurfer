/**
 * @file  hmem.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.6 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef HMEM_H
#define HMEM_H

#ifndef _MSDOS
#define hmemset   memset
#define hmemcpy   memcpy
#define hcalloc   calloc
#else

#define hmemset(buf, val, count)   huge_memset((huge char *)buf, (char)val, (long)count)
#define hcalloc(n,size)            huge_alloc((hsize_t)n*size, 1)
#define hmemcpy(Dst, Src, Count)   huge_memcpy((huge char *)Dst, (huge char *)Src, (long)Count)

#endif

#endif
