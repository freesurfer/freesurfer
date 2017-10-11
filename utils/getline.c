/**
   @file getline.c
   @brief Implementation of replacement getline function.

   Copyright (C) 2005 Free Software Foundation, Inc.
*/

/* Written by Simon Josefsson. */

#if defined(Darwin) || defined(SunOS)
// getline does not exist on Mac OS X or Solaris

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include "getdelim.h"
#include "getline.h"

ssize_t getline(char **lineptr, size_t *n, FILE *stream) { return getdelim(lineptr, n, '\n', stream); }

#endif  // Darwin or SunOS
