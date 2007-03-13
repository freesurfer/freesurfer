/* getdelim.h --- Prototype for replacement getdelim function.
   Copyright (C) 2005 Free Software Foundation, Inc.
*/
/* Written by Simon Josefsson. */

/* Get size_t, FILE, ssize_t.  And getdelim, if available.  */
# include <stddef.h>
# include <stdio.h>
# include <sys/types.h>

#if !HAVE_DECL_GETDELIM
ssize_t getdelim (char **lineptr, size_t *n, int delimiter, FILE *stream);
#endif /* !HAVE_GETDELIM */
