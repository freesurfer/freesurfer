/**
 * @file  colortab.c
 * @brief color table utilities
 *
 * An entry in a color table has:
 *   1. string name
 *   2. rgb (in both int and float)
 * An annotation is is an int packed with the int values of the
 * rgb in the first 3 bytes of the annotation int.
 */
/*
 * Original Authors: Kevin Teich, Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/03/28 19:15:30 $
 *    $Revision: 1.22 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "error.h"
#include "colortab.h"
#include "fio.h"
#include "proto.h"

#define CTAB_VERSION_TO_WRITE 2

/* Different binary i/o versions. */
static int         CTABwriteIntoBinaryV1(COLOR_TABLE *ct, FILE *fp);
static int         CTABwriteIntoBinaryV2(COLOR_TABLE *ct, FILE *fp);
static COLOR_TABLE *CTABreadFromBinaryV1(FILE *fp, int nentries);
static COLOR_TABLE *CTABreadFromBinaryV2(FILE *fp);


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *CTABreadASCII(char *fname)
{
  COLOR_TABLE *ct;
  char        line[STRLEN], *cp;
  int         max_structure;
  FILE        *fp;
  int         structure;
  char        name[STRLEN];
  int         r, g, b, t;
  int         line_num;

  /* Try to open the file. */
  fp = fopen(fname, "r");
  if (fp == NULL)
    ErrorReturn
      (NULL,
       (ERROR_NOFILE, "CTABreadASCII(%s): could not open file", fname));

  /* Scan through the file and see what our max entry number is. */
  max_structure = -1;
  while ((cp = fgetl(line, STRLEN, fp)) != NULL)
  {
    /* See if this line is in the right format. If not, it's
    probably just a comment and we can ignore it. */
    if (sscanf(line, "%d %s %d %d %d %d",
               &structure, name, &r, &g, &b, &t) == 6)
    {
      if (structure > max_structure)
        max_structure = structure;
    }
  }

  /* If no structures, bail. */
  if (max_structure <= 0)
  {
    fclose(fp);
    ErrorReturn
      (NULL,
       (ERROR_NOFILE, "CTABreadASCII(%s): badly formed file", fname));
  }

  /* Allocate our table. */
  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (ct == NULL)
    ErrorReturn
      (NULL, 
       (ERROR_NO_MEMORY, 
        "CTABreadASCII(%s): could not allocate table", fname));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = max_structure + 1;
  ct->entries = (COLOR_TABLE_ENTRY**)
                calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY*));
  if (ct->entries == NULL)
    ErrorReturn
      (NULL, 
       (ERROR_NO_MEMORY, 
        "CTABreadASCII(%s): could not allocate %d entries", 
        fname, ct->nentries));

  /* Copy in the file name. */
  strncpy(ct->fname, fname, sizeof(ct->fname));

  /* We'll write this version if we write to binary. */
  ct->version = CTAB_VERSION_TO_WRITE;

  /* Rewind the file and go through it again. For each entry we find,
     allocate a CTE. This will leave the items in the array for which
     we don't have entries NULL. */
  line_num = 1;
  rewind(fp);
  while ((cp = fgets(line, STRLEN, fp)) != NULL)
  {
    if (sscanf (line, "%d %s %d %d %d %d",
                &structure, name, &r, &g, &b, &t) == 6)
    {

      /* If this entry already exists, there's a duplicate entry
         in the file. Warn, but then continue on.*/
      if (ct->entries[structure] != NULL)
      {
        printf ("CTABreadASCII(%s): Line %d: Duplicate structure "
                "index %d, was %s %d %d %d %d\n",
                fname, line_num, structure,
                ct->entries[structure]->name,
                ct->entries[structure]->ri,
                ct->entries[structure]->gi,
                ct->entries[structure]->bi,
                ct->entries[structure]->ai);
      }
      else
      {
        /* Try to create a new entry.*/
        ct->entries[structure] = (CTE*) malloc(sizeof(CTE));
        if (NULL == ct->entries[structure])
        {
          fclose(fp);
          CTABfree(&ct);
          ErrorReturn
            (NULL, 
             (ERROR_NO_MEMORY, 
              "CTABreadASCII(%s): could not allocate entry for structure %d", 
              fname, structure));
        }

        /* Fill out the entry. */
        strncpy (ct->entries[structure]->name, name,
                 sizeof(ct->entries[structure]->name));
        ct->entries[structure]->ri = r;
        ct->entries[structure]->gi = g;
        ct->entries[structure]->bi = b;
        ct->entries[structure]->ai = (255-t); /* alpha = 255-trans */

        /* Now calculate the float versions. */
        ct->entries[structure]->rf =
          (float)ct->entries[structure]->ri / 255.0;
        ct->entries[structure]->gf =
          (float)ct->entries[structure]->gi / 255.0;
        ct->entries[structure]->bf =
          (float)ct->entries[structure]->bi / 255.0;
        ct->entries[structure]->af =
          (float)ct->entries[structure]->ai / 255.0;
      }
    }
    line_num++;
  }

  fclose(fp);

  /* Return the new color table. */
  return(ct);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABfree(COLOR_TABLE **pct)
{
  int i;
  COLOR_TABLE *ct;

  if (NULL==pct)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "CTABfree: pct was NULL"));
  if (NULL==*pct)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "CTABfree: *pct was NULL"));

  ct = *pct;

  /* Free all the entries. */
  for (i = 0; i < ct->nentries; i++)
    if (NULL != ct->entries[i])
      free (ct->entries[i]);

  free (ct->entries);
  free (ct);

  /* Set argument to null */
  *pct = NULL;

  return(NO_ERROR);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *CTABdeepCopy(COLOR_TABLE *ct)
{
  COLOR_TABLE* copy;
  int structure;

  if (NULL==ct)
    ErrorReturn(NULL,(ERROR_BADPARM, "CTABdeepCopy: ct was NULL"));

  /* Make a copy of the table. Allocate our table. */
  copy = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (copy == NULL)
    ErrorReturn(NULL,
                (ERROR_NO_MEMORY, "CTABdeepCopy: could not allocate table"));

  /* Allocate the array of NULL CTE ptrs. */
  copy->nentries = ct->nentries;
  copy->entries = (COLOR_TABLE_ENTRY**)
                  calloc(copy->nentries, sizeof(COLOR_TABLE_ENTRY*));
  if (copy->entries == NULL)
    ErrorReturn(NULL,(ERROR_NO_MEMORY, "CTABdeepCopy: could not allocate "
                      "%d entries", copy->nentries));

  /* Copy in the file name and verion. */
  strncpy(copy->fname, ct->fname, sizeof(copy->fname));
  copy->version = ct->version;

  /* Go through the table. For each entry we find, allocate a CTE in
     the copy.*/
  for (structure = 0; structure < ct->nentries; structure++)
  {
    if (ct->entries[structure] != NULL)
    {
      /* Try to create a new entry.*/
      copy->entries[structure] = (CTE*) malloc(sizeof(CTE));
      if (NULL == copy->entries[structure])
      {
        CTABfree(&copy);
        ErrorReturn(NULL,(ERROR_NO_MEMORY, "CTABdeepCopy: could not "
                          "allocate entry for structure %d", structure));
      }

      /* Copy the entry. */
      strncpy (copy->entries[structure]->name,
               ct->entries[structure]->name,
               sizeof(copy->entries[structure]->name));
      copy->entries[structure]->ri = ct->entries[structure]->ri;
      copy->entries[structure]->gi = ct->entries[structure]->gi;
      copy->entries[structure]->bi = ct->entries[structure]->bi;
      copy->entries[structure]->ai = ct->entries[structure]->ai;
      copy->entries[structure]->rf = ct->entries[structure]->rf;
      copy->entries[structure]->gf = ct->entries[structure]->gf;
      copy->entries[structure]->bf = ct->entries[structure]->bf;
      copy->entries[structure]->af = ct->entries[structure]->af;
    }
  }

  /* Return the new copy. */
  return copy;
}



/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *CTABreadFromBinary(FILE *fp)
{
  COLOR_TABLE   *ct;
  int           version;

  if (NULL==fp)
    ErrorReturn(NULL,
                (ERROR_BADPARM, "CTABreadFromBinary: fp was NULL"));

  /* Our goal here is to see what vesion we're reading/writing. Look
     at the first int in the stream. If it's > 0, it's the old format,
     and this int is the number of bins to have. If it's < 0, it's the
     new format, and is the negative version number. */
  version = freadInt(fp);

  ct = NULL;
  if (version > 0 )
  {
    /* With v1, we pass in the "version" number we just got, as it's
    the number of entries. */
    ct = CTABreadFromBinaryV1 (fp, version);
  }
  else
  {
    /* Take the negative to get the real version number. */
    version = -version;

    /* Read the right version. */
    if (version == 2)
    {
      ct = CTABreadFromBinaryV2 (fp);
    }
    else
    {
      /* Unsupported version. */
      ErrorReturn(NULL,
                  (ERROR_BADFILE, "CTABreadFromBinary: unknown version"));
    }
  }

  return(ct);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABwriteIntoBinary(COLOR_TABLE *ct, FILE *fp)
{
  int result;

  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABwriteIntoBinary: ct was NULL"));
  if (NULL==fp)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABwriteIntoBinary: fp was NULL"));

  /* Just switch on the version we want to write. */
  switch (ct->version)
  {
  case 1:
    result = CTABwriteIntoBinaryV1 (ct, fp);
    break;
  case 2:
    result = CTABwriteIntoBinaryV2 (ct, fp);
    break;
  default:
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABwriteIntoBinary: unknown version"));
  }

  return(result);
}



/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *CTABalloc(int nentries)
{
  COLOR_TABLE        *ct;
  int                structure;

  if (nentries < 0)
    ErrorReturn(NULL,
                (ERROR_BADPARM, "CTABalloc: nentries was %d", nentries));

  /* Allocate our table. */
  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (ct == NULL)
    ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABalloc(%d): could not "
                       "allocate table", nentries));


  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = nentries;
  ct->entries = (COLOR_TABLE_ENTRY**)
                calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY*));
  if (ct->entries == NULL)
  {
    CTABfree (&ct);
    ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABalloc: could not "
                       "allocate %d entries", ct->nentries));
  }

  /* Copy in default name. */
  strcpy (ct->fname, "none");

  /* We'll write this version if we write to binary. */
  ct->version = CTAB_VERSION_TO_WRITE;

  /* Init all entries to random colors. */
  for (structure = 0; structure < ct->nentries; structure++)
  {

    /* Try to create a new entry.*/
    ct->entries[structure] = (CTE*) malloc(sizeof(CTE));
    if (NULL == ct->entries[structure])
    {
      CTABfree(&ct);
      ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABalloc: could not allocate "
                         "entry for structure %d", structure));
    }

    /* Random colors. */
    ct->entries[structure]->ri = nint(randomNumber(0, 255));
    ct->entries[structure]->gi = nint(randomNumber(0, 255));
    ct->entries[structure]->bi = nint(randomNumber(0, 255));
    ct->entries[structure]->ai = 255;

    /* Now calculate the float versions. */
    ct->entries[structure]->rf =
      (float)ct->entries[structure]->ri / 255.0;
    ct->entries[structure]->gf =
      (float)ct->entries[structure]->gi / 255.0;
    ct->entries[structure]->bf =
      (float)ct->entries[structure]->bi / 255.0;
    ct->entries[structure]->af =
      (float)ct->entries[structure]->ai / 255.0;

    /* Make a fake name. */
    sprintf (ct->entries[structure]->name, "cluster%d", structure);
  }

  return(ct);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *
CTABreadFromBinaryV1(FILE *fp, int nentries)
{
  COLOR_TABLE        *ct;
  int                structure, len;
  char               *name;
  int                t;

  if (nentries < 0)
    ErrorReturn
      (NULL,
       (ERROR_BADPARM, "CTABreadFromBinaryV1: nentries was %d", nentries));

  /* Allocate our table. */
  ct = (COLOR_TABLE *)calloc (1, sizeof(COLOR_TABLE));
  if (ct == NULL)
    ErrorReturn
      (NULL,
       (ERROR_NO_MEMORY, "CTABreadFromBinaryV1: could not allocate table"));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = nentries;
  ct->entries = (COLOR_TABLE_ENTRY**)
                calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY*));
  if (ct->entries == NULL)
  {
    CTABfree (&ct);
    ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABreadFromBinaryV1: could not "
                       "allocate %d entries", ct->nentries));
  }

  /* We'll write the same version to binary. */
  ct->version = 1;

  /* Read the file name. Read it into a temp buffer first to avoid
     overflow. */
  len = freadInt (fp);
  if (len < 0)
  {
    CTABfree (&ct);
    ErrorReturn(NULL, (ERROR_BADFILE, "CTABreadFromBinaryV1: file name "
                       "length was %d", len));
  }
  name = (char*) malloc (len+1);
  fread (name, sizeof(char), len, fp);
  strncpy (ct->fname, name, sizeof(ct->fname));
  free (name);

  /* For each entry, read in the info. We assume these have sequential
     structure indices. */
  for (structure = 0; structure < ct->nentries; structure++)
  {
    ct->entries[structure] = (CTE*) malloc(sizeof(CTE));
    if (NULL == ct->entries[structure])
    {
      CTABfree (&ct);
      ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABreadFromBinaryV1: could "
                         "not allocate entry for structure %d", structure));
    }

    /* Read the structure name. */
    len = freadInt(fp);
    if (len < 0)
    {
      CTABfree (&ct);
      ErrorReturn(NULL, (ERROR_BADFILE, "CTABreadFromBinaryV1: structure "
                         "%d name length was %d", structure, len));
    }
    name = (char*) malloc (len+1);
    fread(name, sizeof(char), len, fp);
    strncpy (ct->entries[structure]->name, name,
             sizeof(ct->entries[structure]->name));

    ct->entries[structure]->ri = freadInt(fp);
    ct->entries[structure]->gi = freadInt(fp);
    ct->entries[structure]->bi = freadInt(fp);
    t = freadInt(fp);
    ct->entries[structure]->ai = 255-t; /* alpha = 255-trans */

    /* Now calculate the float versions. */
    ct->entries[structure]->rf =
      (float)ct->entries[structure]->ri / 255.0;
    ct->entries[structure]->gf =
      (float)ct->entries[structure]->gi / 255.0;
    ct->entries[structure]->bf =
      (float)ct->entries[structure]->bi / 255.0;
    ct->entries[structure]->af =
      (float)ct->entries[structure]->ai / 255.0;
  }

  return(ct);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int
CTABwriteIntoBinaryV1(COLOR_TABLE *ct, FILE *fp)
{
  int  i;
  int  t;

  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABwriteIntoBinaryV1: ct was NULL"));
  if (NULL==fp)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABwriteIntoBinaryV1: fp was NULL"));

  /* This is the old version which didn't have a version number
     written at the beginning of the stream, so just start
     writing. First we write the number of entries we have. */
  fwriteInt (ct->nentries, fp);

  /* Now the length of the filename and the filename with an extra
     character for the null terminator. */
  fwriteInt (strlen(ct->fname)+1, fp);
  fwrite (ct->fname, sizeof(char), strlen (ct->fname)+1, fp);

  /* Now for each bine, if it's not null, write it to the stream. We
     don't write the structure number as it's implicit in the
     order. Note that this doesn't save structure indecies if there
     are skipped entries properly, but that's a feature of v1. */
  for (i = 0; i < ct->nentries; i++)
  {
    if (NULL != ct->entries[i])
    {
      fwriteInt (strlen(ct->entries[i]->name)+1, fp);
      fwrite (ct->entries[i]->name,
              sizeof(char), strlen (ct->entries[i]->name)+1, fp);
      fwriteInt (ct->entries[i]->ri, fp);
      fwriteInt (ct->entries[i]->gi, fp);
      fwriteInt (ct->entries[i]->bi, fp);
      t = 255 - ct->entries[i]->ai; /* alpha = 255-trans */
      fwriteInt (t, fp);
    }
  }

  return(NO_ERROR);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *
CTABreadFromBinaryV2(FILE *fp)
{
  COLOR_TABLE        *ct;
  int                nentries, num_entries_to_read, i;
  int                structure, len;
  char               *name;
  int                t;

  /* Read the number of entries from the stream. Note that this is
     really the max structure index; some of these entries could be
     blank. */
  nentries = freadInt (fp);
  if (nentries < 0)
    ErrorReturn
      (NULL,
       (ERROR_BADFILE, "CTABreadFromBinaryV2: nentries was %d", nentries));

  /* Allocate our table. */
  ct = (COLOR_TABLE *)calloc (1, sizeof(COLOR_TABLE));
  if (ct == NULL)
    ErrorReturn
      (NULL,
       (ERROR_NO_MEMORY, "CTABreadFromBinaryV2: could not allocate table"));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = nentries;
  ct->entries = (COLOR_TABLE_ENTRY**)
                calloc (ct->nentries, sizeof(COLOR_TABLE_ENTRY*));
  if (ct->entries == NULL)
    ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABreadFromBinaryV2: could not "
                       "allocate %d entries", ct->nentries));

  /* We'll write the same version to binary. */
  ct->version = 2;

  /* Read the file name. Read it into a temp buffer first to avoid
     overflow. */
  len = freadInt (fp);
  if (len < 0)
    ErrorReturn(NULL, (ERROR_BADFILE, "CTABreadFromBinaryV2: file name length "
                       "was %d", len));
  name = (char*) malloc (len+1);
  fread (name, sizeof(char), len, fp);
  strncpy (ct->fname, name, sizeof(ct->fname));
  free (name);

  /* Read the number of entries to read. */
  num_entries_to_read = freadInt (fp);

  /* For each entry, read in the info. */
  for (i = 0; i < num_entries_to_read; i++)
  {
    /* Read a structure number first. */
    structure = freadInt(fp);
    if (structure < 0)
    {
      CTABfree (&ct);
      ErrorReturn(NULL, (ERROR_BADFILE, "CTABreadFromBinaryV2: read entry "
                         "index %d", structure));
    }

    /* See if we already have an entry here. */
    if (NULL != ct->entries[structure])
    {
      CTABfree (&ct);
      ErrorReturn(NULL,(ERROR_BADFILE, "CTABreadFromBinaryV2: Duplicate "
                        "structure %d", structure));
    }

    /* Create the entry */
    ct->entries[structure] = (CTE*) malloc (sizeof(CTE));
    if (NULL == ct->entries[structure])
      ErrorReturn
        (NULL, 
         (ERROR_NO_MEMORY, 
          "CTABreadFromBinaryV2: could not allocate entry for structure %d", 
          structure));
    
    /* Read the structure name. */
    len = freadInt (fp);
    if (len < 0)
    {
      CTABfree (&ct);
      ErrorReturn(NULL, (ERROR_BADFILE, "CTABreadFromBinaryV2: structure "
                         "%d name length was %d", structure, len));
    }
    name = (char*) malloc (len+1);
    fread(name, sizeof(char), len, fp);
    strncpy (ct->entries[structure]->name, name,
             sizeof(ct->entries[structure]->name));

    /* Read in the color. */
    ct->entries[structure]->ri = freadInt(fp);
    ct->entries[structure]->gi = freadInt(fp);
    ct->entries[structure]->bi = freadInt(fp);
    t = freadInt(fp);
    ct->entries[structure]->ai = 255-t; /* alpha = 255-trans */

    /* Now calculate the float versions. */
    ct->entries[structure]->rf =
      (float)ct->entries[structure]->ri / 255.0;
    ct->entries[structure]->gf =
      (float)ct->entries[structure]->gi / 255.0;
    ct->entries[structure]->bf =
      (float)ct->entries[structure]->bi / 255.0;
    ct->entries[structure]->af =
      (float)ct->entries[structure]->ai / 255.0;
  }

  return(ct);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int
CTABwriteIntoBinaryV2(COLOR_TABLE *ct, FILE *fp)
{
  int structure;
  int i, t;
  int num_entries_to_write;

  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABwriteIntoBinaryV2: ct was NULL"));
  if (NULL==fp)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABwriteIntoBinaryV2: fp was NULL"));

  /* Write our negative version number. */
  fwriteInt (-2, fp);

  /* First we write the number of entries we have. Note that this is
     really the max structure index; some of these entries could be
     blank. */
  fwriteInt (ct->nentries, fp);

  /* Now the length of the filename and the filename with an extra
     character for the null terminator. */
  fwriteInt (strlen(ct->fname)+1, fp);
  fwrite (ct->fname, sizeof(char), strlen (ct->fname)+1, fp);

  /* We have to run through our table and count our non-null
     entries. */
  num_entries_to_write = 0;
  for (i = 0; i < ct->nentries; i++)
    if (NULL != ct->entries[i])
      num_entries_to_write++;
  fwriteInt (num_entries_to_write, fp);

  /* Now for each bin, if it's not null, write it to the stream. */
  for (structure = 0; structure < ct->nentries; structure++)
  {
    if (NULL != ct->entries[structure])
    {
      /* Write the structure number, then name, then color
         info. */
      fwriteInt (structure, fp);
      fwriteInt (strlen(ct->entries[structure]->name)+1, fp);
      fwrite (ct->entries[structure]->name,
              sizeof(char), strlen (ct->entries[structure]->name)+1, fp);
      fwriteInt (ct->entries[structure]->ri, fp);
      fwriteInt (ct->entries[structure]->gi, fp);
      fwriteInt (ct->entries[structure]->bi, fp);
      t = 255 - ct->entries[structure]->ai; /* alpha = 255-trans */
      fwriteInt (t, fp);
    }
  }

  return(NO_ERROR);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int
CTABcopyFileName(COLOR_TABLE *ct, char *name, size_t name_len)
{
  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABcopyName: ct was NULL"));
  if (NULL==name)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABcopyName: output parameter was NULL"));

  strncpy (name, ct->fname, name_len);

  return(NO_ERROR);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABgetNumberOfValidEntries(COLOR_TABLE *ct, int *num)
{
  int valid_entries;
  int structure;

  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABgetNumberOfValidEntries: ct was NULL"));
  if (NULL==num)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABgetNumberOfValidEntries: num was NULL"));

  /* Count the non-NULL entries. */
  valid_entries = 0;
  for (structure = 0; structure < ct->nentries; structure++)
    if (NULL != ct->entries[structure])
      valid_entries++;

  *num = valid_entries;
  return (NO_ERROR);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABgetNumberOfTotalEntries(COLOR_TABLE *ct, int *num)
{
  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABgetNumberOfTotalEntries: ct was NULL"));
  if (NULL==num)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABgetNumberOfTotalEntries: num was NULL"));

  *num = ct->nentries;
  return (NO_ERROR);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABisEntryValid(COLOR_TABLE *ct, int index, int *valid)
{
  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABisEntryValid: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABisEntryValid: index %d was OOB", index));
  if (NULL==valid)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABisEntryValid: valid was NULL"));

  /* Return whether or not this entry is not NULL. */
  *valid = (NULL != ct->entries[index]);

  return (NO_ERROR);

}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int
CTABrgbAtIndexi(COLOR_TABLE *ct, int index, int*r, int*g, int*b)
{
  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABrgbAtIndexi: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABrgbAtIndexi: index %d was OOB", index));
  if (NULL==r || NULL==g || NULL==b)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABrgbAtIndexi: output parameter was NULL"));

  if (NULL == ct->entries[index])
    return (ERROR_BADPARM);

  *r = ct->entries[index]->ri;
  *g = ct->entries[index]->gi;
  *b = ct->entries[index]->bi;

  return(NO_ERROR);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int
CTABrgbAtIndexf(COLOR_TABLE *ct, int index, float*r, float*g, float*b)
{
  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABrgbAtIndexf: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABrgbAtIndexf: index %d was OOB", index));
  if (NULL==r || NULL==g || NULL==b)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "CTABrgbAtIndexf: output parameter was NULL"));

  if (NULL == ct->entries[index])
    return (ERROR_BADPARM);

  *r = ct->entries[index]->rf;
  *g = ct->entries[index]->gf;
  *b = ct->entries[index]->bf;

  return(NO_ERROR);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int
CTABrgbaAtIndexi(COLOR_TABLE *ct, int index, int*r, int*g, int*b, int*a)
{
  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABrgbaAtIndexi: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABrgbaAtIndexi: index %d was OOB", index));
  if (NULL==r || NULL==g || NULL==b || NULL==a)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABrgbaAtIndexi: output parameter was NULL"));

  if (NULL == ct->entries[index])
    return (ERROR_BADPARM);

  *r = ct->entries[index]->ri;
  *g = ct->entries[index]->gi;
  *b = ct->entries[index]->bi;
  *a = ct->entries[index]->ai;

  return(NO_ERROR);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABrgbaAtIndexf(COLOR_TABLE *ct, int index,
                     float*r, float*g, float*b, float*a)
{
  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABrgbaAtIndexf: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABrgbaAtIndexf: index %d was OOB", index));
  if (NULL==r || NULL==g || NULL==b || NULL==a)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABrgbaAtIndexf: output parameter was NULL"));

  if (NULL == ct->entries[index])
    return (ERROR_BADPARM);

  *r = ct->entries[index]->rf;
  *g = ct->entries[index]->gf;
  *b = ct->entries[index]->bf;
  *a = ct->entries[index]->af;

  return(NO_ERROR);
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABcopyName(COLOR_TABLE *ct, int index, char *name, size_t name_len)
{
  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABcopyName: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABcopyName: index %d was OOB", index));
  if (NULL==name)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABcopyName: output parameter was NULL"));

  if (NULL == ct->entries[index])
    return (ERROR_BADPARM);

  strncpy (name, ct->entries[index]->name, name_len);

  return(NO_ERROR);
}


/*-------------------------------------------------------------------
  CTABannotationAtIndex() - given the index into the ctab, compute
  the annotation from the rgb values.
  ----------------------------------------------------------------*/
int CTABannotationAtIndex(COLOR_TABLE *ct, int index, int *annot)
{
  CTE *e;
  int  annotation;

  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABannotationAtIndex: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABannotationAtIndex: index %d was OOB", index));
  if (NULL==annot)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABannotationAtIndex: output parameter was NULL"));

  /* Shift the values over into a single integer. */
  e = ct->entries[index];
  if (NULL == ct->entries[index])
    return (ERROR_BADPARM);

  annotation = (e->bi << 16) + (e->gi << 8) + e->ri;
  *annot = annotation;

  return(NO_ERROR);
}


/*------------------------------------------------------------------
  CTABfindAnnotation() - given the annotation of a structure, return
  the index into the color table that matches the annotation (yes,
  this function is named incorrectly -- it should be something like
  CTABannot2Index).
  -----------------------------------------------------------------*/
int CTABfindAnnotation(COLOR_TABLE *ct, int annotation, int *index)
{
  int r, g, b;
  int result;

  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABfindAnnotation: ct was NULL"));
  if (NULL==index)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABfindAnnotation: output parameter was NULL"));

  /* Separate the annotation into colors and then find the color. */
  r = annotation & 0x0000ff;
  g = (annotation >> 8) & 0x0000ff;
  b = (annotation >> 16) & 0x0000ff;

  result = CTABfindRGBi (ct, r, g, b, index);

  return (result);
}


/*------------------------------------------------------------------
  CTABfindRGBi() - given the RGB of a structure, return the
  index into the color table that matches the RGB (yes, this function
  is named incorrectly -- it should be something like CTABrgb2Index).
  -----------------------------------------------------------------*/
int CTABfindRGBi(COLOR_TABLE *ct, int r, int g, int b, int *index)
{
  int structure;

  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABfindRGBi: ct was NULL"));
  if (r<0 || r>=256 || g<0 || g>=256 || b<0 || b>=256)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABfindRGBi: rgb was invalid"));
  if (NULL==index)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABfindRGBi: output parameter was NULL"));

  for (structure = 0; structure < ct->nentries; structure++)
  {
    if (NULL != ct->entries[structure])
    {
      if (ct->entries[structure]->ri == r &&
          ct->entries[structure]->gi == g &&
          ct->entries[structure]->bi == b)
      {
        *index = structure;
        return (NO_ERROR);
      }
    }
  }

  /* If we got here, we didn't find a match. */
  *index = -1;
  return (NO_ERROR);
}


/*------------------------------------------------------------------
  CTABfindName() - given the string name of a structure, return the
  index into the color table that matches the name (yes, this function
  is named incorrectly -- it should be something like CTABname2Index).
  -----------------------------------------------------------------*/
int CTABfindName(COLOR_TABLE *ct, char *name, int *index)
{
  int structure;

  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABfindName: ct was NULL"));
  if (NULL==name)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABfindName: name was NULL"));
  if (NULL==index)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABfindName: output parameter was NULL"));

  for (structure = 0; structure < ct->nentries; structure++)
  {
    if (NULL != ct->entries[structure])
    {
      if (stricmp(name, ct->entries[structure]->name) == 0)
      {
        *index = structure;
        return (NO_ERROR);
      }
    }
  }

  /* If we got here, we didn't find a match. */
  *index = -1;
  return (NO_ERROR);
}



/*--------------------------------------------------------------*/
int CTABprintASCII(COLOR_TABLE *ct, FILE *fp)
{
  int structure;

  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABprintASCII: ct was NULL"));
  if (NULL==fp)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABprintASCII: fp was NULL"));

  for (structure = 0; structure < ct->nentries; structure++)
  {
    if (NULL != ct->entries[structure])
    {
      fprintf (fp, "%3d  %-30s  %3d %3d %3d  %3d\n",
               structure, ct->entries[structure]->name,
               ct->entries[structure]->ri,
               ct->entries[structure]->gi,
               ct->entries[structure]->bi,
               255 - ct->entries[structure]->ai); /* alpha = 255-trans */
    }
  }

  return (ERROR_NONE);
}


/*--------------------------------------------------------------*/
int CTABwriteFileASCII(COLOR_TABLE *ct, char *fname)
{
  FILE *fp;
  int result;

  if (NULL==ct)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABwriteFileASCII: ct was NULL"));
  if (NULL==fname)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "CTABwriteFileASCII: fname was NULL"));


  fp = fopen(fname,"w");
  if (fp == NULL)
  {
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, 
                 "CTABwriteFileASCII(%s): could not open for writing\n", 
                 fname));
  }

  result = CTABprintASCII (ct, fp);
  fclose (fp);

  return (result);
}
