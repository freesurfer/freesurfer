/**
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cma.h"
#include "colortab.h"
#include "error.h"
#include "fio.h"
#include "fsenv.h"
#include "proto.h"
#include "utils.h"

/* Different binary i/o versions. */
static int CTABwriteIntoBinaryV1(COLOR_TABLE *ct, FILE *fp);
static int CTABwriteIntoBinaryV2(COLOR_TABLE *ct, FILE *fp);
static COLOR_TABLE *CTABreadFromBinaryV1(FILE *fp, int nentries);
static COLOR_TABLE *CTABreadFromBinaryV2(FILE *fp);

static int znzCTABwriteIntoBinaryV1(COLOR_TABLE *ct, znzFile fp);
static int znzCTABwriteIntoBinaryV2(COLOR_TABLE *ct, znzFile fp);
static COLOR_TABLE *znzCTABreadFromBinaryV1(znzFile fp, int nentries);
static COLOR_TABLE *znzCTABreadFromBinaryV2(znzFile fp);

int ctabDuplicates;

/*-------------------------------------------------------------------
COLOR_TABLE *CTABreadASCII(const char *fname)
Reads in 7th column as tissue type
  ----------------------------------------------------------------*/

COLOR_TABLE *CTABreadASCII(const char *fname) { return CTABreadASCII2(fname, 1); }

/* Discussions of Label File, LUT, and Annotation:
 *   https://surfer.nmr.mgh.harvard.edu/fswiki/LabelsClutsAnnotationFiles
 *
 * LUT is read into COLOR_TABLE:
 *   a. first colum in LUT becomes COLOR_TABLE index 
 *   b. number of COLOR_TABLE entries to create is max(1st LUT column) + 1
 *   c. COLOR_TABLE index is 0 .. n
 *   d. w/ or w/o '0  Unknown 0    0      0      0' in the first line, index 0 will be there
 *      the difference: w/  the line, CTABfindAnnotation() returns 0  for annotation=0;
 *                      w/o the line, CTABfindAnnotation() returns -1 for annotation=0
 *   e. skipped numbers in LUT will leave holes (empty entries) in COLOR_TABLE,
 *      their corresponding index will have unknown annotations
 */
COLOR_TABLE *CTABreadASCII2(const char *fname, int checkDuplicateNames)
{
  COLOR_TABLE *ct;
  char line[STRLEN], *cp;
  int max_structure;
  FILE *fp;
  int structure;
  char name[STRLEN];
  int r, g, b, t, tt = 0;
  int line_num, nscan;

  /* Try to open the file. */
  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(NULL, (ERROR_NOFILE, "CTABreadASCII(%s): could not open file", fname));

  /* Scan through the file and see what our max entry number is. */
  max_structure = -1;
  while ((cp = fgetl(line, STRLEN, fp)) != NULL) {
    /* See if this line is in the right format. If not, it's
    probably just a comment and we can ignore it. */
    if (sscanf(line, "%d %s %d %d %d %d", &structure, name, &r, &g, &b, &t) == 6) {
      if (structure > max_structure) max_structure = structure;
    }
  }

  /* If no structures, bail. */
  if (max_structure <= 0) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_NOFILE, "CTABreadASCII(%s): badly formed file", fname));
  }

  /* Allocate our table. */
  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (ct == NULL) ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABreadASCII(%s): could not allocate table", fname));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = max_structure + 1;
  ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY *));
  if (ct->entries == NULL)
    ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABreadASCII(%s): could not allocate %d entries", fname, ct->nentries));

  /* Copy in the file name. */
  strncpy(ct->fname, fname, sizeof(ct->fname)-1);

  /* We'll write this version if we write to binary. */
  ct->version = CTAB_VERSION_TO_WRITE;

  /* Rewind the file and go through it again. For each entry we find,
     allocate a CTE. This will leave the items in the array for which
     we don't have entries NULL. */
  line_num = 1;
  ctabDuplicates = 0;
  rewind(fp);
  while ((cp = fgets(line, STRLEN, fp)) != NULL) {
    nscan = sscanf(line, "%d %s %d %d %d %d %d", &structure, name, &r, &g, &b, &t, &tt);
    if (nscan != 7) tt = -1;  // set tissue type to -1 if not there
    if (nscan >= 6) {
      /* If this entry already exists, there's a duplicate entry
         in the file. Warn, but then continue on.*/
      if (ct->entries[structure] != NULL) {
        printf(
            "CTABreadASCII(%s): Line %d: Duplicate structure "
            "index %d, was %s %d %d %d %d\n",
            fname,
            line_num,
            structure,
            ct->entries[structure]->name,
            ct->entries[structure]->ri,
            ct->entries[structure]->gi,
            ct->entries[structure]->bi,
            ct->entries[structure]->ai);
        ctabDuplicates++;
      }
      else {
        /* Try to create a new entry.*/
        ct->entries[structure] = (CTE *)malloc(sizeof(CTE));
        if (NULL == ct->entries[structure]) {
          fclose(fp);
          CTABfree(&ct);
          ErrorReturn(
              NULL,
              (ERROR_NO_MEMORY, "CTABreadASCII(%s): could not allocate entry for structure %d", fname, structure));
        }

        /* Fill out the entry. */
        strncpy(ct->entries[structure]->name, name, sizeof(ct->entries[structure]->name));
        ct->entries[structure]->ri = r;
        ct->entries[structure]->gi = g;
        ct->entries[structure]->bi = b;
        ct->entries[structure]->ai = (255 - t); /* alpha = 255-trans */

        /* Now calculate the float versions. */
        ct->entries[structure]->rf = (float)ct->entries[structure]->ri / 255.0;
        ct->entries[structure]->gf = (float)ct->entries[structure]->gi / 255.0;
        ct->entries[structure]->bf = (float)ct->entries[structure]->bi / 255.0;
        ct->entries[structure]->af = (float)ct->entries[structure]->ai / 255.0;
        ct->entries[structure]->TissueType = tt;
        ct->entries[structure]->count = 0;
      }
    }
    line_num++;
  }

  fclose(fp);

  if (checkDuplicateNames) CTABfindDuplicateNames(ct);

  // CTABfindDuplicateAnnotations(ct);

  // This will read the tissue type ctab if it is there
  ct->ctabTissueType = CTABreadASCIIttHeader(fname);

  /* Return the new color table. */
  return (ct);
}

/*
  \fn COLOR_TABLE *CTABreadASCIIttHeader(const char *fname)
  \brief reads tissue type from header of a ctab (line that start with
  #ctTType) as written by CTABwriteFileASCIItt().
 */
COLOR_TABLE *CTABreadASCIIttHeader(const char *fname)
{
  FILE *fp;
  COLOR_TABLE *ct = NULL;
  int nct, ni, segid, segidmax;
  char *cp, *str, line[5000];
  int structure;
  char name[STRLEN];
  int r, g, b, t;

  fp = fopen(fname, "r");
  if (fp == NULL) {
    printf("ERROR: CTABreadASCIItt(): cannot open %s\n", fname);
    return (NULL);
  }

  segidmax = -1;
  nct = 0;
  cp = (char *)1;
  while (cp != NULL) {
    cp = fgets(line, STRLEN, fp);
    ni = CountItemsInString(line);
    if (ni == 0) continue;
    str = GetNthItemFromString(line, 0);
    if (strcmp(str, "#ctTType") != 0) {
      free(str);
      continue;
    }
    free(str);
    str = GetNthItemFromString(line, 1);
    sscanf(str, "%d", &segid);
    free(str);
    if (segidmax < segid) segidmax = segid;
    nct++;
  }
  // printf("nct = %d, segidmax %d\n",nct,segidmax);
  if (nct == 0) return (NULL);

  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  ct->nentries = segidmax + 1;
  ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY *));

  rewind(fp);

  nct = 0;
  cp = (char *)1;
  while (cp != NULL) {
    cp = fgets(line, STRLEN, fp);
    ni = CountItemsInString(line);
    if (ni == 0) continue;
    str = GetNthItemFromString(line, 0);
    if (strcmp(str, "#ctTType") != 0) {
      free(str);
      continue;
    }
    free(str);
    sscanf(line, "%*s %d %s  %d %d %d  %d", &structure, name, &r, &g, &b, &t);
    ct->entries[structure] = (COLOR_TABLE_ENTRY *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
    strncpy(ct->entries[structure]->name, name, STRLEN);
    ct->entries[structure]->ri = r;
    ct->entries[structure]->gi = g;
    ct->entries[structure]->bi = b;
    ct->entries[structure]->ai = (255 - t); /* alpha = 255-trans */
    ct->entries[structure]->TissueType = structure;
    nct++;
  }

  fclose(fp);
  return (ct);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABfree(COLOR_TABLE **pct)
{
  int i;
  COLOR_TABLE *ct;

  if (NULL == pct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfree: pct was NULL"));
  if (NULL == *pct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfree: *pct was NULL"));

  ct = *pct;

  /* Free all the entries. */
  for (i = 0; i < ct->nentries; i++)
    if (NULL != ct->entries[i]) free(ct->entries[i]);

  if (ct->ctabTissueType) CTABfree(&ct->ctabTissueType);

  free(ct->entries);
  free(ct);

  /* Set argument to null */
  *pct = NULL;

  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *CTABdeepCopy(COLOR_TABLE *ct)
{
  COLOR_TABLE *copy;
  int structure;

  if (NULL == ct) ErrorReturn(NULL, (ERROR_BADPARM, "CTABdeepCopy: ct was NULL"));

  /* Make a copy of the table. Allocate our table. */
  copy = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (copy == NULL) ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABdeepCopy: could not allocate table"));

  /* Allocate the array of NULL CTE ptrs. */
  copy->nentries = ct->nentries;
  copy->entries = (COLOR_TABLE_ENTRY **)calloc(copy->nentries, sizeof(COLOR_TABLE_ENTRY *));
  if (copy->entries == NULL)
    ErrorReturn(NULL,
                (ERROR_NO_MEMORY,
                 "CTABdeepCopy: could not allocate "
                 "%d entries",
                 copy->nentries));

  /* Copy in the file name and verion. */
  strncpy(copy->fname, ct->fname, sizeof(copy->fname));
  copy->version = ct->version;

  /* Go through the table. For each entry we find, allocate a CTE in
     the copy.*/
  for (structure = 0; structure < ct->nentries; structure++) {
    if (ct->entries[structure] != NULL) {
      /* Try to create a new entry.*/
      copy->entries[structure] = (CTE *)malloc(sizeof(CTE));
      if (NULL == copy->entries[structure]) {
        CTABfree(&copy);
        ErrorReturn(NULL,
                    (ERROR_NO_MEMORY,
                     "CTABdeepCopy: could not "
                     "allocate entry for structure %d",
                     structure));
      }

      /* Copy the entry. */
      strncpy(copy->entries[structure]->name, ct->entries[structure]->name, sizeof(copy->entries[structure]->name));
      copy->entries[structure]->ri = ct->entries[structure]->ri;
      copy->entries[structure]->gi = ct->entries[structure]->gi;
      copy->entries[structure]->bi = ct->entries[structure]->bi;
      copy->entries[structure]->ai = ct->entries[structure]->ai;
      copy->entries[structure]->rf = ct->entries[structure]->rf;
      copy->entries[structure]->gf = ct->entries[structure]->gf;
      copy->entries[structure]->bf = ct->entries[structure]->bf;
      copy->entries[structure]->af = ct->entries[structure]->af;
      copy->entries[structure]->TissueType = ct->entries[structure]->TissueType;
    }
  }
  if (ct->ctabTissueType) copy->ctabTissueType = CTABdeepCopy(ct->ctabTissueType);
  strcpy(copy->TissueTypeSchema,ct->TissueTypeSchema);
  /* Return the new copy. */
  return copy;
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *CTABreadFromBinary(FILE *fp)
{
  COLOR_TABLE *ct;
  int version;

  if (NULL == fp) ErrorReturn(NULL, (ERROR_BADPARM, "CTABreadFromBinary: fp was NULL"));

  /* Our goal here is to see what vesion we're reading/writing. Look
     at the first int in the stream. If it's > 0, it's the old format,
     and this int is the number of bins to have. If it's < 0, it's the
     new format, and is the negative version number. */
  version = freadInt(fp);

  ct = NULL;
  if (version > 0) {
    /* With v1, we pass in the "version" number we just got, as it's
    the number of entries. */
    ct = CTABreadFromBinaryV1(fp, version);
  }
  else {
    /* Take the negative to get the real version number. */
    version = -version;

    /* Read the right version. */
    if (version == 2) {
      ct = CTABreadFromBinaryV2(fp);
    }
    else {
      /* Unsupported version. */
      ErrorReturn(NULL, (ERROR_BADFILE, "CTABreadFromBinary: unknown version"));
    }
  }

  return (ct);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABwriteIntoBinary(COLOR_TABLE *ct, FILE *fp)
{
  int result;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinary: ct was NULL"));
  if (NULL == fp) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinary: fp was NULL"));

  /* Just switch on the version we want to write. */
  switch (ct->version) {
    case 1:
      result = CTABwriteIntoBinaryV1(ct, fp);
      break;
    case 2:
      result = CTABwriteIntoBinaryV2(ct, fp);
      break;
    default:
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinary: unknown version"));
  }

  return (result);
}
/*------------------------------------------------------------------------*/
/*!
  \fn COLOR_TABLE *CTABalloc(int nentries)
  \brief Allocates a color table and fills it with random
  entries. Note: if the RGBs are not unique, then all hell breaks
  loose for surface annotations, so there is a search at the end of
  this function which uses a brute force search to over 10 max
  iterations to find a unique set.
*/
COLOR_TABLE *CTABalloc(int nentries)
{
  COLOR_TABLE *ct;
  int structure;

  if (nentries < 0) ErrorReturn(NULL, (ERROR_BADPARM, "CTABalloc: nentries was %d", nentries));

  /* Allocate our table. */
  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (ct == NULL)
    ErrorReturn(NULL,
                (ERROR_NO_MEMORY,
                 "CTABalloc(%d): could not "
                 "allocate table",
                 nentries));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = nentries;
  ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY *));
  if (ct->entries == NULL) {
    CTABfree(&ct);
    ErrorReturn(NULL,
                (ERROR_NO_MEMORY,
                 "CTABalloc: could not "
                 "allocate %d entries",
                 ct->nentries));
  }

  /* Copy in default name. */
  strcpy(ct->fname, "none");

  /* We'll write this version if we write to binary. */
  ct->version = CTAB_VERSION_TO_WRITE;

  /* Init all entries to random colors. */
  for (structure = 0; structure < ct->nentries; structure++) {
    /* Try to create a new entry.*/
    ct->entries[structure] = (CTE *)malloc(sizeof(CTE));
    if (NULL == ct->entries[structure]) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_NO_MEMORY,
                   "CTABalloc: could not allocate "
                   "entry for structure %d",
                   structure));
    }

    /* Random colors. */
    ct->entries[structure]->ri = nint(randomNumber(0, 255));
    ct->entries[structure]->gi = nint(randomNumber(0, 255));
    ct->entries[structure]->bi = nint(randomNumber(0, 255));
    ct->entries[structure]->ai = 255;
      
    /* Now calculate the float versions. */
    ct->entries[structure]->rf = (float)ct->entries[structure]->ri / 255.0;
    ct->entries[structure]->gf = (float)ct->entries[structure]->gi / 255.0;
    ct->entries[structure]->bf = (float)ct->entries[structure]->bi / 255.0;
    ct->entries[structure]->af = (float)ct->entries[structure]->ai / 255.0;

    /* Make a fake name. */
    sprintf(ct->entries[structure]->name, "cluster%d", structure);
  }
  ct->idbase = 0;
  CTABunique(ct, 100);

  return (ct);
}
/*------------------------------------------------------------------------*/
/*!
  \fn int CTABunique(COLOR_TABLE *ct, int nmax)
  \brief Fill an already allocated color table with unique, random
  entries.  First entry set to unknown with RBG=0. This is just a
  brute force search over at most nmax tries to find a set that are
  unique.
*/
int CTABunique(COLOR_TABLE *ct, int nmax)
{
  int n;
  n = 0;
  while (n < nmax) {
    n++;
    CTABrandom(ct);
    if(CTABcountRepeats(ct,1)==0) break;
  }
  // Set first entry to unknown with RBG=0
  ct->entries[0]->ri = 0;
  ct->entries[0]->gi = 0;
  ct->entries[0]->bi = 0;
  ct->entries[0]->ai = 0;
  ct->entries[0]->rf = 0;
  ct->entries[0]->gf = 0;
  ct->entries[0]->bf = 0;
  ct->entries[0]->af = 0;
  sprintf(ct->entries[0]->name, "unknown");
  if (n == nmax) {
    printf("INFO: CTABunique() could not find a unique set in %d tries\n", nmax);
    return (-1);
  }
  return (0);
}
/*------------------------------------------------------------------------*/
/*!
  \fn int CTABcountRepeats(COLOR_TABLE *ct, int break_after_found)
  \brief Returns the number of pairs of entries that have the same RGB
  \param if break_after_found !=0 will stop after finding the first
         repeat (as in call from CTABunique to prevent N^2 search
*/
int CTABcountRepeats(COLOR_TABLE *ct, int break_after_found)
{
  int i,j,nrepeats;
  nrepeats=0;
  for(i = 0; i < ct->nentries; i++) {
    for(j = i+1; j < ct->nentries; j++) {
      if(i==j) continue;
      if(ct->entries[i]->rf == ct->entries[j]->rf &&
	 ct->entries[i]->gf == ct->entries[j]->gf &&
	 ct->entries[i]->bf == ct->entries[j]->bf){
	//printf("Entries %d and %d have the same RGB\n",i,j);
	nrepeats++;
	if (break_after_found)
	  break ;
      }
    }
  }
  // printf("Found %d nrepeatss\n",nrepeats);
  return (nrepeats);
}
/*------------------------------------------------------------------------*/
/*!
  \fn int CTABrandom(COLOR_TABLE *ct)
  \brief Fills an already allocated color table with random colors
*/
int CTABrandom(COLOR_TABLE *ct)
{
  int structure;
  /* Set all entries to random colors. */
  for (structure = 0; structure < ct->nentries; structure++) {
    /* Random colors. */
    ct->entries[structure]->ri = nint(randomNumber(0, 255));
    ct->entries[structure]->gi = nint(randomNumber(0, 255));
    ct->entries[structure]->bi = nint(randomNumber(0, 255));
    ct->entries[structure]->ai = 255;
    /* Now calculate the float versions. */
    ct->entries[structure]->rf = (float)ct->entries[structure]->ri / 255.0;
    ct->entries[structure]->gf = (float)ct->entries[structure]->gi / 255.0;
    ct->entries[structure]->bf = (float)ct->entries[structure]->bi / 255.0;
    ct->entries[structure]->af = (float)ct->entries[structure]->ai / 255.0;
  }
  return (0);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *CTABreadFromBinaryV1(FILE *fp, int nentries)
{
  COLOR_TABLE *ct;
  int structure, len;
  char *name;
  int t;

  if (nentries < 0) ErrorReturn(NULL, (ERROR_BADPARM, "CTABreadFromBinaryV1: nentries was %d", nentries));

  /* Allocate our table. */
  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (ct == NULL) ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABreadFromBinaryV1: could not allocate table"));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = nentries;
  ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY *));
  if (ct->entries == NULL) {
    CTABfree(&ct);
    ErrorReturn(NULL,
                (ERROR_NO_MEMORY,
                 "CTABreadFromBinaryV1: could not "
                 "allocate %d entries",
                 ct->nentries));
  }

  /* We'll write the same version to binary. */
  ct->version = 1;

  /* Read the file name. Read it into a temp buffer first to avoid
     overflow. */
  len = freadInt(fp);
  if (len < 0) {
    CTABfree(&ct);
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "CTABreadFromBinaryV1: file name "
                 "length was %d",
                 len));
  }
  name = (char *)malloc(len + 1);
  if (fread(name, sizeof(char), len, fp) != (unsigned)len) {
    ErrorPrintf(ERROR_BADFILE, "CTABreadFromBinaryV1: could not read parameter(s)");
  }
  strncpy(ct->fname, name, STRLEN-1);
  free(name);

  /* For each entry, read in the info. We assume these have sequential
     structure indices. */
  for (structure = 0; structure < ct->nentries; structure++) {
    ct->entries[structure] = (CTE *)malloc(sizeof(CTE));
    if (NULL == ct->entries[structure]) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_NO_MEMORY,
                   "CTABreadFromBinaryV1: could "
                   "not allocate entry for structure %d",
                   structure));
    }

    /* Read the structure name. */
    len = freadInt(fp);
    if (len < 0) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "CTABreadFromBinaryV1: structure "
                   "%d name length was %d",
                   structure,
                   len));
    }
    name = (char *)malloc(len + 1);
    if (fread(name, sizeof(char), len, fp) != (unsigned)len) {
      ErrorPrintf(ERROR_BADFILE, "CTABreadFromBinaryV1: could not read parameter(s)");
    }
    strncpy(ct->entries[structure]->name, name, STRLEN-1);
    ct->entries[structure]->name[len] = 0;

    ct->entries[structure]->ri = freadInt(fp);
    ct->entries[structure]->gi = freadInt(fp);
    ct->entries[structure]->bi = freadInt(fp);
    t = freadInt(fp);
    ct->entries[structure]->ai = 255 - t; /* alpha = 255-trans */

    /* Now calculate the float versions. */
    ct->entries[structure]->rf = (float)ct->entries[structure]->ri / 255.0;
    ct->entries[structure]->gf = (float)ct->entries[structure]->gi / 255.0;
    ct->entries[structure]->bf = (float)ct->entries[structure]->bi / 255.0;
    ct->entries[structure]->af = (float)ct->entries[structure]->ai / 255.0;
  }

  return (ct);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABwriteIntoBinaryV1(COLOR_TABLE *ct, FILE *fp)
{
  int i;
  int t;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinaryV1: ct was NULL"));
  if (NULL == fp) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinaryV1: fp was NULL"));

  /* This is the old version which didn't have a version number
     written at the beginning of the stream, so just start
     writing. First we write the number of entries we have. */
  fwriteInt(ct->nentries, fp);

  /* Now the length of the filename and the filename with an extra
     character for the null terminator. */
  fwriteInt(strlen(ct->fname) + 1, fp);
  fwrite(ct->fname, sizeof(char), strlen(ct->fname) + 1, fp);

  /* Now for each bine, if it's not null, write it to the stream. We
     don't write the structure number as it's implicit in the
     order. Note that this doesn't save structure indecies if there
     are skipped entries properly, but that's a feature of v1. */
  for (i = 0; i < ct->nentries; i++) {
    if (NULL != ct->entries[i]) {
      fwriteInt(strlen(ct->entries[i]->name) + 1, fp);
      fwrite(ct->entries[i]->name, sizeof(char), strlen(ct->entries[i]->name) + 1, fp);
      fwriteInt(ct->entries[i]->ri, fp);
      fwriteInt(ct->entries[i]->gi, fp);
      fwriteInt(ct->entries[i]->bi, fp);
      t = 255 - ct->entries[i]->ai; /* alpha = 255-trans */
      fwriteInt(t, fp);
    }
  }

  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *CTABreadFromBinaryV2(FILE *fp)
{
  COLOR_TABLE *ct;
  int nentries, num_entries_to_read, i;
  int structure, len;
  char *name;
  int t;

  /* Read the number of entries from the stream. Note that this is
     really the max structure index; some of these entries could be
     blank. */
  nentries = freadInt(fp);
  if (nentries < 0) ErrorReturn(NULL, (ERROR_BADFILE, "CTABreadFromBinaryV2: nentries was %d", nentries));

  /* Allocate our table. */
  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (ct == NULL) ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABreadFromBinaryV2: could not allocate table"));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = nentries;
  ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY *));
  if (ct->entries == NULL)
    ErrorReturn(NULL,
                (ERROR_NO_MEMORY,
                 "CTABreadFromBinaryV2: could not "
                 "allocate %d entries",
                 ct->nentries));

  /* We'll write the same version to binary. */
  ct->version = 2;

  /* Read the file name. Read it into a temp buffer first to avoid
     overflow. */
  len = freadInt(fp);
  if (len < 0)
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "CTABreadFromBinaryV2: file name length "
                 "was %d",
                 len));
  name = (char *)malloc(len + 1);
  if (fread(name, sizeof(char), len, fp) != (unsigned)len) {
    ErrorPrintf(ERROR_BADFILE, "CTABreadFromBinaryV1: could not read parameter(s)");
  }
  strncpy(ct->fname, name, STRLEN-1);
  free(name);

  /* Read the number of entries to read. */
  num_entries_to_read = freadInt(fp);

  if (Gdiag & DIAG_SHOW)
    printf("[DEBUG] CTABreadFromBinaryV2(): ct->nentries=%d, num_entries_to_read=%d\n", ct->nentries, num_entries_to_read);
 
  /* For each entry, read in the info. */
  for (i = 0; i < num_entries_to_read; i++) {
    /* Read a structure number first. */
    structure = freadInt(fp);
    if (structure < 0) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "CTABreadFromBinaryV2: read entry "
                   "index %d",
                   structure));
    }

    /* See if we already have an entry here. */
    if (NULL != ct->entries[structure]) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "CTABreadFromBinaryV2: Duplicate "
                   "structure %d",
                   structure));
    }

    /* Create the entry */
    ct->entries[structure] = (CTE *)malloc(sizeof(CTE));
    if (NULL == ct->entries[structure])
      ErrorReturn(NULL,
                  (ERROR_NO_MEMORY, "CTABreadFromBinaryV2: could not allocate entry for structure %d", structure));

    /* Read the structure name. */
    len = freadInt(fp);
    if (len < 0) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "CTABreadFromBinaryV2: structure "
                   "%d name length was %d",
                   structure,
                   len));
    }
    name = (char *)malloc(len + 1);
    if (fread(name, sizeof(char), len, fp) != (unsigned)len) {
      ErrorPrintf(ERROR_BADFILE, "CTABreadFromBinaryV1: could not read parameter(s)");
    }
    strncpy(ct->entries[structure]->name, name, STRLEN-1);

    /* Read in the color. */
    ct->entries[structure]->ri = freadInt(fp);
    ct->entries[structure]->gi = freadInt(fp);
    ct->entries[structure]->bi = freadInt(fp);
    t = freadInt(fp);
    ct->entries[structure]->ai = 255 - t; /* alpha = 255-trans */

    /* Now calculate the float versions. */
    ct->entries[structure]->rf = (float)ct->entries[structure]->ri / 255.0;
    ct->entries[structure]->gf = (float)ct->entries[structure]->gi / 255.0;
    ct->entries[structure]->bf = (float)ct->entries[structure]->bi / 255.0;
    ct->entries[structure]->af = (float)ct->entries[structure]->ai / 255.0;
  }

  return (ct);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABwriteIntoBinaryV2(COLOR_TABLE *ct, FILE *fp)
{
  int structure;
  int i, t;
  int num_entries_to_write;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinaryV2: ct was NULL"));
  if (NULL == fp) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinaryV2: fp was NULL"));

  /* CTAB binary format:
   *   negative version number (-2)
   *   number of CTAB entries include empty (null) ones (ct->nentries)
   *   length of file name string follows (strlen(ct->fname) + 1)
   *   file name (ct->fname)
   *   number of non-null entries
   *   (COLOR_TABLE_ENTRY) x (number of non-null entries)
   */

  /* Write our negative version number. */
  fwriteInt(-2, fp);

  /* First we write the number of entries we have. Note that this is
     really the max structure index; some of these entries could be
     blank. */
  fwriteInt(ct->nentries, fp);

  /* Now the length of the filename and the filename with an extra
     character for the null terminator. */
  fwriteInt(strlen(ct->fname) + 1, fp);
  fwrite(ct->fname, sizeof(char), strlen(ct->fname) + 1, fp);

  /* We have to run through our table and count our non-null
     entries. */
  num_entries_to_write = 0;
  for (i = 0; i < ct->nentries; i++)
    if (NULL != ct->entries[i]) num_entries_to_write++;
  fwriteInt(num_entries_to_write, fp);

  if (Gdiag & DIAG_SHOW)
    printf("[DEBUG] CTABwriteIntoBinaryV2(): ct->entries=%d, to_write=%d\n", ct->nentries, num_entries_to_write);

  /* Now for each bin, if it's not null, write it to the stream. */
  for (structure = 0; structure < ct->nentries; structure++) {
    if (NULL != ct->entries[structure]) {
      /* Write the structure number, then name, then color
         info. */
      fwriteInt(structure, fp);
      fwriteInt(strlen(ct->entries[structure]->name) + 1, fp);
      fwrite(ct->entries[structure]->name, sizeof(char), strlen(ct->entries[structure]->name) + 1, fp);
      fwriteInt(ct->entries[structure]->ri, fp);
      fwriteInt(ct->entries[structure]->gi, fp);
      fwriteInt(ct->entries[structure]->bi, fp);
      t = 255 - ct->entries[structure]->ai; /* alpha = 255-trans */
      fwriteInt(t, fp);
    }
  }

  return (NO_ERROR);
}

/*------------------------------------------------------------------
znzlib support
-------------------------------------------------------------------*/

COLOR_TABLE *znzCTABreadFromBinary(znzFile fp)
{
  COLOR_TABLE *ct;
  int version;

  if (znz_isnull(fp)) ErrorReturn(NULL, (ERROR_BADPARM, "CTABreadFromBinary: fp was NULL"));

  /* Our goal here is to see what vesion we're reading/writing. Look
  at the first int in the stream. If it's > 0, it's the old format,
  and this int is the number of bins to have. If it's < 0, it's the
  new format, and is the negative version number. */
  version = znzreadInt(fp);

  ct = NULL;
  if (version > 0) {
    /* With v1, we pass in the "version" number we just got, as it's
    the number of entries. */
    ct = znzCTABreadFromBinaryV1(fp, version);
  }
  else {
    /* Take the negative to get the real version number. */
    version = -version;

    /* Read the right version. */
    if (version == 2) {
      ct = znzCTABreadFromBinaryV2(fp);
    }
    else {
      /* Unsupported version. */
      ErrorReturn(NULL, (ERROR_BADFILE, "CTABreadFromBinary: unknown version"));
    }
  }

  return (ct);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int znzCTABwriteIntoBinary(COLOR_TABLE *ct, znzFile fp)
{
  int result;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinary: ct was NULL"));
  if (znz_isnull(fp)) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinary: fp was NULL"));

  /* Just switch on the version we want to write. */
  switch (ct->version) {
    case 1:
      result = znzCTABwriteIntoBinaryV1(ct, fp);
      break;
    case 2:
      result = znzCTABwriteIntoBinaryV2(ct, fp);
      break;
    default:
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinary: unknown version"));
  }

  return (result);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *znzCTABreadFromBinaryV1(znzFile fp, int nentries)
{
  COLOR_TABLE *ct;
  int structure, len;
  char *name;
  int t;

  if (nentries < 0) ErrorReturn(NULL, (ERROR_BADPARM, "CTABreadFromBinaryV1: nentries was %d", nentries));

  /* Allocate our table. */
  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (ct == NULL) ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABreadFromBinaryV1: could not allocate table"));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = nentries;
  ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY *));
  if (ct->entries == NULL) {
    CTABfree(&ct);
    ErrorReturn(NULL,
                (ERROR_NO_MEMORY,
                 "CTABreadFromBinaryV1: could not "
                 "allocate %d entries",
                 ct->nentries));
  }

  /* We'll write the same version to binary. */
  ct->version = 1;

  /* Read the file name. Read it into a temp buffer first to avoid
  overflow. */
  len = znzreadInt(fp);
  if (len < 0) {
    CTABfree(&ct);
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "CTABreadFromBinaryV1: file name "
                 "length was %d",
                 len));
  }
  name = (char *)malloc(len + 1);
  znzread(name, sizeof(char), len, fp);
  strncpy(ct->fname, name, sizeof(ct->fname)-1);
  free(name);

  /* For each entry, read in the info. We assume these have sequential
  structure indices. */
  for (structure = 0; structure < ct->nentries; structure++) {
    ct->entries[structure] = (CTE *)malloc(sizeof(CTE));
    if (NULL == ct->entries[structure]) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_NO_MEMORY,
                   "CTABreadFromBinaryV1: could "
                   "not allocate entry for structure %d",
                   structure));
    }

    /* Read the structure name. */
    len = znzreadInt(fp);
    if (len < 0) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "CTABreadFromBinaryV1: structure "
                   "%d name length was %d",
                   structure,
                   len));
    }
    name = (char *)malloc(len + 1);
    znzread(name, sizeof(char), len, fp);
    strncpy(ct->entries[structure]->name, name, sizeof(ct->entries[structure]->name)-1);
    ct->entries[structure]->name[len] = 0;

    ct->entries[structure]->ri = znzreadInt(fp);
    ct->entries[structure]->gi = znzreadInt(fp);
    ct->entries[structure]->bi = znzreadInt(fp);
    t = znzreadInt(fp);
    ct->entries[structure]->ai = 255 - t; /* alpha = 255-trans */

    /* Now calculate the float versions. */
    ct->entries[structure]->rf = (float)ct->entries[structure]->ri / 255.0;
    ct->entries[structure]->gf = (float)ct->entries[structure]->gi / 255.0;
    ct->entries[structure]->bf = (float)ct->entries[structure]->bi / 255.0;
    ct->entries[structure]->af = (float)ct->entries[structure]->ai / 255.0;
  }

  return (ct);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int znzCTABwriteIntoBinaryV1(COLOR_TABLE *ct, znzFile fp)
{
  int i;
  int t;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinaryV1: ct was NULL"));
  if (znz_isnull(fp)) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinaryV1: fp was NULL"));

  /* This is the old version which didn't have a version number
  written at the beginning of the stream, so just start
  writing. First we write the number of entries we have. */
  znzwriteInt(ct->nentries, fp);

  /* Now the length of the filename and the filename with an extra
  character for the null terminator. */
  znzwriteInt(strlen(ct->fname) + 1, fp);
  znzwrite(ct->fname, sizeof(char), strlen(ct->fname) + 1, fp);

  /* Now for each bine, if it's not null, write it to the stream. We
  don't write the structure number as it's implicit in the
  order. Note that this doesn't save structure indecies if there
  are skipped entries properly, but that's a feature of v1. */
  for (i = 0; i < ct->nentries; i++) {
    if (NULL != ct->entries[i]) {
      znzwriteInt(strlen(ct->entries[i]->name) + 1, fp);
      znzwrite(ct->entries[i]->name, sizeof(char), strlen(ct->entries[i]->name) + 1, fp);
      znzwriteInt(ct->entries[i]->ri, fp);
      znzwriteInt(ct->entries[i]->gi, fp);
      znzwriteInt(ct->entries[i]->bi, fp);
      t = 255 - ct->entries[i]->ai; /* alpha = 255-trans */
      znzwriteInt(t, fp);
    }
  }

  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
COLOR_TABLE *znzCTABreadFromBinaryV2(znzFile fp)
{
  COLOR_TABLE *ct;
  int nentries, num_entries_to_read, i;
  int structure, len;
  char *name;
  int t;

  /* Read the number of entries from the stream. Note that this is
  really the max structure index; some of these entries could be
  blank. */
  nentries = znzreadInt(fp);
  if (nentries < 0) ErrorReturn(NULL, (ERROR_BADFILE, "CTABreadFromBinaryV2: nentries was %d", nentries));

  /* Allocate our table. */
  ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (ct == NULL) ErrorReturn(NULL, (ERROR_NO_MEMORY, "CTABreadFromBinaryV2: could not allocate table"));

  /* Allocate the array of NULL CTE ptrs. */
  ct->nentries = nentries;
  ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY *));
  if (ct->entries == NULL)
    ErrorReturn(NULL,
                (ERROR_NO_MEMORY,
                 "CTABreadFromBinaryV2: could not "
                 "allocate %d entries",
                 ct->nentries));

  /* We'll write the same version to binary. */
  ct->version = 2;

  /* Read the file name. Read it into a temp buffer first to avoid
  overflow. */
  len = znzreadInt(fp);
  if (len < 0)
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "CTABreadFromBinaryV2: file name length "
                 "was %d",
                 len));
  name = (char *)malloc(len + 1);
  /* 
   * if the file comes from surfa.io.fsio.write_binary_lookup_table(), len = 0.
   * if the file comes from freesurfer/utils/colortab.cpp::znzCTABwriteIntoBinaryV2(), len > 0.
   */
  znzread(name, sizeof(char), len, fp);
  strncpy(ct->fname, name, STRLEN-1);
  free(name);

  /* Read the number of entries to read. */
  num_entries_to_read = znzreadInt(fp);

  /* For each entry, read in the info. */
  for (i = 0; i < num_entries_to_read; i++) {
    /* Read a structure number first. */
    structure = znzreadInt(fp);
    if (structure < 0) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "CTABreadFromBinaryV2: read entry "
                   "index %d",
                   structure));
    }

    /* See if we already have an entry here. */
    if (NULL != ct->entries[structure]) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "CTABreadFromBinaryV2: Duplicate "
                   "structure %d",
                   structure));
    }

    /* Create the entry */
    ct->entries[structure] = (CTE *)malloc(sizeof(CTE));
    if (NULL == ct->entries[structure])
      ErrorReturn(NULL,
                  (ERROR_NO_MEMORY, "CTABreadFromBinaryV2: could not allocate entry for structure %d", structure));

    /* Read the structure name. */
    len = znzreadInt(fp);
    if (len < 0) {
      CTABfree(&ct);
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "CTABreadFromBinaryV2: structure "
                   "%d name length was %d",
                   structure,
                   len));
    }
    name = (char *)malloc(len + 1);
    znzread(name, sizeof(char), len, fp);
    strncpy(ct->entries[structure]->name, name, STRLEN-1);

    /* Read in the color. */
    ct->entries[structure]->ri = znzreadInt(fp);
    ct->entries[structure]->gi = znzreadInt(fp);
    ct->entries[structure]->bi = znzreadInt(fp);
    t = znzreadInt(fp);
    ct->entries[structure]->ai = 255 - t; /* alpha = 255-trans */

    /* Now calculate the float versions. */
    ct->entries[structure]->rf = (float)ct->entries[structure]->ri / 255.0;
    ct->entries[structure]->gf = (float)ct->entries[structure]->gi / 255.0;
    ct->entries[structure]->bf = (float)ct->entries[structure]->bi / 255.0;
    ct->entries[structure]->af = (float)ct->entries[structure]->ai / 255.0;
  }

  return (ct);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int znzCTABwriteIntoBinaryV2(COLOR_TABLE *ct, znzFile fp)
{
  int structure;
  int i, t;
  int num_entries_to_write;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinaryV2: ct was NULL"));
  if (znz_isnull(fp)) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteIntoBinaryV2: fp was NULL"));

  /* Write our negative version number. */
  znzwriteInt(-2, fp);

  /* First we write the number of entries we have. Note that this is
  really the max structure index; some of these entries could be
  blank. */
  znzwriteInt(ct->nentries, fp);

  /* Now the length of the filename and the filename with an extra
  character for the null terminator. */
  znzwriteInt(strlen(ct->fname) + 1, fp);
  znzwrite(ct->fname, sizeof(char), strlen(ct->fname) + 1, fp);

  /* We have to run through our table and count our non-null
  entries. */
  num_entries_to_write = 0;
  for (i = 0; i < ct->nentries; i++)
    if (NULL != ct->entries[i]) num_entries_to_write++;
  znzwriteInt(num_entries_to_write, fp);

  /* Now for each bin, if it's not null, write it to the stream. */
  for (structure = 0; structure < ct->nentries; structure++) {
    if (NULL != ct->entries[structure]) {
      /* Write the structure number, then name, then color
      info. */
      znzwriteInt(structure, fp);
      znzwriteInt(strlen(ct->entries[structure]->name) + 1, fp);
      znzwrite(ct->entries[structure]->name, sizeof(char), strlen(ct->entries[structure]->name) + 1, fp);
      znzwriteInt(ct->entries[structure]->ri, fp);
      znzwriteInt(ct->entries[structure]->gi, fp);
      znzwriteInt(ct->entries[structure]->bi, fp);
      t = 255 - ct->entries[structure]->ai; /* alpha = 255-trans */
      znzwriteInt(t, fp);
    }
  }

  return (NO_ERROR);
}


// Reads the default color table from $FREESURFER_HOME/FreeSurferColorLUT.txt
COLOR_TABLE* CTABreadDefault() {
  FSENV *fsenv = FSENVgetenv();
  std::string filename = std::string(fsenv->FREESURFER_HOME) + "/FreeSurferColorLUT.txt";
  FSENVfree(&fsenv);
  return CTABreadASCII(filename.c_str());
}


/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABcopyFileName(COLOR_TABLE *ct, char *name, size_t name_len)
{
  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABcopyName: ct was NULL"));
  if (NULL == name) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABcopyName: output parameter was NULL"));

  strncpy(name, ct->fname, name_len);

  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABgetNumberOfValidEntries(COLOR_TABLE *ct, int *num)
{
  int valid_entries;
  int structure;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABgetNumberOfValidEntries: ct was NULL"));
  if (NULL == num) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABgetNumberOfValidEntries: num was NULL"));

  /* Count the non-NULL entries. */
  valid_entries = 0;
  for (structure = 0; structure < ct->nentries; structure++)
    if (NULL != ct->entries[structure]) valid_entries++;

  *num = valid_entries;
  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABgetNumberOfTotalEntries(COLOR_TABLE *ct, int *num)
{
  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABgetNumberOfTotalEntries: ct was NULL"));
  if (NULL == num) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABgetNumberOfTotalEntries: num was NULL"));

  *num = ct->nentries;
  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABisEntryValid(COLOR_TABLE *ct, int index, int *valid)
{
  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABisEntryValid: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABisEntryValid: index %d was OOB", index));
  if (NULL == valid) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABisEntryValid: valid was NULL"));

  /* Return whether or not this entry is not NULL. */
  *valid = (NULL != ct->entries[index]);

  return (NO_ERROR);
}

/*-------------------------------------------------------------------
This is a safer version of CTABisEntryValid that actually returns the
correct answer. I'm not fixing CTABisEntryValid since it's used in a lot
of places and I don't want to break any logic.
  ----------------------------------------------------------------*/
bool CTABhasEntry(COLOR_TABLE *ct, int index)
{
  if (NULL == ct) return false;
  if (index < 0 || index >= ct->nentries) return false;
  return (NULL != ct->entries[index]);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABrgbAtIndexi(COLOR_TABLE *ct, int index, int *r, int *g, int *b)
{
  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbAtIndexi: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbAtIndexi: index %d was OOB", index));
  if (NULL == r || NULL == g || NULL == b)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbAtIndexi: output parameter was NULL"));

  if (NULL == ct->entries[index]) return (ERROR_BADPARM);

  *r = ct->entries[index]->ri;
  *g = ct->entries[index]->gi;
  *b = ct->entries[index]->bi;

  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABrgbAtIndexf(COLOR_TABLE *ct, int index, float *r, float *g, float *b)
{
  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbAtIndexf: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbAtIndexf: index %d was OOB", index));
  if (NULL == r || NULL == g || NULL == b)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbAtIndexf: output parameter was NULL"));

  if (NULL == ct->entries[index]) return (ERROR_BADPARM);

  *r = ct->entries[index]->rf;
  *g = ct->entries[index]->gf;
  *b = ct->entries[index]->bf;

  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABrgbaAtIndexi(COLOR_TABLE *ct, int index, int *r, int *g, int *b, int *a)
{
  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbaAtIndexi: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbaAtIndexi: index %d was OOB", index));
  if (NULL == r || NULL == g || NULL == b || NULL == a)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbaAtIndexi: output parameter was NULL"));

  if (NULL == ct->entries[index]) return (ERROR_BADPARM);

  *r = ct->entries[index]->ri;
  *g = ct->entries[index]->gi;
  *b = ct->entries[index]->bi;
  *a = ct->entries[index]->ai;

  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABrgbaAtIndexf(COLOR_TABLE *ct, int index, float *r, float *g, float *b, float *a)
{
  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbaAtIndexf: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbaAtIndexf: index %d was OOB", index));
  if (NULL == r || NULL == g || NULL == b || NULL == a)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABrgbaAtIndexf: output parameter was NULL"));

  if (NULL == ct->entries[index]) return (ERROR_BADPARM);

  *r = ct->entries[index]->rf;
  *g = ct->entries[index]->gf;
  *b = ct->entries[index]->bf;
  *a = ct->entries[index]->af;

  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABcopyName(COLOR_TABLE *ct, int index, char *name, size_t name_len)
{
  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABcopyName: ct was NULL"));
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABcopyName: index %d was OOB", index));
  if (NULL == name) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABcopyName: output parameter was NULL"));

  if (NULL == ct->entries[index]) return (ERROR_BADPARM);

  strncpy(name, ct->entries[index]->name, name_len);

  return (NO_ERROR);
}

/*-------------------------------------------------------------------
  CTABrgb2Annotation(int r, int g, int b)
  Converts an rgb triplet into an annotation value.
  ----------------------------------------------------------------*/
int CTABrgb2Annotation(int r, int g, int b)
{
  int annotation;
  annotation = (b << 16) + (g << 8) + r;
  return (annotation);
}

/*-------------------------------------------------------------------
  CTABentryNameToIndex(char *EntryName, COLOR_TABLE *ct)
  Return the color table index given the name of the entry.
  ----------------------------------------------------------------*/
int CTABentryNameToIndex(const char *EntryName, COLOR_TABLE *ct)
{
  CTE *cte;
  int i;

  for (i = 0; i < ct->nentries; i++) {
    cte = ct->entries[i];
    // cte might be NULL, so this check is essential.
    if (cte != NULL) {
      if (!stricmp(cte->name, EntryName)) return (i);
    }
  }
  return (-1);  // error
}
/*-------------------------------------------------------------------
  CTABentryNameToIndex(char *EntryName, COLOR_TABLE *ct)
  Return the color table annotation given the name of the entry.
  ----------------------------------------------------------------*/
int CTABentryNameToAnnotation(const char *EntryName, COLOR_TABLE *ct)
{
  CTE *cte;
  int index, annotation;
  index = CTABentryNameToIndex(EntryName, ct);
  if (index == -1) return (-1);  // error
  cte = ct->entries[index];
  annotation = CTABrgb2Annotation(cte->ri, cte->gi, cte->bi);
  return (annotation);
}
/*-------------------------------------------------------------------
  CTABannotationAtIndex() - given the index into the ctab, compute
  the annotation from the rgb values.
  ----------------------------------------------------------------*/
int CTABannotationAtIndex(COLOR_TABLE *ct, int index, int *annot)
{
  CTE *e;
  int annotation;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABannotationAtIndex: ct was NULL"));
  // invalid index -1 is treated as error. This is different from index_to_annotation().
  if (index < 0 || index >= ct->nentries)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABannotationAtIndex: index %d was OOB", index));
  if (NULL == annot) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABannotationAtIndex: output parameter was NULL"));

  /* Shift the values over into a single integer. */
  e = ct->entries[index];
  if (NULL == ct->entries[index]) return (ERROR_BADPARM);

  // This should give the same, but have not tested it
  // annotation = CTABrgb2Annotation(e->ri, e->gi, e->bi);
  annotation = (e->bi << 16) + (e->gi << 8) + e->ri;
  *annot = annotation;

  return (NO_ERROR);
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

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindAnnotation: ct was NULL"));
  if (NULL == index) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindAnnotation: output parameter was NULL"));

  /* Separate the annotation into colors and then find the color. */
  r = annotation & 0x0000ff;
  g = (annotation >> 8) & 0x0000ff;
  b = (annotation >> 16) & 0x0000ff;

  result = CTABfindRGBi(ct, r, g, b, index);

  return (result);
}

/*------------------------------------------------------------------
  CTABgetAnnotationName() - given the annotation of a structure,
  return the name of this annotation label.  returns "NOT_FOUND" if
  annotation not found in colortable.
  -----------------------------------------------------------------*/
const char *CTABgetAnnotationName(COLOR_TABLE *ct, int annotation)
{
  int r, g, b;
  int index = -1;

  if (NULL == ct) ErrorExit(ERROR_BADPARM, "CTABfindAnnotationName: ct was NULL");

  /* Separate the annotation into colors and then find the color. */
  r = annotation & 0x0000ff;
  g = (annotation >> 8) & 0x0000ff;
  b = (annotation >> 16) & 0x0000ff;

  CTABfindRGBi(ct, r, g, b, &index);

  if ((index < 0) || (index >= ct->nentries)) return "NOT_FOUND";

  return ct->entries[index]->name;
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABfindDuplicateAnnotations(COLOR_TABLE *ct)
{
  int idx1;
  int idx2;
  int dupCount = 0;
  for (idx1 = 0; idx1 < ct->nentries; idx1++) {
    if (NULL != ct->entries[idx1]) {
      int annot1 = 0;
      CTABannotationAtIndex(ct, idx1, &annot1);
      for (idx2 = 0; idx2 < ct->nentries; idx2++) {
        if ((NULL != ct->entries[idx2]) && (idx1 != idx2)) {
          int annot2 = 0;
          CTABannotationAtIndex(ct, idx2, &annot2);
          if (annot1 == annot2) {
            dupCount++;
            printf(
                "ERROR: labels '%s'\n"
                "          and '%s'\n"
                "       have duplicate annotations!\n",
                ct->entries[idx1]->name,
                ct->entries[idx2]->name);
          }
        }
      }
    }
  }
  if (dupCount > 0) {
    printf(
        "CTABfindDuplicateAnnotations: %d duplicate annotation "
        "values found in:\n  %s",
        dupCount / 2,
        ct->fname);
  }
  return dupCount / 2;
}

/*-------------------------------------------------------------------
  ----------------------------------------------------------------*/
int CTABfindDuplicateNames(COLOR_TABLE *ct)
{
  int idx1;
  int idx2;
  int dupCount = 0;
  for (idx1 = 0; idx1 < ct->nentries; idx1++) {
    if (NULL != ct->entries[idx1]) {
      for (idx2 = 0; idx2 < ct->nentries; idx2++) {
        if ((NULL != ct->entries[idx2]) && (idx1 != idx2)) {
          char *name1 = ct->entries[idx1]->name;
          char *name2 = ct->entries[idx2]->name;
          if (strcmp(name1, name2) == 0) {
            dupCount++;
            printf("ERROR: label name '%s' is duplicated!\n", name1);
          }
        }
      }
    }
  }
  if (dupCount > 0) {
    printf("CTABfindDuplicateNames: %d duplicate names found in:\n  %s\n", dupCount / 2, ct->fname);
  }
  return dupCount / 2;
}

int CTABfindIndexFromAnnotation(COLOR_TABLE *ct, int annot, int *index)
{
  int r, g, b;

  AnnotToRGB(annot, r, g, b);  // macro
  return (CTABfindRGBi(ct, r, g, b, index));
}
/*------------------------------------------------------------------
  CTABfindRGBi() - given the RGB of a structure, return the
  index into the color table that matches the RGB (yes, this function
  is named incorrectly -- it should be something like CTABrgb2Index).
  -----------------------------------------------------------------*/
int CTABfindRGBi(COLOR_TABLE *ct, int r, int g, int b, int *index)
{
  int structure;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindRGBi: ct was NULL"));
  if (r < 0 || r >= 256 || g < 0 || g >= 256 || b < 0 || b >= 256)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindRGBi: rgb was invalid"));
  if (NULL == index) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindRGBi: output parameter was NULL"));

  for (structure = 0; structure < ct->nentries; structure++) {
    if (NULL != ct->entries[structure]) {
      if (ct->entries[structure]->ri == r && ct->entries[structure]->gi == g && ct->entries[structure]->bi == b) {
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
int CTABfindName(COLOR_TABLE *ct, const char *name, int *index)
{
  int structure;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindName: ct was NULL"));
  if (NULL == name) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindName: name was NULL"));
  if (NULL == index) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindName: output parameter was NULL"));

  for (structure = 0; structure < ct->nentries; structure++) {
    if (NULL != ct->entries[structure]) {
      if (stricmp(name, ct->entries[structure]->name) == 0) {
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
  CTABfindEntryByName() - given the string name of a structure, return the
  entry number( NOT index!) in the color table that matches the name
  -----------------------------------------------------------------*/
int CTABfindEntryByName(COLOR_TABLE *ct, const char *name, int *nEntry)
{
  int structure;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindName: ct was NULL"));
  if (NULL == name) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindName: name was NULL"));
  if (NULL == nEntry) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABfindName: output parameter was NULL"));

  int n = 0;
  for (structure = 0; structure < ct->nentries; structure++) {
    if (NULL != ct->entries[structure]) {
      if (stricmp(name, ct->entries[structure]->name) == 0) {
        *nEntry = n;
        return (NO_ERROR);
      }
      n++;
    }
  }

  /* If we got here, we didn't find a match. */
  *nEntry = -1;
  return (NO_ERROR);
}

/*--------------------------------------------------------------*/
int CTABprintASCII(COLOR_TABLE *ct, FILE *fp)
{
  int structure;
  char *tmpstr;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABprintASCII: ct was NULL"));
  if (NULL == fp) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABprintASCII: fp was NULL"));

  for (structure = 0; structure < ct->nentries; structure++) {
    if (NULL != ct->entries[structure]) {
      tmpstr = deblank(ct->entries[structure]->name);
      fprintf(fp,
              "%3d  %-30s  %3d %3d %3d  %3d\n",
              structure + ct->idbase,
              tmpstr,
              ct->entries[structure]->ri,
              ct->entries[structure]->gi,
              ct->entries[structure]->bi,
              255 - ct->entries[structure]->ai); /* alpha = 255-trans */
      free(tmpstr);
    }
  }

  return (ERROR_NONE);
}

/*--------------------------------------------------------------*/
int CTABwriteFileASCII(COLOR_TABLE *ct, const char *fname)
{
  FILE *fp;
  int result;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteFileASCII: ct was NULL"));
  if (NULL == fname) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABwriteFileASCII: fname was NULL"));

  fp = fopen(fname, "w");
  if (fp == NULL) {
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "CTABwriteFileASCII(%s): could not open for writing\n", fname));
  }

  result = CTABprintASCII(ct, fp);
  fclose(fp);

  return (result);
}

/*-------------------------------------------------------
  COLOR_TABLE *CTABaddEntry(COLOR_TABLE *ctold, char *name)
  Creates a new color table by adding name to input table.
  Colors are assigned randomly.
  -------------------------------------------------------*/
COLOR_TABLE *CTABaddEntry(COLOR_TABLE *ctold, const char *name)
{
  COLOR_TABLE *ct;
  COLOR_TABLE_ENTRY *cte;
  int nentries, i;

  nentries = ctold->nentries;
  ct = CTABalloc(nentries + 1);

  for (i = 0; i < nentries; i++) memmove(ct->entries[i], ctold->entries[i], sizeof(COLOR_TABLE_ENTRY));

  //    *(ct->entries[i]) = *(ctold->entries[i]) ;

  cte = ct->entries[nentries];
  sprintf(cte->name, "%s", name);
  cte->ri = floor(drand48() * 256);
  cte->gi = floor(drand48() * 256);
  cte->bi = floor(drand48() * 256);
  cte->rf = (float)cte->ri / 255.0f;
  cte->gf = (float)cte->gi / 255.0f;
  cte->bf = (float)cte->bi / 255.0f;
  printf("RGB %d %d %d\n", cte->ri, cte->gi, cte->bi);

  return (ct);
}

/*!
\fn COLOR_TABLE *TissueTypeSchema(COLOR_TABLE *ct, char *schema)
\brief Adds tissue type information to a color table using the
named schema. The added information includes (1) the tissue type
for each/most entries, (2) a second embedded color table that
describes the tissue classes, and (3) the schema name.
\param ct - color table with tissue type info (may or may not be null
depending upon the schema)
*/
COLOR_TABLE *TissueTypeSchema(COLOR_TABLE *ct, const char *schema)
{
  if(strcmp(schema, "default-jan-2014")==0){
    ct = TissueTypeSchemaDefault(ct);
    return (ct);
  }
  if(strcmp(schema, "default-jan-2014+head")==0 || strcmp(schema, "default-apr-2019+head")==0) {
    ct = TissueTypeSchemaDefaultHead(ct);
    return (ct);
  }
  printf("ERROR: tissue type schema %s unrecognized\n", schema);
  return (NULL);
}

/*!
\fn COLOR_TABLE *TissueTypeSchemaDefault(COLOR_TABLE *ct)
\brief Adds tissue type information to a color table using the
default FreeSurfer schema.
\param ct - color table with tissue type info (if null then
uses FreeSurferColorLUT.txt)
*/
COLOR_TABLE *TissueTypeSchemaDefault(COLOR_TABLE *ct)
{
  COLOR_TABLE_ENTRY *cte;
  FSENV *fsenv;
  char tmpstr[2000];
  int n, TT, TTUnknown, TTCtxGM, TTSubCtxGM, TTWM, TTCSF;
  TTUnknown = 0;
  TTCtxGM = 1;
  TTSubCtxGM = 2;
  TTWM = 3;
  TTCSF = 4;

  if (ct == NULL) {
    fsenv = FSENVgetenv();
    sprintf(tmpstr, "%s/FreeSurferColorLUT.txt", fsenv->FREESURFER_HOME);
    ct = CTABreadASCII(tmpstr);
  }

  sprintf(ct->TissueTypeSchema, "default-jan-2014");
  ct->ctabTissueType = CTABalloc(5);
  cte = ct->ctabTissueType->entries[0];
  sprintf(cte->name, "unknown");
  cte->ri = 0;
  cte->gi = 0;
  cte->bi = 0;
  cte = ct->ctabTissueType->entries[1];
  sprintf(cte->name, "cortex");
  cte->ri = 205;
  cte->gi = 62;
  cte->bi = 78;
  cte = ct->ctabTissueType->entries[2];
  sprintf(cte->name, "subcort_gm");
  cte->ri = 230;
  cte->gi = 148;
  cte->bi = 34;
  cte = ct->ctabTissueType->entries[3];
  sprintf(cte->name, "wm");
  cte->ri = 0;
  cte->gi = 255;
  cte->bi = 0;
  cte = ct->ctabTissueType->entries[4];
  sprintf(cte->name, "csf");
  cte->ri = 120;
  cte->gi = 18;
  cte->bi = 134;

  for (n = 0; n < ct->nentries; n++) {
    cte = ct->entries[n];
    if (cte == NULL) continue;

    TT = -1;
    switch (n) {
      case 0:  // unknown
        TT = TTUnknown;
        break;

      case Left_Cerebral_Cortex:
      case Right_Cerebral_Cortex:
        TT = TTCtxGM;
        break;

      case Left_Cerebellum_Cortex:
      case Right_Cerebellum_Cortex:
      case Left_Hippocampus:
      case Right_Hippocampus:
      case Left_Amygdala:
      case Right_Amygdala:
      case Left_Pallidum:
      case Right_Pallidum:
      case Left_Thalamus:
      case Right_Thalamus:
      case Right_Putamen:
      case Left_Putamen:
      case Right_Caudate:
      case Left_Caudate:
      case Left_Accumbens_area:
      case Right_Accumbens_area:
      case Left_VentralDC:
      case Right_VentralDC:
      case Brain_Stem: // was WM
      case Spinal_Cord: //( same as brainstem)
      case 174:  // Pons,  was WM
      case Left_undetermined:
      case Right_undetermined:
      case non_WM_hypointensities:  // not sure
        TT = TTSubCtxGM;
        break;

      case Left_Cerebral_White_Matter:
      case Right_Cerebral_White_Matter:
      case Left_Cerebellum_White_Matter:
      case Right_Cerebellum_White_Matter:
      case 690:  // Cerebellum CbmWM_Gyri_Left
      case 691:  // Cerebellum CbmWM_Gyri_Right
      case WM_hypointensities:
      case Left_WM_hypointensities:
      case Right_WM_hypointensities:
      case Optic_Chiasm:
      case 930: // Optic-Nerve
      case Corpus_Callosum:
      case CC_Posterior:
      case CC_Mid_Posterior:
      case CC_Central:
      case CC_Mid_Anterior:
      case CC_Anterior:
      case 34: // wmcrowns, lh
      case 66: // wmcrowns, rh
        TT = TTWM;
        break;

      case Third_Ventricle:
      case Fourth_Ventricle:
      case CSF:
      case CSF_ExtraCerebral:
      case Left_Lateral_Ventricle:
      case Right_Lateral_Ventricle:
      case Left_Inf_Lat_Vent:
      case Right_Inf_Lat_Vent:
      case Left_choroid_plexus:
      case Right_choroid_plexus:
      case Fifth_Ventricle:
      case Left_vessel:
      case Right_vessel:
      case 914: // Vein -- putting it here as it is typically labeled as CSF
        TT = TTCSF;
        break;
    }

    if (TT == -1) {
      // still not assigned
      if (n >= 1000 && n <= 1035) TT = TTCtxGM;
      if (n >= 2000 && n <= 2035) TT = TTCtxGM;
      if (n >= 3000 && n <= 3035) TT = TTWM;
      if (n >= 4000 && n <= 4035) TT = TTWM;
      if (n >= 1100 && n <= 1181) TT = TTCtxGM;
      if (n >= 2100 && n <= 2181) TT = TTCtxGM;
      if (n >= 3100 && n <= 3181) TT = TTWM;
      if (n >= 4100 && n <= 4181) TT = TTWM;
      if (n == 5001 || n == 5002) TT = TTWM;
      if (n >= 11100 && n <= 11300) TT = TTCtxGM;
      if (n >= 12100 && n <= 12300) TT = TTCtxGM;
    }

    cte->TissueType = TT;
  }

  return (ct);
}

/*!
\fn COLOR_TABLE *TissueTypeSchemaDefaultHead(COLOR_TABLE *ct)
\brief Adds tissue type information to a color table using the
default FreeSurfer schema and includes a Head tissue type.
Apr 2019, changed ventraldc, brainstem, and pons to be subcort GM
\param ct - color table with tissue type info (if null then
uses FreeSurferColorLUT.txt)
*/
COLOR_TABLE *TissueTypeSchemaDefaultHead(COLOR_TABLE *ct)
{
  COLOR_TABLE_ENTRY *cte;
  FSENV *fsenv;
  char tmpstr[2000];
  int n, TT, TTUnknown, TTCtxGM, TTSubCtxGM, TTWM, TTCSF, TTHead;
  TTUnknown = 0;
  TTCtxGM = 1;
  TTSubCtxGM = 2;
  TTWM = 3;
  TTCSF = 4;
  TTHead = 5;

  printf("Entering TissueTypeSchemaDefaultHead()\n");

  if (ct == NULL) {
    fsenv = FSENVgetenv();
    sprintf(tmpstr, "%s/FreeSurferColorLUT.txt", fsenv->FREESURFER_HOME);
    ct = CTABreadASCII(tmpstr);
  }

  // changed schema name to help indicate that ventraldc, brainstem,
  // and pons are now subcortical GM
  sprintf(ct->TissueTypeSchema, "default-apr-2019+head"); 
  printf("schema %s\n",ct->TissueTypeSchema);
  ct->ctabTissueType = CTABalloc(6);
  cte = ct->ctabTissueType->entries[0];
  sprintf(cte->name, "unknown");
  cte->ri = 0;
  cte->gi = 0;
  cte->bi = 0;
  cte = ct->ctabTissueType->entries[1];
  sprintf(cte->name, "cortex");
  cte->ri = 205;
  cte->gi = 62;
  cte->bi = 78;
  cte = ct->ctabTissueType->entries[2];
  sprintf(cte->name, "subcort_gm");
  cte->ri = 230;
  cte->gi = 148;
  cte->bi = 34;
  cte = ct->ctabTissueType->entries[3];
  sprintf(cte->name, "wm");
  cte->ri = 0;
  cte->gi = 255;
  cte->bi = 0;
  cte = ct->ctabTissueType->entries[4];
  sprintf(cte->name, "csf");
  cte->ri = 120;
  cte->gi = 18;
  cte->bi = 134;
  cte = ct->ctabTissueType->entries[5];
  sprintf(cte->name, "head");
  cte->ri = 150;
  cte->gi = 150;
  cte->bi = 200;

  for (n = 0; n < ct->nentries; n++) {
    cte = ct->entries[n];
    if (cte == NULL) continue;

    TT = -1;
    switch (n) {
      case 0:  // unknown
        TT = TTUnknown;
        break;

      case Left_Cerebral_Cortex:
      case Right_Cerebral_Cortex:
        TT = TTCtxGM;
        break;

      case Left_Cerebellum_Cortex:
      case Right_Cerebellum_Cortex:
      case 172:  // Vermis
      case 179:  // Floculus
      case Left_Hippocampus:
      case Right_Hippocampus:
      case Left_Amygdala:
      case Right_Amygdala:
      case Left_Pallidum:
      case Right_Pallidum:
      case Left_Thalamus:
      case Right_Thalamus:
      case Right_Putamen:
      case Left_Putamen:
      case Right_Caudate:
      case Left_Caudate:
      case Left_Accumbens_area:
      case Right_Accumbens_area:
      case Left_VentralDC:
      case Right_VentralDC:
      case Brain_Stem: // was WM
      case Spinal_Cord: //( same as brainstem)
      case 174:  // Pons,  was WM
      case 267:  // Pons-Belly-Area
      case non_WM_hypointensities:  // not sure
      case Left_undetermined:
      case Right_undetermined:
        TT = TTSubCtxGM;
        break;

      case Left_Cerebral_White_Matter:
      case Right_Cerebral_White_Matter:
      case Left_Cerebellum_White_Matter:
      case Right_Cerebellum_White_Matter:
      case 690:  // Cerebellum CbmWM_Gyri_Left
      case 691:  // Cerebellum CbmWM_Gyri_Right
      case WM_hypointensities:
      case Left_WM_hypointensities:
      case Right_WM_hypointensities:
      case Optic_Chiasm:
      case 930: // Optic-Nerve
      case Corpus_Callosum:
      case CC_Posterior:
      case CC_Mid_Posterior:
      case CC_Central:
      case CC_Mid_Anterior:
      case CC_Anterior:
      case 34: // wmcrowns, lh
      case 66: // wmcrowns, rh
        TT = TTWM;
        break;

      case Third_Ventricle:
      case Fourth_Ventricle:
      case CSF:
      case CSF_ExtraCerebral:
      case Left_Lateral_Ventricle:
      case 75: // was Left_Lateral_Ventricle
      case Right_Lateral_Ventricle:
      case 76: // was Right_Lateral_Ventricle
      case Left_Inf_Lat_Vent:
      case Right_Inf_Lat_Vent:
      case Left_choroid_plexus:
      case Right_choroid_plexus:
      case Fifth_Ventricle:
      case Left_vessel:
      case Right_vessel:
      case 914: // Vein -- putting it here as it is typically labeled as CSF
        TT = TTCSF;
        break;

      case Head_ExtraCerebral:
      case 165:  // Skull
      case 259:  // PossibleSkull
      case 118:  // Epidermis
      case 119:  // Conn-Tissue
      case 120:  // SC-Fat-Muscle
      case 121:  // Cranium
      case 122:  // CSF-SA
      case 123:  // Muscle
      case 124:  // Ear
      case 127:  // Soft-Tissue
      case 129:  // Bone
      case 130:  // Air
      case 131:  // Oribital-Fat
      case 132:  // Tongue
      case 133:  // Nasal-Structures
      case 134:  // Globe
      case 135:  // Teeth
      case 143:  // Vitreous-Humor
      case 144:  // Lens
      case 145:  // Atrieous-Humor
      case 167:  // Scalp
      case 262:  // Sinus
      case 263:  // Left-Eustachian
      case 264:  // Right-Eustachian
      case 902:  //Artery
      case 907:  //Other Tissue
      case 908:  //Eye muscles
      case 909:  //Mucosa
      case 911:  //Skin
      case 915:  //Bone-Cortical
      case 916:  //Bone-Cancellous
        TT = TTHead;
        break;
    }

    if (TT == -1) {
      // still not assigned
      if (n >= 1000 && n <= 1035) TT = TTCtxGM;
      if (n >= 2000 && n <= 2035) TT = TTCtxGM;
      if (n >= 3000 && n <= 3035) TT = TTWM;
      if (n >= 4000 && n <= 4035) TT = TTWM;
      if (n >= 3201 && n <= 3207) TT = TTWM;
      if (n >= 4201 && n <= 4207) TT = TTWM;
      if (n >= 1100 && n <= 1181) TT = TTCtxGM;
      if (n >= 2100 && n <= 2181) TT = TTCtxGM;
      if (n >= 3100 && n <= 3181) TT = TTWM;
      if (n >= 4100 && n <= 4181) TT = TTWM;
      if (n == 5001 || n == 5002) TT = TTWM;
      if (n >= 11100 && n <= 11300) TT = TTCtxGM;
      if (n >= 12100 && n <= 12300) TT = TTCtxGM;
    }

    cte->TissueType = TT;
  }

  return (ct);
}

/*!
\fn COLOR_TABLE *TissueTypeSchemaLat(COLOR_TABLE *ct)
\brief Adds tissue type information to a color table using the
default FreeSurfer schema including a Head tissue type. 
Creates lateralized tissue types.
\param ct - color table with tissue type info (if null then
uses FreeSurferColorLUT.txt)
*/
COLOR_TABLE *TissueTypeSchemaLat(COLOR_TABLE *ct)
{
  COLOR_TABLE_ENTRY *cte;
  FSENV *fsenv;
  char tmpstr[2000];
  int n, TT, TTUnknown, TTCtxGMlh, TTSubCtxGMlh, TTSubCtxGMmid, TTCtxGMrh, TTSubCtxGMrh, TTWM, TTCSF, TTHead;
  TTUnknown = 0;
  TTCtxGMlh = 1;
  TTCtxGMrh = 2;
  TTSubCtxGMlh = 3;
  TTSubCtxGMrh = 4;
  TTSubCtxGMmid = 5;
  TTWM = 6;
  TTCSF = 7;
  TTHead = 8;

  if (ct == NULL) {
    fsenv = FSENVgetenv();
    sprintf(tmpstr, "%s/FreeSurferColorLUT.txt", fsenv->FREESURFER_HOME);
    ct = CTABreadASCII(tmpstr);
  }

  sprintf(ct->TissueTypeSchema, "default-apr-2019+head+lat");
  ct->ctabTissueType = CTABalloc(9); // should dealloc existing?
  cte = ct->ctabTissueType->entries[0];
  sprintf(cte->name, "unknown");
  cte->ri = 0;
  cte->gi = 0;
  cte->bi = 0;
  cte = ct->ctabTissueType->entries[1];
  sprintf(cte->name, "cortex-lh");
  cte->ri = 205;
  cte->gi = 62;
  cte->bi = 78;
  cte = ct->ctabTissueType->entries[2];
  sprintf(cte->name, "cortex-rh");
  cte->ri = 205;
  cte->gi = 62;
  cte->bi = 78;
  cte = ct->ctabTissueType->entries[3];
  sprintf(cte->name, "subcort_gm-lh");
  cte->ri = 230;
  cte->gi = 148;
  cte->bi = 34;
  cte = ct->ctabTissueType->entries[4];
  sprintf(cte->name, "subcort_gm-rh");
  cte->ri = 230;
  cte->gi = 148;
  cte->bi = 34;
  cte = ct->ctabTissueType->entries[5];
  sprintf(cte->name, "subcort_gm-mid");
  cte->ri = 230;
  cte->gi = 148;
  cte->bi = 34;
  cte = ct->ctabTissueType->entries[6];
  sprintf(cte->name, "wm");
  cte->ri = 0;
  cte->gi = 255;
  cte->bi = 0;
  cte = ct->ctabTissueType->entries[7];
  sprintf(cte->name, "csf");
  cte->ri = 120;
  cte->gi = 18;
  cte->bi = 134;
  cte = ct->ctabTissueType->entries[8];
  sprintf(cte->name, "head");
  cte->ri = 150;
  cte->gi = 150;
  cte->bi = 200;

  for (n = 0; n < ct->nentries; n++) {
    cte = ct->entries[n];
    if (cte == NULL) continue;

    TT = -1;
    switch (n) {
      case 0:  // unknown
        TT = TTUnknown;
        break;

      case Left_Cerebral_Cortex:
        TT = TTCtxGMlh;
        break;
      case Right_Cerebral_Cortex:
        TT = TTCtxGMrh;
        break;

      case Left_Cerebellum_Cortex:
      case Left_Hippocampus:
      case Left_Amygdala:
      case Left_Thalamus:
      case Left_Putamen:
      case Left_Pallidum:
      case Left_Caudate:
      case Left_Accumbens_area:
      case Left_VentralDC:
      case 179:  // Floculus
      case 183:  // Left-Vermis
      case non_WM_hypointensities:  // not sure
      case Left_undetermined:
        TT = TTSubCtxGMlh;
        break;

      case Right_Cerebellum_Cortex:
      case Right_Hippocampus:
      case Right_Amygdala:
      case Right_Thalamus:
      case Right_Putamen:
      case Right_Pallidum:
      case Right_Caudate:
      case Right_Accumbens_area:
      case Right_VentralDC:
      case Right_undetermined:
      case 184:  // Right-Vermis
        TT = TTSubCtxGMrh;
        break;

      case 172:  // Vermis // not clear how to laterlize
      case 174:  // Pons (now considered GM)
      case 267:  // Pons-Belly-Area
      case Brain_Stem: //(now considered GM)
      case Spinal_Cord: //( same as brainstem)
        TT = TTSubCtxGMmid;
        break;

      case Left_Cerebral_White_Matter:
      case Right_Cerebral_White_Matter:
      case Left_Cerebellum_White_Matter:
      case Right_Cerebellum_White_Matter:
      case 690:  // Cerebellum CbmWM_Gyri_Left
      case 691:  // Cerebellum CbmWM_Gyri_Right
      case WM_hypointensities:
      case Left_WM_hypointensities:
      case Right_WM_hypointensities:
      case Optic_Chiasm:
      case 930: // Optic-Nerve
      case Corpus_Callosum:
      case CC_Posterior:
      case CC_Mid_Posterior:
      case CC_Central:
      case CC_Mid_Anterior:
      case CC_Anterior:
      case 34: // wmcrowns, lh
      case 66: // wmcrowns, rh
        TT = TTWM;
        break;

      case Third_Ventricle:
      case Fourth_Ventricle:
      case CSF:
      case CSF_ExtraCerebral:
      case Left_Lateral_Ventricle:
      case Right_Lateral_Ventricle:
      case Left_Inf_Lat_Vent:
      case Right_Inf_Lat_Vent:
      case Left_choroid_plexus:
      case Right_choroid_plexus:
      case Fifth_Ventricle:
      case Left_vessel:
      case Right_vessel:
        TT = TTCSF;
        break;

      case Head_ExtraCerebral:
      case 165:  // Skull
      case 259:  // PossibleSkull
      case 118:  // Epidermis
      case 119:  // Conn-Tissue
      case 120:  // SC-Fat-Muscle
      case 121:  // Cranium
      case 122:  // CSF-SA
      case 123:  // Muscle
      case 124:  // Ear
      case 127:  // Soft-Tissue
      case 129:  // Bone
      case 130:  // Air
      case 131:  // Oribital-Fat
      case 132:  // Tongue
      case 133:  // Nasal-Structures
      case 134:  // Globe
      case 135:  // Teeth
      case 143:  // Vitreous-Humor
      case 144:  // Lens
      case 145:  // Atrieous-Humor
      case 167:  // Scalp
      case 262:  // Sinus
      case 263:  // Left-Eustachian
      case 264:  // Right-Eustachian
      case 902:  //Artery
      case 907:  //Other Tissue
      case 908:  //Eye muscles
      case 909:  //Mucosa
      case 911:  //Skin
      case 915:  //Bone-Cortical
      case 916:  //Bone-Cancellous
        TT = TTHead;
        break;
    }

    if (TT == -1) {
      // still not assigned
      if (n >= 1000 && n <= 1035) TT = TTCtxGMlh;
      if (n >= 2000 && n <= 2035) TT = TTCtxGMrh;
      if (n >= 3000 && n <= 3035) TT = TTWM;
      if (n >= 4000 && n <= 4035) TT = TTWM;
      if (n >= 3201 && n <= 3207) TT = TTWM;
      if (n >= 4201 && n <= 4207) TT = TTWM;
      if (n >= 1100 && n <= 1181) TT = TTCtxGMlh;
      if (n >= 2100 && n <= 2181) TT = TTCtxGMrh;
      if (n >= 3100 && n <= 3181) TT = TTWM;
      if (n >= 4100 && n <= 4181) TT = TTWM;
      if (n == 5001 || n == 5002) TT = TTWM;
      if (n >= 11100 && n <= 11300) TT = TTCtxGMlh;
      if (n >= 12100 && n <= 12300) TT = TTCtxGMrh;
    }

    cte->TissueType = TT;
  }

  return (ct);
}

/*--------------------------------------------------------------*/
/*!
\fn int CTABprintASCIItt(COLOR_TABLE *ct, FILE *fp)
\brief Prints the color table including tissue type information.
with ribbon values if the aseg is CtxGM or CtxWM or unknown.
\param ct - color table with tissue type info
*/
int CTABprintASCIItt(COLOR_TABLE *ct, FILE *fp)
{
  int structure;
  char *tmpstr;
  COLOR_TABLE_ENTRY *cte;

  if (NULL == ct) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABprintASCIItt: ct was NULL"));
  if (NULL == fp) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABprintASCIItt: fp was NULL"));

  if (ct->ctabTissueType == NULL)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "CTABprintASCIItt: tissue type ctab was NULL"));

  fprintf(fp, "# TissueTypeSchema %s\n", ct->TissueTypeSchema);
  for (structure = 0; structure < ct->ctabTissueType->nentries; structure++) {
    cte = ct->ctabTissueType->entries[structure];
    if (cte == NULL) continue;
    tmpstr = deblank(cte->name);
    fprintf(fp,
            "#ctTType %3d  %-30s  %3d %3d %3d  %3d\n",
            structure + ct->idbase,
            tmpstr,
            cte->ri,
            cte->gi,
            cte->bi,
            255 - cte->ai);
    free(tmpstr);
  }

  for (structure = 0; structure < ct->nentries; structure++) {
    cte = ct->entries[structure];
    if (cte == NULL) continue;
    // if(cte->TissueType == -1) continue;
    tmpstr = deblank(cte->name);
    fprintf(fp,
            "%3d  %-30s  %3d %3d %3d  %3d  %2d\n",
            structure + ct->idbase,
            tmpstr,
            cte->ri,
            cte->gi,
            cte->bi,
            255 - cte->ai,
            cte->TissueType);
    free(tmpstr);
  }

  return (ERROR_NONE);
}
/*--------------------------------------------------------------*/
/*!
\fn int CTABwriteFileASCIItt(COLOR_TABLE *ct, const char *fname)
\brief Writes the color table including tissue type information.
with ribbon values if the aseg is CtxGM or CtxWM or unknown.
\param ct - color table with tissue type info
*/
int CTABwriteFileASCIItt(COLOR_TABLE *ct, const char *fname)
{
  FILE *fp = fopen(fname, "w");
  if (fp == NULL) {
    printf("ERROR: could not open %s for writing\n", fname);
    return (1);
  }
  CTABprintASCIItt(ct, fp);
  fclose(fp);
  return (0);
}

/*--------------------------------------------------------------*/
/*!
\fn int CTABmerge(COLOR_TABLE *ct, const COLOR_TABLE *merge)
\brief Takes items in merge and copies them into ct. If ct
already has an item for that structure, it is deleted and
the new one is used to overwrite it.
*/
int CTABmerge(COLOR_TABLE *ct, const COLOR_TABLE *merge)
{
  int n;
  CTE *cte, *cte0;

  for (n = 0; n < merge->nentries; n++) {
    cte = merge->entries[n];
    if (cte == NULL) continue;
    cte0 = ct->entries[n];
    if (cte0 == NULL) {
      cte0 = (CTE *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
      ct->entries[n] = cte0;
    }
    memcpy(cte0, cte, sizeof(CTE));
  }
  return (0);
}

static int ctabMinDist(COLOR_TABLE *ct, int r, int g, int b)
{
  int i, dist, min_dist = 3 * 256;

  for (i = 0; i < ct->nentries; i++) {
    if (ct->entries[i] == NULL) continue;
    dist = abs(ct->entries[i]->ri - r) + abs(ct->entries[i]->gi - g) + abs(ct->entries[i]->bi - b);
    if (dist < min_dist) min_dist = dist;
  }
  return (min_dist);
}

/*--------------------------------------------------------------*/
/*!
\fn int CTABaddUniqueEntry(COLOR_TABLE *ct, char *name, int min_dist)
\brief Creates and adds a unique entry into the ct (increasing ct->nentries)
The new entry will be at least min_dist from any existing rgb value so it can be
visually distinguished (rdist+gdist+bdist)
*/
int CTABaddUniqueEntry(COLOR_TABLE *ct, char *name, int min_dist)
{
  int dist, i, r, g, b;
  COLOR_TABLE_ENTRY *cte, **pcte;

  while (min_dist > 0) {
    for (i = 0; i < 1000; i++) {
      r = nint(randomNumber(0, 255));
      g = nint(randomNumber(0, 255));
      b = nint(randomNumber(0, 255));
      dist = ctabMinDist(ct, r, g, b);
      if (dist <= min_dist) break;
    }
    if (dist > min_dist)
      min_dist--;
    else
      break;
  }

  if (min_dist <= 0) return (-1);

  for (i = 0; i < ct->nentries; i++)  // see if there are any unused slots
    if (ct->entries[i] == NULL) {
      ct->entries[i] = (COLOR_TABLE_ENTRY *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
      break;
    }

  if (i >= ct->nentries)  // allocate and copy over new table
  {
    pcte = ct->entries;
#if GCC_VERSION > 80000
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Walloc-size-larger-than="
#endif
    ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY *));
#if GCC_VERSION > 80000
#pragma GCC diagnostic pop
#endif
    for (i = 0; i < ct->nentries; i++) {
      ct->entries[i] = (COLOR_TABLE_ENTRY *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
      memmove(ct->entries[i], pcte[i], sizeof(COLOR_TABLE_ENTRY));
    }
    i = ct->nentries;
    ct->nentries++;
  }

  cte = ct->entries[i];
  cte->ri = r;
  cte->gi = g;
  cte->bi = b;
  cte->rf = (float)cte->ri / 255.0f;
  cte->gf = (float)cte->gi / 255.0f;
  cte->bf = (float)cte->bi / 255.0f;
  cte->ai = 255;
  strcpy(cte->name, name);

  return (i);
}


/* add an unique entry at the end of COLOR_TABLE
 * the ctab index of newly added entry is return in *ctabindex
 */
int CTABaddUniqueEntryAtEnd(COLOR_TABLE *ct, char *name, int *ctabindex)
{
  *ctabindex = -1;

  COLOR_TABLE_ENTRY **pcte = ct->entries;

  // allocate one more COLOR_TABLE_ENTRY
  ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries+1, sizeof(COLOR_TABLE_ENTRY *));

  for (int i = 0; i < ct->nentries; i++)
  {
    // point to old memory for each COLOR_TABLE_ENTRY
    ct->entries[i] = pcte[i];
#if 0
    // copy over the data, free the old memory for each COLOR_TABLE_ENTRY
    if (pcte[i] != NULL)
    {
      ct->entries[i] = (COLOR_TABLE_ENTRY *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
      memmove(ct->entries[i], pcte[i], sizeof(COLOR_TABLE_ENTRY));
      free(pcte[i]);
    }
#endif
  }
  
  int r, g, b;
  int rgbi = -1;
  do
  { 
    r = nint(randomNumber(0, 255));
    g = nint(randomNumber(0, 255));
    b = nint(randomNumber(0, 255));

    CTABfindRGBi(ct, r, g, b, &rgbi);
  } while (rgbi >= 0);

  COLOR_TABLE_ENTRY *cte = (COLOR_TABLE_ENTRY *)calloc(1, sizeof(COLOR_TABLE_ENTRY));
  cte->ri = r;
  cte->gi = g;
  cte->bi = b;
  cte->ai = 255;
  cte->rf = (float)cte->ri / 255.0f;
  cte->gf = (float)cte->gi / 255.0f;
  cte->bf = (float)cte->bi / 255.0f;
  cte->af = (float)cte->ai / 255.0f;
  strcpy(cte->name, name);

  ct->entries[ct->nentries] = cte;
  *ctabindex = ct->nentries;
  ct->nentries++;

  return NO_ERROR;
}


// generate unique annotation the given COLOR_TABLE
int CTABgenUniqueAnnotation(const COLOR_TABLE *ctab1, const COLOR_TABLE *ctab2)
{
  int ri, gi, bi;
  int rgbi_1 = -1, rgbi_2 = -1;
  do
  { 
    ri = nint(randomNumber(0, 255));
    gi = nint(randomNumber(0, 255));
    bi = nint(randomNumber(0, 255));

    CTABfindRGBi((COLOR_TABLE *)ctab1, ri, gi, bi, &rgbi_1);
    if (ctab2 != NULL)
      CTABfindRGBi((COLOR_TABLE *)ctab1, ri, gi, bi, &rgbi_2);
  } while (rgbi_1 >= 0 || rgbi_2 >= 0);

  int annotation;
  RGBToAnnot(ri, gi, bi, annotation);

  return annotation;
}

// if suggestAnnot is true,
// duplicated annotation assignment in ctab will be replaced with newly calculated suggestion
int CTABprintAnnotationAssignment(COLOR_TABLE *ctab, bool suggestAnnot, FILE *fp, const char *exception)
{
  int dupCount = 0;
  
  if (fp == NULL)
    fp = stdout;

  if (suggestAnnot)
    setRandomSeed(12);

  // unique_annot_map is sorted
  //   key   = annotation
  //   value = label-id vector  
  std::map<int, std::vector<int>> unique_annot_map;
  for (int numentry = 0; numentry < ctab->nentries; numentry++)
  {
    if (ctab->entries[numentry] != NULL)
    {
      int r = ctab->entries[numentry]->ri;
      int g = ctab->entries[numentry]->gi;
      int b = ctab->entries[numentry]->bi;

      /* Duplicated macros are defined 
       *   in mrisurf.h as MRISAnnotToRGB and MRISRGBToAnnot, and
       *   in colortable.h as AnnotToRGB and RGBToAnnot.
       * annotation => RGB are also calculated in CTABfindAnnotation(), CTABgetAnnotationName()
       * RGB => annotation are also calculated in CTABrgb2Annotation(), CTABannotationAtIndex()
       */
      int annotation;  //int annotation = CTABrgb2Annotation(r, g, b);
      RGBToAnnot(r, g, b, annotation);
      unique_annot_map[annotation].push_back(numentry);
    }
  }

  // print annotation and all label-ids assinged
  // suggest new annotation if it is requested
  std::map<int, std::vector<int>>::const_iterator annotIt = unique_annot_map.begin();
  for (annotIt = unique_annot_map.begin(); annotIt != unique_annot_map.end(); annotIt++)
  {
    int r, g, b;

    int annot = annotIt->first; 
    AnnotToRGB(annot, r, g, b);
    int nlabels = annotIt->second.size();
    fprintf(fp, "%8d (%-3d) :", annot, nlabels); // left-justified

    int first = 1;
    std::vector<int>::const_iterator labelIt = annotIt->second.begin();
    for (; labelIt != annotIt->second.end(); labelIt++)
    {
      const char *labelname = ctab->entries[*labelIt]->name;
      if (exception != NULL && strncmp(labelname, exception, strlen(exception)) == 0 && !first)
      {
	fprintf(fp, "%16s %-6d (%-3d %-3d %-3d) (%8d) %s\n", " ", *labelIt, r, g, b, annot, labelname);
	continue;
      }
      
      if (first || !suggestAnnot)
      {
	if (first)
	{
          first = 0;
          fprintf(fp, " %-6d (%-3d %-3d %-3d) (%8d) %s\n", *labelIt, r, g, b, annot, labelname);
	}
	else
	  fprintf(fp, "%16s %-6d (%-3d %-3d %-3d) (%8d) %s\n", " ", *labelIt, r, g, b, annot, labelname);
      }
      else
      {
	dupCount++;

        // these are label ids assigned to the same annotation
        // suggest a different annotation
	int ri, gi, bi;
	annot = CTABgenUniqueAnnotation(ctab);
	AnnotToRGB(annot, ri, gi, bi);
        fprintf(fp, "%15s *%-6d (%-3d %-3d %-3d) (%8d) %s\n", " ", *labelIt, ri, gi, bi, annot, labelname);

	// replace the entry with newly calculated annotation
	ctab->entries[*labelIt]->ri = ri;
	ctab->entries[*labelIt]->gi = gi;
	ctab->entries[*labelIt]->bi = bi;

	// the float versions
	ctab->entries[*labelIt]->rf = (float)ri / 255.0;
	ctab->entries[*labelIt]->gf = (float)gi / 255.0;
	ctab->entries[*labelIt]->bf = (float)bi / 255.0;
      }
    }
  }

  return dupCount;
}
