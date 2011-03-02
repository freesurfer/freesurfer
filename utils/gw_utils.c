/**
 * @file  gw_utils.c
 * @brief miscellaneous utility functions contributed by Graham Wideman
 *
 */
/*
 * Original Author: Graham Wideman
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
 *    $Revision: 1.5 $
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


#ifndef Darwin
#include <values.h> // MAXSHORT
#endif
#ifndef MAXSHORT
#define MAXSHORT 32767
#endif
#include <time.h>

#include "error.h"
#include "gw_utils.h"

/*-----------------------------------------------
  GWU_make_surface_from_lists
  Adapted from ic642_make_surface()

  Assumes:
  -- Zero-based vertex and face numbering

  -------------------------------------------------*/
MRI_SURFACE *GWU_make_surface_from_lists(GWUTILS_VERTEX * vertices,
                                         int vertexcount,
                                         GWUTILS_FACE   * faces,
                                         int facecount )
{
  MRI_SURFACE *mris ;
  int         vno, fno, n, vn, n1, n2 ;
  VERTEX      *v = NULL;
  FACE        *f ;

  mris = MRISoverAlloc(0,0, vertexcount, facecount) ;

  //-----------------------------------------
  // Read vertex data into mris
  //-----------------------------------------
  for (vno = 0 ; vno < vertexcount ; vno++)
  {
    v = &mris->vertices[vno] ;

    v->x = vertices[vno].x ;
    v->y = vertices[vno].y ;
    v->z = vertices[vno].z ;
  }

  //-----------------------------------------
  // Read face data into mris, and count
  // # of faces each vertex is part of
  //-----------------------------------------
  for (fno = 0 ; fno < facecount ; fno++)
  {
    f = &mris->faces[fno] ;

    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      f->v[n] = faces[fno].vno[n] ;  // already zero-based
      v = &mris->vertices[f->v[n]] ;
      v->num++ ;
      v->vnum += 2 ;   // Not sure what the extra +2 is for...
    }
  }

  //-----------------------------------------
  // Allocate space for neighbor vertices
  // v->v is neighbor vertex list; v->vnum is count
  //-----------------------------------------
  for (vno = 0 ; vno < vertexcount ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->v = (int *)calloc(v->vnum/2, sizeof(int)) ;
    if (!v->v)
      ErrorExit(ERROR_NOMEMORY,
                "%s: could not allocate %dth vertex list.", __func__, vno) ;
    v->vnum = 0 ;
  }

  //-----------------------------------------
  // For each face,
  //   for each vertex
  //     Tell vertex that other vertices in face are its neighbors
  //-----------------------------------------
  for (fno = 0 ; fno < facecount ; fno++)
  {
    f = &mris->faces[fno] ;

    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      v = &mris->vertices[f->v[n]] ;

      // [sic] now add an edge to other 2 vertices if not already in list
      // [GW] Ie: tell each vertex about its neighbors from this face
      for (n1 = 0 ; n1 < VERTICES_PER_FACE ; n1++)
      {
        if (n1 == n)               // don't connect vertex to itself
          continue ;
        vn = faces[fno].vno[n1] ;  // already zero-based

        // now check to make sure it's not a duplicate
        for (n2 = 0 ; n2 < v->vnum ; n2++)
        {
          if (v->v[n2] == vn)
          {
            vn = -1 ; // mark it as a duplicate
            break ;
          }
        }

        if (vn >= 0)   // add only non-duplicates
          v->v[v->vnum++] = vn ;
      }
    }
  }

  //--------------------------------------------
  // In each vertex, allocate face array
  //--------------------------------------------
  for (vno = 0 ; vno < vertexcount ; vno++)
  {
    v->vtotal = v->vnum ;
    v = &mris->vertices[vno] ;
    v->f = (int *)calloc(v->num, sizeof(int)) ;
    if (!v->f)
      ErrorExit(ERROR_NO_MEMORY,
                "%s: could not allocate %d faces",__func__,v->num);
    v->n = (uchar *)calloc(v->num,sizeof(uchar));
    if (!v->n)
      ErrorExit(ERROR_NO_MEMORY,
                "%s: could not allocate %d nbrs",__func__, v->n);
    v->num = 0 ;   /* for use as counter in next section */
    v->dist = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist)
      ErrorExit(ERROR_NOMEMORY,
                "%s: could not allocate list of %d ""dists at v=%d",
                __func__, v->vnum, vno) ;
    v->dist_orig = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist_orig)
      ErrorExit(ERROR_NOMEMORY,
                "%s: could not allocate list of %d ""dists at v=%d",
                __func__, v->vnum, vno) ;
    v->nsize = 1 ;
    v->vtotal = v->vnum ;
  }

  //---------------------------------------------
  // Tell each vertex what faces it is part of
  //---------------------------------------------
  for (fno = 0 ; fno < facecount ; fno++)
  {
    f = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      v = &mris->vertices[f->v[n]] ;
      v->n[v->num] = n ;
      v->f[v->num++] = fno ;
    }
  }

  //---------------------------------------------
  // Final housekeeping
  //---------------------------------------------
  MRIScomputeMetricProperties(mris) ;

  mris->type = MRIS_ICO_SURFACE ;
  MRISsetNeighborhoodSize(mris, 1) ;
  return(mris) ;

}

/*-------------------------------------------------------
  Translate MHT data to an MRI volume in various different ways so
  that it can be visualized.
  Does not attempt to be a proper spatial transform to MRI space,
  (because what to do when mht->vres is large?)
  just a way to read out the MHT data in 3D
  ---------------------------------------------------------*/
MRI * MRIFromMHTandMRIS(MHT * mht, MRIS * mris, MFMM_Option_t mfmm_option)
{
  MRI * amri;
  int mhtvx, mhtvy, mhtvz;       // MHT "voxels"
  int mrix , mriy , mriz ;       // MRI voxels
  MHBT * bucket;
  MHB  * bin;
  int binnum;
  int fno_usage;   //
  int outval;

#define HALFMHTFOV 200
#define HALFMRIFOV 128

  fno_usage = mht->fno_usage;

  amri = MRIalloc(256, 256, 256, MRI_SHORT);
  for (mriz = 0; mriz < 255; mriz++) {
    for (mriy = 0; mriy < 255; mriy++) {
      for (mrix = 0; mrix < 255; mrix++) {
        MRISvox(amri, mrix, mriy, mriz) = 100;
      }
    }
  }
//  goto done;

  for (mriz = 0; mriz < 255; mriz++) {
    mhtvz = HALFMHTFOV - HALFMRIFOV + mriz;
    for (mriy = 0; mriy < 255; mriy++) {
      mhtvy = HALFMHTFOV - HALFMRIFOV + mriy;
      for (mrix = 0; mrix < 255; mrix++) {
        mhtvx = HALFMHTFOV - HALFMRIFOV + mrix;

        bucket = MHTgetBucketAtVoxIx(mht, mhtvx, mhtvy, mhtvz);

        outval = 0;

        if (bucket) {
          if (MFMM_None == mfmm_option) {
            outval = 1;
            goto outval_done;
          }
          for (binnum = 0; binnum < bucket->nused; binnum++ ){
            bin = &(bucket->bins[binnum]);

            switch (mfmm_option) {
            case MFMM_None: break;
            case MFMM_Num:
              outval = bin->fno;
              goto outval_done;
            case MFMM_NumDiv16:
              outval = bin->fno >> 4;
              goto outval_done;
            case MFMM_Count:
              outval++;
              break;
            }
          }
        }
      outval_done:
        if (outval > MAXSHORT)
          outval = MAXSHORT;

        // MRI?vox is a type-specific macro!
        MRISvox(amri, mrix, mriy, mriz) = outval;

      } // for mrix
    } // for mriy
  } // for mriz
//done:
  return amri;
}

// Some simple log functions for test programs.

static char local_Progname[500]     = "uninitialized";
static char local_Progversion[100]  = "uninitialized";
static char local_Logfilepath[1000] = "uninitialized";

//----------------------------------------
int gw_log_init(char * AProgname,
                char * AProgversion,
                char * ALogfilepath,
                int newfile) { // 0 for OK
//----------------------------------------
  FILE * afile;
  int rslt = 0;
  strcpy(local_Progname    , AProgname    );
  strcpy(local_Progversion , AProgversion );
  strcpy(local_Logfilepath , ALogfilepath );
  if (newfile) {
    afile = fopen(local_Logfilepath, "w");
  } else {
    afile = fopen(local_Logfilepath, "a");
  }

  if (afile) {
    fclose(afile);
  } else {
    rslt = 1;
  }
  return rslt;
}

//------------------------------
void gw_log_message(char * msg) {
//------------------------------
  FILE * afile;
  afile = fopen(local_Logfilepath, "a");
  fprintf(afile, "%s\n", msg);
  fclose(afile);
}

//------------------------------
static void nowstr(char * buf) {
//------------------------------
  time_t tim;
  struct tm *tmr;
  int rslt;

  time(&tim);
  tmr = localtime(&tim);
  rslt = strftime(buf, 100, "%Y-%m-%d %H:%M:%S", tmr);
}

//------------------------------
void gw_log_timestamp(char * label) {
//------------------------------
  char datestr[100];
  char msg[200];

  nowstr(datestr);
  sprintf(msg, "---[%s]--- %s version %s at %s",
          label, local_Progname, local_Progversion, datestr);
  gw_log_message(msg);
}

//------------------------------
void gw_log_begin(void) {
//------------------------------
  gw_log_timestamp("Begin");
}

//------------------------------
void gw_log_end(void) {
//------------------------------
  gw_log_timestamp("End");
}

