/**
 * @brief API for dijkstra related processing.
 *
 * Provides an API for dijkstra search through freesurfer structures.
 */
/*
 * Original Author: Rudolph Pienaar / Christian Haselgrove
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


#include <stdlib.h>
#include <assert.h>

#include "dijkstra.h"

struct d_node *d_list = NULL;

int addToList(MRIS *surf, int vno) {

  struct d_node *dn_new, *dn_this;

  dn_new = (struct d_node *)malloc(sizeof(struct d_node));
  if (dn_new == NULL) {
    ErrorReturn(ERROR_NO_MEMORY, (ERROR_NO_MEMORY, "list(): error in malloc()"));
  }

  dn_new->vno = vno;
  dn_new->val = surf->vertices[vno].val;
  dn_new->next = NULL;

  /* - case 1: no elements in the list yet - */

  if (d_list == NULL) {
    d_list = dn_new;
    return(NO_ERROR);
  }

  /* - case 2: the dn_new element goes at the front of the list - */

  if (dn_new->val < d_list->val) {
    dn_new->next = d_list;
    d_list = dn_new;
    return(NO_ERROR);
  }

  /* - place "dn_new" in the list - */

  /* "dn_this" will always precede "dn_new" */

  dn_this = d_list;

  while (dn_this->next != NULL) {
    if (dn_this->next->val > dn_new->val) {
      dn_new->next = dn_this->next;
      dn_this->next = dn_new;
      return(NO_ERROR);
    }
    dn_this = dn_this->next;
  }

  dn_this->next = dn_new;

  return(NO_ERROR);

} /* end list() */

void unlist(int vno) {

  struct d_node *last, *dn;

  last = NULL;
  for (dn = d_list;dn != NULL;dn = dn->next) {
    if (dn->vno == vno) {
      if (last == NULL)
        d_list = dn->next;
      else
        last->next = dn->next;
      free(dn);
      break;
    }
    last = dn;
  }

  return;

} /* end unlist() */

int
mark(MRIS *surf, int vno, int m, bool b_overwrite) {
  // HISTORY
  // 18 November 2004
  // o Added b_overwrite to force marking irrespective
  //  of current value. Needed in the case of a "HUP"
  //   signal.
  //

  VERTEX *v;
  int rv;

  v = &surf->vertices[vno];

  if (!b_overwrite) {
    if (v->marked == DIJK_IN_PLAY) {
      unlist(vno);
    }
    if (m == DIJK_IN_PLAY) {
      rv = addToList(surf, vno);
      if (rv != NO_ERROR)
        return(rv);
    }
  }
  v->marked = m;
  return(NO_ERROR);
} /* end mark() */

int dijkstra(
    s_env&          st_env,
    float           af_maxAllowedCost,
    bool            ab_surfaceCostVoid)
{
    int             i, j;
    int             vno_c   = -1;
    int             vno_n   = -1;
    int             vno_i, vno_f;
    float           cost, f_pathCost;
    struct d_node   *dn, *dn_next;
    int             rv;
//    s_iterInfo      st_iterInfo;
    MRIS*           surf                = st_env.pMS_active;
//    bool            b_relNextReference  = true;

  // If we aren't going to preserve cost history in the environment, then we
  // will by default always be able to write path costs
  bool              b_canWriteCostVal   = !st_env.b_costHistoryPreserve;
  static int        calls               = 0;
  int               marked              = 0;
  int               totalLoops          = -1;

  /* --- sanity checks --- */
  vno_i  = st_env.startVertex;
  vno_f  = st_env.endVertex;

  assert(vno_i >= 0);
  assert(vno_i < surf->nvertices);
  assert(vno_f >= 0);
  assert(vno_f < surf->nvertices);

  if (!st_env.b_costHistoryPreserve) {
    assert(!surf->vertices[vno_i].ripflag);
    assert(!surf->vertices[vno_f].ripflag);
  }

  /* --- initialize --- */
  for (i = 0; i < surf->nvertices; i++) {
    bool b_overwrite = true;
    if (mark(surf, i, DIJK_VIRGIN, b_overwrite) != NO_ERROR)
      goto error;
    // Set all vertex values to -1 - only the very first time
    // that this function is called, or if explicitly
    // specified in the calling parameters.
    if (!calls || ab_surfaceCostVoid) surf->vertices[i].val = -1;
  }
  calls++;

  surf->vertices[vno_i].val = 0.0;
  surf->vertices[vno_i].old_undefval = vno_f;
  if (mark(surf, vno_i, DIJK_IN_PLAY) != NO_ERROR)
    goto error;

  // If the start and end vertices are coincident in the problem environment,
  // we should loop through at least once, and ignore the vno_c!=vno_f
  // condition, otherwise the while() will terminate after one loop.
  while ( vno_c!=vno_f || !totalLoops) {
    totalLoops++;
    if(totalLoops >= st_env.pMS_primary->nvertices-1) {
      // If this condition is true, we have processed all available vertices
      // in the mesh -- typically only occurs if the startVertex == endVertex
      // and is used for 'autodijk' type calculations.
      rv = TRUE;
      goto clean;
    }
    
    /* set vno_c (find min) */
    if (d_list == NULL) {
//      ErrorPrintf(ERROR_BADPARM, "dijkstra(): out of vertices");
//      colprintf(st_env.lw, st_env.rw, "start:stop", "[ %d:%d ]\n",
//                st_env.startVertex,
//                st_env.endVertex);
      goto error;
    }

    vno_c  = d_list->vno;
    
    
    VERTEX_TOPOLOGY const * const v_ct = &surf->vertices_topology[vno_c];
    VERTEX          const * const v_c  = &surf->vertices         [vno_c];

    /* mark it */
    if (mark(surf, vno_c, DIJK_DONE) != NO_ERROR)
      goto error;

    /* update neighbors */
    //cout << "neighbors = " << (int) v_c->num << endl;
    for (j = 0; j < (int) v_ct->vnum; j++) {
      //cout << "neighbor = " << j << endl;
      vno_n = v_ct->v[j];
      VERTEX * const v_n  = &surf->vertices[vno_n];

      //if(v_n->ripflag) continue;

      if (v_n->marked == DIJK_DONE) continue; // for "circular" path searches

//      cost = st_env.costFunc_do(st_env, &st_iterInfo, vno_c, j, 
//          			b_relNextReference);
      cost = s_env_edgeCostFind(st_env, vno_c, vno_n);
      f_pathCost = v_c->val + cost;

      // Break out of while if af_maxAllowedCost is violated.
      if (af_maxAllowedCost && (f_pathCost > af_maxAllowedCost)) continue;

      if ( (v_n->marked == DIJK_VIRGIN) || (f_pathCost < v_n->val) ) {
        // The check in the <if> statement preserves pathCost values in the MRIS
        // from previous calls to this function. This history is important
        // in determing ply distances from a given target path. A pathCost
        // is only written to a new vertex iff that vertex has never been
        // visited, or if the cost is less than an older value.
        if (st_env.b_costHistoryPreserve) {
          b_canWriteCostVal = (f_pathCost<v_n->val||v_n->val==-1);
        }
        if (b_canWriteCostVal) {
          marked++;
          v_n->val   = f_pathCost;
          v_n->old_undefval  = vno_c;
          //cout << vno_c << "<---" << vno_n << endl;
        }
      }
//     cout << "v->marked in dijkstra " << v_n->marked << endl;
      if (v_n->marked == DIJK_VIRGIN)
        if (mark(surf, vno_n, DIJK_IN_PLAY) != NO_ERROR) goto error;
    }
  }

  rv = TRUE;
  goto clean;

error:
  rv = FALSE;

clean:
  for (dn = d_list;dn != NULL;dn = dn_next) {
    dn_next = dn->next;
    free(dn);
  }
  d_list = NULL;

  return(rv);

} /* end dijkstra() */

/* eof */
