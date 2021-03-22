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

#ifndef __DIJKSTRA_H__
#define __DIJKSTRA_H__

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "error.h"

#include "general.h"
#include "env.h"
#include "pstream.h"

#define DIJK_VIRGIN  0
#define DIJK_IN_PLAY  1
#define DIJK_DONE  2

typedef enum {
  eDIJK_VIRGIN = 0,
  eDIJK_IN_PLAY = 1,
  eDIJK_DONE = 2
}
e_DIJKMARK;


/// d_node is the core structure that maintains a linked list of ordered
/// vertex nodes (with cost).
struct d_node {
  int    vno;
  float    val;
  struct d_node*  next;
};

/// \fn int list(MRIS *surf, int vno)
/// \brief Add another node to the linked list of minimum paths. This node is added in the correct position to maintain a sorted cost order.
/// \param  surf  The surface containing all the vertices
/// \param  vno The vertex index to add to the list
/// \return  (int) Error condition
int   addToList(  MRIS   *surf,
                  int   vno);

/// \fn void unlist(int vno)
/// \brief Remove a node from the list
/// \param  vno The vertex index to remove
/// \return  (none)
void   unlist(   int   vno);

/// \fn int mark(MRIS *surf,int vno, int m, bool b_overwrite = false)
/// \brief Write the value 'm' to the node vno in the surface surf.
/// \param  surf The surface structure
/// \param  vno The current vertex index
/// \param  m The value to write to the 'marked' field of the vertex
/// \param b_overwrite If true, ignore the current marked value and overwrite
/// \return    ERROR | NO_ERROR
int   mark(   MRIS   *surf,
              int   vno,
              int   m,
              bool  b_overwrite = false);

/// \fn int dijkstra(s_env& st_env, float *cost_fct)
/// \fn int dijkstra(s_env& st_env, float *cost_fct)
/*/// \fn int dijkstra */
/// \brief The core entry point into the dijkstra engine. The start and end vertex indices are contained in
///  the environment structure, st_env. This function determines the shortest path between these
/// vertices. If the optional af_maxAllowedCost is non-zero, if any path cost is found to exceed
/// this value, the function will return.
/// \param  st_env The problem environment
/// \param  *cost_fct A cost function to use on the weighted graph
/// \param af_maxAllowedCost An optional termination condition on the dijkstra function.
/// \return    ERROR | NO_ERROR
int   dijkstra(  s_env&  st_env,
                 float  af_maxAllowedCost = 0.0,
                 bool  ab_surfaceCostVoid = false);


typedef struct _node s_node;
typedef struct _node {
  int    id; // typically an MRIS vertex index
  float    f_val; // the value of this MRIS vertex
  s_node*   p_next; // pointer to next node
}
s_node;

class C_dlist {
  // data structures

protected:

  s_node*       mpnode_ordered;
  MRIS*         mpMS_surface;

  // methods

public:
  //
  // constructor / destructor block
  //
  C_dlist(
    MRIS*  apMS_surface
  );
  ~C_dlist();

  //
  // Add / remove / mark nodes
  //
  int node_addOrdered(
    int  avertex
  );
  int node_remove(
    int  avertex
  );
  int node_mark(
    int  avertex,
    e_DIJKMARK ae_dijkMark,
    bool  ab_forceMark = false
  );

};

#endif //__DIJKSTRA_H__

