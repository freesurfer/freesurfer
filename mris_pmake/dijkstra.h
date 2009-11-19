/***************************************************************************
 *   Copyright (C) 2004 by Rudolph Pienaar / Christian Haselgrove          *
 *    Center for Morphometric Analysis       *
 *    Massachusetts General Hospital        *
 *    Building 149, 13th St.         *
 *    Charlestown, MA 02129         *
 *    {ch|rudolph}@nmr.mgh.harvard.edu      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/// \file dijkstra.h
///
/// \brief
/// API for dijkstra related processing.
///
///
/// \b DESCRIPTION
///
/// Provides an API for dijkstra search through freesurfer structures.
///
/// \b HISTORY
///
///  Week of 20 September 2004 - kdevelop integration / cvs setup
/// 09 November 2004
/// o Added s_iterInfo
///

using namespace std;

#ifndef __DIJKSTRA_H__
#define __DIJKSTRA_H__

#include <string>
#include <iostream>
#include <fstream>

#ifdef __cplusplus
extern  "C" {
#endif

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "error.h"

#ifdef __cplusplus
}
#endif

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

