/***************************************************************************
 *   Copyright (C) 2004 by Rudolph Pienaar / Christian Haselgrove          *
 *    Center for Morphometric Analysis       *
 *    Massachusetts General Hospital        *
 * Building 149, 13th St.         *
 *  Charlestown, MA 02129         *
 *     {ch|rudolph}@nmr.mgh.harvard.edu      *
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
/// \file c_label.h
///
/// \brief Brief description
/// The label related object API.
///
/// \b DESCRIPTION
/// Label type functions include saving / loading surface structures
/// using label-type files.
///
/// \b HISTORY
/// 15 March 2005 - Initial consolidation from several other sources.
/// $Id: c_label.h,v 1.1 2009/09/08 22:39:27 nicks Exp $
///
///

#ifndef __C_LABEL_H__
#define __C_LABEL_H__

#include "general.h"
#include "env.h"
#include <deque>

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

#include <string>
using namespace std;


/// \fn void label_terminalsFind(s_env &ast_env)
/// \brief Find the "terminals" associated with a label. These terminals are
///  the "end points" of a "line" label. More formally, terminals
///  are vertices that link to a single neighbour.
/// \param apmris    surface to process
/// \param astr_fileName  filename containing the label to load
/// \param deque<int>  STL queue containing the terminal
///     vertex numbers
/// \return bool - TRUE if terminals were found, FALSE otherwise
bool
label_terminalsFind(
  MRIS*   apmris,
  string   astr_fileName,
  deque<int>&  aque_terminal

);

/// \fn void label_ply_do(s_env &ast_env)
/// \brief Perform a ply "spread" about any vertices that have their ripflag
///  fields set to TRUE
/// \param ast_env   the program environment
/// \return    (void)
void
label_ply_do(
  s_env   &ast_env
);

/// \fn void singleVertexSet(MRIS* apmris, int avertex, void (*vertex_labelMark) (VERTEX* pvertex, void* marker), void* apv_marker )
/// \fn void singleVertexSet(MRIS* apmris, int avertex, void (*vertex_labelMark) (VERTEX* pvertex, void* marker), void* apv_marker )
/*/// \fn label_coreLoad */
/// \brief For a passed surface, set the specified vertex index using the callback function, vertex_labelMark().
/// \param apmris    surface to process
/// \param avertex   index on the surface to a target vertex
/// \param vertex_labelMark a callback function that operates on
///     each loaded vertex
/// \param apv_marker  a pointer to void that embodies the
///     argument list to (*vertex_labelMark)()
/// \return    (void)
void
label_singleVertexSet(
  MRIS*   apmris,
  int   avertex,
  void   (*vertex_labelMark)
  (VERTEX* pvertex,
   void*  marker),
  void*   apv_marker
);

/// \fn void label_coreLoad(MRIS& amris, string astr_fileName, bool (*vertex_labelMark)(VERTEX* pvertex,void* marker), void* apv_marker )
/// \fn void label_coreLoad(MRIS& amris, string astr_fileName, bool (*vertex_labelMark)(VERTEX* pvertex,void* marker), void* apv_marker )
/*/// \fn label_coreLoad */
/// \brief The core label function. Typically called by wrapper functions within a specific context. A label file is read from disk, and each vertex of this label is marked on the passed surface by vertex_labelMark().
/// \param amris    surface on to which the label is loaded
/// \param astr_fileName  filename containing the label to load
/// \param vertex_labelMark a callback function that operates on
///     each loaded vertex
/// \param apv_marker  a pointer to void that embodies the
///     argument list to (*vertex_labelMark)()
/// \return    (void)
void
label_coreLoad(
  MRIS*   apmris,
  string   astr_fileName,
  void   (*vertex_labelMark)
  (VERTEX* pvertex,
   void*  marker),
  void*   apv_marker
);

/// \fn void workingSurface_loadFrom(s_env& st_env, void (*vertex_labelMark)(VERTEX* pvertex, void* marker), void* apv_marker);
/// \fn void workingSurface_loadFrom(s_env& st_env, void (*vertex_labelMark)(VERTEX* pvertex, void* marker), void* apv_marker);
/*/// \fn label_coreSave*/
/// \brief Load a surface from a label filename as specified in the program
///  environment. Each label file is processed by
///  (*vertex_labelMark)().
/// \param ast_env   program core environment
/// \param vertex_labelMark a callback function that operates on
///     each loaded vertex
/// \param apv_marker  a pointer to void that embodies the
///     argument list to (*vertex_labelMark)()
/// \return    (void)
void
label_workingSurface_loadFrom(
  s_env&   st_env,
  void   (*vertex_labelMark)
  (VERTEX* pvertex,
   void*  marker),
  void*   apv_marker
);

/// \fn void label_coreSave(MRIS& amris, string astr_fileName, bool *vertex_satisfyTestCondition)
/// \fn void label_coreSave(MRIS& amris, string astr_fileName, bool *vertex_satisfyTestCondition)
/*/// \fn label_coreSave*/
/// \brief The core label save function. Typically called by wrapper functions within a specific context. This parses a processed surface, and writes a label file to disk. This label file contains each vertex flagged by the callback function vertex_satisfyTestCondition().
/// \param amris     surface from which to save the
///       label file
/// \param astr_fileName   filename to contain the saved
///      label
/// \param  *vertex_satisfyTestCondition  a call back function used on
///      each vertex to determine
///       whether it should be saved or
///      not.
/// \return    (void)
void
label_coreSave(
  MRIS*   apmris,
  string   astr_fileName,
  bool   (*vertex_satisfyTestCondition)
  (VERTEX* apvertex,
   void*  apv_void),
  void*   apv_fromCaller
);

/// \fn void label_ply_save(s_env &st_env, string astr_filePrefix, bool b_staggered = true)
/// \brief Save concentric ply areas.
/// \param st_env   the program environment
/// \param astr_fileName  the filename prefix to save ply labels
///     If b_staggered is true, each integer
///     increase in a ply distance will be
///     saved to a different file.
/// \param b_staggered  a boolean flag that toggles "staggered"
///     saves ON | OFF.
/// \return    (void)
void
label_ply_save(
  s_env&  st_env,
  string  astr_filePrefix,
  bool  b_staggered  = true
);

/// \fn void label_auxSurface_saveTo(s_env& st_env, bool (*vertex_satisfyTestCondition) (VERTEX* pvertex, void* apv_void), void* apv_fromCaller)
/// \fn void label_auxSurface_saveTo(s_env& st_env, bool (*vertex_satisfyTestCondition) (VERTEX* pvertex, void* apv_void), void* apv_fromCaller)
/*/// \fn label_coreLoad */
/// \brief Save the auxillary surface to the label filename specified in the
///  program environment. Each label file is processed by
///  (*vertex_satisfyTestCondition)() before being added to the
///  save label structure.
/// \param st_env    the program environment
/// \param vertex_satisfyTestCondition a callback function that
///      operates on each vertex and
///      test if this vertex is to be
///      saved.
/// \param apv_fromCaller   a pointer to void that embodies
///       the argument list to
///      (*vertex_satisfyTestCondition)()
/// \return    (void)
void
label_auxSurface_saveTo(
  s_env&   st_env,
  bool   (*vertex_satisfyTestCondition)
  (VERTEX* pvertex,
   void*  apv_void),
  void*   apv_fromCaller
);

/// \fn void label_workingSurface_saveTo(s_env& st_env, bool (*vertex_satisfyTestCondition) (VERTEX* pvertex, void* apv_void), void* apv_fromCaller)
/// \fn void label_workingSurface_saveTo(s_env& st_env, bool (*vertex_satisfyTestCondition) (VERTEX* pvertex, void* apv_void), void* apv_fromCaller)
/*/// \fn label_coreLoad */
/// \brief Save the working surface to the label filename specified in the
///  program environment. Each label file is processed by
///  (*vertex_satisfyTestCondition)() before being added to the
///  save label structure.
/// \param st_env    the program environment
/// \param vertex_satisfyTestCondition a callback function that
///      operates on each vertex and
///      test if this vertex is to be
///      saved.
/// \param apv_fromCaller   a pointer to void that embodies
///       the argument list to
///      (*vertex_satisfyTestCondition)()
/// \return    (void)
void
label_workingSurface_saveTo(
  s_env&   st_env,
  bool   (*vertex_satisfyTestCondition)
  (VERTEX* pvertex,
   void*  apv_void),
  void*   apv_fromCaller
);


#endif //__C_LABEL_H__


