/**
 * @brief The label related object API.
 *
 * Label type functions include saving / loading surface structures
 * using label-type files.
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


#include "c_label.h"
#include "dijkstra.h"

#include "c_surface.h"
#include "c_vertex.h"

#include "C_mpmProg.h"
#include "asynch.h"

#include <sstream>

void
label_ply_do(
    s_env       &ast_env
) {
    //
    // PRECONDITIONS
    // o The plyDepth field of ast_env is set.
    //  o The ply function has been assigned.
    // o A pattern of TRUE ripflags across the vertex space corresponding
    //   to the label to ply has been defined.
    //
    // POSTCONDITIONS
    // o Using the dijkstra engine, each vertex of the label is surrounded by
    //   concentric "circles". The overlap of these circles defines a unidistant
    //   delta region about the path defined by the labelled vertices.
    //

    // Should any cost values remaining in the surface be zeroed? Yes.

    bool        b_origHistoryFlag       = ast_env.b_costHistoryPreserve;
    bool        b_surfaceCostVoid       = false;
    int         i                       = 0;
    int         j                       = 0;
    float       f_plyDepth              = s_env_plyDepth_get(ast_env);
    ast_env.b_costHistoryPreserve       = true;

    //s_env_costFctSet(&ast_env, costFunc_unityReturn, e_unity);

    for (i=0; i<ast_env.pMS_active->nvertices; i++) {
        if (ast_env.pMS_active->vertices[i].ripflag == TRUE) {
        b_surfaceCostVoid  = !j++;
        ast_env.startVertex = i;
        ast_env.endVertex = i;
        dijkstra(ast_env, f_plyDepth, b_surfaceCostVoid);
        }
    }
    ast_env.b_costHistoryPreserve = b_origHistoryFlag;
}

void
label_singleVertexSet(
    MRIS*       apmris,
    int         avertex,
    void        (*vertex_labelMark)
    (VERTEX*    pvertex,
        void*   marker),
    void*       apv_marker
) {
    //
    // POSTCONDITIONS
    // o The existing pattern of TRUE ripflags is cleared across the entire
    //  surface.
    // o The single vertex index has is passed to the vertex function.
    // o Any Dijkstra "weights" that might already exist across the surface
    //  are left intact.
    //
    // HISTORY
    // 09 March 2005
    // o Initial design and coding.
    //

    VERTEX*  pvertex;
    int   i;

    for (i = 0; i < apmris->nvertices; i++)
        apmris->vertices[i].ripflag = FALSE;

    pvertex = &apmris->vertices[avertex];
    vertex_labelMark(pvertex, apv_marker);
}

bool
label_terminalsFind(
  MRIS*   apmris,
  string   astr_fileName,
  deque<int>&  aque_terminal

) {
  //
  // PRECONDITIONS
  // o No checks are performed on <astr_fileName>.
  //
  // POSTCONDITIONS
  // o This function does not modify *any* internal data structures -
  //  it is for most practical purposes, standalone.
  // o <aque_terminal> contains a list of "terminal" vertex numbers.
  //
  // HISTORY
  // 22 July 2005
  // o Initial design and coding.
  //

  LABEL*  pLBL;
  int   vno_i   = 0;
  int   vno_j  = 0;
  int   vno_k  = 0;
  int   i, j, k;
  short  inLabel  = 0;
  bool  b_ret  = false;

  pLBL = LabelRead((char*)"", (char*) astr_fileName.c_str());
  if (pLBL == NULL)
    error_exit("allocating a pLBL structure", "some error occurred", 1);

  aque_terminal.clear();
  for (i=0; i <pLBL->n_points; i++) {
    vno_i  = pLBL->lv[i].vno;
    VERTEX_TOPOLOGY const * const pvertex = &apmris->vertices_topology[vno_i];
    inLabel = 0;
    for (j=0; j<pvertex->vnum; j++) {
      vno_j = pvertex->v[j];
      for (k=0; k< pLBL->n_points; k++) {
        vno_k = pLBL->lv[k].vno;
        if (vno_j==vno_k)
          inLabel++;
      }
    }
    if (inLabel==1) {
      aque_terminal.push_back(vno_i);
      b_ret = true;
    }
  }
  LabelFree(&pLBL);
  return b_ret;
}


void
label_coreLoad(
    MRIS*       apmris,
    string      astr_fileName,
    void        (*vertex_labelMark)
    (VERTEX*    pvertex,
        void*   marker),
    void*       apv_marker
) {
  //
  // POSTCONDITIONS
  // o The existing pattern of TRUE ripflags is cleared across the entire
  //  surface.
  // o A label is read from disk, and vertices on the passed surface are
  //  marked according to (*vertex_labelMark).
  // o No ordered "connection" profile is created between nodes.
  // o Any Dijkstra "weights" that might already exist across the surface
  //  are left intact.
  //
  // HISTORY
  // 16 February 2005
  // o Initial design and coding.
  //
  // 15 February 2006
  // o Convert <astr_fileName> to absolute directory spec.
  //

  LABEL*    pLBL;
  VERTEX*   pvertex;
  int       vno         = 0;
  int       i;

  pLBL = LabelRead((char*)"", (char*) astr_fileName.c_str());
  if (pLBL == NULL)
    error_exit("allocating a pLBL structure", "some error occurred", 1);

  for (i = 0; i < apmris->nvertices; i++)
    apmris->vertices[i].ripflag = FALSE;

  for (i = 0; i < pLBL->n_points; i++) {
    vno  = pLBL->lv[i].vno;
    pvertex = &apmris->vertices[vno];
    vertex_labelMark(pvertex, apv_marker);
  }
  LabelFree(&pLBL);
}

void
label_workingSurface_loadFrom(
    s_env&      st_env,
    void        (*vertex_labelMark)
    (VERTEX*    pvertex,
        void*   marker),
    void*       apv_marker
) {
    //
    // ARGS
    // st_env   in   environment data
    //
    // DESCRIPTION
    //  This is thin wrapper that calls label_coreLoad(...) correctly.
    // It accesses the working surface label file as defined in the
    // environment and loads this into the internal working surface
    // structure. Vertices are marked with the vertex_labelMark function.
    //
    // HISTORY
    // 07 February 2005
    // o Split from "label_save()".
    //

    string      str_labelFileName;
    MRIS*       pMS_primary;

    bool        b_clearWholeSurface = true;
    surface_ripClear(st_env, b_clearWholeSurface);

    str_labelFileName   = st_env.str_workingDir + st_env.str_labelFileName;
    pMS_primary         = st_env.pMS_primary;
    ULOUT(str_labelFileName);
    label_coreLoad(pMS_primary, str_labelFileName, vertex_labelMark, apv_marker);
    nULOUT("\t\t\t\t[ ok ]");
}

void
label_coreSave(
        MRIS*                   apmris,
        string                  astr_fileName,
        bool   (*vertex_satisfyTestCondition)
                    (VERTEX* apvertex,
        void*  apv_void),
        void*   apv_fromCaller
) {
      //
      // ARGS
      // amris              in              surface to examine for label
      //                                    "rips"
      // astr_fileName      in              filename to contain the
      //                                    saved label
      //
      // DESCRIPTION
      // Saves a label file onto the passed surface
      //
      // POSTCONDITIONS
      // o Label file is written to disk.
      //
      // HISTORY
      // 07 February 2005
      // o Split from "label_save()".
      //

      LABEL*        pLBL;
      VERTEX*       pvertex;
      int           n = 0;
      int           i;

      for (i = 0;i < apmris->nvertices;i++) {
    // if(apmris->vertices[i].ripflag == TRUE)
        pvertex = &apmris->vertices[i];
        if (vertex_satisfyTestCondition(pvertex, apv_fromCaller))
          n++;
      }
      pLBL = LabelAlloc(n, (char*)"", (char*) (astr_fileName).c_str());
      if (pLBL == NULL)
        error_exit("allocating a pLBL structure", "some error occurred", 1);

      for (i = 0;i < apmris->nvertices;i++) {
        pvertex = &apmris->vertices[i];
        if (vertex_satisfyTestCondition(pvertex, apv_fromCaller)) {
          pLBL->lv[pLBL->n_points].vno = i;
          pLBL->lv[pLBL->n_points].x = apmris->vertices[i].x;
          pLBL->lv[pLBL->n_points].y = apmris->vertices[i].y;
          pLBL->lv[pLBL->n_points].z = apmris->vertices[i].z;
          pLBL->n_points++;
        }
      }
      LabelWrite(pLBL, (char*) (astr_fileName).c_str());
      LabelFree(&pLBL);
}

void
label_ply_save(
        s_env&        st_env,
        string        astr_filePrefix,
        bool          b_staggered
) {
  //
  // PRECONDITIONS
  //  o ply labels should be defined on the internal working surface.
  //
  // POSTCONDITIONS
  // o A series of label files are saved to the current working directory.
  //   These are labelled <astr_filePrefix>.<N>.label where
  //   N = [1... <plyDepth>] and denote a label containing vertices that
  //   are less than <N> cost function evaluations distant from a
  //   reference label curve.
  //

  string                str_labelFileName;
  MRIS*                 pMS_surface;
  stringstream          sout("");

  float                 f_plyStart;
  float                 f_plyIncrement;
  float                 f_plyDepth;
  float*                pf_plyDepth;
  void*                 pv_fromCaller;

  pMS_surface  = st_env.pMS_active;

  f_plyIncrement        = s_env_plyIncrement_get(st_env);
  f_plyDepth            = s_env_plyDepth_get(st_env);
  f_plyStart            = f_plyDepth;
  if (b_staggered)
    f_plyStart = 0.0;

  while(f_plyStart <= f_plyDepth) {
    pf_plyDepth = &f_plyStart;
    pv_fromCaller = (void*) pf_plyDepth;
    sout << st_env.str_workingDir + astr_filePrefix;
    if (b_staggered) sout << "-ply" << f_plyStart;
    sout << ".label";
    ULOUT(sout.str());
    label_coreSave(pMS_surface, sout.str(),
                   vertex_valLTE, pv_fromCaller);
    nULOUT("\t\t\t\t\t\t\t[ ok ]\n");
    sout.str("");
    f_plyStart += f_plyIncrement;
  }
}

void
label_secondarySurface_saveTo(
    s_env&      st_env,
    bool        (*vertex_satisfyTestCondition)
    (VERTEX*    pvertex,
        void*   apv_void),
    void*       apv_fromCaller
) {
    //
    // ARGS
    // st_env   in   environment data
    //
    // DESCRIPTION
    // This is thin wrapper that calls label_coreSave(...) correctly.
    //  It access marked vertices on the original surface and saves these to
    // a label file as specified in the environment structure.
    //
    // HISTORY
    // 07 February 2005
    // o Split from "label_save()".
    //

    string        str_labelFileName;
    MRIS*         pMS_secondarySurface;

    str_labelFileName   = st_env.str_workingDir + st_env.str_labelFileNameOS;
    pMS_secondarySurface= st_env.pMS_secondary;
    label_coreSave(pMS_secondarySurface, str_labelFileName,
                 vertex_satisfyTestCondition, apv_fromCaller);
}

void
label_workingSurface_saveTo(
    s_env&      st_env,
    bool        (*vertex_satisfyTestCondition)
    (VERTEX*    pvertex,
        void*   apv_void),
    void*       apv_fromCaller
) {
    //
    // ARGS
    // st_env   in   environment data
    //
    // DESCRIPTION
    //  This is thin wrapper that calls label_coreSave(...) correctly.
    //  It access marked vertices on the working surface and saves these to
    // a label file as specified in the environment structure.
    //
    // HISTORY
    // 07 February 2005
    // o Split from "label_save()".
    //

    string      str_labelFileName;
    MRIS*       pMS_primary;

    str_labelFileName   = st_env.str_workingDir + st_env.str_labelFileName;
    pMS_primary         = st_env.pMS_primary;
    label_coreSave(pMS_primary, str_labelFileName,
                 vertex_satisfyTestCondition, apv_fromCaller);
}


/* eof */
