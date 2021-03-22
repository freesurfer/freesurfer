/**
 * @brief The surface related object API.
 *
 * Surface type functions setting whole surface function pointers, and
 * function definitions.
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


#include "c_surface.h"
#include "c_vertex.h"
#include "c_label.h"

#include <sstream>

void
surface_vertexFunction_do(
      s_env&            st_env,
      bool              (*vertex_satisfyTestCondition)
                              (VERTEX*  pvertex,
                               void*    pv_extra),
      void*             apv_conditional,
      void              (*vertex_function)
                              (VERTEX*  pvertex,
                               void*    pv_extra),
      void*             apv_functional
) {

    int         i;
    VERTEX*     pvertex;

  for (i = 0;i < st_env.pMS_active->nvertices;i++) {
    pvertex = &st_env.pMS_active->vertices[i];
    if (vertex_satisfyTestCondition(pvertex, apv_conditional)) {
      vertex_function(pvertex, apv_functional);
    }
  }
}

void
surface_annotation_do(
  s_env&   st_env
) {
  int  i;
  VERTEX* pvertex;

  for (i = 0;i < st_env.pMS_active->nvertices;i++) {
    pvertex = &st_env.pMS_active->vertices[i];
    pvertex->annotation = i;
  }
}

void
surface_vertexPatternCopy(
    s_env&      st_env,
    MRIS*       apMS_source,
    MRIS*       apMS_target,
    bool        (*vertex_satisfyTestCondition)
    (VERTEX* pvertex,
    void*       pv_extra),
    void*       apvsource_extra,
    void        (*vertex_modify)
    (VERTEX*    pvertex,
    void*       pv_extra),
    void*       apvtarget_extra
) {

  //
  // PRECONDITIONS
  // o The source surface has a information spread across its vertices that
  //   is defined by (*vertex_pattern).
  //  o The source and target surfaces are not checked explicitly for
  //   compatibility!
  //
  // POSTCONDITIONS
  // o The target surface will have the same vertex pattern as the source.
  // o Any other information in the target surface that is not specified
  //   by (*vertex_pattern) is left untouched.
  //

  int  i;
  VERTEX* pvertex_source;
  VERTEX* pvertex_target;

  for (i = 0;i < apMS_source->nvertices;i++) {
    pvertex_source = &apMS_source->vertices[i];
    if (vertex_satisfyTestCondition(pvertex_source, apvsource_extra)) {
      pvertex_target = &apMS_target->vertices[i];
      vertex_modify(pvertex_target, apvtarget_extra);
    }
  }
}

void
surface_rawCurveMinMax_do(
  s_env&  st_env
) {

  s_minMax  sMM;
  void*  pv_sMM = (void*) &sMM;
  void*  pv_dummy = NULL;

  stringstream sout("");

  sMM.f_min  = 1;
  sMM.f_max  = -1;

  surface_vertexFunction_do( st_env,
                             vertex_alwaysTrue,
                             pv_dummy,
                             vertex_rawCurveMinMax,
                             pv_sMM
                           );

  sout << "Min:\t" << sMM.f_min << " (at vertex) " << sMM.minVertex << endl;
  sout << "Max:\t" << sMM.f_max << " (at vertex) " << sMM.maxVertex << endl;

  nRLOUT(sout.str());
}

void
surface_rawCurveSum_do(
  s_env&  st_env
) {

  s_rawCurve  sRC;
  void*  pv_sRC = (void*) &sRC;
  void*  pv_dummy = NULL;

  stringstream sout("");

  sRC.f_posCount = 0;
  sRC.f_negCount = 0;

  surface_vertexFunction_do( st_env,
                             vertex_alwaysTrue,
                             pv_dummy,
                             vertex_rawCurveSum,
                             pv_sRC
                           );

  sout << "Positive sum: " << sRC.f_posCount << endl;
  sout << "Negative sum: " << sRC.f_negCount << endl;

  nRLOUT(sout.str());
}

void
surface_signumFunction_do(
  s_env&  st_env
) {

  s_signumFunctional sgnf;
  void*  pv_sgnf = (void*) &sgnf;
  void*  pv_dummy = NULL;

  stringstream sout("");

  sgnf.posCount = 0;
  sgnf.negCount = 0;
  sgnf.zeroCount = 0;

  surface_vertexFunction_do( st_env,
                             vertex_alwaysTrue,
                             pv_dummy,
                             vertex_signumFunctional,
                             pv_sgnf
                           );

  sout << "Positive count: " << sgnf.posCount << endl;
  sout << "Negative count: " << sgnf.negCount << endl;
  sout << "Zero     count: " << sgnf.zeroCount << endl;

  nRLOUT(sout.str());
}


void
surface_averageIntegratedCurveArea_do(
  s_env&  st_env,
  e_CURVATURE ae_curvature
) {
  //
  // PRECONDITIONS
  //  o Make sure that pMS_active is set to the correct surface!
  //
  // POSTCONDITIONS
  //  o Rectified integrated curvature measurement is transmitted to the
  //   result log channel.
  //  o NB!! Will fail if working/aux surfaces are re-defined while engine
  //   is active.
  //

  s_integratedCurve sIC;
  void*  pv_sIC = (void*) &sIC;
  void*  pv_dummy = NULL;
  static int   calls = 0;

  stringstream sout("");

  sIC.e_curvature = ae_curvature;
  sIC.f_sum  = 0.0;

  if (!calls) {
    ULOUT("Computing Second Fundamental Form on primary\t\t\t");
    MRISsetNeighborhoodSizeAndDist(st_env.pMS_primary, 2) ;
    MRIScomputeSecondFundamentalForm(st_env.pMS_primary);
    nULOUT("[ ok ]\n");
    ULOUT("Computing Second Fundamental Form on secondary\t\t\t");
    MRISsetNeighborhoodSizeAndDist(st_env.pMS_secondary, 2) ;
    MRIScomputeSecondFundamentalForm(st_env.pMS_secondary);
    nULOUT("[ ok ]\n");
  }

  surface_vertexFunction_do( st_env,
                             vertex_alwaysTrue,
                             pv_dummy,
                             vertex_curveAreaSum,
                             pv_sIC
                           );

  float f_avIC = 0.0;
  f_avIC = sIC.f_sum / st_env.pMS_active->nvertices;

  string str_curvature;
  switch (ae_curvature) {
  case e_gaussian:
    str_curvature  = "Gaussian";
    break;
  case e_mean:
    str_curvature = "mean";
    break;
  }
  sout << "Average integrated rectified " << str_curvature << " curvature:";
  sout << endl << f_avIC << endl;

  nRLOUT(sout.str());
  calls++;
}

void
surface_correlationFunction_do(
  s_env&  st_env
) {
  //
  //PRECONDITIONS
  // o The "target" surface (usually the working curvature) should
  //   have a rip pattern denoting a path of interest. Typically this
  //   implies that after a standard path search, the output label needs
  //   to have been re-loaded in order to restore the rip pattern.
  //
  // POSTCONDITIONS
  //  o A float correlation is dumped on the results channel. If this
  //   correlation is '-1', then *no* rips were found on the working
  //   surface. This denotes an error state.
  //

  float  f_correl  = 0.0;
  int  denom  = 0;
  int  membership = 0;
  VERTEX* pvertex  = NULL;
  void* pv_void  = NULL;

  for (int i = 0;i < st_env.pMS_primary->nvertices;i++) {
    pvertex = &st_env.pMS_primary->vertices[i];
    if (vertex_ripFlagIsTrue(pvertex, pv_void)) {
      denom++;
      pvertex = &st_env.pMS_secondary->vertices[i];
      if (vertex_ripFlagIsTrue(pvertex, pv_void)) {
        //cout << i << endl;
        membership++;
      }
    }
  }

  //cout << "number overlap: " << membership << endl;
  //cout << "size of source: " << denom << endl;

  if (!denom)
    f_correl = -1;
  else
    f_correl = (float)membership / (float)denom;

  cout << "Correlation: " << f_correl << endl;

  stringstream sout("");
  sout << f_correl;
  nRLOUT(sout.str());
}

void
surface_primaryToSecondary_ripTrueCopy(
  s_env&  st_env
) {

  char      ch_rip;
  char*     pch_rip;
  void*     pv_rip;

  ch_rip    = (char) TRUE;
  pch_rip   = &ch_rip;
  pv_rip    = (void*) pch_rip;

  surface_vertexPatternCopy( st_env,
                             st_env.pMS_primary,
                             st_env.pMS_secondary,
                             vertex_ripFlagIsTrue,
                             pv_rip,
                             vertex_ripFlagMark,
                             pv_rip
                           );
}

float
surface_ripMark(
    s_env&              st_env
) {
    //
    // ARGS
    // st_env   in   environment data
    //
    // DESCRIPTION
    // Adds a TRUE to the ripflag of each vertex that is part of the
    // Dijskstra path. Also saves the route itself and cost information
    // to a text file.
    //
    // PRECONDITIONS
    // o Rips on the "active" surface are marked.
    // o The cumulative cost of traveling along the marked path is returned.
    //
    // HISTORY
    // 18 November 2004
    //  o Initial design and coding.
    //
    // 10 March 2005
    // o Added "active" surface.
    //
    // 30 October 2009
    // o Added cumulative cost value to result channel and stdout.
    //

//    s_iterInfo  st_iterInfo;
    string      str_costFile            = st_env.str_workingDir +
                                          st_env.str_costFileName;
    int         i;
    int         ii, jj;
    float       f_cost                  = 0.;
    float       f_costSum               = 0.;
//    bool        b_relNextReference    = false;

    ofstream    ofs(str_costFile.c_str(), ios::out);
    ofs.flags(ios::fixed );

    for (i = st_env.endVertex; i != st_env.startVertex;
        i = st_env.pMS_active->vertices[i].old_undefval) {
        ii = st_env.pMS_active->vertices[i].old_undefval;
        jj = i;
//        f_cost = st_env.costFunc_do(st_env, &st_iterInfo, ii, jj,
//                                    b_relNextReference);
	f_cost = s_env_edgeCostFind(st_env, ii, jj);
        st_env.pMS_active->vertices[i].ripflag = TRUE;
        f_costSum += f_cost;
        if(st_env.b_costPathSave) {
            ofs << ii << "\t" << jj;
            ofs << "\t"  << f_cost;
            ofs << "\t"  << st_env.pst_iterInfo->iter;
            ofs << "\t"  << st_env.pst_iterInfo->f_distance;
            ofs << "\t"  << st_env.pst_iterInfo->f_curvature;
            ofs << "\t"  << st_env.pst_iterInfo->f_sulcalHeight;
            ofs << "\t"  << st_env.pst_iterInfo->f_dir;
            ofs << endl;
        }
    }
    st_env.pMS_active->vertices[i].ripflag = TRUE;
    if(st_env.b_costPathSave) ofs.close();
    return f_costSum;
}

void
surface_ripClear(
  s_env&   st_env,
  bool   b_wholeSurfaceForce
) {
  //
  // ARGS
  // st_env                     in              environment data
  // b_wholeSurfaceForce        in              if true, force a "clear"
  //                                            across the whole
  //                                            surface
  //
  // DESCRIPTION
  //  Clears all ripflags. Necessary for "HUP" considerations.
  //
  // PRECONDITIONS
  // o Assumes that patch_ripMark() has been called and a rip path
  //   marked.
  //
  // POSTCONDITIONS
  // o If b_wholeSurfaceForce is true, then every vertex node
  //   will have its 'ripflag' set to 'FALSE'.
  //
  // HISTORY
  // 18 November 2004
  //  o Initial design and coding.
  //
  // 17 February 2005
  //  o Added b_wholeSurfaceForce
  //
  // 10 March 2005
  // o Added active surface
  //
  // 7 September 2006
  //  o Test / retest experiments sometimes fail in certain cases. It seems
  //   as though interference from a previous run might be the cause. To
  //   "fix" this, the default behaviour of this function has been changed
  //   to also reset the cost value of a vertex.
  //
  // 22 January 2010
  // o Fixed logic error with wholeSurface iteration! Did not have
  //   {} around for loop internals!
  //

  int   i;

  if (!b_wholeSurfaceForce) {
    for (i = st_env.endVertex; i != st_env.startVertex;
         i = st_env.pMS_active->vertices[i].old_undefval) {
      st_env.pMS_active->vertices[i].ripflag  = FALSE;
      st_env.pMS_active->vertices[i].val  = -1.0;
    }
    st_env.pMS_active->vertices[i].ripflag   = FALSE;
    st_env.pMS_active->vertices[i].val  = -1.0;
  } else {
    for (i = 0; i < st_env.pMS_active->nvertices; i++) {
      st_env.pMS_active->vertices[i].ripflag  = FALSE;
      st_env.pMS_active->vertices[i].val  = -1.0;
    }
  }
}

/* eof */
