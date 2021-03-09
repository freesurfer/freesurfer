/**
 * @brief The environment object API.
 *
 * This file contains old legacy code, consolidated here during testing
 * to facilitate its eventual removal from the mainline code.
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

#include "legacy.h"
#include "general.h"

extern string G_SELF;

void
s_weights_scan(
  s_weights&  st_costWeight,
  C_scanopt&  cso_options
) {
  //
  // ARGS
  // cso_options  in  scanopt structure to be parsed
  // st_costWeight  in/out  weight structure to be filled
  //
  // DESCRIPTION
  // Scans the options file structure for Dijkstra weights
  //
  // HISTORY
  // 17 November 2004
  // o Initial design and coding.
  //

  string str_value;

  if (cso_options.scanFor("wd", &str_value))
    st_costWeight.wd = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a wd weight.",
               100);
  if (cso_options.scanFor("wc", &str_value))
    st_costWeight.wc = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a wc weight.",
               101);
  if (cso_options.scanFor("wh", &str_value))
    st_costWeight.wh = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a wh weight.",
               102);
  if (cso_options.scanFor("wdc", &str_value))
    st_costWeight.wdc = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a wdc weight.",
               103);
  if (cso_options.scanFor("wdh", &str_value))
    st_costWeight.wdh = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a wdh weight.",
               104);
  if (cso_options.scanFor("wch", &str_value))
    st_costWeight.wch = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a wch weight.",
               105);
  if (cso_options.scanFor("wdch", &str_value))
    st_costWeight.wdch = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a wdch weight.",
               106);
  if (cso_options.scanFor("wdir", &str_value))
    st_costWeight.wdir = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a wdir weight.",
               107);
}

void
s_Dweights_scan(
  s_Dweights&  st_DcostWeight,
  C_scanopt&  cso_options
) {
  //
  // ARGS
  // cso_options  in  scanopt structure to be parsed
  // st_DcostWeight  in/out  Dweight structure to be filled
  //
  // DESCRIPTION
  // Scans the options file structure for Dijkstra weights
  //
  // HISTORY
  // 17 November 2004
  // o Initial design and coding.
  //

  string str_value;

  if (cso_options.scanFor("Dwd", &str_value))
    st_DcostWeight.Dwd = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a Dwd weight.",
               200);
  if (cso_options.scanFor("Dwc", &str_value))
    st_DcostWeight.Dwc = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a Dwc weight.",
               201);
  if (cso_options.scanFor("Dwh", &str_value))
    st_DcostWeight.Dwh = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a Dwh weight.",
               202);
  if (cso_options.scanFor("Dwdc", &str_value))
    st_DcostWeight.Dwdc = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a Dwdc weight.",
               203);
  if (cso_options.scanFor("Dwdh", &str_value))
    st_DcostWeight.Dwdh = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a Dwdh weight.",
               204);
  if (cso_options.scanFor("Dwch", &str_value))
    st_DcostWeight.Dwch = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a Dwch weight.",
               205);
  if (cso_options.scanFor("Dwdch", &str_value))
    st_DcostWeight.Dwdch = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a Dwdch weight.",
               206);
  if (cso_options.scanFor("Dwdir", &str_value))
    st_DcostWeight.Dwdir = atof(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a Dwdir weight.",
               207);
}

#if 1
void
s_weights_print(
  s_weights& asw
) {
  int lw     = 20;
  int rw     = 20;
  //std::_Ios_Fmtflags origFlags;

  //origFlags  = cout.flags();
  cout.setf(ios::left);

  cout << "_________________________"  << endl;
  CW(18, "weight");
  CWn(rw, "value");
  cout << "_________________________"  << endl;
  cout << endl;
  CW(lw, "wd");
  CWn(rw, asw.wd);
  CW(lw, "wc");
  CWn(rw, asw.wc);
  CW(lw, "wh");
  CWn(rw, asw.wh);
  CW(lw, "wdc");
  CWn(rw, asw.wdc);
  CW(lw, "wdh");
  CWn(rw, asw.wdh);
  CW(lw, "wch");
  CWn(rw, asw.wch);
  CW(lw, "wdch");
  CWn(rw, asw.wdch);
  CW(lw, "wdir");
  CWn(rw, asw.wdir);

  //cout.flags(origFlags);
}
#endif

void
s_weights_setAll(
    s_weights&  asw,
    float       af
) {
    asw.wd      = af;
    asw.wc      = af;
    asw.wh      = af;
    asw.wdc     = af;
    asw.wdh     = af;
    asw.wch     = af;
    asw.wdch    = af;
    asw.wdir    = af;
}

#if 1
void
s_Dweights_print(
  s_Dweights& asw
) {
  int lw     = 20;
  int rw     = 20;
  //std::_Ios_Fmtflags origFlags;

  //origFlags  = cout.flags();
  cout.setf(ios::left);

  cout << "_________________________"  << endl;
  CW(18, "Dweight");
  CWn(rw, "value");
  cout << "_________________________"  << endl;
  cout << endl;
  CW(lw, "Dwd");
  CWn(rw, asw.Dwd);
  CW(lw, "Dwc");
  CWn(rw, asw.Dwc);
  CW(lw, "Dwh");
  CWn(rw, asw.Dwh);
  CW(lw, "Dwdc");
  CWn(rw, asw.Dwdc);
  CW(lw, "Dwdh");
  CWn(rw, asw.Dwdh);
  CW(lw, "Dwch");
  CWn(rw, asw.Dwch);
  CW(lw, "Dwdch");
  CWn(rw, asw.Dwdch);
  CW(lw, "Dwdir");
  CWn(rw, asw.Dwdir);

  //cout.flags(origFlags);
}
#endif

void
s_Dweights_setAll(
    s_Dweights&     asw,
    float           af
) {
    asw.Dwd     = af;
    asw.Dwc     = af;
    asw.Dwh     = af;
    asw.Dwdc    = af;
    asw.Dwdh    = af;
    asw.Dwch    = af;
    asw.Dwdch   = af;
    asw.Dwdir   = af;
}

void s_weights_copy(
    s_weights&		sw_target,
    s_weights&		sw_source
) {
    //
    // ARGS
    //	sw_target	in/out		target struct 
    //	sw_source	in/out		source struct
    //
    // DESC
    //	"Deep" copy the component of sw_source to sw_target.
    //

    sw_target.wd	= sw_source.wd;
    sw_target.wc	= sw_source.wc;
    sw_target.wh	= sw_source.wh;
    sw_target.wdc	= sw_source.wdc;
    sw_target.wdh	= sw_source.wdh;
    sw_target.wch	= sw_source.wch;
    sw_target.wdch	= sw_source.wdch;
    sw_target.wdir	= sw_source.wdir;
}

void s_Dweights_copy(
    s_Dweights&		sw_target,
    s_Dweights&		sw_source
) {
    //
    // ARGS
    //	sw_target	in/out		target struct 
    //	sw_source	in/out		source struct
    //
    // DESC
    //	"Deep" copy the component of sw_source to sw_target.
    //

    sw_target.Dwd	= sw_source.Dwd;
    sw_target.Dwc	= sw_source.Dwc;
    sw_target.Dwh	= sw_source.Dwh;
    sw_target.Dwdc	= sw_source.Dwdc;
    sw_target.Dwdh	= sw_source.Dwdh;
    sw_target.Dwch	= sw_source.Dwch;
    sw_target.Dwdch	= sw_source.Dwdch;
    sw_target.Dwdir	= sw_source.Dwdir;
}

void
s_env_costFctList(
  s_env& ast_env
) {
  int  lw = 30;
  int  rw = 20;

  CW(lw, "Current cost function:");
  CWn(rw, ast_env.pstr_functionName[ast_env.ecf_current]);

  cout << "All available cost functions:" << endl;
  for (int i=0; i<ast_env.totalNumFunctions; i++) {
    cout << "Function index: " << i << ": ";
    cout << ast_env.pstr_functionName[i] << endl;
  }
}

int
s_env_costFctSetIndex(
    s_env*      apst_env,
    int         aindex
) {
  int   ret = -1;
  switch (aindex) {
  case 0:
    s_env_costFctSet(  apst_env,
                       costFunc_defaultDetermine,
                       (e_COSTFUNCTION) 0);
    break;
  case 1:
    s_env_costFctSet(  apst_env,
                       costFunc_unityReturn,
                       (e_COSTFUNCTION) 1);
    break;
  case 2:
    s_env_costFctSet(  apst_env,
                       costFunc_EuclideanReturn,
                       (e_COSTFUNCTION) 2);
    break;
  case 3:
    s_env_costFctSet(  apst_env,
                       costFunc_distanceReturn,
                       (e_COSTFUNCTION) 3);
    break;
  default:
    s_env_costFctSet(  apst_env,
                       costFunc_defaultDetermine,
                       (e_COSTFUNCTION) 0);
    break;
  }
  if(aindex < 0 || aindex > 3)
      ret       = -1;
  else
      ret       = aindex;
  return ret;
}

void
s_env_costFctSet(
    s_env*          pst_env,
    float               (*acost_fct) (
        s_env&              st_env,
        s_iterInfo*         pst_iterInfo,
        int                 vno_c,
        int                 j,
        bool                b_relNextReference
    ),
    e_COSTFUNCTION  aecf_new
) {
    pst_env->costFunc_do = acost_fct;
    pst_env->ecf_current = aecf_new;
};

float
costFunc_defaultDetermine(
    s_env&          st_env,
    s_iterInfo*     pst_iterInfo,
    int             vno_c,
    int             j,
    bool            b_relNextReference) {
    //
    // HISTORY
    // 09 November 2004
    // o Added st_iterInfo
    //

    int           vno_n;
    VERTEX*       v_n;
    float         dist, ave_curv, curv, max_height, max_curv, cost;
    s_weights*    pSTw = st_env.pSTw;
    MRIS*         surf = st_env.pMS_primary;

    VERTEX_TOPOLOGY const * const v_ct = &surf->vertices_topology[vno_c];
    VERTEX          const * const v_c  = &surf->vertices         [vno_c];
    if (b_relNextReference) {
        vno_n = v_ct->v[j];
        v_n = &surf->vertices[vno_n];
    } else {
        v_n = &surf->vertices[j];
    }

    float wd  = pSTw->wd;
    float wdc = pSTw->wdc;
    float wc  = pSTw->wc;
    float wdh = pSTw->wdh;
    float wh  = pSTw->wh;
    float wch = pSTw->wch;

    float wdch  = pSTw->wdch;
    float wdir  = pSTw->wdir;

    float f_height      = 0.;
    float f_dir         = 0.;

    st_V3D        V3_c;           // current point
    st_V3D        V3_n;           // next points
    st_V3D        V3_cn;          // vector from current to next
    static st_V3D V3_e;           // end point
    st_V3D        V3_ce;          // vector from current to end
    static int    calls = 0;

    if (!calls) {
        V3_e.f_x = st_env.pMS_primary->vertices[st_env.endVertex].x;
        V3_e.f_y = st_env.pMS_primary->vertices[st_env.endVertex].y;
        V3_e.f_z = st_env.pMS_primary->vertices[st_env.endVertex].z;
    }
    calls++;

    // Cartesian points "current" and "next"
    V3_c.f_x = v_c->x;
    V3_n.f_x = v_n->x;
    V3_c.f_y = v_c->y;
    V3_n.f_y = v_n->y;
    V3_c.f_z = v_c->z;
    V3_n.f_z = v_n->z;

    // direction vectors
    V3D_normalizedDirection_find(V3_c, V3_n, &V3_cn);
    V3D_normalizedDirection_find(V3_c, V3_e, &V3_ce);

    float f_dot = V3D_dot(V3_cn, V3_ce);        // dot product is "1"
                                                //+ for colinear
                                                //+ vectors
    f_dir  = 1 - f_dot;

    if (b_relNextReference)
        dist = v_c->dist[j];
    else
        dist = V3D_distance(V3_c, V3_n);

    // Arguably two possibilities for thinking about the
    // curvature connecting the current vertex to its
    // neighbour. Initially I thought to take a simple
    // average between the current and neighbour vertex
    // curvatures, but this has 'artifacts' across zero
    // crossings when the curvature changes sign. I then
    // simply just took the neighbour curvature as the
    // curvature between 'here' and 'there' to avoid
    // all these sign-issues.

    // ave_curv  = (v_c->curv + v_n->curv) / 2.0;
    ave_curv                     = v_n->curv;

    pst_iterInfo->iter           = calls;
    pst_iterInfo->f_distance     = dist;
    pst_iterInfo->f_curvature    = ave_curv;
    pst_iterInfo->f_sulcalHeight = st_env.pMS_secondary->vertices[vno_c].curv;
    pst_iterInfo->f_dir          = f_dir;

    // Initial testing revealed that 'wdch' was particularly sensitive to *=10,
    //  and resulted in considerable (negative) impact on overall
    // trajectory

    if (st_env.b_transitionPenalties && ave_curv<0) {
        wc      *= st_env.pSTDw->Dwc;
        wdc     *= st_env.pSTDw->Dwdc;
        wch     *= st_env.pSTDw->Dwch;
        wdch    *= st_env.pSTDw->Dwdch;
    }
    if (st_env.b_transitionPenalties &&  f_height<0) {
        wh      *= st_env.pSTDw->Dwh;
        wdh     *= st_env.pSTDw->Dwdh;
        wch     *= st_env.pSTDw->Dwch;
        wdch    *= st_env.pSTDw->Dwdch;
    }

    max_height  = (st_env.pMS_secondary->max_curv);
    max_curv    = (st_env.pMS_primary->max_curv);
    if(st_env.b_useAbsCurvs) {
        f_height        = fabs(f_height);
        curv            = fabs(ave_curv);
    } else {
        f_height        = max_height - st_env.pMS_secondary->vertices[vno_c].curv;
        curv            = max_curv   - ave_curv;
    }
    // cost   = dist + 20.0 * dist * curv;
    cost  = wd*dist                     + wc*curv               +
            wh*f_height                 + wdc*dist*curv         +
            wdh*dist*f_height           + wch*curv*f_height     +
            wdch*dist*curv*f_height     + wdir*f_dir;

    return(cost);
}

float
costFunc_unityReturn(
    s_env&          st_env,
    s_iterInfo*     pst_iterInfo,
    int             vno_c,
    int             j,
    bool            b_relNextReference) {
    //
    // POSTCONDITIONS
    //  o Will always return a 1.0 as the transition cost. This is used 
    //    primarily in determining logical distances between nodes in the 
    //    vertex.
    //
    //   Most of the function arguments are superfluous in this case.
    //
    // HISTORY
    // 15 February 2005
    // o Initial development
    //

    float   cost = 1.0;
    return(cost);
}

float
costFunc_distanceReturn(
    s_env&          st_env,
    s_iterInfo*     pst_iterInfo,
    int             vno_c,
    int             j,
    bool            b_relNextReference)
{
    //
    // PRECONDITIONS
    // o Distance concept is only valid for "relative" neighbours.
    //
    // POSTCONDITIONS
    //  o Returns the weighted distance as stored in the MRIS.
    //
    //   Most of the function arguments are superfluous in this case.
    //
    // HISTORY
    // 09 March 2005
    // o Initial development
    //
    // 24 February 2011
    // o For "absolute" references, determine the relative neighbor.
    //

    float       	f_cost      = 0.0;
    float       	f_distance  = 0.0;
    s_weights*  	pSTw        = st_env.pSTw;
    float       	wd          = pSTw->wd;
    MRIS*    		surf        = st_env.pMS_primary;
    const char*       	pch_proc    = "costFunc_distanceReturn(...)";
    char        	pch_txt[65536];
    static bool 	b_warned    = false;

    VERTEX_TOPOLOGY const * const v_ct = &surf->vertices_topology[vno_c];
    VERTEX          const * const v_c  = &surf->vertices         [vno_c];

    if(!b_relNextReference) {
	int 	jrelcount;
	int 	jrel = 0;
	for(jrelcount=0; jrelcount< v_ct->vnum; jrelcount++) {
	    if(v_ct->v[jrelcount] == j) {
		jrel = jrelcount;
		break;
	    }
	}
	j = jrel;
    }
    f_distance  = v_c->dist[j];

    f_cost  = f_distance * wd;
    if (wd <= 0.) {
        sprintf(pch_txt, "calculating cost in %s", pch_proc);
        error_exit(pch_txt, "wd must be greater than zero.", 1);
    }
    if (wd != 1. && !b_warned) {
        sprintf(pch_txt, "calculating cost in %s", pch_proc);
        warn(pch_txt, "wd is not equal to 1. Distances will be skewed", 1);
        b_warned = true;
    }

    return(f_cost);
}

float
costFunc_EuclideanReturn(
    s_env&          st_env,
    s_iterInfo*     pst_iterInfo,
    int             vno_c,
    int             j,
    bool            b_relNextReference) {
    //
    // POSTCONDITIONS
    //  o Returns the Euclidean distance between a vertex and its neighbour.
    //   This distance is weighted by 'w_d'.
    //
    //   Most of the function arguments are superfluous in this case.
    //
    // HISTORY
    // 09 March 2005
    // o Initial development
    //

    float       f_cost      = 0.0;
    float       f_distance  = 0.0;
    int         vno_n       = 0;

    VERTEX*     v_n         = NULL;
    MRIS*       surf        = st_env.pMS_primary;
    static int  calls       = 0;
    const char*       pch_proc    = "costFunc_EuclideanReturn(...)";
    char        pch_txt[65536];
    static bool b_warned = false;

    VERTEX_TOPOLOGY const * const v_ct = &surf->vertices_topology[vno_c];
    VERTEX          const * const v_c  = &surf->vertices         [vno_c];

    if (b_relNextReference) {
        vno_n = v_ct->v[j];
        v_n = &surf->vertices[vno_n];
    } else {
        v_n = &surf->vertices[j];
    }

    s_weights*  pSTw  = st_env.pSTw;
    float  wd  = pSTw->wd;

    if (wd <= 0.) {
        sprintf(pch_txt, "calculating cost in %s", pch_proc);
        error_exit(pch_txt, "wd must be greater than zero.", 1);
    }
    if (wd != 1. && !b_warned) {
        sprintf(pch_txt, "calculating cost in %s", pch_proc);
        warn(pch_txt, "wd is not equal to 1. Distances will be skewed", 1);
        b_warned = true;
    }

    st_V3D   V3_c; // current point
    st_V3D   V3_n; // next points

    calls++;

    // Cartesian points "current" and "next"
    V3_c.f_x    = v_c->x;
    V3_n.f_x    = v_n->x;
    V3_c.f_y    = v_c->y;
    V3_n.f_y    = v_n->y;
    V3_c.f_z    = v_c->z;
    V3_n.f_z    = v_n->z;

    f_distance  = V3D_distance(V3_c, V3_n);

    f_cost      = f_distance * wd;

    return(f_cost);
}

#if 0
bool
asynchEvent_processWGHT(
  s_env&    ast_env,
  string    astr_comms
) {

  int lw     = ast_env.lw;
  int rw     = ast_env.rw;

  string str_errorAct = "checking <WGHT>";

  string str_object = "";
  string str_verb = "";
  string str_modifier = "";
  string str_sep  = " ";
  float f_val  = 0.0;

  //std::_Ios_Fmtflags origFlags;
  //origFlags  = cout.flags();
  cout.setf(ios::left);

  if (!str_3parse( astr_comms, str_object, str_verb, str_modifier))
    warn(str_errorAct, "Some error occurred in the 3parse.", 1);

  if (str_object == "all") {
    if (str_verb == "get") {
      s_weights_print(*ast_env.pSTw);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      s_weights_setAll(*ast_env.pSTw, f_val);
    }
  }
  if (str_object == "wd") {
    if (str_verb == "get") {
      CW(lw, "wd");
      CWn(rw, ast_env.pSTw->wd);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTw->wd =  f_val;
    }
  }
  if (str_object == "wc") {
    if (str_verb == "get") {
      CW(lw, "wc");
      CWn(rw, ast_env.pSTw->wc);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTw->wc =  f_val;
    }
  }
  if (str_object == "wh") {
    if (str_verb == "get") {
      CW(lw, "wh");
      CWn(rw, ast_env.pSTw->wh);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTw->wh =  f_val;
    }
  }
  if (str_object == "wdc") {
    if (str_verb == "get") {
      CW(lw, "wdc");
      CWn(rw, ast_env.pSTw->wdc);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTw->wdc =  f_val;
    }
  }
  if (str_object == "wdh") {
    if (str_verb == "get") {
      CW(lw, "wdh");
      CWn(rw, ast_env.pSTw->wdh);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTw->wdh =  f_val;
    }
  }
  if (str_object == "wch") {
    if (str_verb == "get") {
      CW(lw, "wch");
      CWn(rw, ast_env.pSTw->wch);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTw->wch =  f_val;
    }
  }
  if (str_object == "wdch") {
    if (str_verb == "get") {
      CW(lw, "wdch");
      CWn(rw, ast_env.pSTw->wdch);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTw->wdch =  f_val;
    }
  }
  if (str_object == "wdir") {
    if (str_verb == "get") {
      CW(lw, "wdir");
      CWn(rw, ast_env.pSTw->wdir);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTw->wdir =  f_val;
    }
  }

  //cout.flags(origFlags);
  return true;
}
#endif

#if 0
bool
asynchEvent_processDWGHT(
  s_env&    ast_env,
  string    astr_comms
) {

  int lw     = ast_env.lw;
  int rw     = ast_env.rw;

  string str_errorAct = "checking <DWGHT>";

  string str_object = "";
  string str_verb = "";
  string str_modifier = "";
  string str_sep  = " ";
  float f_val  = 0.0;

  //std::_Ios_Fmtflags origFlags;
  //origFlags  = cout.flags();
  cout.setf(ios::left);

  if (!str_3parse( astr_comms, str_object, str_verb, str_modifier))
    warn(str_errorAct, "Some error occurred in the 3parse.", 1);

  if (str_object == "all") {
    if (str_verb == "get") {
      s_Dweights_print(*ast_env.pSTDw);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      s_Dweights_setAll(*ast_env.pSTDw, f_val);
    }
  }
  if (str_object == "Dwd") {
    if (str_verb == "get") {
      CW(lw, "Dwd");
      CWn(rw, ast_env.pSTDw->Dwd);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTDw->Dwd =  f_val;
    }
  }
  if (str_object == "Dwc") {
    if (str_verb == "get") {
      CW(lw, "Dwc");
      CWn(rw, ast_env.pSTDw->Dwc);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTDw->Dwc =  f_val;
    }
  }
  if (str_object == "Dwh") {
    if (str_verb == "get") {
      CW(lw, "Dwh");
      CWn(rw, ast_env.pSTDw->Dwh);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTDw->Dwh =  f_val;
    }
  }
  if (str_object == "Dwdc") {
    if (str_verb == "get") {
      CW(lw, "Dwdc");
      CWn(rw, ast_env.pSTDw->Dwdc);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTDw->Dwdc =  f_val;
    }
  }
  if (str_object == "Dwdh") {
    if (str_verb == "get") {
      CW(lw, "Dwdh");
      CWn(rw, ast_env.pSTDw->Dwdh);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTDw->Dwdh =  f_val;
    }
  }
  if (str_object == "Dwch") {
    if (str_verb == "get") {
      CW(lw, "Dwch");
      CWn(rw, ast_env.pSTDw->Dwch);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTDw->Dwch =  f_val;
    }
  }
  if (str_object == "Dwdch") {
    if (str_verb == "get") {
      CW(lw, "Dwdch");
      CWn(rw, ast_env.pSTDw->Dwdch);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTDw->Dwdch =  f_val;
    }
  }
  if (str_object == "Dwdir") {
    if (str_verb == "get") {
      CW(lw, "Dwdir");
      CWn(rw, ast_env.pSTDw->Dwdir);
    } else if (str_verb == "set") {
      if (!str_modifier.length())
        return false;
      f_val = atof(str_modifier.c_str());
      ast_env.pSTDw->Dwdir =  f_val;
    }
  }

  //cout.flags(origFlags);
  return true;
}
#endif

/* eof */
