/**
 * @brief The environment object API.
 *
 * The environment acts as the main access / interface between the core
 * components of the program and the external system. Several structures
 * and objects, taken together, constitute this environment and define
 * amongst others, operational parameters, cost functions, weight
 * structures, etc.
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

#ifndef __ENV_H__
#define __ENV_H__

#include <sys/stat.h>

//#include "legacy.h"
//#include "general.h"

#include "c_SMessage.h"
#include "c_SSocket.h"
#include "scanopt.h"

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "error.h"

#include <string>
#include <vector>
using namespace std;

// Forward declaration
class C_mpmProg;
class C_mpmOverlay;

// For legacy handling...
typedef struct _weights 	s_weights;
typedef struct _Dweights	s_Dweights;

/// s_iterInfo contains information pertaining to a particular iteration
/// of the main program loop. It is accessed during each call to the
/// cost function, and is populated with per-call information.
///
typedef struct _iterInfo {
  int       iter;
  float     f_distance;
  float     f_curvature;
  float     f_sulcalHeight;
  float     f_dir;
}
s_iterInfo;

typedef enum {
  e_default, e_unity, e_euclid, e_distance
} e_COSTFUNCTION;

#if 0
/// Weights and values for the main cost function polynomial.
typedef struct _weights {
  float wd;         // distance
  float wc;         // curve
  float wh;         // sulcal height
  float wdc;        // distance X curve
  float wdh;        // distance X height
  float wch;        // curve X height
  float wdch;       // distance X curve X height
  float wdir;       // direction vector
}
s_weights;

void
s_weights_scan(
  s_weights&  st_costWeight,
  C_scanopt&  cso_options
);

/// \fn void s_weights_print(s_weights& asw)
/// \brief Print the current internal weight structure to stdout.
/// \param  asw The weight structure to print
/// \return    (void)
void s_weights_print( s_weights& asw);

/// \fn void s_weights_setAll(s_weights& asw, float af)
/// \brief Sets all the weights in the current internal weight structure to the passed value.
/// \param  asw The weight structure to access
/// \param af Value assigned to each member
/// \return    (void)
void  s_weights_setAll( s_weights& asw,
                        float  af);
/// Del weights and values.
///
/// These are used when zero-crossings
/// of curvature and sulcal height occur. By incorporating these
/// transition weights, we can force the path search algorithm
/// to choose to strongly avoid areas of zero-crossing.
typedef struct _Dweights {
  float Dwd;            // distance
  float Dwc;            // curve
  float Dwh;            // sulcal height
  float Dwdc;           // distance X curve
  float Dwdh;           // distance X height
  float Dwch;           // curve X height
  float Dwdch;          // distance X curve X height
  float Dwdir;          // direction vector
}
s_Dweights;

void
s_Dweights_scan(
  s_Dweights&  st_DcostWeight,
  C_scanopt&  cso_options
);

/// \fn void s_Dweights_print(s_Dweights& asw)
/// \brief Print the current internal Dweight structure to stdout.
/// \param  asw The weight structure to print
/// \return    (void)
void s_Dweights_print( s_Dweights& asw);

/// \fn void s_Dweights_setAll(s_Dweights& asw, float af)
/// \brief Sets all the delta weights in the current internal Dweight structure to the passed value.
/// \param  asw The Dweight structure to access
/// \param af Value assigned to each member
/// \return    (void)
void  s_Dweights_setAll( s_Dweights& asw,
                         float  af);

/// s_iterInfo contains information pertaining to a particular iteration
/// of the main program loop. It is accessed during each call to the
/// cost function, and is populated with per-call information.
///
typedef struct _iterInfo {
  int       iter;
  float     f_distance;
  float     f_curvature;
  float     f_sulcalHeight;
  float     f_dir;
}
s_iterInfo;

typedef enum {
  e_default, e_unity, e_euclid, e_distance
} e_COSTFUNCTION;
#endif

typedef enum {
  e_workingCurvature, e_workingSulcal, e_auxillary
} e_SURFACE;

typedef enum {
  e_user, e_sys, e_result
} e_LOG;

/// The main environment structure. This structure records important
/// variables, mostly interpreted from the process <b>options</b> file.
typedef struct _env s_env;

/// 
/// mpm MODULE enums
///

typedef enum _e_MODULE {
	e_mpmProg, e_mpmOverlay, emodule
} e_MODULE;

// String names for these are defined in env.cpp using 'push_back()'
// To add new modules:
//      o Create the index here
//      o Edit the s_env_mpmProgSetIndex(...)
//
// Internal module name check on passed command line args:
//      o help::mpmProg_check()

typedef enum _e_mpmProg {
    emp_NULL 		= 0, 
    emp_NOP 		= 1,
    emp_pathFind	= 2,
    emp_autodijk 	= 3, 
    emp_autodijk_fast 	= 4,
    emp_ROI             = 5,
    emp_externalMesh	= 6,
    empmprog
} e_MPMPROG;

// enum typedef for mpmOverlays
// ordering should be same as e_COSTFUNCTION
//
// String names for these are defined in env.cpp using 'push_back()'
//      o env::s_env_nullify()
//
// Add case handling to:
//      o asynch:asynchEvent_processMRMPROG()
//      o env::s_env_mpmProgSetIndex()
//      o help::mpmOverlay_check()
//
typedef enum _e_mpmOverlay {
    emo_LEGACY          = 0,    // LEGACY overlay -- toggles cost calc to legacy
                                // engine.
    emo_NULL 		= 1,	// NULL overlay -- for debugging
    emo_NOP		= 2,	// NOP overlay -- for debugging 
    emo_unity		= 3,	// returns '1' for each internode distance	
    emo_euclidean	= 4,	// returns distance between nodes (calculated)
    emo_distance	= 5,	// returns distance between nodes (read)
    emo_fscurvs		= 6,	// returns weighted cost function of curvs
    emo_curvature       = 7,    // returns curvature between nodes
    empmoverlay
} e_MPMOVERLAY;


///
/// Main environment structure
///

typedef struct _env {
    int           argc;                     // original number of command line
                                            // args
    char**        ppch_argv;                // char array of command line args
    
    C_scanopt*	  pcso_options;		    // class that houses the parsed
    					    //+ options file for the environment
    C_SMessage*   pcsm_optionsFile;         // message wrapper for options file
    					    //+ used to create a new options
    					    //+ file.
    
    int           timeoutSec;               // listen timeout
    int           port;                     // port on which to listen
                                            //+ for async control

    int           lw;                       // left width (for stdout format)
    int           rw;                       // right width (for stdout format)

    bool          b_syslogPrepend;          // prepend syslog style
    string        str_stdout;
    string        str_userMsgLog;
    string        str_sysMsgLog;
    string        str_resultMsgLog;
    C_SMessage*   pcsm_stdout;              // stdout C_SMessage object
    C_SMessage*   pcsm_syslog;              // log file for "sys" events
    C_SMessage*   pcsm_userlog;             // log file for "user" events
    C_SMessage*   pcsm_resultlog;           // log file for "result" event
    int           serverControlPort;        // port on which internal server 
                                            //+ listens

    bool          b_labelFile_save;         // flag: save label file
    bool          b_patchFile_save;         // flag: save patch file?
    bool          b_transitionPenalties;    // flag: apply transition?
                                            // penalties
    bool          b_useAbsCurvs;            // flag: if true, use abs(curv)
                                            //+ for cost function calculations
    string        str_workingDir;           // directory containing input
                                            //+ and output files
    string        str_patchFileStem;
    string        str_labelFileStem;
    string        str_patchFileName;        // file to contain path "patch"
    string        str_labelFileName;        // file to contain path "label"
    string        str_labelFileNameOS;      // aux file to contain path "label"
                                            //+ projected back onto
                                            //+ the "Original Surface"
    bool          b_optionsFileUse;         // Initialize from options file
                                            //+ or command line arguments
    string        str_optionsFileName;      // file containing meta options
    string        str_costFileName;         // file containing final costs
    bool          b_costPathSave;           // toggle for creating the costFile
                                            //+ by turning OFF, extra speed
                                            //+ is available for cases when
                                            //+ system running in tight
                                            //+ loops, typically when executing
                                            //+ mpmProgs

    int           startVertex;              // index of the start vertex
    int           endVertex;                // index of the end vertex
    float         plyDepth;                 // ply depth around a core path

    e_SURFACE     esf_active;               // the current active surface
    string*       pstr_activeName;          // names of each surface
    int           totalNumSurfaces;         // total number of surfaces
    MRIS*         pMS_active;               // a pointer (used by several
                                            // functions) to specifiy
                                            // a particular surface to
                                            // process
    string        str_subject;
    string        str_hemi;
    string        str_surface;

    MRIS*         pMS_primary;              // primary surface
    string        str_primarySurfaceFileName;
    string        str_primaryCurvatureFileName;
    bool          b_primaryCurvature;

    MRIS*         pMS_secondary;            // secondary (optional) surface
    string        str_secondarySurfaceFileName;
    bool          b_secondarySurface;
    string        str_secondaryCurvatureFileName;
    bool          b_secondaryCurvature;

    MRIS*         pMS_auxillary;            // auxillary surface

    bool          b_surfacesKeepInSync;     // flag: behavioural /
                                            //+ conditional. Evaluate
                                            //+ if wanting to merge
                                            //+ information from one
                                            //+ surface to another.

    bool          b_surfacesClear;          // flag: behavioural /
                                            //+ conditional. Should be
                                            //+ true for most cases.
                                            //+ If false, prevents the
                                            //+ clearing of rips
                                            //+ after a path search.
                                            //+ Useful when keeping
                                            //+ a pattern for cases
                                            //+ when correlation needs
                                            //+ to be calculated.

    bool          b_costHistoryPreserve;    // flag: preserve cost history
                                            //+ on surface between
                                            //+ successive calls to
                                            //+ dijkstra function.
                                            //+ Set to TRUE if finding
                                            //+ ply distances from
                                            //+ existing path

    //
    // LEGACY CODE
    s_weights*    pSTw;                     // weight structure
    s_Dweights*   pSTDw;                    // Del weight structure

    int           totalNumFunctions;        // total number of cost
                                            // functions
    s_iterInfo*	  pst_iterInfo;	            // structure that houses
    					    //+ per-iteration information
    e_COSTFUNCTION ecf_current;             // the current cost function
    string*       pstr_functionName;        // names of each cost function
    float  (*costFunc_do)                   // a cost function to operate
    (                                       //+ on this environment
        s_env&          st_env,
        s_iterInfo*     pst_iterInfo,
        int             vno_c,
        int             j,
        bool            b_relNextReference
    );
    // LEGACY CODE
    //

    //
    // Modules
    //

    vector<string>	vstr_mpm;	    // contains list of all module
    					    // type names
    
    // mpmProgs
    int                 totalmpmProgs;      // total number of mpmProgs
    vector<string>	vstr_mpmProgName;   // names of each mpmProg
    bool                b_mpmProgUse;       // flag toggle on using mpm's
    e_MPMPROG           empmProg_current;   // mpm program index to run
    C_mpmProg*          pCmpmProg;          // handle to mpmProg object to run
    string              str_mpmArgs;        // User spec'd, semi-colon delimited
                                            //+ arg string
    bool                b_exitOnDone;       // if true, terminate the main
                                            //+ mris_pmake process when an
                                            //+ mpmProg is finished.
    // autodijk options
    string              str_costCurvFile;   // file containing per vertex costs
                                            //+ for 'autodijk'

    // mpmOverlays
    bool		b_mpmOverlayUse;    // debugging flag to maintain
    					    //+ legacy compatibility. If true,
    					    //+ use overlay engine, else use
    					    //+ legacy engine. This flag will
    					    //+ probably go away (along with the
    					    //+ old engine code!)
    int                 totalmpmOverlays;   // total number of mpmOverlays
    vector<string>	vstr_mpmOverlayName;// names of each mpmOverlay
    e_MPMOVERLAY	empmOverlay_current;// mpmOverlay program index to use
    C_mpmOverlay*       pCmpmOverlay;       // handle to mpmProg object to use
    string              str_mpmOverlayArgs; // User spec'd, semi-colon delimited
                                            //+ arg string
} s_env;

float 
s_env_edgeCostFind(
    s_env&		ast_env,
    int			avertexi,
    int			avertexj
);

void
s_env_mpmPrint(
    s_env&		ast_env,
    string 		astr_msg	= "",
    e_MODULE		ae_module	= e_mpmProg
); 

void
s_env_defaultsSet(
    s_env&              st_env
);

void
s_env_optionsFile_write(
    s_env&              st_env,
    bool                ab_setToDefaults        = false                        
);

void
s_env_scan(
    s_env&              st_env
);

string 
s_env_HUP(
    s_env&			st_env,
    c_SSocket_UDP_receive**    	pCSSocketReceive
);

/// \fn void s_env_nullify( s_env& st_env);
/// \brief Initialise (nullify, i.e. primitive constructor) a passed environment.
/// \param  st_env The environment structure
/// \return    (void)
void s_env_nullify( s_env& st_env);

bool s_env_b_surfacesKeepInSync_get(
    s_env&      ast_env
);

void s_env_b_surfacesKeepInSync_set(
    s_env&      ast_env,
    bool        b_val
);

bool s_env_b_surfacesClear_get(
    s_env&   ast_env
);

void s_env_b_surfacesClear_set(
    s_env&      ast_env,
    bool        b_val
);

float s_env_plyDepth_get(
    s_env&   ast_env
);

void s_env_plyDepth_set(
    s_env&      ast_env,
    float       af_val
);

float s_env_plyIncrement_get(
    s_env&   ast_env
);

void s_env_plyIncrement_set(
    s_env&      ast_env,
    float       af_val
);


bool s_env_surfaceFile_set(
    s_env&      st_env,
    string      astr_fileName
);

bool s_env_surfaceCurvature_set(
    s_env&      st_env,
    string      astr_fileName
);

bool s_env_secondarySurface_setCurvature(
    s_env&      st_env,
    string      astr_fileName
);

bool s_env_auxSurfaceFile_set(
    s_env&      st_env,
    string      astr_fileName
);

bool s_env_auxSurfaceCurvature_set(
    s_env&      st_env,
    string      astr_fileName
);

void  s_env_log_file_changeTo(
    s_env&      ast_env,
    e_LOG       ae_log,
    string      astr_newName
);

void s_env_activeSurfaceList(
    s_env&      ast_env
);

void s_env_activeSurfaceSetIndex(
    s_env*      apst_env,
    int         aindex
);

int s_env_mpmProgSetIndex(
    s_env*      apst_env,
    int         aindex
);

int
s_env_mpmOverlaySetIndex(
    s_env*      apst_env,
    int         aindex
);

#if 0
void s_env_costFctList(
    s_env&      ast_env
);

int s_env_costFctSetIndex(
    s_env*      apst_env,
    int         aindex
);

void  s_env_costFctSet(
    s_env*      pst_env,
    float       (*cost_fct) (
        s_env&          st_env,
        s_iterInfo*     pst_iterInfo,
        int             vno_c,
        int             j,
        bool            b_relNextReference
    ),
    e_COSTFUNCTION  aecf_new = e_default
);

/// \fn float costFunc_defaultDetermine(s_env& st_env, s_iterInfo* pst_iterInfo, int vno_c, int j, bool b_relNextReference = true);
/// \brief Calculate the cost in moving from one vertex to another.
/// \param  st_env The environment structure
/// \param pst_iterInfo A structure that records vertex transition information (cost, distance, curvature, etc.)
/// \param  vno The current vertex index
/// \param  j The neighbour vertex index - either relative or absolute
/// \param b_relNextReference If true, implies that 'j' is a reference relative to vno, otherwise j is an absolute reference
/// \return    (float) cost of the transition
float   costFunc_defaultDetermine(
  s_env&  st_env,
  s_iterInfo* pst_iterInfo,
  int   vno_c,
  int   j,
  bool  b_relNextReference = true
);

/// \fn float costFunc_unityReturn(s_env& st_env, s_iterInfo* pst_iterInfo, int vno_c, int j, bool b_relNextReference = true);
/// \brief Always returns a '1' as the transition cost between a node and its neighbour.
/// \param  st_env The environment structure
/// \param pst_iterInfo A structure that records vertex transition information (cost, distance, curvature, etc.)
/// \param  vno The current vertex index
/// \param  j The neighbour vertex index - either relative or absolute
/// \param b_relNextReference If true, implies that 'j' is a reference relative to vno, otherwise j is an absolute reference
/// \return    (float) cost of the transition - always 1
float   costFunc_unityReturn(
  s_env&  st_env,
  s_iterInfo* pst_iterInfo,
  int   vno_c,
  int   j,
  bool  b_relNextReference = true
);

/// \fn float costFunc_EuclideanReturn(s_env& st_env, s_iterInfo* pst_iterInfo, int vno_c, int j, bool b_relNextReference = true);
/// \brief Returns the Euclidean distance between a vertex and its neighbour.
/// \param  st_env The environment structure
/// \param pst_iterInfo A structure that records vertex transition information (cost, distance, curvature, etc.)
/// \param  vno The current vertex index
/// \param  j The neighbour vertex index - either relative or absolute
/// \param b_relNextReference If true, implies that 'j' is a reference relative to vno, otherwise j is an absolute reference
/// \return    (float) cost of the transition - always 1
float   costFunc_EuclideanReturn(
  s_env&  st_env,
  s_iterInfo* pst_iterInfo,
  int   vno_c,
  int   j,
  bool  b_relNextReference = true
);

/// \fn float costFunc_distanceReturn(s_env& st_env, s_iterInfo* pst_iterInfo, int vno_c, int j, bool b_relNextReference = true);
/// \brief Returns the distance (as encoded in the MRIS itself) between a vertex and its neighbour. For a given vertex, the
///  distance between itself and neighbouring vertices is recorded in the MRIS data structure itself.
/// \param  st_env The environment structure
/// \param pst_iterInfo A structure that records vertex transition information (cost, distance, curvature, etc.)
/// \param  vno The current vertex index
/// \param  j The neighbour vertex index - either relative or absolute
/// \param b_relNextReference If true, implies that 'j' is a reference relative to vno, otherwise j is an absolute reference
/// \return    (float) cost of the transition - always 1
float   costFunc_distanceReturn(
  s_env&  st_env,
  s_iterInfo* pst_iterInfo,
  int   vno_c,
  int   j,
  bool  b_relNextReference = true
);
#endif

#endif //__ENV_H__


