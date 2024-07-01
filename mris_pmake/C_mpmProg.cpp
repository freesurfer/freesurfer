/*
 * @brief The internal 'program' API.
 *
 *  C_mpmProgs are overloaded classes that perform specific functions in the
 *  context of the Dijkstra system. The default contol system allows an
 *  external network-based program to drive the core dijkstra search. This
 *  process can be very slow on large-scale problems. To alleviate that, 
 *  mpm_programs can be written. These map directly to external dsh scripts
 *  but without any of the network overhead.
 */
/*
 * Original Author: Rudolph Pienaar
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


#include "C_mpmProg.h"
#include "dijkstra.h"

#include "c_surface.h"
#include "c_label.h"
#include "c_vertex.h"
#include "fio.h"

#include "unistd.h"
#include <libgen.h>

#include <sstream>

extern 	bool     	Gb_stdout;
extern	stringstream	Gsout;

//
//\\\---
// Base class issues --->>>
/////---
//

//-----------------------------------------------------------
// Base socket constructors/destructor, output dump routines
//-----------------------------------------------------------

void
C_mpmProg::debug_push(
    string      astr_currentProc)
{
    if (stackDepth_get() >= SSTACKDEPTH-1) {
        cout << "Current stackDepth:\t" << stackDepth_get()         << endl;
        cout << "stackdepth limit:\t"   << SSTACKDEPTH              << endl;
        for (int i=0; i<SSTACKDEPTH; i++)
            cout << "Stack depth: " << i << "\t" << str_proc_get(i) << endl;
        error("Out of str_proc stack depth");
    }
    stackDepth_set(stackDepth_get()+1);
    str_proc_set(stackDepth_get(), astr_currentProc);
}

void C_mpmProg::debug_pop() {
    mstackDepth--;
}

void
C_mpmProg::error(
    string      astr_msg        /* = ""  Error message */,
    int         code            /* = 1  Error ID  */)
{

    // comment this line out - it is causing compiler errors for rocky8 gcc10 
    //extern int errno;

    cerr << "\nFatal error encountered.\n";
    cerr << "\tC_mpmProg `" << mstr_name << "' (id: " << mid << ")\n";
    cerr << "\tCurrent function: " << mstr_obj << "::" << str_proc_get() << "\n";
    cerr << "\t" << astr_msg << "\n";
    perror("Error returned from system");
    print();
    cerr << "Throwing exception to (this) with code " << code << "\n\n";
    throw(this);
}

void
C_mpmProg::warn(
    string      astr_msg        /* = ""  Warning message */,
    int         code            /* = 1  Warning code  */)
{
    if (mwarnings) {
        cerr << "\nWarning.\n";
        cerr << "\tC_mpmProg `" << mstr_name << "' (id: " << mid << ")\n";
        cerr << "\tCurrent function: " << mstr_obj << "::" << str_proc_get() << "\n";
        cerr << "\t" << astr_msg << " (warning code: " << code << ")\n";
    }
}

void
C_mpmProg::function_trace(
    string      astr_msg        /* = ""  Trace message */,
    string      astr_separator  /* = ""  Separator message */)
{
    int                 i;
    string              str_tab                 = "";
    static string       str_objectName          = "";
    static string       str_funcName            = "";

    if (mverbosity>=mstackDepth) {
        cerr << astr_separator;
        for (i=0; i<mstackDepth; i++)
        str_tab += "    ";
        if (str_objectName != mstr_name) {
	    cerr << "\n" << mstr_obj << " `";
            cerr << str_name_get() << "' (id: " << mid << ")" << endl;
	}
        if (str_funcName != str_proc_get()) {
            cerr << "\n" << str_tab << "Current function: " << mstr_obj << "::";
            cerr << str_proc_get();
        }
        cerr << "\n" << str_tab << astr_msg << endl;
    }
    str_objectName      = str_name_get();
    str_funcName        = str_proc_get();
}

void
C_mpmProg::core_construct(
    string      astr_name       /* = "unnamed" */,
    int         a_id            /* = -1  */,
    int         a_verbosity     /* = 0  */,
    int         a_warnings      /* = 0  */,
    int         a_stackDepth    /* = 0  */,
    string      astr_proc       /* = "noproc" */)
{
    //
    // ARGS
    //  astr_name       in      name of object
    //  a_id            in      internal id number of object
    //  a_verbosity     in      verbosity level
    // a_warnings       in      warnings level
    // a_stackDepth     in      stack depth (debugging)
    // astr_proc        in      currently executing proc (debugging)
    //
    // DESC
    //  Common core statements for all constructors.
    //

    mstr_name   = astr_name;
    mid         = a_id;
    mverbosity  = a_verbosity;
    mwarnings   = a_warnings;
    mstackDepth = a_stackDepth;

    str_proc_set(stackDepth_get(), "no name");
    mstr_obj = "C_mpmProg";
}

C_mpmProg::C_mpmProg(
    s_env*      aps_env
){
    string      str_name        =       "Unnamed C_mpmProg";
    int         id              = 0;
    core_construct(str_name, id);
    mps_env                     = aps_env;
}

C_mpmProg::C_mpmProg(const C_mpmProg &C_mpmProg) {
    //
    // Copy constructor
    //

    debug_push("C_mpmProg (copy constructor)");
    error("Copy constructor not yet defined");
    debug_pop();
}

C_mpmProg & C_mpmProg::operator=(const C_mpmProg & C_mpmProg) {
    //
    // Overloaded (=) operator
    //

    debug_push("operator=");
    error("Overloaded operator= not yet defined");
    return *this;       // Never executed --
    debug_pop();
}

C_mpmProg::~C_mpmProg(void) {}

void
C_mpmProg::print() {
    cout << "object name:\t"    << mstr_name    << endl;
    cout << "object id:\t"      << mid          << endl;
    cout << "object type:\t"    << mstr_obj     << endl;
}

//
//\\\***
// C_mpmProg_NOP definitions ****>>>>
/////***
//

C_mpmProg_NOP::C_mpmProg_NOP(
    s_env*      aps_env) : C_mpmProg(aps_env)
{
    //
    // ARGS
    //
    // DESC
    // Basically a thin "fall-through" constructor to the base
    // class.
    //
    // PRECONDITIONS
    // o aps_env must be fully instantiated.
    //
    // HISTORY
    // 14 December 2009
    // o Initial design and coding.
    //

    debug_push("C_mpmProg_NOP");
    mstr_obj	= "C_mpmProg_NOP";

    msleepSeconds       = 5;

    debug_pop();
}

C_mpmProg_NOP::~C_mpmProg_NOP() {
    //
    // Destructor
    //

}

int
C_mpmProg_NOP::run() {
    //
    // DESC
    // Main entry to the actual 'run' core of the mpmProg
    //

    int         ret     = 1;

    debug_push("run");

    // Sleep for 'msleepSeconds'
    mps_env->pcsm_stdout->colprintf("Sleeping for interval (seconds)", "[ %d ]\n",
                                     msleepSeconds);
    sleep(msleepSeconds);

    debug_pop();
    return ret;
}

//
//\\\***
// C_mpmProg_pathFind definitions ****>>>>
/////***
//

C_mpmProg_pathFind::C_mpmProg_pathFind(
    s_env*      aps_env,
    int		amvertex_start,
    int		amvertex_end) : C_mpmProg(aps_env)
{
    //
    // ARGS
    //	aps_env			in		parent env structure
    //	amvertex_start		in		default start vertex
    //	amvertex_end		in		default end vertex -- if set to
    //						-1 the end vertex is assumed to
    //						mean the total number of
    //						vertices.
    //
    // DESC
    // Basically a thin "fall-through" constructor to the base
    // class.
    //
    // PRECONDITIONS
    // o aps_env must be fully instantiated.
    //
    // HISTORY
    // 08 July 2010
    // o Initial design and coding.
    //

    debug_push("C_mpmProg_pathFind");

    mstr_obj	= "C_mpmProg_pathFind";
    vertexStart_set(amvertex_start);
    mvertex_end                 = 0;
    mb_surfaceRipClear          = false;

    s_env_activeSurfaceSetIndex(mps_env, 0);
    mvertex_total		= mps_env->pMS_primary->nvertices;

    if(	amvertex_start >= mvertex_total	||
        amvertex_start < 0 )
	vertexStart_set(0);
    else
	vertexStart_set(amvertex_start);

    if(	amvertex_end == -1 		|| 
        amvertex_end >= mvertex_total	|| 
        amvertex_end < 0)
	vertexEnd_set(mvertex_total-1);
    else
	vertexEnd_set(amvertex_end);

    debug_pop();
}

C_mpmProg_pathFind::~C_mpmProg_pathFind() {
    //
    // Destructor
    //

}

void
C_mpmProg_pathFind::print() {
  //
  // DESC
  // Simple info print method
  //

  C_mpmProg::print();
}

float
C_mpmProg_pathFind::cost_compute(
    int         a_start,
    int         a_end
) {
    //
    // ARGS
    // a_start          in              start vertex index
    // a_end            in              end vertex index
    //
    // DESC
    // Sets the main environment and then calls a dijkstra
    // computation.
    //

    int         ok;
    float       f_cost  = 0.0;
    static int  calls   = -1;
    char 	pch_buffer[65536];

    mps_env->startVertex = a_start;
    mps_env->endVertex   = a_end;

    calls++;
    ok = dijkstra(*mps_env);  // costs are written to vertex elements along
                              //+ the path in the MRIS structure.

    if(!ok) {
        fprintf(stderr, " fail ]\n");
        fprintf(stderr, "dijkstra failure, returning to system.\n");
        exit(1);
    }
    pULOUT(colsprintf(-50, 20, pch_buffer,
           "Marking (rip) path along vertices", "[ ok ]\n"));
//    s_env_activeSurfaceSetIndex(&st_env, (int) e_workingCurvature);
    f_cost      = surface_ripMark(*mps_env);
    if(mb_surfaceRipClear) {
        surface_ripClear(*mps_env, mb_surfaceRipClear);
    }
    return f_cost;
}

int
C_mpmProg_pathFind::run() {
    //
    // DESC
    // Main entry to the actual 'run' core of the pathFind module
    // -- for the most part, these innards comprise what used to be
    // the core of the original 'mris_pmake'.
    //

    float       f_cost  = 0.0;
    int         ret     = 1;
    string      str_optionsFQName       = "";
    string	str_patchFQName         = "";
    char 	pch_buffer[65536];

    debug_push("run");

    pSLOUT("PROCESSING: path\n");
    colsprintf(	mps_env->lw, mps_env->rw, pch_buffer, 
        	"Start->End vertices", "[ %d->%d ]\n",
        	mps_env->startVertex, mps_env->endVertex);

    // gcc11 complains about the indentation
    // not sure the author's intention. put {} around the first statement to be safe
    if(Gb_stdout) { printf("%s", pch_buffer); } pULOUT(pch_buffer);

    f_cost            = cost_compute(mvertex_start, mvertex_end);

    colsprintf(	mps_env->lw, mps_env->rw, pch_buffer, 
        	"Total path cost", " [ %f ]\n", f_cost);
    
    // gcc11 complains about the indentation
    // not sure the author's intention. put {} around the first statement to be safe
    if(Gb_stdout) { printf("%s", pch_buffer); } pULOUT(pch_buffer);
    pnRLOUT(lsprintf(mps_env->lw, pch_buffer, "%f", f_cost));

    if (mps_env->b_patchFile_save) {
        str_patchFQName =  mps_env->str_workingDir +
        		   mps_env->str_patchFileName;
        if (MRISwritePatch(mps_env->pMS_primary,
                               (char*) str_patchFQName.c_str()) != NO_ERROR)
        	exit(1);
	pULOUT(colsprintf(mps_env->lw, mps_env->rw, pch_buffer, 
        	"Saving patch file...", " [ ok ]\n"));
    }

    if (mps_env->b_labelFile_save) {
        //label_save(st_env);
        void* pv_void = NULL;
        label_workingSurface_saveTo(*mps_env, vertex_ripFlagIsTrue, pv_void);
        if (mps_env->b_surfacesKeepInSync) {
            surface_primaryToSecondary_ripTrueCopy(*mps_env);
            label_secondarySurface_saveTo(*mps_env, vertex_ripFlagIsTrue, pv_void);
        }
	pULOUT(colsprintf(mps_env->lw, mps_env->rw, pch_buffer, 
        	"Labeling and saving all target vertices...", " [ ok ]\n"));
    }

    if (mps_env->b_surfacesClear) {
        s_env_activeSurfaceSetIndex(mps_env, (int) e_workingCurvature);
        surface_ripClear(*mps_env, true);
        if (mps_env->b_surfacesKeepInSync) {
            s_env_activeSurfaceSetIndex(mps_env, (int) e_auxillary);
            surface_ripClear(*mps_env, true);
            // NB!! Remember to set the "active" surface back
            // to the "working" surface. The dijkstra()
            // function operates on this "active" surface.
            s_env_activeSurfaceSetIndex(mps_env, (int) e_workingCurvature);
        }
	pULOUT(colsprintf(mps_env->lw, mps_env->rw, pch_buffer, 
        	"Clearing (rip) path along vertices...", " [ ok ]\n"));
    }
    fflush(stdout);
    debug_pop();
    return ret;
}

//
//\\\***
// C_mpmProg_autodijk definitions ****>>>>
/////***
//

C_mpmProg_autodijk::C_mpmProg_autodijk(
    s_env*      aps_env) : C_mpmProg(aps_env)
{
    //
    // ARGS
    //
    // DESC
    // Basically a thin "fall-through" constructor to the base
    // class.
    //
    // PRECONDITIONS
    // o aps_env must be fully instantiated.
    //
    // HISTORY
    // 17 November 2009
    // o Initial design and coding.
    //

    debug_push("C_mpmProg_autodijk");

    mstr_obj	= "C_mpmProg_autodijk";
    mvertex_polar               = 0;
    mvertex_start               = 0;
    mvertex_step                = 1;
    mvertex_end                 = 0;
    m_costFunctionIndex         = 0;
    mb_performExhaustive        = false;
    mb_surfaceRipClear          = false;
    mb_worldMap                 = false;
    mb_simpleStatsShow          = true;
    mprogressIter               = 100;

    mpOverlayDistance           = NULL;
    mpOverlayOrig               = mps_env->pCmpmOverlay;

    s_env_activeSurfaceSetIndex(mps_env, 0);
    mstr_costFileName   = mps_env->str_costCurvFile;
    mstr_costFullPath   = mps_env->str_workingDir + "/" + mstr_costFileName;
    mvertex_end         = mps_env->pMS_primary->nvertices - 1;
    mvertex_total       = mps_env->pMS_primary->nvertices;
    mpf_cost            = new float[mvertex_total];
    mpf_persistent      = new float[mvertex_total];
    for(int i=0; i<mvertex_total; i++) {
        mpf_cost[i]             = 0.0;
        mpf_persistent[i]       = 0.0;
    }

    // By default, the mpf_fileSaveData points to the mpf_cost
    mpf_fileSaveData = mpf_cost;
    
    debug_pop();
}

C_mpmProg_autodijk::~C_mpmProg_autodijk() {
    //
    // Destructor
    //

    delete [] mpf_cost;
    delete [] mpf_persistent;
    if(mpOverlayDistance) delete mpOverlayDistance;
}

void
C_mpmProg_autodijk::worldMap_set(int avalue) {
    /*
     * This method controls the runtime behavior of the autodijkstra
     * mpmProg.
     * 
     */
    mb_worldMap         = avalue;
    if(mb_worldMap) {
        // Delete any existing distance overlays
        if(mpOverlayDistance) delete mpOverlayDistance;
        // and create a new distance overlay
        mpOverlayDistance = new C_mpmOverlay_distance(mps_env);
    }
}

bool 
C_mpmProg_autodijk::worldMap_shouldCreate() {
   return(mb_worldMap);
}


/*!
  \fn C_mpmProg_autodijk::CURV_fileWrite()
  \brief Write a FreeSurfer curvature array to a file
  \return If curvature file is successfully written, return e_OK, else return e_WRITEACCESSERROR.
*/
e_FILEACCESS
C_mpmProg_autodijk::CURV_fileWrite()
{
  FILE* FP_curv;
  int   i;
  char  pch_readMessage[65536];

  if((FP_curv = fopen(mstr_costFullPath.c_str(), "w")) == NULL)
    return(e_WRITEACCESSERROR);
  fwrite3(NEW_VERSION_MAGIC_NUMBER, FP_curv);
  fwriteInt(mvertex_total, FP_curv);
  fwriteInt(mps_env->pMS_primary->nfaces, FP_curv);
  fwriteInt((int)1, FP_curv);
  sprintf(pch_readMessage, "Writing %s", mstr_costFileName.c_str());
  for(i=0; i<mvertex_total; i++) {
    CURV_arrayProgress_print(mvertex_total, i, pch_readMessage);
    fwriteFloat(mpf_fileSaveData[i], FP_curv);
  }
  fclose(FP_curv);
  return(e_OK);
}


void
C_mpmProg_autodijk::print() {
  //
  // DESC
  // Simple info print method
  //

  C_mpmProg::print();
}

float
C_mpmProg_autodijk::cost_compute(
    int         a_start,
    int         a_end
) {
    //
    // ARGS
    // a_start          in              start vertex index
    // a_end            in              end vertex index
    //
    // DESC
    // Sets the main environment and then calls a dijkstra
    // computation.
    //
    // HISTORY
    // 19 January 2012
    // o Removed call to 'surface_ripMark()' -- autodijk's have no "paths".
    //   

    int         ok;
    float       f_cost  = 0.0;
    static int  calls   = -1;
    int         ret;
    
    mps_env->startVertex = a_start;
    mps_env->endVertex   = a_end;

    if(mb_worldMap && !mb_worldMapDistanceCalc) calls++;
    if(!mb_worldMap)                            calls++;
    if(!(calls % (mprogressIter))) {
        if(mb_worldMap && !mb_worldMapDistanceCalc) {
            mpf_fileSaveData = mpf_persistent;
            if((ret=CURV_fileWrite()) != e_OK)
                error("I could not save the cost curv file.");
            mps_env->pcsm_stdout->colprintf(
                            "here/end->opposite/total = cost",
                            "[ %6d/%6d -> %6d/%6d = ",
                            a_start,
                            mvertex_end,
                            ms_stats.indexMax,
                            mvertex_total);
        }
        if (!mb_worldMap) {
            mps_env->pcsm_stdout->colprintf(
                            "start->end = cost", "[ %6d -> %6d/%6d = 0.0]\n",
                            mps_env->startVertex,
                            mps_env->endVertex,
                            mvertex_total);
        }
    }
    ok = dijkstra(*mps_env);  // costs are written to vertex elements along
                              //+ the path in the MRIS structure.

    if(!ok) {
        mps_env->pcsm_stdout->printf(" fail ]\n");
        mps_env->pcsm_stdout->printf("dijkstra failure, returning to system.\n");
        exit(1);
    }
    f_cost = mps_env->pMS_active->vertices[ms_stats.indexMax].val;
    if(!(calls % (mprogressIter))) {
        if(mb_worldMap && !mb_worldMapDistanceCalc)
        mps_env->pcsm_stdout->printf(" %f ]\n", f_cost);
    }
    if(mb_surfaceRipClear) {
        // WARNING! This destroys the cost values and paths stored in the mesh!
        surface_ripClear(*mps_env, mb_surfaceRipClear);
    }
    return f_cost;
}

int 
C_mpmProg_autodijk::vertexCosts_pack(e_stats& a_stats) {
    /*
     * Pack the cost values stored in the FreeSurfer mesh into 
     * an array.
     */
    a_stats.f_max               =  0.0;
    a_stats.indexMax            = -1;
    for(int v = 0; v < mvertex_total; v++) {
        mpf_cost[v]     = mps_env->pMS_active->vertices[v].val;
        if(a_stats.f_max < mpf_cost[v]) {
            a_stats.f_max   = mpf_cost[v];
            a_stats.indexMax    = v;
        }
    }
    return true;
}

int
C_mpmProg_autodijk::run() {
    //
    // DESC
    // Main entry to the actual 'run' core of the mpmProg
    //

    float       f_cost                  =  0.0;
    float       f_valueAtFurthest       = 0.0;
    int         ret                     =  1;

    e_MPMPROG	        e_prog          = mps_env->empmProg_current;
    e_MPMOVERLAY        e_overlay       = mps_env->empmOverlay_current;

    debug_push("run");

    ms_stats.f_max        = 0.0;
    ms_stats.indexMax     = -1;

    mps_env->pcsm_stdout->colprintf("mpmProg (ID)", "[ %s (%d) ]\n",
                            mps_env->vstr_mpmProgName[e_prog].c_str(),
                            mps_env->empmProg_current);
    mps_env->pcsm_stdout->colprintf("mpmArgs", "[ %s ]\n",
                            mps_env->str_mpmArgs.c_str());
    mps_env->pcsm_stdout->colprintf("mpmOverlay bool", "[ %d ]\n",
                            mps_env->b_mpmOverlayUse);
    mps_env->pcsm_stdout->colprintf("mpmOverlay (ID)", "[ %s (%d) ]\n",
                            mps_env->vstr_mpmOverlayName[e_overlay].c_str(),
                            mps_env->empmOverlay_current);
    mps_env->pcsm_stdout->colprintf("mpmOverlayArgs", "[ %s ]\n",
                            mps_env->str_mpmOverlayArgs.c_str());

    if(mb_performExhaustive) {
      // Calculate the costs from polar to every other vertex in
      // mesh *explicitly*. DO NOT USE -- ONLY FOR DEBUGGING!! A
      // single sweep of the dijkstra from polar->polar will
      // have the same result in seconds... the exhaustive search
      // can take multiple hours.
      for(int v = mvertex_start; v <= mvertex_end; v+=mvertex_step) {
          f_cost        = cost_compute(mvertex_polar, v);
          mpf_cost[v]   = f_cost;
      }
    } else if(mb_worldMap) {
        mps_env->pcsm_stdout->lprintf("Computing World Map...\n");
        mps_env->pcsm_stdout->colprintf(
                              "Progress iteration interval",
                              "[ %d ]\n", mprogressIter);
        for(int v = mvertex_start; v <= mvertex_end; v+=mvertex_step) {
            // First we need to find the anti-pole for current 'v'
            // by setting the env overlay to the distance object, 
            // performing a single sweep from 'v' to 'v' and then
            // finding the highest "cost" which corresponds to the vertex at
            // furthest distance from 'v' on the mesh.
            mps_env->pCmpmOverlay = mpOverlayDistance;
            mb_worldMapDistanceCalc = true;
            f_cost      = cost_compute(v, v);
            vertexCosts_pack(ms_stats);
            // Now set the env overlay back to the original, and recompute
            mps_env->pCmpmOverlay = mpOverlayOrig;
            mb_worldMapDistanceCalc = false;
            f_cost      = cost_compute(v, v);
            // For this run, we only store the single cost value from the 
            // maxIndex vertex from the distance overlay
            f_valueAtFurthest = mps_env->pMS_active->vertices[ms_stats.indexMax].val;
            mpf_persistent[v] = f_valueAtFurthest;
        }
        // Now, set the save pointer to the correct data to save
        mpf_fileSaveData = mpf_persistent;
    } else {
      f_cost            = cost_compute(mvertex_polar, mvertex_polar);
      vertexCosts_pack(ms_stats);
    }
    if(mb_simpleStatsShow) {
        mps_env->pcsm_stdout->colprintf("max(cost) @ index",
                                        "[ %f : %d ]\n",
                                        ms_stats.f_max,
                                        ms_stats.indexMax);
    }
    // Write the cost curv file
    // NOTE: This saves the data pointed to by mpf_fileSaveData!
    if((ret=CURV_fileWrite()) != e_OK)
        error("I could not save the cost curv file.");

    debug_pop();
    return ret;
}

//
//\\\***
// C_mpmProg_autodijk_fast definitions ****>>>>
/////***
//

///
/// Constructor
/// \param apps_env Freesurfer environment
///
C_mpmProg_autodijk_fast::C_mpmProg_autodijk_fast(s_env* aps_env) :
    C_mpmProg_autodijk(aps_env)
{
    debug_push("C_mpmProg_autodijk_fast");
    mstr_obj	= "C_mpmProg_autodijk_fast";

    debug_pop();
}

///
/// Destructor
///
C_mpmProg_autodijk_fast::~C_mpmProg_autodijk_fast()
{
}

///
/// Convert the freesurfer environment MRI surface representation into
/// a representation that is formatted for consumption by the OpenCL
/// algorithm
/// \param graph Pointer to the graph data to generate
///
void C_mpmProg_autodijk_fast::genOpenCLGraphRepresentation(GraphData *graph)
{
    MRIS*           surf = mps_env->pMS_active;
//    s_iterInfo      st_iterInfo;
//    bool            b_relNextReference  = true;

    cout << "Converting graph to fast representation... " << endl;
    // Allocate memory for each of the vertices
    graph->vertexCount = surf->nvertices;
    graph->vertexArray = new int[surf->nvertices];

    // Determine the size of the edge and weight buffers
    int numEdges = 0;
    for(int i = 0; i < surf->nvertices; i++)
    {
        numEdges += surf->vertices_topology[i].vnum;
    }

    // Allocate memory for the edge and weight arrays
    graph->edgeCount = numEdges;
    graph->edgeArray = new int[numEdges];
    graph->weightArray = new float[numEdges];

    // Now create the edge list and compute all of the edge weights
    int curEdgeIndex = 0;
    for(int i = 0; i < surf->nvertices; i++)
    {
        VERTEX_TOPOLOGY const * const curVertex = &surf->vertices_topology[i];

        // Assign the edge index for this vertex
        graph->vertexArray[i] = curEdgeIndex;

        for(int j = 0; j < curVertex->vnum; j++)
        {
	    int ij;
            // index for the neighbor
            graph->edgeArray[curEdgeIndex] = curVertex->v[j];
	    ij = curVertex->v[j];

            // Compute the weight for this edge
//            float cost = mps_env->costFunc_do(*mps_env, &st_iterInfo, i, j, b_relNextReference);
	    float cost = s_env_edgeCostFind(*mps_env, i, ij);	    
            graph->weightArray[curEdgeIndex] = cost;

            curEdgeIndex++;
        }
    }
}
///
/// Run method - the main entry to the 'run' core of the mpmProg
///
int C_mpmProg_autodijk_fast::run()
{
    int         ret     = 1;

    debug_push("run");

    GraphData graph;
    genOpenCLGraphRepresentation(&graph);

    // Create the list of source indices, here we only ever do one source
    int sourceVertices = mvertex_polar;

    // Allocate array for results (the algorithm will write out a result
    // for each vertex in the graph)
    float *results = new float[graph.vertexCount];

    // Perform a Dijkstra run
    cout << "Running Dijkstra's algorithm..." << endl;

#ifdef FS_OPENCL
    // We have implemented various flavors of the Dijkstra algorithm that runs with
    // OpenCL using either CPU-only, GPU+CPU, Multi GPU, or Multi GPU + CPU.  Based
    // on what set of devices is available, this function will choose which implementation
    // to use.  This version of Dijkstra selects which devices to use automatically.

    // If compiled with OpenCL support, run the OpenCL version of the algorithm
    runDijkstraOpenCL(&graph, &sourceVertices, results, 1);
#else
    // If not compiled with OpenCL, run the reference version of the algorithm
    runDijkstraRef(&graph, &sourceVertices, results, 1);
#endif

    cout << "Done." << endl;

    // Save back the resulting costs
    for(int v = mvertex_start; v <= mvertex_end; v+=mvertex_step)
    {
        mpf_cost[v] = results[v];
    }

    // Free temporary results buffer
    delete [] results;

    // Free up memory from allocation of graph data
    delete [] graph.vertexArray;
    delete [] graph.edgeArray;
    delete [] graph.weightArray;

    // Write the cost curv file
    if((ret=CURV_fileWrite()) != e_OK)
        error("I could not save the cost curv file.");

    debug_pop();
    return ret;
}

//
//\\\***
// C_mpmProg_ROI definitions ****>>>>
/////***
//

C_mpmProg_ROI::C_mpmProg_ROI(
    s_env*      aps_env, 
    string      astr_vertexFile,
    float 	af_radius) : C_mpmProg(aps_env)
{
    //
    // ARGS
    //
    // DESC
    // Basically a thin "fall-through" constructor to the base
    // class.
    //
    // PRECONDITIONS
    // o aps_env must be fully instantiated.
    //
    // HISTORY
    // 05 May 2011
    // o Initial design and coding.
    //

    debug_push("C_mpmProg_ROI");
    mstr_obj	        = "C_mpmProg_ROI";

    if(astr_vertexFile.length())
        vertexFile_load(astr_vertexFile);

    mf_radius           = af_radius;
    mf_plyIncrement     = 1.0;
    mb_surfaceRipClear  = true;
    m_borderSize	= 0;

    debug_pop();
}

int
C_mpmProg_ROI::labelFile_load(
        string                  astr_fileName
) {
    //
    // ARGS
    // astr_fileName            string          label file to read
    //
    // DESC
    //  This method reads the contents of the <astr_fileName>, ripping
    //  each vertex index in the internal mesh -- in so doing defining
    //  a 'label' on the mesh.
    //
    // Loading a label file typically implies that *all* the labelled
    // vertices together create a single ROI, and a single ROI generated
    // output file is created (cf vertexFile_load()).
    //
    // PRECONDITIONS
    //  o <astr_fileName> is a FreeSurfer Label file.
    //
    // POSTCONDITIONS
    //  o For each vertex index in the label file, the mesh vertex
    //    has its ripflag set to TRUE.
    //  o Vertex indices are also stored in the internal vertex vector.
    //

    MRIS*       pmesh           = NULL;
    LABEL*      pLBL            = NULL;
    string      lstr_ext        = ".";
    int         i               = 0;
    char        ch_mark         = TRUE;
    char*       pch_mark        = &ch_mark;
    void*       pv_mark         = (void*) pch_mark;

    mstr_labelFile              = astr_fileName;
    lstr_ext                    += fio_extension(mstr_labelFile.c_str());
    mstr_outputStem             = fio_basename(mstr_labelFile.c_str(),
                                  lstr_ext.c_str());
    pmesh                       = mps_env->pMS_primary;
    // Mark the mesh
    label_coreLoad(pmesh, astr_fileName, vertex_ripFlagMark, pv_mark);

    // and store the "rip"ed vertices -- this actually opens
    // the label file again...
    if(mv_vertex.size()) {
        mv_vertex.clear();
    }
    pLBL = LabelRead((char*)"", (char*) astr_fileName.c_str());
    for (i = 0; i < pLBL->n_points; i++) {
      mv_vertex.push_back(pLBL->lv[i].vno);
    }
    LabelFree(&pLBL);

    mb_ROIsInSeparateLabels = false;
    return true;
}

int
C_mpmProg_ROI::label_savePly(
        string                  astr_filePrefix,
        bool                    ab_staggered,
        float                   af_plyIncrement
) {
    mf_plyIncrement     = af_plyIncrement;
    mb_saveStaggered    = ab_staggered;
    label_ply_save(*mps_env, astr_filePrefix, mb_saveStaggered);
    return true;
}

int
C_mpmProg_ROI::labelFile_save(
        string                  astr_fileName
) {

    MRIS*       pmesh           = NULL;

    pmesh       = mps_env->pMS_primary;
    label_coreSave(pmesh, astr_fileName, vertex_ripFlagIsTrue,
            (void*) (char)TRUE);
    return true;
}

int
C_mpmProg_ROI::border_mark(void)
{
/*
 * ARGS
 * 	(void)
 *
 * DESC
 * Redefines the pattern of RIPFLAGS on a mesh so that contiguous
 * regions only have outer border vertices ripped.
 *
 * PRECONDITIONS
 * o An existing pattern of rips must already be present on the surface.
 * o The "primary" surface mesh is processed.
 *
 * POSTCONDITIONS
 * o Contiguous regions will have border vertices marked.
 * o Number of border vertices is returned.
 *
 */
    unsigned int	i		= 0;
    int         	vertex          = -1;
    int			neighbor	= -1;
    int			neighborCount	= -1;
    int			markedCount	= 0;
    MRIS*       	mesh            = NULL;
    bool		b_innerVertex 	= true;

    mesh = mps_env->pMS_primary;

    // Start by clearing the mv_vertex array. This array will be re-purposed
    // to contain the indices of the border vertices.
    if(mv_vertex.size()) {
        mv_vertex.clear();
    }

    m_borderSize = 0;
    for(vertex = 0; vertex < mesh->nvertices; vertex++) {
	if(mesh->vertices[vertex].ripflag) {
	    markedCount++;
	    b_innerVertex = true;
	    VERTEX_TOPOLOGY const * const SVertex = &mesh->vertices_topology[vertex];
	    for(neighborCount = 0; neighborCount < SVertex->vnum; neighborCount++) {
		neighbor = SVertex->v[neighborCount];
		b_innerVertex &= mesh->vertices[neighbor].ripflag;
	    }
	    if(!b_innerVertex) {
		m_borderSize++;
		mv_vertex.push_back(vertex);
	    }
	}
    }

    if(m_borderSize) {
        // At this point, the mv_vertex array contains all the border vertices.
        // Now, clear the current ripflag pattern, and then assign rips
        // according the mv_vertex array
        surface_ripClear(*mps_env, true);
        for(i=0; i<mv_vertex.size(); i++) {
            vertex = mv_vertex[i];
            mesh->vertices[vertex].ripflag = true;
        }
    }
    return m_borderSize;
}

int
C_mpmProg_ROI::vertexFile_load(string astr_fileName)
{
    //
    // ARGS
    // astr_fileName            string          file to read
    //
    // DESC
    //  This method reads the contents of the <astr_fileName> into the
    //  internal list vertex vector list.
    //
    //  The text filename a simple line-delimited list of vertex indices:
    //
    //          v1
    //          v2
    //          ...
    //          vn
    //
    // Passing a vertex file typically implies that each vertex
    // is the seed of a single ROI that should be saved to its
    // own output file.
    //
    // PRECONDITIONS
    //  o <astr_fileName> must exist and MUST be Nx1
    //
    // POSTCONDITIONS
    //  o The internal mv_vertex is initialized.
    //  o Returns the actual number of indices read from file.
    //  o Vertex ripflag is NOT set.
    //

    int         vertex          = -1;
    int         readCount       = 0;
    string      lstr_ext        = ".";
    MRIS*       mesh            = NULL;


    mesh = mps_env->pMS_primary;

    // First turn off any "rip" marks on the surface
    for(vertex = 0; vertex < mesh->nvertices; vertex++)
        mesh->vertices[vertex].ripflag = FALSE;

    if(!astr_fileName.length())
        astr_fileName           = mstr_vertexFile;
    else
        mstr_vertexFile         = astr_fileName;

    lstr_ext                    += fio_extension(mstr_vertexFile.c_str());
    mstr_outputStem             = fio_basename(mstr_vertexFile.c_str(),
                                  lstr_ext.c_str());

    ifstream    ifs_vertex(astr_fileName.c_str());
    if(!ifs_vertex) return false;

    if(mv_vertex.size()) {
        mv_vertex.clear();
    }

    while(!ifs_vertex.eof()) {
        readCount++;
        ifs_vertex >> vertex;
        mv_vertex.push_back(vertex);
    }

    mb_ROIsInSeparateLabels = true;

    return readCount;
}

C_mpmProg_ROI::~C_mpmProg_ROI() {
    //
    // Destructor
    //

}

int
C_mpmProg_ROI::run() {
    //
    // DESC
    // Main entry to the actual 'run' core of the mpmProg
    //

    int          ret                    = 1;
    bool         b_origHistoryFlag      = mps_env->b_costHistoryPreserve;
    bool         b_surfaceCostVoid      = false;
    unsigned int i			= 0;
    int          j                      = 0;
    mps_env->b_costHistoryPreserve      = true;
    MRIS*        pmesh                  = mps_env->pMS_active;

    debug_push("run");

    if(mb_boundaryOnly) border_mark();

    if(!mb_ROIsInSeparateLabels) {
        for (i=0; i<(unsigned int)pmesh->nvertices; i++) {
            if (pmesh->vertices[i].ripflag == TRUE) {
                b_surfaceCostVoid       = !j++;
                mps_env->startVertex    = i;
                mps_env->endVertex      = i;
                ret = dijkstra(*mps_env, mf_radius, b_surfaceCostVoid);
            }
        }
        mps_env->b_costHistoryPreserve = b_origHistoryFlag;
        Gsout.str(std::string());
        Gsout << "-r" << mf_radius;
        label_savePly(mstr_outputStem + Gsout.str(),
                      mb_saveStaggered, mf_plyIncrement);
        // For the "single" label file output, the surface clearing
        // is dependent on the mb_surfaceRipClear flag.
        if(mb_surfaceRipClear) {
            // WARNING! This destroys the cost values and paths
            // stored in the mesh!
            surface_ripClear(*mps_env, mb_surfaceRipClear);
        }
    } else {
        for(i=0; i<mv_vertex.size(); i++) {
            mps_env->startVertex        = mv_vertex[i];
            mps_env->endVertex          = mv_vertex[i];
            ret = dijkstra(*mps_env, mf_radius, b_surfaceCostVoid);
            Gsout.str(std::string());
            Gsout << "-v" << mv_vertex[i] << "-r" << mf_radius << ".label";
            label_savePly(mstr_outputStem + Gsout.str(),
                          mb_saveStaggered, mf_plyIncrement);
            // In the case of separate label files, clear the surface
            // of rips after each save.
            surface_ripClear(*mps_env, true);
        }
    }

    debug_pop();
    return ret;
}

//
//\\\***
// C_mpmProg_externalMesh definitions ****>>>>
/////***
//

C_mpmProg_externalMesh::C_mpmProg_externalMesh(
    s_env*      aps_env,
    string      astr_meshFile) : C_mpmProg(aps_env)
{
    //
    // ARGS
    //
    // DESC
    // Basically a thin "fall-through" constructor to the base
    // class.
    //
    // PRECONDITIONS
    // o aps_env must be fully instantiated.
    //
    // HISTORY
    // 05 May 2011
    // o Initial design and coding.
    //

    debug_push("C_mpmProg_externalMesh");
    mstr_obj	        = "C_mpmProg_externalMesh";

    mb_surfaceRipClear  = true;

    debug_pop();
}

C_mpmProg_externalMesh::~C_mpmProg_externalMesh() {
    //
    // Destructor
    //

}

int
C_mpmProg_externalMesh::run() {
    //
    // DESC
    // Main entry to the actual 'run' core of the mpmProg
    //

    int         ret                     = 1;
    debug_push("run");

    debug_pop();
    return ret;
}


/* eof */
