/**
 * @file  C_mpmOverlay.cpp
 * @brief The internal 'program' API
 *
 *  The C_mpmOverlay class (and derivatives) presents an encapsulated
 *  interface to 'overlays' which are the 'curv' field component of
 *  of FreeSurfer meshes.
 * 
 *  The class also maintains cost function evaluation and weight vector
 *  structures.
 * 
 *  Since the class contains a pointed to the mris_pmake 'env' structure,
 *  its name is prefixed by 'mpm'.
 */
/*
 * Original Author: Rudolph Pienaar
 * CVS Revision Info:
 *    $Author: rudolph $
 *    $Date: 2011/05/31 18:18:49 $
 *    $Revision: 1.3 $
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


#include "C_mpmOverlay.h"
#include "dijkstra.h"
#include "general.h"

#include "c_surface.h"
#include "c_vertex.h"
#include "unistd.h"

#include <sstream>
#include <iostream>
#include <fstream>

extern 	bool     	Gb_stdout;
extern  stringstream 	Gsout;

//
//\\\---
// Base class issues --->>>
/////---
//

//-----------------------------------------------------------
// Base socket constructors/destructor, output dump routines
//-----------------------------------------------------------

void
C_mpmOverlay::debug_push(
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

void C_mpmOverlay::debug_pop() {
    mstackDepth--;
}

void
C_mpmOverlay::error(
    string      astr_msg        /* = ""  Error message */,
    int         code            /* = 1  Error ID  */)
{

    extern int errno;

    cerr << "\nFatal error encountered.\n";
    cerr << "\tC_mpmOverlay `" << mstr_name << "' (id: " << mid << ")\n";
    cerr << "\tCurrent function: " << mstr_obj << "::" << str_proc_get() << "\n";
    cerr << "\t" << astr_msg << "\n";
    perror("Error returned from system");
    print();
    cerr << "Throwing exception to (this) with code " << code << "\n\n";
    throw(this);
}

void
C_mpmOverlay::warn(
    string      astr_msg        /* = ""  Warning message */,
    int         code            /* = 1  Warning code  */)
{
    if (mwarnings) {
        cerr << "\nWarning.\n";
        cerr << "\tC_mpmOverlay `" << mstr_name << "' (id: " << mid << ")\n";
        cerr << "\tCurrent function: " << mstr_obj << "::" << str_proc_get() << "\n";
        cerr << "\t" << astr_msg << " (warning code: " << code << ")\n";
    }
}

void
C_mpmOverlay::function_trace(
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
C_mpmOverlay::core_construct(
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
    mstr_obj = "C_mpmOverlay";
}

C_mpmOverlay::C_mpmOverlay(
    s_env*      aps_env
){
    string      str_name        =       "Unnamed C_mpmOverlay";
    int         id              = 0;
    core_construct(str_name, id);
    mps_env                     = aps_env;
}

C_mpmOverlay::C_mpmOverlay(const C_mpmOverlay &C_mpmOverlay) {
    //
    // Copy constructor
    //

    debug_push("C_mpmOverlay (copy constructor)");
    error("Copy constructor not yet defined");
    debug_pop();
}

C_mpmOverlay & C_mpmOverlay::operator=(const C_mpmOverlay & C_mpmOverlay) {
    //
    // Overloaded (=) operator
    //

    debug_push("operator=");
    error("Overloaded operator= not yet defined");
    return *this;       // Never executed --
    debug_pop();
}

C_mpmOverlay::~C_mpmOverlay(void) {}

void
C_mpmOverlay::print() {
    cout << "object name:\t"    << mstr_name    << endl;
    cout << "object id:\t"      << mid          << endl;
    cout << "object type:\t"    << mstr_obj     << endl;
}

string
C_mpmOverlay::curvFileName_get(
    EOVERLAY		ae_overlay
    ) 
    //
    // ARGS
    // ae_overlay		enum 		overlay enumeration
    //
    // DESC
    //  This method returns the filename for a given overlay, based 
    //  on current hemisphere and env surface choice.
    //
{
    string	str_hemi	= mps_env->str_hemi;
    string	str_prefix("");
    string	str_suffix("");
    string	str_curvFileName("");

    str_prefix		= mstr_curvPrefixTemplate[ae_overlay];
    str_findAndReplace(str_prefix, "HEMI.", str_hemi);
    str_findAndReplace(str_prefix, "SURF.", mstr_surface);
    str_suffix		= mstr_curvSuffix[ae_overlay];
    str_curvFileName	= str_prefix + str_suffix;

    return str_curvFileName;
}

bool
C_mpmOverlay::costVector_read(string astr_fileName)
{
    //
    // ARGS
    // astr_fileName		string 		file to read
    //
    // DESC
    // 	This method reads the contents of the <astr_fileName> into the
    //	vectors for cost function weights and delta weights.
    //
    //	The filename contains:
    //	w1	dw1
    //	w2	dw2
    //	...	...
    //	wn	dwn
    //
    // PRECONDITIONS
    // 	o <astr_fileName> must exist and MUST be Nx2
    //
    // POSTCONDITIONS
    // 	o The internal mpv_costWeight and mpv_costWeightDel are initialized.
    //

    float 	f_valcost	= 0.;
    float 	f_valcostdel	= 0.;

    if(!astr_fileName.length())
	astr_fileName		= mstr_costWeightFile;
    else
	mstr_costWeightFile	= astr_fileName;

    ifstream	ifs_costVector(astr_fileName.c_str());
    if(!ifs_costVector) return false;

    if(mv_costWeight.size()) {
	mv_costWeight.clear();
	mv_costWeight.clear();
    }	    

    while(!ifs_costVector.eof()) {
	ifs_costVector >> f_valcost;
	ifs_costVector >> f_valcostdel;
	mv_costWeight.push_back(f_valcost);
	mv_costWeightDel.push_back(f_valcostdel);
    }
    return true;
}

bool
C_mpmOverlay::costVector_write(string astr_fileName)
{
    //
    // ARGS
    // astr_fileName		string 		file to write
    //
    // DESC
    // 	This method writes the cost weights and delta weights to
    //  <astr_fileName>.
    //
    //	The filename contains:
    //	w1	dw1
    //	w2	dw2
    //	...	...
    //	wn	dwn
    //
    // PRECONDITIONS
    // 	o The internal mpv_costWeight and mpv_costWeightDel must be
    //    initialized.
    //
    // POSTCONDITIONS
    // 	o <astr_fileName> contains the cost vector and delta weights.
    //

    int	row	= 0;

    if(!astr_fileName.length())
	astr_fileName		= mstr_costWeightFile;
    else
	mstr_costWeightFile	= astr_fileName;

    ofstream	ofs_costVector(astr_fileName.c_str());
    if(!ofs_costVector) return false;

    for(row = 0; row < (int) mv_costWeight.size(); row++) {
	ofs_costVector << mv_costWeight[row] 		<< " ";
	ofs_costVector << mv_costWeightDel[row] 	<< endl;
    }

    ofs_costVector.close();
    return true;
}


//
//\\\***
// C_mpmOverlay_NOP definitions ****>>>>
/////***
//

void
C_mpmOverlay_NOP::costWeightVector_init(void) {
    //
    // ARGS
    //
    // DESC
    // Initialize the cost weight vector for this class
    //
    // POSTCONDITIONS
    // 	o Both weight vectors are initialized to '1'.
    //
    // HISTORY
    // Late June 2010
    // o Initial design and coding.
    //

    // Use the 'eoverlay' ord value to determine size of full 
    // cost vector
    mv_costWeight.clear();
    mv_costWeightDel.clear();
    mv_costWeight.push_back(0.);
    mv_costWeightDel.push_back(0.);
}

C_mpmOverlay_NOP::C_mpmOverlay_NOP(
    s_env*      aps_env) : C_mpmOverlay(aps_env)
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

    debug_push("C_mpmOverlay_NOP");
    mstr_obj	= "C_mpmProg_NOP";
    
    mstr_costWeightFile	= "M_weights_NOP.mat";

    if(!costVector_read()) {
	costWeightVector_init();
	costVector_write();
    }

    mb_created	= true;
    debug_pop();
}

C_mpmOverlay_NOP::~C_mpmOverlay_NOP() {
    //
    // Destructor
    //

}

float
C_mpmOverlay_NOP::costEdge_calc(int i, int j) {

	return 0.0;
}

//
//\\\***
// C_mpmOverlay_unity definitions ****>>>>
/////***
//

void
C_mpmOverlay_unity::costWeightVector_init(void) {
    //
    // ARGS
    //
    // DESC
    // Initialize the cost weight vector for this class
    //
    // POSTCONDITIONS
    // 	o Both weight vectors are initialized to '1'.
    //
    // HISTORY
    // Late June 2010
    // o Initial design and coding.
    //

    // Use the 'eoverlay' ord value to determine size of full 
    // cost vector
    mv_costWeight.clear();
    mv_costWeightDel.clear();
    mv_costWeight.push_back(1.);
    mv_costWeightDel.push_back(1.);
}

C_mpmOverlay_unity::C_mpmOverlay_unity(
    s_env*      aps_env) : C_mpmOverlay(aps_env)
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

    debug_push("C_mpmOverlay_unity");
    mstr_obj	= "C_mpmProg_unity";
    mstr_costWeightFile	= "M_weights_unity.mat";

    if(!costVector_read()) {
	costWeightVector_init();
	costVector_write();
    }

    mb_created	= true;
    debug_pop();
}

C_mpmOverlay_unity::~C_mpmOverlay_unity() {
    //
    // Destructor
    //

}

float
C_mpmOverlay_unity::costEdge_calc(int i, int j) {
    return 1.0 * mv_costWeight[0];
}

//
//\\\***
// C_mpmOverlay_FScurvs definitions ****>>>>
/////***
//

C_mpmOverlay_FScurvs::C_mpmOverlay_FScurvs(
    s_env*      aps_env,
    string	astr_costWeightFile) : C_mpmOverlay(aps_env)
{
    //
    // ARGS
    //
    // DESC
    // For each EOVERLAY "curvature", read the corresponding FreeSurfer
    // file into local structures.
    //
    // PRECONDITIONS
    // o aps_env must be fully instantiated.
    //
    // HISTORY
    // Late June 2010
    // o Initial design and coding.
    //

    debug_push("C_mpmOverlay_FScurvs");
    mstr_obj	= "C_mpmProg_FScurvs";
    
    if(astr_costWeightFile.length())
	mstr_costWeightFile	= astr_costWeightFile;
    else
	mstr_costWeightFile	= "M_weights_FScurvs.mat";

    if(!costVector_read()) {
	costWeightVector_init();
	costVector_write();
    }

    // Size of overlay arrays
    mv_size	= aps_env->pMS_curvature->nvertices;
    mb_created	= true;

    debug_pop();
}

C_mpmOverlay_FScurvs::~C_mpmOverlay_FScurvs() {
    //
    // Destructor
    //

}

void
C_mpmOverlay_FScurvs::costWeightVector_init(void) {
    //
    // ARGS
    //
    // DESC
    // Initialize the cost weight vector for this class
    //
    // POSTCONDITIONS
    // 	o Both weight vectors are initialized to '1'.
    //
    // HISTORY
    // Late June 2010
    // o Initial design and coding.
    //

    // Use the 'eoverlay' ord value to determine size of full 
    // cost vector
    mv_costWeight.clear();
    mv_costWeightDel.clear();
    for(int i=0; i<(int)eoverlay; i++) {
	mv_costWeight.push_back(1.);
	mv_costWeightDel.push_back(1.);
    }
}

float
C_mpmOverlay_FScurvs::costEdge_calc(int i, int j) {

	return 0.0;
}

//
//\\\***
// C_   mpmOverlay_distance definitions ****>>>>
/////***
//

void
C_mpmOverlay_distance::costWeightVector_init(void) {
    //
    // ARGS
    //
    // DESC
    // Initialize the cost weight vector for this class
    //
    // POSTCONDITIONS
    // 	o Both weight vectors are initialized to '1'.
    //
    // HISTORY
    // April 2011
    // o Initial design and coding.
    //

    mv_costWeight.clear();
    mv_costWeightDel.clear();
    mv_costWeight.push_back(1.);
    mv_costWeightDel.push_back(1.);
}

C_mpmOverlay_distance::C_mpmOverlay_distance(
    s_env*      aps_env) : C_mpmOverlay(aps_env)
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

    debug_push("C_mpmOverlay_distance");
    mstr_obj	= "C_mpmOverlay_distance";
    mstr_costWeightFile	= "M_weights_distance.mat";

    if(!costVector_read()) {
	costWeightVector_init();
	costVector_write();
    }

    mb_created	= true;
    debug_pop();
}

C_mpmOverlay_distance::~C_mpmOverlay_distance() {
    //
    // Destructor
    //

}

float
C_mpmOverlay_distance::costEdge_calc(int i, int j) {
    //
    // Return the cost in moving from vertex 'i' to vertex 'j'.
    // In this overlay, this is simply the edge length between
    // the two vertices, as defined in the FreeSurfer mesh
    // structure.
    //

    VERTEX*     pVrtx_i         = &mps_env->pMS_active->vertices[i];
    float       wd              = mv_costWeight[0];
    float       f_distance      = 0.;
    float       f_cost          = 0.;
    int 	jrel            = 0;
    int 	jrelcount;

    debug_push ("costEdge_calc (...)");
    
    for(jrelcount=0; jrelcount< pVrtx_i->vnum; jrelcount++) {
        if(pVrtx_i->v[jrelcount] == j) {
            jrel = jrelcount;
            break;
        }
    }
    j = jrel;
    f_distance  = pVrtx_i->dist[j];

    f_cost  = f_distance * wd;
    if (wd <= 0.) error("negative edge length detected!", 1);

    debug_pop();
    return(f_cost);
}

/* eof */
