/***************************************************************************
 *   Copyright (C) 2009 by Rudolph Pienaar                                 *
 *   Childrens Hospital Boston                                             *
 *   rudolph.pienaar@childrens.harvard.edu                                 *
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
// $Id: C_mpmProg.cpp,v 1.2 2009/11/20 22:05:10 rudolph Exp $

#include "C_mpmProg.h"
#include "dijkstra.h"

#include "c_surface.h"
#include "c_vertex.h"

#include <sstream>

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

    extern int errno;

    cerr << "\nFatal error encountered.\n";
    cerr << "\tSSocket `" << mstr_name << "' (id: " << mid << ")\n";
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
        cerr << "\tSSocket `" << mstr_name << "' (id: " << mid << ")\n";
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
        if (str_objectName != mstr_name)
        cerr << "\nSSocket `" << str_name_get() << "' (id: " << mid << ")" << endl;
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
    mstr_obj = "C_SSocket";
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
    // HISTORY
    // 17 November 2009
    // o Initial design and coding.
    //

    debug_push("C_mpmProg_autodijk");

    mvertex_polar       = 0;
    mvertex_start       = 0;
    mvertex_step        = 1;
    mvertex_end         = 0;
    m_costFunctionIndex = 0;
    mb_surfaceRipClear  = true;

    if(s_env_costFctSetIndex(mps_env, m_costFunctionIndex) == -1)
        error_exit("setting costFunctionIndex", "Could not set index", 1);
    s_env_activeSurfaceSetIndex(mps_env, 0);
    mvertex_end         = mps_env->pMS_curvature->nvertices;
    
    debug_pop();
}

C_mpmProg_autodijk::~C_mpmProg_autodijk() {
  //
  // Destructor
  //
}

void
C_mpmProg_autodijk::print() {
  //
  // DESC
  // Simple info print method
  //

  C_mpmProg::print();
}

int
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

    int         ret;

    mps_env->startVertex = a_start;
    mps_env->endVertex   = a_end;

    ret = dijkstra(*mps_env); // costs are written to vertex elements along
                              //+ the path in the MRIS structure.

    if(mb_surfaceRipClear) {
        surface_ripClear(*mps_env, mb_surfaceRipClear);
    }                    
    return ret;
}

int
C_mpmProg_autodijk::run() {
    //
    // DESC
    // Main entry to the actual 'run' core of the mpmProg
    //

    int         ret = 0;

    for(int v = mvertex_start; v < mvertex_end; v+=mvertex_step) {
        ret = cost_compute(mvertex_polar, v);
    }
    return ret;
}


/* eof */
