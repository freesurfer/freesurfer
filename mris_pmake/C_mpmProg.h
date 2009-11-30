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
/// \file C_mpmProg.h
///
/// \brief Brief description
/// The internal 'program' API.
///
/// \b DESCRIPTION
/// C_mpmProgs are overloaded classes that perform specific functions in the
/// context of the Dijkstra system. The default contol system allows an
/// external network-based program to drive the core dijkstra search. This
/// process can be very slow on large-scale problems. To alleviate that, 
/// mpm_programs can be written. These map directly to external dsh scripts
/// but without any of the network overhead.
///
/// \b HISTORY
/// 16 November 2009 - Initial consolidation from several other sources.
/// $Id: C_mpmProg.h,v 1.4 2009/11/30 20:14:44 rudolph Exp $
///
///

#ifndef __C_MPM_PROG_H__
#define __C_MPM_PROG_H__

#ifdef __cplusplus
extern  "C" {
#endif

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "error.h"
#include "fio.h"

#ifdef __cplusplus
}
#endif

#include "env.h"

#include <string>
using namespace std;

const int       MPMSTACKDEPTH     = 64;

class C_mpmProg {

    //
    // Data members
    //

  public:

    // type and var info
    string      mstr_obj;                       // name of object class
    string      mstr_name;                      // name of object variable
    int         mid;                            // id of socket
    int         mverbosity;                     // Debug related value 
    int         mwarnings;                      // Show warnings
    int         mstackDepth;                    // Current procedure stackDepth
    string      mstr_proc[MPMSTACKDEPTH];       // Used to track the current
                                                //+ procedure being executed
  protected:

    // base class info

    s_env*      mps_env;                        // Pointer to the main
                                                //+ environment
    bool        b_created;                                                                   

    //
    // Method members
    //

  public:

    //
    // Constructor / destructor block
    //
    void core_construct(
        string          astr_name       = "unnamed",
        int             a_id            = -1,
        int             a_verbosity     = 0,
        int             a_warnings      = 0,
        int             a_stackDepth    = 0,
        string          astr_proc       = "noproc"
    );
    C_mpmProg(s_env* aps_env);
    virtual ~C_mpmProg(void);
    C_mpmProg(
        const C_mpmProg & C_mpmProg
    );
    C_mpmProg & operator=(const C_mpmProg & C_mpmProg);

    //
    // Error / warn /  print block
    //
    void        debug_push(     string  astr_currentProc);
    void        debug_pop();
    void        error(          string  astr_msg        = "",
                                int     code            = 1);
    void        warn(           string  astr_msg        = "",
                                int     code            = 1);
    void        function_trace( string  astr_msg        = "",
                                string  astr_separator  = "");
    void        print();

    //
    // Access "housekeeping" state info
    //

    const string        str_obj_get()           const {
        return mstr_obj;
    };
    const string        str_name_get()          const {
        return mstr_name;
    };
    const int           id_get()                const {
        return mid;
    };
    const int           verbosity_get()         const {
        return mverbosity;
    };
    const int           warnings_get()          const {
        return mwarnings;
    };
    const int           stackDepth_get()        const {
        return mstackDepth;
    };

    void                str_obj_set(string astr_val) {
        mstr_obj        = astr_val;
    };
    void                str_name_set(string astr_val) {
        mstr_name       = astr_val;
    };
    void                str_proc_set(int depth, string astr_proc) {
        mstr_proc[depth] = astr_proc;
    } ;
    void                id_set(int value) {
        mid             = value ;
    } ;
    void                verbosity_set(int value) {
        mverbosity      = value ;
    } ;
    void                warnings_set(int value) {
        mwarnings       = value ;
    } ;
    void                stackDepth_set(int value) {
        mstackDepth     = value ;
    } ;
    const string        str_proc_get()          const {
        return mstr_proc[stackDepth_get()];
    };
    const string        str_proc_get(int i) {
        return mstr_proc[i];
    };

    //
    // Core class methods
    //

    virtual int         run(void)                = 0;

};

class C_mpmProg_autodijk : public C_mpmProg {

  protected:

    int         mvertex_polar;
    int         mvertex_start;
    int         mvertex_step;
    int         mvertex_end;
    int         mvertex_total;
    int         m_costFunctionIndex;
    bool        mb_surfaceRipClear;
    int         mprogressIter;                  // Number of iterations to
                                                //+ loop before showing
                                                //+ progress to stdout
    float*      mpf_cost;                       // Cost as calculated by
                                                //+ autodijk
    string      mstr_costFileName;              // Parsed from the environment
                                                //+ structure
    string      mstr_costFullPath;              // Full path to cost file                                            
      
  public:
    C_mpmProg_autodijk(s_env* aps_env);
    ~C_mpmProg_autodijk(void);

    //
    // Access block
    //
    void        surfaceRipClear_set(bool avalue) {
            mb_surfaceRipClear  = avalue;
    };
    int         surfaceRipClear_get() {
            return(mb_surfaceRipClear);
    };
    void        vertexPolar_set(int avalue) {
            mvertex_polar       = avalue;
    };
    int         vertexPolar_get() {
            return(mvertex_polar);
    };
    void        vertexStart_set(int avalue) {
            mvertex_start       = avalue;
    };
    int         vertexStart_get() {
            return(mvertex_start);
    };
    void        vertexStep_set(int avalue) {
            mvertex_step        = avalue;
    };
    int         vertexStep_get() {
            return(mvertex_step);
    };
    void        vertexEnd_set(int avalue) {
            mvertex_end         = avalue;
    };
    int         vertexEnd_get() {
            return(mvertex_end);
    };
    void        progressIter_set(int avalue) {
            mprogressIter       = avalue;
    };
    int         progressIter_get() {
            return(mprogressIter);
    };
    void        print(void);

    //
    // Functional block
    //

    virtual int         run(void);
    float               cost_compute(int start, int end);
    e_FILEACCESS        CURV_fileWrite();
};


#endif //__C_MPM_PROG_H__


