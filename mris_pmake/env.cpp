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

#include "env.h"
#include "pathconvert.h"
#include "C_mpmProg.h"
#include "C_mpmOverlay.h"
#include "asynch.h"

#include "legacy.h"

#include <stdlib.h>
#include <assert.h>

extern string G_SELF;

#if 0
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

void
s_weights_print(
  s_weights& asw
) {
  int lw     = 20;
  int rw     = 20;
  std::_Ios_Fmtflags origFlags;

  origFlags  = cout.flags();
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

  cout.flags(origFlags);
}

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

void
s_Dweights_print(
  s_Dweights& asw
) {
  int lw     = 20;
  int rw     = 20;
  std::_Ios_Fmtflags origFlags;

  origFlags  = cout.flags();
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

  cout.flags(origFlags);
}

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
#endif

float
s_env_plyDepth_get(
    s_env&   ast_env
) {
    //
    // DESC
    //  The 'plyDepth' is a parameter used by the ROI
    //  mpmProg, corresponding to the 'f_radius' member
    //  variable. Getting the plyDepth via a function
    //  call is safer since the env.plyDepth will be
    //  depreciated. Moreover, the env.plyDepth is not
    //  synch'ed with the mpmProg.
    //

    C_mpmProg_ROI*      pC_ROI;
    if(!pC_ROI_cast(ast_env.pCmpmProg, pC_ROI))
        error_exit( "setting plyDepth",
                    "internal mpmProg is not of type ROI",
                    20);
    return pC_ROI->radius_get();
}

void
s_env_plyDepth_set(
    s_env&      ast_env,
    float       af_val
) {
    //
    // DESC
    //  The 'plyDepth' is a parameter used by the ROI
    //  mpmProg, corresponding to the 'f_radius' member
    //  variable. Setting the plyDepth via a function
    //  call is safer since the env.plyDepth will be
    //  depreciated. Moreover, the env.plyDepth is not
    //  synch'ed with the mpmProg.
    //

    C_mpmProg_ROI*      pC_ROI;
    if(!pC_ROI_cast(ast_env.pCmpmProg, pC_ROI))
        error_exit( "setting plyDepth",
                    "internal mpmProg is not of type ROI",
                    20);
    pC_ROI->radius_set(af_val);
}

float
s_env_plyIncrement_get(
    s_env&   ast_env
) {
    //
    // DESC
    //  The 'plyIncrement' is a parameter used by the ROI
    //  mpmProg, and used when creating "staggered" label
    //  saves. Each successive label save includes points
    //  corresponding to an additional <af_val> shell about
    //  a previous label.
    //

    C_mpmProg_ROI*      pC_ROI;
    if(!pC_ROI_cast(ast_env.pCmpmProg, pC_ROI))
        error_exit( "setting plyDepth",
                    "internal mpmProg is not of type ROI",
                    20);
    return pC_ROI->plyIncrement_get();
}

void
s_env_plyIncrement_set(
    s_env&      ast_env,
    float       af_val
) {
    //
    // DESC
    //  The 'plyIncrement' is a parameter used by the ROI
    //  mpmProg, and used when creating "staggered" label
    //  saves. Each successive label save includes points
    //  corresponding to an additional <af_val> shell about
    //  a previous label.
    //

    C_mpmProg_ROI*      pC_ROI;
    if(!pC_ROI_cast(ast_env.pCmpmProg, pC_ROI))
        error_exit( "setting plyDepth",
                    "internal mpmProg is not of type ROI",
                    20);
    pC_ROI->plyIncrement_set(af_val);
}

void
s_env_nullify(
    s_env&  st_env) {
    //
    // ARGS
    // st_env  			in  			environment to nullify
    //
    // DESC
    // "nullify", i.e. set relevant records files to NULL / zero / ""
    // This is something of a pseudo constructor.
    //
    // PRECONDITIONS
    // o Use *directly* after declaring a s_env struct
    //
    // HISTORY
    // 15 December 2004
    // o Initial design and coding.
    //

    st_env.pSTw                     = NULL;
    st_env.pSTDw                    = NULL;

    st_env.timeoutSec               = 0;
    st_env.port                     = 0;

    st_env.lw                       = -50;
    st_env.rw                       =  20;

    st_env.b_syslogPrepend          = false;
    st_env.b_exitOnDone             = false;
    st_env.pcsm_optionsFile         = NULL;
    st_env.pcsm_stdout              = NULL;
    st_env.pcsm_syslog              = NULL;
    st_env.pcsm_userlog             = NULL;
    st_env.pcsm_resultlog           = NULL;

    st_env.b_labelFile_save         = false;
    st_env.b_patchFile_save         = false;
    st_env.b_transitionPenalties    = false;

    st_env.str_workingDir           = "";
    st_env.str_patchFileName        = "";
    st_env.str_labelFileNameOS      = "";
    st_env.str_labelFileName        = "";
    st_env.str_costFileName         = "";
    st_env.b_costPathSave           = true;

    st_env.startVertex              = 0;
    st_env.endVertex                = 0;
    st_env.plyDepth                 = 0;

    st_env.pMS_active               = NULL;
    st_env.pMS_secondary            = NULL;
    st_env.pMS_primary              = NULL;
    st_env.pMS_auxillary            = NULL;

    st_env.b_useAbsCurvs            = false;
    st_env.b_surfacesKeepInSync     = false;
    st_env.b_surfacesClear          = true;
    st_env.b_costHistoryPreserve    = false;
    st_env.costFunc_do              = NULL;

    // Define the legacy cost functions for human readable setting / getting
    st_env.totalNumFunctions        = 4;
    st_env.ecf_current              = e_default;
    st_env.pstr_functionName        = new string[st_env.totalNumFunctions];
    st_env.pstr_functionName[0]     = "default";
    st_env.pstr_functionName[1]     = "unity";
    st_env.pstr_functionName[2]     = "euclid";
    st_env.pstr_functionName[3]     = "distance";

    // Module names
    st_env.vstr_mpm.clear();
    st_env.vstr_mpm.push_back("mpmProg");
    st_env.vstr_mpm.push_back("mpmOverlay");
	
    // Define the internal mpmProg modules for human readable setting / getting
    st_env.totalmpmProgs            = (int) empmprog;
    st_env.b_mpmProgUse             = false;
    st_env.empmProg_current         = emp_NOP;
    st_env.vstr_mpmProgName.clear();
    st_env.vstr_mpmProgName.push_back("NULL");
    st_env.vstr_mpmProgName.push_back("NOP");
    st_env.vstr_mpmProgName.push_back("pathFind");
    st_env.vstr_mpmProgName.push_back("autodijk");
    st_env.vstr_mpmProgName.push_back("autodijk_fast");
    st_env.vstr_mpmProgName.push_back("ROI");
    st_env.vstr_mpmProgName.push_back("externalMesh");
    st_env.totalmpmProgs	    = st_env.vstr_mpmProgName.size();
    st_env.pCmpmProg                = NULL; // Not yet created!
    // autodijk
    st_env.str_costCurvFile         = "autodijk.cost.crv";

    // Define the internal mpmOverlay modules 
    st_env.totalmpmOverlays	    = (int) empmoverlay;
    st_env.empmOverlay_current	    = emo_NOP;
    st_env.vstr_mpmOverlayName.clear();
    st_env.vstr_mpmOverlayName.push_back("LEGACY");
    st_env.vstr_mpmOverlayName.push_back("NULL");
    st_env.vstr_mpmOverlayName.push_back("NOP");
    st_env.vstr_mpmOverlayName.push_back("unity");
    st_env.vstr_mpmOverlayName.push_back("euclidean");
    st_env.vstr_mpmOverlayName.push_back("distance");
    st_env.vstr_mpmOverlayName.push_back("fscurvs");
    st_env.vstr_mpmOverlayName.push_back("curvature");
    st_env.totalmpmOverlays	    = st_env.vstr_mpmOverlayName.size();
    st_env.pCmpmOverlay		    = NULL; // Not yet created!

    // Define the active surface tracker
    st_env.totalNumSurfaces         = 3;
    st_env.esf_active               = e_workingCurvature;
    st_env.pstr_activeName          = new string[st_env.totalNumSurfaces];
    st_env.pstr_activeName[0]       = "workingCurvature";
    st_env.pstr_activeName[1]       = "workingSulcal";
    st_env.pstr_activeName[2]       = "auxillary";

}

#if 0
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
#endif

string 
s_env_HUP(
    s_env&			st_env,
    c_SSocket_UDP_receive**	pCSSocketReceive
) {
    // ARGS
    //	st_env			in/out		reference to environment struct
    //	pCSSocketReceive	in/out		pointer to "server" that 
    //						receives asynch UDP comms
    //
    // DESC 
    // 	This "method" handles a HUP "event" and reparses/rebuilds its core 
    // 	environment.
    //
    //	It returns an initialized pointer to the UDP comms handler.
    //
    
    //
    // If the system receives a "HUP" comms (which is the first-run 
    // default), then it will (re)parse the environment file (typically 
    // 'options.txt') in the working directory.
    //
    // NOTE: the working dir and optionsFile name are set to startup
    // defaults on very first run.
    //

    string		str_asynchComms         = "RUNPROG";
    string		str_optionsFQName;
    struct stat		st_fileInfo;
#if 0
    s_weights           st_costWeight;
    s_Dweights		st_DcostWeight;
#endif
    static int		oldport;

    str_optionsFQName = 	st_env.str_workingDir + 
	      			st_env.str_optionsFileName;
    if(stat(str_optionsFQName.c_str(), &st_fileInfo))
        error_exit("checking on the options file,",
                   "I couldn't access the options file. Does it exist?",
                    40);
    if(st_env.pcso_options) delete st_env.pcso_options;
    st_env.pcso_options = new C_scanopt(str_optionsFQName, e_EquLink);
	  
    // Parse the options file
    s_env_scan(st_env);
    // Initialize the mpmOverlay and mpmProg modules:
    // The mpmOverlay module determines edge costs
    s_env_mpmOverlaySetIndex(&st_env, st_env.empmOverlay_current);
    // The mpmProg module executes "programs" on the mesh based on 
    // costs
    s_env_mpmProgSetIndex(&st_env, st_env.empmProg_current);

#if 0
    // start legacy
    // WARNING!!
    s_env_costFctSetIndex (&st_env, st_env.empmOverlay_current);
    s_weights_scan(   st_costWeight,  *(st_env.pcso_options));
    s_Dweights_scan(  st_DcostWeight, *(st_env.pcso_options));
    st_env.pSTw                         = new s_weights;
    st_env.pSTDw                        = new s_Dweights;
    s_weights_copy(*(st_env.pSTw), 	st_costWeight);
    s_Dweights_copy(*(st_env.pSTDw),	st_DcostWeight);
    // LEGACY DEBUGGING!!
    // st_env.b_mpmOverlayUse		= false;
    // end legacy
#endif
    
    if(st_env.port != oldport) {
	if(*pCSSocketReceive) {
	    delete *pCSSocketReceive;
	    *pCSSocketReceive	= NULL;
	}
	oldport	= st_env.port;
    }

    if(!(*pCSSocketReceive))
	*pCSSocketReceive	= new c_SSocket_UDP_receive(
                                	st_env.port, st_env.timeoutSec);

    // Finally, if the system is set to run in "server" mode, return
    // a NULL str_asynchComms
    if(!st_env.b_exitOnDone) str_asynchComms = "NULL";
    
    return str_asynchComms;
}

void
s_env_mpmPrint(
    s_env&		ast_env,
    string 		astr_msg,
    e_MODULE		ae_module) 
    //
    // ARGS
    // 	ast_env			in		reference to environment
    // 	astr_msg		in		intro message
    //	ae_module		in		module enum to print
    //
    // DESCRIPTION
    //	Default informational printing of internal mpm module data.
    //
    // PRECONDITIONS
    // 	o ast_env must be non-NULL. 
    //   
    // POSTCONDITIONS
    //	o Prints various internal information relating to given mpm module 
    //	  type.
    //   
    // HISTORY
    // Late June 2010
    // o Initial design and coding.
    //
{

    int  		lw       	= ast_env.lw;
    int  		rw       	= ast_env.rw;
    int  		moduleIndex	= 0;
    vector<string>*	p_vstr;
    string		str_moduleType;
    string		str_moduleName;
    void*		p_moduleAddress	= NULL;

    p_vstr 		= &ast_env.vstr_mpmProgName;
    str_moduleType	= ast_env.vstr_mpm[ae_module];
    switch(ae_module) {
	case e_mpmProg:
	    p_vstr		= &ast_env.vstr_mpmProgName;
	    moduleIndex		= ast_env.empmProg_current;	
	    p_moduleAddress	= ast_env.pCmpmProg;
	    break;
	case e_mpmOverlay:
	    p_vstr		= &ast_env.vstr_mpmOverlayName;
	    moduleIndex		= ast_env.empmOverlay_current;
	    p_moduleAddress	= ast_env.pCmpmOverlay;
	    break;
	case emodule:
	    break;
    }

    cout << astr_msg	<< endl;

    colprintf(lw, rw, "Module type", "[ %s ]\n", str_moduleType.c_str());
    colprintf(lw, rw, "Module index:name" , "[ %d:%s ]\n",
            (int) moduleIndex,
            (*p_vstr)[moduleIndex].c_str());
    colprintf(lw, rw, "mpm module pointer:", "[ 0x%x ]\n", p_moduleAddress);
    if(p_moduleAddress)
	colprintf(lw, rw, "mpm module initialized:", "[ ok ]\n");
    else 
	colprintf(lw, rw, "mpm module initialized:", "[ no ]\n");

    lprintf(lw, "\nAll available %s modules (index:name):-\n", 
        str_moduleType.c_str());
    for (unsigned int i=0; i<p_vstr->size(); i++) {
        lprintf(lw, "%d: %s\n",
                i, (*p_vstr)[i].c_str());
    }
}

void
s_env_defaultsSet(
    s_env&              st_env
) {
    //
    // ARGS
    //  st_env          in/out          environment structure to be filled
    //
    // DESCRIPTION
    // Sets environment structure with basic (non-NULL) defaults. Its main
    // purpose is to set a basic number of values such that a complete
    // options file can be created.
    //
    // PRECONDITIONS
    // o st_env must be non-NULL. To be safe, this function also calls 
    //   's_env_nullify' on the env structure.
    //   
    // POSTCONDITIONS
    // o defaults are all set.
    // o complex structures are *not* initialized!
    // o data is stored to internal pcsm_optionsFile class.
    //   
    // HISTORY
    // 04 December, 2009
    // o Initial design and coding.
    //

    s_env_nullify(st_env);
    st_env.str_workingDir               = "./";
    st_env.str_optionsFileName          = "./options.txt";
    st_env.pcsm_optionsFile             = new C_SMessage( "",
                                                eSM_raw,
                                                st_env.str_optionsFileName,
                                                eSM_c,
                                                eOverwrite);
    st_env.pcsm_optionsFile->lw_set(25);
    st_env.pcsm_optionsFile->rw_set(55);

    st_env.startVertex                  = 0;
    st_env.endVertex                    = 0;

    st_env.str_primarySurfaceFileName           = "";
    st_env.str_secondarySurfaceFileName         = "";
    st_env.str_primaryCurvatureFileName         = "";
    st_env.str_secondaryCurvatureFileName       = "";

    st_env.b_patchFile_save             = false;
    st_env.b_labelFile_save             = true;

    st_env.b_useAbsCurvs                = false;
    st_env.str_costFileName             = "cost.txt";
    st_env.str_patchFileStem            = "dijk";
    st_env.str_labelFileStem            = "dijk";
    st_env.str_labelFileNameOS          = "dijkAuxSurface";

    st_env.b_syslogPrepend              = true;
    st_env.str_stdout                   = "./stdout.log";
    st_env.str_userMsgLog               = "./user_msg.log";
    st_env.str_sysMsgLog                = "localhost:1834";
    st_env.str_resultMsgLog             = "localhost:1835";
    st_env.serverControlPort            = 1701;
    st_env.timeoutSec                   = 60;

    //
    // LEGACY CODE
    st_env.pSTw                         = new s_weights;
    st_env.pSTDw                        = new s_Dweights;

    s_weights_setAll(*st_env.pSTw, 0.0);
    st_env.pSTw->wc                    	= 1.0;

    st_env.b_transitionPenalties        = false;
    s_Dweights_setAll(*st_env.pSTDw, 1.0);
    // LEGACY CODE
    //

    st_env.str_costCurvFile             = st_env.str_hemi+st_env.str_surface+\
                                            ".autodijk.crv";
    st_env.b_exitOnDone                 = true;
    st_env.b_costPathSave               = false;

}

void
s_env_optionsFile_write(
    s_env&              st_env,
    bool                ab_setToDefaults
) {
    //
    // ARGS
    //  st_env                  in      environment structure to be saved to
    //                                  + an optionsFile
    //  ab_setToDefaults        in      if true, call s_env_defaultsSet(...)
    //
    // DESCRIPTION
    // Writes a "quick-and-dirty" options file to disk.
    // 
    // POSTCONDITIONS
    // o the old options files is closed.
    // o a new options file is written to disk (overwriting the old).
    //
    // HISTORY
    // 04 December 2009
    // o Initial design and coding.
    //
    // 17 February 2011
    // o Closed/open of options file (to force overwrite).
    //

    C_SMessage*                 O;
    char                        pch_commandLine[65536];
    int                         i;

    if(ab_setToDefaults)        s_env_defaultsSet(st_env);

    O                           = st_env.pcsm_optionsFile;
    strcpy(pch_commandLine, "");
    for(i=0; i<st_env.argc; i++) {
        strcat(pch_commandLine, st_env.ppch_argv[i]);
        strcat(pch_commandLine, " ");
    }
    if(O) {
	delete O;
	O             		= new C_SMessage( "",
                                                eSM_raw,
                                                st_env.str_optionsFileName,
                                                eSM_c,
                                                eOverwrite);

        O->pprintf("\n#\n# auto-generated optionsFile\n#\n\n");

        O->pprintf("\n# Input surfaces and embedded curvature overlays\n");
        O->pcolprintf("surfaceFile",        " = %s\n",
                        st_env.str_primarySurfaceFileName.c_str());
        if(st_env.b_secondarySurface)
            O->pcolprintf("secondarySurfaceFile",     " = %s\n",
                        st_env.str_secondarySurfaceFileName.c_str());
        if(st_env.b_primaryCurvature)
            O->pcolprintf("primaryCurvature",      " = %s\n",
                        st_env.str_primaryCurvatureFileName.c_str());
        if(st_env.b_secondaryCurvature)
            O->pcolprintf("secondaryCurvature",   " = %s\n",
                        st_env.str_secondaryCurvatureFileName.c_str());
#if 0
        O->pprintf("\n# Start and End vertices\n");
        O->pcolprintf("startVertex",        " = %d\n",
                        st_env.startVertex);
        O->pcolprintf("endVertex",          " = %d\n",
                        st_env.endVertex);
#endif
        O->pprintf("\n# Control flags and settings\n");
        O->pcolprintf("controlPort",         " = %d\n",
                        st_env.serverControlPort);
        O->pcolprintf("timeoutSec",          " = %d\n",
                        st_env.timeoutSec);
        O->pcolprintf("b_transitionPenalties", " = %d\n",
                        st_env.b_transitionPenalties);
        O->pcolprintf("b_syslogPrepend",     " = %d\n",
                        st_env.b_syslogPrepend);
        O->pcolprintf("b_patchFile_save",   " = %d\n",
                        st_env.b_patchFile_save);
        O->pcolprintf("b_labelFile_save",   " = %d\n",
                        st_env.b_labelFile_save);
        O->pcolprintf("b_useAbsCurvs",      " = %d\n",
                        st_env.b_useAbsCurvs);
        O->pcolprintf("b_exitOnDone",       " = %d\n",
                        st_env.b_exitOnDone);
        O->pcolprintf("b_costPathSave",     " = %d\n",
                        st_env.b_costPathSave);
        O->pprintf("\n# Output file names and file stems\n");
        O->pcolprintf("costFile",           " = %s\n",
                        st_env.str_costFileName.c_str());
        O->pcolprintf("patchFile",          " = %s\n",
                        st_env.str_patchFileStem.c_str());
        O->pcolprintf("labelFile",          " = %s\n",
                        st_env.str_labelFileStem.c_str());
        O->pcolprintf("labelFileAuxSurface"," = %s\n",
                        st_env.str_labelFileNameOS.c_str());
        O->pprintf("\n# Messaging channels\n");
        O->pcolprintf("stdout",              " = %s\n",
                        st_env.str_stdout.c_str());
        O->pcolprintf("userMessages",        " = %s\n",
                        st_env.str_userMsgLog.c_str());
        O->pcolprintf("sysMessages",         " = %s\n",
                        st_env.str_sysMsgLog.c_str());
        O->pcolprintf("resultMessages",      " = %s\n",
                        st_env.str_resultMsgLog.c_str());
#if 0
        O->pprintf("\n# Weights\n");
        O->pcolprintf("wd",                 " = %f\n",
                        st_env.pSTw->wd);
        O->pcolprintf("wc",                 " = %f\n",
                        st_env.pSTw->wc);
        O->pcolprintf("wh",                 " = %f\n",
                        st_env.pSTw->wh);
        O->pcolprintf("wdc",                " = %f\n",
                        st_env.pSTw->wdc);
        O->pcolprintf("wdh",                " = %f\n",
                        st_env.pSTw->wdh);
        O->pcolprintf("wch",                " = %f\n",
                        st_env.pSTw->wch);
        O->pcolprintf("wdch",               " = %f\n",
                        st_env.pSTw->wdch);
        O->pcolprintf("wdir",               " = %f\n",
                        st_env.pSTw->wdir);
        O->pprintf("\n# Transitional penality weights\n");
        O->pcolprintf("Dwd",                " = %f\n",
                        st_env.pSTDw->Dwd);
        O->pcolprintf("Dwc",                " = %f\n",
                        st_env.pSTDw->Dwc);
        O->pcolprintf("Dwh",                " = %f\n",
                        st_env.pSTDw->Dwh);
        O->pcolprintf("Dwdc",               " = %f\n",
                        st_env.pSTDw->Dwdc);
        O->pcolprintf("Dwdh",               " = %f\n",
                        st_env.pSTDw->Dwdh);
        O->pcolprintf("Dwch",               " = %f\n",
                        st_env.pSTDw->Dwch);
        O->pcolprintf("Dwdch",              " = %f\n",
                        st_env.pSTDw->Dwdch);
        O->pcolprintf("Dwdir",              " = %f\n",
                        st_env.pSTDw->Dwdir);
#endif
        O->pprintf("\n# mpmProg\n");
        O->pcolprintf("mpmProgID",          " = %d\n",
                        st_env.empmProg_current);
        O->pcolprintf("mpmArgs",            " = %s\n",
                        st_env.str_mpmArgs.c_str());
        O->pcolprintf("costCurvFile",       " = %s\n",
                        st_env.str_costCurvFile.c_str());
        O->pprintf("\n# mpmOverlay\n");
        O->pcolprintf("mpmOverlayID",          " = %d\n",
                        st_env.empmOverlay_current);
        O->pcolprintf("mpmOverlayArgs",        " = %s\n",
                        st_env.str_mpmOverlayArgs.c_str());
        O->pprintf("\n# Debugging -- change at your own risk!\n");
        O->pcolprintf("b_mpmProgUse",		" = %d\n",
                        st_env.b_mpmProgUse);
        O->pcolprintf("b_mpmOverlayUse",	" = %d\n",
                        st_env.b_mpmOverlayUse);
        O->pcolprintf("argc", " = %d\n", st_env.argc);
        O->pcolprintf("argv", " = %s\n", pch_commandLine);

        O->dump();
    }
}

void
s_env_scan(
    s_env&          st_env
) {
  //
  // ARGS
  //    st_env          in/out          environment structure to be filled
  //
  // DESCRIPTION
  // Scans the options file structure for environment data.
  //
  // HISTORY
  // 17 November 2004
  // o Initial design and coding.
  //
  // 12 December 2004
  //  o Additional parse terms added.
  //
  // 15 December 2004
  // o C_SMessage
  //
  // 27 January 2005
  // o Added 'patchFileOS" (Orig Surface) processing.
  //
  // 08 March 2005
  // o Changed "orig" surface to "aux" surface.
  //
  // 05 April 2005
  //  o Added separate IO channel for result communications - comms
  //   is dependant on type of (user event-driven) function applied
  //   to the system. The need arose initially during genetic algorithm
  //   optimisation of cost weights - the system needed to communicate
  //   the result of a correlation operation to some external process.
  //
  // 29 October 2009
  // o Added 'b_useAbsCurvs'.
  // 
  // 23 January 2012
  // o Secondary surface and curvature reads filtered...
  //

  static int    calls                   = 0;

  static MRIS*  pMS_primary             = NULL;
  static MRIS*  pMS_secondary           = NULL;
  string        str_value               = "";
  string        str_surfaceFileName     = "";
  string        str_auxSurfaceFileName  = "";
  string        str_curvatureFileName   = "";
  string        str_secondaryCurvatureFile      = "";
  string        str_patchFileName       = "";
  string        str_labelFileName       = "";
  string        str_labelFileNameOS     = "";
  string        str_costFileName        = "";
  string        str_stdout              = "";
  string        str_userMsgFileName     = "";
  string        str_sysMsgFileName      = "";
  string        str_resultMsgFileName   = "";
  string        str_mpmArgs             = "-x";
  string        str_mpmOverlayArgs      = "-x";

  bool          b_useAbsCurvs           = true;
  bool          b_labelFile_save        = true;
  bool          b_transitionPenalties   = true;
  bool          b_patchFile_save        = true;
  bool          b_syslogPrepend         = true;
  bool          b_mpmProgUse            = false;
  bool          b_mpmOverlayUse		= false;
  bool          b_exitOnDone            = false;
  bool          b_costPathSave          = true;

  int           startVertex             = 0;
  int           endVertex               = 0;
  int           port                    = 0;
  int           timeoutSec              = 0;
  int           mpmProgID               = -1;
  int		mpmOverlayID		= -1;

  // These are used to re-read possibly new files if a HUP
  // is sent to the process with changed options file.
  static string str_surfaceFileNameOld          = "";
  static string str_auxSurfaceFileNameOld       = "";
  static string str_curvatureFileNameOld        = "";
  static string str_secondaryCurvatureFileOld   = "";
  static string str_stdoutOld                   = "";
  static string str_userMsgFileNameOld          = "";
  static string str_sysMsgFileNameOld           = "";
  static string str_resultMsgFileNameOld        = "";

  // mpmProg options
  static string str_costCurvFile        = "autodijk.cost.crv";
  C_scanopt	cso_options		= *st_env.pcso_options;

#if 0
  if (cso_options.scanFor("startVertex", &str_value))
    startVertex  = atoi(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a startVertex index.",
               10);
  if (cso_options.scanFor("endVertex", &str_value))
    endVertex  = atoi(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find a endVertex index.",
               11);
#endif

  if (cso_options.scanFor("surfaceFile", &str_value))
    str_surfaceFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find surfaceFile.",
               20);
  if (cso_options.scanFor("secondarySurfaceFile", &str_value)) {
    str_auxSurfaceFileName      =  str_value;
    st_env.b_secondarySurface   = true;
  } else
      st_env.b_secondarySurface = false;
  if (cso_options.scanFor("curvatureFile", &str_value)) {
      str_curvatureFileName     =  str_value;
      st_env.b_primaryCurvature = true;
  } else
      st_env.b_primaryCurvature = false;
  if (cso_options.scanFor("secondaryCurvature", &str_value)) {
    str_secondaryCurvatureFile          =  str_value;
    st_env.b_secondaryCurvature = true;
  } else
    st_env.b_secondaryCurvature = false;
  if (cso_options.scanFor("patchFile", &str_value))
    str_patchFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find patchFile.",
               24);
  if (cso_options.scanFor("labelFile", &str_value))
    str_labelFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find labelFile.",
               25);
  if (cso_options.scanFor("labelFileAuxSurface", &str_value))
    str_labelFileNameOS =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find labelFileAuxSurface.",
               26);
  if (cso_options.scanFor("costFile", &str_value))
    str_costFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find costFile.",
               27);

  if (cso_options.scanFor("b_labelFile_save", &str_value))
    b_labelFile_save =  atoi(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find b_labelFile_save.",
               30);
  if (cso_options.scanFor("b_patchFile_save", &str_value))
    b_patchFile_save =  atoi(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find b_patchFile_save.",
               31);
  if (cso_options.scanFor("b_transitionPenalties", &str_value))
    b_transitionPenalties =  atoi(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find b_transitionPenalties.",
               32);

  if (cso_options.scanFor("b_useAbsCurvs", &str_value))
    b_useAbsCurvs =  atoi(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find b_useAbsCurvs.",
               33);

  if (cso_options.scanFor("controlPort", &str_value))
    port   =  atoi(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find controlPort.",
               40);
  if (cso_options.scanFor("timeoutSec", &str_value))
    timeoutSec  = atoi(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find timeoutSec.",
               41);

  if (cso_options.scanFor("b_syslogPrepend", &str_value))
    b_syslogPrepend  =  atoi(str_value.c_str());
  else
    error_exit("scanning user options",
               "I couldn't find b_syslogPrepend.",
               50);
  if (cso_options.scanFor("userMessages", &str_value))
    str_userMsgFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find userMessages.",
               51);
  if (cso_options.scanFor("sysMessages", &str_value))
    str_sysMsgFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find sysMessages.",
               52);
  if (cso_options.scanFor("resultMessages", &str_value))
    str_resultMsgFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find resultMessages.",
               53);
  if (cso_options.scanFor("stdout", &str_value))
    str_stdout =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find stdout.",
               54);

  if (cso_options.scanFor("costCurvFile",       &str_value))
      str_costCurvFile  =  str_value;
  if (cso_options.scanFor("b_mpmProgUse",       &str_value))
      b_mpmProgUse      = atoi(str_value.c_str());
  if (cso_options.scanFor("b_mpmOverlayUse",	&str_value))
      b_mpmOverlayUse   = atoi(str_value.c_str());
  if (cso_options.scanFor("mpmProgID",          &str_value))
      mpmProgID         = atoi(str_value.c_str());
  if (cso_options.scanFor("mpmArgs",            &str_value))
      str_mpmArgs       = str_value;
  if (cso_options.scanFor("mpmOverlayArgs",     &str_value))
      str_mpmOverlayArgs= str_value;
  if (cso_options.scanFor("b_exitOnDone",       &str_value))
      b_exitOnDone      = atoi(str_value.c_str());
  if (cso_options.scanFor("b_costPathSave",     &str_value))
      b_costPathSave    = atoi(str_value.c_str());

  if (cso_options.scanFor("mpmOverlayID",	&str_value))
      mpmOverlayID	= atoi(str_value.c_str());
  else
      error_exit("scanning user options",
	    	   "I couldn't find mpmOverlayID",
	    	   54);
    
  st_env.b_syslogPrepend = b_syslogPrepend;
  e_SMessageIO esm_io  = eSM_cpp;
  int pos   = 0;

  if (str_stdout != str_stdoutOld) {
    if (st_env.pcsm_stdout)
      delete st_env.pcsm_stdout;
    pos  = str_userMsgFileName.find_first_of(":");
    // C-style IO for the internal "printf"ing...
    esm_io = ((unsigned) pos != (unsigned) string::npos) ? eSS : eSM_c;
    st_env.pcsm_stdout  = new C_SMessage( "",
                                           eSM_raw,
                                           str_stdout,
                                           esm_io);
    st_env.pcsm_stdout->str_syslogID_set(G_SELF);
    st_env.pcsm_stdout->b_syslogPrepend_set(true);
    st_env.pcsm_stdout->lw_set(st_env.lw);
    st_env.pcsm_stdout->rw_set(st_env.rw);
    st_env.pcsm_stdout->b_canPrint_set(true);
  }

  if (str_userMsgFileName != str_userMsgFileNameOld) {
    if (st_env.pcsm_userlog)
      delete st_env.pcsm_userlog;
    pos  = str_userMsgFileName.find_first_of(":");
    esm_io = ((unsigned) pos != (unsigned) string::npos) ? eSS : eSM_cpp;
    st_env.pcsm_userlog  = new C_SMessage( "",
                                           eSM_raw,
                                           str_userMsgFileName,
                                           esm_io);
    st_env.pcsm_userlog->str_syslogID_set(G_SELF);
  }

  if (str_sysMsgFileName != str_sysMsgFileNameOld) {
    pos  = str_sysMsgFileName.find_first_of(":");
    esm_io = ((unsigned) pos !=  (unsigned) string::npos) ? eSS : eSM_cpp;
    if (st_env.pcsm_syslog)
      delete st_env.pcsm_syslog;
    st_env.pcsm_syslog  = new C_SMessage( "",
                                          eSM_raw,
                                          str_sysMsgFileName,
                                          esm_io);
    st_env.pcsm_syslog->str_syslogID_set(G_SELF);
  }
  if (str_resultMsgFileName != str_resultMsgFileNameOld) {
    pos  = str_resultMsgFileName.find_first_of(":");
    esm_io = ((unsigned) pos !=  (unsigned) string::npos) ? eSS : eSM_cpp;
    if (st_env.pcsm_resultlog)
      delete st_env.pcsm_resultlog;
    st_env.pcsm_resultlog  = new C_SMessage( "",
                             eSM_raw,
                             str_resultMsgFileName,
                             esm_io);
    st_env.pcsm_resultlog->str_syslogID_set(G_SELF);
  }

  if (!calls)
    st_env.pcsm_syslog->timer(eSM_start);

  SLOUT("PARSING: env");
  string str_surfaceFileNameAbs;
  string str_auxSurfaceFileNameAbs;
  string str_curvatureFileNameAbs;
  string str_secondaryCurvatureFileAbs;

  if (str_surfaceFileName  != str_surfaceFileNameOld) {
    str_rel2absDirSpec_change(str_surfaceFileName, str_surfaceFileNameAbs);
    str_rel2absDirSpec_change(str_auxSurfaceFileName, str_auxSurfaceFileNameAbs);
    //cout << "-->" << str_surfaceFileNameAbs << endl;
    //cout << "-->" << str_auxSurfaceFileNameAbs << endl;
    ULOUT("Reading primary surface mesh...");
    pMS_primary  = MRISread( (char*)str_surfaceFileNameAbs.c_str());
    if(!pMS_primary) {
        nULOUT("\t\t\t[ failure ]\n");
        error_exit("reading the primary surface,",
                "I couldn't access the file. Does it exist? Are permissions OK?",
                 30);
    }
    nULOUT("\t\t\t[ ok ]\n");
    if(st_env.b_secondarySurface) {
        ULOUT("Reading secondary surface mesh...");
        pMS_secondary   = MRISread( (char*)str_surfaceFileNameAbs.c_str());
        if(!pMS_secondary) {
            nULOUT("\t\t\t[ failure ]\n");
            error_exit("reading the secondary surface,",
                    "I couldn't access the file. Does it exist? Are permissions OK?",
                     30);
        }
        nULOUT("\t\t\t\t[ ok ]\n");
    }
  }
  if (str_curvatureFileName  != str_curvatureFileNameOld &&
          st_env.b_primaryCurvature) {
    str_rel2absDirSpec_change(str_curvatureFileName,
            str_curvatureFileNameAbs);
    //cout << "-->" << str_curvatureFileNameAbs << endl;
    ULOUT("Mapping curvature texture on primary surface...");
    if(MRISreadCurvature(pMS_primary,
        (char*)str_curvatureFileNameAbs.c_str()) != NO_ERROR) {
        nULOUT("\t\t\t[ failure ]\n");
        error_exit("reading the primary curvature file,",
                "I couldn't access the file. Does it exist? Are permissions OK?",
                 30);
    }
    nULOUT("\t\t\t[ ok ]\n");
  }
  if (str_secondaryCurvatureFile  != str_secondaryCurvatureFileOld &&
          st_env.b_secondaryCurvature) {
    str_rel2absDirSpec_change(str_secondaryCurvatureFile,
            str_secondaryCurvatureFileAbs);
    //cout << "-->" << str_secondaryCurvatureFileAbs << endl;
    ULOUT("Mapping secondary texture on secondary surface...");
    if(MRISreadCurvature(pMS_secondary,
        (char*)str_secondaryCurvatureFileAbs.c_str()) != NO_ERROR) {
        nULOUT("\t\t\t[ failure ]\n");
        error_exit("reading the secondary curvature file,",
                "I couldn't access the file. Does it exist? Are permissions OK?",
                 30);
    }
    nULOUT("\t\t\t[ ok ]\n");
  }

    st_env.b_useAbsCurvs                = b_useAbsCurvs;
    st_env.b_labelFile_save             = b_labelFile_save;
    st_env.b_transitionPenalties        = b_transitionPenalties;
    st_env.b_patchFile_save             = b_patchFile_save;
    st_env.startVertex                  = startVertex;
    st_env.endVertex                    = endVertex;
    st_env.port                         = port;
    st_env.timeoutSec                   = timeoutSec;
    st_env.str_patchFileName            = str_patchFileName + ".patch";
    st_env.str_labelFileName            = str_labelFileName;
    st_env.str_labelFileNameOS          = str_labelFileNameOS;
    st_env.str_costFileName             = str_costFileName;
    st_env.b_costPathSave               = b_costPathSave;

    //    if(!calls) {
    st_env.pMS_primary                  = pMS_primary;
    st_env.pMS_active                   = pMS_primary;

    str_surfaceFileNameOld              = str_surfaceFileName;
    str_curvatureFileNameOld            = str_curvatureFileName;
    str_secondaryCurvatureFileOld       = str_secondaryCurvatureFile;
    str_userMsgFileNameOld              = str_userMsgFileName;
    str_sysMsgFileNameOld               = str_sysMsgFileName;
    //    }

    // mpmProg
    st_env.str_costCurvFile             = str_costCurvFile;
    st_env.b_mpmProgUse                 = b_mpmProgUse;
    st_env.str_mpmArgs                  = str_mpmArgs;
    st_env.empmProg_current             = (e_MPMPROG) mpmProgID;
    st_env.b_exitOnDone                 = b_exitOnDone;

    // mpmOverlay
    st_env.empmOverlay_current		= (e_MPMOVERLAY) mpmOverlayID;
    st_env.b_mpmOverlayUse              = b_mpmOverlayUse;
    st_env.str_mpmOverlayArgs           = str_mpmOverlayArgs;
    
    if(    !st_env.str_hemi.length()
        || !st_env.str_subject.length()
        || !st_env.str_primarySurfaceFileName.length()) {
        // Parse the surface text to extract the hemisphere 
        //+ and subject name
        vector<string>                  v_dir;
        vector<string>                  v_surface;
        vector<string>::iterator        i;
        string                          str_surfaceFile = "";
        int                             tokens  = 0;

        tokens = str_tokenize(str_surfaceFileName, v_dir, "/");
        if(!tokens)
            str_surfaceFile = str_surfaceFileName;
        else
            str_surfaceFile = v_dir.at(tokens-1);
        str_tokenize(str_surfaceFile, v_surface, ".");
        st_env.str_hemi                 = v_surface.at(0);
        st_env.str_surface              = v_surface.at(1);
        //st_env.str_mainSurfaceFileName  = v_surface.at(1);
        st_env.str_subject              = v_dir.at(tokens-3);
    }

    nSLOUT("\t\t\t\t\t\t\t[ ok ]\n");

    calls++;
}

bool
s_env_surfaceFile_set(
  s_env&   st_env,
  string   astr_fileName
) {
  //
  // PRECONDITIONS
  // o Internal MRIS structure should ideally speaking
  //   already exist.
  //
  // POSTCONDITIONS
  // o A new surfaceFile is read from file and housed
  //   internally, i.e. new pMS_curvature and pMS_sulcal
  //   are created.
  // o Curvature and sulcal height textures will need to be
  //    explicitly remapped to new surface.
  //
  // HISTORY
  // 26 April 2005
  // o Initial design and coding.
  //

  if (st_env.pMS_primary) {
    MRISfree(&st_env.pMS_primary);
  }
  ULOUT("Reading primary surface...");
  st_env.pMS_primary  = MRISread( (char*)astr_fileName.c_str());
  nULOUT("\t\t\t[ ok ]\n");
  st_env.pMS_active   = st_env.pMS_primary;

  return true;
}

bool
s_env_surfaceCurvature_set(
  s_env&   st_env,
  string   astr_fileName
) {
  //
  // PRECONDITIONS
  // o Internal MRIS structure pMS_curvature should
  //   already exist.
  //
  // POSTCONDITIONS
  // o The curvatures described in astr_filenName are mapped
  //   onto pMS_curvature.
  //
  // HISTORY
  // 26 April 2005
  // o Initial design and coding.
  //

  if (!st_env.pMS_primary)
    return false;

  ULOUT("Mapping curvature texture on primary surface...");
  MRISreadCurvature(st_env.pMS_primary,
                    (char*)astr_fileName.c_str());
  nULOUT("\t\t\t[ ok ]\n");

  return true;
}

bool
s_env_secondarySurface_setCurvature(
  s_env&   st_env,
  string   astr_fileName
) {
  //
  // PRECONDITIONS
  // o Internal MRIS structure pMS_sulcal should
  //   already exist.
  //
  // POSTCONDITIONS
  // o The sulcal heights described in astr_filenName are mapped
  //   onto pMS_sulcal.
  //
  // HISTORY
  // 26 April 2005
  // o Initial design and coding.
  //

  if (!st_env.pMS_secondary)
    return false;

  ULOUT("Mapping curvature texture on secondary surface...");
  MRISreadCurvature(st_env.pMS_secondary,
                    (char*)astr_fileName.c_str());
  nULOUT("\t\t\t[ ok ]\n");

  return true;
}

bool
s_env_auxSurfaceCurvature_set(
  s_env&   st_env,
  string   astr_fileName
) {
  //
  // PRECONDITIONS
  // o Internal MRIS structure pMS_auxSurface should
  //   already exist.
  //
  // POSTCONDITIONS
  // o The curvatures described in astr_filenName are mapped
  //   onto pMS_auxSurface.
  //
  // HISTORY
  // 26 April 2005
  // o Initial design and coding.
  //

  if (!st_env.pMS_secondary)
    return false;

  ULOUT("Mapping curvature texture on auxillary surface...");
  MRISreadCurvature(st_env.pMS_secondary,
                    (char*)astr_fileName.c_str());
  nULOUT("\t\t\t[ ok ]\n");

  return true;
}

bool
s_env_auxSurfaceFile_set(
  s_env&   st_env,
  string   astr_fileName
) {
  //
  // PRECONDITIONS
  // o Internal MRIS structure should ideally speaking
  //   already exist.
  //
  // POSTCONDITIONS
  // o A new auxSurfaceFile is read from file and housed
  //   internally, ie new pMS_auxSurface.
  // o Curvature textures will need to be
  //    explicitly remapped to new surface.
  //
  // HISTORY
  // 26 April 2005
  // o Initial design and coding.
  //

  if (st_env.pMS_secondary)
    MRISfree(&st_env.pMS_secondary);
  ULOUT("Reading surface for auxillary curvature...");
  st_env.pMS_secondary = MRISread( (char*)astr_fileName.c_str());
  nULOUT("\t\t\t[ ok ]\n");

  return true;
}

void
s_env_log_file_changeTo(
  s_env&  ast_env,
  e_LOG  ae_log,
  string  astr_newName
) {
  //
  // PRECONDITIONS
  //  o None during normal operation.
  //
  // POSTCONDITIONS
  //  o The log file spec'd by <ae_log> has its file_changeTo() method
  //   called.
  //
  // HISTORY
  // 06 April 2005
  //  o Initial design and coding.
  //

  C_SMessage*  pcsm = NULL;
  int   pos;
  e_SMessageIO esm_io;

  pos  = astr_newName.find_first_of(":");
  esm_io  = ((unsigned) pos !=  (unsigned) string::npos) ? eSS : eSM_cpp;

  switch (ae_log) {
  case e_user:
    pcsm = ast_env.pcsm_userlog;
    break;
  case e_sys:
    pcsm = ast_env.pcsm_syslog;
    break;
  case e_result:
    pcsm = ast_env.pcsm_resultlog;
    break;
  }

  pcsm->file_changeTo(astr_newName, esm_io);
}

bool
s_env_b_surfacesKeepInSync_get(
  s_env&   ast_env
) {
  return ast_env.b_surfacesKeepInSync;
};

void
s_env_b_surfacesKeepInSync_set(
  s_env&   ast_env,
  bool   b_val
) {
  ast_env.b_surfacesKeepInSync = b_val;
};

bool
s_env_b_surfacesClear_get(
  s_env&   ast_env
) {
  return ast_env.b_surfacesClear;
};

void
s_env_b_surfacesClear_set(
  s_env&   ast_env,
  bool   b_val
) {
  ast_env.b_surfacesClear = b_val;
};

void
s_env_activeSurfaceList(
  s_env& ast_env
) {
  int  lw  = 30;
  int  rw = 20;

  CW(lw, "Current active surface:");
  CWn(rw, ast_env.pstr_activeName[ast_env.esf_active]);

  cout << "All available active surfaces:" << endl;
  for (int i=0; i<ast_env.totalNumSurfaces; i++) {
    cout << "Surface index: " << i << ": ";
    cout << ast_env.pstr_activeName[i] << endl;
  }
}

void
s_env_activeSurfaceSetIndex(
  s_env*  apst_env,
  int  aindex
) {
  switch (aindex) {
  case 0:
    apst_env->pMS_active = apst_env->pMS_primary;
    apst_env->esf_active = e_workingCurvature;
    break;
  case 1:
    apst_env->pMS_active = apst_env->pMS_auxillary;
    apst_env->esf_active = e_workingSulcal;
    break;
  case 2:
    apst_env->pMS_active = apst_env->pMS_secondary;
    apst_env->esf_active = e_auxillary;
    break;
  default:
    apst_env->pMS_active = apst_env->pMS_primary;
    apst_env->esf_active = e_workingCurvature;
    break;
  }
}

int
s_env_mpmProgSetIndex(
    s_env*      apst_env,
    int         aindex
) {
    int         ret             = -1;
    int         lw              = apst_env->lw;
    int         rw              = apst_env->rw;

    apst_env->b_mpmProgUse	= true;
    if(apst_env->pCmpmProg) {
	e_MPMPROG	e_prog = apst_env->empmProg_current;
        lprintf(lw, "Non-NULL mpmProg pointer detected.\n");
        lprintf(lw, "Deleting existing mpmProg '%s'...", 
	        	apst_env->vstr_mpmProgName[e_prog].c_str());
        delete apst_env->pCmpmProg;
        lprintf(rw, "[ ok ]\n");
    }

    switch (aindex) {
      case emp_NOP:
        apst_env->pCmpmProg         	= new C_mpmProg_NOP(apst_env);
        break;
      case emp_pathFind:
	apst_env->pCmpmProg		= new C_mpmProg_pathFind(apst_env,
	    						apst_env->startVertex,
	    						apst_env->endVertex);
          // Check for any command-line spec'd args for this mpmProg:
          if(apst_env->str_mpmArgs != "-x") {
            C_scanopt                   cso_mpm(apst_env->str_mpmArgs, ",",
                                                e_EquLink, "", ":");
            string      str_vertexStart = "0";
	    string	str_vertexEnd	= "0";
            int         vertexStart     = 0;
	    int		vertexEnd	= apst_env->pMS_primary->nvertices;
            C_mpmProg_pathFind* pC_pathFind	= NULL;
            pC_pathFind_cast(apst_env->pCmpmProg, pC_pathFind);
            if(cso_mpm.scanFor("vertexStart", &str_vertexStart)) {
                vertexStart     = atoi(str_vertexStart.c_str());
                pC_pathFind->vertexStart_set(vertexStart);
            }
            if(cso_mpm.scanFor("vertexEnd", &str_vertexEnd)) {
                vertexEnd     = atoi(str_vertexEnd.c_str());
                pC_pathFind->vertexEnd_set(vertexEnd);
            }
	    s_env_optionsFile_write(*apst_env);
          }
	break;
      case emp_autodijk: 
      case emp_autodijk_fast:
	if(aindex == emp_autodijk)
	  apst_env->pCmpmProg	    	= new C_mpmProg_autodijk(apst_env);
	if(aindex == emp_autodijk_fast)
          apst_env->pCmpmProg		= new C_mpmProg_autodijk_fast(apst_env);
        // Check for any command-line spec'd args for the 'autodijk' mpmProg:
        if(apst_env->str_mpmArgs != "-x") {
            C_scanopt                   cso_mpm(apst_env->str_mpmArgs, ",",
                                                e_EquLink, "", ":");
            string      str_vertexPolar         = "0";
            string      str_vertexStart         = "0";
            string      str_vertexStep          = "1";
            string      str_vertexEnd           = "0";
            string      str_costCurvStem        = "";
            string      str_worldMapCreate      = "";
            bool        b_worldMapCreate        = false;
            int         vertexPolar             = 0;
            int         vertexStart             = 0;
            int         vertexEnd               = 0;
            int         vertexStep              = 1;
            C_mpmProg_autodijk* pC_autodijk     = NULL;
            pC_autodijk_cast(apst_env->pCmpmProg, pC_autodijk);
            if(cso_mpm.scanFor("vertexPolar", &str_vertexPolar)) {
                vertexPolar     = atoi(str_vertexPolar.c_str());
                pC_autodijk->vertexPolar_set(vertexPolar);
            }
            if(cso_mpm.scanFor("vertexStart", &str_vertexStart)) {
                vertexStart     = atoi(str_vertexStart.c_str());
                pC_autodijk->vertexStart_set(vertexStart);
            }
            if(cso_mpm.scanFor("vertexStep", &str_vertexStep)) {
                vertexStep      = atoi(str_vertexStep.c_str());
                pC_autodijk->vertexStep_set(vertexStep);
            }
            if(cso_mpm.scanFor("vertexEnd", &str_vertexEnd)) {
                vertexEnd       = atoi(str_vertexEnd.c_str());
                pC_autodijk->vertexEnd_set(vertexEnd);
            }
            if(cso_mpm.scanFor("worldMapCreate", &str_worldMapCreate)) {
                b_worldMapCreate = atoi(str_worldMapCreate.c_str());
                pC_autodijk->worldMap_set(b_worldMapCreate);
            }
            if(cso_mpm.scanFor("costCurvStem", &str_costCurvStem)) {
                apst_env->str_costCurvFile =
                        apst_env->str_hemi + "."                +
                        apst_env->str_surface + "."             +
                        "autodijk-"                             +
                        str_costCurvStem                        + 
                        ".crv";
                pC_autodijk->costFile_set(apst_env->str_costCurvFile);
                s_env_optionsFile_write(*apst_env);
            }
          }
        break;
      case emp_ROI:
          apst_env->pCmpmProg             = new C_mpmProg_ROI(apst_env);
            // Check for any command-line spec'd args for this mpmProg:
            if(apst_env->str_mpmArgs != "-x") {
              C_scanopt                   cso_mpm(apst_env->str_mpmArgs, ",",
                                                  e_EquLink, "", ":");
              float       f_radius              = 0.0;
              bool        b_saveStaggered       = false;
              bool	  b_borderRegion	= false;
              string      str_option            = "";
              C_mpmProg_ROI* pC_ROI             = NULL;
              pC_ROI_cast(apst_env->pCmpmProg, pC_ROI);
              if(cso_mpm.scanFor("radius", &str_option)) {
                  f_radius              = atof(str_option.c_str());
                  pC_ROI->radius_set(f_radius);
              }
              if(cso_mpm.scanFor("vertexFile", &str_option)) {
                  if(!pC_ROI->vertexFile_load(str_option))
                      error_exit("reading the vertexFile",
                                  "a file access error occurred", 10);
              }
              if(cso_mpm.scanFor("labelFile", &str_option)) {
                  if(!pC_ROI->labelFile_load(str_option))
                      error_exit("reading the labelFile",
                                  "a file access error occurred", 10);
              }
              if(cso_mpm.scanFor("borderOnly", &str_option)) {
        	  b_borderRegion 	= atoi(str_option.c_str());
        	  pC_ROI->boundaryOnly(b_borderRegion);
              }
              if(cso_mpm.scanFor("plySaveStaggered", &str_option)) {
                  b_saveStaggered       = atoi(str_option.c_str());
                  pC_ROI->plySaveStaggered_set(b_saveStaggered);
              }
              s_env_optionsFile_write(*apst_env);
            }
          break;
      case emp_externalMesh:
          apst_env->pCmpmProg             = new C_mpmProg_externalMesh(apst_env);
            // Check for any command-line spec'd args for this mpmProg:
            if(apst_env->str_mpmArgs != "-x") {
              C_scanopt                   cso_mpm(apst_env->str_mpmArgs, ",",
                                                  e_EquLink, "", ":");
              float       f_radius              = 0.0;
              bool        b_saveStaggered       = false;
              string      str_option            = "";
              C_mpmProg_ROI* pC_ROI             = NULL;
              pC_ROI_cast(apst_env->pCmpmProg, pC_ROI);
              if(cso_mpm.scanFor("radius", &str_option)) {
                  f_radius              = atof(str_option.c_str());
                  pC_ROI->radius_set(f_radius);
              }
              if(cso_mpm.scanFor("vertexFile", &str_option)) {
                  if(!pC_ROI->vertexFile_load(str_option))
                      error_exit("reading the vertexFile",
                                  "a file access error occurred", 10);
              }
              if(cso_mpm.scanFor("labelFile", &str_option)) {
                  if(!pC_ROI->labelFile_load(str_option))
                      error_exit("reading the labelFile",
                                  "a file access error occurred", 10);
              }
              if(cso_mpm.scanFor("plySaveStaggered", &str_option)) {
                  b_saveStaggered       = atoi(str_option.c_str());
                  pC_ROI->plySaveStaggered_set(b_saveStaggered);
              }
              s_env_optionsFile_write(*apst_env);
            }
          break;
    default:
        apst_env->empmProg_current      = (e_MPMPROG) 0;
	apst_env->pCmpmProg		= new C_mpmProg_NOP(apst_env);
        break;
    }

    if(aindex < 0 || aindex >= apst_env->totalmpmProgs)
        ret       = -1;
    else
        ret       = aindex;
    apst_env->empmProg_current	= (e_MPMPROG) ret;
    return ret;
}

int
s_env_mpmOverlaySetIndex(
    s_env*      apst_env,
    int         aindex
) {
    int         ret                     = -1;
    int         lw                      = apst_env->lw;
    int         rw                      = apst_env->rw;
    string      str_curvatureFile       = "-x";

    if(apst_env->pCmpmOverlay) {
	e_MPMOVERLAY e_overlay = apst_env->empmOverlay_current;
    	lprintf(lw, "Non-NULL mpmOverlay pointer detected.\n");
        lprintf(lw, "Deleting existing mpmOverlay '%s'...", 
	    apst_env->vstr_mpmOverlayName[e_overlay].c_str());
        delete apst_env->pCmpmOverlay;
        lprintf(rw, "[ ok ]\n");
    }
    
    switch ((e_MPMOVERLAY) aindex) {
      case emo_LEGACY:
        lprintf(lw, "Forcing overlay engine to LEGACY mode.\n");
        apst_env->b_mpmOverlayUse       = false;
        break;
      case emo_NULL:
	break;
      case emo_NOP:
        apst_env->pCmpmOverlay	= new C_mpmOverlay_NOP(apst_env);
        break;
      case emo_unity:
	apst_env->pCmpmOverlay	= new C_mpmOverlay_unity(apst_env);
	break;
      case emo_distance:
	apst_env->pCmpmOverlay	= new C_mpmOverlay_distance(apst_env);
	break;
      case emo_euclidean:
	apst_env->pCmpmOverlay	= new C_mpmOverlay_euclidean(apst_env);
        break;
      case emo_fscurvs:
        break;
      case emo_curvature:
        // Check for any overlay args, specifically the name of the curv file
        if(apst_env->str_mpmOverlayArgs != "-x") {
            C_scanopt                   cso_mpm(apst_env->str_mpmOverlayArgs, 
                                                ",",
                                                e_EquLink, "", ":");
            if(!cso_mpm.scanFor("primaryCurvature", &str_curvatureFile)) {
                error_exit ("checking for curvature file to read",
                            "it seems no file was specified",
                            10);
            }
        }
        apst_env->pCmpmOverlay  = new C_mpmOverlay_curvature(apst_env,
                                                             str_curvatureFile);
        break;
      default:
        apst_env->empmOverlay_current   = (e_MPMOVERLAY) 0;
        apst_env->pCmpmOverlay		= new C_mpmOverlay_NOP(apst_env);
        break;
  }
  if(aindex < -2 || aindex >= apst_env->totalmpmOverlays)
      ret       = -2;
  else
      ret       = aindex;
  apst_env->empmOverlay_current  = (e_MPMOVERLAY) ret;
  return ret;
}

#if 0
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
#endif

float 
s_env_edgeCostFind(
    s_env&		ast_env,
    int			avertexi,
    int			avertexj
) {
    //
    // ARGS
    //	ast_env			in		environment structure
    //	avertexi		in		"start" vertex
    //	avertexj		in		"end" vertex
    //
    // DESC 
    //	This "method" determines the cost of an edge between vertices
    //	<avertexi> and <avertexj>.
    //
    // PRECONDITIONS
    //	o The vertices are assumed to have only one edge between them, i.e.
    //	  there are no vertices between <avertexi> and <avertexj>.
    //  o The vertex indices are "absolute" indices.
    //
    // POSTCONDITIONS
    //	o The cost of the edge is returned.
    //	o The ast_env st_iterInfo structure is updated.
    //
    // HISTORY
    // 08 July 2010
    // 	o Initial design and coding.
    //

    float 	f_cost			= 0.0;
    bool 	b_relNextReference	= false;

    if(!ast_env.b_mpmOverlayUse)
    	f_cost = ast_env.costFunc_do(	ast_env, &*(ast_env.pst_iterInfo),
        		    		avertexi, avertexj,
        		    		b_relNextReference);
    else
	f_cost = ast_env.pCmpmOverlay->costEdge_calc(avertexi, avertexj);

    return f_cost;
}

#if 0
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
    VERTEX*       v_c;
    VERTEX*       v_n;
    float         dist, ave_curv, curv, max_height, max_curv, cost;
    s_weights*    pSTw = st_env.pSTw;
    MRIS*         surf = st_env.pMS_primary;

    v_c = &surf->vertices[vno_c];
    if (b_relNextReference) {
        vno_n = v_c->v[j];
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
    VERTEX*     	v_c         = NULL;
    MRIS*    		surf        = st_env.pMS_primary;
    const char*       	pch_proc    = "costFunc_distanceReturn(...)";
    char        	pch_txt[65536];
    static bool 	b_warned    = false;

    v_c         = &surf->vertices[vno_c];

    if(!b_relNextReference) {
	int 	jrelcount;
	int 	jrel = 0;
	for(jrelcount=0; jrelcount< v_c->vnum; jrelcount++) {
	    if(v_c->v[jrelcount] == j) {
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

    VERTEX*     v_c         = NULL;
    VERTEX*     v_n         = NULL;
    MRIS*       surf        = st_env.pMS_primary;
    static int  calls       = 0;
    const char*       pch_proc    = "costFunc_EuclideanReturn(...)";
    char        pch_txt[65536];
    static bool b_warned = false;

    v_c = &surf->vertices[vno_c];
    if (b_relNextReference) {
        vno_n = v_c->v[j];
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
#endif

/* eof */
