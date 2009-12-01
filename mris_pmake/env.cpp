/***************************************************************************
 *   Copyright (C) 2004 by Rudolph Pienaar / Christian Haselgrove          *
 *   {ch|rudolph}@nmr.mgh.harvard.edu                                      *
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
// $Id: env.cpp,v 1.13 2009/12/01 19:46:09 rudolph Exp $

#include "env.h"
#include "pathconvert.h"
#include "C_mpmProg.h"

#include <stdlib.h>
#include <assert.h>

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

void
s_env_nullify(
    s_env&  st_env) {
    //
    // ARGS
    // st_env  in  environment to nullify
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

    st_env.lw                       = 40;
    st_env.rw                       = 20;

    st_env.b_syslogPrepend          = false;
    st_env.stdout                   = new C_SMessage( "",
                                           eSM_raw,
                                           "stdout",
                                           eSM_c);
    st_env.stdout->b_syslogPrepend_set(true);
    st_env.stdout->lw_set(st_env.lw);
    st_env.stdout->rw_set(st_env.rw);
    st_env.stdout->b_canPrint_set(true);
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
    st_env.pMS_auxSurface           = NULL;
    st_env.pMS_curvature            = NULL;
    st_env.pMS_sulcal               = NULL;

    st_env.b_useAbsCurvs            = false;
    st_env.b_surfacesKeepInSync     = false;
    st_env.b_surfacesClear          = true;
    st_env.b_costHistoryPreserve    = false;
    st_env.costFunc_do              = NULL;

    // Define the cost functions for human readable setting / getting
    st_env.totalNumFunctions        = 4;
    st_env.ecf_current              = e_default;
    st_env.pstr_functionName        = new string[st_env.totalNumFunctions];
    st_env.pstr_functionName[0]     = "default";
    st_env.pstr_functionName[1]     = "unity";
    st_env.pstr_functionName[2]     = "euclid";
    st_env.pstr_functionName[3]     = "distance";

    // Define the internal mpmProg modules for human readable setting / getting
    st_env.totalmpmProgs            = 1;
    st_env.b_mpmProgUse             = false;
    st_env.empm_current             = e_autodijk;
    st_env.pstr_mpmProgName         = new string[st_env.totalmpmProgs];
    st_env.pstr_mpmProgName[0]      = "autodijk";
    st_env.pCmpmProg                = NULL; // Not yet created!
    // autodijk
    st_env.str_costCurvFile         = "autodijk.cost.crv";

    // Define the active surface tracker
    st_env.totalNumSurfaces         = 3;
    st_env.esf_active               = e_workingCurvature;
    st_env.pstr_activeName          = new string[st_env.totalNumSurfaces];
    st_env.pstr_activeName[0]       = "workingCurvature";
    st_env.pstr_activeName[1]       = "workingSulcal";
    st_env.pstr_activeName[2]       = "auxillary";

}

void
s_env_scan(
    s_env&          st_env,
    C_scanopt&      cso_options
) {
  //
  // ARGS
  //    cso_options     in              scanopt structure to be parsed
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

  static int    calls                   = 0;

  static MRIS*  pMS_curvature           = NULL;
  static MRIS*  pMS_auxSurface          = NULL;
  static MRIS*  pMS_sulcal              = NULL;
  string        str_value               = "";
  string        str_surfaceFileName     = "";
  string        str_auxSurfaceFileName  = "";
  string        str_curvatureFileName   = "";
  string        str_sulcalFileName      = "";
  string        str_patchFileName       = "";
  string        str_labelFileName       = "";
  string        str_labelFileNameOS     = "";
  string        str_costFileName        = "";
  string        str_userMsgFileName     = "";
  string        str_sysMsgFileName      = "";
  string        str_resultMsgFileName   = "";

  bool          b_useAbsCurvs           = true;
  bool          b_labelFile_save        = true;
  bool          b_transitionPenalties   = true;
  bool          b_patchFile_save        = true;
  bool          b_syslogPrepend         = true;

  int           startVertex             = 0;
  int           endVertex               = 0;
  int           port                    = 0;
  int           timeoutSec              = 0;

  // These are used to re-read possibly new files if a HUP
  // is sent to the process with changed options file.
  static string str_surfaceFileNameOld      = "";
  static string str_auxSurfaceFileNameOld   = "";
  static string str_curvatureFileNameOld    = "";
  static string str_sulcalFileNameOld       = "";
  static string str_userMsgFileNameOld      = "";
  static string str_sysMsgFileNameOld       = "";
  static string str_resultMsgFileNameOld    = "";

  // mpmProg options
  static string str_costCurvFile            = "autodijk.cost.crv";
  
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

  if (cso_options.scanFor("surfaceFile", &str_value))
    str_surfaceFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find surfaceFile.",
               20);
  if (cso_options.scanFor("auxSurfaceFile", &str_value))
    str_auxSurfaceFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find auxSurfaceFile.",
               21);
  if (cso_options.scanFor("curvatureFile", &str_value))
    str_curvatureFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find curvatureFile.",
               22);
  if (cso_options.scanFor("sulcalHeightFile", &str_value))
    str_sulcalFileName =  str_value;
  else
    error_exit("scanning user options",
               "I couldn't find sulcalHeightFile.",
               23);
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

  if (cso_options.scanFor("costCurvFile", &str_value))
    str_costCurvFile    =  str_value;

  st_env.b_syslogPrepend = b_syslogPrepend;
  e_SMessageIO esm_io  = eSM_cpp;
  int pos   = 0;
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
  string  str_surfaceFileNameAbs;
  string str_auxSurfaceFileNameAbs;
  string  str_curvatureFileNameAbs;
  string str_sulcalFileNameAbs;

  if (str_surfaceFileName  != str_surfaceFileNameOld) {
    str_rel2absDirSpec_change(str_surfaceFileName, str_surfaceFileNameAbs);
    str_rel2absDirSpec_change(str_auxSurfaceFileName, str_auxSurfaceFileNameAbs);
    //cout << "-->" << str_surfaceFileNameAbs << endl;
    //cout << "-->" << str_auxSurfaceFileNameAbs << endl;
    ULOUT("Reading surface for primary curvature...");
    pMS_curvature  = MRISread( (char*)str_surfaceFileNameAbs.c_str());
    nULOUT("\t\t\t[ ok ]\n");
    ULOUT("Reading surface for primary sulcal...");
    pMS_sulcal   = MRISread( (char*)str_surfaceFileNameAbs.c_str());
    nULOUT("\t\t\t\t[ ok ]\n");
    ULOUT("Reading surface for auxillary curvature...");
    pMS_auxSurface  = MRISread( (char*)str_auxSurfaceFileNameAbs.c_str());
    nULOUT("\t\t\t[ ok ]\n");
  }
  if (str_curvatureFileName  != str_curvatureFileNameOld) {
    str_rel2absDirSpec_change(str_curvatureFileName, str_curvatureFileNameAbs);
    //cout << "-->" << str_curvatureFileNameAbs << endl;
    ULOUT("Mapping curvature texture on primary surface...");
    MRISreadCurvature(pMS_curvature,  (char*)str_curvatureFileNameAbs.c_str());
    nULOUT("\t\t\t[ ok ]\n");
    ULOUT("Mapping curvature texture on auxillary surface...");
    MRISreadCurvature(pMS_auxSurface,  (char*)str_curvatureFileNameAbs.c_str());
    nULOUT("\t\t[ ok ]\n");
  }
  if (str_sulcalFileName  != str_sulcalFileNameOld) {
    str_rel2absDirSpec_change(str_sulcalFileName, str_sulcalFileNameAbs);
    //cout << "-->" << str_sulcalFileNameAbs << endl;
    ULOUT("Mapping sulcal texture on primary surface...");
    MRISreadCurvature(pMS_sulcal, (char*)str_sulcalFileNameAbs.c_str());
    nULOUT("\t\t\t[ ok ]\n");
  }

    st_env.b_useAbsCurvs              = b_useAbsCurvs;
    st_env.b_labelFile_save           = b_labelFile_save;
    st_env.b_transitionPenalties      = b_transitionPenalties;
    st_env.b_patchFile_save           = b_patchFile_save;
    st_env.startVertex                = startVertex;
    st_env.endVertex                  = endVertex;
    st_env.port                       = port;
    st_env.timeoutSec                 = timeoutSec;
    st_env.str_patchFileName          = str_patchFileName + ".patch";
    st_env.str_labelFileName          = str_labelFileName;
    st_env.str_labelFileNameOS        = str_labelFileNameOS;
    st_env.str_costFileName           = str_costFileName;
    st_env.b_costPathSave             = true;

    //    if(!calls) {
    st_env.pMS_active                 = pMS_curvature;
    st_env.pMS_auxSurface             = pMS_auxSurface;
    st_env.pMS_sulcal                 = pMS_sulcal;
    st_env.pMS_curvature              = pMS_curvature;

    str_surfaceFileNameOld            = str_surfaceFileName;
    str_curvatureFileNameOld          = str_curvatureFileName;
    str_sulcalFileNameOld             = str_sulcalFileName;
    str_userMsgFileNameOld            = str_userMsgFileName;
    str_sysMsgFileNameOld             = str_sysMsgFileName;
    //    }

    // mpmProg
    st_env.str_costCurvFile           = str_costCurvFile;

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

  if (st_env.pMS_curvature) {
    MRISfree(&st_env.pMS_curvature);
    MRISfree(&st_env.pMS_sulcal);
  }
  ULOUT("Reading surface for primary curvature...");
  st_env.pMS_curvature  = MRISread( (char*)astr_fileName.c_str());
  nULOUT("\t\t\t[ ok ]\n");
  ULOUT("Reading surface for primary sulcal...");
  st_env.pMS_sulcal  = MRISread( (char*)astr_fileName.c_str());
  nULOUT("\t\t\t\t[ ok ]\n");
  st_env.pMS_active   = st_env.pMS_curvature;

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

  if (!st_env.pMS_curvature)
    return false;

  ULOUT("Mapping curvature texture on primary surface...");
  MRISreadCurvature(st_env.pMS_curvature,
                    (char*)astr_fileName.c_str());
  nULOUT("\t\t\t[ ok ]\n");

  return true;
}

bool
s_env_surfaceSulcal_set(
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

  if (!st_env.pMS_sulcal)
    return false;

  ULOUT("Mapping sulcal texture on primary surface...");
  MRISreadCurvature(st_env.pMS_sulcal,
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

  if (!st_env.pMS_auxSurface)
    return false;

  ULOUT("Mapping curvature texture on auxillary surface...");
  MRISreadCurvature(st_env.pMS_auxSurface,
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

  if (st_env.pMS_auxSurface)
    MRISfree(&st_env.pMS_auxSurface);
  ULOUT("Reading surface for auxillary curvature...");
  st_env.pMS_auxSurface = MRISread( (char*)astr_fileName.c_str());
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
    apst_env->pMS_active = apst_env->pMS_curvature;
    apst_env->esf_active = e_workingCurvature;
    break;
  case 1:
    apst_env->pMS_active = apst_env->pMS_sulcal;
    apst_env->esf_active = e_workingSulcal;
    break;
  case 2:
    apst_env->pMS_active = apst_env->pMS_auxSurface;
    apst_env->esf_active = e_auxillary;
    break;
  default:
    apst_env->pMS_active = apst_env->pMS_curvature;
    apst_env->esf_active = e_workingCurvature;
    break;
  }
}

void
s_env_mpmProgList(
  s_env& ast_env
) {
  int  lw       = ast_env.lw;
  int  rw       = ast_env.rw;

  colprintf(lw, rw, "Current mpmProg index:name" , "[ %d:%s ]\n",
            (int) ast_env.empm_current,
            ast_env.pstr_mpmProgName[ast_env.empm_current].c_str());
  colprintf(lw, rw, "ENV mpmProg use flag:", "[ %d ]\n", ast_env.b_mpmProgUse);
  colprintf(lw, rw, "mpmProg pointer:", "[ %d ]\n", ast_env.pCmpmProg);
  if(ast_env.pCmpmProg) colprintf(lw, rw, "mpmProg initialized:", "[ ok ]\n");
  else colprintf(lw, rw, "mpmProg initialized:", "[ no ]\n");

  lprintf(lw, "\nAll available mpmProgs:-\n");
  for (int i=0; i<ast_env.totalmpmProgs; i++) {
      colprintf(lw, rw, "index:name", "[ %d:%s ]\n",
                i, ast_env.pstr_mpmProgName[i].c_str());
  }
}

int
s_env_mpmProgSetIndex(
    s_env*      apst_env,
    int         aindex
) {
  int   ret     = -1;
  int   lw      = apst_env->lw;
  int   rw      = apst_env->rw;
  switch (aindex) {
    case 0:
        apst_env->b_mpmProgUse      = true;
        apst_env->empm_current      = (e_MPMPROG) aindex;
        if(apst_env->pCmpmProg) {
            lprintf(lw, "Non-NULL mpmProg pointer detected.\n");
            lprintf(lw, "Deleting existing mpmProg '%s'...", apst_env->pstr_mpmProgName[0].c_str());
            delete apst_env->pCmpmProg;
            lprintf(rw, "[ ok ]\n");
        }
        apst_env->pCmpmProg         = new C_mpmProg_autodijk(apst_env);
        break;
    default:
        apst_env->empm_current      = (e_MPMPROG) 0;
        break;
  }
  if(aindex < 0 || aindex >= apst_env->totalmpmProgs)
      ret       = -1;
  else
      ret       = aindex;
  return ret;
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
    //  09 November 2004
    // o Added st_iterInfo
    //

    int           vno_n;
    VERTEX*       v_c;
    VERTEX*       v_n;
    float         dist, ave_curv, curv, max_height, max_curv, cost;
    s_weights*    pSTw = st_env.pSTw;
    MRIS*         surf = st_env.pMS_curvature;

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
        V3_e.f_x = st_env.pMS_curvature->vertices[st_env.endVertex].x;
        V3_e.f_y = st_env.pMS_curvature->vertices[st_env.endVertex].y;
        V3_e.f_z = st_env.pMS_curvature->vertices[st_env.endVertex].z;
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
    pst_iterInfo->f_sulcalHeight = st_env.pMS_sulcal->vertices[vno_c].curv;
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

    max_height  = (st_env.pMS_sulcal->max_curv);
    max_curv    = (st_env.pMS_curvature->max_curv);
    if(st_env.b_useAbsCurvs) {
        f_height        = fabs(f_height);
        curv            = fabs(ave_curv);
    } else {
        f_height        = max_height - st_env.pMS_sulcal->vertices[vno_c].curv;
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

    float       f_cost      = 0.0;
    float       f_distance  = 0.0;
    s_weights*  pSTw        = st_env.pSTw;
    float       wd          = pSTw->wd;
    VERTEX*     v_c         = NULL;
    MRIS*       surf        = st_env.pMS_curvature;
    const char*       pch_proc    = "costFunc_distanceReturn(...)";
    char        pch_txt[65536];
    static bool b_warned    = false;

    v_c         = &surf->vertices[vno_c];
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
    MRIS*       surf        = st_env.pMS_curvature;
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


/* eof */
