/**
 * @brief help text utils
 *
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

#include "help.h"
#include "utils.h"
#include "mris_pmake.help.xml.h"

void
synopsis_show(void)
{
  outputHelpXml(mris_pmake_help_xml,
                mris_pmake_help_xml_len);
}

void
version_show(void)
{
  //
  // DESC
  //  Show program version.
  //
  // HISTORY
  // 21 October 2004
  //  o Initial design and coding.
  //

  cout << endl << "\t\t" << G_VERSION;
  cout << endl << "";
  exit(0);
}

bool
mpmProg_check(
  s_env&    st_env,
  string    str_mpmProg
)
{

  bool b_validMpmProg = false;

  if(!str_mpmProg.length())
  {
    s_env_mpmPrint
      (st_env,
       "\nNo mpmProg specified! You must specify an mpmProg to use.\n",
       e_mpmProg);
    error_exit("processing command line options,",
               "you must also specify '--mpmProg <mpmProg>'.",
               20);
  }
  if(str_mpmProg == "NULL")
  {
    b_validMpmProg        	= true;
    st_env.empmProg_current 	= emp_NULL;
  }
  if(str_mpmProg == "NOP")
  {
    b_validMpmProg        	= true;
    st_env.empmProg_current 	= emp_NOP;
  }
  if(str_mpmProg == "pathFind")
  {
    b_validMpmProg        	= true;
    st_env.empmProg_current 	= emp_pathFind;
  }
  if(str_mpmProg == "autodijk")
  {
    b_validMpmProg        	= true;
    st_env.empmProg_current 	= emp_autodijk;
  }
  if(str_mpmProg == "autodijk_fast")
  {
    b_validMpmProg        	= true;
    st_env.empmProg_current 	= emp_autodijk_fast;
  }
  if(str_mpmProg == "ROI")
  {
    b_validMpmProg              = true;
    st_env.empmProg_current     = emp_ROI;
  }
  if(str_mpmProg == "externalMesh")
  {
    b_validMpmProg              = true;
    st_env.empmProg_current     = emp_externalMesh;
  }
  return b_validMpmProg;
}

bool
mpmOverlay_check(
  s_env&    st_env,
  string    str_mpmOverlay
)
{

  bool b_validMpmOverlay  = false;

  if(!str_mpmOverlay.length())
  {
    s_env_mpmPrint
      (st_env,
       "\nNo overlay specified! You must specify an mpmOverlay to use.\n",
       e_mpmOverlay);
    error_exit("processing command line options,",
               "you must also specify '--mpmOverlay <mpmOverlayID>'.",
               20);
  }
  if(str_mpmOverlay == "legacy")
  {
    b_validMpmOverlay           = true;
    st_env.empmOverlay_current  = emo_LEGACY;
  }
  if(str_mpmOverlay == "NULL")
  {
    b_validMpmOverlay           = true;
    st_env.empmOverlay_current  = emo_NULL;
  }
  if(str_mpmOverlay == "NOP")
  {
    b_validMpmOverlay           = true;
    st_env.empmOverlay_current  = emo_NOP;
  }
  if(str_mpmOverlay == "unity")
  {
    b_validMpmOverlay           = true;
    st_env.empmOverlay_current  = emo_unity;
  }
  if(str_mpmOverlay == "distance")
  {
    b_validMpmOverlay           = true;
    st_env.empmOverlay_current  = emo_distance;
  }
  if(str_mpmOverlay == "euclidean")
  {
    b_validMpmOverlay           = true;
    st_env.empmOverlay_current  = emo_euclidean;
  }
  if(str_mpmOverlay == "fscurvs")
  {
    b_validMpmOverlay           = true;
    st_env.empmOverlay_current  = emo_fscurvs;
  }
  if(str_mpmOverlay == "curvature")
  {
    b_validMpmOverlay           = true;
    st_env.empmOverlay_current  = emo_curvature;
  }

  return b_validMpmOverlay;
}

string
commandLineOptions_process(
  int         argc,
  char**      ppch_argv,
  s_env&      st_env
)
{

    //
    // PRECONDITIONS
    // o An 'options.txt' file if no command line arguments are
    //	 specified.
    //
    // POSTCONDITIONS
    // o If any "trigger" command line argument is passed, the 'options.txt'
    //   file is regenerated (and overwritten if it already exists). 
    //
    
    
  bool        b_optionsFileUse                = true;
  bool        b_useAbsCurvs                   = false;

  string      str_asynchComms                 = "HUP";
  string      str_subjectsDir                 = "";
  char*       pch_subjectsDir;
  string      str_subject                     = "";
  string      str_hemi                        = "";
  int	      port                            = 1701;

  string      str_mpmProg                     = "";
  string      str_mpmArgs                     = "-x";
  string      str_mpmOverlay         	      = "";
  string      str_mpmOverlayArgs              = "-x";

  string      str_mainSurfaceFileName         = "inflated";
  string      str_auxSurfaceFileName          = "smoothwm";
  string      str_mainCurvatureFileName       = "smoothwm.H.crv";
  string      str_auxCurvatureFileName        = "sulc";

  st_env.b_primaryCurvature                   = false;
  st_env.b_secondarySurface                   = false;
  st_env.b_secondaryCurvature                 = false;


  if( (pch_subjectsDir = getenv("SUBJECTS_DIR")) == NULL)
    error_exit("processing environment,",
               "it seems that the SUBJECTS_DIR env variable is not set.", 10);
  str_subjectsDir       = pch_subjectsDir;
  st_env.argc           = argc;
  st_env.ppch_argv      = ppch_argv;
  while (1)
  {
    int opt;
    int optionIndex = 0;
    opt = getopt_long(argc, ppch_argv, "", longopts, &optionIndex);
    if ( opt == -1)
    {
      break;
    }

    switch (opt)
    {
    case 'o':
      st_env.str_optionsFileName.assign(optarg, strlen(optarg));
      break;
    case 'D':
      st_env.str_workingDir.assign(optarg, strlen(optarg));
      break;
    case '?':
    case 'u':
      synopsis_show();
      exit(1);
      break;
    case 'v':
      version_show();
      break;
	    
    // All the following options trigger a regeneration of the 'options.txt'
    // file. If it already exsts, it will be overwritten.
    case 'a':
      b_useAbsCurvs           	      = true;
      b_optionsFileUse        	      = false;
      break;
    case 'h':
      str_hemi                        = optarg;
      b_optionsFileUse                = false;
      break;
    case 'S':
      str_subject                     = optarg;
      b_optionsFileUse                = false;
      break;
    case 's':
      str_mainSurfaceFileName         = optarg;
      b_optionsFileUse                = false;
      break;
    case 't':
      str_auxSurfaceFileName          = optarg;
      b_optionsFileUse                = false;
      st_env.b_secondarySurface       = true;
      break;
    case 'c':
      str_mainCurvatureFileName       = optarg;
      b_optionsFileUse                = false;
      st_env.b_primaryCurvature       = true;
      break;
    case 'd':
      str_auxCurvatureFileName        = optarg;
      b_optionsFileUse                = false;
      st_env.b_secondaryCurvature     = true;
      break;
    case 'm':
      str_mpmProg                     = optarg;
      st_env.b_mpmProgUse             = true;
      b_optionsFileUse                = false;
      break;
    case 'M':
      str_mpmArgs                     = optarg;
      b_optionsFileUse                = false;
      break;
    case 'O':
      str_mpmOverlay                    = optarg;
      st_env.b_mpmOverlayUse            = true;
      b_optionsFileUse                  = false;
      break;
    case 'p':
      port                              = atoi(optarg);
      b_optionsFileUse                  = false;
      break;
    case 'V':
      str_mpmOverlayArgs                = optarg;
      b_optionsFileUse                  = false;
      break;
    default:
      cout << "?? getopt returned character code " << opt << endl;
      break;
    }
  }
  st_env.str_hemi             = str_hemi;
  string str_p                = str_subjectsDir + "/" + str_subject + "/surf/";
  st_env.b_optionsFileUse     = b_optionsFileUse;
  while(!st_env.b_optionsFileUse)
  {
    s_env_defaultsSet(st_env);
    st_env.serverControlPort            = port;
    st_env.b_useAbsCurvs                = b_useAbsCurvs;
    st_env.str_primarySurfaceFileName   = str_p + str_hemi + "." +
            str_mainSurfaceFileName;
    st_env.str_secondarySurfaceFileName = str_p + str_hemi + "." +
            str_auxSurfaceFileName;
    st_env.str_primaryCurvatureFileName = str_p + str_hemi + "." +
            str_mainCurvatureFileName;
    st_env.str_secondaryCurvatureFileName = str_p+str_hemi + "." +
            str_auxCurvatureFileName;
    if(!mpmProg_check (st_env, str_mpmProg))
    {
      s_env_mpmPrint(st_env,
                     "\nInvalid mpmProg! You must specify a correct mpmProg to use.\n",
                     e_mpmProg);
      error_exit("processing command line options,",
                 "I didn't recognize the <mpmProg>.",
                 20);
    }
    st_env.b_mpmProgUse                     = true;
    st_env.str_mpmArgs                      = str_mpmArgs;
    st_env.b_optionsFileUse                 = true;
    st_env.b_exitOnDone                     = true;

    if(!mpmOverlay_check (st_env, str_mpmOverlay))
    {
      s_env_mpmPrint(st_env,
                     "\nInvalid mpmOverlay! You must specify an mpmOverlay to use.\n",
                     e_mpmOverlay);
      error_exit("processing command line options,",
                 "I didn't recognize the <mpmOverlay>.",
                 20);
    }
    st_env.str_mpmOverlayArgs               = str_mpmOverlayArgs;
    s_env_optionsFile_write(st_env);

    str_asynchComms                         = "HUP";
  }
  return str_asynchComms;
}

#if 1
void    asynchEvent_processHELP(
  s_env&   st_env,
  string   str_event
)
{
  //
  // PRECONDITIONS
  //    o Valid environment
  //   o The str_event has its "primary" event string removed.
  //
  // POSTCONDITIONS
  //    o Help information is printed on stdout.
  //
  // HISTORY
  // 15 March 2005
  //    o Initial design and coding.
  //

  int lc         = 30;
  int rc         = 50;
  //std::_Ios_Fmtflags    origFlags;

  //origFlags    = cout.flags();
  cout.setf(ios::left);

  cout << endl;
  cout.fill('_');
  CWn(rc+lc, "");
  CCn(rc+lc, "Currently understood commands (case-sensitive):");
  CWn(rc+lc, "");
  cout.fill(' ');
  //  cout.flags(origFlags);
  cout.fill('_');
  CWn(rc+lc, "");
  CCn(rc+lc, "No argument commands:");
  CWn(rc+lc, "");
  cout.fill(' ');
  cout.setf(ios::left);
  CW(lc, "HUP");
  CWn(rc, "Re-read options file and re-process.");
  CW(lc, "RUN");
  CWn(rc, "Re-run with current weight values.");
  CW(lc, "HELP");
  CWn(rc, "This message.");
  CW(lc, "TERM");
  CWn(rc, "Terminate and exit to system.");
  CW(lc, "LISTEN");
  CWn(rc, "Listen on server port *and* re-read");
  CW(lc, "");
  CWn(rc, "\toptions file, i.e. create");
  CW(lc, "");
  CWn(rc, "\tenvironment.");
  CW(lc, "LISTENPORT");
  CWn(rc, "Listen on server port, but do *not*");
  CW(lc, "");
  CWn(rc, "\tread options file, or create");
  CW(lc, "");
  CWn(rc, "\tenvironment.");
  cout.fill('_');
  CWn(rc+lc, "");
  CCn(rc+lc, "Single argument commands:");
  CWn(rc+lc, "");
  cout.fill(' ');
  cout.setf(ios::left);
  CW(lc, "SYSECHO <str>");
  CWn(rc, "Write <str> to sys_msg channel.");
  CW(lc, "USERECHO <str>");
  CWn(rc, "Write <str> to user_msg channel.");
  CW(lc, "CWD <path>");
  CWn(rc, "Change working directory to <path>.");
  CW(lc, "OPT <file>");
  CWn(rc, "Change options file to <file>. If a");
  CW(lc, "");
  CWn(rc, "\tpath is prefixed, the working");
  CW(lc, "");
  CWn(rc, "\tdirectory is also changed.");
  cout.fill('_');
  CWn(rc+lc, "");
  CCn(rc+lc, "Multiple argument commands:");
  CWn(rc+lc, "");
  cout.fill(' ');
  cout.setf(ios::left);
  CW(lc,  "<OBJECT> \\");
  CWn(rc, "Perform operations on <OBJECT>.");
  cout.setf(ios::right);
  CW(lc/2, "<qualifier> \\");
  //cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "depends on <OBJECT> -- see below");
  cout.setf(ios::right);
  CW(lc/2, "<verb> \\");
  //cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "one of {'get', 'set', 'do', ...");
  CW(lc, "");
  CWn(rc, "\t'saveFrom', 'loadTo'}");
  cout.setf(ios::right);
  CW(lc/2, "[<modifier>]");
  //cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "(value depends on <OBJECT><qualifier><verb>");

  cout.fill('_');
  CWn(rc+lc, "");
  cout.fill(' ');
  //cout.flags(origFlags);
  cout.setf(ios::right);
  CW(lc/2, "<OBJECT>");
  CW(lc/2, "");
  //cout.flags(origFlags);
  cout.setf(ios::left);
  CWn(rc+lc/2, "<qualifier>");
  cout.setf(ios::right);
  CW(lc/2, "WGHT");
  //cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'all' 'wd' 'wc' 'wh' 'wdc' 'wdh'");
  CW(lc, "");
  CWn(rc, "\t'wch' 'wdch' 'wdir'}");
  cout.setf(ios::right);
  CW(lc/2, "dWGHT");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'all'  'Dwd' 'Dwc' 'Dwh' 'Dwdc' 'Dwdh'");
  CW(lc, "");
  CWn(rc, "\t'Dwch' 'Dwdch' 'Dwdir'}");
  cout.setf(ios::right);
  CW(lc/2, "VERTEX");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'start' 'end'}");
  cout.setf(ios::right);
  CW(lc/2, "LABEL");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'ply' 'workingSurface' 'auxSurface'}");
  cout.setf(ios::right);
  CW(lc/2, "SURFACE");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'active'}");

  cout.fill('_');
  CWn(rc+lc, "");
  cout.fill(' ');
  //  cout.flags(origFlags);
  cout.setf(ios::right);
  CW(lc/2, "<qualifier>");
  CW(lc/2, "");
  //  cout.flags(origFlags);
  cout.setf(ios::left);
  CWn(rc+lc/2, "<verb> (object specific)");
  cout.setf(ios::right);
  CW(lc/2, "all");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'get' 'set <value>'}");
  cout.setf(ios::right);
  CW(lc/2, "start");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'get' 'set <value>'}");
  cout.setf(ios::right);
  CW(lc/2, "end");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'get' 'set <value>'}");
  cout.setf(ios::right);
  CW(lc/2, "ply");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'get' 'set <depth>' 'do' 'save <prefix>'");
  cout.setf(ios::right);
  CW(lc/2, "active");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'get' 'set aux|working' ");
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "'functionList' 'functionAssign <index>'}");
  cout.setf(ios::right);
  CW(lc/2, "workingSurface");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'loadTo <fileToLoadToSurface>'");
  CW(lc, "");
  CWn(rc, "'saveFrom <fileToContainSurface>'}");
  CW(lc, "");
  CWn(rc, "'singleVertexSet <vertexIndex>'}");
  cout.setf(ios::right);
  CW(lc/2, "auxSurface");
  //  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'loadTo <fileToLoadToSurface>'");
  CW(lc, "");
  CWn(rc, "'saveFrom <fileToContainSurface>'}");
  CW(lc, "");
  CWn(rc, "'singleVertexSet <vertexIndex>'}");

  cout.fill('_');
  CWn(rc+lc, "");
  CCn(rc+lc, "Examples:");
  CWn(rc+lc, "");
  cout.fill(' ');
  cout.setf(ios::left);
  CW(lc, "[D]WGHT all get");
  CWn(rc,"Print all [delta] polynomial weights.");
  CW(lc, "[D]WGHT all set 0.5");
  CWn(rc,"Set all [delta] polynomial weights to 0.5.");
  cout << endl;
  CW(lc, "VERTEX start get");
  CWn(rc,"Print the value of the start VERTEX.");
  CW(lc, "VERTEX end set 96650");
  CWn(rc,"Set the end VERTEX to 96650.");
  cout << endl;
  CW(lc, "LABEL auxSurface \\");
  CWn(rc,"Load the label file 'dijk.label' to internal core");
  CW(lc, "    loadFrom dijk.label");
  CWn(rc, "and project onto the orig surface.")
  CW(lc, "LABEL workingSurface \\");
  CWn(rc,"Project the internal label trajectory to the");
  CW(lc, "    saveTo dijk.label");
  CWn(rc,"working surface and save to 'dijk.label'.");
  CW(lc, "LABEL workingSurface \\");
  CWn(rc,"Clear any existing labels and set vertex");
  CW(lc, "    singleVertexSet 78832");
  CWn(rc,"78832 as a 'label' for further processing.");
  CW(lc, "LABEL ply set 6");
  CWn(rc, "Set the ply depth around the");
  CW(lc, "");
  CWn(rc, "core label to 6");
  CW(lc, "LABEL ply functionList");
  CWn(rc, "An indexed list of available");
  CW(lc, "");
  CWn(rc, "ply functions");
  CW(lc, "LABEL ply functionAssign 2");
  CWn(rc, "Set the ply function index to");
  CW(lc, "");
  CWn(rc, "'2', i.e. unity cost function");
  CW(lc, "LABEL ply do");
  CWn(rc, "Calculate ply depth around the");
  CW(lc, "");
  CWn(rc, "core label.");
  CW(lc, "LABEL ply surfaceSet \\");
  CWn(rc, "Associate the internal ply labels");
  CW(lc, "    'aux'|'working'");
  CWn(rc,"to either the 'orig' or 'working' surface.");
  CW(lc, "LABEL ply saveStaggered <prefix>");
  CWn(rc, "Save the ply labels to a series");
  CW(lc, "");
  CWn(rc, "of file names starting with <prefix> and");
  CW(lc, "");
  CWn(rc, "ending with '.ply.N.label' where");
  CW(lc, "");
  CWn(rc, "N = 1... <plyDepth>.");
  CW(lc, "LABEL ply save <prefix>");
  CWn(rc, "Save the greatest extent ply label");
  CW(lc, "");
  CWn(rc, "to a file name starting with <prefix> and");
  CW(lc, "");
  CWn(rc, "ending with '.ply.N.label' where");
  CW(lc, "");
  CWn(rc, "N = <plyDepth>.");
  CW(lc, "SURFACE active set aux");
  CWn(rc, "Set the 'active' surface");
  CW(lc, "");
  CWn(rc, "to the auxillary surface");
  cout << endl;

  //cout.flags(origFlags);
}
#endif

/* eof */
