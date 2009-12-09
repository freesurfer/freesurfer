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
// $Id: help.cpp,v 1.2 2009/12/09 22:30:02 rudolph Exp $

#include "help.h"


void
synopsis_show(void) {
  //
  // DESC
  //  Show a simple synopsis of program usage.
  //
  // HISTORY
  // 24 September 2004
  //  o Initial design and coding.
  //

  const char*   pch_synopsis = "\n\
 \n\
 \n\
NAME \n\
 \n\
    mris_pmake \n\
 \n\
SYNOPSIS \n\
 \n\
    mris_pmake          [--optionsFile=<fileName>]              \\ \n\
                        [--dir=<workingDir>]                    \\ \n\
                        [--listen | --listenOnPort <port>] \n\
 \n\
DESCRIPTION \n\
 \n\
    'mris_pmake' generates paths on FreeSurfer surfaces based on an edge cost \n\
    and Dijkstra's algorithm. \n\
 \n\
    In its simplest usage, a <start> and <end> vertex index on the surface is \n\
    specified (typically in the <optionsFile>), and the program calculates the \n\
    shortest path connected the points as well as the path cost. \"Shortest\" \n\
    path in this case only has meaning in the context of the cost function that \n\
    is being evaluated (see COST FUNCTION for details). \n\
 \n\
    The program can also be used to calculate regions, specifying a <center> \n\
    vertex and tagging all vertices that satisfy some cost constraint away from \n\
    this vertex. Usually this is used to tag all vertices a certain distance \n\
    away from the <start> vertex. \n\
 \n\
    An interactive mode of operation is also available through a companion \n\
    Python script called 'dsh' that allows asychronous setting of <start> and \n\
    <end> vertices, changes in the cost function weights, etc. This 'dsh' \n\
    script is probably the best and easiest way to run 'mris_pmake'. \n\
 \n\
COST FUNCTION \n\
 \n\
    The cost function is currently a multi-dimensional weight vector of \n\
    following form: \n\
 \n\
       p = w  d +  w c + w h + w  dc + w  dh + w  ch + w   dch + w   (dir) \n\
            d       c     h     dc      dh      ch      dch       dir \n\
 \n\
       where \n\
 \n\
       w_d     : weighting factor for distance, d \n\
       w_c     : weighting factor for curvature, c \n\
       w_h     : weighting factor for sulcal height, h \n\
       w_dc    : weighting factor for product of distance and curve \n\
       w_dh    : weighting factor for product of distance and height \n\
       w_ch    : weighting factor for product of curve and height \n\
       w_dch   : weighting factor for product of distance, curve, and height \n\
       w_dir   : weighting factor for direction \n\
 \n\
    The curvature, c, is specified in the <optionsFile> with the 'curvatureFile' \n\
    option, and the height, h, is specified in the <optionsFile> with the \n\
    'sulcalHeightFile'. These names are somewhat historical, and in theory any \n\
    valid FreeSurfer overlay can be used for 'c' and 'h'. \n\
 \n\
    An additional non-linear penalty is also available, and if \n\
    'b_transitionPenalties' is TRUE, will be applied to the cost function, by \n\
    an index-to-index multiplication of the cost vector. It currently triggered \n\
    if the original 'c' value undergoes a zero-crossing along a trajectory \n\
    path. \n\
 \n\
OPTIONS \n\
 \n\
    --optionsFile=<fileName> \n\
    The main configuration file that specifies the startup run-time behaviour \n\
    of the program, including cost function variables, start and terminal \n\
    vertex indices, cost function, input files, output files, etc. \n\
 \n\
    If the <fileName> contains a directory prefix, this directory will be \n\
    assumed to be the working directory. \n\
 \n\
    Default is 'options.txt' \n\
 \n\
    --dir=<workingDir> \n\
    The working directory. This will override a working directory that might \n\
    have been specified in the <fileName> prefix. \n\
 \n\
    Defaults is current directory. \n\
 \n\
    --listen \n\
    Start in LISTEN mode, i.e. initialize the program and read the default \n\
    'options.txt' file parsing surfaces and curvatures, but do not actually \n\
    calculate a path. Instead, once ready, start listening on the embedded \n\
    server socket for instructions. Send a 'RUN' string in a UDP packet to \n\
    the port specified in <optionsFile> to perform the path search. \n\
 \n\
    --listenOnPort <port> \n\
    Similar to above, but do not interpret the <optionsFile> environment. \n\
    Essentially create the server port on <port> and do nothing else. In this \n\
    mode, the program requires an explicit 'HUP' text string to be sent as \n\
    a UDP packet to <port> to read the default enviroment, or an options file \n\
    can be spec'd by sending a UDP string 'OPT <optionsFile>'. \n\
 \n\
EXAMPLE USE \n\
 \n\
    The best mechanism to run a 'mris_pmake' process is from a companion \n\
    'shell' called 'dsh'. The use of 'dsh' is beyond the scope of this help, \n\
    but in the simplest case (and assuming a valid <optionsFile>), simply \n\
    run 'dsh' and wait for it to provide a prompt. At the prompt type 'RUN' \n\
    and wait for the next prompt, at which simply type 'quit'. \n\
 \n\
    Direct running of 'mris_pmake' can be somewhat cumbersome, since by default \n\
    the process was designed to parse a surface, calculate a cost, and then \n\
    stay resident for additional instructions. These instructions are delivered \n\
    using UDP sockets communication. In this manner, an external 'controller' \n\
    could initialze a 'mris_pmake', read surfaces and curvatures, and then \n\
    RUN different path searches with different vertices and/or modified cost \n\
    function vectors. \n\
 \n\
    'mris_pmake' communicates on three different channels. Most operational \n\
    data is sent to a channel called <userMessages> (in the <optionsFile>). \n\
    System-type messages are sent to a channel called <sysMessages> and results \n\
    are sent to <resultMessages>. If these are defined as files, then the \n\
    channel simply logs to the file. If these are specifed as 'localhost:XXXX' \n\
    then these channels are created as UDP sockets to which 'mris_pmake' \n\
    transmits data. \n\
 \n\
    That all having been said, how do I do a quick-and-dirty run for \n\
    'mris_pmake'? \n\
 \n\
        o Make sure you have a valid <optionsFile> in the working directory \n\
        o run 'mris_pmake' \n\
        o Monitor the <userMessages> channel and wait for: \n\
                \"Listening for socket comms...\" \n\
        o Send a UDP string \"TERM\" to the <controlPort> of 'mris_pmake' \n\
 \n\
    OR, just run 'mris_pmake', wait until a file defined by <costFile> in \n\
    the <optionsFile> appears on the file system (typically a few seconds) \n\
    and then hit <ctrl>C. \n\
 \n\
\n";

   cout << pch_synopsis;
}

void
version_show(void) {
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

string
commandLineOptions_process(
    int         argc,
    char**      ppch_argv,
    s_env&      st_env
) {

    bool        b_optionsFileUse                = true;
    
    string      str_asynchComms                 = "NOP";
    string      str_subjectsDir                 = "";
    string      str_subject                     = "";
    string      str_hemi                        = "";

    string      str_mpmProg                     = "autodijk";
    string      str_mpmArgs                     = "";

    string      str_mainSurfaceFileName         = "inflated";
    string      str_auxSurfaceFileName          = "smoothwm";
    string      str_mainCurvatureFileName       = "smoothwm.H.crv";
    string      str_auxCurvatureFileName        = "sulc";

    try {
        str_subjectsDir = getenv("SUBJECTS_DIR");
    } catch (...) {
        error_exit("processing environment",
                    "SUBJECTS_DIR env variable is not set.", 10);
    }
    while (1) {
        int opt;
        int optionIndex = 0;
        opt = getopt_long(argc, ppch_argv, "", longopts, &optionIndex);
        if ( opt == -1)
        break;

        switch (opt) {
            case 'o':
                st_env.str_optionsFileName.assign(optarg, strlen(optarg));
            break;
            case 'D':
                st_env.str_workingDir.assign(optarg, strlen(optarg));
            break;
            case 'l':
                str_asynchComms         = "LISTEN";
            break;
            case 'L':
                str_asynchComms         = "LISTENPORT";
                st_env.port             = atoi(optarg);
                st_env.timeoutSec       = 60;
            break;
            case '?':
                synopsis_show();
                exit(1);
            break;
            case 'v':
                version_show();
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
            break;
            case 'c':
                str_mainCurvatureFileName       = optarg;
                b_optionsFileUse                = false;
            break;
            case 'd':
                str_auxCurvatureFileName        = optarg;
                b_optionsFileUse                = false;
            break;
            case 'm':
                str_mpmProg                     = optarg;
                st_env.b_mpmProgUse             = true;
                
                b_optionsFileUse                = false;
            break;
            case 'M':
                str_mpmArgs                     = optarg;
                b_optionsFileUse                = false;
            default:
                cout << "?? getopt returned character code " << opt << endl;
        }
    }
    st_env.str_hemi             = str_hemi;
    string str_p                = str_subjectsDir + "/" + str_subject + "/surf/s";
    st_env.b_optionsFileUse     = b_optionsFileUse;
    while(!st_env.b_optionsFileUse) {
        s_env_defaultsSet(st_env);
        st_env.str_mainSurfaceFileName      = str_p + str_mainSurfaceFileName;
        st_env.str_auxSurfaceFileName       = str_p + str_auxSurfaceFileName;
        st_env.str_mainCurvatureFileName    = str_p + str_mainCurvatureFileName;
        st_env.str_auxCurvatureFileName     = str_p + str_auxCurvatureFileName;
        s_env_optionsFile_write(st_env);
        st_env.b_optionsFileUse             = true;
    }
    return str_asynchComms;
}

void    asynchEvent_processHELP(
  s_env&   st_env,
  string   str_event
) {
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
  std::_Ios_Fmtflags    origFlags;

  origFlags    = cout.flags();
  cout.setf(ios::left);

  cout << endl;
  cout.fill('_');
  CWn(rc+lc, "");
  CCn(rc+lc, "Currently understood commands (case-sensitive):");
  CWn(rc+lc, "");
  cout.fill(' ');
  cout.flags(origFlags);
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
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "depends on <OBJECT> -- see below");
  cout.setf(ios::right);
  CW(lc/2, "<verb> \\");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "one of {'get', 'set', 'do', ...");
  CW(lc, "");
  CWn(rc, "\t'saveFrom', 'loadTo'}");
  cout.setf(ios::right);
  CW(lc/2, "[<modifier>]");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "(value depends on <OBJECT><qualifier><verb>");

  cout.fill('_');
  CWn(rc+lc, "");
  cout.fill(' ');
  cout.flags(origFlags);
  cout.setf(ios::right);
  CW(lc/2, "<OBJECT>");
  CW(lc/2, "");
  cout.flags(origFlags);
  cout.setf(ios::left);
  CWn(rc+lc/2, "<qualifier>");
  cout.setf(ios::right);
  CW(lc/2, "WGHT");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'all' 'wd' 'wc' 'wh' 'wdc' 'wdh'");
  CW(lc, "");
  CWn(rc, "\t'wch' 'wdch' 'wdir'}");
  cout.setf(ios::right);
  CW(lc/2, "dWGHT");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'all'  'Dwd' 'Dwc' 'Dwh' 'Dwdc' 'Dwdh'");
  CW(lc, "");
  CWn(rc, "\t'Dwch' 'Dwdch' 'Dwdir'}");
  cout.setf(ios::right);
  CW(lc/2, "VERTEX");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'start' 'end'}");
  cout.setf(ios::right);
  CW(lc/2, "LABEL");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'ply' 'workingSurface' 'auxSurface'}");
  cout.setf(ios::right);
  CW(lc/2, "SURFACE");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'active'}");

  cout.fill('_');
  CWn(rc+lc, "");
  cout.fill(' ');
  cout.flags(origFlags);
  cout.setf(ios::right);
  CW(lc/2, "<qualifier>");
  CW(lc/2, "");
  cout.flags(origFlags);
  cout.setf(ios::left);
  CWn(rc+lc/2, "<verb> (object specific)");
  cout.setf(ios::right);
  CW(lc/2, "all");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'get' 'set <value>'}");
  cout.setf(ios::right);
  CW(lc/2, "start");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'get' 'set <value>'}");
  cout.setf(ios::right);
  CW(lc/2, "end");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'get' 'set <value>'}");
  cout.setf(ios::right);
  CW(lc/2, "ply");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'get' 'set <depth>' 'do' 'save <prefix>'");
  cout.setf(ios::right);
  CW(lc/2, "active");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'get' 'set aux|working' ");
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "'functionList' 'functionAssign <index>'}");
  cout.setf(ios::right);
  CW(lc/2, "workingSurface");
  cout.flags(origFlags);
  CW(lc/2, "");
  cout.setf(ios::left);
  CWn(rc, "{'loadTo <fileToLoadToSurface>'");
  CW(lc, "");
  CWn(rc, "'saveFrom <fileToContainSurface>'}");
  CW(lc, "");
  CWn(rc, "'singleVertexSet <vertexIndex>'}");
  cout.setf(ios::right);
  CW(lc/2, "auxSurface");
  cout.flags(origFlags);
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

  cout.flags(origFlags);
}


/* eof */
