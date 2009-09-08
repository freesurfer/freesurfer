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
// $Id: help.cpp,v 1.1 2009/09/08 22:39:27 nicks Exp $

#include "help.h"

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
