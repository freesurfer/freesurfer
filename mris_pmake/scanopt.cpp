/***************************************************************************
                          scanopt.cpp  -  description
                             -------------------
    begin                : Thu Sep 7 2000
    copyright            : (C) 2000 by Rudolph Pienaar
    email                : pienaar@bme.ri.ccf.org
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
//
// NAME
//
//      scanopt.cpp
//
// VERSION
//
// $Id: scanopt.cpp,v 1.1 2009/09/08 22:39:27 nicks Exp $
//
// DESC
//
//      See header file
//
// HISTORY
//      See header file for generic history. Specific methods have their
//      own history log.
//

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iterator>
using namespace std;

#include <scanopt.h>
#include <string.h>

//
//\\\***
// C_scanopt definitions ****>>>>
/////***
//

void
C_scanopt::debug_push(
  string                          astr_currentProc) {
  //
  // ARGS
  //  astr_currentProc        in      method name to
  //                                          "push" on the "stack"
  //
  // DESC
  //  This attempts to keep a simple record of methods that
  //  are called. Note that this "stack" is severely crippled in
  //  that it has no "memory" - names pushed on overwrite those
  //  currently there.
  //

  if (stackDepth_get() >= C_scanopt_STACKDEPTH-1)
    error(  "Out of str_proc stack depth");
  stackDepth_set(stackDepth_get()+1);
  str_proc_set(stackDepth_get(), astr_currentProc);
}

void
C_scanopt::debug_pop() {
  //
  // DESC
  //  "pop" the stack. Since the previous name has been
  //  overwritten, there is no restoration, per se. The
  //  only important parameter really is the stackDepth.
  //

  stackDepth_set(stackDepth_get()-1);
}

void
C_scanopt::error(
  string          astr_msg        /*= "Some error has occured"    */,
  int             code            /*= -1                          */) {
  //
  // ARGS
  //  astr_class              in              name of the current class
  //  atr_msg                 in              message to dump to stderr
  //  code                    in              error code
  //
  // DESC
  //  Print error related information. This routine throws an exception
  //  to the class itself, allowing for coarse grained, but simple
  //  error flagging.
  //

  cerr << "\nFatal error encountered.\n";
  cerr << "\tscanopt object `" << str_name << "' (id: " << id << ")\n";
  cerr << "\tCurrent function: " << str_obj << "::" << str_proc_get() << "\n";
  cerr << "\t" << astr_msg << "\n";
  cerr << "Throwing an exception to (this) with code " << code << "\n\n";
  throw(this);
}

void
C_scanopt::warn(
  string          astr_class      /*= "C_scanopt::"       */,
  string          astr_msg,
  int             code            /*= -1                  */
) {
  //
  // ARGS
  //  astr_class       in              name of the current class
  //  atr_msg          in              message to dump to stderr
  //  code             in              error code
  //
  // DESC
  //  Print error related information. Conceptually identical to
  //  the `error' method, but no expection is thrown.
  //

  cerr << "\nWarning.\n";
  cerr << "\tWorld `" << str_name << "' (id: " << id << ")\n";
  cerr << "\tCurrent function: " << astr_class << str_proc_get() << "\n";
  cerr << "\t" << astr_msg << "(code: " << code << ")\n";
}

void
C_scanopt::function_trace(
  string          astr_class,
  string          astr_msg,
  string          astr_separator) {
  //
  // ARGS
  //  astr_class       in      current class (or derivative)
  //  astr_msg         in      message to print
  //  astr_separator   in      separator char between successive calls
  //
  // DESC
  //  This method allows for run-time debugging-related information
  //  in a particular class and method to be displayed to stderr.
  //

  string str_tab                      = "";
  static string str_objectName        = "";
  static string str_funcName          = "";

  if (verbosity_get() >= stackDepth_get()) {
    cerr << astr_separator;
    for (int i = 0; i < stackDepth_get(); i++)
      str_tab += "\t";
    if (str_objectName != str_name_get() )
      cerr << "\nStochastic World `" << str_name_get();
    cerr << "' (id: " << id_get() << ")\n";
    if (str_funcName != str_proc_get()) {
      cerr << "\n" << str_tab << "Current function: " << astr_class;
      cerr << "::" << str_proc_get() << endl;
      cerr << "\tverbosity = "    << verbosity_get();
      cerr << ", stackDepth = "   << stackDepth_get() << endl;
    }
    cerr << "\n" << str_tab << astr_msg;
  }
  str_objectName      = str_name_get();
  str_funcName        = str_proc_get();
}

void
C_scanopt::core_construct(
  string          astr_name       /*= "unnamed"           */,
  int             a_id            /*= -1                  */,
  int             a_iter          /*= 0                   */,
  int             a_verbosity     /*= 0                   */,
  int             a_warnings      /*= 0                   */,
  int             a_stackDepth    /*= 0                   */,
  string          astr_proc       /*= "noproc"            */
) {
  //
  // ARGS
  //  astr_name        in              name of object
  //  a_id             in              id of object
  //  a_iter           in              current iteration in arbitrary scheme
  //  a_verbosity      in              verbosity of object
  //  a_stackDepth     in              stackDepth
  //  astr_proc        in              current that has been "debug_push"ed
  //
  // DESC
  //  Simply fill in the core values of the object with some defaults
  //
  // HISTORY
  // 07 September 2000
  //  o Initial design and coding
  //

  str_name    = astr_name;
  id          = a_id;
  iter        = a_iter;
  verbosity   = a_verbosity;
  warnings    = a_warnings;
  stackDepth  = a_stackDepth;
  str_proc[stackDepth]    = astr_proc;

  str_obj     = "C_scanopt";

  str_optDes  = "--";
  str_optSep  = "";

}

void
C_scanopt::map_opt_build(
  e_SCANOPT_tokType       e_tokType       /*= e_DesTag*/

) {
  //
  // ARGS
  //  e_tokType       in      type of token, either the
  //                                  command-line
  //                                  --<tok> <val> format    (e_DesTag)
  //                                  or a more conventional
  //                                  <tok> = <val> format    (e_EquLink)
  // DESC
  //  This method creates the fundamental data structure
  //  of the class.
  //
  //  Using the internal pstr_body vector, process this and
  //  create the option map.
  //
  // HISTORY
  // 08 September 2000
  //  o Budded off initial constructor
  //
  // 02 May 2003
  //  o Expanded to accommodate e_EquLink format.
  //

  int         i;
  string      str_argvA;
  string      str_argvB;
  string      str_Equ;
  string      str_opt;
  string      str_val;
  const bool  b_found = false;                // if a search string
  //      is found at the
  //      beginning of a
  //      search, index is
  //      zero, which means
  //      that the bool of
  //      a positive find
  //      is false :-P

  switch (e_tokType) {
  case e_DesTag:
    // Check all the passed strings for couplets (up until the next
    // to last string)
    Lstr_iter   = Lstr_body.begin();
    for (i=0; i<argc_get()-1; i++) {
      str_argvA       = *(Lstr_iter);        // read in a pair of options
      str_argvB       = *(++Lstr_iter);      // that may be a couplet
      //cout << i << "\t" << argc << "\t" << str_argvA << "\t" << str_argvB << endl;
      if (str_argvA.find(str_optDes) == b_found) {
        str_argvA.erase(0, str_optDes.length());
        if (str_argvB.find(str_optDes) == b_found)
          map_opt.insert(pair<string, string>(str_argvA, str_NONCOUPLET));
        else
          map_opt.insert(pair<string, string>(str_argvA, str_argvB));
      }
    }
    // Now we still have the very last appch_argv left, which may be
    // a non-couplet switch, otherwise it is assumed to be noise
    str_argvA   = *(Lstr_iter);
    if (str_argvA.find(str_optDes) == b_found) {
      str_argvA.erase(0, str_optDes.length());
      map_opt.insert(pair<string, string>(str_argvA, str_NONCOUPLET));
    }
    break;
  case e_EquLink:
    Lstr_iter   = Lstr_body.begin();
    for (i=0; i<argc_get()-2; i++) {
      str_argvA       = *(Lstr_iter);         // read in three options
      str_Equ         = *(++Lstr_iter);       // that may be a couplet
      str_argvB       = *(++Lstr_iter);       // linked by an equal sign
      //cout << i << "\t" << argc << "\t" << str_argvA << "\t" << str_Equ << "\t" << str_argvB << endl;
      if (str_Equ.find("=") == b_found) {
        map_opt.insert(pair<string, string>(str_argvA, str_argvB));
      }
      --Lstr_iter;
    }
    break;
  }
}

C_scanopt::C_scanopt(
  int             a_argc,
  char**          appch_argv,
  string          astr_optDes             /*= "--"                */,
  string          astr_optSep             /*= " "                 */
) {
  //
  // ARGS
  //  a_argc           in              number of char* "strings"
  //  appch_argv       in              array of char* "strings"
  //
  // DESC
  //  Constructor
  //  Copies the passed "strings" into the internal str_body
  //  string for future processing.
  //
  // NOTE
  //  The a_argc and appch_argv that are passed to this routine should
  //  *NOT* be the same as the program argc and argv! It is assumed that
  //  appch_argv[0] is the first of the command line options that were
  //  passed to the main program (usually argv[0] contains the program
  //  name and argv[1] is the first option that was passed).
  //
  //  This method will still work, though, even if the external argc and
  //  argv are passed :-)
  //
  // HISTORY
  // 07 September 2000
  //  o Initial design and coding.
  //
  // 08 September 2000
  //  o Fleshing out - budded off map_opt_build method
  //

  core_construct();
  str_optDes_set(     astr_optDes);
  str_optSep_set(     astr_optSep);

  argc_set(   a_argc);
  for (int i=0; i<a_argc; i++)
    Lstr_body.push_back(appch_argv[i]);

  map_opt_build();
}

C_scanopt::C_scanopt(
  string                  astr_filename,
  e_SCANOPT_tokType       e_tokType               /*= e_DesTag    */,
  string                  astr_optDes             /*= "--"        */,
  string                  astr_optSep             /*= " "         */
) {
  //
  // ARGS
  //  astr_filename            in              file containing options
  //  astr_optDes              in              option DESignator
  //  astr_optSep              in              option separator
  //
  // DESC
  //  Constructor
  //  Reads a file into pstr_body payload and build the option map.
  //
  // HISTORY
  // 08 September 2000
  //  o Initial design and coding.
  //

  core_construct();

  debug_push("C_scanopt");
  int         size = 0;
  string      str_word;


  str_optDes_set(     astr_optDes);
  str_optSep_set(     astr_optSep);

  ifstream    istream_optFile(astr_filename.c_str());
  if (!istream_optFile) {
    error("Cannot find input file: " + astr_filename);
  }

  // run through file to find the number of strings
  while (istream_optFile         >> str_word) {
    Lstr_body.push_back(str_word);
    size++;
  }

  argc = size;

  map_opt_build(e_tokType);

  istream_optFile.close();
  debug_pop();
}

C_scanopt::~C_scanopt() {
//    delete [] pstr_body;
}

void
C_scanopt::print() {
  //
  // DESC
  //  Print the current state of the object
  //
  // HISTORY
  // 07 September 2000
  //  o Initial design and coding
  //

  cout << " **** "    << str_obj_get()        << " data dump **** " << endl;
  cout << "\tname\t"  << str_name_get()               << endl;
  cout << "\tid\t"    << id_get()                     << endl;
  cout                                                << endl;
  cout << "\tRaw argument list\n\t\t";
  Lstr_iter   = Lstr_body.begin();
  for (int i=0; i<argc_get(); i++)
    cout << ">" << *(Lstr_iter++) << "< ";
  cout                                                << endl;

  cout << "\tList length\t" << argc                     << endl;

  cout << "\tSTL option map"                          << endl;
  map_iter    = map_opt.begin();
  while (map_iter != map_opt.end()) {
    cout << "\t\t" << map_iter->first;
    cout << "\t\t\t" << map_iter->second              << endl;
    map_iter++;
  }
  cout                                                << endl;
  cout << "\tstr_optDes\t"    << str_optDes_get()     << endl;
  cout << "\tstr_optSep\t"    << str_optSep_get()     << endl;
}

bool
C_scanopt::scanFor(
  string                  astr_target,
  string*                 apstr_value) {
  //
  // ARGS
  //  astr_target              in              map key to search for
  //  apstr_value              out             value of key (if found)
  //
  // RETURN
  //  0        key NOT found
  //  1        key found
  //
  // DESC
  //  The main entry point for methods using this class. After
  //  an object is constructed on a possible set of option pairs,
  //  this method is used to scan for key values.
  //
  //  If a key is found, `true' is returned, and the value
  //  of the key in pstr_value.
  //
  // HISTORY
  // 08 September 2000
  //  o Initial design and coding.
  //

  bool b_ret;

  map_iter    = map_opt.find(astr_target);
  if (map_iter != map_opt.end()) {
    apstr_value->assign(map_iter->second);
    b_ret = true;
  } else
    b_ret = false;
  return b_ret;
}

