/**
 * @brief This class does option processing
 *
 *
 *   `scanopt.h' provides the header definition for the C_scanopt class.
 *   This class does option processing, similar (at least conceptually)
 *   to the getopt family of functions bundled with the C-library.
 *
 *   Although the original design of the class was to process command
 *   line arguments, it has been used increasingly to parse "option"
 *   files, extracting specific tokens and their corresponding values.
 *
 *   The class is surprising adept at this parsing, mostly due to its dogged
 *   "single mindedness" (i.e. stupidity). It understands only the most
 *   basic syntax, caring little for structure or semantics of the files
 *   it parses. Everything that does not conform to its basic syntax is
 *   simply ignored as noise.
 *
 *   Token:value pairs can be stipulated in one of two ways. The first
 *   harkens back to the class's origin as a command line parser and
 *   consists of the following syntax:
 *
 *           [tokenPrefix][token]<whitespace>[value]
 *
 *   where [tokenPrefix] is by default '--' (however I have often used
 *   the ":" character as well). The following example illustrates
 *   this syntax:
 *
 *           --colourValue1          Red
 *           --colourValue2          Blue
 *
 *   Note that ANY (7-bit ASCII )token and ANY (7-bit ASCII) value pair
 *   can be specified with this syntax. The class triggers on the presence
 *   of the "--" tokenPrefix.
 *
 *   The second token:value paring can be presented as:
 *
 *           [token]<[whitespace]>=<[whitespace]>[value]
 *
 *   for example
 *
 *           colourValue1    =       Red
 *           colourValue2    =       Blue
 *
 *   Here, the presence of the equals character "=" links a token:value
 *   pair. This link character defaults to "=" but can be set to 
 *   meaningful character.
 *
 *   If constructed with a filename argument (and assuming a valid file),
 *   the class creates a bi-directional linked list of strings containing
 *   the file's contents. A linked list is used simply because it allows
 *   for an easy means of dynamically recording a file's contents
 *   internally.
 *
 *   This list is then parsed to construct a STL-type
 *   map<string, string> which can be rapidly searched for target tokens
 *   and their values.
 *
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iterator>
#include <algorithm>
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
  string          astr_msg        /*= "Some error has occurred"    */,
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
      cerr << "\nC_scanopt `" << str_name_get();
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
  str_equ     = "=";

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
                                              //+ is found at the
                                              //+ beginning of a
                                              //+ search, index is
                                              //+ zero, which means
                                              //+ that the bool of
                                              //+ a positive find
                                              //+ is false :-P

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
      map_opt.insert(pair<string, string>(str_trim(str_argvA), str_NONCOUPLET));
    }
    break;
  case e_EquLink:
    Lstr_iter   = Lstr_body.begin();
    for (i=0; i<argc_get()-2; i++) {
      str_argvA       = *(Lstr_iter);         // read in three options
      str_Equ         = *(++Lstr_iter);       // that may be a couplet
      str_argvB       = *(++Lstr_iter);       // linked by an equal sign
      //cout << i << "\t" << argc << "\t" << str_argvA << "\t" << str_Equ << "\t" << str_argvB << endl;
      if (str_Equ.find(str_equ) == b_found) {
	// Strip the str_optDes tag from the first string (if found)
        if (str_argvA.find(str_optDes) == b_found) 
	    str_argvA.erase(0, str_optDes.length());
        map_opt.insert(pair<string, string>(str_trim(str_argvA), str_trim(str_argvB)));
      }
      --Lstr_iter;
    }
    break;
  }
}

C_scanopt::C_scanopt(
  int             	a_argc,
  char**          	appch_argv,
  string          	astr_optDes             /*= "--"                */
) {
  //
  // ARGS
  //  a_argc            in		numbe of char* "strings"
  //  appch_argv        in              array of char* "strings"
  //  astr_optDes	in		option designator string
  //
  // DESC
  //  Constructor
  //  Copies the passed "strings" into the internal str_body
  //  string for future processing.
  //
  // NOTE
  //  This method can process command line options in the e_DesTag
  //  format.
  //
  // HISTORY
  // 07 September 2000
  //  o Initial design and coding.
  //
  // 08 September 2000
  //  o Fleshing out - budded off map_opt_build method
  //
  // 22 February 2011
  //  o Explicit limiting to e_DesTag type specification.
  //

    core_construct();
    str_optDes_set(     astr_optDes);

    argc_set(   a_argc);
    for (int i=0; i<a_argc; i++)
        Lstr_body.push_back(appch_argv[i]);

    map_opt_build(e_DesTag);
}

C_scanopt::C_scanopt(
  string                  astr_filename,
  e_SCANOPT_tokType       e_tokType               /*= e_DesTag    */,
  string                  astr_optDes             /*= "--"        */,
  string                  astr_equ                /*= "="         */
) {
  //
  // ARGS
  //  astr_filename             in              file containing options
  //  e_tokType                 in              token link type
  //  astr_optDes               in              option DESignator
  //  astr_equ                  in              the "equal" string for 
  //                                            e_EquLink
  //
  // DESC
  //  Constructor
  //  Reads a file into pstr_body payload and build the option map.
  //
  // HISTORY
  // 08 September 2000
  //  o Initial design and coding.
  //  
  // 18 February 2011
  //  o Cleanup.
  //

  core_construct();

  debug_push("C_scanopt");
  int         size = 0;
  string      str_word;


  str_optDes_set(     astr_optDes);
  str_equ_set(        astr_equ);

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

int
str_tokenize(
    const string&       str,
    vector<string>&     tokens,
    const string&       delimiters)
{
    int tokenCount      = 0;
    
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
        tokenCount++;
    }
    return tokenCount;
}

string 
str_trim(
    string 			astr_toTrim) {
 
    size_t found;
    found = astr_toTrim.find(" ");
 
    while(astr_toTrim.length() != 0 && found != string::npos) {
	astr_toTrim.erase(found,1);		
	found = astr_toTrim.find(" ");
    }
    return astr_toTrim;
}

C_scanopt::C_scanopt(
    string                      astr_options,
    string                      astr_delimiter          /*= ";"         */,
    e_SCANOPT_tokType           e_tokType               /*= e_EquLink   */,
    string                      astr_optDes             /*= "--"        */,
    string                      astr_equ                /*= "="         */
) {
    //
    // ARGS
    //  astr_options            in                      string with all options
    //  astr_delimiter          in                      delimited string
    //  e_tokType               in                      token type
    //  astr_optDes             in                      option DESignator
    //  astr_equ                in                      "equal" string for 
    //                                                  e_EqulLink.
    //
    // DESC
    //  Constructor
    //  Parses the <astr_options> into tokens delimited by <astr_delimiter>, and
    //  builds internal map.
    //
    // PRECONDITIONS
    // o Valid <tag> <value> pairs in astr_options.
    //
    // POSTCONDITIONS
    // o <tag> <value> pairs are pushed to internal map:
    // 	+ if <astr_optDes> is specified and non-zero, only <tag>s that start
    //	  with <astr_optDes> are pushed.
    //  + if e_EquLink is specified, then the astr_equ is also pushed:
    //	   	<astr_optDes><tag> = <val>
    //
    // HISTORY
    // 21 December 2009
    //  o Initial design and coding.
    //
    // February 2011
    // 	o Updated / corrected
    //

    core_construct();

    debug_push("C_scanopt");
    int                 pairLines       = 0;
    int                 size            = 0;
    vector<string>      v_lines;
    vector<string>      v_option;
    string 		str_token("");

    str_optDes_set(     astr_optDes);
    str_equ_set(        astr_equ);

    v_lines.clear();
    if((pairLines=str_tokenize(astr_options, v_lines, astr_delimiter))!=0) {
        // First loop over each "line" of options... 
        //+ each line is separated by <astr_delimiter>
        for(vector<string>::iterator i = v_lines.begin();
            i != v_lines.end(); i++) {
            // Now tokenizing each line on the <astr_eq>
            // cout << "Tokenizing " << *i << endl;
            v_option.clear();
	    str_token = e_tokType == e_EquLink ? astr_equ : " ";
  	    if( str_tokenize(*i, v_option, str_token) ) {
                // Push the option and its value into the internal
                // Lstr_body...
		bool b_canPush = true;
		//printf("length optDes: %d\n", (int)astr_optDes.length());
		if(astr_optDes.length()) {
		    if(v_option.begin()->find(astr_optDes)!=0)
		        b_canPush = false;
		}
		if(b_canPush) {
                    for(vector<string>::iterator j = v_option.begin();
                        j != v_option.end(); j++) {
                        // cout << "\tPushing " << *j << endl;
                        Lstr_body.push_back(*j);
                        if(j==v_option.begin() && e_tokType == e_EquLink) {
			    // Need to explicitly push the astr_equ into the
			    // body for the e_EquLink case.
                            size++;
			    Lstr_body.push_back(astr_equ);
                        }
                        size++;
                    }
                }
	    }
        }
    }

    argc = size;

    map_opt_build(e_tokType);

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

  cout << "\tList length\t" << argc                   << endl;

  cout << "\tSTL option map"                          << endl;
  map_iter    = map_opt.begin();
  while (map_iter != map_opt.end()) {
    cout << "\t\t" << map_iter->first;
    cout << "\t\t\t" << map_iter->second              << endl;
    map_iter++;
  }
  cout                                                << endl;
  cout << "\tstr_optDes\t"    << "->" << str_optDes_get() << "<-" << endl;
  cout << "\tstr_equ\t\t"     << "->" << str_equ_get()    << "<-" << endl;
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

