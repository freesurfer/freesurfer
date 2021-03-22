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
#include <string>
#include <map>
#include <list>

#include <vector>

using namespace std;

#ifndef __SCANOPT_H__
#define __SCANOPT_H__

const string    str_NONCOUPLET  = "< nonCouplet >";
const int       C_scanopt_STACKDEPTH      = 64;

typedef enum {
  e_DesTag,               // --<token>    <val>
  e_EquLink               // <token>  =   <val>
}
e_SCANOPT_tokType;

class C_scanopt {

  // data structures

protected:
  //
  // generic object structures - used for internal bookkeeping
  // and debugging / automated tracing methods. The stackDepth
  // and str_proc[] variables are maintained by the debug_push|pop
  // methods
  //
  string  str_obj;                  // name of object class
  string  str_name;                 // name of object variable
  int     id;                       // id of agent
  int     iter;                     // current iteration in an
                                    //+ arbitrary processing scheme
  int     verbosity;                // debug related value for object
  int     warnings;                 // show warnings (and warnings level)
  int     stackDepth;               // current pseudo stack depth
  string  str_proc[C_scanopt_STACKDEPTH];   // execution procedure stack


protected:
  int                             argc;             // size of string vector
  list<string>                    Lstr_body;        // string list containing
                                                    //+ options and values
  list<string>::iterator          Lstr_iter;        // an iterator to move over 
                                                    //+ the list

  map<string, string>             map_opt;          // a map associating 
                                                    //+ possible optName and 
                                                    //+ optVal pairs
  map<string, string>::iterator   map_iter;         // an iterator to move over the map

  string                          str_optDes;       // the "pattern" that defines
                                                    //+ a string as an option.
                                                    //+ Typically this will be
                                                    //+ either "-" or "--".
  string                          str_equ;          // The string for the
                                                    //+ EquLink. Defaults to
                                                    //+ "="
  


    // methods

public:
  //
  // constructor / destructor block
  //

  // conventional constructor using passed command line args
  C_scanopt(      int                     a_argc,
                  char**                  appch_argv,
                  string                  astr_optDes     = "--");

  // constructor reading options from file
  C_scanopt(      string                  astr_filename,
                  e_SCANOPT_tokType       e_tokType       = e_DesTag,
                  string                  astr_optDes     = "--",
                  string                  astr_equ        = "=");

  C_scanopt(      string                  astr_options,
                  string                  astr_delimiter  = ";",
                  e_SCANOPT_tokType       e_tokType       = e_EquLink,
                  string                  astr_optDes     = "--",
                  string                  astr_equ        = "=");
                  
  void    core_construct( string astr_name    = "unnamed",
                          int a_id            = -1,
                          int a_iter          = 0,
                          int a_verbosity     = 0,
                          int a_warnings      = 0,
                          int a_stackDepth    = 0,
                          string astr_proc    = "noproc");
  ~C_scanopt();

  //
  // error / warn / print block
  //
  void        debug_push(         string astr_currentProc);
  void        debug_pop();

  void        error(              string  astr_msg        = "Some error has occured",
                                  int     code            = -1);
  void        warn(               string astr_class = "C_scanopt::",
                                  string astr_msg = "", int code = -1);
  void        function_trace(     string astr_class,
                                  string astr_msg,
                                  string astr_separator);
  void        function_trace(     string astr_msg) {
    function_trace("C_scanopt", astr_msg, "");
  };

    //
    // access block
    //
    void    print();                        // print object
    int     stackDepth_get()const {return stackDepth;};
    void    stackDepth_set(int anum) {stackDepth = anum;};
    int     iter_get()      const   {return iter;};
    void    iter_set(int anum)      {iter = anum;};
    int     id_get()        const   {return id;};
    void    id_set(int anum)        {id = anum;};
    int     verbosity_get() const   {return verbosity;};
    void    verbosity_set(int anum) {verbosity = anum;};
    int     warnings_get()  const   {return warnings;};
    void    warnings_set(int anum)  {warnings = anum;};
    string  str_obj_get()   const   {return str_obj;};
    void    str_obj_set(string astr) {str_obj = astr;};
    string  str_name_get()  const   {return str_name;};
    void    str_name_set(string astr) {str_name = astr;};
    string  str_proc_get()  const   {return str_proc[stackDepth_get()];};
    void    str_proc_set(int depth, string astr) {str_proc[depth] = astr;};
    int     argc_get()      const   {return argc;};
    void    argc_set(int anum)      {argc = anum;};
    string  str_optDes_get() const  {return str_optDes;};
    void    str_optDes_set(string astr) {str_optDes = astr;};
    string  str_equ_get() const  {return str_equ;};
    void    str_equ_set(string astr) {str_equ = astr;};

  //
  // miscellaneous block
  //
private:
  void    map_opt_build(
    e_SCANOPT_tokType e_tokType = e_DesTag
  );                // build option map
protected:
  void    nullify();                      // "reset" structure to 0
public:
  bool    scanFor(string  astr_target,    // the main purpose of the class
                  string* apstr_value);
};

string 	str_trim(string		astr_toTrim);
int     str_tokenize(
            const string&       str,
            vector<string>&     tokens,
            const string&       delimiters = " " );


#endif //__SCANOPT_H__
