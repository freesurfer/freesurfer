#include <sstream>
#include <iostream>
#include <cassert>
#include <algorithm>

#include "argparse.h"
#include "version.h"
#include "utils.h"


/// Makes sure option key is valid (i.e it isn't empty and
/// has the correct number of leading dashes)
static std::string verifyOption(const std::string& name)
{
  if (name.empty())
    fs::fatal() << "invalid argument configuration. Argument names must not be empty";
  if ((name.size() == 2 && name[0] != '-') || name.size() == 3)
    fs::fatal() << "invalid argument configuration for '" << name << "'. Short names must begin with '-'";
  if (name.size() > 3 && (name[0] != '-' || name[1] != '-'))
    fs::fatal() << "invalid argument configuration for '" << name << "'. Multi-character names must begin with '--'";
  return name;
}


/// Strips an option flag 'name' of its leading dashes
static std::string strip(const std::string& name)
{
  size_t begin = 0;
  begin += name.size() > 0 ? name[0] == '-' : 0;
  begin += name.size() > 3 ? name[1] == '-' : 0;
  return name.substr(begin);
}


/// Converts a string to boolean
static bool stob(const std::string& str)
{
  std::string upper = std::string(str);
  std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
  if ((upper == "TRUE") || (upper == "YES") || (upper == "1")) return true;
  if ((upper == "FALSE") || (upper == "NO") || (upper == "0")) return false;
  throw std::invalid_argument("invalid boolean");
}


/// This creates an invalid, empty argument 
ArgumentParser::Argument::Argument()
  : short_name(""),
    name(""),
    argtype(ArgType::Unknown),
    consumed(0),
    fixed(true),
    required(false),
    positional(false), 
    valid(false),
    fixed_nargs(0)
{}


ArgumentParser::Argument::Argument(const ArgumentParser::String& _short_name, const ArgumentParser::String& _name,
                                   char nargs, ArgType _argtype, bool _required)
  : short_name(_short_name),
    name(_name),
    argtype(_argtype),
    consumed(0),
    required(_required),
    valid(true)
{
  // check whether this argument is a positional arg or option 
  if (name[0] != '-') {
    // positional
    positional = true;
    // positionals should always be required
    required = true;
    // positionals must require more than 0 inputs
    if (nargs == 0) nargs = 1;
    if (nargs == '*') nargs = '+';
  } else {
    // flagged
    positional = false;
  }

  // set the required number of inputs
  if (nargs == '+' || nargs == '*') {
    // variable args
    variable_nargs = nargs;
    if (variable_nargs == '+') min_args = 1;
    if (variable_nargs == '*') min_args = 0;
    fixed = false;
  } else {
    // fixed args
    fixed_nargs = min_args = nargs;
    fixed = true;
  }

  // set type defaults if not set (default to string except for option flags)
  if (argtype == ArgType::Unknown) {
    if (min_args == 0) argtype = ArgType::Bool;
    else argtype = ArgType::String;
  }

  // check for illogical option flags
  if ((min_args == 0) && (argtype != ArgType::Bool)) {
    fs::fatal() << "invalid argument configuration for '" << canonicalName() << "'. "
                << "Option flags that accept no input must be of type ArgType::Bool";
  }
  if ((min_args == 0) && required) {
    fs::fatal() << "invalid argument configuration for '" << canonicalName() << "'. "
                << "Required flags must accept at least one input";
  }
}


/// Checks if the correct amount of inputs have been parsed for this argument
void ArgumentParser::Argument::validate()
{
  if (positional && consumed < min_args)
    fs::fatal(2) << "not enough positional arguments supplied";
  if (fixed && fixed_nargs != consumed)
    fs::fatal(2) << "not enough inputs passed to option '" << canonicalName() << "' (expected " << fixed_nargs << ")";
  if (!fixed && variable_nargs == '+' && consumed < 1)
    fs::fatal(2) << "option '" << canonicalName() << "' requires at least one input";
}


/// Returns a valid name for the configured argument
ArgumentParser::String ArgumentParser::Argument::canonicalName() const
{
    return (name.empty()) ? short_name : name;
}


/// Returns the argument type as a string
ArgumentParser::String ArgumentParser::Argument::typeName() const
{
  switch (argtype) {
    case ArgType::Int     : return "int";
    case ArgType::Float   : return "float";
    case ArgType::Bool    : return "bool";
    case ArgType::String  : return "std::string";
    default               : return "unknown";
  }
}


/// Configures a command line argument. This can either be a positional argument or an
/// option (if 'name' is preceded with one or two dashes). 
void ArgumentParser::addArgument(const ArgumentParser::String& name, char nargs, ArgType argtype, bool required)
{
  // check if key represents a positional or flagged argument
  if (name[0] != '-') {
    // positional argument (must be required)
    Argument arg("", name, nargs, argtype, true);
    insertArgument(arg);
  } else {
    // flagged option
    if (name.size() > 2) {
      Argument arg("", verifyOption(name), nargs, argtype, required);
      insertArgument(arg);
    } else {
      Argument arg(verifyOption(name), "", nargs, argtype, required);
      insertArgument(arg);
    }
  }
}


/// Configures a command line option. 'short_name' must be prefaced with a single dash,
/// and 'name' must be prefaced with two
void ArgumentParser::addArgument(const ArgumentParser::String& short_name, const ArgumentParser::String& name, char nargs, ArgType argtype, bool required)
{
  Argument arg(verifyOption(short_name), verifyOption(name), nargs, argtype, required);
  insertArgument(arg);
}


/// Configures a --help option that prints the supplied xml help text
void ArgumentParser::addHelp(const unsigned char *text, unsigned int size)
{
  helptext = text;
  helptextsize = size;
  addArgument("-h", "--help", 0, Bool, false);
}


/// The main parsing routine. This will error out if the command line input
/// does not match the argument configuration
void ArgumentParser::parse(size_t ac, char** av)
{
  // TODO progname should be an optionally-set member variable, but extracting from argv for now
  std::string progname;
  std::stringstream progss(av[0]);
  while (std::getline(progss, progname, '/'));  // do nothing in loop, just get basename

  // create argv vector
  StringVector argv = StringVector(av, av + ac);

  // name the app
  if (!argv.empty()) app_name = argv[0];

  // first do a quick and dirty sweep of the options, making sure the minimum
  // amount of arguments have been provided
  for (StringVector::const_iterator in = argv.begin() + 1; in < argv.end(); ++in) {
    String element = *in;
    if (index.count(element) != 0) {
      // count the number of input args following this option
      unsigned int args_following = 0;
      for (StringVector::const_iterator fin = in + 1 ; fin < argv.end() ; fin++) {
        String future = *fin;
        if (index.count(future) == 0) args_following++;
      }
      if (arguments[index[element]].min_args > args_following) {
        fs::fatal(2) << "not enough inputs supplied to '" << element << "'";
      }
    }
  }

  // init our cursors
  unsigned int posidx = -1;
  Argument active;

  // now do the real parsing, and iterate over each element in the input array
  for (StringVector::const_iterator in = argv.begin() + 1; in < argv.end(); ++in) {
    String element = *in;
    // check if the element is a known flag
    if (index.count(element) == 0) {
      // it's not, it must be a positional argument
      if (!active.valid && !positionals.empty()) active = positionals[++posidx];

      // has the current active argument reached its required number of inputs?
      if (active.fixed && (active.fixed_nargs == active.consumed)) {
        // if so, let's get the next set of positional arguments
        posidx++;
        if (posidx >= positionals.size()) {
          fs::fatal(2) << "unexpected argument '" << element << "'";
        } else {
          active = positionals[posidx];
        }
      } else if (!active.fixed && (active.consumed >= active.min_args)) {
        // if we're parsing a variable number of inputs, we need to make sure that we've left
        // enough for any required positional arguments at the end. To do this, we need to find
        // out the minimum possible number of positional arguments left to parse
        int minimum = 0;
        for (unsigned int idx = posidx + 1 ; idx < positionals.size() ; idx++) minimum += positionals[idx].min_args;
        // now we jump ahead to find the actual number of positional arguments remaining
        if (minimum > 0) {
          int positionals_remaining = 0;
          for (StringVector::const_iterator fin = in + 1 ; fin < argv.end() ; fin++) {
            String future = *fin;
            if (future[0] == '-') {
              // future argument is just a flag, so we ignore it as well as it's potential inputs
              positionals_remaining -= arguments[index[future]].min_args;
            } else {
              positionals_remaining++;
            }
          }
          // move on to the next set of positionals if we've reached our limit
          if (minimum > positionals_remaining) active = positionals[++posidx];
        }
      }

      // here we add the input element, but first test if it can be converted to the proper arg type
      size_t N = index[active.canonicalName()];
      try {
        if (active.fixed && active.fixed_nargs == 1) {
          switch(active.argtype) {
            case ArgType::Int   : variables[N].castTo<int>() = std::stoi(element); break;
            case ArgType::Float : variables[N].castTo<float>() = std::stof(element); break;
            case ArgType::Bool  : variables[N].castTo<bool>() = stob(element); break;
            default             : variables[N].castTo<String>() = element; break;
          }
        } else {
          switch(active.argtype) {
            case ArgType::Int   : variables[N].castTo<IntVector>().push_back(std::stoi(element)); break;
            case ArgType::Float : variables[N].castTo<FloatVector>().push_back(std::stof(element)); break;
            case ArgType::Bool  : variables[N].castTo<std::vector<bool>>().push_back(stob(element)); break;
            default             : variables[N].castTo<StringVector>().push_back(element); break;
          }
        }
      } catch (...) {
        fs::fatal(2) << "input '" << element << "' cannot be converted to expected type (" << active.typeName() << ")";
      }
      variables[N].exists = true;
      active.consumed++;
    } else {
      // we've found a new flag, now check whether the active argument has parsed enough elements
      if (active.valid) active.validate();
      active = arguments[index[element]];

      // if argument is a flag that accept no inputs, set to true
      if (ArgType::Bool && active.fixed && (active.fixed_nargs == 0)) {
        size_t N = index[active.canonicalName()];
        variables[N].castTo<bool>() = true;
        variables[N].exists = true;
      }
    }
  }
  // validate the final argument
  if (active.valid) active.validate();
  
  // check if the default --version or --all-info flags were provided
  int numInfoFlags = 0;
  if (exists("version"))  {
    std::cout << progname << " freesurfer " << getVersion() << std::endl;
    numInfoFlags += 1;
  }
  if (exists("all-info")) {
    std::cout << getAllInfo(ac, av, progname) << std::endl;
    numInfoFlags += 1;
  }
  // exit cleanly if only --version or --all-info commands were used
  if ((numInfoFlags > 0) && (ac - numInfoFlags) == 1) exit(0);

  // check for the help flag
  if ((helptextsize > 0) && (exists("help"))) {
    outputHelpXml(helptext, helptextsize);
    exit(0);
  }

  // check that all of the required arguments have been provided
  for (ArgumentVector::const_iterator it = arguments.begin(); it != arguments.end(); ++it) {
    Argument arg = *it;
    if (arg.required && !exists(arg.canonicalName())) {
      fs::fatal(2) << "missing required input '" << arg.canonicalName() << "'";
    }
  }
}


/// Inserts dashes to the front of a stripped key name if the key
/// exists (has been configured with addArgument("--name"))
ArgumentParser::String ArgumentParser::unstrip(const String& name)
{
  for (IndexMap::iterator it = index.begin(); it != index.end(); it++) {
    if (strip(it->first) == name) return String(it->first);
  }
  return String(name);
}


/// Returns true if inputs for a particular argument key (specified by 'name')
/// were provided on the command line. The 'name' parameter doesn't have
/// to be preceded by option dashes in order for the key to be found
bool ArgumentParser::exists(const String& name)
{
  // first check if name is a valid argument key
  String unstripped = unstrip(name);
  if (index.count(unstripped) == 0) fs::fatal() << "'" << unstripped << "' is not a known argument";
  return variables[index[unstripped]].exists;
}


/// Used by addArgument() to actually insert a typed argument into
/// the ArgumentParser index
void ArgumentParser::insertArgument(const ArgumentParser::Argument& arg)
{
  size_t N = arguments.size();
  arguments.push_back(arg);
  if (arg.positional) positionals.push_back(arg);
  if (arg.fixed && arg.fixed_nargs <= 1) {
    switch(arg.argtype) {
      case ArgType::Int    : variables.push_back(int(0)); break;
      case ArgType::Float  : variables.push_back(float(0)); break;
      case ArgType::Bool   : variables.push_back(false); break;
      case ArgType::String : variables.push_back(String()); break;
      default : fs::fatal() << "unknown argument type for '" << arg.canonicalName() << "'";
    }
  } else {
    switch(arg.argtype) {
      case ArgType::Int    : variables.push_back(IntVector()); break;
      case ArgType::Float  : variables.push_back(FloatVector()); break;
      case ArgType::Bool   : variables.push_back(std::vector<bool>()); break;
      case ArgType::String : variables.push_back(StringVector()); break;
      default : fs::fatal() << "unknown argument type for '" << arg.canonicalName() << "'";
    }
  }

  // make sure name doesn't already exist
  for (IndexMap::iterator it = index.begin(); it != index.end(); it++) {
    String stripped = strip(it->first);
    if (stripped == strip(arg.short_name) || stripped == strip(arg.name)) {
      fs::fatal() << "invalid argument configuration. '" << arg.canonicalName() << "' is used twice";
    }
  }

  // make sure there are no illogical combinations of arguments with a variable number of inputs 
  if (!arg.fixed) {
    if (arg.positional) {
      if (variable_positional) {
        fs::fatal() << "invalid argument configuration for '" << arg.canonicalName() << "'. "
                    << "Two positional arguments cannot both have a variable amount of inputs, "
                    << "as this could lead to an undefined boundary between the two";
      }
      variable_positional = true;
    } else {
      variable_flag = true;
    }
    if (variable_positional && variable_flag) {
      fs::fatal() << "invalid argument configuration for '" << arg.canonicalName() << "'. "
                  << "A positional argument and flagged argument cannot both have a variable "
                  << "amount of inputs, as this could lead to an undefined boundary between the two";
    }
  }

  // add argument name to the index
  if (!arg.short_name.empty()) index[arg.short_name] = N;
  if (!arg.name.empty()) index[arg.name] = N;
}
