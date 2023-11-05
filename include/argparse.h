// Note: this file contains code that was derived from the ArgumentParser project
// by Hilton Bristow, which is developed with the BSD 3-Clause License:
// https://github.com/hbristow/argparse

#ifndef ARGPARSE_H
#define ARGPARSE_H
#include <string>

#if __cplusplus >= 201103L
#include <unordered_map>
typedef std::unordered_map<std::string, size_t> IndexMap;
#else
#include <map>
typedef std::map<std::string, size_t> IndexMap;
#endif
#include <vector>
#include <typeinfo>
#include <stdexcept>

#include "log.h"

/*!

  \class ArgumentParser
  \brief A simple command-line argument parser
  
  ArgumentParser is a class that can parse arguments from
  the command-line (or from any array of strings). The syntax is
  similar to python's argparse. Here's an example for a program
  that requires specifying one output file and at least one input file
  via flagged arguments:
  
  int main(int argc, const char **argv)
  {
    // create a parser and add the options
    ArgumentParser parser;
    parser.addArgument("-i", "--inputs", '+', String, true);
    parser.addArgument("-o", "--output", 1, String, true);

    // parse the command-line
    parser.parse(argc, argv);

    // retrieve the args
    vector<std::string> inputs = parser.retrieve<vector<std::string>>("inputs");
    std::string output = parser.retrieve<std::string>("output");
  }

  To configure a parser, an ArgumentParser object must first be instantiated. Valid
  positional and flagged arguments can then be configured with the addArgument() member function,
  which accepts the following parameters. At the absolute minimum, 'name' must be specified:

    addArgument(shortname, name, nargs, argtype, required)

  details:

    std::string  shortname  Shortened flag name, such as "-i". This is optional and only
                            allowed for flagged arguments.
    
    std::string  name       Positional argument name or full flagged argument name,
                            such as "--input". Names with leading dashes will be registered
                            as flagged arguments, while those without will be registered as
                            positional arguments.

    optional parameters:

    char         nargs      Number of required inputs. '*' allows for a variable number of
                            args, and '+' specifies that at least 1 arg is required. The
                            default is 0 for flagged arguments and 1 for positional arguments.

    ArgType      argtype    Specifies the expected input type. Valid types are:
                            Bool, String, Int, or Float. If the user provides an incorrect
                            argument type on the command-line, a usage error will be
                            thrown. The default type is String.

    bool         required   Species whether the argument is required. The default is false for
                            flagged arguments, but always true for positionals.

  After configuring, the parse() function (as used in the example above) must be called
  once to validate the input. Then, the retrieve function can be used to return what has
  been parsed. The result must be typecasted to the type that was specified in addArgument.
  Either the 'shortname' or 'name' can be used as a reference to retrieve an input. For example:

    std::string output = parser.retrieve<std::string>("output");

  Flagged arguments that do not accept any inputs are usually meant to turn things on/off. So
  the easiest way to check whether an option has simply been provided on the command-line is
  with the exists() function, which returns true if the argument exists. For example:
  
    bool do_thing = parser.exists("thing")

*/


enum ArgType {Unknown, Bool, String, Int, Float};


class ArgumentParser {
private:
  class Any;
  class PlaceHolder;
  class Holder;
  struct Argument;
  typedef std::string String;
  typedef std::vector<Any> AnyVector;
  typedef std::vector<String> StringVector;
  typedef std::vector<int> IntVector;
  typedef std::vector<float> FloatVector;
  typedef std::vector<Argument> ArgumentVector;

  // --------------------- type-erasure internal storage ----------------------

  class Any {
  public:
    // constructor
    Any() : exists(false), content(0) {}
    // destructor
    ~Any() { delete content; }
    // inward conversions
    Any(const Any& other) : exists(false), content(other.content ? other.content->clone() : 0) {}
    template <typename ValueType>
    Any(const ValueType& other) : exists(false), content(new Holder<ValueType>(other)) {}
    Any& swap(Any& other) {
      std::swap(content, other.content);
      return *this;
    }
    Any& operator=(const Any& rhs) {
      Any tmp(rhs);
      return swap(tmp);
    }
    template <typename ValueType>
    Any& operator=(const ValueType& rhs) {
      Any tmp(rhs);
      return swap(tmp);
    }
    // outward conversions
    template <typename ValueType>
    ValueType* toPtr() const {
      return content->type_info() == typeid(ValueType) ? &static_cast<Holder<ValueType>*>(content)->held : 0;
    }
    template <typename ValueType>
    ValueType& castTo() {
      if (!toPtr<ValueType>()) throw std::bad_cast();
      return *toPtr<ValueType>();
    }
    template <typename ValueType>
    const ValueType& castTo() const {
      if (!toPtr<ValueType>()) throw std::bad_cast();
      return *toPtr<ValueType>();
    }
    bool exists;

  private:
    // inner placeholder interface
    class PlaceHolder {
    public:
      virtual ~PlaceHolder() {}
      virtual const std::type_info& type_info() const = 0;
      virtual PlaceHolder* clone() const = 0;
    };
    // Inner template concrete instantiation of PlaceHolder
    template <typename ValueType>
    class Holder : public PlaceHolder {
    public:
      ValueType held;
      Holder(const ValueType& value) : held(value) {}
      virtual const std::type_info& type_info() const { return typeid(ValueType); }
      virtual PlaceHolder* clone() const { return new Holder(held); }
    };
    PlaceHolder* content;
  };

  // -------------------------------- argument --------------------------------

  struct Argument {
    Argument();
    Argument(const String& short_name, const String& name, char nargs, ArgType argtype, bool required);

    String short_name;
    String name;
    ArgType argtype;
    unsigned int min_args;
    unsigned int consumed;
    bool fixed;
    bool required;
    bool positional;
    bool valid;
    union { size_t fixed_nargs; char variable_nargs; };

    void validate();
    String canonicalName() const;
    String typeName() const;
  };

public:
  
  // ----------------------------- argument parser ----------------------------

  ArgumentParser() : variable_positional(false), variable_flag(false) {
    // add default --all-info and --version flags
    addArgument("--all-info", 0, Bool, false);
    addArgument("--version", 0, Bool, false);
  }

  void addArgument(const String& name, char nargs = 0, ArgType argtype = Unknown, bool required = false);
  void addArgument(const String& short_name, const String& name, char nargs = 0, ArgType argtype = Unknown, bool required = false);

  void addHelp(const unsigned char *text, unsigned int size);

  void parse(size_t argc, char** argv);

  bool exists(const String& name);

  /// Returns the parsed inputs for a given argument key (specified by 'name').
  /// The output must be correctly typecasted based on the configured arg type,
  /// but if it's not, the error message should hopefully provide the correct fix
  template <typename T>
  T retrieve(const String& name) {
    String unstripped = unstrip(name);
    if (index.count(unstripped) == 0) fs::fatal() << "'" << unstripped << "' is not a known argument";
    size_t N = index[unstripped];
    T retrieved{};
    // try to cast the arguments
    try {
      retrieved = variables[N].castTo<T>();
    } catch (std::bad_cast&) {
      // if casting fails, print out a VERY detailed debug message
      String fulltype, sentence_starter;
      if (arguments[N].fixed && (arguments[N].fixed_nargs <= 1)) {
        fulltype = arguments[N].typeName();
        sentence_starter = "This input is";
      } else {
        fulltype = "std::vector<" + arguments[N].typeName() + ">";
        sentence_starter = "These inputs are";
      }
      fs::fatal() << "invalid cast of argument '" << name << "'. " << sentence_starter << " of type '"
                  << arguments[N].typeName() << "' and should be retrieved via " << term::dim()
                  << "retrieve<" << fulltype << ">(\"" << name << "\")" << term::reset() << ". " 
                  << "To change the expected type, modify the call to "
                  << term::dim() << "addArgument()" << term::reset();
    }
    return retrieved;
  }

private:
  // utility
  String unstrip(const String& name);
  void insertArgument(const Argument& arg);

  // member variables
  IndexMap index;
  bool variable_positional;
  bool variable_flag;
  String app_name;
  ArgumentVector arguments;
  ArgumentVector positionals;
  AnyVector variables;
  const unsigned char *helptext;
  unsigned int helptextsize = 0;
};

#endif
