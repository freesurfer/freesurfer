// Note: this file contains code that was derived from the ArgumentParser project
// by Hilton Bristow, which is developed with the BSD 3-Clause License:
// https://github.com/hbristow/argparse

#ifndef ARGPARSE_HPP
#define ARGPARSE_HPP

#if __cplusplus >= 201103L
#include <unordered_map>
typedef std::unordered_map<std::string, size_t> IndexMap;
#else
#include <map>
typedef std::map<std::string, size_t> IndexMap;
#endif
#include <string>
#include <vector>
#include <typeinfo>
#include <stdexcept>

#include "log.hpp"


enum ArgType {Unknown, Bool, String, Int, Float};

/// \class ArgumentParser
/// \brief A simple command-line argument parser
///
/// ArgumentParser is a simple class that can parse arguments from
/// the command-line or any array of strings. The syntax is familiar to
/// anyone who has used python's ArgumentParser:
/// \code
///   // create a parser and add the options
///   ArgumentParser parser;
///   parser.addArgument("-i", "--inputs", '+');
///   parser.addArgument("-o", "--output", 1);
///   // parse the command-line arguments
///   parser.parse(argc, argv);
///   // get the inputs and iterate over them
///   std::string output = parser.retrieve<std::string>("output");
///   vector<std::string> inputs = parser.retrieve<vector<std::string>>("inputs");
/// \endcode
class ArgumentParser {
private:
  class Any;
  class Argument;
  class PlaceHolder;
  class Holder;
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
    Any() : content(0) {}
    // destructor
    ~Any() { delete content; }
    // inward conversions
    Any(const Any& other) : content(other.content ? other.content->clone() : 0) {}
    template <typename ValueType>
    Any(const ValueType& other) : content(new Holder<ValueType>(other)) {}
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
    bool exists = false;

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

  ArgumentParser() : variable_positional(false), variable_flag(false) {}

  void addArgument(const String& name, char nargs = 0, ArgType argtype = Unknown, bool required = false);
  void addArgument(const String& short_name, const String& name, char nargs = 0, ArgType argtype = Unknown, bool required = false);

  void parse(size_t argc, const char** argv);
  void parse(const StringVector& argv);

  bool exists(const String& name);


  /// Returns the parsed inputs for a given argument key (specified by 'name').
  /// The output must be correctly typecasted based on the configured arg type,
  /// but if it's not, the error message should hopefully provide the correct fix
  template <typename T>
  T& retrieve(const String& name) {
    String unstripped = unstrip(name);
    if (index.count(unstripped) == 0) errExit(1) << "'" << unstripped << "' is not a known argument";
    size_t N = index[unstripped];
    // try to cast the arguments
    try {
      return variables[N].castTo<T>();
    } catch (std::bad_cast) {
      // if casting fails, print out a VERY detailed debug message
      String fulltype, sentence_starter;
      if (arguments[N].fixed && (arguments[N].fixed_nargs <= 1)) {
        fulltype = arguments[N].typeName();
        sentence_starter = "This input is";
      } else {
        fulltype = "std::vector<" + arguments[N].typeName() + ">";
        sentence_starter = "These inputs are";
      }
      errExit(1) << "invalid cast of argument '" << name << "'. " << sentence_starter << " of type '"
                 << arguments[N].typeName() << "' and should be retrieved via " << log::dim()
                 << "retrieve<" << fulltype << ">(\"" << name << "\")" << log::reset() << ". " 
                 << "To change the expected type, modify the call to "
                 << log::dim() << "addArgument()" << log::reset();
    }
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
};

#endif