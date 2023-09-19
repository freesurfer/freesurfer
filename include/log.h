#pragma once

#include <iostream>
#include <sstream>

#include "diag.h"
#include "fsinit.h"


// terminal colors 
namespace term {
  const char* black();
  const char* red();
  const char* green();
  const char* yellow();
  const char* blue();
  const char* purple();
  const char* cyan();
  const char* white();
  const char* bold();
  const char* dim();
  const char* underline();
  const char* reset();
}


// global settings
void throwExceptions(bool setting);
void setErrorLog(const std::string& filename);

namespace fs {
  namespace detail {

    void writeToErrorLog(const std::string& message);
    void errorExit(int code);

    struct logger
    {
      template<typename T> logger& operator << (const T& t) { ss << t; return *this; }
      logger& operator << (std::ostream&(*f)(std::ostream&)) { f(ss); return *this; }
      std::ostringstream ss;
    };
  }
  
  struct fatal : public detail::logger
  {
    int ret;
    fatal(int err = 1) : ret(err) {}
    ~fatal() {
      std::cerr << term::red() << "error: " << term::reset() << this->ss.str() << "\n";
      detail::writeToErrorLog(this->ss.str());
      detail::errorExit(this->ret);
    }
  };

  struct error : public detail::logger
  {
    ~error() {
      std::cerr << term::red() << "error: " << term::reset() << this->ss.str() << "\n";
      detail::writeToErrorLog(this->ss.str());
    }
  };

  struct warning : public detail::logger
  {
    ~warning() { std::cerr << term::yellow() << "warning: " << term::reset() << this->ss.str() << "\n"; }
  };

  struct debug : public detail::logger
  {
    ~debug() { if (DIAG_VERBOSE_ON) std::cerr << term::cyan() << "debug: " << term::reset() << this->ss.str() << "\n"; }
  };

}

// old c-style error functions
void    ErrorExit(int ecode, const char *fmt, ...);
void    ErrorPrintf(int ecode, const char *fmt, ...);
#define ErrorReturn(ret, args)  { ErrorPrintf args; return(ret); }
#define ErrorInit(a, b, c) FSinit();

// global error
extern int Gerror;

// error codes
#define NO_ERROR              0
#define ERROR_NONE            NO_ERROR
#define ERROR_NO_FILE         -1
#define ERROR_NOFILE          ERROR_NO_FILE
#define ERROR_NO_MEMORY       -2
#define ERROR_NOMEMORY        ERROR_NO_MEMORY
#define ERROR_UNSUPPORTED     -3
#define ERROR_BADPARM         -4
#define ERROR_BAD_PARM        ERROR_BADPARM
#define ERROR_BADFILE         -5
#define ERROR_BAD_FILE        ERROR_BADFILE
#define ERROR_SIZE            -6
#define ERROR_BADLOOP         -7
#define ERROR_OUT_OF_BOUNDS   -8
