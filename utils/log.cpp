#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cstdarg>
#include <vector>
#include <unistd.h>

#include "log.h"


namespace term {

  // Checks to make sure that stdout and stderr are ttys that accept colors.
#if 0
  static bool termAllowsColor() {
    if (!isatty(fileno(stderr))) return false;
    if (!isatty(fileno(stdout))) return false;
    if (const char* term = getenv("TERM")) {
      return 0 == strcmp(term, "cygwin")
          || 0 == strcmp(term, "linux")
          || 0 == strcmp(term, "rxvt-unicode-256color")
          || 0 == strcmp(term, "screen")
          || 0 == strcmp(term, "screen-256color")
          || 0 == strcmp(term, "tmux-256color")
          || 0 == strcmp(term, "xterm")
          || 0 == strcmp(term, "xterm-256color")
          || 0 == strcmp(term, "xterm-termite")
          || 0 == strcmp(term, "xterm-color");
    }
    else return false;
  }
#else
  // colors aren't fully tested yet, so let's always return false until confident
  static bool termAllowsColor() { return false; }
#endif

  const char* black()     { return termAllowsColor() ? "\e[30m" : ""; }
  const char* red()       { return termAllowsColor() ? "\e[31m" : ""; }
  const char* green()     { return termAllowsColor() ? "\e[32m" : ""; }
  const char* yellow()    { return termAllowsColor() ? "\e[33m" : ""; }
  const char* blue()      { return termAllowsColor() ? "\e[34m" : ""; }
  const char* purple()    { return termAllowsColor() ? "\e[35m" : ""; }
  const char* cyan()      { return termAllowsColor() ? "\e[36m" : ""; }
  const char* white()     { return termAllowsColor() ? "\e[37m" : ""; }
  const char* bold()      { return termAllowsColor() ? "\e[1m"  : ""; }
  const char* dim()       { return termAllowsColor() ? "\e[2m"  : ""; }
  const char* underline() { return termAllowsColor() ? "\e[4m"  : ""; }
  const char* reset()     { return termAllowsColor() ? "\e[0m"  : ""; }

}


// static logging options
static std::string errorlog = std::string();
static bool exceptions = false;

int Gerror = NO_ERROR;


/**
  Configures the global error log file. All messages sent to `fs::fatal()` and
  `fs::error()` are written to this file.
*/
void setErrorLog(const std::string& filename)
{
  errorlog = filename;
}


/**
  Configures a global setting to throw exceptions instead of exiting
  when `fs::fatal()` and `ErrorExit()` are called.
*/
void throwExceptions(bool setting)
{
  exceptions = setting;
}


/**
  Writes a message the global error log if configured with `setErrorLog(filename)`.
  **This is for internal use only** - use `fs::error()` to throw errors.
*/
void fs::detail::writeToErrorLog(const std::string& message)
{
  if (errorlog.empty()) return;
  
  std::ofstream file;
  file.open(errorlog, std::ios_base::app);
  file << message << std::endl; 
}


/**
  Exits with a given return code, or throws exceptions if enabled via `throwExceptions(true)`.
  **This is for internal use only** - use `fs::fatal(code)` to throw a fatal error.
*/
void fs::detail::errorExit(int code)
{
  if (exceptions) {
    if (code != 0) throw code;
  } else {
    exit(code);
  }
}


// a rather ugly but safe way to convert a va_list to a string
#define vargsToString(format, message)  \
  va_list vargs;  \
  va_start(vargs, format);  \
  va_list vaCopy;  \
  va_copy(vaCopy, vargs);  \
  const int iLen = std::vsnprintf(NULL, 0, format, vaCopy);  \
  va_end(vaCopy);  \
  std::vector<char> zc(iLen + 1);  \
  std::vsnprintf(zc.data(), zc.size(), format, vargs);  \
  va_end(vargs);  \
  const std::string message = std::string(zc.data(), zc.size());


/**
  A c-style error function that wraps the cpp stream logging utilities. Ideally, `fs::error() << message`
  should be used instead.
*/
void ErrorPrintf(int errcode, const char *format, ...)
{
  Gerror = errcode;
  if (errno) fs::error() << strerror(errno);
  vargsToString(format, message);
  fs::error() << message;
}


/**
  A c-style fatal error that wraps the cpp stream logging utilities. Ideally, `fs::fatal(code) << message`
  should be used instead.
*/
void ErrorExit(int errcode, const char *format, ...)
{
  Gerror = errcode;
  if (errno) fs::error() << strerror(errno);
  vargsToString(format, message);
  fs::fatal(errcode) << message;
}
