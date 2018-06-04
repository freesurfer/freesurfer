#include <iostream>
#include <cstring>
#include <unistd.h>
#include <string>

#include "log.hpp"


namespace Log
{
  Logger::Logger(Status status) : msg_status(status), exitout(false) {}
  Logger::Logger(Status status, int exitcode) : msg_status(status), retcode(exitcode), exitout(true) {}

  // destructor calls the print routine
  Logger::~Logger()
  {
    log(ss.str());
    if (exitout) exit(retcode);
  }

  /// prints the actual prefaced message
  void Logger::log(std::string content)
  {
    switch(msg_status) {
      case WARNING  : fprintf(stderr, "%swarning:%s %s\n", yellow(), reset(), content.c_str()); break;
      case ARGERROR : fprintf(stderr, "%sargument error:%s %s\n", red(), reset(), content.c_str()); break;
      case ERROR    : fprintf(stderr, "%serror:%s %s\n", red(), reset(), content.c_str()); break;
      default       : fprintf(stdout, "%s\n", content.c_str()); break;
    }
  }

  /// checks to make sure the output is a tty that accepts colors
  static bool term_allows_color()
  {
    if (!isatty(fileno(stderr))) return false;
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

  // colors
  const char* black()      { return term_allows_color() ? "\e[30m" : ""; }
  const char* red()        { return term_allows_color() ? "\e[31m" : ""; }
  const char* green()      { return term_allows_color() ? "\e[32m" : ""; }
  const char* yellow()     { return term_allows_color() ? "\e[33m" : ""; }
  const char* blue()       { return term_allows_color() ? "\e[34m" : ""; }
  const char* purple()     { return term_allows_color() ? "\e[35m" : ""; }
  const char* cyan()       { return term_allows_color() ? "\e[36m" : ""; }
  const char* light_gray() { return term_allows_color() ? "\e[37m" : ""; }
  const char* white()      { return term_allows_color() ? "\e[37m" : ""; }
  const char* light_red()  { return term_allows_color() ? "\e[91m" : ""; }
  const char* dim()        { return term_allows_color() ? "\e[2m"  : ""; }

  // formating
  const char* bold()       { return term_allows_color() ? "\e[1m" : ""; }
  const char* underline()  { return term_allows_color() ? "\e[4m" : ""; }

  // reset
  const char* reset()      { return term_allows_color() ? "\e[0m" : ""; }
}