#include <iostream>
#include <cstring>
#include <unistd.h>
#include <string>

#include "log.hpp"


namespace fs
{
  Logger::Logger(MessageStatus status) : msg_status(status), exitout(false) {}
  Logger::Logger(MessageStatus status, int exitcode) : msg_status(status), retcode(exitcode), exitout(true) {}

  // destructor calls the print routine
  Logger::~Logger()
  {
    print(ss.str());
    if (exitout) exit(retcode);
  }

  /// prints the actual prefaced message
  void Logger::print(std::string content)
  {
    switch(msg_status) {
      case Warning  : fprintf(stderr, "%swarning:%s %s\n", term::yellow(), term::reset(), content.c_str()); break;
      case Error    : fprintf(stderr, "%serror:%s %s\n", term::red(), term::reset(), content.c_str()); break;
      default       : fprintf(stdout, "%s\n", content.c_str()); break;
    }
  }

  /// checks to make sure the output is a tty that accepts colors
  static bool termAllowsColor()
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

  namespace term
  {
    // colors
    const char* black()      { return termAllowsColor() ? "\e[30m" : ""; }
    const char* red()        { return termAllowsColor() ? "\e[31m" : ""; }
    const char* green()      { return termAllowsColor() ? "\e[32m" : ""; }
    const char* yellow()     { return termAllowsColor() ? "\e[33m" : ""; }
    const char* blue()       { return termAllowsColor() ? "\e[34m" : ""; }
    const char* purple()     { return termAllowsColor() ? "\e[35m" : ""; }
    const char* cyan()       { return termAllowsColor() ? "\e[36m" : ""; }
    const char* light_gray() { return termAllowsColor() ? "\e[37m" : ""; }
    const char* white()      { return termAllowsColor() ? "\e[37m" : ""; }
    const char* light_red()  { return termAllowsColor() ? "\e[91m" : ""; }
    const char* dim()        { return termAllowsColor() ? "\e[2m"  : ""; }

    // formating
    const char* bold()       { return termAllowsColor() ? "\e[1m" : ""; }
    const char* underline()  { return termAllowsColor() ? "\e[4m" : ""; }

    // reset
    const char* reset()      { return termAllowsColor() ? "\e[0m" : ""; }
  }
}