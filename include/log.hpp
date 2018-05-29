// Note: this file contains code that was derived from the Loguru logging
// library (by Emil Ernerfeldt), which is developed in the public domain:
// https://github.com/emilk/loguru

#ifndef LOG_HPP
#define LOG_HPP

#include <sstream>


namespace log
{
  // message status types
  enum Status {
    WARNING,
    ARGERROR,
    ERROR,
    DEBUG
  };


  /// \class Logger
  /// A stream-style logging object
  class Logger
  {
  public:
    Logger(Status status);
    Logger(Status status, int exitcode);
    ~Logger() noexcept(false);

    template<typename T>
    Logger& operator << (const T& t) {
      ss << t;
      return *this;
    }

    // for std::endl
    Logger& operator << (std::ostream&(*f)(std::ostream&)) {
      f(ss);
      return *this;
    }

  private:
    void log(std::string content);
    Status msg_status;
    int retcode;
    bool exitout;
    std::ostringstream ss;
  };

  // terminal output colors 
  const char* black();
  const char* red();
  const char* green();
  const char* yellow();
  const char* blue();
  const char* purple();
  const char* cyan();
  const char* light_gray();
  const char* white();
  const char* light_red();
  const char* dim();
  // formating
  const char* bold();
  const char* underline();
  // reset
  const char* reset();
}

// macros for easy logging of standard message types
#define warning       log::Logger(log::WARNING)
#define error         log::Logger(log::ERROR)
#define errExit(ret)  log::Logger(log::ERROR, ret)
#define argError      log::Logger(log::ARGERROR, 2)
#define debug         log::Logger(log::DEBUG)

#endif
