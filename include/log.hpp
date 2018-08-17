// Note: this file contains code that was derived from the Loguru logging
// library (by Emil Ernerfeldt), which is developed in the public domain:
// https://github.com/emilk/loguru

#ifndef LOG_HPP
#define LOG_HPP

#include <sstream>


namespace fs
{
  // message status types
  enum MessageStatus {Debug, Warning, Error};

  /// \class Logger
  /// A stream-style logging object
  class Logger
  {
  public:
    Logger(MessageStatus status);
    Logger(MessageStatus status, int exitcode);
    ~Logger();

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
    void print(std::string content);
    MessageStatus msg_status;
    int retcode;
    bool exitout;
    std::ostringstream ss;
  };

  // terminal output colors 
  namespace term
  {
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
}

// macros for easy logging of standard message types
#define fs_debug      fs::Logger(fs::Debug)
#define fs_warning    fs::Logger(fs::Warning)
#define fs_error      fs::Logger(fs::Error)
#define fs_fatal(ret) fs::Logger(fs::Error, ret)

#endif
