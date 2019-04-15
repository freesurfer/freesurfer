#ifndef LOG_H
#define LOG_H

#include <sstream>


/// \class StreamLogger
///
/// A stream-style logging object for printing status messages.
///
/// In general, this class shouldn't used directly. Status messages should
/// be printed via the standard logging macros. For example:
///
///     logWarning << "exceeded standard iteration count";
///     logError << "vertex has no neighbors";
///     logFatal(1) << "could not find file";

class StreamLogger {
public:
  enum MessageStatus {Warning, Error, Debug};

  StreamLogger(MessageStatus status) : status(status), exitout(false) {};
  StreamLogger(MessageStatus status, int exitcode) : status(status), retcode(exitcode), exitout(true) {};
  ~StreamLogger();

  template<typename T>
  StreamLogger& operator << (const T& t) {
    ss << t;
    return *this;
  }

  // for std::endl
  StreamLogger& operator << (std::ostream&(*f)(std::ostream&)) {
    f(ss);
    return *this;
  }

private:
  MessageStatus status;
  int retcode;
  bool exitout;
  std::ostringstream ss;
};

// terminal output colors 
namespace term {
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
#define logDebug      StreamLogger(StreamLogger::Debug)
#define logWarning    StreamLogger(StreamLogger::Warning)
#define logError      StreamLogger(StreamLogger::Error)
#define logFatal(ret) StreamLogger(StreamLogger::Error, ret)

#endif
