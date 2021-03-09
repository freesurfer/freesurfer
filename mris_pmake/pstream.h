/**
 * @brief PStreams - POSIX Process I/O for C++
 *
 * Defines classes redi::ipstream, redi::opstream, redi::pstream
 * and redi::rpstream.
 */
/*
 * Original Author: Jonathan Wakely
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef REDI_PSTREAM_H_SEEN
#define REDI_PSTREAM_H_SEEN

#include <ios>
#include <streambuf>
#include <istream>
#include <ostream>
#include <string>
#include <vector>
#include <algorithm>    // for min()
#include <cstring>      // for memcpy(), memmove() etc.
#include <cerrno>       // for errno
#include <cstddef>      // for size_t
#include <cstdlib>      // for exit()
#include <sys/types.h>  // for pid_t
#include <sys/wait.h>   // for waitpid()
#include <unistd.h>     // for pipe() fork() exec() and filedes functions
#include <signal.h>     // for kill()
#if REDI_EVISCERATE_PSTREAMS
#include <stdio.h>       // for FILE, fdopen()
#endif


/// The library version.
#define PSTREAMS_VERSION 0x0050   // 0.5.0

/**
 *  @namespace redi
 *  @brief  All PStreams classes are declared in namespace redi.
 *
 *  Like the standard IOStreams, PStreams is a set of class templates,
 *  taking a character type and traits type. As with the standard streams
 *  they are most likely to be used with @c char and the default
 *  traits type, so typedefs for this most common case are provided.
 *
 *  The @c pstream_common class template is not intended to be used directly,
 *  it is used internally to provide the common functionality for the
 *  other stream classes.
 */
namespace redi {
/// Common base class providing constants and typenames.
struct pstreams {
  /// Type used to specify how to connect to the process.
  typedef std::ios_base::openmode           pmode;

  /// Type used to hold the arguments for a command.
  typedef std::vector<std::string>          argv_type;

  static const pmode pstdin  = std::ios_base::out; ///< Write to stdin
  static const pmode pstdout = std::ios_base::in;  ///< Read from stdout
  static const pmode pstderr = std::ios_base::app; ///< Read from stderr

protected:
  enum {
    bufsz = 32
  }
  ;  ///< Size of pstreambuf buffers.
  enum {
    pbsz  = 2
  }
  ;   ///< Number of putback characters kept.
};

/// Class template for stream buffer.
template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_pstreambuf
      : public std::basic_streambuf<CharT, Traits>
      , public pstreams {
public:
  // Type definitions for dependent types
  typedef CharT                             char_type;
  typedef Traits                            traits_type;
  typedef typename traits_type::int_type    int_type;
  typedef typename traits_type::off_type    off_type;
  typedef typename traits_type::pos_type    pos_type;
  /// Type used for file descriptors.
  typedef int                               fd_type;
  /** @deprecated use fd_type instead. */
  typedef fd_type                           fd_t;

  /// Default constructor.
  basic_pstreambuf();

  /// Constructor that initialises the buffer with @a command.
  basic_pstreambuf(const std::string& command, pmode mode);

  /// Constructor that initialises the buffer with @a file and @a argv.
  basic_pstreambuf( const std::string& file,
                    const argv_type& argv,
                    pmode mode );

  /// Destructor.
  ~basic_pstreambuf();

  /// Initialise the stream buffer with @a command.
  basic_pstreambuf*
  open(const std::string& command, pmode mode);

  /// Initialise the stream buffer with @a file and @a argv.
  basic_pstreambuf*
  open(const std::string& file, const argv_type& argv, pmode mode);

  /// Close the stream buffer and wait for the process to exit.
  basic_pstreambuf*
  close();

  /// Send a signal to the process.
  basic_pstreambuf*
  kill(int signal = SIGTERM);

  /// Close the pipe connected to the process' stdin.
  void
  peof();

  /// Change active input source.
  bool
  read_err(bool readerr = true);

  /// Report whether the stream buffer has been initialised.
  bool
  is_open() const;

  /// Report whether the process has exited.
  bool
  exited();

#if REDI_EVISCERATE_PSTREAMS
  /// Obtain FILE pointers for each of the process' standard streams.
  std::size_t
  fopen(std::FILE*& in, std::FILE*& out, std::FILE*& err);
#endif

  /// Return the exit status of the process.
  int
  status() const;

  /// Return the error number for the most recent failed operation.
  int
  error() const;

protected:
  /// Transfer characters to the pipe when character buffer overflows.
  int_type
  overflow(int_type c);

  /// Transfer characters from the pipe when the character buffer is empty.
  int_type
  underflow();

  /// Make a character available to be returned by the next extraction.
  int_type
  pbackfail(int_type c = traits_type::eof());

  /// Write any buffered characters to the stream.
  int
  sync();

  /// Insert multiple characters into the pipe.
  std::streamsize
  xsputn(const char_type* s, std::streamsize n);

  /// Insert a sequence of characters into the pipe.
  std::streamsize
  write(char_type* s, std::streamsize n);

  /// Extract a sequence of characters from the pipe.
  std::streamsize
  read(char_type* s, std::streamsize n);

protected:
  /// Enumerated type to indicate whether stdout or stderr is to be read.
  enum buf_read_src {
    rsrc_out = 0, rsrc_err = 1
  };

  /// Initialise pipes and fork process.
  pid_t
  fork(pmode mode);

  /// Wait for the child process to exit.
  int
  wait(bool nohang = false);

  /// Return the file descriptor for the output pipe.
  fd_type&
  wpipe();

  /// Return the file descriptor for the active input pipe.
  fd_type&
  rpipe();

  /// Return the file descriptor for the specified input pipe.
  fd_type&
  rpipe(buf_read_src which);

  void
  create_buffers(pmode mode);

  void
  destroy_buffers(pmode mode);

  /// Writes buffered characters to the process' stdin pipe.
  bool
  empty_buffer();

  bool
  fill_buffer();

  /// Return the active input buffer.
  char_type*
  rbuffer();

  buf_read_src
  switch_read_buffer(buf_read_src);

private:
  basic_pstreambuf(const basic_pstreambuf&);
  basic_pstreambuf& operator=(const basic_pstreambuf&);

  void
  init_rbuffers();

  pid_t         ppid_;        // pid of process
  fd_type       wpipe_;       // pipe used to write to process' stdin
  fd_type       rpipe_[2];    // two pipes to read from, stdout and stderr
  char_type*    wbuffer_;
  char_type*    rbuffer_[2];
  char_type*    rbufstate_[3];
  /// Index into rpipe_[] to indicate active source for read operations.
  buf_read_src  rsrc_;
  int           status_;      // hold exit status of child process
  int           error_;       // hold errno if fork() or exec() fails
};

/// Class template for common base class.
template <typename CharT, typename Traits = std::char_traits<CharT> >
class pstream_common
      : virtual public std::basic_ios<CharT, Traits>
      , virtual public pstreams {
protected:
  typedef basic_pstreambuf<CharT, Traits>       streambuf_type;

  /// Default constructor.
  pstream_common();

  /// Constructor that initialises the stream by starting a process.
  pstream_common(const std::string& command, pmode mode);

  /// Constructor that initialises the stream by starting a process.
  pstream_common( const std::string& file,
                  const argv_type& argv,
                  pmode mode );

  /// Pure virtual destructor.
  virtual
  ~pstream_common() = 0;

  /// Start a process.
  void
  do_open(const std::string& command, pmode mode);

  /// Start a process.
  void
  do_open(const std::string& file, const argv_type& argv, pmode mode);

public:
  /// Close the pipe.
  void
  close();

  /// Report whether the stream's buffer has been initialised.
  bool
  is_open() const;

  /// Return the command used to initialise the stream.
  const std::string&
  command() const;

  /// Return a pointer to the stream buffer.
  streambuf_type*
  rdbuf() const;

#if REDI_EVISCERATE_PSTREAMS
  /// Obtain FILE pointers for each of the process' standard streams.
  std::size_t
  fopen(std::FILE*& in, std::FILE*& out, std::FILE*& err);
#endif

protected:
  std::string       command_; ///< The command used to start the process.
  streambuf_type    buf_;     ///< The stream buffer.
};


/**
 * @class basic_ipstream
 * @brief Class template for Input PStreams.
 *
 * Reading from an ipstream reads the command's standard output and/or
 * standard error (depending on how the ipstream is opened)
 * and the command's standard input is the same as that of the process
 * that created the object, unless altered by the command itself.
 */

template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_ipstream
      : public std::basic_istream<CharT, Traits>
      , public pstream_common<CharT, Traits>
      , virtual public pstreams {
  typedef std::basic_istream<CharT, Traits>     istream_type;
  typedef pstream_common<CharT, Traits>         pbase_type;

  using pbase_type::buf_;  // declare name in this scope

public:
  /// Type used to specify how to connect to the process.
  typedef typename pbase_type::pmode            pmode;

  /// Type used to hold the arguments for a command.
  typedef typename pbase_type::argv_type        argv_type;

  /// Default constructor, creates an uninitialised stream.
  basic_ipstream()
      : istream_type(NULL), pbase_type() { }

  /**
   * @brief Constructor that initialises the stream by starting a process.
   *
   * Initialises the stream buffer by calling do_open() with the supplied
   * arguments.
   *
   * @param command  a string containing a shell command.
   * @param mode     the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, pmode)
   */
  basic_ipstream(const std::string& command, pmode mode = pstdout)
      : istream_type(NULL), pbase_type(command, mode|pstdout) { }

  /**
   * @brief Constructor that initialises the stream by starting a process.
   *
   * Initialises the stream buffer by calling do_open() with the supplied
   * arguments.
   *
   * @param file  a string containing the pathname of a program to execute.
   * @param argv  a vector of argument strings passed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, const argv_type&, pmode)
   */
  basic_ipstream( const std::string& file,
                  const argv_type& argv,
                  pmode mode = pstdout )
      : istream_type(NULL), pbase_type(file, argv, mode|pstdout) { }

  /**
   * @brief Destructor.
   *
   * Closes the stream and waits for the child to exit.
   */
  ~basic_ipstream() { }

  /**
   * @brief Start a process.
   *
   * Calls do_open( @a %command , @a mode|pstdout ).
   *
   * @param command  a string containing a shell command.
   * @param mode     the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, pmode)
   */
  void
  open(const std::string& command, pmode mode = pstdout) {
    this->do_open(command, mode|pstdout);
  }

  /**
   * @brief Start a process.
   *
   * Calls do_open( @a file , @a argv , @a mode|pstdout ).
   *
   * @param file  a string containing the pathname of a program to execute.
   * @param argv  a vector of argument strings passed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, const argv_type&, pmode)
   */
  void
  open( const std::string& file,
        const argv_type& argv,
        pmode mode = pstdout ) {
    this->do_open(file, argv, mode|pstdout);
  }

  /**
   * @brief Set streambuf to read from process' @c stdout.
   * @return  @c *this
   */
  basic_ipstream&
  out() {
    this->buf_.read_err(false);
    return *this;
  }

  /**
   * @brief Set streambuf to read from process' @c stderr.
   * @return  @c *this
   */
  basic_ipstream&
  err() {
    this->buf_.read_err(true);
    return *this;
  }
};


/**
 * @class basic_opstream
 * @brief Class template for Output PStreams.
 *
 * Writing to an open opstream writes to the standard input of the command;
 * the command's standard output is the same as that of the process that
 * created the pstream object, unless altered by the command itself.
 */

template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_opstream
      : public std::basic_ostream<CharT, Traits>
      , public pstream_common<CharT, Traits>
      , virtual public pstreams {
  typedef std::basic_ostream<CharT, Traits>     ostream_type;
  typedef pstream_common<CharT, Traits>         pbase_type;

  using pbase_type::buf_;  // declare name in this scope

public:
  /// Type used to specify how to connect to the process.
  typedef typename pbase_type::pmode            pmode;

  /// Type used to hold the arguments for a command.
  typedef typename pbase_type::argv_type        argv_type;

  /// Default constructor, creates an uninitialised stream.
  basic_opstream()
      : ostream_type(NULL), pbase_type() { }

  /**
   * @brief Constructor that initialises the stream by starting a process.
   *
   * Initialises the stream buffer by calling do_open() with the supplied
   * arguments.
   *
   * @param command  a string containing a shell command.
   * @param mode     the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, pmode)
   */
  basic_opstream(const std::string& command, pmode mode = pstdin)
      : ostream_type(NULL), pbase_type(command, mode|pstdin) { }

  /**
   * @brief Constructor that initialises the stream by starting a process.
   *
   * Initialises the stream buffer by calling do_open() with the supplied
   * arguments.
   *
   * @param file  a string containing the pathname of a program to execute.
   * @param argv  a vector of argument strings passed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, const argv_type&, pmode)
   */
  basic_opstream( const std::string& file,
                  const argv_type& argv,
                  pmode mode = pstdin )
      : ostream_type(NULL), pbase_type(file, argv, mode|pstdin) { }

  /**
   * @brief Destructor
   *
   * Closes the stream and waits for the child to exit.
   */
  ~basic_opstream() { }

  /**
   * @brief Start a process.
   *
   * Calls do_open( @a %command , @a mode|pstdin ).
   *
   * @param command  a string containing a shell command.
   * @param mode     the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, pmode)
   */
  void
  open(const std::string& command, pmode mode = pstdin) {
    this->do_open(command, mode|pstdin);
  }

  /**
   * @brief Start a process.
   *
   * Calls do_open( @a file , @a argv , @a mode|pstdin ).
   *
   * @param file  a string containing the pathname of a program to execute.
   * @param argv  a vector of argument strings passed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, const argv_type&, pmode)
   */
  void
  open( const std::string& file,
        const argv_type& argv,
        pmode mode = pstdin) {
    this->do_open(file, argv, mode|pstdin);
  }
};


/**
 * @class basic_pstream
 * @brief Class template for Bidirectional PStreams.
 *
 * Writing to a pstream opened with @c pmode @c pstdin writes to the
 * standard input of the command.
 * Reading from a pstream opened with @c pmode @c pstdout and/or @c pstderr
 * reads the command's standard output and/or standard error.
 * Any of the process' @c stdin, @c stdout or @c stderr that is not
 * connected to the pstream (as specified by the @c pmode)
 * will be the same as the process that created the pstream object,
 * unless altered by the command itself.
 */
template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_pstream
      : public std::basic_iostream<CharT, Traits>
      , public pstream_common<CharT, Traits>
      , virtual public pstreams {
  typedef std::basic_iostream<CharT, Traits>    iostream_type;
  typedef pstream_common<CharT, Traits>         pbase_type;

  using pbase_type::buf_;  // declare name in this scope

public:
  /// Type used to specify how to connect to the process.
  typedef typename pbase_type::pmode            pmode;

  /// Type used to hold the arguments for a command.
  typedef typename pbase_type::argv_type        argv_type;

  /// Default constructor, creates an uninitialised stream.
  basic_pstream()
      : iostream_type(NULL), pbase_type() { }

  /**
   * @brief Constructor that initialises the stream by starting a process.
   *
   * Initialises the stream buffer by calling do_open() with the supplied
   * arguments.
   *
   * @param command  a string containing a shell command.
   * @param mode     the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, pmode)
   */
  basic_pstream( const std::string& command,
                 pmode mode = pstdout|pstdin )
      : iostream_type(NULL), pbase_type(command, mode) { }

  /**
   * @brief Constructor that initialises the stream by starting a process.
   *
   * Initialises the stream buffer by calling do_open() with the supplied
   * arguments.
   *
   * @param file  a string containing the pathname of a program to execute.
   * @param argv  a vector of argument strings passed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, const argv_type&, pmode)
   */
  basic_pstream( const std::string& file,
                 const argv_type& argv,
                 pmode mode = pstdout|pstdin )
      : iostream_type(NULL), pbase_type(file, argv, mode) { }

  /**
   * @brief Destructor
   *
   * Closes the stream and waits for the child to exit.
   */
  ~basic_pstream() { }

  /**
   * @brief Start a process.
   *
   * Calls do_open( @a %command , @a mode ).
   *
   * @param command  a string containing a shell command.
   * @param mode     the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, pmode)
   */
  void
  open(const std::string& command, pmode mode = pstdout|pstdin) {
    this->do_open(command, mode);
  }

  /**
   * @brief Start a process.
   *
   * Calls do_open( @a file , @a argv , @a mode ).
   *
   * @param file  a string containing the pathname of a program to execute.
   * @param argv  a vector of argument strings passed to the new program.
   * @param mode  the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, const argv_type&, pmode)
   */
  void
  open( const std::string& file,
        const argv_type& argv,
        pmode mode = pstdout|pstdin ) {
    this->do_open(file, argv, mode);
  }

  /**
   * @brief Set streambuf to read from process' @c stdout.
   * @return  @c *this
   */
  basic_pstream&
  out() {
    this->buf_.read_err(false);
    return *this;
  }

  /**
   * @brief Set streambuf to read from process' @c stderr.
   * @return  @c *this
   */
  basic_pstream&
  err() {
    this->buf_.read_err(true);
    return *this;
  }
};


/**
 * @class basic_rpstream
 * @brief template for Restricted PStreams.
 *
 * Writing to an rpstream opened with @c pmode @c pstdin writes to the
 * standard input of the command.
 * It is not possible to read directly from an rpstream object, to use
 * an rpstream as in istream you must call either basic_rpstream::out()
 * or basic_rpstream::err(). This is to prevent accidental reads from
 * the wrong input source. If the rpstream was not opened with @c pmode
 * @c pstderr then the class cannot read the process' @c stderr, and
 * basic_rpstream::err() will return an istream that reads from the
 * process' @c stdout, and vice versa.
 * Reading from an rpstream opened with @c pmode @c pstdout and/or
 * @c pstderr reads the command's standard output and/or standard error.
 * Any of the process' @c stdin, @c stdout or @c stderr that is not
 * connected to the pstream (as specified by the @c pmode)
 * will be the same as the process that created the pstream object,
 * unless altered by the command itself.
 */

template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_rpstream
      : public std::basic_ostream<CharT, Traits>
      , private std::basic_istream<CharT, Traits>
      , private pstream_common<CharT, Traits>
      , virtual public pstreams {
  typedef std::basic_ostream<CharT, Traits>     ostream_type;
  typedef std::basic_istream<CharT, Traits>     istream_type;
  typedef pstream_common<CharT, Traits>         pbase_type;

  using pbase_type::buf_;  // declare name in this scope

public:
  /// Type used to specify how to connect to the process.
  typedef typename pbase_type::pmode            pmode;

  /// Type used to hold the arguments for a command.
  typedef typename pbase_type::argv_type        argv_type;

  /// Default constructor, creates an uninitialised stream.
  basic_rpstream()
      : ostream_type(NULL) , istream_type(NULL) , pbase_type() { }

  /**
   * @brief  Constructor that initialises the stream by starting a process.
   *
   * Initialises the stream buffer by calling do_open() with the supplied
   * arguments.
   *
   * @param command a string containing a shell command.
   * @param mode the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, pmode)
   */
  basic_rpstream(const std::string& command, pmode mode = pstdout|pstdin)
      : ostream_type(NULL) , istream_type(NULL) , pbase_type(command, mode) { }

  /**
   * @brief  Constructor that initialises the stream by starting a process.
   *
   * Initialises the stream buffer by calling do_open() with the supplied
   * arguments.
   *
   * @param file a string containing the pathname of a program to execute.
   * @param argv a vector of argument strings passed to the new program.
   * @param mode the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, const argv_type&, pmode)
   */
  basic_rpstream( const std::string& file,
                  const argv_type& argv,
                  pmode mode = pstdout|pstdin )
      : ostream_type(NULL), istream_type(NULL), pbase_type(file, argv, mode) { }

  /// Destructor
  ~basic_rpstream() { }

  /**
   * @brief  Start a process.
   *
   * Calls do_open( @a %command , @a mode ).
   *
   * @param command a string containing a shell command.
   * @param mode the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, pmode)
   */
  void
  open(const std::string& command, pmode mode = pstdout|pstdin) {
    this->do_open(command, mode);
  }

  /**
   * @brief  Start a process.
   *
   * Calls do_open( @a file , @a argv , @a mode ).
   *
   * @param file a string containing the pathname of a program to execute.
   * @param argv a vector of argument strings passed to the new program.
   * @param mode the I/O mode to use when opening the pipe.
   * @see   do_open(const std::string&, const argv_type&, pmode)
   */
  void
  open( const std::string& file,
        const argv_type& argv,
        pmode mode = pstdout|pstdin ) {
    this->do_open(file, argv, mode);
  }

  /**
   * @brief  Obtain a reference to the istream that reads
   *         the process' @c stdout.
   * @return @c *this
   */
  istream_type&
  out() {
    this->buf_.read_err(false);
    return *this;
  }

  /**
   * @brief  Obtain a reference to the istream that reads
   *         the process' @c stderr.
   * @return @c *this
   */
  istream_type&
  err() {
    this->buf_.read_err(true);
    return *this;
  }
};


/// Type definition for common template specialisation.
typedef basic_pstreambuf<char> pstreambuf;
/// Type definition for common template specialisation.
typedef basic_ipstream<char> ipstream;
/// Type definition for common template specialisation.
typedef basic_opstream<char> opstream;
/// Type definition for common template specialisation.
typedef basic_pstream<char> pstream;
/// Type definition for common template specialisation.
typedef basic_rpstream<char> rpstream;


/**
 * When inserted into an output pstream the manipulator calls
 * basic_pstreambuf<C,T>::peof() to close the output pipe,
 * causing the child process to receive the End Of File indicator
 * on subsequent reads from its @c stdin stream.
 *
 * @brief   Manipulator to close the pipe connected to the process' stdin.
 * @param   s  An output PStream class.
 * @return  The stream object the manipulator was invoked on.
 * @see     basic_pstreambuf<C,T>::peof()
 * @relates basic_opstream basic_pstream basic_rpstream
 */
template <typename C, typename T>
inline std::basic_ostream<C,T>&
peof(std::basic_ostream<C,T>& s) {
  typedef basic_pstreambuf<C,T> pstreambuf;
  if (pstreambuf* p = dynamic_cast<pstreambuf*>(s.rdbuf()))
    p->peof();
  return s;
}


/*
 * member definitions for pstreambuf
 */


/**
 * @class basic_pstreambuf
 * Provides underlying streambuf functionality for the PStreams classes.
 */

/** Creates an uninitialised stream buffer. */
template <typename C, typename T>
inline
basic_pstreambuf<C,T>::basic_pstreambuf()
    : ppid_(-1)   // initialise to -1 to indicate no process run yet.
    , wpipe_(-1)
    , wbuffer_(NULL)
    , rsrc_(rsrc_out)
    , status_(-1)
    , error_(0) {
  init_rbuffers();
}

/**
 * Initialises the stream buffer by calling open() with the supplied
 * arguments.
 *
 * @param command a string containing a shell command.
 * @param mode    the I/O mode to use when opening the pipe.
 * @see   open()
 */
template <typename C, typename T>
inline
basic_pstreambuf<C,T>::basic_pstreambuf(const std::string& command, pmode mode)
    : ppid_(-1)   // initialise to -1 to indicate no process run yet.
    , wpipe_(-1)
    , wbuffer_(NULL)
    , rsrc_(rsrc_out)
    , status_(-1)
    , error_(0) {
  init_rbuffers();
  open(command, mode);
}

/**
 * Initialises the stream buffer by calling open() with the supplied
 * arguments.
 *
 * @param file  a string containing the name of a program to execute.
 * @param argv  a vector of argument strings passsed to the new program.
 * @param mode  the I/O mode to use when opening the pipe.
 * @see   open()
 */
template <typename C, typename T>
inline
basic_pstreambuf<C,T>::basic_pstreambuf( const std::string& file,
    const argv_type& argv,
    pmode mode )
    : ppid_(-1)   // initialise to -1 to indicate no process run yet.
    , wpipe_(-1)
    , wbuffer_(NULL)
    , rsrc_(rsrc_out)
    , status_(-1)
    , error_(0) {
  init_rbuffers();
  open(file, argv, mode);
}

/**
 * Closes the stream by calling close().
 * @see close()
 */
template <typename C, typename T>
inline
basic_pstreambuf<C,T>::~basic_pstreambuf() {
  close();
}

/**
 * Starts a new process by passing @a command to the shell
 * and opens pipes to the process with the specified @a mode.
 *
 * There is no way to tell whether the shell command succeeded, this
 * function will always succeed unless resource limits (such as
 * memory usage, or number of processes or open files) are exceeded.
 *
 * @param   command  a string containing a shell command.
 * @param   mode     a bitwise OR of one or more of @c out, @c in, @c err.
 * @return  NULL if the shell could not be started or the
 *          pipes could not be opened, @c this otherwise.
 * @see     execlp(3)
 */
template <typename C, typename T>
basic_pstreambuf<C,T>*
basic_pstreambuf<C,T>::open(const std::string& command, pmode mode) {
  basic_pstreambuf<C,T>* ret = NULL;

  if (!is_open()) {
    switch (fork(mode)) {
    case 0 : {
      // this is the new process, exec command
      ::execlp("sh", "sh", "-c", command.c_str(), (void*)NULL);

      // can only reach this point if exec() failed

      // parent can get exit code from waitpid()
      ::_exit(errno);
      // using std::exit() would make static dtors run twice
    }
    case -1 : {
      // couldn't fork, error already handled in pstreambuf::fork()
      break;
    }
    default : {
      // this is the parent process
      // activate buffers
      create_buffers(mode);
      ret = this;
    }
    }
  }
  return ret;
}

/**
 * Starts a new process by executing @a file with the arguments in
 * @a argv and opens pipes to the process with the specified @a mode.
 *
 * By convention argv[0] should be the file name of the file being executed.
 * Will duplicate the actions of  the  shell  in searching for an
 * executable file if the specified file name does not contain a slash (/)
 * character.
 *
 * @param   file  a string containing the pathname of a program to execute.
 * @param   argv  a vector of argument strings passed to the new program.
 * @param   mode  a bitwise OR of one or more of @c out, @c in and @c err.
 * @return  NULL if a pipe could not be opened or if the program could
 *          not be executed, @c this otherwise.
 * @see     execvp(3)
 */
template <typename C, typename T>
basic_pstreambuf<C,T>*
basic_pstreambuf<C,T>::open( const std::string& file,
                             const argv_type& argv,
                             pmode mode ) {
  basic_pstreambuf<C,T>* ret = NULL;

  if (!is_open()) {
    switch (fork(mode)) {
    case 0 : {
      // this is the new process, exec command

      char** arg_v = new char*[argv.size()+1];
      for (std::size_t i = 0; i < argv.size(); ++i) {
        const std::string& src = argv[i];
        char*& dest = arg_v[i];
        dest = new char[src.size()+1];
        dest[ src.copy(dest, src.size()) ] = '\0';
      }
      arg_v[argv.size()] = NULL;

      ::execvp(file.c_str(), arg_v);

      // can only reach this point if exec() failed

      // parent can get exit code from waitpid()
      ::_exit(errno);
      // using std::exit() would make static dtors run twice
    }
    case -1 : {
      // couldn't fork, error already handled in pstreambuf::fork()
      break;
    }
    default : {
      // this is the parent process
      // activate buffers
      create_buffers(mode);
      ret = this;
    }
    }
  }
  return ret;
}

/**
 * @brief  Helper function to close an array of file descriptors.
 *
 * Inspects each of the @a count file descriptors in the array @a filedes
 * and calls @c close() if they have a non-negative value.
 * @param filedes an array of file descriptors
 * @param count size of the array.
 * @relates basic_pstreambuf
 */
inline void
close_fd_array(int* filedes, std::size_t count) {
  for (std::size_t i = 0; i < count; ++i)
    if (filedes[i] >= 0)
      if (::close(filedes[i]) == 0)
        filedes[i] = -1;
}

/**
 * Creates pipes as specified by @a mode and calls @c fork() to create
 * a new process. If the fork is successful the parent process stores
 * the child's PID and the opened pipes and the child process replaces
 * its standard streams with the opened pipes.
 *
 * If an error occurs the error code will be set to one of the possile
 * errors for @c pipe() or @c fork().
 * See your system's documentation for these error codes.
 *
 * @param   mode  an OR of pmodes specifying which of the child's
 *                standard streams to connect to.
 * @return  On success the PID of the child is returned in the parent's
 *          context and zero is returned in the child's context.
 *          On error -1 is returned and the error code is set appropriately.
 */
template <typename C, typename T>
pid_t
basic_pstreambuf<C,T>::fork(pmode mode) {
  pid_t pid = -1;

  // three pairs of file descriptors, for pipes connected to the
  // process' stdin, stdout and stderr
  // (stored in a single array so close_fd_array() can close all at once)
  fd_type fd[6] =  {-1, -1, -1, -1, -1, -1};
  fd_type* const pin = fd;
  fd_type* const pout = fd+2;
  fd_type* const perr = fd+4;

  // constants for read/write ends of pipe
  const int RD = 0;
  const int WR = 1;

  // N.B.
  // For the pstreambuf pin is an output stream and
  // pout and perr are input streams.

  if (!error_ && mode&pstdin && ::pipe(pin))
    error_ = errno;

  if (!error_ && mode&pstdout && ::pipe(pout))
    error_ = errno;

  if (!error_ && mode&pstderr && ::pipe(perr))
    error_ = errno;

  if (!error_) {
    pid = ::fork();
    switch (pid) {
    case 0 : {
      // this is the new process

      // for each open pipe close one end and redirect the
      // respective standard stream to the other end

      if (*pin >= 0) {
        ::close(pin[WR]);
        ::dup2(pin[RD], STDIN_FILENO);
        ::close(pin[RD]);
      }
      if (*pout >= 0) {
        ::close(pout[RD]);
        ::dup2(pout[WR], STDOUT_FILENO);
        ::close(pout[WR]);
      }
      if (*perr >= 0) {
        ::close(perr[RD]);
        ::dup2(perr[WR], STDERR_FILENO);
        ::close(perr[WR]);
      }
      break;
    }
    case -1 : {
      // couldn't fork for some reason
      error_ = errno;
      // close any open pipes
      close_fd_array(fd, 6);
      break;
    }
    default : {
      // this is the parent process, store process' pid
      ppid_ = pid;

      // store one end of open pipes and close other end
      if (*pin >= 0) {
        wpipe_ = pin[WR];
        ::close(pin[RD]);
      }
      if (*pout >= 0) {
        rpipe_[rsrc_out] = pout[RD];
        ::close(pout[WR]);
      }
      if (*perr >= 0) {
        rpipe_[rsrc_err] = perr[RD];
        ::close(perr[WR]);
      }

      if (rpipe_[rsrc_out] == -1 && rpipe_[rsrc_err] >= 0) {
        // reading stderr but not stdout, so use stderr for all reads
        read_err(true);
      }
    }
    }
  } else {
    // close any pipes we opened before failure
    close_fd_array(fd, 6);
  }
  return pid;
}

/**
 * Closes all pipes and calls wait() to wait for the process to finish.
 * If an error occurs the error code will be set to one of the possible
 * errors for @c waitpid().
 * See your system's documentation for these errors.
 *
 * @return  @c this on successful close or @c NULL if there is no
 *          process to close or if an error occurs.
 */
template <typename C, typename T>
basic_pstreambuf<C,T>*
basic_pstreambuf<C,T>::close() {
  basic_pstreambuf<C,T>* ret = NULL;
  if (is_open()) {
    sync();

    destroy_buffers(pstdin|pstdout|pstderr);

    // close pipes before wait() so child gets EOF/SIGPIPE
    close_fd_array(&wpipe_, 1);
    close_fd_array(rpipe_, 2);

    if (wait() == 1) {
      ret = this;
    }
  }
  return ret;
}

/**
 *  Called on construction to initialise the arrays used for reading.
 */
template <typename C, typename T>
inline void
basic_pstreambuf<C,T>::init_rbuffers() {
  rpipe_[rsrc_out] = rpipe_[rsrc_err] = -1;
  rbuffer_[rsrc_out] = rbuffer_[rsrc_err] = NULL;
  rbufstate_[0] = rbufstate_[1] = rbufstate_[2] = NULL;
}

template <typename C, typename T>
void
basic_pstreambuf<C,T>::create_buffers(pmode mode) {
  if (mode & pstdin) {
    delete[] wbuffer_;
    wbuffer_ = new char_type[bufsz];
    this->setp(wbuffer_, wbuffer_ + bufsz);
  }
  if (mode & pstdout) {
    delete[] rbuffer_[rsrc_out];
    rbuffer_[rsrc_out] = new char_type[bufsz];
    if (rsrc_ == rsrc_out)
      this->setg(rbuffer_[rsrc_out] + pbsz, rbuffer_[rsrc_out] + pbsz,
                 rbuffer_[rsrc_out] + pbsz);
  }
  if (mode & pstderr) {
    delete[] rbuffer_[rsrc_err];
    rbuffer_[rsrc_err] = new char_type[bufsz];
    if (rsrc_ == rsrc_err)
      this->setg(rbuffer_[rsrc_err] + pbsz, rbuffer_[rsrc_err] + pbsz,
                 rbuffer_[rsrc_err] + pbsz);
  }
}

template <typename C, typename T>
void
basic_pstreambuf<C,T>::destroy_buffers(pmode mode) {
  if (mode & pstdin) {
    this->setp(NULL, NULL);
    delete[] wbuffer_;
    wbuffer_ = NULL;
  }
  if (mode & pstdout) {
    if (rsrc_ == rsrc_out)
      this->setg(NULL, NULL, NULL);
    delete[] rbuffer_[rsrc_out];
    rbuffer_[rsrc_out] = NULL;
  }
  if (mode & pstderr) {
    if (rsrc_ == rsrc_err)
      this->setg(NULL, NULL, NULL);
    delete[] rbuffer_[rsrc_err];
    rbuffer_[rsrc_err] = NULL;
  }
}

template <typename C, typename T>
typename basic_pstreambuf<C,T>::buf_read_src
basic_pstreambuf<C,T>::switch_read_buffer(buf_read_src src) {
  if (rsrc_ != src) {
    char_type* tmpbufstate[] = {this->eback(), this->gptr(), this->egptr()};
    this->setg(rbufstate_[0], rbufstate_[1], rbufstate_[2]);
    for (std::size_t i = 0; i < 3; ++i)
      rbufstate_[i] = tmpbufstate[i];
    rsrc_ = src;
  }
  return rsrc_;
}

/**
 * Suspends execution and waits for the associated process to exit, or
 * until a signal is delivered whose action is to terminate the current
 * process or to call a signal handling function. If the process has
 * already exited wait() returns immediately.
 *
 * @param   nohang  true to return immediately if the process has not exited.
 * @return  1 if the process has exited.
 *          0 if @a nohang is true and the process has not exited yet.
 *          -1 if no process has been started or if an error occurs,
 *          in which case the error can be found using error().
 */
template <typename C, typename T>
int
basic_pstreambuf<C,T>::wait(bool nohang) {
  int exited = -1;
  if (is_open()) {
    int status;
    switch (::waitpid(ppid_, &status, nohang ? WNOHANG : 0)) {
    case 0 :
      // nohang was true and process has not exited
      exited = 0;
      break;
    case -1 :
      error_ = errno;
      break;
    default :
      // process has exited
      ppid_ = 0;
      status_ = status;
      exited = 1;
      destroy_buffers(pstdin|pstdout|pstderr);
      close_fd_array(&wpipe_, 1);
      close_fd_array(rpipe_, 2);
      break;
    }
  }
  return exited;
}

/**
 * Sends the specified signal to the process.  A signal can be used to
 * terminate a child process that would not exit otherwise.
 *
 * If an error occurs the error code will be set to one of the possible
 * errors for @c kill().  See your system's documentation for these errors.
 *
 * @param   signal  A signal to send to the child process.
 * @return  @c this or @c NULL if @c kill() fails.
 */
template <typename C, typename T>
inline basic_pstreambuf<C,T>*
basic_pstreambuf<C,T>::kill(int signal) {
  basic_pstreambuf<C,T>* ret = NULL;
  if (is_open()) {
    if (::kill(ppid_, signal))
      error_ = errno;
    else {
      // TODO call exited() to check for exit and clean up? leave to user?
      ret = this;
    }
  }
  return ret;
}

/**
 *  @return  True if the associated process has exited, false otherwise.
 *  @see     basic_pstreambuf<C,T>::close()
 */
template <typename C, typename T>
inline bool
basic_pstreambuf<C,T>::exited() {
  return ppid_ == 0 || wait(true)==1;
}


/**
 *  @return  The exit status of the child process, or -1 if close()
 *           has not yet been called to wait for the child to exit.
 *  @see     basic_pstreambuf<C,T>::close()
 */
template <typename C, typename T>
inline int
basic_pstreambuf<C,T>::status() const {
  return status_;
}

/**
 *  @return  The error code of the most recently failed operation, or zero.
 */
template <typename C, typename T>
inline int
basic_pstreambuf<C,T>::error() const {
  return error_;
}

/**
 *  Closes the output pipe, causing the child process to receive the
 *  End Of File indicator on subsequent reads from its @c stdin stream.
 */
template <typename C, typename T>
inline void
basic_pstreambuf<C,T>::peof() {
  sync();
  destroy_buffers(pstdin);
  close_fd_array(&wpipe_, 1);
}

/**
 * @return  true if a previous call to open() succeeded and wait() has
 *          not been called and determined that the process has exited,
 *          false otherwise.
 * @warning This function can not be used to determine whether the
 *          command used to initialise the buffer was successfully
 *          executed or not. If the shell command failed this function
 *          will still return true.
 *          You can use exited() to see if it's still open.
 */
template <typename C, typename T>
inline bool
basic_pstreambuf<C,T>::is_open() const {
  return ppid_ > 0;
}

/**
 * Toggle the stream used for reading. If @a readerr is @c true then the
 * process' @c stderr output will be used for subsequent extractions, if
 * @a readerr is false the the process' stdout will be used.
 * @param   readerr  @c true to read @c stderr, @c false to read @c stdout.
 * @return  @c true if the requested stream is open and will be used for
 *          subsequent extractions, @c false otherwise.
 */
template <typename C, typename T>
inline bool
basic_pstreambuf<C,T>::read_err(bool readerr) {
  buf_read_src src = readerr ? rsrc_err : rsrc_out;
  if (rpipe_[src]>=0) {
    switch_read_buffer(src);
    return true;
  }
  return false;
}

/**
 * Called when the internal character buffer is not present or is full,
 * to transfer the buffer contents to the pipe. For unbuffered streams
 * this is called for every insertion.
 *
 * @param   c  a character to be written to the pipe.
 * @return  @c traits_type::not_eof(c) if @a c is equal to @c
 *          traits_type::eof(). Otherwise returns @a c if @a c can be
 *          written to the pipe, or @c traits_type::eof() if not.
 */
template <typename C, typename T>
typename basic_pstreambuf<C,T>::int_type
basic_pstreambuf<C,T>::overflow(int_type c) {
  if (!empty_buffer())
    return traits_type::eof();
  else if (!traits_type::eq_int_type(c, traits_type::eof()))
    return this->sputc(c);
  else
    return traits_type::not_eof(c);
}


template <typename C, typename T>
int
basic_pstreambuf<C,T>::sync() {
  return !exited() && empty_buffer() ? 0 : -1;
}

/**
 * @param   s  character buffer.
 * @param   n  buffer length.
 * @return  the number of characters written.
 */
template <typename C, typename T>
std::streamsize
basic_pstreambuf<C,T>::xsputn(const char_type* s, std::streamsize n) {
  if (n < this->epptr() - this->pptr()) {
    std::memcpy(this->pptr(), s, n * sizeof(char_type));
    this->pbump(n);
    return n;
  } else {
    for (std::streamsize i = 0; i < n; ++i) {
      if (traits_type::eq_int_type(this->sputc(s[i]), traits_type::eof()))
        return i;
    }
    return n;
  }
}

/**
 * @return  true if the buffer was emptied, false otherwise.
 */
template <typename C, typename T>
bool
basic_pstreambuf<C,T>::empty_buffer() {
  const std::streamsize count = this->pptr() - this->pbase();
  const std::streamsize written = this->write(this->wbuffer_, count);
  if (count > 0 && written == count) {
    this->pbump(-written);
    return true;
  }
  return false;
}

/**
 * Called when the internal character buffer is is empty, to re-fill it
 * from the pipe.
 *
 * @return The first available character in the buffer,
 * or @c traits_type::eof() in case of failure.
 * @see uflow()
 */
template <typename C, typename T>
typename basic_pstreambuf<C,T>::int_type
basic_pstreambuf<C,T>::underflow() {
  if (this->gptr() < this->egptr() || fill_buffer())
    return traits_type::to_int_type(*this->gptr());
  else
    return traits_type::eof();
}

/**
 * Attempts to make @a c available as the next character to be read by
 * @c sgetc().
 *
 * @param   c   a character to make available for extraction.
 * @return  @a c if the character can be made available,
 *          @c traits_type::eof() otherwise.
 */
template <typename C, typename T>
typename basic_pstreambuf<C,T>::int_type
basic_pstreambuf<C,T>::pbackfail(int_type c) {
  if (this->gptr() != this->eback()) {
    this->gbump(-1);
    if (!traits_type::eq_int_type(c, traits_type::eof()))
      *this->gptr() = traits_type::to_char_type(c);
    return traits_type::not_eof(c);
  } else
    return traits_type::eof();
}

/**
 * @return  true if the buffer was filled, false otherwise.
 */
template <typename C, typename T>
bool
basic_pstreambuf<C,T>::fill_buffer() {
  const std::streamsize pb1 = this->gptr() - this->eback();
  const std::streamsize pb2 = pbsz;
  const std::streamsize npb = std::min(pb1, pb2);

  std::memmove( rbuffer() + pbsz - npb,
                this->gptr() - npb,
                npb * sizeof(char_type) );

  const std::streamsize rc = read(rbuffer() + pbsz, bufsz - pbsz);

  if (rc > 0) {
    this->setg( rbuffer() + pbsz - npb,
                rbuffer() + pbsz,
                rbuffer() + pbsz + rc );
    return true;
  } else {
    this->setg(NULL, NULL, NULL);
    return false;
  }
}

/**
 * Writes up to @a n characters to the pipe from the buffer @a s.
 * This currently only works for fixed width character encodings where
 * each character uses sizeof(char_type) bytes.
 *
 * @param   s  character buffer.
 * @param   n  buffer length.
 * @return  the number of characters written.
 */
template <typename C, typename T>
inline std::streamsize
basic_pstreambuf<C,T>::write(char_type* s, std::streamsize n) {
  return wpipe() >= 0 ? ::write(wpipe(), s, n * sizeof(char_type)) : 0;
}

/**
 * Reads up to @a n characters from the pipe to the buffer @a s.
 * This currently only works for fixed width character encodings where
 * each character uses sizeof(char_type) bytes.
 *
 * @param   s  character buffer.
 * @param   n  buffer length.
 * @return  the number of characters read.
 */
template <typename C, typename T>
inline std::streamsize
basic_pstreambuf<C,T>::read(char_type* s, std::streamsize n) {
  return rpipe() >= 0 ? ::read(rpipe(), s, n * sizeof(char_type)) : 0;
}

/** @return a reference to the output file descriptor */
template <typename C, typename T>
inline typename basic_pstreambuf<C,T>::fd_type&
basic_pstreambuf<C,T>::wpipe() {
  return wpipe_;
}

/** @return a reference to the active input file descriptor */
template <typename C, typename T>
inline typename basic_pstreambuf<C,T>::fd_type&
basic_pstreambuf<C,T>::rpipe() {
  return rpipe_[rsrc_];
}

/** @return a reference to the specified input file descriptor */
template <typename C, typename T>
inline typename basic_pstreambuf<C,T>::fd_type&
basic_pstreambuf<C,T>::rpipe(buf_read_src which) {
  return rpipe_[which];
}

/** @return a pointer to the start of the active input buffer area. */
template <typename C, typename T>
inline typename basic_pstreambuf<C,T>::char_type*
basic_pstreambuf<C,T>::rbuffer() {
  return rbuffer_[rsrc_];
}


/*
 * member definitions for pstream_common
 */

/**
 * @class pstream_common
 * Abstract Base Class providing common functionality for basic_ipstream,
 * basic_opstream and basic_pstream.
 * pstream_common manages the basic_pstreambuf stream buffer that is used
 * by the derived classes to initialise an IOStream class.
 */

/** Creates an uninitialised stream. */
template <typename C, typename T>
inline
pstream_common<C,T>::pstream_common()
    : std::basic_ios<C,T>(NULL)
    , command_()
    , buf_() {
  this->init(&buf_);
}

/**
 * Initialises the stream buffer by calling
 * do_open( @a command , @a mode )
 *
 * @param command a string containing a shell command.
 * @param mode    the I/O mode to use when opening the pipe.
 * @see   do_open(const std::string&, pmode)
 */
template <typename C, typename T>
inline
pstream_common<C,T>::pstream_common(const std::string& command, pmode mode)
    : std::basic_ios<C,T>(NULL)
    , command_(command)
    , buf_() {
  this->init(&buf_);
  do_open(command, mode);
}

/**
 * Initialises the stream buffer by calling
 * do_open( @a file , @a argv , @a mode )
 *
 * @param file  a string containing the pathname of a program to execute.
 * @param argv  a vector of argument strings passed to the new program.
 * @param mode  the I/O mode to use when opening the pipe.
 * @see do_open(const std::string&, const argv_type&, pmode)
 */
template <typename C, typename T>
inline
pstream_common<C,T>::pstream_common( const std::string& file,
                                     const argv_type& argv,
                                     pmode mode )
    : std::basic_ios<C,T>(NULL)
    , command_(file)
    , buf_() {
  this->init(&buf_);
  do_open(file, argv, mode);
}

/**
 * This is a pure virtual function to make @c pstream_common abstract.
 * Because it is the destructor it will be called by derived classes
 * and so must be defined.  It is also protected, to discourage use of
 * the PStreams classes through pointers or references to the base class.
 *
 * @sa If defining a pure virtual seems odd you should read
 * http://www.gotw.ca/gotw/031.htm (and the rest of the site as well!)
 */
template <typename C, typename T>
inline
pstream_common<C,T>::~pstream_common() {}

/**
 * Calls rdbuf()->open( @a command , @a mode )
 * and sets @c failbit on error.
 *
 * @param command a string containing a shell command.
 * @param mode    the I/O mode to use when opening the pipe.
 * @see   basic_pstreambuf::open(const std::string&, pmode)
 */
template <typename C, typename T>
inline void
pstream_common<C,T>::do_open(const std::string& command, pmode mode) {
  if (!buf_.open((command_=command), mode))
    this->setstate(std::ios_base::failbit);
}

/**
 * Calls rdbuf()->open( @a file, @a  argv, @a mode )
 * and sets @c failbit on error.
 *
 * @param file  a string containing the pathname of a program to execute.
 * @param argv  a vector of argument strings passed to the new program.
 * @param mode  the I/O mode to use when opening the pipe.
 * @see   basic_pstreambuf::open(const std::string&, const argv_type&, pmode)
 */
template <typename C, typename T>
inline void
pstream_common<C,T>::do_open( const std::string& file,
                              const argv_type& argv,
                              pmode mode ) {
  if (!buf_.open((command_=file), argv, mode))
    this->setstate(std::ios_base::failbit);
}

/** Calls rdbuf->close() and sets @c failbit on error. */
template <typename C, typename T>
inline void
pstream_common<C,T>::close() {
  if (!buf_.close())
    this->setstate(std::ios_base::failbit);
}

/**
 * @return  rdbuf()->is_open().
 * @see     basic_pstreambuf::is_open()
 */
template <typename C, typename T>
inline bool
pstream_common<C,T>::is_open() const {
  return buf_.is_open();
}

/** @return a string containing the command used to initialise the stream. */
template <typename C, typename T>
inline const std::string&
pstream_common<C,T>::command() const {
  return command_;
}

/** @return a pointer to the private stream buffer member. */
// TODO  document behaviour if buffer replaced.
template <typename C, typename T>
inline typename pstream_common<C,T>::streambuf_type*
pstream_common<C,T>::rdbuf() const {
  return const_cast<streambuf_type*>(&buf_);
}


#if REDI_EVISCERATE_PSTREAMS
/**
 * @def REDI_EVISCERATE_PSTREAMS
 * If this macro has a non-zero value then certain internals of the
 * @c basic_pstreambuf template class are exposed. In general this is
 * a Bad Thing, as the internal implementation is largely undocumented
 * and may be subject to change at any time, so this feature is only
 * provided because it might make PStreams useful in situations where
 * it is necessary to do Bad Things.
 */

/**
 * @warning  This function exposes the internals of the stream buffer and
 *           should be used with caution. It is the caller's responsibility
 *           to flush streams etc. in order to clear any buffered data.
 *           The POSIX.1 function @c fdopen(3) is used to obtain the
 *           @c FILE pointers from the streambuf's private file descriptor
 *           members so consult your system's documentation for @c fdopen().
 *
 * @param   in    A FILE* that will refer to the process' stdin.
 * @param   out   A FILE* that will refer to the process' stdout.
 * @param   err   A FILE* that will refer to the process' stderr.
 * @return  An OR of zero or more of @c pstdin, @c pstdout, @c pstderr.
 *
 * For each open stream shared with the child process a @c FILE* is
 * obtained and assigned to the corresponding parameter. For closed
 * streams @c NULL is assigned to the parameter.
 * The return value can be tested to see which parameters should be
 * @c !NULL by masking with the corresponding @c pmode value.
 */
template <typename C, typename T>
std::size_t
basic_pstreambuf<C,T>::fopen(std::FILE*& in, std::FILE*& out, std::FILE*& err) {
  in = out = err = NULL;
  std::size_t open_files = 0;
  if (wpipe() > -1) {
    if ((in = ::fdopen(wpipe(), "w"))) {
      open_files |= pstdin;
    }
  }
  if (rpipe(rsrc_out) > -1) {
    if ((out = ::fdopen(rpipe(rsrc_out), "r"))) {
      open_files |= pstdout;
    }
  }
  if (rpipe(rsrc_err) > -1) {
    if ((err = ::fdopen(rpipe(rsrc_err), "r"))) {
      open_files |= pstderr;
    }
  }
  return open_files;
}

/**
 *  @warning This function exposes the internals of the stream buffer and
 *  should be used with caution.
 *
 *  @param  in   A FILE* that will refer to the process' stdin.
 *  @param  out  A FILE* that will refer to the process' stdout.
 *  @param  err  A FILE* that will refer to the process' stderr.
 *  @return A bitwise-or of zero or more of @c pstdin, @c pstdout, @c pstderr.
 *  @see    basic_pstreambuf::fopen()
 */
template <typename C, typename T>
inline std::size_t
pstream_common<C,T>::fopen(std::FILE*& in, std::FILE*& out, std::FILE*& err) {
  return buf_.fopen(in, out, err);
}

#endif // REDI_EVISCERATE_PSTREAMS


} // namespace redi

/**
 * @mainpage PStreams Reference
 * @htmlinclude mainpage.html
 */

#endif  // REDI_PSTREAM_H_SEEN

// vim: ts=2 sw=2 expandtab

