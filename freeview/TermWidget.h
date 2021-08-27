/*
 * Original Author: Ruopeng Wang
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
#ifndef TERMWIDGET_H
#define TERMWIDGET_H

#include <ios>
#include <QWidget>
#include <QSocketNotifier>
#include <unistd.h> //Provides STDIN_FILENO

namespace Ui
{
class TermWidget;
}

template< class Elem = char, class Tr = std::char_traits< Elem > >
class StdRedirector : public std::basic_streambuf< Elem, Tr >
{
  /**
    * Callback Function.
    */
  typedef void (*pfncb) ( const Elem*, std::streamsize _Count, void* pUsrData );

public:
  /**
    * Constructor.
    * @param a_Stream the stream to redirect
    * @param a_Cb the callback function
    * @param a_pUsrData user data passed to callback
    */
  StdRedirector( std::ostream& a_Stream, pfncb a_Cb, void* a_pUsrData ) :
    m_Stream( a_Stream ),
    m_pCbFunc( a_Cb ),
    m_pUserData( a_pUsrData )
  {
    //redirect stream
    m_pBuf = m_Stream.rdbuf( this );
  };

  /**
    * Destructor.
    * Restores the original stream.
    */
  ~StdRedirector()
  {
    m_Stream.rdbuf( m_pBuf );
  }

  /**
    * Override xsputn and make it forward data to the callback function.
    */
  std::streamsize xsputn( const Elem* _Ptr, std::streamsize _Count )
  {
    m_pCbFunc( _Ptr, _Count, m_pUserData );
    return _Count;
  }

  /**
    * Override overflow and make it forward data to the callback function.
    */
  typename Tr::int_type overflow( typename Tr::int_type v )
  {
    Elem ch = Tr::to_char_type( v );
    m_pCbFunc( &ch, 1, m_pUserData );
    return Tr::not_eof( v );
  }

protected:
  std::basic_ostream<Elem, Tr>& m_Stream;
  std::streambuf*               m_pBuf;
  pfncb                         m_pCbFunc;
  void*                         m_pUserData;
};

class TermWidget : public QWidget
{
  Q_OBJECT

public:
  explicit TermWidget(QWidget *parent = 0);
  ~TermWidget();

  bool GetDarkTheme();

  bool GetRedirectStdOutput()
  {
    return m_bRedirectStdOutput;
  }

  void EnableListeningStdin();

public slots:
  void SetDarkTheme(bool bDark = true);
  void SetRedirectStdOutput(bool bRedir = true);
  void SetDuplicateStdOutput(bool bDup = true)
  {
    m_bDuplicateStdOutput = bDup;
  }

protected slots:
  void OnCommandTriggered(const QString& cmd);
  void OnTimeOut();
  void OnStdinActivated();

private:
  void ScrollToBottom();
  void AppendErrorString(const QString& strg);

  Ui::TermWidget  *ui;
  StdRedirector<>*  m_stdOut;
  StdRedirector<>*  m_stdErr;
  QString           m_bufferStdOut;
  QString           m_bufferStdErr;
  QString           m_strLogColor;
  QString           m_strErrorColor;
  QSocketNotifier*  m_stdinNotifier;

  bool    m_bRedirectStdOutput;
  bool    m_bDuplicateStdOutput;
};

#endif // TERMWIDGET_H
