/**
 * @brief simplifies the command-line parsing interface in a new application
 *
 */
/*
 * Original Author: Gheorghe Postelnicu
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

#ifndef H_CMD_LINE_INTERFACE_H
#define H_CMD_LINE_INTERFACE_H

#include <iostream>
#include <list>
#include <vector>
#include <string>

#include <stdio.h> // printf
#include <stdlib.h>

#define GMP_ISOPTION(c) ((c) == '-')


// generic class that will hold and handle data for a specific type of option
// when the user wants a certain type of option to be handled
//     they will just have to specify
//     - name
//     - add pointers to the addresses that will be affected 
//       when the option is found
//     - help text that gives a generic description
template <class T>
class CCmdLineOption {
protected:
  std::string       m_strOption;
  std::string       m_strHelp;
  std::vector<T *>  m_pars; // this contains parameters assigned to the option
  int               m_iPars;
  int               m_iCur;
  bool              m_bInput;

public:

  CCmdLineOption() {}
  CCmdLineOption(const char *i_strName, int i_iPars, const char *i_strHelp)
    : m_strOption(), m_strHelp() {
    m_strOption = i_strName;
    m_iPars = i_iPars;
    if ( i_strHelp )
      m_strHelp = i_strHelp;
    else
      m_strHelp = " no help available for this option ";
    m_iCur  = 0;
    m_bInput = false;
  }


  void Add(T *i_par) {
    if ( m_iCur < m_iPars && !m_bInput) {
      m_pars.push_back(i_par);
      m_iCur++;
    } else {
      std::cerr << " ERROR - too many parameters for option " << m_strOption
      << std::endl;
      exit(-1);
    }

    if ( m_iCur == m_iPars ) {
      m_bInput = true;
      m_iCur = 0;
    }
  }

  bool Compare(const std::string i_strOption) const {
    return ( m_strOption.compare(i_strOption) == 0 );
  }

  bool Compare(const char *i_strOption) const {
    std::string strBuf = i_strOption;
    return ( m_strOption.compare(strBuf) == 0 );
  }

  void Set(T val) {
    if ( !m_bInput ) {
      std::cout << " can not set values for option " << m_strOption
      << " until all the pars have been set\n";
    }

    if ( (unsigned)m_iCur < m_pars.size() ) {
      *(m_pars[m_iCur]) = val;
      m_iCur++;
    }
  }


  int  GetArgs() const {
    return m_iPars;
  }

  std::string GetHelp() const {
    char chBuf[10];
    sprintf(chBuf, "%d", m_iPars);
    std::string strBuf = std::string("--") 
      + m_strOption + "\t pars = " + chBuf + "\t" + m_strHelp;
    return strBuf;
  }

  void outStatus() const {
    std::cout << m_strOption ;
    for ( typename std::vector<T*>::const_iterator cit = m_pars.begin();
          cit != m_pars.end(); cit++) {
      std::cout << "\t" << *(*cit) ;
    }
    std::cout << "\n";
  }

};

// specializations to be able to do strong-typing
typedef CCmdLineOption<int> CCmdLineOptionInt;
typedef CCmdLineOption<float> CCmdLineOptionFloat;
typedef CCmdLineOption<std::string> CCmdLineOptionString;

//
// a partial specialization is dedicated to the bool class
// in the case of a boolean - this will act like a flag
// the sole presence of the option will be enough to activate the flag
//
class CCmdLineOptionBool : public CCmdLineOption<bool> {

public:
  CCmdLineOptionBool(const char *i_strName, const char *i_strHelp)
      : CCmdLineOption<bool>(i_strName, 1, i_strHelp) {}

  void Set() {
    *(m_pars[0]) = true;
  }

  int GetArgs() const {
    return 0; // there is no actual parameter involved 
    // - the presence of the option IS  the flag
  }

  std::string GetHelp() const {
    char chBuf[10];
    sprintf(chBuf, "%d", m_iPars);
    std::string strBuf = std::string("--") + m_strOption + "\t" + m_strHelp;
    return strBuf;
  }

  void outStatus() const {
    std::cout << m_strOption << "\t" << *(m_pars[0]) << "\n";
  }

};

class CCmdLineIo {
  std::string *m_pStrItem;
  std::string  m_strHelp;
  bool    m_bOk;
public:
  CCmdLineIo() {
    m_pStrItem = NULL;
    m_bOk = false;
  }
  CCmdLineIo(std::string *i_pStrItem, const char* i_strHelp=NULL) {
    m_pStrItem = i_pStrItem;
    if ( i_strHelp )
      m_strHelp  = i_strHelp;
    else
      m_strHelp = " no help available for IO item";
    m_bOk = false;
  }

  void Set(std::string i_strItem) {
    *m_pStrItem = i_strItem;
    m_bOk = true;
  }

  std::string printHelp() {
    return m_strHelp;
  }

  bool Ok() const {
    return m_bOk;
  }

  std::string GetHelp() const {
    return m_strHelp;
  }

  void outStatus() const {
    std::cout << *m_pStrItem << "\n";
  }
};

//
// class CCmdLineInterface
//
//   this is a placeholder for a list of potential options
//   the redundance in the code is present here because I did not 
//   know how to handle the typing in a generic way
//
//   one must notice that another limitation is that all options 
//   are only allowed to have parameters that have the same type

class CCmdLineInterface {
  std::string   m_strProgName;

  std::list<CCmdLineOptionFloat>   m_lstFloatOption;
  std::list<CCmdLineOptionInt>     m_lstIntOption;
  std::list<CCmdLineOptionString>  m_lstStringOption;
  std::list<CCmdLineOptionBool>    m_lstBoolOption;

  std::list<CCmdLineIo>            m_lstIo;

public:

  CCmdLineInterface(char *i_strProgName) {
    m_strProgName = i_strProgName;
  }

  // the pointers in the following functions HAVE to be 
  //    initialized at the time of the call
  //    this code has been designed so that these pointers 
  //    hold variables addresses
  //    this simplifies the writing and the passage of parameters....
  void AddOptionInt
    (const char *i_cstrName, int *piVal, const char *help=NULL);
  void AddOptionInt
    (const char *i_cstrName, int *piVal1,int *piVal2, const char *help=NULL);
  void AddOptionFloat
    (const char *i_cstrName, float *pfVal, char *help=NULL);
  void AddOptionFloat
    (const char *i_cstrName, float *pfVal1, float *pfVal2, char *help=NULL);
  void AddOptionString
    (const char *i_cstrName, std::string *pstrVal, const char *help=NULL);
  void AddOptionString
    (const char *i_cstrName, std::string *pstrVal_1,
     std::string *pstrVal_2, const char* help=NULL);
  // a bool option is basically just a flag
  // no parameters at all
  void AddOptionBool
    (const char *i_cstrName, bool *pbVal, const char *help=NULL);
  void AddIoItem(std::string *pStrItem, const char *help=NULL) {
    CCmdLineIo item(pStrItem, help);

    m_lstIo.push_back(item);
  }

  //
  // this will do the actual parsing
  //   all the options should of course be registered 
  //   using one of the preceding functions
  //   at the time of this call
  //
  // any option name should be preceded by a '-'
  //
  // the return value is the number of read items in the 
  // command-line string array
  bool Parse(int argc, char *argv[]);

  void PrintHelp() const;

  void print() const;

};

#endif // H_CMD_LINE_INTERFACE_H
