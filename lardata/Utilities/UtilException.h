/**
 * \file UtilException.h
 *
 * \ingroup Utilities
 *
 * \brief Class def header for exception classes used in GeometryUtilities
 *
 * @author kazuhiro
 */

/** \addtogroup Util

    @{*/
#ifndef UTILEXCEPTION_HH
#define UTILEXCEPTION_HH

#include <exception>
#include <string>

namespace util {
  /**
     \class UtilException
     Generic (base) exception class
  */
  class UtilException : public std::exception {

  public:
    UtilException(std::string msg = "") : std::exception()
    {
      _msg = "\033[93m";
      _msg += msg;
      _msg += "\033[00m";
    }

    virtual ~UtilException() throw(){};
    virtual const char* what() const throw() { return _msg.c_str(); }

  private:
    std::string _msg;
  };

}
#endif
/** @} */ // end of doxygen group
