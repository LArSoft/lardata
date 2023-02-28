/**
 * @file   DumpAssociations.h
 * @brief  Dumps on a stream the content of associations
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 25th, 2015
 */

#ifndef UTIL_DUMPASSOCIATIONS_H
#define UTIL_DUMPASSOCIATIONS_H 1

// C//C++ standard libraries
#include <type_traits> // std::is_same<>
#include <typeinfo>

// framework and supporting libraries
#include "canvas/Persistency/Common/Assns.h"
#include "cetlib_except/demangle.h"

namespace util {

  /**
   * @brief Dumps a short introduction about specified association
   * @tparam Stream type of output stream
   * @tparam Left first type in the association
   * @tparam Right second type in the association
   * @tparam Data metadata type in the association
   * @param out output stream
   * @param assns the associations to be dumped
   */
  template <typename Stream, typename Left, typename Right, typename Data>
  void DumpAssociationsIntro(Stream&& out, art::Assns<Left, Right, Data> const& assns)
  {
    out << "Association between '" << cet::demangle_symbol(typeid(Left).name()) << "' and '"
        << cet::demangle_symbol(typeid(Right).name()) << "'";
    if (std::is_same<Data, void>()) {
      out << " with '" << cet::demangle_symbol(typeid(Data).name()) << "' metadata";
    }
    if (assns.size() > 0) { out << " contains " << assns.size() << " relations"; }
    else {
      out << " is empty";
    }
  } // DumpAssociationsIntro<Data>()

} // namespace util

#endif // UTIL_DUMPASSOCIATIONS_H
