/**
 * @file   SpacePointDumpers.h
 * @brief  Functions dumping space points
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 18th, 2015
 * @see    SpacePointDumpers.cc DumpSpacePoints_module.cc
 */

#ifndef LARDATA_RECOBASE_DUMPERS_SPACEPOINTDUMPERS_H
#define LARDATA_RECOBASE_DUMPERS_SPACEPOINTDUMPERS_H 1

// LArSoft libraries
#include "lardataobj/RecoBase/SpacePoint.h"

// C/C++ standard libraries
#include <type_traits> // std::decay<>
#include <utility>     // std::forward<>()

// --- for the implementation ---
// LArSoft libraries
#include "lardata/ArtDataHelper/Dumpers/NewLine.h"
#include "lardata/ArtDataHelper/Dumpers/hexfloat.h"

namespace recob {
  namespace dumper {

    /// Collection of available printing style options
    struct SpacePointPrintOptions_t {
      IndentOptions_t indent; ///< indentation string
      bool hexFloats = false; ///< print all floating point numbers in base 16

      /**
       * @brief Default constructor
       *
       * By default, the options are:
       *
       *  * no indentation
       *  * same indentation for the first and the following lines
       *  * real numbers printed in base 10
       *
       */
      SpacePointPrintOptions_t() = default;

      SpacePointPrintOptions_t(IndentOptions_t indentOptions, bool bHexFloats)
        : indent(indentOptions), hexFloats(bHexFloats)
      {}

    }; // SpacePointPrintOptions_t

    /**
     * @brief Dumps the content of the specified space point into a stream
     * @tparam Stream the type of the output stream
     * @tparam NewLineRef NewLine reference type (to get a universal reference)
     * @param out the output stream
     * @param sp the space point to be dumped
     * @param options indentation and formatting options
     */
    template <typename Stream, typename NewLineRef = recob::dumper::NewLine<Stream>>
    auto DumpSpacePoint(Stream&& out,
                        recob::SpacePoint const& sp,
                        SpacePointPrintOptions_t const& options = {})
      -> std::enable_if_t<
        std::is_same<NewLine<std::decay_t<Stream>>, std::decay_t<NewLineRef>>::value>;

  } // namespace dumper
} // namespace recob

//==============================================================================
//=== template implementation
//===
//------------------------------------------------------------------------------
//--- recob::dumper::DumpSpacePoint
//---
template <typename Stream, typename NewLineRef>
auto recob::dumper::DumpSpacePoint(Stream&& out,
                                   recob::SpacePoint const& sp,
                                   SpacePointPrintOptions_t const& options /* = {} */
                                   )
  -> std::enable_if_t<std::is_same<NewLine<std::decay_t<Stream>>, std::decay_t<NewLineRef>>::value>
{
  double const* pos = sp.XYZ();
  double const* err = sp.ErrXYZ();

  NewLineRef nl(out, options.indent);
  lar::OptionalHexFloat hexfloat(options.hexFloats);

  nl() << "ID=" << sp.ID() << " at (" << hexfloat(pos[0]) << ", " << hexfloat(pos[1]) << ", "
       << hexfloat(pos[2]) << ") cm, chi^2/NDF=" << hexfloat(sp.Chisq());

  nl() << "variances { x^2=" << hexfloat(err[0]) << " y^2=" << hexfloat(err[2])
       << " z^2=" << hexfloat(err[5]) << " xy=" << hexfloat(err[1]) << " xz=" << hexfloat(err[3])
       << " yz=" << hexfloat(err[4]) << " }";

} // recob::dumper::DumpSpacePoint()

//------------------------------------------------------------------------------

#endif // LARDATA_RECOBASE_DUMPERS_SPACEPOINTDUMPERS_H
