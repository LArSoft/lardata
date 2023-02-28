/**
 * @file   PCAxisDumpers.h
 * @brief  Functions dumping principal component axis objects
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 18th, 2015
 * @see    PCAxisDumpers.cc DumpPCAxis_module.cc
 */

#ifndef LARDATA_RECOBASE_DUMPERS_PCAXISDUMPERS_H
#define LARDATA_RECOBASE_DUMPERS_PCAXISDUMPERS_H 1

// LArSoft libraries
#include "lardataobj/RecoBase/PCAxis.h"

// C/C++ standard libraries
#include <iomanip>
#include <string>
#include <type_traits> // std::decay<>

// --- for the implementation ---
// LArSoft libraries
#include "lardata/ArtDataHelper/Dumpers/NewLine.h"

namespace recob {
  namespace dumper {

    /// Dumps the content of the specified PCA axis (indentation info in nl)
    /// @tparam Stream the type of the output stream
    /// @tparam NewLineRef NewLine reference type (to get a universal reference)
    template <typename Stream, typename NewLineRef>
    std::enable_if_t<
      std::is_same<recob::dumper::NewLine<std::decay_t<Stream>>, std::decay_t<NewLineRef>>::value>
    DumpPCAxis(Stream&& out, recob::PCAxis const& pca, NewLineRef&& nl);

    /** ************************************************************************
     * @brief Dumps the content of the specified PCA axis into a stream
     * @tparam Stream the type of output stream
     * @param out the output stream
     * @param pca the principal component axis to be dumped
     * @param indent indentation string (none by default)
     * @param indentFirst whether to indent the first line (yes by default)
     *
     * Insertion operators are required that insert into Stream basic types.
     *
     * This function does not insert a end-of-line after its output.
     */
    template <typename Stream>
    void DumpPCAxis(Stream&& out,
                    recob::PCAxis const& pca,
                    std::string indent = "",
                    bool indentFirst = true)
    {
      DumpPCAxis(std::forward<Stream>(out), pca, makeNewLine(out, indent, !indentFirst));
    }

  } // namespace dumper
} // namespace lar

//==============================================================================
//=== template implementation
//===
//------------------------------------------------------------------------------
//--- recob::dumper::DumpPCAxis
//---
template <typename Stream, typename NewLineRef>
std::enable_if_t<
  std::is_same<recob::dumper::NewLine<std::decay_t<Stream>>, std::decay_t<NewLineRef>>::value>
recob::dumper::DumpPCAxis(Stream&& out, recob::PCAxis const& pca, NewLineRef&& nl)
{

  if (!pca.getSvdOK()) {
    nl() << "<not valid>";
    return;
  }

  nl() << std::setiosflags(std::ios::fixed) << std::setprecision(2) << " ID " << pca.getID()
       << " run on " << pca.getNumHitsUsed() << " space points";
  nl() << "  - center position: " << std::setw(6) << pca.getAvePosition()[0] << ", "
       << pca.getAvePosition()[1] << ", " << pca.getAvePosition()[2];
  nl() << "  - eigen values: " << std::setw(8) << std::right << pca.getEigenValues()[0] << ", "
       << pca.getEigenValues()[1] << ", " << pca.getEigenValues()[2];
  nl() << "  - average doca: " << pca.getAveHitDoca();
  nl() << "  - principle axis: " << std::setw(7) << std::setprecision(4)
       << pca.getEigenVectors()[0][0] << ", " << pca.getEigenVectors()[0][1] << ", "
       << pca.getEigenVectors()[0][2];
  nl() << "  - second axis: " << std::setw(7) << std::setprecision(4) << pca.getEigenVectors()[1][0]
       << ", " << pca.getEigenVectors()[1][1] << ", " << pca.getEigenVectors()[1][2];
  nl() << "  - third axis: " << std::setw(7) << std::setprecision(4) << pca.getEigenVectors()[2][0]
       << ", " << pca.getEigenVectors()[2][1] << ", " << pca.getEigenVectors()[2][2];

} // recob::dumper::DumpPCAxis()

//------------------------------------------------------------------------------

#endif // LARDATA_RECOBASE_DUMPERS_PCAXISDUMPERS_H
