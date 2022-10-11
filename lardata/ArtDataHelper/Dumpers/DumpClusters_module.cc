/**
 * @file   DumpClusters_module.cc
 * @brief  Dumps on screen the content of the clusters
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   August 13th, 2014
 */

// LArSoft includes
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C//C++ standard libraries
#include <algorithm> // std::min(), std::sort()
#include <iomanip>   // std::setw()
#include <sstream>
#include <string>

namespace recob {

  /**
   * @brief Prints the content of all the clusters on screen
   *
   * This analyser prints the content of all the clusters into the
   * LogInfo/LogVerbatim stream.
   *
   * Configuration parameters
   * -------------------------
   *
   * - *ClusterModuleLabel* (string, _required_): input tag from the
   *   producer used to create the recob::Cluster collection to be dumped
   * - *OutputCategory* (string, default: `"DumpClusters"`): the category
   *   used for the output (useful for filtering)
   * - *HitsPerLine* (integer, default: `20`): the dump of hits
   *   will put this many of them for each line
   *
   */
  class DumpClusters : public art::EDAnalyzer {
  public:
    /// Configuration object
    struct Config {
      using Comment = fhicl::Comment;
      using Name = fhicl::Name;

      fhicl::Atom<art::InputTag> ClusterModuleLabel{
        Name("ClusterModuleLabel"),
        Comment("input tag for the clusters to be dumped")};
      fhicl::Atom<std::string> OutputCategory{
        Name("OutputCategory"),
        Comment("name of the category used for message facility output"),
        "DumpClusters"};
      fhicl::Atom<unsigned int> HitsPerLine{
        Name("HitsPerLine"),
        Comment("number of hits per line (0 suppresses hit dumping)"),
        20U};

    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;

    /// Default constructor
    explicit DumpClusters(Parameters const& config);

    /// Does the printing
    void analyze(const art::Event& evt);

  private:
    art::InputTag fClusterModuleLabel; ///< tag of the cluster data product
    std::string fOutputCategory;       ///< category for LogInfo output
    unsigned int fHitsPerLine;         ///< hits per line in the output

  }; // class DumpWires

} // namespace recob

//------------------------------------------------------------------------------
namespace {

  /// Returns the length of the string representation of the specified object
  template <typename T>
  size_t StringLength(const T& value)
  {
    std::ostringstream sstr;
    sstr << value;
    return sstr.str().length();
  } // StringLength()

} // local namespace

namespace recob {

  //-------------------------------------------------
  DumpClusters::DumpClusters(Parameters const& config)
    : EDAnalyzer(config)
    , fClusterModuleLabel(config().ClusterModuleLabel())
    , fOutputCategory(config().OutputCategory())
    , fHitsPerLine(config().HitsPerLine())
  {}

  //-------------------------------------------------
  void DumpClusters::analyze(const art::Event& evt)
  {

    // fetch the data to be dumped on screen
    art::InputTag ClusterInputTag(fClusterModuleLabel);

    auto Clusters = evt.getValidHandle<std::vector<recob::Cluster>>(ClusterInputTag);

    // get cluster-hit associations
    art::FindManyP<recob::Hit> HitAssn(Clusters, evt, ClusterInputTag);

    mf::LogInfo(fOutputCategory) << "The event contains " << Clusters->size() << " '"
                                 << ClusterInputTag.encode() << "' clusters";

    unsigned int iCluster = 0;
    std::vector<size_t> HitBuffer(fHitsPerLine), LastBuffer;
    for (const recob::Cluster& cluster : *Clusters) {
      decltype(auto) ClusterHits = HitAssn.at(iCluster);

      // print a header for the cluster
      mf::LogVerbatim(fOutputCategory)
        << "Cluster #" << (iCluster++) << " from " << ClusterHits.size() << " hits: " << cluster;

      // print the hits of the cluster
      if ((fHitsPerLine > 0) && !ClusterHits.empty()) {
        std::vector<size_t> HitIndices;
        for (art::Ptr<recob::Hit> pHit : ClusterHits)
          HitIndices.push_back(pHit.key());
        std::sort(HitIndices.begin(), HitIndices.end());

        unsigned int Padding = ::StringLength(HitIndices.back());

        mf::LogVerbatim(fOutputCategory) << "  hit indices:";

        std::vector<size_t>::const_iterator iHit = HitIndices.begin(), hend = HitIndices.end();
        size_t RangeStart = *iHit, RangeStop = RangeStart;
        std::ostringstream output_line;
        size_t nItemsInLine = 0;
        while (++iHit != hend) {

          if (*iHit == RangeStop + 1) { ++RangeStop; }
          else {
            // the new item does not belong to the current range:
            // - print the current range
            if (RangeStart == RangeStop) {
              output_line << "  " << std::setw(Padding) << RangeStart;
              ++nItemsInLine;
            }
            else {
              char fill = (RangeStart + 1 == RangeStop) ? ' ' : '-';
              output_line << "  " << std::setw(Padding) << RangeStart << fill << fill
                          << std::setw(Padding) << std::setfill(fill) << RangeStop
                          << std::setfill(' ');
              nItemsInLine += 2;
            }
            // - start a new one
            RangeStart = RangeStop = *iHit;
          } // if ... else

          // if we have enough stuff in the buffer, let's print it
          if (nItemsInLine >= fHitsPerLine) {
            nItemsInLine = 0;
            mf::LogVerbatim(fOutputCategory) << " " << output_line.str();
            output_line.str("");
          }

        } // while

        mf::LogVerbatim line_out(fOutputCategory);
        line_out << " " << output_line.str();
        if (RangeStart == RangeStop)
          line_out << "  " << std::setw(Padding) << RangeStart;
        else {
          char fill = (RangeStart + 1 == RangeStop) ? ' ' : '-';
          line_out << "  " << std::setw(Padding) << RangeStart << fill << fill << std::setw(Padding)
                   << std::setfill(fill) << RangeStop;
        }
      } // if dumping the hits

    } // for clusters

  } // DumpClusters::analyze()

  DEFINE_ART_MODULE(DumpClusters)

} // namespace recob
