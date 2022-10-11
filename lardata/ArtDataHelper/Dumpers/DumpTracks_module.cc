/**
 * @file   DumpTracks_module.cc
 * @brief  Dumps on screen the content of the tracks
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 12th, 2014
 */

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"

// C//C++ standard libraries
#include <algorithm>  // std::max(), std::sort(), std::transform()
#include <functional> // std::mem_fn()
#include <iomanip>    // std::setw()
#include <iterator>   // std::back_inserter()
#include <memory>     // std::unique_ptr()
#include <sstream>
#include <string>

namespace {

  /// Returns the length of the string representation of the specified object
  template <typename T>
  size_t StringLength(const T& value)
  {
    std::ostringstream sstr;
    sstr << value;
    return sstr.str().length();
  } // StringLength()

  /**
   * @brief Prints a table with compact integers
   * @param log stream to send the output to
   * @param Indices sorted container of the indices to be printed
   * @param IndicesPerLine number of indices printed in one line
   * @param IndentStr string to be printed at the beginning of each new line
   *
   */
  template <typename STREAM, typename CONT>
  void PrintCompactIndexTable(STREAM& log,
                              const CONT& Indices,
                              unsigned int IndicesPerLine,
                              std::string IndentStr = "")
  {
    unsigned int Padding = StringLength(Indices.back());

    typename CONT::const_iterator iIndex = Indices.begin(), iend = Indices.end();
    size_t RangeStart = *iIndex, RangeStop = RangeStart;
    std::ostringstream output_line;
    size_t nItemsInLine = 0;
    while (++iIndex != iend) {

      if (*iIndex == RangeStop + 1) { ++RangeStop; }
      else {
        // the new item does not belong to the current range:
        // - print the current range
        if (nItemsInLine) output_line << "  ";
        if (RangeStart == RangeStop) {
          output_line << std::setw(Padding) << RangeStart;
          ++nItemsInLine;
        }
        else {
          char fill = (RangeStart + 1 == RangeStop) ? ' ' : '-';
          output_line << std::setw(Padding) << RangeStart << fill << fill << std::setw(Padding)
                      << std::setfill(fill) << RangeStop << std::setfill(' ');
          nItemsInLine += 2;
        }
        // - start a new one
        RangeStart = RangeStop = *iIndex;
      } // if ... else

      // if we have enough stuff in the buffer, let's print it
      if (nItemsInLine >= IndicesPerLine) {
        nItemsInLine = 0;
        log << IndentStr << output_line.str() << "\n";
        output_line.str("");
      }

    } // while

    // print what we have accumulated so far of the last line
    log << IndentStr << output_line.str();
    // print the last range (or single element)
    if (nItemsInLine) log << "  ";
    if (RangeStart == RangeStop)
      log << std::setw(Padding) << RangeStart;
    else {
      char fill = (RangeStart + 1 == RangeStop) ? ' ' : '-';
      log << std::setw(Padding) << RangeStart << fill << fill << std::setw(Padding)
          << std::setfill(fill) << RangeStop << std::setfill(' ');
    }
    log << std::endl;
  } // PrintCompactIndexTable()

  /**
   * @brief Prints a table with keys from a container of objects
   * @param log stream to send the output to
   * @param Indices sorted container of the indices to be printed
   * @param IndicesPerLine number of indices printed in one line
   * @param IndexExtractor functor extracting the index from a element
   * @param IndentStr string to be printed at the beginning of each new line
   *
   * The key extractor must oterate on elements of Indices and returning
   * something that can be converted to a size_t.
   */
  template <typename STREAM, typename CONT, typename GETINDEX>
  void PrintCompactIndexTable(STREAM& log,
                              const CONT& Objects,
                              unsigned int IndicesPerLine,
                              GETINDEX IndexExtractor,
                              std::string IndentStr)
  {
    if ((IndicesPerLine == 0) || Objects.empty()) return;

    std::vector<size_t> Indices;
    Indices.reserve(Objects.size());
    std::transform(Objects.begin(), Objects.end(), std::back_inserter(Indices), IndexExtractor);
    std::sort(Indices.begin(), Indices.end());
    PrintCompactIndexTable(log, Indices, IndicesPerLine, IndentStr);

  } // PrintCompactIndexTable()

  /**
   * @brief Prints a table with indices from a container of objects
   * @param log stream to send the output to
   * @param Objects container of art::Ptr to objects to be printed
   * @param IndicesPerLine number of indices printed in one line
   * @param IndentStr string to be printed at the beginning of each new line
   *
   * The class in the container must have a key() member returning something
   * that can be converted to a size_t.
   * This is designed for containers like std::vector<art::Ptr<T>>.
   */
  template <typename STREAM, typename T>
  inline void PrintAssociatedIndexTable(STREAM& log,
                                        const std::vector<art::Ptr<T>>& Objects,
                                        unsigned int IndicesPerLine,
                                        std::string IndentStr = "")
  {
    PrintCompactIndexTable(log, Objects, IndicesPerLine, std::mem_fn(&art::Ptr<T>::key), IndentStr);
  } // PrintAssociatedIndexTable()

  /**
   * @brief Prints a table with indices from a container of objects
   * @param log stream to send the output to
   * @param Objects container of art::Ptr to objects to be printed
   * @param IndicesPerLine number of indices printed in one line
   * @param IndentStr string to be printed at the beginning of each new line
   *
   * The class in the container must have a key() member returning something
   * that can be converted to a size_t.
   * This is designed for containers like std::vector<art::Ptr<T>>.
   */
  template <typename STREAM, typename T>
  inline void PrintAssociatedIDTable(STREAM& log,
                                     const std::vector<art::Ptr<T>>& Objects,
                                     unsigned int IndicesPerLine,
                                     std::string IndentStr = "")
  {
    PrintCompactIndexTable(
      log, Objects, IndicesPerLine, [](const art::Ptr<T>& ptr) { return ptr->ID(); }, IndentStr);
  } // PrintAssociatedIDTable()

} // local namespace

namespace recob {

  /**
   * @brief Prints the content of all the tracks on screen
   *
   * This analyser prints the content of all the tracks into the
   * LogInfo/LogVerbatim stream.
   * The associated objects are printed only if they were produced with the
   * same input tag as the tracks.
   *
   * Configuration parameters
   * -------------------------
   *
   * - *TrackModuleLabel* (string, _required_): label of the
   *   producer used to create the recob::Track collection to be dumped
   * - *OutputCategory* (string, default: `"DumpTracks"`): the category
   *   used for the output (useful for filtering)
   * - *WayPoints* (unsigned integer, default: `10`): approximate number
   *   of way points printed in the output
   * - *SpacePointAssociations* (boolean, default: `true`): prints the number
   *   of space points associated with the tracks
   * - *PrintSpacePoints* (boolean, default: `false`): print the index of all hits
   *   associated with the tracks
   * - *HitAssociations* (boolean, default: `true`): prints the number of
   *   hits associated with the tracks
   * - *PrintHits* (boolean, default: `false`): print the index of all hits
   *   associated with the tracks
   * - *ParticleAssociations* (boolean, default: `true`): prints the number
   *   of particle-flow particles associated with the tracks
   *
   */
  class DumpTracks : public art::EDAnalyzer {
  public:
    /// Configuration object
    struct Config {
      using Comment = fhicl::Comment;
      using Name = fhicl::Name;

      fhicl::Atom<art::InputTag> TrackModuleLabel{Name("TrackModuleLabel"),
                                                  Comment("input tag for the tracks to be dumped")};
      fhicl::Atom<std::string> OutputCategory{
        Name("OutputCategory"),
        Comment("name of the category used for message facility output"),
        "DumpTracks"};
      fhicl::Atom<unsigned int> WayPoints{Name("WayPoints"),
                                          Comment("number of points along the trajectory printed"),
                                          10U};
      fhicl::Atom<bool> SpacePointAssociations{
        Name("SpacePointAssociations"),
        Comment("prints the number of space points associated to the track"),
        true};
      fhicl::Atom<bool> PrintSpacePoints{
        Name("PrintSpacePoints"),
        Comment("prints the index of all space points associated to the track"),
        false};
      fhicl::Atom<bool> HitAssociations{
        Name("HitAssociations"),
        Comment("prints the number of hits associated to the track"),
        true};
      fhicl::Atom<bool> PrintHits{Name("PrintHits"),
                                  Comment("prints the index of all hits associated to the track"),
                                  false};
      fhicl::Atom<bool> ParticleAssociations{
        Name("ParticleAssociations"),
        Comment("prints the number of PF particles associated to the track"),
        true};

    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;

    /// Default constructor
    explicit DumpTracks(Parameters const& config);

    /// Does the printing
    void analyze(const art::Event& evt);

  private:
    art::InputTag fTrackModuleLabel; ///< name of module that produced the tracks
    std::string fOutputCategory;     ///< category for LogInfo output
    unsigned int fPrintWayPoints;    ///< number of printed way points

    bool fPrintNHits;        ///< prints the number of associated hits
    bool fPrintNSpacePoints; ///< prints the number of associated space points
    bool fPrintNParticles;   ///< prints the number of associated PFParticles
    bool fPrintHits;         ///< prints the index of associated hits
    bool fPrintSpacePoints;  ///< prints the index of associated space points
    bool fPrintParticles;    ///< prints the index of associated PFParticles

    /// Dumps information about the specified track
    void DumpTrack(unsigned int iTrack, recob::Track const& track) const;

  }; // class DumpTracks

} // namespace recob

//------------------------------------------------------------------------------

namespace recob {

  //-------------------------------------------------
  DumpTracks::DumpTracks(Parameters const& config)
    : EDAnalyzer(config)
    , fTrackModuleLabel(config().TrackModuleLabel())
    , fOutputCategory(config().OutputCategory())
    , fPrintWayPoints(config().WayPoints())
    , fPrintNHits(config().HitAssociations())
    , fPrintNSpacePoints(config().SpacePointAssociations())
    , fPrintNParticles(config().ParticleAssociations())
    , fPrintHits(config().PrintHits())
    , fPrintSpacePoints(config().PrintSpacePoints())
    , fPrintParticles(config().ParticleAssociations())
  {}

  //-------------------------------------------------
  void DumpTracks::analyze(const art::Event& evt)
  {

    // fetch the data to be dumped on screen
    auto Tracks = evt.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);

    mf::LogInfo(fOutputCategory) << "The event contains " << Tracks->size() << " '"
                                 << fTrackModuleLabel.encode() << "'tracks";

    std::unique_ptr<art::FindManyP<recob::Hit>> pHits(
      fPrintNHits ? new art::FindManyP<recob::Hit>(Tracks, evt, fTrackModuleLabel) : nullptr);
    if (pHits && !pHits->isValid()) {
      throw art::Exception(art::errors::ProductNotFound)
        << "No hit associated with '" << fTrackModuleLabel.encode() << "' tracks.\n";
    }

    std::unique_ptr<art::FindManyP<recob::SpacePoint>> pSpacePoints(
      fPrintNSpacePoints ? new art::FindManyP<recob::SpacePoint>(Tracks, evt, fTrackModuleLabel) :
                           nullptr);
    if (pSpacePoints && !pSpacePoints->isValid()) {
      throw art::Exception(art::errors::ProductNotFound)
        << "No space point associated with '" << fTrackModuleLabel.encode() << "' tracks.\n";
    }

    std::unique_ptr<art::FindManyP<recob::PFParticle>> pPFParticles(
      fPrintNParticles ? new art::FindManyP<recob::PFParticle>(Tracks, evt, fTrackModuleLabel) :
                         nullptr);
    if (pPFParticles && !pPFParticles->isValid()) {
      throw art::Exception(art::errors::ProductNotFound)
        << "No particle-flow particle associated with '" << fTrackModuleLabel.encode()
        << "' tracks.\n";
    }

    for (unsigned int iTrack = 0; iTrack < Tracks->size(); ++iTrack) {
      const recob::Track& track = Tracks->at(iTrack);

      // print track information
      DumpTrack(iTrack, track);

      mf::LogVerbatim log(fOutputCategory);
      if (pHits || pSpacePoints || pPFParticles) {
        log << "\n  associated with:";
        if (pHits) log << " " << pHits->at(iTrack).size() << " hits;";
        if (pSpacePoints) log << " " << pSpacePoints->at(iTrack).size() << " space points;";
        if (pPFParticles) log << " " << pPFParticles->at(iTrack).size() << " PF particles;";
      } // if we have any association

      if (pHits && fPrintHits) {
        const auto& Hits = pHits->at(iTrack);
        log << "\n  hit indices (" << Hits.size() << "):\n";
        PrintAssociatedIndexTable(log, Hits, 10 /* 10 hits per line */, "    ");
      } // if print individual hits

      if (pSpacePoints && fPrintSpacePoints) {
        const auto& SpacePoints = pSpacePoints->at(iTrack);
        log << "\n  space point IDs (" << SpacePoints.size() << "):\n";
        PrintAssociatedIDTable(log, SpacePoints, 10 /* 10 hits per line */, "    ");
      } // if print individual space points

      if (pPFParticles && fPrintParticles) {
        const auto& PFParticles = pPFParticles->at(iTrack);
        log << "\n  particle indices (" << PFParticles.size() << "):\n";
        // currently a particle has no ID
        PrintAssociatedIndexTable(log, PFParticles, 10 /* 10 hits per line */, "    ");
      } // if print individual particles
    }   // for tracks
  }     // DumpTracks::analyze()

  //---------------------------------------------------------------------------
  void DumpTracks::DumpTrack(unsigned int iTrack, recob::Track const& track) const
  {
    // print a header for the track
    const unsigned int nPoints = track.NumberTrajectoryPoints();
    mf::LogVerbatim log(fOutputCategory);
    log << "Track #" << iTrack << " ID: " << track.ID() << std::fixed << std::setprecision(3)
        << " theta: " << track.Theta() << " rad, phi: " << track.Phi()
        << " rad, length: " << track.Length() << " cm"
        << "\n  start at: ( " << track.Vertex().X() << " ; " << track.Vertex().Y() << " ; "
        << track.Vertex().Z() << " ), direction: ( " << track.VertexDirection().X() << " ; "
        << track.VertexDirection().Y() << " ; " << track.VertexDirection().Z() << " )"
        << "\n  end at:   ( " << track.End().X() << " ; " << track.End().Y() << " ; "
        << track.End().Z() << " ), direction: ( " << track.EndDirection().X() << " ; "
        << track.EndDirection().Y() << " ; " << track.EndDirection().Z() << " )"
        << "\n  with " << nPoints << " trajectory points";

    if (fPrintWayPoints > 0) {
      // print up to 10 (actually, 8 or 9) way points
      log << "\n  passes through:";
      unsigned int skip = std::max(nPoints / fPrintWayPoints, 1U);
      unsigned int iPoint = 0;
      while ((iPoint += skip) < nPoints) {
        const auto& point = track.LocationAtPoint(iPoint);
        log << "\n    [#" << iPoint << "] (" << point.X() << ", " << point.Y() << ", " << point.Z()
            << ")";
      } // while (iPoint)
    }   // if print way points
  }     // DumpTracks::DumpTrack()

  //---------------------------------------------------------------------------
  DEFINE_ART_MODULE(DumpTracks)

} // namespace recob
