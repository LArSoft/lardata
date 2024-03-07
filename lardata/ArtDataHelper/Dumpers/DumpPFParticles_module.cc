/**
 * @file   DumpPFParticles_module.cc
 * @brief  Dumps on screen the content of ParticleFlow particles
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 25th, 2015
 */

// LArSoft includes
#include "lardata/ArtDataHelper/Dumpers/hexfloat.h"
#include "lardataobj/RecoBase/PFParticle.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Utilities/InputTag.h"

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// C//C++ standard libraries
#include <string>

// ... and more in the implementation part

namespace recob {

  /**
   * @brief Prints the content of all the ParticleFlow particles on screen
   *
   * This analyser prints the content of all the ParticleFlow particles into the
   * LogInfo/LogVerbatim stream.
   *
   * Configuration parameters
   * =========================
   *
   * - *PFModuleLabel* (art::InputTag, _required_): label of the
   *   producer used to create the recob::PFParticle collection to be dumped
   * - *OutputCategory* (string, default: `"DumpPFParticles"`): the category
   *   used for the output (useful for filtering)
   * - *PrintHexFloats* (boolean, default: `false`): print all the floating
   *   point numbers in base 16
   * - *MaxDepth* (unsigned int, optional): if specified, at most this number of
   *   particle generations will be printed; 1 means printing only primaries and
   *   their daughters, 0 only primaries. If not specified, no limit will be
   *   applied. This is useful for buggy PFParticles with circular references.
   * - *MakeParticleGraphs* (boolean, default: `false`): creates a DOT file for
   *   each event, with a graph of PFParticle relations; each file is named as:
   *   `ProcessName_ModuleLabel_InstanceName_Run#_Subrun#_Event#_particles.dot`,
   *   where the the input label elements refer to the data product being
   *   plotted.
   *
   *
   * Particle connection graphs
   * ---------------------------
   *
   * When _MakeParticleGraphs_ configuration option is activated, a file is
   * created for each event, that contains the particle flow tree in GraphViz
   * format. The GraphViz `dot` command can be used to render it into a PDF,
   * SVG, EPS or one of the many supported bitmap formats.
   * The typical command to use is:
   *
   *     dot -Tpdf -oPMTrk.pdf PMTrk.dot
   *
   * A `bash` command to convert all files into a `OutputFormat` format:
   *
   *     OutputFormat='pdf'
   *     for DotFile in *.dot ; do
   *       OutputFile="${DotFile%.dot}.${OutputFormat}"
   *       [[ "$OutputFile" -ot "$DotFile" ]] || continue # up to date already
   *       echo "${DotFile} => ${OutputFile} ..."
   *       dot -T"$OutputFormat" -o"$OutputFile" "$DotFile" || break
   *     done
   *
   * which will also skip files already converted.
   *
   * The output shows one cell ("node") per particle. The format of the node
   * follows these prescriptions:
   *
   * * the label has the particle ID number prepended by a hash character (`#`)
   * * if the particle has a PDG ID, that also appears in the label (either the
   *   name of the corresponding particle, or, if unknown, just the PDG ID
   *   number)
   * * if the particle is primary, it is rendered in bold font
   * * if the particle is referred by other particles, but it is not present
   *   ("ghost particle"), its border is red and dashed
   *
   * The relations between particles in the flow are represented by connecting
   * lines ("edges"). Connection information is redundant: the parent particle
   * should have the daughter in the daughter list, and the daughter should have
   * the parent particle referenced as such. Since the connection is usually
   * from two sources, there are usually two arrow heads, each one close to the
   * particle that provides information on that connection; all arrow heads
   * point from parent to daughter.
   *
   * * when the information of daughter and parent is consistent, a black line
   *   with two arrow heads both pointing to the daughter is shown
   * * when the parent is ghost, the arrow head close to the daughter is hollow;
   *   ghost particles have no arrow heads close to them
   * * when the daughter is ghost, the arrow head close to the parent is hollow;
   *   ghost particles have no arrow heads close to them
   *
   *
   * If you are trying to interpret an existing diagram, the following list is
   * more direct to the point.
   * Nodes: represent particles (see above for the label content)
   *
   *  * bold label: primary particle
   *  * red, dashed border: "ghost particle" (missing but referenced by others)
   *  * other: just a particle
   *
   * Connecting lines ("edges"):
   *  * all arrow heads point from parent to daughter
   *  * black with two full arrow heads: regular parent to daughter
   *  * black with a single inward empty arrow head: the particle close to the
   *    arrow claims the particle pointed by the arrow as a daughter, but there
   *    is no information on that daughter (ghost daughter)
   *  * black with a single outward empty arrow head: the particle at the tip of
   *    the arrow claims to be daughter of the other particle, but there is no
   *    information on that parent (ghost parent)
   *  * red, outward arrow: the daughter (at the tip of the only arrow) claims
   *    the other particle as parent, but that parent does not recognise it as
   *    daughter
   *  * orange, inward arrow: the parent (close to the only arrow head) claims
   *    the other particle as daughter, but that daighter does not recognise it
   *    as parent
   *
   */
  class DumpPFParticles : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> PFModuleLabel{
        Name("PFModuleLabel"),
        Comment("label of producer of the recob::PFParticle to be dumped")};

      fhicl::Atom<std::string> OutputCategory{
        Name("OutputCategory"),
        Comment("message facility category used for output (for filtering)"),
        "DumpPFParticles"};

      fhicl::Atom<bool> PrintHexFloats{Name("PrintHexFloats"),
                                       Comment("print all the floating point numbers in base 16"),
                                       false};

      fhicl::OptionalAtom<unsigned int> MaxDepth{
        Name("MaxDepth"),
        Comment("at most this number of particle generations will be printed")};

      fhicl::Atom<bool> MakeParticleGraphs{
        Name("MakeParticleGraphs"),
        Comment("creates a DOT file with particle information for each event"),
        false};

    }; // struct Config

    using Parameters = art::EDAnalyzer::Table<Config>;

    /// Default constructor
    explicit DumpPFParticles(Parameters const& config);

    /// Does the printing
    virtual void analyze(const art::Event& evt) override;

  private:
    art::InputTag fInputTag;     ///< input tag of the PFParticle product
    std::string fOutputCategory; ///< category for LogInfo output
    bool fPrintHexFloats;        ///< whether to print floats in base 16
    unsigned int fMaxDepth;      ///< maximum generation to print (0: only primaries)
    bool fMakeEventGraphs;       ///< whether to create one DOT file per event

    static std::string DotFileName(art::EventID const& evtID, art::Provenance const& prodInfo);

    void MakePFParticleGraph(art::Event const& event,
                             art::ValidHandle<std::vector<recob::PFParticle>> const& handle) const;

  }; // class DumpPFParticles

} // namespace recob

//==============================================================================
//===  Implementation section
//==============================================================================
// LArSoft includes
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

// art libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C//C++ standard libraries
#include <algorithm> // std::count(), std::find()
#include <fstream>
#include <limits>  // std::numeric_limits<>
#include <utility> // std::swap()

namespace {
  /**
   * @brief A container keyed by integer key (size_t).
   * @tparam T type of contained element
   */
  template <typename T>
  class int_map : private std::vector<T> {
    using map_type = std::vector<T>;

  public:
    using this_type = int_map<T>;

    using typename map_type::const_iterator;
    using typename map_type::const_reference;
    using typename map_type::iterator;
    using typename map_type::reference;
    using typename map_type::value_type;

    using typename map_type::size_type;

    /// Default constructor: "invalid data" is default-constructed
    int_map() = default;

    /// Constructor: sets the value of invalid data
    int_map(value_type invalid_value) : invalid(invalid_value) {}

    /// Resizes the map to accomodate this many items
    void resize(size_t new_size) { get_map().resize(new_size, invalid_value()); }

    /// Creates the item at specified position with invalid value and returns it
    reference operator[](size_t pos)
    {
      if (pos >= size()) resize(pos + 1);
      return get_map()[pos];
    }

    /// Returns the item at specified position, or invalid value if not existing
    const_reference operator[](size_t pos) const
    {
      return (pos >= size()) ? invalid : get_map()[pos];
    }

    using map_type::empty;
    using map_type::size;

    using map_type::begin;
    using map_type::cbegin;
    using map_type::cend;
    using map_type::crbegin;
    using map_type::crend;
    using map_type::end;
    using map_type::rbegin;
    using map_type::rend;

    /// Swaps the content of data and of this map
    void swap(this_type& data)
    {
      swap(data.get_map());
      std::swap(data.invalid, invalid);
    }

    /// Swaps the content of data and of this map
    void swap(std::vector<T>& data) { map_type::swap(data); }

    // non-standard calls follow
    /// Returns the number of invalid elements
    size_type n_invalid() const { return std::count(begin(), end(), invalid); }

    /// Returns the number of valid elements
    size_type n_valid() const { return size() - n_invalid(); }

    /// Returns whether the element at specified position is valid
    bool is_valid(size_type pos) const
    {
      return (pos < size()) ? is_valid_value(get_map()[pos]) : false;
    }

    /// Returns whether the specified value is valid
    bool is_valid_value(value_type const& v) const { return v != invalid; }

    /// Returns the invalid value
    value_type const& invalid_value() const { return invalid; }

  private:
    value_type const invalid; ///< value of invalid data

    map_type& get_map() { return static_cast<map_type&>(*this); }
    map_type const& get_map() const { return static_cast<map_type const&>(*this); }

  }; // int_map<>

  //----------------------------------------------------------------------------
  int_map<size_t> CreateMap(std::vector<recob::PFParticle> const& particles)
  {

    int_map<size_t> pmap(std::numeric_limits<size_t>::max());
    pmap.resize(particles.size());

    size_t const nParticles = particles.size();
    for (size_t iPart = 0; iPart < nParticles; ++iPart)
      pmap[particles[iPart].Self()] = iPart;

    return pmap;
  } // CreateMap()

  //----------------------------------------------------------------------------
  bool hasDaughter(recob::PFParticle const& particle, size_t partID)
  {
    auto const& daughters = particle.Daughters();
    return daughters.end() != std::find(daughters.begin(), daughters.end(), partID);
  } // hasDaughter()

  //----------------------------------------------------------------------------
  class ParticleDumper {
  public:
    /// Collection of available printing style options
    struct PrintOptions_t {
      bool hexFloats = false; ///< print all floating point numbers in base 16
      /// number of particle generations to descent into (0: only primaries)
      unsigned int maxDepth = std::numeric_limits<unsigned int>::max();
      /// name of the output stream
      std::string streamName;
    }; // PrintOptions_t

    /// Constructor; will dump particles from the specified list
    ParticleDumper(std::vector<recob::PFParticle> const& particle_list)
      : ParticleDumper(particle_list, {})
    {}

    /// Constructor; will dump particles from the specified list
    ParticleDumper(std::vector<recob::PFParticle> const& particle_list,
                   PrintOptions_t print_options)
      : particles(particle_list)
      , options(print_options)
      , visited(particles.size(), 0U)
      , particle_map(CreateMap(particles))
    {}

    /// Sets the vertices associated to each particle
    void SetVertices(art::FindOne<recob::Vertex> const* vtx_query) { vertices = vtx_query; }

    /// Sets the vertices associated to each particle
    void SetTracks(art::FindMany<recob::Track> const* trk_query) { tracks = trk_query; }

    /// Sets the clusters associated to each particle
    void SetClusters(art::FindMany<recob::Cluster> const* cls_query) { clusters = cls_query; }

    /// Sets the seeds associated to each particle
    void SetSeeds(art::FindMany<recob::Seed> const* seed_query) { seeds = seed_query; }

    /// Sets the 3D points associated to each particle
    void SetSpacePoints(art::FindMany<recob::SpacePoint> const* sp_query)
    {
      spacepoints = sp_query;
    }

    /// Sets the 3D points associated to each particle
    void SetPCAxes(art::FindMany<recob::PCAxis> const* pca_query) { pcaxes = pca_query; }

    /// Dump a particle specified by its index in the input particle list
    /// @param gen max generations to print
    template <typename Stream>
    void DumpParticle(Stream&& out,
                      size_t iPart,
                      std::string indentstr = "",
                      unsigned int gen = 0) const;

    /// Dump a particle specified by its ID
    /// @param gen max generations to print
    template <typename Stream>
    void DumpParticleWithID(Stream&& out,
                            size_t pID,
                            std::string indentstr = "",
                            unsigned int gen = 0) const;

    /// Dumps all primary particles
    void DumpAllPrimaries(std::string indentstr = "") const;

    /// Dumps all particles in the input list
    void DumpAllParticles(std::string indentstr = "") const;

    template <typename Stream>
    static void DumpPDGID(Stream&& out, int ID);

  protected:
    std::vector<recob::PFParticle> const& particles; ///< input list

    PrintOptions_t options; ///< printing and formatting options

    /// Associated vertices (expected same order as for particles)
    art::FindOne<recob::Vertex> const* vertices = nullptr;

    /// Associated tracks (expected same order as for particles)
    art::FindMany<recob::Track> const* tracks = nullptr;

    /// Associated clusters (expected same order as for particles)
    art::FindMany<recob::Cluster> const* clusters = nullptr;

    /// Associated seeds (expected same order as for particles)
    art::FindMany<recob::Seed> const* seeds = nullptr;

    /// Associated space points (expected same order as for particles)
    art::FindMany<recob::SpacePoint> const* spacepoints = nullptr;

    /// Associated principal component axes (expected same order as particles)
    art::FindMany<recob::PCAxis> const* pcaxes = nullptr;

    /// Number of dumps on each particle
    mutable std::vector<unsigned int> visited;

    int_map<size_t> const particle_map; ///< fast lookup index by particle ID

    template <typename Stream>
    void DumpPFParticleInfo(Stream&& out,
                            recob::PFParticle const& part,
                            unsigned int iPart,
                            std::string indentstr) const;

    template <typename Stream>
    void DumpVertex(Stream&& out,
                    cet::maybe_ref<recob::Vertex const> VertexRef,
                    std::string indentstr) const;

    template <typename Stream>
    void DumpPCAxisDirection(Stream&& out, recob::PCAxis const& axis) const;

    template <typename Stream>
    void DumpPCAxes(Stream&& out,
                    std::vector<recob::PCAxis const*> const& myAxes,
                    std::string indentstr) const;

    template <typename Stream>
    void DumpTrack(Stream&& out, recob::Track const& track) const;

    template <typename Stream>
    void DumpTracks(Stream&& out,
                    std::vector<recob::Track const*> const& myTracks,
                    std::string indentstr) const;

    template <typename Stream>
    void DumpSeed(Stream&& out, recob::Seed const& seed, std::string indentstr) const;

    template <typename Stream>
    void DumpSeeds(Stream&& out,
                   std::vector<recob::Seed const*> const& mySeeds,
                   std::string indentstr) const;

    template <typename Stream>
    void DumpSpacePoint(Stream&& out, recob::SpacePoint const& sp) const;

    template <typename Stream>
    void DumpSpacePoints(Stream&& out,
                         std::vector<recob::SpacePoint const*> const& mySpacePoints,
                         std::string indentstr) const;

    template <typename Stream>
    void DumpClusterSummary(Stream&& out, recob::Cluster const& track) const;

    template <typename Stream>
    void DumpClusters(Stream&& out,
                      std::vector<recob::Cluster const*> const& myClusters,
                      std::string indentstr) const;

  }; // class ParticleDumper

  //----------------------------------------------------------------------------
  class PFParticleGraphMaker {
  public:
    PFParticleGraphMaker() = default;

    template <typename Stream>
    void MakeGraph(Stream&& out, std::vector<recob::PFParticle> const& particles) const;

    template <typename Stream>
    void WriteGraphHeader(Stream&& out) const;

    template <typename Stream>
    void WriteParticleRelations(Stream&& out,
                                std::vector<recob::PFParticle> const& particles) const;

    template <typename Stream>
    void WriteParticleNodes(Stream&& out, std::vector<recob::PFParticle> const& particles) const;

    template <typename Stream>
    void WriteParticleEdges(Stream&& out, std::vector<recob::PFParticle> const& particles) const;

    template <typename Stream>
    void WriteGraphFooter(Stream&& out) const;

  }; // class PFParticleGraphMaker

  //----------------------------------------------------------------------------
  //--- ParticleDumper implementation
  //---
  template <typename Stream>
  void ParticleDumper::DumpParticle(Stream&& out,
                                    size_t iPart,
                                    std::string indentstr /* = "" */,
                                    unsigned int gen /* = 0 */
  ) const
  {
    lar::OptionalHexFloat hexfloat(options.hexFloats);

    recob::PFParticle const& part = particles.at(iPart);
    ++visited[iPart];

    if (visited[iPart] > 1) {
      out << indentstr << "particle " << part.Self() << " already printed!!!";
      return;
    }

    //
    // intro
    //
    DumpPFParticleInfo(std::forward<Stream>(out), part, iPart, indentstr);

    //
    // vertex information
    //
    if (vertices) DumpVertex(std::forward<Stream>(out), vertices->at(iPart), indentstr);

    // daughters tag
    if (part.NumDaughters() > 0)
      out << " with " << part.NumDaughters() << " daughters";
    else
      out << " with no daughter";

    //
    // axis
    //
    if (pcaxes) DumpPCAxes(std::forward<Stream>(out), pcaxes->at(iPart), indentstr);

    //
    // tracks
    //
    if (tracks) DumpTracks(std::forward<Stream>(out), tracks->at(iPart), indentstr);

    //
    // seeds
    //
    if (seeds) DumpSeeds(std::forward<Stream>(out), seeds->at(iPart), indentstr);

    //
    // space points
    //
    if (spacepoints) {
      DumpSpacePoints(std::forward<Stream>(out), spacepoints->at(iPart), indentstr);
    }

    //
    // clusters
    //
    if (clusters) { DumpClusters(std::forward<Stream>(out), clusters->at(iPart), indentstr); }

    //
    // daughters
    //
    auto const PartID = part.Self();
    if (part.NumDaughters() > 0) {
      std::vector<size_t> const& daughters = part.Daughters();
      out << "\n" << indentstr << "  " << part.NumDaughters() << " particle daughters";
      if (gen > 0) {
        out << ":";
        for (size_t DaughterID : daughters) {
          if (DaughterID == PartID) {
            out << "\n" << indentstr << "    oh, just great: this particle is its own daughter!";
          }
          else {
            out << '\n';
            DumpParticleWithID(out, DaughterID, indentstr + "  ", gen - 1);
          }
        }
      } // if descending
      else {
        out << " (further descend suppressed)";
      }
    } // if has daughters

    //
    // warnings
    //
    if (visited[iPart] == 2) {
      out << "\n"
          << indentstr << "  WARNING: particle ID=" << PartID << " connected more than once!";
    }

    //
    // done
    //

  } // ParticleDumper::DumpParticle()

  //----------------------------------------------------------------------------
  template <typename Stream>
  void ParticleDumper::DumpParticleWithID(Stream&& out,
                                          size_t pID,
                                          std::string indentstr /* = "" */,
                                          unsigned int gen /* = 0 */
  ) const
  {
    size_t const pos = particle_map[pID];
    if (particle_map.is_valid_value(pos)) { DumpParticle(out, pos, indentstr, gen); }
    else {
      out /* << "\n" */ << indentstr << "<ID=" << pID << " not found>";
    }
  } // ParticleDumper::DumpParticleWithID()

  //----------------------------------------------------------------------------
  void ParticleDumper::DumpAllPrimaries(std::string indentstr /* = "" */) const
  {
    indentstr += "  ";
    size_t const nParticles = particles.size();
    unsigned int nPrimaries = 0;
    for (size_t iPart = 0; iPart < nParticles; ++iPart) {
      if (!particles[iPart].IsPrimary()) continue;
      DumpParticle(mf::LogVerbatim(options.streamName), iPart, indentstr, options.maxDepth);
    } // for
    if (nPrimaries == 0) {
      mf::LogVerbatim(options.streamName) << indentstr << "No primary particle found";
    }
  } // ParticleDumper::DumpAllPrimaries()

  //----------------------------------------------------------------------------
  void ParticleDumper::DumpAllParticles(std::string indentstr /* = "" */) const
  {
    // first print all the primary particles
    DumpAllPrimaries(indentstr);
    // then find out if there are any that are "disconnected"
    unsigned int const nDisconnected = std::count(visited.begin(), visited.end(), 0U);
    if (nDisconnected) {
      mf::LogVerbatim(options.streamName)
        << indentstr << nDisconnected << " particles not coming from primary ones:";
      size_t const nParticles = visited.size();
      for (size_t iPart = 0; iPart < nParticles; ++iPart) {
        if (visited[iPart] > 0) continue;
        DumpParticle(
          mf::LogVerbatim(options.streamName), iPart, indentstr + "  ", options.maxDepth);
      } // for unvisited particles
      mf::LogVerbatim(options.streamName)
        << indentstr << "(end of " << nDisconnected << " particles not from primaries)";
    } // if there are disconnected particles
    // TODO finally, note if there are multiply-connected particles

  } // ParticleDumper::DumpAllParticles()

  //----------------------------------------------------------------------------
  template <typename Stream>
  void ParticleDumper::DumpPDGID(Stream&& out, int ID)
  {
    switch (ID) {
    case -11: out << "e+"; break;
    case 11: out << "e-"; break;
    case -13: out << "mu+"; break;
    case 13: out << "mu-"; break;
    default: out << "MCID=" << ID; break;
    } // switch
  }   // DumpPDGID()

  //----------------------------------------------------------------------------
  template <typename Stream>
  void ParticleDumper::DumpPFParticleInfo(Stream&& out,
                                          recob::PFParticle const& part,
                                          unsigned int iPart,
                                          std::string indentstr) const
  {
    auto const PartID = part.Self();
    out << indentstr << "ID=" << PartID;
    if (iPart != PartID) out << " [#" << iPart << "]";
    out << " (type: ";
    DumpPDGID(out, part.PdgCode());
    out << ")";
    if (part.IsPrimary())
      out << " (primary)";
    else
      out << " from ID=" << part.Parent();
  } // ParticleDumper::DumpPFParticleInfo()

  //----------------------------------------------------------------------------
  template <typename Stream>
  void ParticleDumper::DumpVertex(Stream&& out,
                                  cet::maybe_ref<recob::Vertex const> VertexRef,
                                  std::string indentstr) const
  {
    if (!VertexRef.isValid()) {
      out << " [no vertex found!]";
      return;
    }
    recob::Vertex const& vertex = VertexRef.ref();
    std::array<double, 3> vtx_pos;
    vertex.XYZ(vtx_pos.data());
    lar::OptionalHexFloat hexfloat(options.hexFloats);
    out << " [decay at (" << hexfloat(vtx_pos[0]) << "," << hexfloat(vtx_pos[1]) << ","
        << hexfloat(vtx_pos[2]) << "), ID=" << vertex.ID() << "]";

  } // ParticleDumper::DumpVertex()

  //----------------------------------------------------------------------------
  template <typename Stream>
  void ParticleDumper::DumpPCAxisDirection(Stream&& out, recob::PCAxis const& axis) const
  {
    if (!axis.getSvdOK()) {
      out << "axis is invalid";
      return;
    }
    lar::OptionalHexFloat hexfloat(options.hexFloats);
    out << "axis ID=" << axis.getID() << ", principal: (" << hexfloat(axis.getEigenVectors()[0][0])
        << ", " << hexfloat(axis.getEigenVectors()[0][1]) << ", "
        << hexfloat(axis.getEigenVectors()[0][2]) << ")";
  } // ParticleDumper::DumpPCAxisDirection()

  template <typename Stream>
  void ParticleDumper::DumpPCAxes(Stream&& out,
                                  std::vector<recob::PCAxis const*> const& myAxes,
                                  std::string indentstr) const
  {
    out << "\n" << indentstr;
    if (myAxes.empty()) {
      out << " [no axis found!]";
      return;
    }
    if (myAxes.size() == 1) { DumpPCAxisDirection(std::forward<Stream>(out), *(myAxes.front())); }
    else {
      out << "  " << myAxes.size() << " axes present:";
      for (recob::PCAxis const* axis : myAxes) {
        out << "\n" << indentstr << "    ";
        DumpPCAxisDirection(std::forward<Stream>(out), *axis);
      } // for axes
    }   // else
  }     // ParticleDumper::DumpPCAxes()

  //----------------------------------------------------------------------------
  template <typename Stream>
  void ParticleDumper::DumpSeed(Stream&& out, recob::Seed const& seed, std::string indentstr) const
  {
    if (!seed.IsValid()) {
      out << "  <invalid>";
      return;
    }
    std::array<double, 3> start, dir;
    seed.GetDirection(dir.data(), nullptr);
    seed.GetPoint(start.data(), nullptr);
    lar::OptionalHexFloat hexfloat(options.hexFloats);
    out << "\n"
        << indentstr << "    (" << hexfloat(start[0]) << "," << hexfloat(start[1]) << ","
        << hexfloat(start[2]) << ")=>(" << hexfloat(dir[0]) << "," << hexfloat(dir[1]) << ","
        << hexfloat(dir[2]) << "), " << hexfloat(seed.GetLength()) << " cm";
  } // ParticleDumper::DumpSeed()

  template <typename Stream>
  void ParticleDumper::DumpSeeds(Stream&& out,
                                 std::vector<recob::Seed const*> const& mySeeds,
                                 std::string indentstr) const
  {
    if (mySeeds.empty()) return;
    out << "\n" << indentstr << "  " << mySeeds.size() << " seeds:";
    for (recob::Seed const* seed : mySeeds)
      DumpSeed(std::forward<Stream>(out), *seed, indentstr);
  } // ParticleDumper::DumpSeeds()

  //----------------------------------------------------------------------------
  template <typename Stream>
  void ParticleDumper::DumpSpacePoint(Stream&& out, recob::SpacePoint const& sp) const
  {
    const double* pos = sp.XYZ();
    lar::OptionalHexFloat hexfloat(options.hexFloats);
    out << "  ID=" << sp.ID() << " at (" << hexfloat(pos[0]) << "," << hexfloat(pos[1]) << ","
        << hexfloat(pos[2]) << ") cm";
  } // ParticleDumper::DumpSpacePoints()

  template <typename Stream>
  void ParticleDumper::DumpSpacePoints(Stream&& out,
                                       std::vector<recob::SpacePoint const*> const& mySpacePoints,
                                       std::string indentstr) const
  {
    out << "\n" << indentstr << "  ";
    if (mySpacePoints.empty()) {
      out << "no space points";
      return;
    }
    constexpr unsigned int PointsPerLine = 2;
    out << mySpacePoints.size() << " space points:";
    for (size_t iSP = 0; iSP < mySpacePoints.size(); ++iSP) {
      if ((iSP % PointsPerLine) == 0) out << "\n" << indentstr << "  ";

      DumpSpacePoint(std::forward<Stream>(out), *(mySpacePoints[iSP]));
    } // for points
  }   // ParticleDumper::DumpSpacePoints()

  //----------------------------------------------------------------------------
  template <typename Stream>
  void ParticleDumper::DumpTrack(Stream&& out, recob::Track const& track) const
  {
    lar::OptionalHexFloat hexfloat(options.hexFloats);
    out << " length " << hexfloat(track.Length()) << "cm from (" << hexfloat(track.Vertex().X())
        << ";" << hexfloat(track.Vertex().Y()) << ";" << hexfloat(track.Vertex().Z()) << ") to ("
        << hexfloat(track.End().X()) << ";" << hexfloat(track.End().Y()) << ";"
        << hexfloat(track.End().Z()) << ") (ID=" << track.ID() << ")";
  } // ParticleDumper::DumpTrack()

  template <typename Stream>
  void ParticleDumper::DumpTracks(Stream&& out,
                                  std::vector<recob::Track const*> const& myTracks,
                                  std::string indentstr) const
  {
    if (myTracks.empty()) return;
    out << "\n" << indentstr << "  " << myTracks.size() << " tracks:";
    for (recob::Track const* track : myTracks) {
      if (myTracks.size() > 1) out << "\n" << indentstr << "   ";
      DumpTrack(std::forward<Stream>(out), *track);
    } // for tracks
  }   // ParticleDumper::DumpTracks()

  //----------------------------------------------------------------------------
  template <typename Stream>
  void ParticleDumper::DumpClusterSummary(Stream&& out, recob::Cluster const& cluster) const
  {
    out << "  " << cluster.NHits() << " hits on " << cluster.Plane() << " (ID=" << cluster.ID()
        << ")";
  } // ParticleDumper::DumpClusterSummary()

  template <typename Stream>
  void ParticleDumper::DumpClusters(Stream&& out,
                                    std::vector<recob::Cluster const*> const& myClusters,
                                    std::string indentstr) const
  {
    if (myClusters.empty()) return;
    out << "\n" << indentstr << "  " << myClusters.size() << " clusters:";
    for (recob::Cluster const* cluster : myClusters) {
      if (myClusters.size() > 1) out << "\n" << indentstr << "   ";
      DumpClusterSummary(std::forward<Stream>(out), *cluster);
    }
  } // ParticleDumper::DumpClusters()

  //----------------------------------------------------------------------------
  //---  PFParticleGraphMaker
  //---
  template <typename Stream>
  void PFParticleGraphMaker::MakeGraph(Stream&& out,
                                       std::vector<recob::PFParticle> const& particles) const
  {

    WriteGraphHeader(std::forward<Stream>(out));
    WriteParticleRelations(std::forward<Stream>(out), particles);
    WriteGraphFooter(std::forward<Stream>(out));

  } // PFParticleGraphMaker::MakeGraph()

  //----------------------------------------------------------------------------
  template <typename Stream>
  void PFParticleGraphMaker::WriteGraphHeader(Stream&& out) const
  {
    out << "\ndigraph {"
        << "\n";
  } // PFParticleGraphMaker::WriteGraphHeader()

  template <typename Stream>
  void PFParticleGraphMaker::WriteParticleNodes(
    Stream&& out,
    std::vector<recob::PFParticle> const& particles) const
  {

    for (auto const& particle : particles) {

      std::ostringstream label;
      label << "#" << particle.Self();
      if (particle.PdgCode() != 0) {
        label << ", ";
        ParticleDumper::DumpPDGID(label, particle.PdgCode());
      }

      out << "\n  P" << particle.Self() << " ["
          << " label = \"" << label.str() << "\"";
      if (particle.IsPrimary()) out << ", style = bold";
      out << " ]";
    } // for

  } // PFParticleGraphMaker::WriteParticleNodes()

  template <typename Stream>
  void PFParticleGraphMaker::WriteParticleEdges(
    Stream&& out,
    std::vector<recob::PFParticle> const& particles) const
  {
    auto particle_map = CreateMap(particles);

    out << "\n  "
        << "\n  // relations"
        << "\n  // "
        << "\n  // the arrow is close to the information provider,;"
        << "\n  // and it points from parent to daughter"
        << "\n  // "
        << "\n  // "
        << "\n  ";

    for (auto const& particle : particles) {
      auto const partID = particle.Self();

      // draw parent line
      if (!particle.IsPrimary()) {
        auto const parentID = particle.Parent();
        size_t const iParent = particle_map[parentID];
        if (!particle_map.is_valid_value(iParent)) {
          // parent is ghost
          out << "\nP" << parentID << " [ style = dashed, color = red"
              << ", label = \"(#" << parentID << ")\" ] // ghost particle"
              << "\nP" << parentID << " -> P" << partID
              << " [ dir = both, arrowhead = empty, arrowtail = none ]";
        }
        else {
          // parent is a known particle
          recob::PFParticle const& parent = particles[iParent];

          // is the relation bidirectional?
          if (hasDaughter(parent, partID)) {
            out << "\nP" << parentID << " -> P" << partID << " [ dir = both, arrowtail = inv ]";
          } // if bidirectional
          else {
            out << "\nP" << parentID << " -> P" << partID << " [ dir = forward, color = red ]"
                << " // claimed parent";
          } // if parent disowns

        } // if ghost ... else
      }   // if not primary

      // print daughter relationship only if daughters do not recognise us
      for (auto daughterID : particle.Daughters()) {
        size_t const iDaughter = particle_map[daughterID];
        if (!particle_map.is_valid_value(iDaughter)) {
          // daughter is ghost
          out << "\nP" << daughterID << " [ style = dashed, color = red"
              << ", label = \"(#" << daughterID << ")\" ] // ghost daughter"
              << "\nP" << partID << " -> P" << daughterID
              << " [ dir = both, arrowhead = none, arrowtail = invempty ]";
        }
        else {
          // daughter is a known particle
          recob::PFParticle const& daughter = particles[iDaughter];

          // is the relation bidirectional? (if so, the daughter will draw)
          if (daughter.Parent() != partID) {
            out << "\nP" << partID << " -> P" << daughterID
                << " [ dir = both, arrowhead = none, arrowtail = inv, color = orange ]"
                << " // claimed daughter";
          } // if parent disowns

        } // if ghost ... else

      } // for all daughters

    } // for all particles

  } // PFParticleGraphMaker::WriteParticleNodes()

  template <typename Stream>
  void PFParticleGraphMaker::WriteParticleRelations(
    Stream&& out,
    std::vector<recob::PFParticle> const& particles) const
  {
    out << "\n  // " << particles.size() << " particles (nodes)";
    WriteParticleNodes(std::forward<Stream>(out), particles);

    out << "\n"
        << "\n  // parent/children relations";
    WriteParticleEdges(std::forward<Stream>(out), particles);

  } // PFParticleGraphMaker::WriteParticleRelations()

  template <typename Stream>
  void PFParticleGraphMaker::WriteGraphFooter(Stream&& out) const
  {
    out << "\n"
        << "\n} // digraph" << std::endl;
  } // PFParticleGraphMaker::WriteGraphFooter()

  //----------------------------------------------------------------------------

} // local namespace

namespace recob {

  //----------------------------------------------------------------------------
  DumpPFParticles::DumpPFParticles(Parameters const& config)
    : EDAnalyzer(config)
    , fInputTag(config().PFModuleLabel())
    , fOutputCategory(config().OutputCategory())
    , fPrintHexFloats(config().PrintHexFloats())
    , fMaxDepth(std::numeric_limits<unsigned int>::max())
    , fMakeEventGraphs(config().MakeParticleGraphs())
  {
    // here we are handling the optional configuration key as it had just a
    // default value
    if (!config().MaxDepth(fMaxDepth)) fMaxDepth = std::numeric_limits<unsigned int>::max();
  }

  //----------------------------------------------------------------------------
  std::string DumpPFParticles::DotFileName(art::EventID const& evtID,
                                           art::Provenance const& prodInfo)
  {
    return prodInfo.processName() + '_' + prodInfo.moduleLabel() + '_' +
           prodInfo.productInstanceName() + "_Run" + std::to_string(evtID.run()) + "_Subrun" +
           std::to_string(evtID.subRun()) + "_Event" + std::to_string(evtID.event()) +
           "_particles.dot";
  } // DumpPFParticles::DotFileName()

  //----------------------------------------------------------------------------
  void DumpPFParticles::MakePFParticleGraph(
    art::Event const& event,
    art::ValidHandle<std::vector<recob::PFParticle>> const& handle) const
  {
    art::EventID const eventID = event.id();
    std::string fileName = DotFileName(eventID, *(handle.provenance()));
    std::ofstream outFile(fileName); // overwrite by default

    outFile << "// " << fileName << "\n// "
            << "\n// Created for run " << eventID.run() << " subrun " << eventID.subRun()
            << " event " << eventID.event() << "\n// "
            << "\n// dump of " << handle->size() << " particles"
            << "\n// " << std::endl;

    PFParticleGraphMaker graphMaker;
    graphMaker.MakeGraph(outFile, *handle);

  } // DumpPFParticles::MakePFParticleGraph()

  //----------------------------------------------------------------------------
  void DumpPFParticles::analyze(const art::Event& evt)
  {

    //
    // collect all the available information
    //
    // fetch the data to be dumped on screen
    art::ValidHandle<std::vector<recob::PFParticle>> PFParticles =
      evt.getValidHandle<std::vector<recob::PFParticle>>(fInputTag);

    if (fMakeEventGraphs) MakePFParticleGraph(evt, PFParticles);

    art::FindOne<recob::Vertex> const ParticleVertices(PFParticles, evt, fInputTag);
    art::FindMany<recob::Track> const ParticleTracks(PFParticles, evt, fInputTag);
    art::FindMany<recob::Cluster> const ParticleClusters(PFParticles, evt, fInputTag);
    art::FindMany<recob::Seed> const ParticleSeeds(PFParticles, evt, fInputTag);
    art::FindMany<recob::SpacePoint> const ParticleSpacePoints(PFParticles, evt, fInputTag);
    art::FindMany<recob::PCAxis> const ParticlePCAxes(PFParticles, evt, fInputTag);

    size_t const nParticles = PFParticles->size();
    mf::LogVerbatim(fOutputCategory) << "Event " << evt.id() << " contains " << nParticles
                                     << " particles from '" << fInputTag.encode() << "'";

    // prepare the dumper
    ParticleDumper::PrintOptions_t options;
    options.hexFloats = fPrintHexFloats;
    options.maxDepth = fMaxDepth;
    options.streamName = fOutputCategory;
    ParticleDumper dumper(*PFParticles, options);
    if (ParticleVertices.isValid())
      dumper.SetVertices(&ParticleVertices);
    else
      mf::LogPrint("DumpPFParticles") << "WARNING: vertex information not available";
    if (ParticleTracks.isValid())
      dumper.SetTracks(&ParticleTracks);
    else
      mf::LogPrint("DumpPFParticles") << "WARNING: track information not available";
    if (ParticleClusters.isValid())
      dumper.SetClusters(&ParticleClusters);
    else
      mf::LogPrint("DumpPFParticles") << "WARNING: cluster information not available";
    if (ParticleSeeds.isValid())
      dumper.SetSeeds(&ParticleSeeds);
    else
      mf::LogPrint("DumpPFParticles") << "WARNING: seed information not avaialble";
    if (ParticleSpacePoints.isValid())
      dumper.SetSpacePoints(&ParticleSpacePoints);
    else {
      mf::LogPrint("DumpPFParticles") << "WARNING: space point information not available";
    }
    if (ParticlePCAxes.isValid())
      dumper.SetPCAxes(&ParticlePCAxes);
    else {
      mf::LogPrint("DumpPFParticles") << "WARNING: principal component axis not available";
    }
    dumper.DumpAllParticles("  ");

    mf::LogVerbatim(fOutputCategory) << "\n"; // two empty lines

  } // DumpPFParticles::analyze()

  DEFINE_ART_MODULE(DumpPFParticles)

} // namespace recob
