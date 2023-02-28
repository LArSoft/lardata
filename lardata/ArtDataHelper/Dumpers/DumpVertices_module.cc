/**
 * @file   DumpVertices_module.cc
 * @brief  Dumps on screen the content of vertices
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 29th, 2015
 */

// LArSoft includes
#include "lardataobj/RecoBase/Vertex.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

// support libraries
#include "fhiclcpp/ParameterSet.h"

// C//C++ standard libraries
#include <string>

// ... and more in the implementation part

namespace recob {

  /**
   * @brief Prints the content of all the vertices on screen
   *
   * This analyser prints the content of all the vertices into the
   * LogInfo/LogVerbatim stream.
   *
   * Configuration parameters
   * =========================
   *
   * - *VertexModuleLabel* (art::InputTag, mandatory): label of the
   *   producer used to create the recob::Vertex collection to be dumped
   * - *OutputCategory* (string, default: `"DumpVertices"`): the category used
   *   for the output (useful for filtering)
   * - *PrintHexFloats* (boolean, default: `false`): print all the floating
   *   point numbers in base 16
   *
   */
  class DumpVertices : public art::EDAnalyzer {
  public:
    /// Default constructor
    explicit DumpVertices(fhicl::ParameterSet const& pset);

    /// Does the printing
    virtual void analyze(const art::Event& evt) override;

  private:
    art::InputTag fInputTag;     ///< input tag of the Vertex product
    std::string fOutputCategory; ///< category for LogInfo output
    bool fPrintHexFloats;        ///< whether to print floats in base 16

  }; // class DumpVertices

} // namespace recob

//==============================================================================
//===  Implementation section
//==============================================================================

// LArSoft includes
#include "lardata/ArtDataHelper/Dumpers/hexfloat.h"
#include "lardataobj/RecoBase/Vertex.h"

// art libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C//C++ standard libraries

namespace {

  //----------------------------------------------------------------------------
  class VertexDumper {
  public:
    /// Collection of available printing style options
    struct PrintOptions_t {
      bool hexFloats = false; ///< print all floating point numbers in base 16
    };                        // PrintOptions_t

    /// Constructor; will dump vertices from the specified list
    VertexDumper(std::vector<recob::Vertex> const& vertex_list) : VertexDumper(vertex_list, {}) {}

    /// Constructor; will dump vertices from the specified list, using options.
    VertexDumper(std::vector<recob::Vertex> const& vertex_list, PrintOptions_t print_options)
      : vertices(vertex_list), options(print_options)
    {}

    /// Dump a vertex specified by its index in the input list
    template <typename Stream>
    void DumpVertex(Stream&& out, size_t iVertex, std::string indentstr = "") const
    {
      lar::OptionalHexFloat hexfloat(options.hexFloats);

      recob::Vertex const& vertex = vertices.at(iVertex);

      //
      // intro
      //
      out << "\n" << indentstr << "[#" << iVertex << "]";

      std::array<double, 3> vtx_pos;
      vertex.XYZ(vtx_pos.data());
      out << " ID=" << vertex.ID() << " at (" << hexfloat(vtx_pos[0]) << "," << hexfloat(vtx_pos[1])
          << "," << hexfloat(vtx_pos[2]) << ")";

      //
      // done
      //

    } // DumpVertex()

    /// Dumps all vertices in the input list
    template <typename Stream>
    void DumpAllVertices(Stream&& out, std::string indentstr = "") const
    {
      indentstr += "  ";
      size_t const nVertices = vertices.size();
      for (size_t iVertex = 0; iVertex < nVertices; ++iVertex)
        DumpVertex(out, iVertex, indentstr);
    } // DumpAllVertices()

  protected:
    std::vector<recob::Vertex> const& vertices; ///< input list

    PrintOptions_t options; ///< printing and formatting options

  }; // VertexDumper

  //----------------------------------------------------------------------------

} // local namespace

namespace recob {

  //----------------------------------------------------------------------------
  DumpVertices::DumpVertices(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fInputTag(pset.get<art::InputTag>("VertexModuleLabel"))
    , fOutputCategory(pset.get<std::string>("OutputCategory", "DumpVertices"))
    , fPrintHexFloats(pset.get<bool>("PrintHexFloats", false))
  {}

  //----------------------------------------------------------------------------
  void DumpVertices::analyze(const art::Event& evt)
  {

    //
    // collect all the available information
    //
    // fetch the data to be dumped on screen
    auto Vertices = evt.getValidHandle<std::vector<recob::Vertex>>(fInputTag);

    size_t const nVertices = Vertices->size();
    mf::LogVerbatim(fOutputCategory) << "Event " << evt.id() << " contains " << nVertices
                                     << " vertices from '" << fInputTag.encode() << "'";

    // prepare the dumper
    VertexDumper::PrintOptions_t options;
    options.hexFloats = fPrintHexFloats;
    VertexDumper dumper(*Vertices, options);

    dumper.DumpAllVertices(mf::LogVerbatim(fOutputCategory), "  ");

    mf::LogVerbatim(fOutputCategory) << "\n"; // two empty lines

  } // DumpVertices::analyze()

  DEFINE_ART_MODULE(DumpVertices)

} // namespace recob
