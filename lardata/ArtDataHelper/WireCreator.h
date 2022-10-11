/** ****************************************************************************
 * @file   WireCreator.h
 * @brief  Helper functions to create a wire
 * @date   December 11, 2014
 * @author petrillo@fnal.gov
 * @see    Wire.h WireCreator.cxx
 *
 * ****************************************************************************/

#ifndef WIRECREATOR_H
#define WIRECREATOR_H

// C/C++ standard library
#include <utility> // std::move()

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"  // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::View_t
#include "lardataobj/RecoBase/Wire.h"

namespace raw {
  class RawDigit;
}

/// Reconstruction base classes
namespace recob {

  /**
   * @brief Class managing the creation of a new recob::Wire object
   *
   * In order to be as simple as possible (Plain Old Data), data products like
   * recob::Wire need to be stripped of most of their functions, including the
   * ability to communicate whether a value we try to store is invalid
   * (that would require a art::Exception -- art -- or at least a message on the
   * screen -- MessageFacility) and the ability to read things from event,
   * services (e.g. geometry) etc.
   *
   * A Creator is a class that creates a temporary data product, and at the
   * end it yields it to the caller for storage.
   * This last step should be by move construction, although a copy method is
   * also provided.
   *
   * An example of creating a Wire object:
   *
   *     // let RoIsignal be a recob::Wire::RegionsOfInterest_t already filled
   *     // with the signal regions, and rawdigit the raw::RawDigit of the
   *     // channel; RoIsignal will become empty
   *     recob::WireCreator wire(std::move(RoIsignal), rawdigit);
   *     wires.push_back(wire.move()); // wire content is not valid any more
   *
   * This is a one-step creation object: the wire is constructed at the same
   * time the WireCreator is, and no facility is offered to modify the
   * constructed wire, or to create another one.
   */
  class WireCreator {
  public:
    /// Alias for the type of regions of interest
    using RegionsOfInterest_t = Wire::RegionsOfInterest_t;

    // destructor, copy and move constructor and assignment as default

    /**
       * @brief Constructor: uses specified signal in regions of interest
       * @param sigROIlist signal organized in regions of interest
       * @param rawdigit the raw digit this channel is associated to
       *
       * The information used from the raw digit are the channel ID and the
       * length in samples (TDC ticks) of the original readout window.
       */
    WireCreator(const RegionsOfInterest_t& sigROIlist, const raw::RawDigit& rawdigit);

    /**
       * @brief Constructor: uses specified signal in regions of interest
       * @param sigROIlist signal organized in regions of interest
       * @param rawdigit the raw digit this channel is associated to
       *
       * The information used from the raw digit are the channel ID and the
       * length in samples (TDC ticks) of the original readout window.
       *
       * Signal information is moved from sigROIlist, that becomes empty.
       */
    WireCreator(RegionsOfInterest_t&& sigROIlist, const raw::RawDigit& rawdigit);

    /**
       * @brief Constructor: uses specified signal in regions of interest
       * @param sigROIlist signal organized in regions of interest
       * @param channel the ID of the channel
       * @param view the view the channel belongs to
       *
       * The information used from the raw digit are the channel ID and the
       * length in samples (TDC ticks) of the original readout window.
       */
    WireCreator(RegionsOfInterest_t const& sigROIlist, raw::ChannelID_t channel, geo::View_t view);

    /**
       * @brief Constructor: uses specified signal in regions of interest
       * @param sigROIlist signal organized in regions of interest
       * @param channel the ID of the channel
       * @param view the view the channel belongs to
       *
       * The information used from the raw digit are the channel ID and the
       * length in samples (TDC ticks) of the original readout window.
       *
       * Signal information is moved from sigROIlist, that becomes empty.
       */
    WireCreator(RegionsOfInterest_t&& sigROIlist, raw::ChannelID_t channel, geo::View_t view);

    /**
       * @brief Prepares the constructed wire to be moved away
       * @return a right-value reference to the constructed wire
       *
       * Despite the name, no move happens in this function.
       * Move takes place in the caller code as proper; for example:
       *
       *     // be wire a WireCreator instance:
       *     std::vector<recob::Wire> Wires;
       *     wire.move();                          // nothing happens
       *     Wires.push_back(wire.move());         // here the copy happens
       *     recob::Wire single_wire(wire.move()); // wrong! wire is empty now
       *
       */
    Wire&& move() { return std::move(wire); }

    /**
       * @brief Returns the constructed wire
       * @return a constant reference to the constructed wire
       *
       * Despite the name, no copy happens in this function.
       * Copy takes place in the caller code as proper; for example:
       *
       *     // be wire a WireCreator instance:
       *     std::vector<recob::Wire> Wires;
       *     wire.copy();                          // nothing happens
       *     Wires.push_back(wire.copy());         // here a copy happens
       *     recob::Wire single_wire(wire.copy()); // wire is copied again
       *
       */
    const Wire& copy() const { return wire; }

  protected:
    Wire wire; ///< local instance of the wire being constructed

  }; // class WireCreator

} // namespace recob

#endif // WIRECREATOR_H
