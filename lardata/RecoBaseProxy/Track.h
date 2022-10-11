/**
 * @file   lardata/RecoBaseProxy/Track.h
 * @brief  Offers `proxy::Tracks` and `proxy::Track` class for
 *         `recob::Track` access.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 *
 */

/**
 * @defgroup LArSoftProxyTracks proxy::Tracks (recob::Track proxy)
 * @brief Proxy for a `recob::Track` collection.
 * @ingroup LArSoftProxyReco
 *
 *
 * Track proxies is a way to facilitate the navigation of `recob::Track` data
 * objects.
 * The ensemble of fundamental data of a track collection includes:
 * * the tracks themselves, in a `std::vector<recob::Track>` collection
 * * the associated hits, in a `art::Assns<recob::Track, recob::Hit>` data
 *   product
 *
 * Special customisations are provided for:
 * * the associated hits (automatically pulled in): the
 *   information is provided with `recob::Hit` tag, with the
 *   dedicated accessor `hits()`, `nHits()` and `hitAtPoint()` of the track
 *    proxy and with `hitPtr()` and `hit()` when accessing a single point
 * * the track fit information: include it with `withFitHitInfo()`; the
 *   information is provided with `recob::TrackFitHitInfo` tag, with the
 *   dedicated accessor `fitInfoAtPoint()` of the track proxy and with
 *   `fitInfoPtr()` when accessing a single point
 *
 * LArSoft prescribes conventions to be followed, which include:
 * * a track has at least two trajectory points
 * * for each track, there is one hit per trajectory point
 * * the association between tracks and hits is created with tracks as first
 *   ("left") element, and hits as second one ("right")
 * * hits in the association are in a well-defined order: first are the hits of
 *   the first track, sorted in the same way as their trajectory points; then
 *   the second track hits come, likewise; and all tracks follow in order;
 *   this is called the "one-to-many sequential association" requirement
 *   documented in `ProxyBase.h`.
 *
 * For track data products respecting this convention, a track proxy provides
 * an interface to navigate the track information.
 *
 * @note The interface is experimental, and it is likely not to include all the
 *       features you may need. If you find a missing feature, or find a use
 *       case violating an assumption, or you find the interface cumbersome to
 *       use, please contact the author for a discussion on how to improve this
 *       utility.
 *
 *
 * Obtaining a track proxy
 * ========================
 *
 * Track proxies are created by specifying the tag of the tracks, and the event
 * to read them from.
 *
 * To create a track proxy:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto tracks = proxy::getCollection<proxy::Tracks>(event, tracksTag);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Here we ask `getCollection()` to create a proxy of category `proxy::Tracks`.
 * Each proxy is a different beast which needs to be explicitly supported:
 * here support for the proxy to `recob::Track` is described.
 *
 * The additional customizations for track proxy are described above.
 * For example, if the module with label `fitTag` stored the fit information for
 * the track, that information can be merged as:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto tracks = proxy::getCollection<proxy::Tracks>
 *   (event, tracksTag, withFitHitInfo(fitTag));
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * while in the more likely case where that information was produced together
 * with the tracks,
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto tracks = proxy::getCollection<proxy::Tracks>
 *   (event, tracksTag, withFitHitInfo());
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will suffice. After this, the interface specific to `recob::TrackFitHitInfo`
 * will be available. Otherwise, any attempt to use it will cause (complicate)
 * compilation errors.
 *
 * In addition, any type of data can be associated using the generic collection
 * proxy interface (`withAssociated()`, `withParallelData()`, etc.).
 *
 *
 * The C++ type of `tracks` object should not matter to the user, and it depends
 * on which additional data are merged in.
 *
 *
 * Types of proxies, and what to do with them
 * ===========================================
 *
 * Currently there are three different type of proxy-like objects for tracks.
 * Each one supports a specific concept:
 *
 * * `proxy::Tracks` represents the whole collection of tracks; it covers
 *   the tracks themselves, and their associated hits. It is obtained by calling
 *   `getCollection()` as described above.
 * * `proxy::Track` represents a single track; the list of hits and points,
 *   and of course the `recob::Track` object itself, can be accessed through it.
 *   Track proxies are obtained from a track collection proxy (`proxy::Tracks`).
 * * `proxy::TrackPoint` represents a single trajectory point in a track.
 *   It provide access to position and momentum of the track at that point,
 *   associated hit and point flags. Track point proxies are obtained from a
 *   track proxy (`proxy::Track`).
 *
 * For the details of the interface and the information that is exposed by each
 * of these proxy classes, please refer to each class documentation. In
 * particular, see `proxy::Tracks` documentation for more usage examples.
 *
 * @note The interface allows by deliberate design only read-only access to the
 *       underlying data.
 *
 *
 * @section LArSoftProxyTracksTech Technical details
 *
 * @subsection LArSoftProxyTracksTechDetails Track collection proxy
 *
 * The track collection proxy object is derived from
 * `proxy::CollectionProxyBase`, which points to the original (track) data
 * product.
 * In addition, it contains a `proxy::details::AssociatedData` object for
 * the `recob::Track`--`recob::Hit` association list.
 *
 * The `proxy::Tracks` interface is currently quite limited: it only allows to
 * access a track by index, or to iterate through all of them, in addition to
 * know how many tracks are available (`size()`).
 *
 * The object returned when accessing the single track information,
 * `proxy::Track`, actually contains enough information so that it is
 * independent of the track collection proxy. This is derived from the generic
 * collection element proxy object (`proxy::CollectionProxyElement`).
 *
 * The object describing the information of a single point has the interface of
 * `TrackPointWrapper`. The underlying storage includes pointers to the track,
 * the associated hit, and the index of the point in the track. This is
 * track-specific and it is not part of the generic proxy infrastructure.
 *
 *
 * @subsection LArSoftProxyTracksCustom Track proxy as an example of proxy customization
 *
 * The proxy utilities provide the basic functionality of the track proxy,
 * including the track collection proxy, representing all the tracks in the
 * event, and the track proxy, representing a single track. Access to a third
 * tier, the single trajectory point, as proxy is not covered by the basic
 * framework, and it has been implemented here from scratch.
 *
 * The most relevant customization pertains proxy to the single track. The proxy
 * is derived from `proxy::details::CollectionProxyElement` for basic
 * functionality, which is enriched by a custom interface that does not in fact
 * add functionality, except for the trajectory point proxy described below.
 * The customizations are fairly trivial, overlaying user-friendly names on the
 * existing functionality. A step beyond that is the access to data that might
 * not be present, that is the fit information. The proxy takes the necessary
 * steps to determine (statically) whether that information is present, and to
 * provide a null pointer to the information if that is not the case. This is
 * a functionality that can't be completely implemented by the basic proxy code,
 * since that generic code has no clue about what type of null pointer to return
 * when the tag `Tracks::TrackFitHitInfoTag` is unknown. This information is
 * provided by the track proxy via the `getIf()` call.
 * For convenience, `proxy::Track` is defined as alias of this single track
 * proxy.
 *
 * The collection proxy customization is way more straightforward.
 * A proxy "tag" is defined, `proxy::Tracks`, with the only purpose of
 * identifying the collection proxy for tracks. For convenience, it hosts a
 * definition used internally to identify one of the optional auxiliary data,
 * and the type of the main data product, but both definitions are contingent
 * and their presence in there is only to centralize some customization in a
 * single place.
 * Traits (`CollectionProxyMakerTraits<Tracks>`) are specialized to inform that
 * this `proxy::Tracks` proxy relies on `std::vector<recob::Track>` as main
 * data product collection type (`std::vector<recob::Track>` is learned from
 * `proxy::Tracks`, but again, this is a contingent detail). The only other
 * customization we need is to have for our proxy our element class above: since
 * we can use the standard collection base, just with the custom element, we
 * define `collection_proxy_impl_t` in that way.
 * Finally, the creation of the collection proxy is customised by specializing
 * `CollectionProxyMaker` (`CollectionProxyMaker<Tracks>`). That class normally
 * takes care of creating the whole proxy, and our purpose is to have it always
 * add the associated hits as auxiliary data, so that the caller does not have
 * to explicitly use `withAssociated<recob::Hit>()` in `getCollection()`.
 * The simple customization does exactly that, under the hood.
 * As candy, some customized functions may be provided for convenience, like
 * `withFitHitInfo()` as alias of `withAssociated<recob::TrackFitHitInfo>()`.
 *
 *
 * The trajectory point proxy has been implemented from scratch here, and it is
 * much less refined than the generic proxy code. It is based on a data
 * structure with a selected list of pointers to the actual data. This structure
 * is implemented as a `std::tuple`. On top of it, a wrapper provides the
 * interface (by interface substitution). This choice is non-essential and has
 * been taken to stress the separation between data storage and interface.
 * Point data structures are created by the proxy on demand, and they are
 * designed so that they don't become invalid when the original proxies fall out
 * of scope, at the price of added memory usage (the minimal information would
 * be a pointer to the track proxy and a point index).
 *
 * Summary of the customization procedure:
 * -# define the tag to identify the proxy (`proxy::Tracks`)
 * -# choose what type of object that will be (`proxy::CollectionProxyBase`
 *      should do for most)
 * -# (_optional_) customize the element type; deriving it from
 *      `proxy::CollectionProxyElement` is recommended
 * -# define the main data product type (`std::vector<recob::Track>`) and set
 *     the traits of the collection proxy, often deriving them from
 *     `proxy::CollectionProxyMakerTraits` with the main data product type as
 *     template argument is enough; in this example, we specified a collection
 *     proxy object with customized element though
 * -# customize the creation of the proxy collection, if special logic or
 *     default components are specified for the proxy (here,
 *     `withAssociated<recob::Hit>()`);
 *     in this case, `proxy::CollectionProxyMaker` must be specialized (for
 *     `proxy::Tracks`), and a starting point may be to derive the
 *     specialization from `proxy::CollectionProxyMakerBase` and redefine its
 *     `make()` member
 *
 *
 * @subsection LArSoftProxyTracksOverhead Overhead
 *
 * See the notes on @ref LArSoftProxyOverhead "overhead" in `ProxyBase.h`.
 *
 *
 */

#ifndef LARDATA_RECOBASEPROXY_TRACK_H
#define LARDATA_RECOBASEPROXY_TRACK_H

// LArSoft libraries
#include "lardata/RecoBaseProxy/ProxyBase.h" // proxy namespace
#include "lardata/Utilities/filterRangeFor.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"

// framework libraries
#include "canvas/Persistency/Common/Ptr.h"

#include <limits>
#include <tuple>
#include <vector>

namespace proxy {

  //----------------------------------------------------------------------------
  // forward declarations
  template <typename TrackProxy>
  class TrackPointIterator;

  template <typename Data>
  class TrackPointWrapper;

  //----------------------------------------------------------------------------
  namespace details {

    template <typename CollProxy>
    struct TrackPointIteratorBox;

    template <typename T>
    struct isTrackProxy;

    template <typename Obj>
    struct StaticAsserts;

    template <typename Data>
    struct StaticAsserts<TrackPointWrapper<Data>> : public std::true_type {
      using Wrapper_t = TrackPointWrapper<Data>;

      static_assert(sizeof(Wrapper_t) == 1U, "Wrapper carries data!");

      static_assert(std::is_same<std::decay_t<decltype(std::declval<Wrapper_t>().position())>,
                                 recob::Track::Point_t>(),
                    "position() is not a recob::Track::Point_t");
      static_assert(std::is_same<std::decay_t<decltype(std::declval<Wrapper_t>().momentum())>,
                                 recob::Track::Vector_t>(),
                    "momentum() is not a recob::Track::Vector_t");
      static_assert(std::is_same<std::decay_t<decltype(std::declval<Wrapper_t>().flags())>,
                                 recob::Track::PointFlags_t>(),
                    "flags() is not a recob::Track::PointFlags_t");
      static_assert(std::is_same<std::decay_t<decltype(std::declval<Wrapper_t>().hitPtr())>,
                                 art::Ptr<recob::Hit>>(),
                    "hit() is not a art::Ptr<recob::Hit>");
      static_assert(
        std::is_same<std::decay_t<decltype(std::declval<Wrapper_t>().index())>, std::size_t>(),
        "index() is not a std::size_t");

    }; // StaticAsserts<TrackPointWrapper<Data>>

  } // namespace details

  //----------------------------------------------------------------------------

  /**
   * @brief Proxy tag for a `recob::Track` collection proxy.
   * @see `proxy::TrackCollectionProxyElement`
   * @ingroup LArSoftProxyTracks
   *
   * This type can be used to get a proxy for `recob::Track` collection:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto tracks = proxy::getCollection<proxy::Tracks>(event, tracksTag);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * An example of usage for a simple track processing loop:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void MyAnalyzer::analyze(art::Event const& event) {
   *
   *   auto tracks = proxy::getCollection<proxy::Tracks>
   *     (event, tracksTag, proxy::withFitHitInfo());
   *
   *   if (tracks.empty()) {
   *     mf::LogVerbatim("TrackProxyTest")
   *       << "No tracks in '" << tracksTag.encode() << "'";
   *     return;
   *   }
   *
   *   mf::LogVerbatim("TrackProxyTest") << "Collection '" << tracksTag.encode()
   *     << "' contains " << tracks.size() << " tracks.";
   *
   *   for (auto track: tracks) {
   *
   *     recob::Track const& trackRef = track.track();
   *
   *     mf::LogVerbatim log("TrackProxyTest");
   *     log << "[#" << track.index() << "] track " << trackRef
   *       << "\n  with " << trackRef.NPoints() << " points and "
   *       << track.nHits() << " hits:";
   *
   *     for (auto const& point: track.points()) {
   *       log <<
   *         "\n  [#" << point.index() << "] at " << point.position()
   *           << " (momentum: " << point.momentum() << "), flags: "
   *           << point.flags();
   *
   *       recob::Hit const* hit = point.hit();
   *       if (hit) {
   *         log << " with a Q=" << hit->Integral() << " hit on channel "
   *           << hit->Channel() << " at tick " << hit->PeakTime()
   *           << ", measured: " << point.fitInfoPtr()->hitMeas();
   *       }
   *       else
   *         log << " (no associated hit)";
   *
   *     } // for points in track
   *
   *   } // for track
   *
   * } // MyAnalyzer::analyze()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In this example, the track proxy accesses the track itself, its associated
   * hits (always implicitly present) and the track fit hit information
   * (explicitly requested). Since both those data products are produced by the
   * same module as the track, there is no need to specify their producer module
   * label.
   *
   * Unfortunately, the proxy object (`tracks` in the example) can be of a
   * different class depending on which data is merged into it: a proxy created
   * by `getCollection<proxy::Tracks>(event, tag, proxy::withFitHitInfo())` has
   * different type than e.g. `getCollection<proxy::Tracks>(event, tag)`.
   * This implies than when passing proxies as arguments to functions, template
   * types must be used. For example, the following code is equivalent to the
   * one above, but with methods processing a single track (a track proxy) and
   * a single trajectory point (a track point proxy):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename TrackPoint>
   * void MyAnalyzer::processPoint(TrackPoint const& point) const {
   *
   *   mf::LogVerbatim log("TrackProxyTest");
   *   log <<
   *     "  [#" << point.index() << "] at " << point.position()
   *       << " (momentum: " << point.momentum() << "), flags: "
   *       << point.flags();
   *
   *   recob::Hit const* hit = point.hit();
   *   if (hit) {
   *     log << " with a Q=" << hit->Integral() << " hit on channel "
   *       << hit->Channel() << " at tick " << hit->PeakTime()
   *       << ", measured: " << point.fitInfoPtr()->hitMeas();
   *   }
   *   else
   *     log << " (no associated hit)";
   *
   * } // MyAnalyzer::processPoint()
   *
   *
   * //------------------------------------------------------------------------------
   * template <typename Track>
   * void MyAnalyzer::processTrack(Track const& track) const {
   *
   *   recob::Track const& trackRef = track.track();
   *
   *   mf::LogVerbatim("TrackProxyTest")
   *     << "[#" << track.index() << "] track " << trackRef
   *     << "\n  with " << trackRef.NPoints() << " points and " << track.nHits()
   *     << " hits:";
   *
   *   for (auto point: track.points()) {
   *     processPoint(point);
   *   } // for points in track
   *
   * } // MyAnalyzer::processTrack()
   *
   *
   * //------------------------------------------------------------------------------
   * void MyAnalyzer::proxyUsageExample(art::Event const& event) {
   *
   *   auto tracks = proxy::getCollection<proxy::Tracks>
   *     (event, tracksTag, proxy::withFitHitInfo());
   *
   *   if (tracks.empty()) {
   *     mf::LogVerbatim("TrackProxyTest") << "No tracks in '"
   *       << tracksTag.encode() << "'";
   *     return;
   *   }
   *
   *   mf::LogVerbatim("TrackProxyTest") << "Collection '" << tracksTag.encode()
   *     << "' contains " << tracks.size() << " tracks.";
   *
   *   for (auto track: tracks) {
   *
   *     processTrack(track);
   *
   *   } // for track
   *
   * } // MyAnalyzer::proxyUsageExample()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * As any other proxy object, other data can be merged to the proxy, but no
   * custom interface will be available. For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto tracks = proxy::getCollection<proxy::Tracks>(
   *   event, tracksTag,
   *   proxy::withParallelData<recob::TrackMomentumFit>(momTag)
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will add a data product `std::vector<recob::TrackMomentumFit>` expected to
   * have one element (of type `recob::TrackMomentumFit`) per track, which will
   * be accessed with the generic proxy interface:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * for (auto track: tracks) {
   *   recob::TrackMomentumFit const& momFit
   *     = track.get<recob::TrackMomentumFit>();
   *   // ...
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The other common features of proxy collections are supported, like tagging
   * of different instances of the same data types:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct MCS {};
   * struct Range {};
   * auto tracks = proxy::getCollection<proxy::Tracks>(
   *   event, tracksTag,
   *   proxy::withParallelDataAs<recob::TrackMomentumFit, Range>(rangeMomTag),
   *   proxy::withParallelDataAs<recob::TrackMomentumFit, MCS>(MCStag)
   *   );
   * for (auto track: tracks) {
   *   recob::TrackMomentumFit const& rangeMom = track.get<Range>();
   *   recob::TrackMomentumFit const& MCSmom = track.get<MCS>();
   *   // ...
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * A new, filtered collection of proxies can be created with obvious means and
   * with a less-than-friendly declaration:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<decltype(tracks)::element_proxy_t> longTracks;
   * for (auto track: tracks) {
   *   if (track->Length() >= 30.0) longTracks.push_back(track);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The collection thus created (`longTracks`) is valid also after the
   * collection proxy (`tracks`) has fallen out of scope.
   *
   *
   * @note `proxy::Tracks` is *not* the type of the collection proxy returned
   *       by `getCollection()`.
   */
  struct Tracks {

    /// Type of the main collection.
    using TrackDataProduct_t = std::vector<recob::Track>;

    /// Tag used for the "standard" track trajectory information.
    using TrackTrajectoryTag = recob::TrackTrajectory;

    /// Tag used for the "standard" track fit information.
    using TrackFitHitInfoTag = recob::TrackFitHitInfo;

    /// Tag used for the associated hits.
    using HitTag = recob::Hit;

    /// Types of tracks and trajectories.
    typedef enum {
      Unfitted, ///< Represents a track trajectory before the final fit.
      Fitted,   ///< Represents a track trajectory from the final fit.
      NTypes    ///< Number of supported track types.
    } TrackType_t;

  }; // struct Tracks

  //----------------------------------------------------------------------------
  //---  track point information
  //---
  /**
   * @brief Container of track point information.
   * @see `proxy::Track`, `proxy::TrackPointWrapper`
   * @ingroup LArSoftProxyTracks
   *
   * This class contains some information pertaining a single point of a
   * `recob::Track`.
   *
   * The information is not extensible via the usual proxy mechanisms. The data
   * supported is stored in `proxy::TrackPointData` class, and it currently
   * includes:
   * * the `recob::Track` the point belongs to
   * * the position, momentum and flags of the point
   * * the hit associated to the point
   * * the index of the point in its track
   * * fit information
   *
   * The access interface is determined by `proxy::TrackPointWrapper`.
   *
   * An object of this type is returned when accessing a single track point from
   * a track proxy. While the data itself is not copied, this object owns some
   * pointers to the actual data, and once created it is independent of the
   * track proxy which created it.
   *
   */
  using TrackPointData = std::
    tuple<recob::Track const*, art::Ptr<recob::Hit>, recob::TrackFitHitInfo const*, std::size_t>;

  /**
   * @brief Wrapper for a track data proxy.
   * @tparam Data the point data type; requirements are described below
   * @ingroup LArSoftProxyTracks
   *
   * This class provides a user interface to data pertaining a single trajectory
   * point of a `recob::Track`.
   *
   *
   * Implementation details
   * -----------------------
   *
   * This type wraps a generic data structure with tuple interface.
   * It is expected that the following information is returned:
   *
   * * pointer to point position as `recob::Track::Point_t` from
   *   `std::get<0>(Data const&)`
   * * pointer to momentum at point as `recob::Track::Vector_t` from
   *   `std::get<1>(Data const&)`
   * * pointer to flags as `recob::Track::PointFlags_t` from
   *   `std::get<2>(Data const&)`
   * * pointer to associated hit as `art::Ptr<recob::Hit>` from
   *   `std::get<3>(Data const&)`
   *
   * The "pointers" can be any object returning the required type when
   * dereferenced.
   *
   */
  template <typename Data>
  class TrackPointWrapper {
    using This_t = TrackPointWrapper<Data>;
    using Wrapped_t = std::add_const_t<Data>;

    static constexpr std::size_t TrackIndex = 0;
    static constexpr std::size_t HitIndex = 1;
    static constexpr std::size_t FitHitInfoIndex = 2;
    static constexpr std::size_t IndexIndex = 3;
    static constexpr std::size_t NIndices = 4;

    static_assert(std::tuple_size<Data>::value == NIndices, "Unexpected data size.");

    Wrapped_t const& base() const { return reinterpret_cast<Wrapped_t const&>(*this); }

    template <std::size_t N>
    auto get() const -> decltype(auto)
    {
      return std::get<N>(base());
    }

  protected:
    TrackPointWrapper() = default;
    TrackPointWrapper(TrackPointWrapper const&) = default;
    TrackPointWrapper(TrackPointWrapper&&) = default;
    TrackPointWrapper& operator=(TrackPointWrapper const&) = default;
    TrackPointWrapper& operator=(TrackPointWrapper&&) = default;

  public:
    /// Returns the track this point belongs to.
    recob::Track const& track() const { return *get<TrackIndex>(); }

    /// Returns the position of the trajectory point.
    /// @see `recob::Track::LocationAtPoint()`
    auto position() const -> decltype(auto)
    {
      return track().Trajectory().LocationAtPoint(index());
    }

    /// Returns the momentum vector of the trajectory point.
    /// @see `recob::Track::MomentumVectorAtPoint()`
    auto momentum() const -> decltype(auto)
    {
      return track().Trajectory().MomentumVectorAtPoint(index());
    }

    /**
     * @{
     * @name Flags interface
     */

    /// Returns the flags associated with the trajectory point.
    /// @see `recob::Track::FlagsAtPoint()`
    auto flags() const -> decltype(auto) { return track().Trajectory().FlagsAtPoint(index()); }

    /**
     * @brief Returns whether the trajectory point is valid.
     *
     * Even if the trajectory point (position and momentum) are not valid,
     * the hit is still associated to the track/tracjectory.
     */
    bool isPointValid() const { return flags().isPointValid(); }

    /// @}

    /**
     * @brief Returns the hit associated with the trajectory point.
     * @return an _art_ pointer to the hit associated to this point
     */
    art::Ptr<recob::Hit> hitPtr() const { return get<HitIndex>(); }

    /**
     * @brief Returns fit info associated with the trajectory point.
     * @return a pointer to the fit info, or `nullptr` if not merged in proxy
     *
     * If the track proxy this point comes from had no fit information,
     * `nullptr` is returned.
     * The fit information is extracted using the tag in
     * `proxy::Tracks::TrackFitHitInfoTag`.
     */
    recob::TrackFitHitInfo const* fitInfoPtr() const { return get<FitHitInfoIndex>(); }

    /// Returns the index of this point in the trajectory.
    auto index() const -> decltype(auto) { return get<IndexIndex>(); }

    /// Returns a pointer to the hit on the trajectory point, if any.
    recob::Hit const* hit() const
    {
      decltype(auto) ptr = hitPtr();
      return ptr ? ptr.get() : nullptr;
    }

  }; // TrackPointWrapper<>

  /**
   * @brief Type of track point information.
   * @ingroup LArSoftProxyTracks
   * @see `proxy::TrackPointWrapper`
   *
   * For its interface, see `proxy::TrackPointWrapper`.
   */
  struct TrackPoint : private TrackPointData, public TrackPointWrapper<TrackPointData> {
    using TrackPointData::TrackPointData;
    TrackPoint(TrackPointData const& data) : TrackPointData(data) {}
    TrackPoint(TrackPointData&& data) : TrackPointData(std::move(data)) {}

  private:
    static constexpr bool asserts =
      details::StaticAsserts<TrackPointWrapper<TrackPointData>>::value;
  }; // class TrackPoint

  /**
   * @brief Returns an object with information about the specified track point.
   * @tparam TrackProxy an instance of proxy::Track template
   * @param track the track (proxy) the points belong to
   * @param index the index of the point within the track
   * @return a `TrackPointData` object with information on that point
   *
   * For an interface to the point information, see `TrackPointWrapper`.
   */
  template <typename TrackProxy>
  TrackPointData makeTrackPointData(TrackProxy const& track, std::size_t index)
  {
    static_assert(details::isTrackProxy<TrackProxy>(), "Not a proxy::Track!");
    return {&(track.track()), track.hitAtPoint(index), track.fitInfoAtPoint(index), index};
  } // makeTrackPointData()

  //--------------------------------------------------------------------------
  /**
    * @brief Class for track proxy elements.
    * @tparam CollProxy type of track proxy collection to get data from
    * @see `proxy::TrackPoint`, `proxy::TrackPointWrapper`
    * @ingroup LArSoftProxyTracks
    *
    * For details on the track point interface see `proxy::TrackPoint`.
    */
  template <typename CollProxy>
  struct TrackCollectionProxyElement : public CollectionProxyElement<CollProxy> {
    using base_t = CollectionProxyElement<CollProxy>; ///< Base type.
    using base_t::base_t;                             // inherit constructors

    /// This type.
    using track_proxy_t = TrackCollectionProxyElement<CollProxy>;

  public:
    /// Iterator for trajectory point information.
    using point_iterator = TrackPointIterator<track_proxy_t>;

    /// Returns the pointed track.
    recob::Track const& track() const { return base_t::operator*(); }

    /**
     * @brief Returns the requested trajectory from the proxy.
     * @param type type of the track trajectory to be returned
     * @return a reference to the requested track trajectory
     *
     */
    recob::TrackTrajectory const* operator()(proxy::Tracks::TrackType_t type) const noexcept;

    // --- BEGIN Direct hit interface ------------------------------------------
    /**
     * @name Direct hit interface.
     *
     * The track prescription requires one hit per trajectory point.
     *
     * @note Remember that in particular cases there might be a hit without
     *       point of vice versa. In those cases, the point will have a dummy
     *       value, or the hit pointer will have `isNull()` true. In the former
     *       case, the point flag `isPointValid()` should be unset.
     *
     * The interface at track proxy level allows for both access to the whole
     * sequence of hits, or to the hit of a specific point:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * assert(track.nHits() > 0);
     * auto const& hits = track.hits();
     * art::Ptr<recob::Hit> maxHit = hits[0]; // direct access
     * for (art::Ptr<recob::Hit> const& hit: hits) {
     *   if (hit.isNull()) continue;
     *   if (maxHit.isNull() || maxHit->Charge() < hit->Charge())
     *     maxHit = hit;
     * } // for
     *
     * art::Ptr<recob::Hit> lastHit = track.hitAtPoint(track.nHits() - 1U);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    /// @{

    /**
     * @brief Returns a collection-like range of hits of this track, at point
     *        order.
     * @return a range of _art_ pointers to hits
     *
     * One hit is expected per trajectory point. Hits can be missing, in which
     * case the art pointer will have `isNull()` as `true`.
     */
    auto hits() const -> decltype(auto) { return base_t::template get<Tracks::HitTag>(); }

    /// Returns an art pointer to the hit associated with the specified point.
    auto hitAtPoint(std::size_t index) const -> decltype(auto) { return hits()[index]; }

    /// Returns the number of hits associated with this track.
    std::size_t nHits() const { return hits().size(); }

    /// @}
    // --- END Direct hit interface --------------------------------------------

    /// Returns fit info for the specified point (`nullptr` if not available).
    recob::TrackFitHitInfo const* fitInfoAtPoint(std::size_t index) const;

    // --- BEGIN Direct track trajectory interface -----------------------------
    /**
     * @name Direct track trajectory interface
     * @see `proxy::TrackPoint`
     *
     * The interface allows to check if this track has a trajectory associated
     * with it, and to obtain a reference to it or its _art_ pointer.
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * bool hasTraj = proxy.hasOriginalTrajectory();
     * if (hasTraj) {
     *   recob::TrackTrajectory const& trajectory = proxy.originalTrajectory();
     *   // ...
     * }
     * art::Ptr<recob::TrackTrajectory> const& trajectoryPtr
     *   = proxy.originalTrajectoryPtr();
     * if (!trajectoryPtr.isNull()) {
     *   recob::TrackTrajectory const& trajectory = *trajectoryPtr;
     *   // ...
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @note This interface can't be used if the track trajectory information
     *       has not been merged into the proxy (typically via
     *       `proxy::withOriginalTrajectory()`).
     */
    /// @{

    /// Returns whether this track is associated to a trajectory.
    bool hasOriginalTrajectory() const { return !originalTrajectoryPtr().isNull(); }

    /// Returns an _art_ pointer to the associated trajectory.
    /// @return pointer to the associated trajectory (`isNull()` `true` if none)
    art::Ptr<recob::TrackTrajectory> const& originalTrajectoryPtr() const
    {
      return base_t::template get<Tracks::TrackTrajectoryTag>();
    }

    /**
     * @brief Returns a reference to the associated trajectory.
     * @return the associated trajectory as a constant reference
     * @see `originalTrajectoryPtr()`, `hasOriginalTrajectory()`
     *
     * If the track is not associated to any trajectory, the return value is
     * undefined. This condition should be checked beforehand, e.g. with
     * `hasTrajectory()`.
     */
    recob::TrackTrajectory const& originalTrajectory() const { return *originalTrajectoryPtr(); }

    /// @}
    // --- END Direct track trajectory interface -------------------------------

    // --- BEGIN Point-by-point iteration interface ----------------------------
    /**
     * @name Point-by-point iteration interface
     *
     * The points on track can be accessed individually with a special,
     * non-extensible proxy.
     *
     * In this example, points are accessed via iteration:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * unsigned int nPoints = track.nPoints();
     * unsigned int nValidHits = 0;
     * for (auto point: track.points()) {
     *   if (!point.hit().isNull()) ++nValidHits;
     * }
     * unsigned int nValidPoints = std::count_if
     *   (track.beginPoint(), track.endPoint(), &TrackPoint::isPointValid);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Random (index-based) access is also available:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * if (track.nPoints() > 1) {
     *   auto point = track.point(1); // or track[1]
     *   // ...
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    /// @{

    /**
     * @brief Returns an iterable range with point-by-point information.
     * @see `proxy::TrackPoint`, `proxy::TrackPointWrapper`
     *
     * The interface of the elements is documented in `TrackPointWrapper`.
     * Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (auto const& pointInfo: track.points()) {
     *
     *   if (!pointInfo.flags().isPointValid()) continue;
     *
     *   auto const& pos = pointInfo.position();
     *
     *   // ...
     *
     * } // for point
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will iterate through all points (including the invalid ones, hence the
     * check).
     */
    auto points() const { return details::TrackPointIteratorBox<CollProxy>(*this); }

    /**
     * @brief Returns an iterable range with only points matching the `mask`.
     * @tparam Pred type of predicate to test on points
     * @param pred predicate to be fulfilled by the points
     * @return an object that can be forward-iterated
     * @see `points()`, `pointsWithFlags()`
     *
     * This methods is used in a way similar to `points()`, with the addition of
     * specifying a criterium (predicate) defining the selected points.
     * The iteration will happen only through the points which fulfil the
     * predicate.
     *
     * The interface of the elements is documented in `TrackPointWrapper`.
     * Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto farPoints = [](auto const& pointInfo)
     *   { return pointInfo.isPointValid() && pointInfo.position().Z() > 50.; };
     * for (auto const& pointInfo: track.selectPoints(farPoints)) {
     *
     *   auto const& pos = pointInfo.position();
     *
     *   // ...
     *
     * } // for point
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will iterate through all points which are valid and whose position is
     * at _z_ absolute coordinate larger than 50 centimeters (whatever it
     * means).
     *
     *
     * Requirements
     * -------------
     *
     * * `Pred` is a unary function object which can accept a `TrackPoint`
     *     object as its sole argument and which returns a value convertible to
     *     `bool`
     *
     */
    template <typename Pred>
    auto selectPoints(Pred&& pred) const;

    /**
     * @brief Returns an iterable range with only points matching the `mask`.
     * @param mask point flag mask to be matched
     * @return an object that can be forward-iterated
     * @see `points()`, `util::flags::BitMask::match()`
     *
     * This methods is used in a way similar to `points()`, with the addition of
     * specifying a `mask` of flags. The iteration will happen only through the
     * points which match the mask. that is for which
     * `pointInfo.flags().match(mask)` is `true`.
     *
     * The interface of the elements is documented in `TrackPointWrapper`.
     * Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * for (auto const& pointInfo
     *   : track.pointsWithFlags(-recob::TrajectoryPointFlags::flag::NoPoint)
     *   )
     * {
     *
     *   auto const& pos = pointInfo.position();
     *
     *   // ...
     *
     * } // for point
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will iterate through only the points which _do not_ have the `NoPoint`
     * flag set (which have in fact a valid position).
     */
    auto pointsWithFlags(recob::TrackTrajectory::PointFlags_t::Mask_t mask) const;

    /// Returns the number of trajectory points in the track.
    std::size_t nPoints() const { return track().NPoints(); }

    /// Returns the iterator to the data of the first point.
    point_iterator beginPoint() const { return {*this, 0}; }

    /// Returns the iterator past the last point.
    point_iterator endPoint() const { return {*this, nPoints()}; }

    /// @{
    /// Extracts information from the specified point.
    TrackPoint point(std::size_t index) const { return {makeTrackPointData(track(), index)}; }

    TrackPoint operator[](std::size_t index) const { return point(index); }

    /// @}
    // --- END Point-by-point iteration interface ------------------------------

    // --- BEGIN Additional utilities ------------------------------------------
    /// @name Additional utilities
    /// @{

    /// @}
    // --- END Additional utilities --------------------------------------------

  private:
    recob::TrackTrajectory const* originalTrajectoryCPtr() const noexcept
    {
      return hasOriginalTrajectory() ? &originalTrajectory() : nullptr;
    }

  }; // TrackCollectionProxyElement<>

  /**
   * @brief Proxy to an element of a proxy collection of `recob::Track` objects.
   * @tparam TrackCollProxy type of the track collection proxy
   * @ingroup LArSoftProxyTracks
   *
   * This class is the proxy equivalent of `recob::Track`, which exposes data
   * associated with the track.
   * An object of this type is returned when accessing a single track from a
   * track collection proxy. While the data itself is not copied, this object
   * owns some pointers to the actual data, and once created it is independent
   * of the track collection proxy which created it.
   *
   * The interface is currently defined by
   * `proxy::TrackCollectionProxyElement`.
   */
  template <typename TrackCollProxy>
  using Track = TrackCollectionProxyElement<TrackCollProxy>;

  // --- BEGIN Auxiliary data --------------------------------------------------
  /**
   * @name Auxiliary data
   *
   * These functions may be used as arguments to
   * `proxy::getCollection<proxy::Tracks>()` call to
   * @ref LArSoftProxyDefinitionMerging "merge" of some data associated to the
   * tracks.
   *
   * @{
   */

  /**
   * @brief Adds `recob::TrackTrajectory` information to the proxy.
   * @param inputTag the data product label to read the data from
   * @return an object driving `getCollection()` to use `recob::TrackTrajectory`
   * @ingroup LArSoftProxyTracks
   * @see `proxy::withOriginalTrajectory()`,
   *      `proxy::Tracks`, `proxy::getCollection()`
   *
   * The behaviour of this function is like `withOriginalTrajectory()`, but
   * reading the original trajectories from the association with the specified
   * label rather than the label of the tracks in the proxy.
   */
  inline auto withOriginalTrajectory(art::InputTag const& inputTag)
  {
    return proxy::withZeroOrOneAs<recob::TrackTrajectory, Tracks::TrackTrajectoryTag>(inputTag);
  }

  /**
   * @brief Adds `recob::TrackTrajectory` information to the proxy.
   * @return an object driving `getCollection()` to use `recob::TrackTrajectory`
   * @ingroup LArSoftProxyTracks
   * @see `proxy::withOriginalTrajectory(art::InputTag const&)`,
   *      `proxy::Tracks`,
   *      `proxy::TrackCollectionProxyElement::hasOriginalTrajectory()`,
   *      `proxy::TrackCollectionProxyElement::originalTrajectory()`,
   *      `proxy::TrackCollectionProxyElement::originalTrajectoryPtr()`
   *
   * The information from the associated trajectories is merged in the proxy.
   * That association data product must have the same input tag as the track
   * collection data product. To specify a different one, use
   * `withOriginalTrajectory(art::InputTag const&)`.
   *
   * The data is available through the regular interface via tag
   * `recob::TrackTrajectory`, or via custom interface, e.g.:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto tracks = proxy::getCollection<proxy::Tracks>
   *   (event, tracksTag, proxy::withOriginalTrajectory());
   *
   * for (auto const& trackProxy: tracks) {
   *
   *   if (!trackProxy.hasOriginalTrajectory()) continue;
   *
   *   const auto& track = *trackProxy;
   *   recob::TrackTrajectory const& original = trackProxy.originalTrajectory();
   *   recob::TrackTrajectory const& fitted = track.Trajectory();
   *
   *   // ...
   *
   * } // for tracks
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Note the subtle difference between the two inner lines: `trackProxy` is a
   * track proxy, and the first line is accessing its interface. The second line
   * is talking to the track object (`recob::Track`, actually) directly. The
   * same effect is obtained using directly the proxy, but with the indirection
   * operator (`->`) instead of the member operator (`.`):
   * `trackProxy->Trajectory()`.
   *
   * The `recob::TrackTrajectory` information is required to be from a _art_
   * association with `recob::Track`. The association must fulfil the
   * @ref LArSoftProxyDefinitionOneToZeroOrOneSeqAssn "one-to-(zero-or-one) sequential association"
   * requirements.
   */
  inline auto withOriginalTrajectory()
  {
    return proxy::withZeroOrOneAs<recob::TrackTrajectory, Tracks::TrackTrajectoryTag>();
  }

  //----------------------------------------------------------------------------
  /**
   * @brief Adds `recob::TrackFitHitInfo` information to the proxy.
   * @param inputTag the data product label to read the data from
   * @return an object driving `getCollection()` to use `recob::TrackFitHitInfo`
   * @ingroup LArSoftProxyTracks
   * @see `proxy::Tracks`, `proxy::getCollection()`, `proxy::withFitHitInfo()`
   *
   * This function behaves like `withFitHitInfo()`, but allows to use `inputTag`
   * as input tag, instead of the same label as for the track collection.
   * See `proxy::withFitHitInfo()` for explanations and examples.
   */
  inline auto withFitHitInfo(art::InputTag const& inputTag)
  {
    return proxy::withParallelDataAs<std::vector<recob::TrackFitHitInfo>,
                                     Tracks::TrackFitHitInfoTag>(inputTag);
  }

  /**
   * @brief Adds `recob::TrackFitHitInfo` information to the proxy.
   * @return an object driving `getCollection()` to use `recob::TrackFitHitInfo`
   * @ingroup LArSoftProxyTracks
   * @see `proxy::withFitHitInfo(art::InputTag const&)`,
   *      `proxy::Tracks`, `proxy::getCollection()`
   *
   * A `recob::TrackFitHitInfo`data product is read from the event and merged
   * into the proxy being created by `proxy::getCollection()`.
   * The data product has the same input tag as the track data product; if a
   * different one is needed, use `proxy::withFitHitInfo(art::InputTag const&)`
   * instead.
   *
   * Example of usage (more can be found in `TrackProxyTest::testTracks()`):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto tracks = proxy::getCollection<proxy::Tracks>
   *   (event, tracksTag, proxy::withFitHitInfo());
   *
   * for (auto const& trackInfo: tracks) {
   *
   *   for (auto const& point: track.points()) {
   *
   *     auto const& pos = point.position();
   *     auto const* hit = point.hit();
   *     auto const* fitInfo = point.fitInfoPtr();
   *
   *     // ...
   *
   *   } // for point
   *
   * } // for tracks
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The proxy helps associating the right set of `recob::TrackFitHitInfo` for
   * each `track` in the outer loop (not shown in the example: it might have
   * looked like `auto const& fitInfo = tracks.get<recob::TrackFitHitInfo>()`).
   * It also helps to pick the fit information of the current `point` in the
   * inner loop of the example (`fitInfo` will never be `nullptr` here, since
   * we _did_ merge the fit information).
   *
   * The collection of `recob::TrackFitHitInfo` is required to be a
   * `std::vector<std::vector<recob::TrackFitHitInfo>>`, where the first index
   * addresses which track the information is about, and the second index which
   * point within that track.
   *
   * The data is also available through the regular interface via tag
   * `recob::TrackFitHitInfo`.
   *
   * The data must satisfy the
   * @ref LArSoftProxyDefinitionParallelData "parallel data product"
   * requirement.
   */
  inline auto withFitHitInfo()
  {
    return proxy::withParallelDataAs<std::vector<recob::TrackFitHitInfo>,
                                     Tracks::TrackFitHitInfoTag>();
  }

  /// @}
  // --- END Auxiliary data ----------------------------------------------------

  //----------------------------------------------------------------------------
  /// Define the traits of `proxy::Tracks` proxy.
  template <>
  struct CollectionProxyMakerTraits<Tracks>
    : public CollectionProxyMakerTraits<Tracks::TrackDataProduct_t> {
    // default traits, plus a collection proxy class with a custom element:
    template <typename MainColl, typename... AuxColl>
    using collection_proxy_impl_t = CollectionProxyBase<Track, MainColl, AuxColl...>;
  };

  //----------------------------------------------------------------------------
  /// Specialization to create a proxy for `recob::Track` collection.
  template <>
  struct CollectionProxyMaker<Tracks> : public CollectionProxyMakerBase<Tracks> {

    /// Base class.
    using maker_base_t = CollectionProxyMakerBase<Tracks>;

    /**
     * @brief Creates and returns a collection proxy for `recob::Track` based on
     *        `proxy::Tracks` tag and with the requested associated data.
     * @tparam Event type of the event to read the information from
     * @tparam WithArgs type of arguments for associated data
     * @param event event to read the information from
     * @param tag input tag of the `recob::Track` collection data product
     * @param withArgs optional associated objects to be included
     * @return a collection proxy to `recob::Track` collection with `tag`
     *
     * For each argument in `withArgs`, an action is taken. Usually that is to
     * add an association to the proxy.
     * Associated hits (tag: `recob::Hit`) are automatically added to the proxy
     * and must not be explicitly specified.
     */
    template <typename Event, typename... WithArgs>
    static auto make(Event const& event, art::InputTag const& tag, WithArgs&&... withArgs)
    {
      // automatically add associated hits with the same input tag;
      // IDEA: allow a withAssociated<recob::Hit>() from withArgs to override
      // this one; the pattern may be:
      // - if withArgs contains a withAssociated<recob::Hit>(), produce a new
      //   withArgs with that one pushed first
      // - otherwise, produce a new withArgs with a new
      //   withAssociated<recob::Hit>(tag) as first element
      // In principle there is no need for these hits to be first; code might
      // be simpler when assuming that though.
      return maker_base_t::make(event,
                                tag,
                                withAssociatedAs<recob::Hit, Tracks::HitTag>(),
                                std::forward<WithArgs>(withArgs)...);
    } // make()

  }; // struct CollectionProxyMaker<>

  /// "Converts" point data into a `proxy::TrackPointWrapper`.
  template <typename Data>
  auto wrapTrackPoint(Data const& wrappedData)
  {
    (void)details::StaticAsserts<TrackPointWrapper<Data>>();
    return reinterpret_cast<TrackPointWrapper<Data> const&>(wrappedData);
  }

  /// Iterator for points of a track proxy. Only supports range-for loops.
  /// @ingroup LArSoftProxyReco
  template <typename TrackProxy>
  class TrackPointIterator {

    /*
     * So, let's go through the list of iterator traits from cppreference.com:
     * [x] Iterator
     *   [x] CopyConstructible
     *   [x] CopyAssignable
     *   [x] Destructible
     *   [x] lvalues are Swappable
     *   [x] value_type
     *   [x] difference_type
     *   [x] reference
     *   [x] pointer
     *   [x] iterator_category
     *   [x] operator*()
     *   [x] operator++()
     * [ ] InputIterator
     *   [x] Iterator (above)
     *   [x] EqualityComparable (operator== (A, B))
     *   [x] operator!= ()
     *   [ ] reference operator*() (convertible to value_type)
     *   [ ] operator->()
     *   [x] It& operator++()
     *   [x] operator++(int)
     *   [x] *i++ equivalent to { auto v = *i; ++i; return v; }
     * [ ] Forward Iterator
     *   [ ] InputIterator (above)
     *   [x] DefaultConstructible
     *   [x] multipass guarantee: a == b => ++a == ++b
     *   [ ] reference = value_type const&
     *   [x] It operator++(int)
     *   [ ] *i++ returns reference
     * That's it! :-|
     */

    using track_proxy_t = TrackProxy;

    track_proxy_t const* track = nullptr;
    std::size_t index = std::numeric_limits<std::size_t>::max();

  public:
    /// @name Iterator traits
    /// @{
    using difference_type = std::ptrdiff_t;
    using value_type = TrackPoint;
    using pointer = TrackPoint const*;
    using reference = TrackPoint; // booo!
    // not quite an input iterator (see above)
    using iterator_category = std::input_iterator_tag;
    /// @}

    TrackPointIterator() = default;

    TrackPointIterator(track_proxy_t const& track, std::size_t index) : track(&track), index(index)
    {}

    TrackPointIterator& operator++()
    {
      ++index;
      return *this;
    }

    TrackPointIterator operator++(int)
    {
      auto it = *this;
      this->operator++();
      return it;
    }

    // we make sure the return value is a temporary
    value_type operator*() const
    {
      return static_cast<value_type>(makeTrackPointData(*track, index));
    }

    bool operator==(TrackPointIterator const& other) const
    {
      return (index == other.index) && (track == other.track);
    }

    bool operator!=(TrackPointIterator const& other) const
    {
      return (index != other.index) || (track != other.track);
    }

  }; // class TrackPointIterator

} // namespace proxy

namespace proxy {

  //----------------------------------------------------------------------------
  namespace details {

    //--------------------------------------------------------------------------
    template <typename T>
    struct isTrackProxy : public std::false_type {};

    template <typename TrackCollProxy>
    struct isTrackProxy<Track<TrackCollProxy>> : public std::true_type {};

    //--------------------------------------------------------------------------
    /// Structure for range-for iteration.
    template <typename CollProxy>
    struct TrackPointIteratorBox {
      using track_proxy_t = Track<CollProxy>;
      using const_iterator = typename track_proxy_t::point_iterator;

      TrackPointIteratorBox(track_proxy_t const& track) : track(&track) {}

      const_iterator begin() const { return track->beginPoint(); }

      const_iterator end() const { return track->endPoint(); }

    private:
      track_proxy_t const* track = nullptr;

    }; // TrackPointIteratorBox<>

    //--------------------------------------------------------------------------

  } // namespace details

  //----------------------------------------------------------------------------
  template <typename CollProxy>
  recob::TrackTrajectory const* TrackCollectionProxyElement<CollProxy>::operator()(
    proxy::Tracks::TrackType_t type) const noexcept
  {
    switch (type) {
    case proxy::Tracks::Fitted: return &(track().Trajectory());
    case proxy::Tracks::Unfitted: return originalTrajectoryCPtr();
    default: return nullptr;
    } // switch
  }   // TrackCollectionProxyElement<>::operator()

  //----------------------------------------------------------------------------
  template <typename CollProxy>
  recob::TrackFitHitInfo const* TrackCollectionProxyElement<CollProxy>::fitInfoAtPoint(
    std::size_t index) const
  {
    if constexpr (base_t::template has<Tracks::TrackFitHitInfoTag>()) {
      auto const& fitInfo = base_t::template get<Tracks::TrackFitHitInfoTag>();
      return &(fitInfo[index]);
    }
    else
      return nullptr;
  } // TrackCollectionProxyElement<>::fitInfoAtPoint()

  //----------------------------------------------------------------------------
  template <typename CollProxy>
  template <typename Pred>
  auto TrackCollectionProxyElement<CollProxy>::selectPoints(Pred&& pred) const
  {
    return util::filterRangeFor(points(), std::forward<Pred>(pred));
  }

  //----------------------------------------------------------------------------
  template <typename CollProxy>
  auto TrackCollectionProxyElement<CollProxy>::pointsWithFlags(
    recob::TrackTrajectory::PointFlags_t::Mask_t mask) const
  {
    return selectPoints([mask](auto&& point) { return point.flags().match(mask); });
  } // TrackCollectionProxyElement<>::pointsWithFlags()

  //----------------------------------------------------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_TRACK_H
