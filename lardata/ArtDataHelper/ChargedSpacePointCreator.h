/**
 * @file   lardata/ArtDataHelper/ChargedSpacePointCreator.h
 * @brief  Helpers to create space points with associated charge.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 21, 2017
 * @see    lardata/ArtDataHelper/ChargedSpacePointCreator.cpp
 *
 * The unit test for this utility is part of the `proxy::ChargedSpacePoints`
 * unit test.
 *
 */

#ifndef LARDATA_ARTDATAHELPER_CHARGEDSPACEPOINTCREATOR_H
#define LARDATA_ARTDATAHELPER_CHARGEDSPACEPOINTCREATOR_H

// LArSoft libraries
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// framework libraries
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Ptr.h"

// C/C++ standard libraries
#include <cstdlib>     // std::size_t
#include <memory>      // std::unique_ptr<>
#include <type_traits> // std::enable_if_t, ...
#include <vector>

namespace art {
  class Event;
  class ProducesCollector;
}

namespace recob {

  /**
   * @brief Creates a collection of space points with associated charge.
   *
   * This class facilitates the creation of data products satisfying the
   * requirements the `proxy::ChargedSpacePoints` proxy relies on.
   * It will keep track of space points and reconstructed charge, and will put
   * them into the event at the end.
   *
   * Requirements of the data products
   * ==================================
   *
   * The requirements guaranteed by the output of this collection creator
   * satisfy the `proxy::ChargedSpacePoints` proxy requirements.
   * They are:
   * * space points and charges stored in two separate data products
   * * space points and charge in the same order, so that charge at position
   *   `i` is associated to space point at the same position `i`
   * * one-to-one correspondence between each space point and its charge
   * * association is implicit in the requirements above: no `art::Assns` data
   *   product is produced
   *
   *
   * Usage
   * ======
   *
   * The usage pattern is made of two main parts:
   * # declaration of the data products, at producer construction time
   * # production of the data products, event by event
   *
   * The second part happens within the context of the producer's `produce()`
   * (or `filter()`, or equivalent) method, and it can be split into three
   * stages:
   * # construction of the collection creator, binding it to the current event
   * # filling of the creator, usually in a loop
   * # explicit transfer of the data products into the event
   *
   *
   * Declaration of the data products
   * ---------------------------------
   *
   * In the same fashion as data products must be declared to _art_ with a
   * `produces()` call, the collection creator will have to perform an
   * `equivalent step. This is achieved by calling the static `produces()`
   * method from your module's constructor (see its documentation for an example).
   *
   *
   * Construction of a collection creator object
   * --------------------------------------------
   *
   * Collection creator objects are bound to a specific event and therefore
   * can't be data members of the producer class. In the `produces()` method,
   * we'll have:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void MyProducer::produce(art::Event& event) {
   *
   *   recob::ChargedSpacePointCollectionCreator spacePoints(event);
   *
   *   // ...
   *
   * } // MyProducer::produce()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * If _art_ pointers to the data products are needed (e.g. to create
   * associations), then the named-constructor idiom should be followed:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void MyProducer::produce(art::Event& event) {
   *
   *   auto spacePoints = recob::ChargedSpacePointCollectionCreator::forPtrs(event);
   *
   *   // ...
   *
   * } // MyProducer::produce()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * In both cases, an instance name can be specified which will be used for all
   * the managed data products (space points and reconstructed charge).
   *
   *
   * Populating the collections
   * ---------------------------
   *
   * The core of the usage of this object is feeding the objects to the data
   * product collections. This is done using the `add()` member function.
   * If the data objects already exist, they can be moved in instead of being
   * copied. For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void MyProducer::produce(art::Event& event) {
   *
   *   recob::ChargedSpacePointCollectionCreator spacePoints(event);
   *
   *   auto hitAssns
   *     = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit>>();
   *
   *   // some processing
   *
   *   for (auto const& hitSet: hitSets) {
   *
   *     recob::SpacePoint spacePoint;
   *     recob::PointCharge charge;
   *
   *     fSpacePointChargeAlgo->run(hitSet, spacePoint, charge);
   *
   *     // add the new space point and charge to the collection
   *     spacePoints.add(std::move(spacePoint), std::move(charge));
   *
   *     // associate this space point with all the hits in the source set
   *     auto const& spacePointPtr = spacePoints.lastSpacePointPtr();
   *     for (art::Ptr<recob::Hit> const& hitPtr: hitSet)
   *       hitAssns.addSingle(spacePointPtr, hitPtr);
   *
   *   } // for
   *
   *   // ...
   *
   * } // MyProducer::produce()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * If your algorithm is creating a subcollection of space points and charges
   * which are in the same order, a shortcut to a loop of `add()` is `addAll()`:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void MyProducer::produce(art::Event& event) {
   *
   *   recob::ChargedSpacePointCollectionCreator spacePoints(event);
   *
   *   // some processing
   *
   *   for (auto const& track: tracks) {
   *
   *     std::vector<recob::SpacePoint> trackSpacePoints;
   *     std::vector<recob::PointCharge> trackCharges;
   *
   *     fSpacePointChargeAlgo->run(track, trackSpacePoints, trackCharges);
   *
   *     spacePoints.addAll
   *       (std::move(trackSpacePoints), std::move(trackCharges));
   *
   *   } // for
   *
   *   // ...
   *
   * } // MyProducer::produce()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   *
   *
   * Operations on the collection
   * -----------------------------
   *
   * While the collection creator object is designed to be just a single-use,
   * single-pattern helper, there are a few operations that it allows:
   *
   * * query: the current number of space points (`size()`), whether there are
   *          any (`empty()`), and whether `put()` has been already called
   *          (`spent()`; see below)
   * * deletion: it is possible to remove *all* the data so far collected;
   *             this may be useful in case of error, where the data product
   *             should be set in as default-constructed (which unfortunately
   *             in the case of space points _might_ be, albeit unlikely, a
   *             legal outcome)
   * * if art::Ptr-creation has been enabled (by calling the static 'forPtrs(...)' function), _art_
   *   pointers can be created with `lastSpacePointPtr()`, `lastChargePtr()`, `spacePointPtr()`
   *   and `chargePtr()` to the elements of the future data products
   *
   *
   * Insertion of the data products into the event
   * ----------------------------------------------
   *
   * Again as in the standard producer pattern, where `Event::put()` needs to
   * be called to finally move the data into _art_, the collection creator will
   * need to be told to do the same:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * void MyProducer::produce(art::Event& event) {
   *
   *   recob::ChargedSpacePointCollectionCreator spacePoints(event);
   *
   *   auto hitAssns
   *     = std::make_unique<art::Assns<recob::SpacePoint, recob::Hit>>();
   *
   *   // ... filling here ...
   *
   *   spacePoints.put();
   *   event.put(std::move(hitAssns)); // other data products go the usual way
   *
   * } // MyProducer::produce()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The instance name used in `put()` had been set at construction time.
   *
   * *After `put()` is called*, the object has served its purpose and *can't be
   * used any further*. In this state, `spent()` method will return `true`.
   *
   */
  class ChargedSpacePointCollectionCreator {

  public:
    //--- BEGIN Constructors ---------------------------------------------------
    /// @{
    /// @name Constructors

    /**
     * @brief Constructor binding this object to a specific _art_ event.
     * @param event the _art_ event to bind to
     * @param instanceName _(default: empty)_ instance name for all data
     *                     products
     *
     * When the object is constructed with this constructor, the creation of
     * _art_ pointers will not be enabled (`canMakePointers()` will return
     * `false`).
     */
    ChargedSpacePointCollectionCreator(art::Event& event, std::string const& instanceName = {});

    /**
     * @brief Static function binding a new object to a specific _art_ event.
     * @param event the _art_ event to bind to
     * @param instanceName _(default: empty)_ instance name for all data
     *                     products
     *
     * This static function follows the named-constructor idiom,
     * enabling the creation of _art_ pointers.
     */
    static ChargedSpacePointCollectionCreator forPtrs(art::Event& event,
                                                      std::string const& instanceName = {});

    /// @}
    //--- END Constructors ---------------------------------------------------

    //--- BEGIN Insertion and finish operations --------------------------------
    /// @{
    /// @name Insertion and finish operations

    //@{
    /**
     * @brief Inserts the specified space point and associated data into the
     *        collection.
     * @param spacePoint the space point to be copied into the collection
     * @param charge the charge to be copied into the collection
     *
     * The data is pushed as the new last element of the collection.
     *
     * Data is copied or moved depending on which variant of this method is
     * used.
     */
    void add(recob::SpacePoint const& spacePoint, recob::PointCharge const& charge);
    void add(recob::SpacePoint&& spacePoint, recob::PointCharge&& charge);
    //@}

    //@{
    /**
     * @brief Inserts all the space points and associated data from the vectors
     *        into the collection.
     * @param spacePoints the space point to be copied into the collection
     * @param charges the charges to be copied into the collection
     * @throw cet::exception (category `ChargedSpacePointCollectionCreator`)
     *                       if the input collections are inconsistent
     *
     * The data is pushed as the new last element of the collection.
     *
     * Data is copied or moved depending on which variant of this method is
     * used.
     *
     * No exception safety is offered here.
     */
    void addAll(std::vector<recob::SpacePoint>&& spacePoints,
                std::vector<recob::PointCharge>&& charges);
    void addAll(std::vector<recob::SpacePoint> const& spacePoints,
                std::vector<recob::PointCharge> const& charges);
    //@}

    /**
     * @brief Puts all data products into the event, leaving the creator
     *        `empty()`.
     *
     * The accumulated data is moved into the event.
     *
     * This is the last valid action of the object.
     * After this, only `empty()`, `spent()` and the _art_ pointer makers
     * (if enabled) are guaranteed to work.
     */
    void put();

    /// @}
    //--- END Insertion and finish operations --------------------------------

    //--- BEGIN Queries and operations -----------------------------------------
    /// @{
    /// @name Queries and operations

    /// Returns whether there are currently no space points in the collection.
    bool empty() const { return spent() || fSpacePoints->empty(); }

    /// Returns the number of space points currently in the collection.
    std::size_t size() const { return spent() ? 0U : fSpacePoints->size(); }

    /// Removes all data from the collection, making it `empty()`.
    void clear();

    /// Returns whether `put()` has already been called.
    bool spent() const { return !fSpacePoints; }

    /// Returns whether _art_ pointer making is enabled.
    bool canMakePointers() const { return bool(fSpacePointPtrMaker); }

    ///@}
    //--- END Queries and operations -------------------------------------------

    //--- BEGIN Complimentary unchecked element access -------------------------
    ///@{
    ///@name Complimentary unchecked element access

    /// Returns the specified space point; undefined behaviour if not there.
    recob::SpacePoint const& spacePoint(std::size_t i) const { return fSpacePoints->operator[](i); }

    /// Returns the last inserted space point; undefined behaviour if `empty()`.
    recob::SpacePoint const& lastSpacePoint() const { return spacePoint(lastIndex()); }

    /// Returns an _art_ pointer to the specified space point (no check done!).
    art::Ptr<recob::SpacePoint> spacePointPtr(std::size_t i) const;

    /// Returns an _art_ pointer to the last inserted space point (no check!).
    art::Ptr<recob::SpacePoint> lastSpacePointPtr() const { return spacePointPtr(lastIndex()); }

    /// Returns the last inserted charge; undefined behaviour if `empty()`.
    recob::PointCharge const& charge(std::size_t i) const { return fCharges->operator[](i); }

    /// Returns the last inserted charge; undefined behaviour if `empty()`.
    recob::PointCharge const& lastCharge() const { return charge(lastIndex()); }

    /// Returns an _art_ pointer to the specified charge (no check done!).
    art::Ptr<recob::PointCharge> chargePtr(std::size_t i) const;

    /// Returns an _art_ pointer to the inserted charge (no check!).
    art::Ptr<recob::PointCharge> lastChargePtr() const { return chargePtr(lastIndex()); }

    /// @}
    //--- END Complimentary unchecked element access ---------------------------

    //--- BEGIN Static constructor interface -----------------------------------
    /// @{
    /** @name Static constructor interface
     *
     * These methods take care of initialization that needs to take place on
     * construction of the module.
     */

    /**
     * @brief Declares the data products being produced.
     * @param producer the module producing the data products
     * @param instanceName _(default: empty)_ name of instance of all data
     *                     products
     *
     * Call this method in the constructor of `producer`, e.g.:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * MyProducer::MyProducer(Parameters const& config) {
     *
     *   recob::ChargedSpacePointCollectionCreator::produces
     *     (producesCollector(), config().instanceName());
     *
     * } // MyProducer::MyProducer()
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    static void produces(art::ProducesCollector& producesCollector,
                         std::string const& instanceName = {});

    /// @}
    //--- END Static constructor interface -------------------------------------

  private:
    art::Event& fEvent; ///< The event this object is bound to.

    std::string fInstanceName; ///< Instance name of all the data products.

    /// Space point data.
    std::unique_ptr<std::vector<recob::SpacePoint>> fSpacePoints;
    /// Space point pointer maker.
    std::unique_ptr<art::PtrMaker<recob::SpacePoint>> fSpacePointPtrMaker;
    /// Charge data.
    std::unique_ptr<std::vector<recob::PointCharge>> fCharges;
    /// Charge pointer maker.
    std::unique_ptr<art::PtrMaker<recob::PointCharge>> fChargePtrMaker;

    /// Returns the index of the last element (undefined if empty).
    std::size_t lastIndex() const { return size() - 1U; }

  }; // class ChargedSpacePointCollectionCreator

} // namespace recob

//------------------------------------------------------------------------------

#endif // LARDATA_ARTDATAHELPER_CHARGEDSPACEPOINTCREATOR_H
