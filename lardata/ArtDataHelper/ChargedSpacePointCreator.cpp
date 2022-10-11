/**
 * @file   lardata/ArtDataHelper/ChargedSpacePointCreator.cpp
 * @brief  Helpers to create space points with associated charge.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 21, 2017
 * @see    lardata/ArtDataHelper/ChargedSpacePointCreator.h
 *
 */

// class header
#include "lardata/ArtDataHelper/ChargedSpacePointCreator.h"

// framework libraries
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <cassert>
#include <iterator> // std::back_inserter()
#include <utility>  // std::move()

//------------------------------------------------------------------------------
recob::ChargedSpacePointCollectionCreator::ChargedSpacePointCollectionCreator(
  art::Event& event,
  std::string const& instanceName /* = {} */)
  : fEvent(event)
  , fInstanceName(instanceName)
  , fSpacePoints(std::make_unique<std::vector<recob::SpacePoint>>())
  , fCharges(std::make_unique<std::vector<recob::PointCharge>>())
{}

//------------------------------------------------------------------------------
recob::ChargedSpacePointCollectionCreator recob::ChargedSpacePointCollectionCreator::forPtrs(
  art::Event& event,
  std::string const& instanceName /* = {} */)
{
  ChargedSpacePointCollectionCreator creator(event, instanceName);
  creator.fSpacePointPtrMaker =
    std::make_unique<art::PtrMaker<recob::SpacePoint>>(event, instanceName);
  creator.fChargePtrMaker =
    std::make_unique<art::PtrMaker<recob::PointCharge>>(event, instanceName);
  return creator;
} // ChargedSpacePointCollectionCreator(ProducesCollector)

//------------------------------------------------------------------------------
void recob::ChargedSpacePointCollectionCreator::add(recob::SpacePoint const& spacePoint,
                                                    recob::PointCharge const& charge)
{
  // if these assertion fail, add() is being called after put()
  assert(fSpacePoints);
  assert(fCharges);

  fSpacePoints->push_back(spacePoint);
  fCharges->push_back(charge);

  assert(fSpacePoints->size() == fCharges->size());

} // recob::ChargedSpacePointCollectionCreator::add(copy)

//------------------------------------------------------------------------------
void recob::ChargedSpacePointCollectionCreator::add(recob::SpacePoint&& spacePoint,
                                                    recob::PointCharge&& charge)
{
  // if these assertion fail, add() is being called after put()
  assert(fSpacePoints);
  assert(fCharges);

  fSpacePoints->push_back(std::move(spacePoint));
  fCharges->push_back(std::move(charge));

  assert(fSpacePoints->size() == fCharges->size());

} // recob::ChargedSpacePointCollectionCreator::add()

//------------------------------------------------------------------------------
void recob::ChargedSpacePointCollectionCreator::addAll(std::vector<recob::SpacePoint>&& spacePoints,
                                                       std::vector<recob::PointCharge>&& charges)
{
  // if these assertion fail, addAll() is being called after put()
  assert(fSpacePoints);
  assert(fCharges);

  if (spacePoints.size() != charges.size()) {
    throw cet::exception("ChargedSpacePointCollectionCreator")
      << "Input collections of inconsistent size:"
      << " " << spacePoints.size() << " (space points)"
      << "and " << charges.size() << " (charges)"
      << "\n";
  }
  if (empty()) {
    *fSpacePoints = std::move(spacePoints);
    *fCharges = std::move(charges);
  }
  else {
    fSpacePoints->reserve(fSpacePoints->size() + spacePoints.size());
    for (auto&& obj : spacePoints)
      fSpacePoints->push_back(std::move(obj));
    spacePoints.clear();

    fCharges->reserve(fCharges->size() + charges.size());
    for (auto&& obj : charges)
      fCharges->push_back(std::move(obj));
    charges.clear();
  }

  assert(fSpacePoints->size() == fCharges->size());

} // recob::ChargedSpacePointCollectionCreator::addAll()

//------------------------------------------------------------------------------
void recob::ChargedSpacePointCollectionCreator::addAll(
  std::vector<recob::SpacePoint> const& spacePoints,
  std::vector<recob::PointCharge> const& charges)
{
  // if these assertion fail, addAll() is being called after put()
  assert(fSpacePoints);
  assert(fCharges);

  if (spacePoints.size() != charges.size()) {
    throw cet::exception("ChargedSpacePointCollectionCreator")
      << "Input collections of inconsistent size:"
      << " " << spacePoints.size() << " (space points)"
      << "and " << charges.size() << " (charges)"
      << "\n";
  }
  fSpacePoints->reserve(fSpacePoints->size() + spacePoints.size());
  std::copy(spacePoints.begin(), spacePoints.end(), std::back_inserter(*fSpacePoints));

  fCharges->reserve(fCharges->size() + charges.size());
  std::copy(charges.begin(), charges.end(), std::back_inserter(*fCharges));

  assert(fSpacePoints->size() == fCharges->size());

} // recob::ChargedSpacePointCollectionCreator::addAll()

//------------------------------------------------------------------------------
/// Puts all data products into the event, leaving the creator `empty()`.
void recob::ChargedSpacePointCollectionCreator::put()
{

  fEvent.put(std::move(fSpacePoints), fInstanceName);
  fEvent.put(std::move(fCharges), fInstanceName);

  assert(spent());
  assert(empty());
} // recob::ChargedSpacePointCollectionCreator::put()

//------------------------------------------------------------------------------
void recob::ChargedSpacePointCollectionCreator::clear()
{

  if (fSpacePoints) fSpacePoints->clear();
  if (fCharges) fCharges->clear();

  assert(empty());

} // recob::ChargedSpacePointCollectionCreator::clear()

//------------------------------------------------------------------------------
art::Ptr<recob::SpacePoint> recob::ChargedSpacePointCollectionCreator::spacePointPtr(
  std::size_t i) const
{
  return fSpacePointPtrMaker ? (*fSpacePointPtrMaker)(i) : art::Ptr<recob::SpacePoint>{};
} // recob::ChargedSpacePointCollectionCreator::spacePointPtr()

//------------------------------------------------------------------------------
art::Ptr<recob::PointCharge> recob::ChargedSpacePointCollectionCreator::chargePtr(
  std::size_t i) const
{
  return fChargePtrMaker ? (*fChargePtrMaker)(i) : art::Ptr<recob::PointCharge>{};
} // recob::ChargedSpacePointCollectionCreator::chargePtr()

//------------------------------------------------------------------------------
void recob::ChargedSpacePointCollectionCreator::produces(art::ProducesCollector& producesCollector,
                                                         std::string const& instanceName)
{

  producesCollector.produces<std::vector<recob::SpacePoint>>(instanceName);
  producesCollector.produces<std::vector<recob::PointCharge>>(instanceName);

} // recob::ChargedSpacePointCollectionCreator::produces()

//------------------------------------------------------------------------------
