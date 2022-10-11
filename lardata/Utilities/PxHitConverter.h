////////////////////////////////////////////////////////////////////////
// \file PxHitConverter.h
//
// \brief conversion utulities from recob::Hit to PxHit
//
// \author andrzej.szelc@yale.edu, based on LarLite code by Kazu
//
////////////////////////////////////////////////////////////////////////

#ifndef UTIL_PXHITCONVERTER_H
#define UTIL_PXHITCONVERTER_H

#include "PxUtils.h"
#include "lardata/Utilities/Dereference.h"
#include "lardataobj/RecoBase/Hit.h"

#include "canvas/Persistency/Common/Ptr.h"

#include <algorithm>
#include <type_traits>
#include <vector>

///General LArSoft Utilities
namespace util {
  class GeometryUtilities;

  class PxHitConverter {
  public:
    explicit PxHitConverter(GeometryUtilities const& geomUtils);

    /// Generate: from 1 set of hits => 1 set of PxHits using indexes (association)
    void GeneratePxHit(const std::vector<unsigned int>& hit_index,
                       const std::vector<art::Ptr<recob::Hit>> hits,
                       std::vector<PxHit>& pxhits) const;

    /// Generate: from 1 set of hits => 1 set of PxHits using using all hits
    void GeneratePxHit(std::vector<art::Ptr<recob::Hit>> const& hits,
                       std::vector<PxHit>& pxhits) const;

    void GenerateSinglePxHit(art::Ptr<recob::Hit> const& hit, PxHit& pxhits) const;

    /// Generates and returns a PxHit out of a recob::Hit
    PxHit HitToPxHit(recob::Hit const& hit) const;

    /// Generates and returns a PxHit out of a pointer to recob::Hit
    /// or a hit itself
    template <typename HitObj>
    PxHit ToPxHit(HitObj const& hit) const;

    /// Returns a vector of PxHit out of a vector of hits
    template <typename Cont, typename Hit = typename Cont::value_type>
    std::vector<PxHit> ToPxHitVector(Cont const& hits) const;

  private:
    GeometryUtilities const& fGeomUtils;
  }; // class PxHitConverter

} //namespace util

//******************************************************************************
//***  Template implementation
//***
template <typename HitObj>
util::PxHit util::PxHitConverter::ToPxHit(HitObj const& hit) const
{
  // check that the argument is an object convertible to a recob::Hit,
  // or it is a pointer to such an object
  static_assert(
    std::is_convertible<typename lar::util::dereferenced_type<HitObj>::type, recob::Hit>::value,
    "The argument to PxHitConverter::ToPxHit() does not point to a recob::Hit");
  return HitToPxHit(lar::util::dereference(hit));
} // PxHitConverter::ToPxHit()

template <typename Cont, typename Hit /* = typename Cont::value_type */>
std::vector<util::PxHit> util::PxHitConverter::ToPxHitVector(Cont const& hits) const
{
  std::vector<PxHit> pxhits;
  pxhits.reserve(hits.size());
  std::transform(hits.begin(), hits.end(), std::back_inserter(pxhits), [this](Hit const& hit) {
    return this->ToPxHit(hit);
  });
  return pxhits;
} // util::PxHitConverter::ToPxHitVector()

#endif // UTIL_PXHITCONVERTER_H
