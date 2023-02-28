/**
 * @file   HitUtils.h
 * @brief  Functions and objects interfacing with recob::Hit
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 16th, 2014
 * @see    HitUtils.cxx
 *
 * The utilities hereby provided should supply the functionality that was
 * removed in the simplification of recob::Hit (removal of wire and digit
 * pointers, etc).
 */

#ifndef HITUTILS_H
#define HITUTILS_H

// LArSoft libraries
#include "lardata/ArtDataHelper/FindAllP.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

/// LArSoft-specific namespace
namespace lar {

  /// LArSoft utility namespace
  namespace util {

    /** ************************************************************************
     * @brief Query object connecting a hit to a wire
     *
     * Once upon a time, recob::Hit had a art::Ptr<recob::Wire> in it, and life
     * was easy.
     * When it was discovered that art pointers in data products were evil, they
     * were banned from recob::Hit.
     * As always, evil turns out to be convenient.
     * This query object tries to provide in an efficient way a connection
     * between a hit and the wire that has generated it.
     *
     * This object expects
     *
     * Example of usage: let hit_ptr be a valid art::Ptr<recob::Hit>. Then
     *
     *     HitToWire HtoW(evt);
     *     art::Ptr<recob::Wire> wire_ptr = HtoW[hit_ptr];
     *
     * If the association label is known, it can be used to selectively load
     * that association:
     *
     *     HitToWire HtoW(evt, AssociationInputTag);
     *     art::Ptr<recob::Wire> wire_ptr = HtoW[hit_ptr];
     *
     * that has little advantage (in fact, it is possibly slower) respect to
     * using art::FindOneP.
     */
    using HitToWire = details::FindAllP<recob::Hit, recob::Wire>;

  } // namespace util

} // namespace lar

#endif // HITUTILS_H
