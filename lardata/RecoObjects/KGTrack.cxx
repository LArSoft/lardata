///////////////////////////////////////////////////////////////////////
///
/// \file   KGTrack.cxx
///
/// \brief  A collection of KHitTracks.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iomanip>

#include "lardata/RecoObjects/KGTrack.h"
#include "lardata/RecoObjects/KHitWireLine.h"
#include "lardata/RecoObjects/KHitWireX.h"
#include "lardata/RecoObjects/PropXYZPlane.h"
#include "lardata/RecoObjects/SurfXYZPlane.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RecoBase/TrajectoryPointFlags.h"

#include "cetlib_except/exception.h"

namespace trkf {

  KGTrack::KGTrack(int prefplane) : fPrefPlane(prefplane) {}

  /// Track at start point.
  const KHitTrack& KGTrack::startTrack() const
  {
    if (!isValid()) throw cet::exception("KGTrack") << "Starting track is not valid.\n";

    // Return track.

    return (*fTrackMap.begin()).second;
  }

  /// Track at end point.
  const KHitTrack& KGTrack::endTrack() const
  {
    /// Throw exception if track is not valid.

    if (!isValid()) throw cet::exception("KGTrack") << "Ending track is not valid.\n";

    // Return track.

    return (*fTrackMap.rbegin()).second;
  }

  /// Modifiable track at start point.
  KHitTrack& KGTrack::startTrack()
  {
    /// Throw exception if track is not valid.

    if (!isValid()) throw cet::exception("KGTrack") << "Starting track is not valid.\n";

    // Return track.

    return (*fTrackMap.begin()).second;
  }

  /// Modifiable track at end point.
  KHitTrack& KGTrack::endTrack()
  {
    /// Throw exception if track is not valid.

    if (!isValid()) throw cet::exception("KGTrack") << "Ending track is not valid.\n";

    // Return track.

    return (*fTrackMap.rbegin()).second;
  }

  /// Add track.
  void KGTrack::addTrack(const KHitTrack& trh)
  {
    if (!trh.isValid()) throw cet::exception("KGTrack") << "Adding invalid track to KGTrack.\n";
    fTrackMap.insert(std::make_pair(trh.getPath() + trh.getHit()->getPredDistance(), trh));
  }

  /// Recalibrate track map.
  ///
  /// Loop over contents of track map.  Copy each KHitTrack into a new multimap track map.
  /// Offset the distance stored in the KHitTracks such that the minimum distance is zero.
  /// Also update multimap keys to agree with distance stored in track.
  ///
  void KGTrack::recalibrate()
  {
    std::multimap<double, KHitTrack> newmap;

    // Loop over old track map.

    bool first = true;
    double s0 = 0.;
    for (std::multimap<double, KHitTrack>::iterator i = fTrackMap.begin(); i != fTrackMap.end();
         ++i) {
      KHitTrack& trh = (*i).second;
      if (first) {
        first = false;
        s0 = trh.getPath();
      }
      double s = trh.getPath() - s0;
      trh.setPath(s);
      newmap.insert(std::make_pair(s, trh));
    }

    // Update data member track map.

    fTrackMap.swap(newmap);
  }

  /// Fill a recob::Track.
  ///
  /// Arguments:
  ///
  /// track - Track to fill.
  ///
  void KGTrack::fillTrack(detinfo::DetectorPropertiesData const& detProp,
                          recob::Track& track,
                          int id) const
  {

    // Make propagator for propating to standard track surface.

    PropXYZPlane prop(detProp, 0., false);

    // Fill collections of trajectory points and direction vectors.

    std::vector<recob::tracking::Point_t> xyz;
    std::vector<recob::tracking::Vector_t> pxpypz;
    std::vector<recob::tracking::SMatrixSym55> cov;
    std::vector<recob::TrajectoryPointFlags> outFlags;

    xyz.reserve(fTrackMap.size());
    pxpypz.reserve(fTrackMap.size());
    outFlags.reserve(fTrackMap.size());

    // Loop over KHitTracks.

    int ndof = 0;
    float totChi2 = 0.;
    unsigned int n = 0;
    for (std::multimap<double, KHitTrack>::const_iterator itr = fTrackMap.begin();
         itr != fTrackMap.end();
         ++itr, ++n) {
      const KHitTrack& trh = (*itr).second;

      // Get position.

      double pos[3];
      trh.getPosition(pos);
      xyz.push_back({pos[0], pos[1], pos[2]});

      // Get momentum vector.
      // Fill direction unit vector and momentum.

      double mom[3];
      trh.getMomentum(mom);
      double p = std::sqrt(mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2]);
      if (p == 0.) throw cet::exception("KGTrack") << __func__ << ": null momentum\n";
      pxpypz.push_back({mom[0], mom[1], mom[2]});

      ndof += 1;
      totChi2 += trh.getChisq();
      outFlags.emplace_back(n, recob::TrajectoryPointFlags::makeMask());

      // Fill error matrix.

      recob::tracking::SMatrixSym55 covar;

      // Construct surface perpendicular to track momentun, and
      // propagate track to that surface (zero distance).

      const std::shared_ptr<const Surface> psurf(
        new SurfXYZPlane(pos[0], pos[1], pos[2], mom[0], mom[1], mom[2]));
      KETrack tre(trh);
      std::optional<double> dist = prop.err_prop(tre, psurf, Propagator::UNKNOWN, false);
      if (!dist) throw cet::exception("KGTrack") << __func__ << ": error propagation failed\n";
      for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j)
          covar(i, j) = tre.getError()(i, j);
      }

      // Only save first and last error matrix.

      if (cov.size() < 2)
        cov.push_back(covar);
      else
        cov.back() = covar;
    }

    // Fill track.

    ndof = ndof - 4; //fit measures 4 parameters: position and direction on plane
    if (xyz.size() >= 2) {
      track = recob::Track(std::move(xyz),
                           std::move(pxpypz),
                           std::move(outFlags),
                           true,
                           this->startTrack().PdgCode(),
                           totChi2,
                           ndof,
                           std::move(cov.front()),
                           std::move(cov.back()),
                           id);
    }
  }

  /// Fill a PtrVector of Hits.
  ///
  /// Arguments:
  ///
  /// hits - Hit vector to fill.
  ///
  void KGTrack::fillHits(art::PtrVector<recob::Hit>& hits,
                         std::vector<unsigned int>& hittpindex) const
  {
    hits.reserve(hits.size() + fTrackMap.size());

    // Loop over KHitTracks and fill hits belonging to this track.

    unsigned int counter = 0; //Index of corresponding trajectory point
    for (std::multimap<double, KHitTrack>::const_iterator it = fTrackMap.begin();
         it != fTrackMap.end();
         ++it) {
      const KHitTrack& track = (*it).second;
      ++counter;
      // Extrack Hit from track.
      const std::shared_ptr<const KHitBase>& hit = track.getHit();
      if (const KHitWireX* phit = dynamic_cast<const KHitWireX*>(&*hit)) {
        const art::Ptr<recob::Hit> prhit = phit->getHit();
        if (!prhit.isNull()) {
          hits.push_back(prhit);
          hittpindex.push_back(counter - 1);
        }
      }
      else if (const KHitWireLine* phit = dynamic_cast<const KHitWireLine*>(&*hit)) {
        const art::Ptr<recob::Hit> prhit = phit->getHit();
        if (!prhit.isNull()) {
          hits.push_back(prhit);
          hittpindex.push_back(counter - 1);
        }
      }
    }
  }

  ///
  /// Printout
  ///
  std::ostream& KGTrack::Print(std::ostream& out) const
  {

    int n = 0;

    double oldxyz[3] = {0., 0., 0.};
    double len = 0.;
    bool first = true;
    for (auto const& ele : fTrackMap) {
      double s = ele.first;
      const KHitTrack& trh = ele.second;
      double xyz[3];
      double mom[3];
      trh.getPosition(xyz);
      trh.getMomentum(mom);
      double tmom = std::sqrt(mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2]);
      if (tmom != 0.) {
        mom[0] /= tmom;
        mom[1] /= tmom;
        mom[2] /= tmom;
      }
      if (!first) {
        double dx = xyz[0] - oldxyz[0];
        double dy = xyz[1] - oldxyz[1];
        double dz = xyz[2] - oldxyz[2];
        len += std::sqrt(dx * dx + dy * dy + dz * dz);
      }
      const KHitBase& hit = *(trh.getHit());
      int plane = hit.getMeasPlane();
      std::ios_base::fmtflags f = out.flags();
      out << "State " << std::setw(4) << n << ", path=" << std::setw(8) << std::fixed
          << std::setprecision(2) << s << ", length=" << std::setw(8) << len
          << ", x=" << std::setw(8) << xyz[0] << ", y=" << std::setw(8) << xyz[1]
          << ", z=" << std::setw(8) << xyz[2] << ", dx=" << std::setw(8) << mom[0]
          << ", dy=" << std::setw(8) << mom[1] << ", dz=" << std::setw(8) << mom[2]
          << ", plane=" << std::setw(1) << plane << "\n";
      out.flags(f);

      oldxyz[0] = xyz[0];
      oldxyz[1] = xyz[1];
      oldxyz[2] = xyz[2];

      ++n;
      first = false;
    }
    return out;
  }

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const KGTrack& trg) { return trg.Print(out); }

} // end namespace trkf
