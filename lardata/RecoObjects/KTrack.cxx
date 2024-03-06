///////////////////////////////////////////////////////////////////////
///
/// \file   KTrack.cxx
///
/// \brief  Basic Kalman filter track class, without error.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KTrack.h"
#include "cetlib_except/exception.h"
#include <cmath>
#include <cstdlib>
#include <ostream>

// Define some particle masses here.

namespace {
  const double mumass = 0.105658367; // Muon
  const double pimass = 0.13957;     // Charged pion
  const double kmass = 0.493677;     // Charged kaon
  const double pmass = 0.938272;     // Proton
}

namespace trkf {

  /// Default constructor.
  KTrack::KTrack() : fDir(Surface::UNKNOWN), fPdgCode(0) {}

  /// Constructor - specify surface only.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  ///
  KTrack::KTrack(const std::shared_ptr<const Surface>& psurf)
    : fSurf(psurf), fDir(Surface::UNKNOWN), fPdgCode(0)
  {}

  /// Constructor - surface + track parameters.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  /// vec   - Track state vector.
  /// dir   - Track direction.
  /// pdg   - Pdg code.
  ///
  KTrack::KTrack(std::shared_ptr<const Surface> psurf,
                 const TrackVector& vec,
                 Surface::TrackDirection dir,
                 int pdg)
    : fSurf(psurf), fVec(vec), fDir(dir), fPdgCode(pdg)
  {}

  /// Destructor.
  KTrack::~KTrack() {}

  /// Track direction accessor.
  /// Track direction implied by track parameters has precedence
  /// over track direction attribute.
  /// If the surface pointer is null, return UNKNOWN.
  Surface::TrackDirection KTrack::getDirection() const
  {
    Surface::TrackDirection result = Surface::UNKNOWN;
    if (fSurf.get() != 0) result = fSurf->getDirection(fVec, fDir);
    return result;
  }

  /// Test if track is valid.
  ///
  /// A default-constructed or partially-constructed track, is
  /// invalid by virtue of having an unknown propagation direction
  /// or a null surface pointer.
  ///
  /// Tracks can become invaliddynamically for other reasons.  This
  /// method also does the following checks:
  /// a) Check for invalid floating point values (inf and nan).
  /// b) Surface-dependent checks via virtual method Surface::isTrackValid.
  bool KTrack::isValid() const
  {
    bool result = true;

    // Check for valid direction.

    if (getDirection() == Surface::UNKNOWN) result = false;

    // Check for non-null surface pointer (for safety, should be redundant
    // with previous check.

    if (result && fSurf.get() == 0) result = false;

    // Check for track parameters containing invalid floating point
    // values.

    if (result) {
      for (unsigned int i = 0; i < fVec.size(); ++i) {
        if (!std::isfinite(fVec(i))) {
          result = false;
          break;
        }
      }
    }

    // Surface-dependent check on track validity.

    if (result && !fSurf->isTrackValid(fVec)) result = false;

    // Done.

    return result;
  }

  /// Particle mass based on pdg code.
  double KTrack::Mass() const
  {
    double mass = 0.;
    int apdg = std::abs(fPdgCode);

    // Muon

    if (apdg == 13) mass = mumass;

    // Charged pion

    else if (apdg == 211)
      mass = pimass;

    // Charged kaon

    else if (apdg == 321)
      mass = kmass;

    // (Anti)proton

    else if (apdg == 2212)
      mass = pmass;

    // Anything else throw exception

    else
      throw cet::exception("KTrack") << "Mass requested for invalid pdg id = " << fPdgCode << "\n";

    // Done.

    return mass;
  }

  /// Get position of track.
  /// Throw an exception if track is not valid.
  ///
  /// Arguments:
  ///
  /// xyz - Position vector.
  ///
  void KTrack::getPosition(double xyz[3]) const
  {
    if (!isValid()) throw cet::exception("KTrack") << "Position requested for invalid track.\n";
    fSurf->getPosition(fVec, xyz);
  }

  /// Get x-latitude.
  ///
  /// The x-latitude is the latitude defined with respect to the
  /// x-axis.  The x-latitude is zero of the track is traveling
  /// parallel to the wire planes.
  ///
  double KTrack::XLatitude() const
  {
    double mom[3];
    getMomentum(mom);
    double ptx = std::sqrt(mom[1] * mom[1] + mom[2] * mom[2]);
    double result = 0.;
    if (ptx > 0. || mom[0] > 0.) result = atan2(mom[0], ptx);
    return result;
  }

  /// Get x-longitude.
  ///
  /// The x-longitude is the longitude defined with respect to the y-
  /// and z-axes.  The x-longitude is zero of the track is parallel to
  /// the z-axis in the yz-plane.
  ///
  double KTrack::XLongitude() const
  {
    double mom[3];
    getMomentum(mom);
    double result = 0.;
    if (mom[1] != 0. || mom[2] != 0.) result = atan2(mom[1], mom[2]);
    return result;
  }

  /// Get momentum vector of track.
  /// Throw an exception if track is not valid.
  ///
  /// Arguments:
  ///
  /// mom - Momentum vector of track.
  ///
  void KTrack::getMomentum(double mom[3]) const
  {
    if (!isValid())
      throw cet::exception("KTrack") << "Momentum vector requested for invalid track.\n";
    Surface::TrackDirection dir = getDirection();
    fSurf->getMomentum(fVec, mom, dir);
  }

  /// Printout
  std::ostream& KTrack::Print(std::ostream& out, bool doTitle) const
  {
    if (doTitle) out << "KTrack:\n";
    double xyz[3];
    double dir[3];
    getPosition(xyz);
    getMomentum(dir);
    double p = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    if (p != 0.) {
      dir[0] /= p;
      dir[1] /= p;
      dir[2] /= p;
    }
    out << "  Surface direction = "
        << (fDir == Surface::FORWARD ? "FORWARD" :
                                       (fDir == Surface::BACKWARD ? "BACKWARD" : "UNKNOWN"))
        << "\n"
        << "  Pdg = " << fPdgCode << "\n"
        << "  Surface: " << *fSurf << "\n"
        << "  Track parameters:\n"
        << "  [";
    for (unsigned int i = 0; i < fVec.size(); ++i) {
      if (i != 0) out << ", ";
      out << fVec(i);
    }
    out << "]\n";
    out << "  Position:  [" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "]\n";
    out << "  Direction: [" << dir[0] << ", " << dir[1] << ", " << dir[2] << "]\n";
    out << "  X-Latitude  = " << XLatitude() << "\n";
    out << "  X-Longitude = " << XLongitude() << "\n";
    return out;
  }

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const KTrack& trk)
  {
    return trk.Print(out);
  }

} // end namespace trkf
