////////////////////////////////////////////////////////////////////////
//  \file GeometryUtilities.cxx
//
//  \brief Functions to calculate distances and angles in 3D and 2D
//
// andrzej.szelc@yale.edu
//
////////////////////////////////////////////////////////////////////////

#include "lardata/Utilities/GeometryUtilities.h"
#include "cetlib/pow.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/Utilities/UtilException.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TLorentzVector.h"

#include <cmath>

namespace util {

  GeometryUtilities::GeometryUtilities(geo::GeometryCore const& geom,
                                       detinfo::DetectorClocksData const& clockData,
                                       detinfo::DetectorPropertiesData const& propData)
    : fGeom{geom}, fClocks{clockData}, fDetProp{propData}
  {
    fNPlanes = fGeom.Nplanes();
    vertangle.resize(fNPlanes);
    for (UInt_t ip = 0; ip < fNPlanes; ip++)
      vertangle[ip] = fGeom.Plane(ip).Wire(0).ThetaZ(false) - TMath::Pi() / 2.; // wire angle

    fWirePitch = fGeom.WirePitch();
    fTimeTick = sampling_rate(fClocks) / 1000.;
    fDriftVelocity = fDetProp.DriftVelocity(fDetProp.Efield(), fDetProp.Temperature());

    fWiretoCm = fWirePitch;
    fTimetoCm = fTimeTick * fDriftVelocity;
    fWireTimetoCmCm = fTimetoCm / fWirePitch;
  }

  //-----------------------------------------------------------------------------
  // omega0 and omega1 (calculated by CPAN in degrees):
  // angle based on distances in wires and time - rescaled to cm.
  // tan(angle)*fMean_wire_pitch/(fTimeTick*fDriftVelocity);
  // as those calculated with Get2Dangle
  // writes phi and theta in degrees.
  /////////////////////////////////////
  Int_t
  GeometryUtilities::Get3DaxisN(Int_t iplane0,
                                Int_t iplane1,
                                Double_t omega0,
                                Double_t omega1,
                                Double_t& phi,
                                Double_t& theta) const
  {
    // y, z, x coordinates
    Double_t ln(0), mn(0), nn(0);
    Double_t phis(0), thetan(0);

    // Pretend collection and induction planes.
    // "Collection" is the plane with the vertical angle equal to zero.
    // If both are non-zero, collection is the one with the negative angle.
    UInt_t Cplane = 0, Iplane = 1;

    // angleC and angleI are the respective angles to vertical in C/I
    // planes and slopeC, slopeI are the tangents of the track.
    Double_t angleC, angleI, slopeC, slopeI, omegaC, omegaI;
    omegaC = kINVALID_DOUBLE;
    omegaI = kINVALID_DOUBLE;

    // Don't know how to reconstruct these yet, so exit with error.
    // In
    if (omega0 == 0 || omega1 == 0) {
      phi = 0;
      theta = -999;
      return -1;
    }

    //////insert check for existence of planes.

    // check if backwards going track
    Double_t alt_backwards = 0;

    if (fabs(omega0) > (TMath::Pi() / 2.0) || fabs(omega1) > (TMath::Pi() / 2.0)) {
      alt_backwards = 1;
    }

    if (vertangle[iplane0] == 0) {
      // first plane is at 0 degrees
      Cplane = iplane0;
      Iplane = iplane1;
      omegaC = omega0;
      omegaI = omega1;
    }
    else if (vertangle[iplane1] == 0) {
      // second plane is at 0 degrees
      Cplane = iplane1;
      Iplane = iplane0;
      omegaC = omega1;
      omegaI = omega0;
    }
    else if (vertangle[iplane0] != 0 && vertangle[iplane1] != 0) {
      // both planes are at non zero degree - find the one with deg<0
      if (vertangle[iplane1] < vertangle[iplane0]) {
        Cplane = iplane1;
        Iplane = iplane0;
        omegaC = omega1;
        omegaI = omega0;
      }
      else if (vertangle[iplane1] > vertangle[iplane0]) {
        Cplane = iplane0;
        Iplane = iplane1;
        omegaC = omega0;
        omegaI = omega1;
      }
      else {
        // throw error - same plane.
        return -1;
      }
    }

    slopeC = tan(omegaC);
    slopeI = tan(omegaI);
    angleC = vertangle[Cplane];
    angleI = vertangle[Iplane];

    // 0 -1 factor depending on if one of the planes is vertical.
    bool nfact = !(vertangle[Cplane]);

    // ln represents y, omega is 2d angle -- in first 2 quadrants y is positive.
    if (omegaC < TMath::Pi() && omegaC > 0)
      ln = 1;
    else
      ln = -1;

    // calculate x and z using y ( ln )
    mn = (ln / (2 * sin(angleI))) * ((cos(angleI) / (slopeC * cos(angleC))) - (1 / slopeI) +
                                     nfact * (cos(angleI) / (cos(angleC) * slopeC) - 1 / slopeI));

    nn = (ln / (2 * cos(angleC))) *
         ((1 / slopeC) + (1 / slopeI) + nfact * ((1 / slopeC) - (1 / slopeI)));

    // Direction angles
    if (fabs(omegaC) > 0.01) // catch numeric error values
    {
      // phi=atan(ln/nn);
      phis = asin(ln / TMath::Sqrt(ln * ln + nn * nn));

      if (fabs(slopeC + slopeI) < 0.001)
        phis = 0;
      else if (fabs(omegaI) > 0.01 && (omegaI / fabs(omegaI) == -omegaC / fabs(omegaC)) &&
               (fabs(omegaC) < 1 * TMath::Pi() / 180 ||
                fabs(omegaC) > 179 * TMath::Pi() / 180))           // angles have
        phis = (fabs(omegaC) > TMath::Pi() / 2) ? TMath::Pi() : 0; // angles are

      if (nn < 0 && phis > 0 &&
          !(!alt_backwards &&
            fabs(phis) <
              TMath::Pi() / 4)) // do not go back if track looks forward and phi is forward
        phis = (TMath::Pi()) - phis;
      else if (nn < 0 && phis < 0 && !(!alt_backwards && fabs(phis) < TMath::Pi() / 4))
        phis = (-TMath::Pi()) - phis;

      phi = phis * 180 / TMath::Pi();
    }
    // If plane2 (collection), phi = 2d angle (omegaC in this case)
    else {
      phis = omegaC;
      phi = omegaC;
    }

    thetan = -asin(mn / std::hypot(ln, mn, nn));
    theta = thetan * 180 / TMath::Pi();

    return 0;
  }

  //////////////////////////////
  // Calculate theta in case phi~0
  // returns theta in angles
  ////////////////////////////////
  Double_t
  GeometryUtilities::Get3DSpecialCaseTheta(Int_t iplane0,
                                           Int_t iplane1,
                                           Double_t dw0,
                                           Double_t dw1) const
  {

    Double_t splane, lplane; // plane in which the track is shorter and longer.
    Double_t sdw, ldw;

    if (dw0 == 0 && dw1 == 0) return -1;

    if (dw0 >= dw1) {
      lplane = iplane0;
      splane = iplane1;
      ldw = dw0;
      sdw = dw1;
    }
    else {
      lplane = iplane1;
      splane = iplane0;
      ldw = dw1;
      sdw = dw0;
    }

    Double_t top = (std::cos(vertangle[splane]) - std::cos(vertangle[lplane]) * sdw / ldw);
    Double_t bottom = tan(vertangle[lplane] * std::cos(vertangle[splane]));
    bottom -= tan(vertangle[splane] * std::cos(vertangle[lplane])) * sdw / ldw;

    Double_t tantheta = top / bottom;

    return atan(tantheta) * vertangle[lplane] / std::abs(vertangle[lplane]) * 180. / TMath::Pi();
  }

  /////////////////////////////////////////////////////////
  // Calculate 3D pitch in beam coordinates
  //
  /////////////////////////////////////////////////////////
  Double_t
  GeometryUtilities::CalculatePitch(UInt_t iplane, Double_t phi, Double_t theta) const
  {

    Double_t pitch = -1.;

    if (fGeom.Plane(iplane).View() == geo::kUnknown || fGeom.Plane(iplane).View() == geo::k3D) {
      mf::LogError(Form("Warning :  no Pitch foreseen for view %d", fGeom.Plane(iplane).View()));
      return pitch;
    }
    else {

      Double_t pi = TMath::Pi();
      Double_t fTheta = pi / 2 - theta;
      Double_t fPhi = -(phi + pi / 2);

      for (UInt_t i = 0; i < fGeom.Nplanes(); ++i) {
        if (i == iplane) {
          Double_t wirePitch = fGeom.WirePitch(i);
          Double_t angleToVert =
            0.5 * TMath::Pi() - fGeom.WireAngleToVertical(fGeom.Plane(i).View());
          Double_t cosgamma =
            TMath::Abs(TMath::Sin(angleToVert) * TMath::Cos(fTheta) +
                       TMath::Cos(angleToVert) * TMath::Sin(fTheta) * TMath::Sin(fPhi));

          if (cosgamma > 0) pitch = wirePitch / cosgamma;
        } // end if the correct view
      }   // end loop over planes
    }     // end if a reasonable view

    return pitch;
  }

  /////////////////////////////////////////////////////////
  // Calculate 3D pitch in polar coordinates
  //
  /////////////////////////////////////////////////////////
  Double_t
  GeometryUtilities::CalculatePitchPolar(UInt_t iplane, Double_t phi, Double_t theta) const
  {

    Double_t pitch = -1.;

    if (fGeom.Plane(iplane).View() == geo::kUnknown || fGeom.Plane(iplane).View() == geo::k3D) {
      mf::LogError(Form("Warning :  no Pitch foreseen for view %d", fGeom.Plane(iplane).View()));
      return pitch;
    }
    else {

      Double_t fTheta = theta;
      Double_t fPhi = phi;

      for (UInt_t i = 0; i < fGeom.Nplanes(); ++i) {
        if (i == iplane) {
          Double_t wirePitch = fGeom.WirePitch(i);
          Double_t angleToVert =
            0.5 * TMath::Pi() - fGeom.WireAngleToVertical(fGeom.Plane(i).View());
          Double_t cosgamma =
            TMath::Abs(TMath::Sin(angleToVert) * TMath::Cos(fTheta) +
                       TMath::Cos(angleToVert) * TMath::Sin(fTheta) * TMath::Sin(fPhi));

          if (cosgamma > 0) pitch = wirePitch / cosgamma;
        } // end if the correct view
      }   // end loop over planes
    }     // end if a reasonable view

    return pitch;
  }

  /////////////////////////////////////////////////////////
  // Calculate 2D slope
  // in "cm" "cm" coordinates
  /////////////////////////////////////////////////////////
  Double_t
  GeometryUtilities::Get2Dslope(Double_t wireend,
                                Double_t wirestart,
                                Double_t timeend,
                                Double_t timestart) const
  {

    return GeometryUtilities::Get2Dslope((wireend - wirestart) * fWiretoCm,
                                         (timeend - timestart) * fTimetoCm);
  }

  /////////////////////////////////////////////////////////
  // Calculate 2D slope
  // in "cm" "cm" coordinates
  /////////////////////////////////////////////////////////
  double
  GeometryUtilities::Get2Dslope(const util::PxPoint* endpoint,
                                const util::PxPoint* startpoint) const
  {
    return Get2Dslope(endpoint->w - startpoint->w, endpoint->t - startpoint->t);
  }

  /////////////////////////////////////////////////////////
  // Calculate 2D slope
  // in wire time coordinates coordinates
  //
  /////////////////////////////////////////////////////////
  Double_t
  GeometryUtilities::Get2Dslope(Double_t dwire, Double_t dtime) const
  {
    return tan(Get2Dangle(dwire, dtime)) / fWireTimetoCmCm;
  }

  /////////////////////////////////////////////////////////
  // Calculate 2D angle
  // in "cm" "cm" coordinates
  /////////////////////////////////////////////////////////
  Double_t
  GeometryUtilities::Get2Dangle(Double_t wireend,
                                Double_t wirestart,
                                Double_t timeend,
                                Double_t timestart) const
  {

    return Get2Dangle((wireend - wirestart) * fWiretoCm, (timeend - timestart) * fTimetoCm);
  }

  /////////////////////////////////////////////////////////
  // Calculate 2D angle
  // in "cm" "cm" coordinates, endpoint and startpoint are assumed to be in
  // cm/cm space
  /////////////////////////////////////////////////////////
  Double_t
  GeometryUtilities::Get2Dangle(const util::PxPoint* endpoint,
                                const util::PxPoint* startpoint) const
  {
    return Get2Dangle(endpoint->w - startpoint->w, endpoint->t - startpoint->t);
  }

  ////////////////////////////
  // Calculate 2D angle
  // in "cm" "cm" coordinates
  ////////////////////////////
  Double_t
  GeometryUtilities::Get2Dangle(Double_t dwire, Double_t dtime) const
  {
    Double_t BC, AC;
    Double_t omega;

    BC = ((Double_t)dwire); // in cm
    AC = ((Double_t)dtime); // in cm
    omega = std::asin(AC / std::hypot(AC, BC));
    if (BC < 0) // for the time being. Will check if it works for AC<0
    {
      if (AC > 0) {
        omega = TMath::Pi() - std::abs(omega); //
      }
      else if (AC < 0) {
        omega = -TMath::Pi() + std::abs(omega);
      }
      else {
        omega = TMath::Pi();
      }
    }
    return omega;
  }

  // accepting phi and theta in degrees
  // returning in radians.

  double
  GeometryUtilities::Get2DangleFrom3D(unsigned int plane, double phi, double theta) const
  {
    TVector3 dummyvector(cos(theta * TMath::Pi() / 180.) * sin(phi * TMath::Pi() / 180.),
                         sin(theta * TMath::Pi() / 180.),
                         cos(theta * TMath::Pi() / 180.) * cos(phi * TMath::Pi() / 180.));

    return Get2DangleFrom3D(plane, dummyvector);
  }

  // accepting TVector3
  // returning in radians as is customary,

  double
  GeometryUtilities::Get2DangleFrom3D(unsigned int plane, TVector3 dir_vector) const
  {
    double alpha = 0.5 * TMath::Pi() - fGeom.WireAngleToVertical(fGeom.Plane(plane).View());
    // create dummy  xyz point in middle of detector and another one in unit
    // length. calculate correspoding points in wire-time space and use the
    // differnces between those to return 2D a angle

    TVector3 start(fGeom.DetHalfWidth(), 0., fGeom.DetLength() / 2.);
    TVector3 end = start + dir_vector;

    // the wire coordinate is already in cm. The time needs to be converted.
    util::PxPoint startp(
      plane,
      (fGeom.DetHalfHeight() * sin(fabs(alpha)) + start[2] * cos(alpha) - start[1] * sin(alpha)),
      start[0]);

    util::PxPoint endp(
      plane,
      (fGeom.DetHalfHeight() * sin(fabs(alpha)) + end[2] * cos(alpha) - end[1] * sin(alpha)),
      end[0]);

    double angle = Get2Dangle(&endp, &startp);

    return angle;
  }

  //////////////////////////////////////
  // Calculate 2D distance
  // in "cm" "cm" coordinates
  ////////////////////////////////////////
  Double_t
  GeometryUtilities::Get2DDistance(Double_t wire1,
                                   Double_t time1,
                                   Double_t wire2,
                                   Double_t time2) const
  {

    return std::hypot((wire1 - wire2) * fWiretoCm, (time1 - time2) * fTimetoCm);
  }

  double
  GeometryUtilities::Get2DDistance(const util::PxPoint* point1, const util::PxPoint* point2) const
  {
    return std::hypot(point1->w - point2->w, point1->t - point2->t);
  }

  ////////////////////////////
  // Calculate 2D distance, using 2D angle
  // in "cm" "cm" coordinates
  ////////////////////////////
  Double_t
  GeometryUtilities::Get2DPitchDistance(Double_t angle, Double_t inwire, Double_t wire) const
  {
    Double_t radangle = TMath::Pi() * angle / 180;
    if (std::cos(radangle) == 0) return 9999;
    return std::abs((wire - inwire) / std::cos(radangle)) * fWiretoCm;
  }

  //----------------------------------------------------------------------------
  Double_t
  GeometryUtilities::Get2DPitchDistanceWSlope(Double_t slope, Double_t inwire, Double_t wire) const
  {

    return std::abs(wire - inwire) * TMath::Sqrt(1 + slope * slope) * fWiretoCm;
  }

  ///////////////////////////////////
  // Calculate wire,time coordinates of the Hit projection onto a line
  //
  ///////////////////////////////////
  Int_t
  GeometryUtilities::GetPointOnLine(Double_t slope,
                                    Double_t intercept,
                                    Double_t wire1,
                                    Double_t time1,
                                    Double_t& wireout,
                                    Double_t& timeout) const
  {
    Double_t invslope = 0;

    if (slope) { invslope = -1. / slope; }

    Double_t ort_intercept = time1 - invslope * (Double_t)wire1;

    if ((slope - invslope) != 0)
      wireout = (ort_intercept - intercept) / (slope - invslope);
    else
      wireout = wire1;
    timeout = slope * wireout + intercept;

    return 0;
  }

  //////////////////////////////////
  // Calculate wire,time coordinates of the Hit projection onto a line
  //  all points are assumed to be in cm/cm space.
  ///////////////////////////////////
  int
  GeometryUtilities::GetPointOnLine(Double_t slope,
                                    const util::PxPoint* startpoint,
                                    const util::PxPoint* point1,
                                    util::PxPoint& pointout) const
  {

    double intercept = startpoint->t - slope * startpoint->w;

    return GetPointOnLine(slope, intercept, point1, pointout);
  }

  ///////////////////////////////////
  // Calculate wire,time coordinates of the Hit projection onto a line
  //  all points assumed to be in cm/cm space.
  ///////////////////////////////////
  int
  GeometryUtilities::GetPointOnLine(double slope,
                                    double intercept,
                                    const util::PxPoint* point1,
                                    util::PxPoint& pointout) const
  {
    double invslope = 0;

    if (slope) { invslope = -1. / slope; }

    double ort_intercept = point1->t - invslope * point1->w;

    if ((slope - invslope) != 0)
      pointout.w = (ort_intercept - intercept) / (slope - invslope);
    else
      pointout.w = point1->w;

    pointout.t = slope * pointout.w + intercept;

    return 0;
  }

  ///////////////////////////////////
  // Calculate wire,time coordinates of the Hit projection onto a line
  //
  ///////////////////////////////////
  Int_t
  GeometryUtilities::GetPointOnLine(Double_t slope,
                                    Double_t wirestart,
                                    Double_t timestart,
                                    Double_t wire1,
                                    Double_t time1,
                                    Double_t& wireout,
                                    Double_t& timeout) const
  {
    Double_t intercept = timestart - slope * (Double_t)wirestart;

    return GetPointOnLine(slope, intercept, wire1, time1, wireout, timeout);
  }

  ///////////////////////////////////
  // Calculate wire,time coordinates of the Hit projection onto a line
  //
  ///////////////////////////////////
  Int_t
  GeometryUtilities::GetPointOnLineWSlopes(Double_t slope,
                                           Double_t intercept,
                                           Double_t ort_intercept,
                                           Double_t& wireout,
                                           Double_t& timeout) const
  {
    Double_t invslope = 0;

    if (slope) { invslope = -1. / slope; }

    invslope *= fWireTimetoCmCm * fWireTimetoCmCm;

    wireout = (ort_intercept - intercept) / (slope - invslope);
    timeout = slope * wireout + intercept;

    wireout /= fWiretoCm;
    timeout /= fTimetoCm;

    return 0;
  }

  ///////////////////////////////////
  // Calculate wire,time coordinates of the Hit projection onto a line
  // slope should be in cm/cm space. PxPoint should be in cm/cm space.
  ///////////////////////////////////
  Int_t
  GeometryUtilities::GetPointOnLineWSlopes(double slope,
                                           double intercept,
                                           double ort_intercept,
                                           util::PxPoint& pointonline) const
  {
    Double_t invslope = 0;

    if (slope) invslope = -1. / slope;

    pointonline.w = (ort_intercept - intercept) / (slope - invslope);
    pointonline.t = slope * pointonline.w + intercept;
    return 0;
  }

  //////////////////////////////////////////////////////////
  Int_t
  GeometryUtilities::GetProjectedPoint(const PxPoint* p0, const PxPoint* p1, PxPoint& pN) const
  {

    // determine third plane number
    for (UInt_t i = 0; i < fNPlanes; ++i) {
      if (i == p0->plane || i == p1->plane) continue;
      pN.plane = i;
    }

    // Assuming there is no problem ( and we found the best pair that comes
    // close in time ) we try to get the Y and Z coordinates for the start of
    // the shower.
    UInt_t chan1 = fGeom.PlaneWireToChannel(p0->plane, p0->w);
    UInt_t chan2 = fGeom.PlaneWireToChannel(p1->plane, p1->w);
    const double origin[3] = {0.};
    Double_t pos[3] = {0.};
    fGeom.Plane(p0->plane).LocalToWorld(origin, pos);
    Double_t x = (p0->t - trigger_offset(fClocks)) * fTimetoCm + pos[0];

    Double_t y, z;
    if (!fGeom.ChannelsIntersect(chan1, chan2, y, z)) return -1;

    pos[0] = x;
    pos[1] = y;
    pos[2] = z;

    pN = Get2DPointProjection(pos, pN.plane);

    return 0;
  }

  //////////////////////////////////////////////////////////
  Int_t
  GeometryUtilities::GetYZ(const PxPoint* p0, const PxPoint* p1, Double_t* yz) const
  {
    Double_t y, z;

    // Force to the closest wires if not in the range
    int z0 = p0->w / fWiretoCm;
    int z1 = p1->w / fWiretoCm;
    if (z0 < 0) {
      std::cout << "\033[93mWarning\033[00m "
                   "\033[95m<<GeometryUtilities::GetYZ>>\033[00m"
                << std::endl
                << " 2D wire position " << p0->w << " [cm] corresponds to negative wire number."
                << std::endl
                << " Forcing it to wire=0..." << std::endl
                << "\033[93mWarning ends...\033[00m" << std::endl;
      z0 = 0;
    }
    else if (z0 >= (int)(fGeom.Nwires(p0->plane))) {
      std::cout << "\033[93mWarning\033[00m "
                   "\033[95m<<GeometryUtilities::GetYZ>>\033[00m"
                << std::endl
                << " 2D wire position " << p0->w << " [cm] exceeds max wire number "
                << (fGeom.Nwires(p0->plane) - 1) << std::endl
                << " Forcing it to the max wire number..." << std::endl
                << "\033[93mWarning ends...\033[00m" << std::endl;
      z0 = fGeom.Nwires(p0->plane) - 1;
    }
    if (z1 < 0) {
      std::cout << "\033[93mWarning\033[00m "
                   "\033[95m<<GeometryUtilities::GetYZ>>\033[00m"
                << std::endl
                << " 2D wire position " << p1->w << " [cm] corresponds to negative wire number."
                << std::endl
                << " Forcing it to wire=0..." << std::endl
                << "\033[93mWarning ends...\033[00m" << std::endl;
      z1 = 0;
    }
    if (z1 >= (int)(fGeom.Nwires(p1->plane))) {
      std::cout << "\033[93mWarning\033[00m "
                   "\033[95m<<GeometryUtilities::GetYZ>>\033[00m"
                << std::endl
                << " 2D wire position " << p1->w << " [cm] exceeds max wire number "
                << (fGeom.Nwires(p0->plane) - 1) << std::endl
                << " Forcing it to the max wire number..." << std::endl
                << "\033[93mWarning ends...\033[00m" << std::endl;
      z1 = fGeom.Nwires(p1->plane) - 1;
    }

    UInt_t chan1 = fGeom.PlaneWireToChannel(p0->plane, z0);
    UInt_t chan2 = fGeom.PlaneWireToChannel(p1->plane, z1);

    if (!fGeom.ChannelsIntersect(chan1, chan2, y, z)) return -1;

    yz[0] = y;
    yz[1] = z;

    return 0;
  }

  //////////////////////////////////////////////////////////
  Int_t
  GeometryUtilities::GetXYZ(const PxPoint* p0, const PxPoint* p1, Double_t* xyz) const
  {
    const double origin[3] = {0.};
    Double_t pos[3] = {0.};
    fGeom.Plane(p0->plane).LocalToWorld(origin, pos);
    Double_t x = (p0->t) - trigger_offset(fClocks) * fTimetoCm + pos[0];
    double yz[2];

    GetYZ(p0, p1, yz);

    xyz[0] = x;
    xyz[1] = yz[0];
    xyz[2] = yz[1];

    return 0;
  }

  //////////////////////////////////////////////////////////////

  PxPoint
  GeometryUtilities::Get2DPointProjection(Double_t* xyz, Int_t plane) const
  {

    PxPoint pN(0, 0, 0);
    const double origin[3] = {0.};
    Double_t pos[3];
    fGeom.Plane(plane).LocalToWorld(origin, pos);
    Double_t drifttick = (xyz[0] / fDriftVelocity) * (1. / fTimeTick);

    pos[1] = xyz[1];
    pos[2] = xyz[2];

    ///\todo: this should use the cryostat and tpc as well in the NearestWire
    /// method

    pN.w = fGeom.NearestWire(pos, plane);
    pN.t = drifttick - (pos[0] / fDriftVelocity) * (1. / fTimeTick) + trigger_offset(fClocks);
    pN.plane = plane;

    return pN;
  }

  //////////////////////////////////////////////////////////////
  // for now this returns the vlause in CM/CM space.
  // this will become the default, but don't want to break the code that depends
  // on the previous version. A.S. 03/26/14
  //////////////////////////////////////

  PxPoint
  GeometryUtilities::Get2DPointProjectionCM(std::vector<double> xyz, int plane) const
  {

    PxPoint pN(0, 0, 0);

    Double_t pos[3]{0., xyz[1], xyz[2]};

    ///\todo: this should use the cryostat and tpc as well in the NearestWire
    /// method

    pN.w = fGeom.NearestWire(pos, plane) * fWiretoCm;
    pN.t = xyz[0];
    pN.plane = plane;

    return pN;
  }

  PxPoint
  GeometryUtilities::Get2DPointProjectionCM(double* xyz, int plane) const
  {

    PxPoint pN(0, 0, 0);

    Double_t pos[3]{0., xyz[1], xyz[2]};

    ///\todo: this should use the cryostat and tpc as well in the NearestWire
    /// method

    pN.w = fGeom.NearestWire(pos, plane) * fWiretoCm;
    pN.t = xyz[0];
    pN.plane = plane;

    return pN;
  }

  PxPoint
  GeometryUtilities::Get2DPointProjectionCM(TLorentzVector* xyz, int plane) const
  {
    double xyznew[3] = {(*xyz)[0], (*xyz)[1], (*xyz)[2]};

    return Get2DPointProjectionCM(xyznew, plane);
  }

  Double_t
  GeometryUtilities::GetTimeTicks(Double_t x, Int_t plane) const
  {

    const double origin[3] = {0.};
    Double_t pos[3];
    fGeom.Plane(plane).LocalToWorld(origin, pos);
    Double_t drifttick = (x / fDriftVelocity) * (1. / fTimeTick);

    return drifttick - (pos[0] / fDriftVelocity) * (1. / fTimeTick) + trigger_offset(fClocks);
  }

  //----------------------------------------------------------------------
  // provide projected wire pitch for the view // copied from track.cxx and
  // modified
  Double_t
  GeometryUtilities::PitchInView(UInt_t plane, Double_t phi, Double_t theta) const
  {
    Double_t dirs[3] = {0.};
    GetDirectionCosines(phi, theta, dirs);

    /// \todo switch to using new Geometry::WireAngleToVertical(geo::View_t)
    /// \todo and Geometry::WirePitch(geo::View_t) methods
    Double_t wirePitch = 0.;
    Double_t angleToVert = 0.;

    wirePitch = fGeom.WirePitch(plane);
    angleToVert = fGeom.WireAngleToVertical(fGeom.Plane(plane).View()) - 0.5 * TMath::Pi();

    //(sin(angleToVert),std::cos(angleToVert)) is the direction perpendicular to
    // wire fDir.front() is the direction of the track at the beginning of its
    // trajectory
    Double_t cosgamma =
      TMath::Abs(TMath::Sin(angleToVert) * dirs[1] + TMath::Cos(angleToVert) * dirs[2]);

    if (cosgamma < 1.e-5)
    // throw UtilException("cosgamma is basically 0, that can't be right");
    {
      std::cout << " returning 100" << std::endl;
      return 100;
    }

    return wirePitch / cosgamma;
  }

  //////////////////////////////////////////////////
  void
  GeometryUtilities::GetDirectionCosines(Double_t phi, Double_t theta, Double_t* dirs) const
  {
    theta *= (TMath::Pi() / 180);
    phi *= (TMath::Pi() / 180); // working on copies, it's ok.
    dirs[0] = TMath::Cos(theta) * TMath::Sin(phi);
    dirs[1] = TMath::Sin(theta);
    dirs[2] = TMath::Cos(theta) * TMath::Cos(phi);
  }

  void
  GeometryUtilities::SelectLocalHitlist(const std::vector<util::PxHit>& hitlist,
                                        std::vector<const util::PxHit*>& hitlistlocal,
                                        util::PxPoint& startHit,
                                        Double_t& linearlimit,
                                        Double_t& ortlimit,
                                        Double_t& lineslopetest) const
  {
    util::PxHit testHit;
    SelectLocalHitlist(
      hitlist, hitlistlocal, startHit, linearlimit, ortlimit, lineslopetest, testHit);
  }

  //////////////////////////////////////////////////////////////////////////////////
  ////
  ///////////////////////////////////////////////////////////////////////////////////

  void
  GeometryUtilities::SelectLocalHitlist(const std::vector<util::PxHit>& hitlist,
                                        std::vector<const util::PxHit*>& hitlistlocal,
                                        util::PxPoint& startHit,
                                        Double_t& linearlimit,
                                        Double_t& ortlimit,
                                        Double_t& lineslopetest,
                                        util::PxHit& averageHit) const
  {

    hitlistlocal.clear();
    std::vector<unsigned int> hitlistlocal_index;

    hitlistlocal_index.clear();

    SelectLocalHitlistIndex(
      hitlist, hitlistlocal_index, startHit, linearlimit, ortlimit, lineslopetest);

    double timesum = 0;
    double wiresum = 0;
    for (size_t i = 0; i < hitlistlocal_index.size(); ++i) {

      hitlistlocal.push_back((const util::PxHit*)(&(hitlist.at(hitlistlocal_index.at(i)))));
      timesum += hitlist.at(hitlistlocal_index.at(i)).t;
      wiresum += hitlist.at(hitlistlocal_index.at(i)).w;
    }

    averageHit.plane = startHit.plane;
    if (hitlistlocal.size()) {
      averageHit.w = wiresum / (double)hitlistlocal.size();
      averageHit.t = timesum / ((double)hitlistlocal.size());
    }
  }

  void
  GeometryUtilities::SelectLocalHitlistIndex(const std::vector<util::PxHit>& hitlist,
                                             std::vector<unsigned int>& hitlistlocal_index,
                                             util::PxPoint& startHit,
                                             Double_t& linearlimit,
                                             Double_t& ortlimit,
                                             Double_t& lineslopetest) const
  {

    hitlistlocal_index.clear();
    double locintercept = startHit.t - startHit.w * lineslopetest;

    for (size_t i = 0; i < hitlist.size(); ++i) {

      util::PxPoint hitonline;

      GetPointOnLine(lineslopetest, locintercept, (const util::PxHit*)(&hitlist.at(i)), hitonline);

      // calculate linear distance from start point and orthogonal distance from
      // axis
      Double_t lindist =
        Get2DDistance((const util::PxPoint*)(&hitonline), (const util::PxPoint*)(&startHit));
      Double_t ortdist =
        Get2DDistance((const util::PxPoint*)(&hitlist.at(i)), (const util::PxPoint*)(&hitonline));

      if (lindist < linearlimit && ortdist < ortlimit) { hitlistlocal_index.push_back(i); }
    }
  }

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  void
  GeometryUtilities::SelectPolygonHitList(std::vector<util::PxHit> const& hitlist,
                                          std::vector<util::PxHit const*>& hitlistlocal) const
  {
    if (empty(hitlist)) { throw UtilException("Provided empty hit list!"); }

    hitlistlocal.clear();
    unsigned char plane = hitlist.front().plane;

    // Define subset of hits to define polygon
    std::map<double, const util::PxHit*> hitmap;
    double qtotal = 0;
    for (auto const& h : hitlist) {
      hitmap.try_emplace(h.charge, &h);
      qtotal += h.charge;
    }
    double qintegral = 0;
    std::vector<const util::PxHit*> ordered_hits;
    ordered_hits.reserve(hitlist.size());
    for (auto hiter = hitmap.rbegin(); qintegral < qtotal * 0.95 && hiter != hitmap.rend();
         ++hiter) {

      qintegral += (*hiter).first;
      ordered_hits.push_back((*hiter).second);
    }

    // Define container to hold found polygon corner PxHit index & distance
    std::vector<size_t> hit_index(8, 0);
    std::vector<double> hit_distance(8, 1e9);

    // Loop over hits and find corner points in the plane view
    // Also fill corner edge points
    std::vector<util::PxPoint> edges(4, PxPoint(plane, 0, 0));
    double wire_max = fGeom.Nwires(plane) * fWiretoCm;
    double time_max = fDetProp.NumberTimeSamples() * fTimetoCm;

    for (size_t index = 0; index < ordered_hits.size(); ++index) {

      if (ordered_hits.at(index)->t < 0 || ordered_hits.at(index)->w < 0 ||
          ordered_hits.at(index)->t > time_max || ordered_hits.at(index)->w > wire_max) {

        throw UtilException(Form("Invalid wire/time (%g,%g) ... range is (0=>%g,0=>%g)",
                                 ordered_hits.at(index)->w,
                                 ordered_hits.at(index)->t,
                                 wire_max,
                                 time_max));
      }

      double dist = 0;

      // Comparison w/ (Wire,0)
      dist = ordered_hits.at(index)->t;
      if (dist < hit_distance.at(1)) {
        hit_distance.at(1) = dist;
        hit_index.at(1) = index;
        edges.at(0).t = ordered_hits.at(index)->t;
        edges.at(1).t = ordered_hits.at(index)->t;
      }

      // Comparison w/ (WireMax,Time)
      dist = wire_max - ordered_hits.at(index)->w;
      if (dist < hit_distance.at(3)) {
        hit_distance.at(3) = dist;
        hit_index.at(3) = index;
        edges.at(1).w = ordered_hits.at(index)->w;
        edges.at(2).w = ordered_hits.at(index)->w;
      }

      // Comparison w/ (Wire,TimeMax)
      dist = time_max - ordered_hits.at(index)->t;
      if (dist < hit_distance.at(5)) {
        hit_distance.at(5) = dist;
        hit_index.at(5) = index;
        edges.at(2).t = ordered_hits.at(index)->t;
        edges.at(3).t = ordered_hits.at(index)->t;
      }

      // Comparison w/ (0,Time)
      dist = ordered_hits.at(index)->w;
      if (dist < hit_distance.at(7)) {
        hit_distance.at(7) = dist;
        hit_index.at(7) = index;
        edges.at(0).w = ordered_hits.at(index)->w;
        edges.at(3).w = ordered_hits.at(index)->w;
      }
    }
    for (size_t index = 0; index < ordered_hits.size(); ++index) {

      double dist = 0;
      // Comparison w/ (0,0)
      dist = cet::sum_of_squares(ordered_hits.at(index)->t - edges.at(0).t,
                                 ordered_hits.at(index)->w - edges.at(0).w);
      if (dist < hit_distance.at(0)) {
        hit_distance.at(0) = dist;
        hit_index.at(0) = index;
      }

      // Comparison w/ (WireMax,0)
      dist = cet::sum_of_squares(ordered_hits.at(index)->t - edges.at(1).t,
                                 ordered_hits.at(index)->w - edges.at(1).w);
      if (dist < hit_distance.at(2)) {
        hit_distance.at(2) = dist;
        hit_index.at(2) = index;
      }

      // Comparison w/ (WireMax,TimeMax)
      dist = cet::sum_of_squares(ordered_hits.at(index)->t - edges.at(2).t,
                                 ordered_hits.at(index)->w - edges.at(2).w);
      if (dist < hit_distance.at(4)) {
        hit_distance.at(4) = dist;
        hit_index.at(4) = index;
      }

      // Comparison w/ (0,TimeMax)
      dist = cet::sum_of_squares(ordered_hits.at(index)->t - edges.at(3).t,
                                 ordered_hits.at(index)->w - edges.at(3).w);
      if (dist < hit_distance.at(6)) {
        hit_distance.at(6) = dist;
        hit_index.at(6) = index;
      }
    }

    // Loop over the resulting hit indexes and append unique hits to define the
    // polygon to the return hit list
    std::set<size_t> unique_index;
    std::vector<size_t> candidate_polygon;
    candidate_polygon.reserve(9);
    for (auto& index : hit_index) {
      if (unique_index.find(index) == unique_index.end()) {
        unique_index.insert(index);
        candidate_polygon.push_back(index);
      }
    }
    for (auto& index : hit_index) {
      candidate_polygon.push_back(index);
      break;
    }

    if (unique_index.size() > 8) throw UtilException("Size of the polygon > 8!");

    // Untangle Polygon
    candidate_polygon = PolyOverlap(ordered_hits, candidate_polygon);

    hitlistlocal.clear();
    for (unsigned int i = 0; i < (candidate_polygon.size() - 1); i++) {
      hitlistlocal.push_back((const util::PxHit*)(ordered_hits.at(candidate_polygon.at(i))));
    }
    // check that polygon does not have more than 8 sides
    if (unique_index.size() > 8) throw UtilException("Size of the polygon > 8!");
  }

  std::vector<size_t>
  GeometryUtilities::PolyOverlap(std::vector<const util::PxHit*> ordered_hits,
                                 std::vector<size_t> candidate_polygon) const
  {
    // loop over edges
    for (unsigned int i = 0; i < (candidate_polygon.size() - 1); i++) {
      double Ax = ordered_hits.at(candidate_polygon.at(i))->w;
      double Ay = ordered_hits.at(candidate_polygon.at(i))->t;
      double Bx = ordered_hits.at(candidate_polygon.at(i + 1))->w;
      double By = ordered_hits.at(candidate_polygon.at(i + 1))->t;
      // loop over edges that have not been checked yet...
      // only ones furhter down in polygon
      for (unsigned int j = i + 2; j < (candidate_polygon.size() - 1); j++) {
        // avoid consecutive segments:
        if (candidate_polygon.at(i) == candidate_polygon.at(j + 1))
          continue;
        else {
          double Cx = ordered_hits.at(candidate_polygon.at(j))->w;
          double Cy = ordered_hits.at(candidate_polygon.at(j))->t;
          double Dx = ordered_hits.at(candidate_polygon.at(j + 1))->w;
          double Dy = ordered_hits.at(candidate_polygon.at(j + 1))->t;

          if ((Clockwise(Ax, Ay, Cx, Cy, Dx, Dy) != Clockwise(Bx, By, Cx, Cy, Dx, Dy)) and
              (Clockwise(Ax, Ay, Bx, By, Cx, Cy) != Clockwise(Ax, Ay, Bx, By, Dx, Dy))) {
            size_t tmp = candidate_polygon.at(i + 1);
            candidate_polygon.at(i + 1) = candidate_polygon.at(j);
            candidate_polygon.at(j) = tmp;
            // check that last element is still first (to close circle...)
            candidate_polygon.at(candidate_polygon.size() - 1) = candidate_polygon.at(0);
            // swapped polygon...now do recursion to make sure
            return PolyOverlap(ordered_hits, candidate_polygon);
          } // if crossing
        }
      } // second loop
    }   // first loop
    return candidate_polygon;
  }

  bool
  GeometryUtilities::Clockwise(double const Ax,
                               double const Ay,
                               double const Bx,
                               double const By,
                               double const Cx,
                               double const Cy) const
  {
    return (Cy - Ay) * (Bx - Ax) > (By - Ay) * (Cx - Ax);
  }

  util::PxHit
  GeometryUtilities::FindClosestHit(std::vector<util::PxHit> const& hitlist,
                                    unsigned int const wirein,
                                    double const timein) const
  {
    return hitlist[FindClosestHitIndex(hitlist, wirein, timein)];
  }

  unsigned int
  GeometryUtilities::FindClosestHitIndex(std::vector<util::PxHit> const& hitlist,
                                         unsigned int const wirein,
                                         double const timein) const
  {
    double min_length_from_start = 99999;
    unsigned int ret_ind = 0;

    for (unsigned int ii = 0; ii < hitlist.size(); ii++) {
      util::PxHit const& hit = hitlist[ii];
      double const dist_mod = Get2DDistance(wirein, timein, hit.w, hit.t);
      if (dist_mod < min_length_from_start) {
        min_length_from_start = dist_mod;
        ret_ind = ii;
      }
    }

    return ret_ind;
  }

} // namespace
