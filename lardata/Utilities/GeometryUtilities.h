////////////////////////////////////////////////////////////////////////
// \file GeometryUtilities.h
//
// \brief Functions to calculate distances and angles in 3D and 2D
//
// \author andrzej.szelc@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef UTIL_GEOMETRYUTILITIES_H
#define UTIL_GEOMETRYUTILITIES_H

#include "RtypesCore.h"
#include "TVector3.h"

#include "PxUtils.h"

#include <limits>
#include <vector>

class TLorentzVector;

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}
namespace geo {
  class GeometryCore;
}

/// General LArSoft Utilities
namespace util {

  constexpr double kINVALID_DOUBLE = std::numeric_limits<double>::max();

  class GeometryUtilities {
  public:
    GeometryUtilities(geo::GeometryCore const& geom,
                      detinfo::DetectorClocksData const& clockData,
                      detinfo::DetectorPropertiesData const& propData);

    int Get3DaxisN(int iplane0,
                   int iplane1,
                   double omega0,
                   double omega1,
                   double& phi,
                   double& theta) const;

    double CalculatePitch(unsigned int iplane0, double phi, double theta) const;

    double CalculatePitchPolar(unsigned int iplane0, double phi, double theta) const;

    double Get3DSpecialCaseTheta(int iplane0, int iplane1, double dw0, double dw1) const;

    double Get2Dangle(double deltawire, double deltatime) const;

    double Get2Dangle(double wireend, double wirestart, double timeend, double timestart) const;

    double Get2Dangle(const PxPoint* endpoint, const PxPoint* startpoint) const;

    double Get2DangleFrom3D(unsigned int plane, double phi, double theta) const;

    double Get2DangleFrom3D(unsigned int plane, TVector3 dir_vector) const;

    double Get2Dslope(double deltawire, double deltatime) const;

    double Get2Dslope(double wireend, double wirestart, double timeend, double timestart) const;

    double Get2Dslope(const PxPoint* endpoint, const PxPoint* startpoint) const;

    double Get2DDistance(double wire1, double time1, double wire2, double time2) const;

    double Get2DDistance(const PxPoint* point1, const PxPoint* point2) const;

    double Get2DPitchDistance(double angle, double inwire, double wire) const;

    double Get2DPitchDistanceWSlope(double slope, double inwire, double wire) const;

    int GetPointOnLine(double slope,
                       double intercept,
                       double wire1,
                       double time1,
                       double& wireout,
                       double& timeout) const;

    int GetPointOnLine(double slope,
                       double wirestart,
                       double timestart,
                       double wire1,
                       double time1,
                       double& wireout,
                       double& timeout) const;

    int GetPointOnLine(double slope,
                       const PxPoint* startpoint,
                       const PxPoint* point1,
                       PxPoint& pointout) const;

    int GetPointOnLine(double slope,
                       double intercept,
                       const PxPoint* point1,
                       PxPoint& pointout) const;

    int GetPointOnLineWSlopes(double slope,
                              double intercept,
                              double ort_intercept,
                              double& wireout,
                              double& timeout) const;

    int GetPointOnLineWSlopes(double slope,
                              double intercept,
                              double ort_intercept,
                              PxPoint& pointonline) const;

    PxPoint Get2DPointProjection(double const* xyz, unsigned int plane) const;

    PxPoint Get2DPointProjectionCM(std::vector<double> const& xyz, unsigned int plane) const;

    PxPoint Get2DPointProjectionCM(double const* xyz, unsigned int plane) const;

    PxPoint Get2DPointProjectionCM(TLorentzVector const* xyz, unsigned int plane) const;

    double GetTimeTicks(double x, unsigned int plane) const;

    int GetProjectedPoint(const PxPoint* p0, const PxPoint* p1, PxPoint& pN) const;

    PxHit FindClosestHit(std::vector<PxHit> const& hitlist,
                         unsigned int wirein,
                         double timein) const;

    unsigned int FindClosestHitIndex(std::vector<PxHit> const& hitlist,
                                     unsigned int wirein,
                                     double timein) const;

    int GetYZ(const PxPoint* p0, const PxPoint* p1, double* yz) const;

    int GetXYZ(const PxPoint* p0, const PxPoint* p1, double* xyz) const;

    double PitchInView(unsigned int plane, double phi, double theta) const;

    void GetDirectionCosines(double phi, double theta, double* dirs) const;

    // interface without average Hit
    void SelectLocalHitlist(const std::vector<PxHit>& hitlist,
                            std::vector<const PxHit*>& hitlistlocal,
                            PxPoint& startHit,
                            double& linearlimit,
                            double& ortlimit,
                            double& lineslopetest) const;

    void SelectLocalHitlist(const std::vector<PxHit>& hitlist,
                            std::vector<const PxHit*>& hitlistlocal,
                            PxPoint& startHit,
                            double& linearlimit,
                            double& ortlimit,
                            double& lineslopetest,
                            PxHit& averageHit) const;

    void SelectLocalHitlistIndex(const std::vector<PxHit>& hitlist,
                                 std::vector<unsigned int>& hitlistlocal_index,
                                 PxPoint& startHit,
                                 double& linearlimit,
                                 double& ortlimit,
                                 double& lineslopetest) const;

    void SelectPolygonHitList(const std::vector<PxHit>& hitlist,
                              std::vector<const PxHit*>& hitlistlocal) const;

    std::vector<size_t> PolyOverlap(std::vector<const PxHit*> ordered_hits,
                                    std::vector<size_t> candidate_polygon) const;

    bool Clockwise(double Ax, double Ay, double Bx, double By, double Cx, double Cy) const;

    double TimeToCm() const { return fTimetoCm; }
    double WireToCm() const { return fWiretoCm; }
    double WireTimeToCmCm() const { return fWireTimetoCmCm; }
    unsigned int Nplanes() const { return fNPlanes; }

  private:
    geo::GeometryCore const& fGeom;
    detinfo::DetectorClocksData const& fClocks;
    detinfo::DetectorPropertiesData const& fDetProp;

    std::vector<double> vertangle; // angle wrt to vertical
    double fWirePitch;
    double fTimeTick;
    double fDriftVelocity;
    unsigned int fNPlanes;
    double fWiretoCm;
    double fTimetoCm;
    double fWireTimetoCmCm;
  }; // class GeometryUtilities

} // namespace util
#endif // UTIL_GEOMETRYUTILITIES_H
