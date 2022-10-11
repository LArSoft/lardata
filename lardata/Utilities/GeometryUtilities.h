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

  constexpr double kINVALID_DOUBLE = std::numeric_limits<Double_t>::max();

  class GeometryUtilities {
  public:
    GeometryUtilities(geo::GeometryCore const& geom,
                      detinfo::DetectorClocksData const& clockData,
                      detinfo::DetectorPropertiesData const& propData);

    Int_t Get3DaxisN(Int_t iplane0,
                     Int_t iplane1,
                     Double_t omega0,
                     Double_t omega1,
                     Double_t& phi,
                     Double_t& theta) const;

    Double_t CalculatePitch(UInt_t iplane0, Double_t phi, Double_t theta) const;

    Double_t CalculatePitchPolar(UInt_t iplane0, Double_t phi, Double_t theta) const;

    Double_t Get3DSpecialCaseTheta(Int_t iplane0, Int_t iplane1, Double_t dw0, Double_t dw1) const;

    Double_t Get2Dangle(Double_t deltawire, Double_t deltatime) const;

    Double_t Get2Dangle(Double_t wireend,
                        Double_t wirestart,
                        Double_t timeend,
                        Double_t timestart) const;

    double Get2Dangle(const util::PxPoint* endpoint, const util::PxPoint* startpoint) const;

    double Get2DangleFrom3D(unsigned int plane, double phi, double theta) const;

    double Get2DangleFrom3D(unsigned int plane, TVector3 dir_vector) const;

    Double_t Get2Dslope(Double_t deltawire, Double_t deltatime) const;

    Double_t Get2Dslope(Double_t wireend,
                        Double_t wirestart,
                        Double_t timeend,
                        Double_t timestart) const;

    double Get2Dslope(const util::PxPoint* endpoint, const util::PxPoint* startpoint) const;

    Double_t Get2DDistance(Double_t wire1, Double_t time1, Double_t wire2, Double_t time2) const;

    double Get2DDistance(const util::PxPoint* point1, const util::PxPoint* point2) const;

    Double_t Get2DPitchDistance(Double_t angle, Double_t inwire, Double_t wire) const;

    Double_t Get2DPitchDistanceWSlope(Double_t slope, Double_t inwire, Double_t wire) const;

    Int_t GetPointOnLine(Double_t slope,
                         Double_t intercept,
                         Double_t wire1,
                         Double_t time1,
                         Double_t& wireout,
                         Double_t& timeout) const;

    Int_t GetPointOnLine(Double_t slope,
                         Double_t wirestart,
                         Double_t timestart,
                         Double_t wire1,
                         Double_t time1,
                         Double_t& wireout,
                         Double_t& timeout) const;

    int GetPointOnLine(Double_t slope,
                       const util::PxPoint* startpoint,
                       const util::PxPoint* point1,
                       util::PxPoint& pointout) const;

    int GetPointOnLine(double slope,
                       double intercept,
                       const util::PxPoint* point1,
                       util::PxPoint& pointout) const;

    Int_t GetPointOnLineWSlopes(Double_t slope,
                                Double_t intercept,
                                Double_t ort_intercept,
                                Double_t& wireout,
                                Double_t& timeout) const;

    Int_t GetPointOnLineWSlopes(double slope,
                                double intercept,
                                double ort_intercept,
                                util::PxPoint& pointonline) const;

    PxPoint Get2DPointProjection(Double_t* xyz, Int_t plane) const;

    PxPoint Get2DPointProjectionCM(std::vector<double> xyz, int plane) const;

    PxPoint Get2DPointProjectionCM(double* xyz, int plane) const;

    PxPoint Get2DPointProjectionCM(TLorentzVector* xyz, int plane) const;

    Double_t GetTimeTicks(Double_t x, Int_t plane) const;

    Int_t GetProjectedPoint(const PxPoint* p0, const PxPoint* p1, PxPoint& pN) const;

    util::PxHit FindClosestHit(std::vector<util::PxHit> const& hitlist,
                               unsigned int wirein,
                               double timein) const;

    unsigned int FindClosestHitIndex(std::vector<util::PxHit> const& hitlist,
                                     unsigned int wirein,
                                     double timein) const;

    Int_t GetYZ(const PxPoint* p0, const PxPoint* p1, Double_t* yz) const;

    Int_t GetXYZ(const PxPoint* p0, const PxPoint* p1, Double_t* xyz) const;

    Double_t PitchInView(UInt_t plane, Double_t phi, Double_t theta) const;

    void GetDirectionCosines(Double_t phi, Double_t theta, Double_t* dirs) const;

    // interface without average Hit
    void SelectLocalHitlist(const std::vector<util::PxHit>& hitlist,
                            std::vector<const util::PxHit*>& hitlistlocal,
                            util::PxPoint& startHit,
                            Double_t& linearlimit,
                            Double_t& ortlimit,
                            Double_t& lineslopetest) const;

    void SelectLocalHitlist(const std::vector<util::PxHit>& hitlist,
                            std::vector<const util::PxHit*>& hitlistlocal,
                            util::PxPoint& startHit,
                            Double_t& linearlimit,
                            Double_t& ortlimit,
                            Double_t& lineslopetest,
                            util::PxHit& averageHit) const;

    void SelectLocalHitlistIndex(const std::vector<util::PxHit>& hitlist,
                                 std::vector<unsigned int>& hitlistlocal_index,
                                 util::PxPoint& startHit,
                                 Double_t& linearlimit,
                                 Double_t& ortlimit,
                                 Double_t& lineslopetest) const;

    void SelectPolygonHitList(const std::vector<util::PxHit>& hitlist,
                              std::vector<const util::PxHit*>& hitlistlocal) const;

    std::vector<size_t> PolyOverlap(std::vector<const util::PxHit*> ordered_hits,
                                    std::vector<size_t> candidate_polygon) const;

    bool Clockwise(double Ax, double Ay, double Bx, double By, double Cx, double Cy) const;

    Double_t TimeToCm() const { return fTimetoCm; }
    Double_t WireToCm() const { return fWiretoCm; }
    Double_t WireTimeToCmCm() const { return fWireTimetoCmCm; }
    UInt_t Nplanes() const { return fNPlanes; }

  private:
    geo::GeometryCore const& fGeom;
    detinfo::DetectorClocksData const& fClocks;
    detinfo::DetectorPropertiesData const& fDetProp;

    std::vector<Double_t> vertangle; // angle wrt to vertical
    Double_t fWirePitch;
    Double_t fTimeTick;
    Double_t fDriftVelocity;
    UInt_t fNPlanes;
    Double_t fWiretoCm;
    Double_t fTimetoCm;
    Double_t fWireTimetoCmCm;
  }; // class GeometryUtilities

} // namespace util
#endif // UTIL_GEOMETRYUTILITIES_H
