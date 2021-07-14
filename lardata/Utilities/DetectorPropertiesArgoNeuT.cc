////////////////////////////////////////////////////////////////////////
//
//  \file DetectorPropertiesServiceArgoNeuT_service.cc
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/Utilities/DatabaseUtil.h"
#include "lardata/Utilities/DetectorPropertiesServiceArgoNeuT.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Art includes
#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "canvas/Utilities/Exception.h"

namespace util {

  //--------------------------------------------------------------------
  DetectorPropertiesArgoNeuT::DetectorPropertiesArgoNeuT(fhicl::ParameterSet const& pset)
  {
    fLP = dynamic_cast<util::LArPropertiesServiceArgoNeuT const*>(
      lar::providerFrom<detinfo::LArPropertiesService>());
    if (!fLP) {
      // this legacy service works only coupled to the corresponding
      // LArProperties legacy service:
      throw art::Exception(art::errors::Configuration)
        << "DetectorPropertiesArgoNeuT service requires"
           " LArPropertiesServiceArgoNeuT";
    }

    if (pset.has_key("SamplingRate"))
      throw cet::exception(__FUNCTION__) << "SamplingRate is a deprecated fcl parameter for "
                                            "DetectorPropertiesArgoNeuT!";
    if (pset.has_key("TriggerOffset"))
      throw cet::exception(__FUNCTION__) << "TriggerOffset is a deprecated fcl parameter for "
                                            "DetectorPropertiesArgoNeuT!";
    fElectronsToADC = pset.get<double>("ElectronsToADC");
    fNumberTimeSamples = pset.get<unsigned int>("NumberTimeSamples");
    fReadOutWindowSize = pset.get<unsigned int>("ReadOutWindowSize");
    fTimeOffsetU = pset.get<double>("TimeOffsetU");
    fTimeOffsetV = pset.get<double>("TimeOffsetV");
    fTimeOffsetZ = pset.get<double>("TimeOffsetZ");

    fSimpleBoundary = pset.get<bool>("SimpleBoundaryProcess", true);
  }

  detinfo::DetectorPropertiesData
  DetectorPropertiesArgoNeuT::DataFor(
    detinfo::DetectorClocksData const& clock_data) const
  {
    art::ServiceHandle<geo::Geometry const> geo;
    auto const larProp = dynamic_cast<util::LArPropertiesServiceArgoNeuT const*>(
      lar::providerFrom<detinfo::LArPropertiesService>());
    assert(larProp);

    double samplingRate = sampling_rate(clock_data);
    double efield = Efield();
    double temperature = larProp->Temperature();
    double driftVelocity = DriftVelocity(efield, temperature);
    double const x_ticks_coefficient = 0.001 * driftVelocity * samplingRate;

    double triggerOffset = trigger_offset(clock_data);

    std::vector<std::vector<std::vector<double>>> x_ticks_offsets(geo->Ncryostats());
    std::vector<std::vector<double>> drift_direction(geo->Ncryostats());

    for (size_t cstat = 0; cstat < geo->Ncryostats(); ++cstat) {
      x_ticks_offsets[cstat].resize(geo->Cryostat(cstat).NTPC());
      drift_direction[cstat].resize(geo->Cryostat(cstat).NTPC());

      for (size_t tpc = 0; tpc < geo->Cryostat(cstat).NTPC(); ++tpc) {
        const geo::TPCGeo& tpcgeom = geo->Cryostat(cstat).TPC(tpc);

        const double dir((tpcgeom.DriftDirection() == geo::kNegX) ? +1.0 : -1.0);
        drift_direction[cstat][tpc] = dir;

        int nplane = tpcgeom.Nplanes();
        x_ticks_offsets[cstat][tpc].resize(nplane, 0.);
        for (int plane = 0; plane < nplane; ++plane) {
          const geo::PlaneGeo& pgeom = tpcgeom.Plane(plane);

          // Get field in gap between planes
          double efieldgap[3];
          double driftVelocitygap[3];
          double x_ticks_coeff_gap[3];
          for (int igap = 0; igap < 3; ++igap) {
            efieldgap[igap] = Efield(igap);
            driftVelocitygap[igap] = DriftVelocity(efieldgap[igap], temperature);
            x_ticks_coeff_gap[igap] = 0.001 * driftVelocitygap[igap] * samplingRate;
          }

          // Calculate geometric time offset.
          // only works if xyz[0]<=0
          const double* xyz = tpcgeom.PlaneLocation(0);

          x_ticks_offsets[cstat][tpc][plane] =
            -xyz[0] / (dir * x_ticks_coefficient) + triggerOffset;

          if (nplane == 3) {
            /*
         |    ---------- plane = 2 (collection)
         |                      Coeff[2]
         |    ---------- plane = 1 (2nd induction)
         |                      Coeff[1]
         |    ---------- plane = 0 (1st induction) x = xyz[0]
         |                      Coeff[0]
         |    ---------- x = 0
         V     For plane = 0, t offset is -xyz[0]/Coeff[0]
         x   */
            for (int ip = 0; ip < plane; ++ip) {
              x_ticks_offsets[cstat][tpc][plane] +=
                tpcgeom.PlanePitch(ip, ip + 1) / x_ticks_coeff_gap[ip + 1];
            }
          }
          else if (nplane == 2) { ///< special case for ArgoNeuT
            /*
         |    ---------- plane = 1 (collection)
         |                      Coeff[2]
         |    ---------- plane = 0 (2nd induction) x = xyz[0]
         |    ---------- x = 0, Coeff[1]
         V    ---------- first induction plane
         x                      Coeff[0]
               For plane = 0, t offset is pitch/Coeff[1] - (pitch+xyz[0])/Coeff[0]
                         = -xyz[0]/Coeff[0] - pitch*(1/Coeff[0]-1/Coeff[1])
            */
            for (int ip = 0; ip < plane; ++ip) {
              x_ticks_offsets[cstat][tpc][plane] +=
                tpcgeom.PlanePitch(ip, ip + 1) / x_ticks_coeff_gap[ip + 2];
            }
            x_ticks_offsets[cstat][tpc][plane] -=
              tpcgeom.PlanePitch() * (1 / x_ticks_coefficient - 1 / x_ticks_coeff_gap[1]);
          }

          // Add view dependent offset
          switch (pgeom.View()) {
          case geo::kU: x_ticks_offsets[cstat][tpc][plane] += TimeOffsetU(); break;
          case geo::kV: x_ticks_offsets[cstat][tpc][plane] += TimeOffsetV(); break;
          case geo::kZ: x_ticks_offsets[cstat][tpc][plane] += TimeOffsetZ(); break;
          default:
            throw cet::exception("DetectorPropertiesServiceArgoNeuT")
              << "Bad view = " << pgeom.View() << "\n";
          }
        }
      }
    }

    return detinfo::DetectorPropertiesData{*this,
                                           x_ticks_coefficient,
                                           move(x_ticks_offsets),
                                           move(drift_direction)};
  }
}
