////////////////////////////////////////////////////////////////////////
// \file DetectorPropertiesArgoNeuT.h
//
// \brief service to contain information about detector electronics, etc
//
// \author brebel@fnal.gov
//
// From the original DetectorProperties.h ; this one preserves the dependency on
// DatabaseUtil service and the ability to read information from a database
// with direct DB connection.
// For new experiments, an indirect connection should be used instead.
//
// PLEASE DO NOT take this as a model to develop a service:
// this is just a backward-compatible hack.
//
////////////////////////////////////////////////////////////////////////
#ifndef UTIL_DETECTORPROPERTIESARGONEUT_H
#define UTIL_DETECTORPROPERTIESARGONEUT_H

#include "lardata/Utilities/LArPropertiesServiceArgoNeuT.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataalg/DetectorInfo/ElecClock.h"

#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace detinfo {
  class DetectorClocks;
}

/// General LArSoft Utilities
namespace util {
  class DetectorPropertiesArgoNeuT final : public detinfo::DetectorProperties {
  public:
    DetectorPropertiesArgoNeuT(fhicl::ParameterSet const& pset);

    void
    SetNumberTimeSamples(unsigned int nsamp)
    {
      fNumberTimeSamples = nsamp;
    }

    double Efield(unsigned int planegap = 0) const override ///< kV/cm
    {
      return fLP->Efield(planegap);
    }

    double
    DriftVelocity(double efield = 0., double temperature = 0.) const override
    {
      return fLP->DriftVelocity(efield, temperature);
    }

    /// dQ/dX in electrons/cm, returns dE/dX in MeV/cm.
    double
    BirksCorrection(double dQdX) const override
    {
      return fLP->BirksCorrection(dQdX);
    }
    double
    BirksCorrection(double dQdX, double EField) const override
    {
      return fLP->BirksCorrection(dQdX, EField);
    }

    double
    ModBoxCorrection(double dQdX) const override
    {
      return fLP->ModBoxCorrection(dQdX);
    }
    double
    ModBoxCorrection(double dQdX, double EField) const override
    {
      return fLP->ModBoxCorrection(dQdX, EField);
    }

    double
    ElectronLifetime() const override
    {
      return fLP->ElectronLifetime();
    }

    double
    Temperature() const override
    {
      return fLP->Temperature();
    }

    double
    Density(double temperature) const override
    {
      return fLP->Density(temperature);
    }
    double
    Density() const override
    {
      return Density(Temperature());
    }

    double
    Eloss(double mom, double mass, double tcut) const override
    {
      return fLP->Eloss(mom, mass, tcut);
    }

    double
    ElossVar(double mom, double mass) const override
    {
      return fLP->ElossVar(mom, mass);
    }

    double
    ElectronsToADC() const override
    {
      return fElectronsToADC;
    }
    unsigned int
    NumberTimeSamples() const override
    {
      return fNumberTimeSamples;
    }
    unsigned int
    ReadOutWindowSize() const override
    {
      return fReadOutWindowSize;
    }
    double
    TimeOffsetU() const override
    {
      return fTimeOffsetU;
    }
    double
    TimeOffsetV() const override
    {
      return fTimeOffsetV;
    }
    double
    TimeOffsetZ() const override
    {
      return fTimeOffsetZ;
    }
    double
    TimeOffsetY() const override
    {
      return 0;
    }

    bool
    SimpleBoundary() const override
    {
      return fSimpleBoundary;
    }

    detinfo::DetectorPropertiesData
    DataFor(detinfo::DetectorClocksData const& clock_data) const override;

  private:
    double fElectronsToADC;          ///< conversion factor for # of ionization electrons
                                     ///< to 1 ADC count
    unsigned int fNumberTimeSamples; ///< number of clock ticks per event
    unsigned int fReadOutWindowSize; ///< number of clock ticks per readout window
    double fTimeOffsetU;             ///< time offsets to convert spacepoint
    double fTimeOffsetV;             ///< coordinates to hit times on each
    double fTimeOffsetZ;             ///< view

    /// Pointer to the specific LArPropertiesServiceArgoNeuT service (provider)
    util::LArPropertiesServiceArgoNeuT const* fLP;

    bool fSimpleBoundary;

  }; // class DetectorPropertiesArgoNeuT
} // namespace util

#endif // UTIL_DETECTORPROPERTIESARGONEUT_H
