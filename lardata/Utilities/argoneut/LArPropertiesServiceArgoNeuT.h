////////////////////////////////////////////////////////////////////////
// LArPropertiesServiceArgoNeuT.h
//
// Utility LAr functions
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// msoderbe@syr.edu
// joshua.spitz@yale.edu
//
// Optical Properties:
// bjpjones@mit.edu
//
// From the original LArProperties.h ; this one preserves the dependency on
// DatabaseUtil service and the ability to read information from a database
// with direct DB connection.
// For new experiments, an indirect connection should be used instead.
//
// PLEASE DO NOT take this as a model to develop a service:
// this is just a backward-compatible hack.
//
////////////////////////////////////////////////////////////////////////
#ifndef LARPROPERTIESSERVICEARGONEUT_H
#define LARPROPERTIESSERVICEARGONEUT_H

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/LArProperties.h"

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "string.h"

///General LArSoft Utilities
namespace util {
  /**
     * @brief Properties related to liquid argon environment in the detector
     *
     * This class can access databases via DatabaseUtil service.
     *
     * @note Some of the database connection properties are established before
     * the beginning of the job and if they change this service will not be
     * aware of it. These properties petrain, so far, only the connection mode
     * and not any content of the databases themselves.
     */
  class LArPropertiesServiceArgoNeuT
    : public detinfo::LArProperties // implements provider interface
    ,
      public detinfo::LArPropertiesService // implements service interface
  {
  public:
    //------------------------------------------------------------------------
    LArPropertiesServiceArgoNeuT(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    //------------------------------------------------------------------------
    //--- service interface
    /// Type of service provider
    using provider_type = LArPropertiesServiceArgoNeuT;

    virtual ~LArPropertiesServiceArgoNeuT() = default;

    void reconfigure(fhicl::ParameterSet const& pset);
    virtual const detinfo::LArProperties*
    provider() const override
    {
      return this;
    }

    //------------------------------------------------------------------------
    //--- service provider interface

    /// Ar atomic number
    virtual double
    AtomicNumber() const override
    {
      return fZ;
    }
    /// Ar atomic mass (g/mol)
    virtual double
    AtomicMass() const override
    {
      return fA;
    }
    /// Ar mean excitation energy (eV)
    virtual double
    ExcitationEnergy() const override
    {
      return fI;
    }

    virtual double
    RadiationLength() const override
    {
      return fRadiationLength;
    } ///< g/cm^2
    virtual double
    Argon39DecayRate() const override
    {
      return fArgon39DecayRate;
    } ///< decays per cm^3 per second

    virtual double
    ScintResolutionScale() const override
    {
      return fScintResolutionScale;
    }
    virtual double
    ScintFastTimeConst() const override
    {
      return fScintFastTimeConst;
    }
    virtual double
    ScintSlowTimeConst() const override
    {
      return fScintSlowTimeConst;
    }
    virtual double
    ScintBirksConstant() const override
    {
      return fScintBirksConstant;
    }

    virtual bool
    ScintByParticleType() const override
    {
      return fScintByParticleType;
    }

    virtual double
    ScintYield(bool prescale = false) const override
    {
      return fScintYield * ScintPreScale(prescale);
    }
    virtual double
    ScintPreScale(bool prescale = true) const override
    {
      return (prescale ? fScintPreScale : 1);
    }
    virtual double
    ScintYieldRatio() const override
    {
      return fScintYieldRatio;
    }

    virtual double
    ProtonScintYield(bool prescale = false) const override
    {
      return fProtonScintYield * ScintPreScale(prescale);
    }
    virtual double
    ProtonScintYieldRatio() const override
    {
      return fProtonScintYieldRatio;
    }
    virtual double
    MuonScintYield(bool prescale = false) const override
    {
      return fMuonScintYield * ScintPreScale(prescale);
    }
    virtual double
    MuonScintYieldRatio() const override
    {
      return fMuonScintYieldRatio;
    }
    virtual double
    KaonScintYield(bool prescale = false) const override
    {
      return fKaonScintYield * ScintPreScale(prescale);
    }
    virtual double
    KaonScintYieldRatio() const override
    {
      return fKaonScintYieldRatio;
    }
    virtual double
    PionScintYield(bool prescale = false) const override
    {
      return fPionScintYield * ScintPreScale(prescale);
    }
    virtual double
    PionScintYieldRatio() const override
    {
      return fPionScintYieldRatio;
    }
    virtual double
    ElectronScintYield(bool prescale = false) const override
    {
      return fElectronScintYield * ScintPreScale(prescale);
    }
    virtual double
    ElectronScintYieldRatio() const override
    {
      return fElectronScintYieldRatio;
    }
    virtual double
    AlphaScintYield(bool prescale = false) const override
    {
      return fAlphaScintYield * ScintPreScale(prescale);
    }
    virtual double
    AlphaScintYieldRatio() const override
    {
      return fAlphaScintYieldRatio;
    }

    virtual bool
    CerenkovLightEnabled() const override
    {
      return fEnableCerenkovLight;
    }

    virtual std::map<double, double> SlowScintSpectrum() const override;
    virtual std::map<double, double> FastScintSpectrum() const override;
    virtual std::map<double, double> RIndexSpectrum() const override;
    virtual std::map<double, double> AbsLengthSpectrum() const override;
    virtual std::map<double, double> RayleighSpectrum() const override;

    virtual std::map<std::string, std::map<double, double>> SurfaceReflectances() const override;
    virtual std::map<std::string, std::map<double, double>> SurfaceReflectanceDiffuseFractions()
      const override;

    //------------------------------------------------------------------------

    // this stuff is moved to DetectorProperties
    double DriftVelocity(double efield = 0., double temperature = 0.) const; ///< cm/us
    double Efield(unsigned int planegap = 0) const;                          ///< kV/cm
    double ElectronLifetime() const;                                         ///< microseconds
    double Density(double temperature = 0.) const;                           ///< g/cm^3
    double Temperature() const;                                              ///< kelvin

    /// Restricted mean dE/dx energy loss (MeV/cm).
    double Eloss(double mom, double mass, double tcut) const;

    /// Energy loss fluctuation (sigma_E^2 / length in MeV^2/cm).
    double ElossVar(double mom, double mass) const;

    /// dQ/dX in electrons/cm, returns dE/dX in MeV/cm.
    double BirksCorrection(double dQdX) const;
    double BirksCorrection(double dQdX, double EField) const;
    double ModBoxCorrection(double dQdX) const;
    double ModBoxCorrection(double dQdX, double EField) const;
    virtual bool
    ExtraMatProperties() const override
    {
      return fExtraMatProperties;
    }
    virtual double
    TpbTimeConstant() const override
    {
      return fTpbTimeConstant;
    }

  private:
    void preBeginRun(art::Run const& run);
    void checkDBstatus() const;

    std::vector<double> fEfield; ///< kV/cm
    double fTemperature;         ///< kelvin
    double fElectronlifetime;    ///< microseconds
    double fDefTemperature;      ///< kelvin
    double fDefElectronlifetime; ///< microseconds
    double fRadiationLength;     ///< g/cm^2

    double fArgon39DecayRate; ///<  decays per cm^3 per second

    bool
      fAlreadyReadFromDB; ///< tests whether the values have alread been picked up from the Database

    // Following parameters are for use in Bethe-Bloch formula for dE/dx.

    double fZ;     ///< Ar atomic number
    double fA;     ///< Ar atomic mass (g/mol)
    double fI;     ///< Ar mean excitation energy (eV)
    double fSa;    ///< Sternheimer parameter a
    double fSk;    ///< Sternheimer parameter k
    double fSx0;   ///< Sternheimer parameter x0
    double fSx1;   ///< Sternheimer parameter x1
    double fScbar; ///< Sternheimer parameter Cbar

    // Optical parameters for LAr

    std::vector<double> fFastScintSpectrum;
    std::vector<double> fFastScintEnergies;
    std::vector<double> fSlowScintSpectrum;
    std::vector<double> fSlowScintEnergies;
    std::vector<double> fRIndexSpectrum;
    std::vector<double> fRIndexEnergies;
    std::vector<double> fAbsLengthSpectrum;
    std::vector<double> fAbsLengthEnergies;
    std::vector<double> fRayleighSpectrum;
    std::vector<double> fRayleighEnergies;

    bool fExtraMatProperties;
    virtual std::map<double, double> TpbAbs() const override;
    virtual std::map<double, double> TpbEm() const override;
    double fTpbTimeConstant;

    bool fScintByParticleType;

    double fProtonScintYield;
    double fProtonScintYieldRatio;
    double fMuonScintYield;
    double fMuonScintYieldRatio;
    double fPionScintYield;
    double fPionScintYieldRatio;
    double fKaonScintYield;
    double fKaonScintYieldRatio;
    double fElectronScintYield;
    double fElectronScintYieldRatio;
    double fAlphaScintYield;
    double fAlphaScintYieldRatio;

    double fScintYield;
    double fScintPreScale;
    double fScintResolutionScale;
    double fScintFastTimeConst;
    double fScintSlowTimeConst;
    double fScintYieldRatio;
    double fScintBirksConstant;

    bool fEnableCerenkovLight;

    std::vector<std::string> fReflectiveSurfaceNames;
    std::vector<double> fReflectiveSurfaceEnergies;
    std::vector<std::vector<double>> fReflectiveSurfaceReflectances;
    std::vector<std::vector<double>> fReflectiveSurfaceDiffuseFractions;

    std::vector<double> fTpbEmmisionEnergies;
    std::vector<double> fTpbEmmisionSpectrum;
    std::vector<double> fTpbAbsorptionEnergies;
    std::vector<double> fTpbAbsorptionSpectrum;

    struct DBsettingsClass {
      DBsettingsClass();

      bool ToughErrorTreatment; ///< equivalent parameter in DatabaseUtil
      bool ShouldConnect;       ///< equivalent parameter in DatabaseUtil
    };                          // DBsettingsClass

    DBsettingsClass DBsettings; ///< settings read from DB access

  }; // class LArPropertiesServiceArgoNeuT

  /// type of provider name following LArSoft name convention
  using LArPropertiesArgoNeuT = LArPropertiesServiceArgoNeuT;

} //namespace util
DECLARE_ART_SERVICE_INTERFACE_IMPL(util::LArPropertiesServiceArgoNeuT,
                                   detinfo::LArPropertiesService,
                                   LEGACY)
#endif // LARPROPERTIESSERVICEARGONEUT_H
