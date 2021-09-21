//
// Name:  LArPropTest_module.cc
//
// Purpose: LArPropTest module.  Test some features of LArProperties
//          service.
//
// Created:  5-Apr-2012  H. Greenlee

#include <iomanip>
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataalg/DetectorInfo/LArProperties.h"

namespace util {
  class LArPropTest : public art::EDAnalyzer {
  public:
    explicit LArPropTest(fhicl::ParameterSet const& pset);

  private:
    void beginJob() override;
    void analyze(const art::Event& evt) override;
  };

  DEFINE_ART_MODULE(LArPropTest)

  LArPropTest::LArPropTest(const fhicl::ParameterSet& pset) : EDAnalyzer(pset) {}

  void
  LArPropTest::beginJob()
  {
    // Make sure assert is enabled.

    bool assert_flag = false;
    assert((assert_flag = true, assert_flag));
    if (!assert_flag) {
      std::cerr << "Assert is disabled" << std::endl;
      abort();
    }

    // Get services.

    detinfo::LArProperties const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();
    auto const detprop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();

    // Test (default) accessors.

    std::cout << "Density = " << detprop.Density() << " g/cm^3" << std::endl;
    std::cout << "Drift velocity = " << detprop.DriftVelocity() << " cm/usec" << std::endl;
    std::cout << "Efield = " << detprop.Efield() << " kV/cm" << std::endl;
    std::cout << "Temperature = " << detprop.Temperature() << " Kelvin" << std::endl;
    std::cout << "Electron lifetime = " << detprop.ElectronLifetime() << " usec" << std::endl;
    std::cout << "Radiation Length = " << larprop->RadiationLength() << " g/cm^2" << std::endl;
    std::cout << "Radiation Length = " << larprop->RadiationLength() / detprop.Density() << " cm"
              << std::endl;

    // Make sure default values are acting correctly.

    assert(detprop.Density() == detprop.Density(detprop.Temperature()));
    assert(detprop.Density() != detprop.Density(detprop.Temperature() + 0.1));
    assert(detprop.DriftVelocity() == detprop.DriftVelocity(detprop.Efield()));
    assert(detprop.DriftVelocity() ==
           detprop.DriftVelocity(detprop.Efield(), detprop.Temperature()));

    // Drift velocity vs. electric field.

    std::cout << "\nDrift Velocity vs. Electric Field.\n"
              << "      E (kV/cm)      v (cm/us)" << std::fixed << std::endl;
    for (int i = 0; i < 3; ++i) {
      double e = 0.5;
      if (i == 1)
        e = 0.666667;
      else if (i == 2)
        e = 0.8;
      double v = detprop.DriftVelocity(e);
      std::cout << std::setprecision(3) << std::setw(15) << e << std::setprecision(4)
                << std::setw(15) << v << std::endl;
    }

    // dE/dx.

    std::cout
      << "\nCompare "
         "http://pdg.lbl.gov/2011/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.dat\n"
      << std::endl;

    double mass = 0.10565839; // Muon.
    double fact[16] = {
      1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0};

    std::ios_base::fmtflags flags = std::cout.flags();
    std::cout << "     T         p     Ionization dE/dx|_R\n"
              << "   [MeV]    [MeV/c]  ---[MeV cm^2/g]----\n"
              << std::scientific << std::setprecision(3);

    // Loop over kinetic energy.

    for (double tbase = 1.; tbase <= 1.e9; tbase *= 10.) {
      for (int i = 0; i < 16; ++i) {

        // T in MeV.

        double t = tbase * fact[i];
        if (t > 1.e9) break;

        // Calculate momentum in GeV/C

        double p = std::sqrt(1.e-6 * t * t + 2.e-3 * t * mass);

        // Calculate restricted and unrestricted dE/dx.

        double dedxr = detprop.Eloss(p, mass, 0.05) / detprop.Density(); // Restricted.
        double dedx = detprop.Eloss(p, mass, 0.) / detprop.Density();    // Unrestricted.
        std::cout << std::setw(10) << t << std::setw(10) << 1000. * p << std::setw(10) << dedx
                  << std::setw(10) << dedxr << std::endl;
      }
    }
    std::cout.flags(flags);
  }

  void
  LArPropTest::analyze(const art::Event& /* evt */)
  {}
}
