
cet_make_library(
         SOURCE GeometryUtilities.cxx
                LArFFTW.cxx
                LArFFTWPlan.cxx
                PxHitConverter.cxx
                Range.cxx
                SignalShaping.cxx
                UniqueRangeSet.cxx
                UtilException.cxx
         LIBRARIES PUBLIC
            lardata_Utilities_LArFFT_service
            lardataobj::RecoBase
            larcorealg::Geometry
            larvecutils::MarqFitAlg
            art::Framework_Core
            art::Framework_Principal
            art::Persistency_Provenance
          PRIVATE
            canvas::canvas
            messagefacility::MF_MessageLogger
            cetlib_except::cetlib_except
            FFTW3::FFTW3
          )

cet_build_plugin(DatabaseUtil art::service
         LIBRARIES PUBLIC
              art::Persistency_Provenance
              art::Framework_Core
              art::Framework_Principal
              art_root_io::RootDB
           PRIVATE
              messagefacility::MF_MessageLogger
              PostgreSQL::PostgreSQL)

cet_build_plugin(LArFFT art::service
         LIBRARIES PUBLIC
              lardataalg::DetectorInfo
           PRIVATE
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::FFTW
              ROOT::Hist
              ROOT::MathCore
              ROOT::Physics)

cet_build_plugin(ComputePi art::module
                 LIBRARIES REG messagefacility::MF_MessageLogger)

install_headers()
install_fhicl()
install_source()

add_subdirectory(test)
