
art_make(NO_PLUGINS
         LIB_LIBRARIES
            lardataobj_RecoBase
            larcorealg_Geometry
              art::Framework_Core
              art::Framework_Principal
              art::Persistency_Provenance
            canvas::canvas
            messagefacility::MF_MessageLogger
            cetlib_except::cetlib_except
            FFTW_LIBRARY)

simple_plugin(LArPropertiesServiceArgoNeuT "service"
              lardata_Utilities_DatabaseUtil_service
              messagefacility::MF_MessageLogger)

simple_plugin(DetectorPropertiesServiceArgoNeuT "service"
              lardata_Utilities_argoneut_LArPropertiesServiceArgoNeuT_service
              lardata_Utilities
              lardata_Utilities_DatabaseUtil_service
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              ${ART_ROOT_IO_ROOTDB}
              ${SQLITE3}
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::RIO)


## Problematic with cetbuildtools 8 and does not appear to be used
##add_library (gshf_MarqFitAlg STATIC MarqFitAlg.cxx)
##set_property(TARGET gshf_MarqFitAlg PROPERTY POSITION_INDEPENDENT_CODE ON)

install_headers()
install_fhicl()
install_source()

add_subdirectory(test)
