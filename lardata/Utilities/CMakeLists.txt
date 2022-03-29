
art_make(NO_PLUGINS
         LIB_LIBRARIES
            lardata_Utilities_LArFFT_service
            lardataobj_RecoBase
            larcorealg_Geometry
              art::Framework_Core
              art::Framework_Principal
              art::Persistency_Provenance
            canvas::canvas
            messagefacility::MF_MessageLogger
            cetlib_except::cetlib_except
            FFTW3::FFTW3
          ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(DatabaseUtil "service"
              messagefacility::MF_MessageLogger
              PostgreSQL::PostgreSQL
              art::Framework_Core
              art::Framework_Principal
              art_root_io::RootDB
              art::Persistency_Provenance) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(LArFFT "service"
              lardataalg_DetectorInfo
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::FFTW
              ROOT::Hist
              ROOT::MathCore
              ROOT::Physics) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(ComputePi "module"
              messagefacility::MF_MessageLogger) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

include(FindOpenMP)
if(OPENMP_FOUND)
  # even if OpenMP is found on a SLF6 machine, it cannot be used.
  # see the discussion in https://cdcvs.fnal.gov/redmine/issues/23001
  if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^2.6.32.*" )
    message(STATUS "OpenMP is not supported on SLF6")
  else()
    add_definitions(-DWITH_OPENMP) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - avoid directory-scope functions: use target_compile_definitions() or target_compile_features() whenever possible
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif(OPENMP_FOUND)

## Problematic with cetbuildtools 8 and does not appear to be used
##add_library (gshf_MarqFitAlg STATIC MarqFitAlg.cxx)
##set_property(TARGET gshf_MarqFitAlg PROPERTY POSITION_INDEPENDENT_CODE ON)

install_headers()
install_fhicl()
install_source()

add_subdirectory(test)
