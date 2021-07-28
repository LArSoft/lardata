include_directories(${FFTW_INCLUDE_DIR})
cet_find_library(FFTW_LIBRARY NAMES fftw3 fftw3-3 PATHS $ENV{FFTW_DIR}/$ENV{FFTW_FQ}/lib )
set(FFTW_LIBRARIES ${FFTW_LIBRARY})

art_make(NO_PLUGINS
         LIB_LIBRARIES
            lardata_Utilities_LArFFT_service
            lardata_Utilities_LArPropertiesServiceArgoNeuT_service
            lardataobj_RecoBase
            larcorealg_Geometry
            canvas::canvas
            messagefacility::MF_MessageLogger
            cetlib_except::cetlib_except
            ${FFTW_LIBRARIES})

simple_plugin(DatabaseUtil "service"
              messagefacility::MF_MessageLogger
              ${PQ}
              art::Framework_Core
              art::Framework_Principal
              ${ART_ROOT_IO_ROOTDB}
              art::Persistency_Provenance)

simple_plugin(LArPropertiesServiceArgoNeuT "service"
              lardata_Utilities_DatabaseUtil_service
              messagefacility::MF_MessageLogger)

simple_plugin(DetectorPropertiesServiceArgoNeuT "service"
              lardata_Utilities_LArPropertiesServiceArgoNeuT_service
              lardata_Utilities
              lardata_Utilities_DatabaseUtil_service
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              ${ART_ROOT_IO_ROOTDB}
              ${SQLITE3}
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::RIO)

simple_plugin(LArFFT "service"
              lardataalg_DetectorInfo
              messagefacility::MF_MessageLogger
              ROOT::Core
              ROOT::FFTW
              ROOT::Hist
              ROOT::MathCore
              ROOT::Physics)

simple_plugin(ComputePi "module"
              messagefacility::MF_MessageLogger)

include(FindOpenMP)
if(OPENMP_FOUND)
  # even if OpenMP is found on a SLF6 machine, it cannot be used.
  # see the discussion in https://cdcvs.fnal.gov/redmine/issues/23001
  if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^2.6.32.*" )
    message(STATUS "OpenMP is not supported on SLF6")
  else()
    add_definitions(-DWITH_OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif(OPENMP_FOUND)

add_library (gshf_MarqFitAlg STATIC MarqFitAlg.cxx)
set_property(TARGET gshf_MarqFitAlg PROPERTY POSITION_INDEPENDENT_CODE ON)

install_headers()
install_fhicl()
install_source()

add_subdirectory(test)
