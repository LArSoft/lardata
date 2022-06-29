cet_make_library(LIBRARY_NAME AssociationUtil INTERFACE
  SOURCE AssociationUtil.h
  LIBRARIES INTERFACE
  art::Framework_Principal
  art::Persistency_Common
  canvas::canvas
  messagefacility::MF_MessageLogger
)

cet_make_library(LIBRARY_NAME ForEachAssociatedGroup INTERFACE
  SOURCE ForEachAssociatedGroup.h
  LIBRARIES INTERFACE
  canvas::AssnsAlgorithms
)

cet_make_library(LIBRARY_NAME SimpleFits INTERFACE
  SOURCE SimpleFits.h
  LIBRARIES INTERFACE
  lardataalg::UtilitiesHeaders
)

cet_make_library(SOURCE
  GeometryUtilities.cxx
  LArFFTW.cxx
  LArFFTWPlan.cxx
  PxHitConverter.cxx
  Range.cxx
  SignalShaper.cxx
  SignalShaping.cxx
  UniqueRangeSet.cxx
  UtilException.cxx
  LIBRARIES PUBLIC
  lardata_Utilities_LArFFT_service
  lardataobj::RecoBase
  larcorealg::Geometry
  larvecutils::MarqFitAlg
  larcore::headers
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
  PostgreSQL::PostgreSQL
)

cet_build_plugin(LArFFT art::service
  LIBRARIES
  PUBLIC
  lardataalg::DetectorInfo
  larcore::headers
  PRIVATE
  messagefacility::MF_MessageLogger
  ROOT::Core
  ROOT::FFTW
  ROOT::Hist
  ROOT::MathCore
  ROOT::Physics
)

cet_build_plugin(ComputePi art::module
  LIBRARIES REG
  messagefacility::MF_MessageLogger
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(test)
