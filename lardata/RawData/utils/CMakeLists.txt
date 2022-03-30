cet_make_library(
         SOURCE LArRawInputDriver.cxx
                LArRawInputDriverJP250L.cxx
                LArRawInputDriverLongBo.cxx
                LArRawInputDriverShortBo.cxx
          LIBRARIES PUBLIC
                       lardataobj::RawData
                       larcorealg::Geometry
                       larcoreobj::SummaryData
                       art::Framework_Core
                       art::Framework_IO_Sources
                       art::Framework_Principal
                       art::Persistency_Common
                       art::Persistency_Provenance
                    PRIVATE
                       canvas::canvas
                       cetlib_except::cetlib_except
                       PostgreSQL::PostgreSQL
                       Boost::serialization
                       Boost::date_time
                       ROOT::Core
                       ROOT::RIO
                       ROOT::Tree
        )

cet_build_plugin(LArRawInputSource art::source
  LIBRARIES REG
  lardataobj::RawData
  lardata_RawData_utils
  art::Framework_IO_Sources
)

install_headers()
install_fhicl()
install_source()
