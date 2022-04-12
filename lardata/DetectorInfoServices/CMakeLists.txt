cet_make_library(LIBRARY_NAME ServicePack INTERFACE
  SOURCE ServicePack.h
  LIBRARIES INTERFACE
  larcore::ServiceUtil
  larcorealg::headers
  )

# Service interfaces.
set(builder_dest Modules)
function(make_simple_service_builder STEM)
  # Define an interface library for services and modules that use this
  # service interface.
  cet_make_library(LIBRARY_NAME ${STEM}Service INTERFACE
    SOURCE ${STEM}Service.h
    LIBRARIES INTERFACE ${ARGN})
  # Make a plugin builder for services of this type.
  cet_make_plugin_builder(lar::${STEM}Service art::service ${builder_dest}
    LIBRARIES CONDITIONAL lardata::${STEM}Service)
endfunction()

make_simple_service_builder(DetectorClocks lardataalg::DetectorInfo)
make_simple_service_builder(DetectorProperties lardataalg::DetectorInfo)
make_simple_service_builder(LArProperties lardataalg::DetectorInfo)

cet_collect_plugin_builders(${builder_dest} lar::DetectorInfoServiceBuilders)

include(lar::DetectorInfoServiceBuilders)

cet_build_plugin(DetectorClocksServiceStandard lar::DetectorClocksService
  LIBRARIES PRIVATE
  art_root_io::RootDB
  art_root_io::detail
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  ROOT::Hist
  ROOT::RIO
)

cet_build_plugin(DetectorPropertiesServiceStandard lar::DetectorPropertiesService
  LIBRARIES
  PUBLIC
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  PRIVATE
  lardata::LArPropertiesService
  lardata::ServicePack
  art_root_io::RootDB
  messagefacility::MF_MessageLogger
  ROOT::Tree
  ROOT::RIO
)

cet_build_plugin(LArPropertiesServiceStandard lar::LArPropertiesService
  LIBRARIES PRIVATE art::Framework_Principal
)

install_headers()
install_fhicl()
install_source()
