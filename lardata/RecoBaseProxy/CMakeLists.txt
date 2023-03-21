cet_make_library(SOURCE Track.cxx
  LIBRARIES
  PUBLIC
  lardataobj::RecoBase
  larcorealg::geo_vectors_utils
  larcorealg::CoreUtils
  canvas::canvas
)

install_headers(SUBDIRS ProxyBase)
install_source(SUBDIRS ProxyBase)
