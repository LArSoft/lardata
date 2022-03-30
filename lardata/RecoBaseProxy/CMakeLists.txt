add_subdirectory(ProxyBase)

cet_make_library( SOURCE Track.cxx
   LIBRARIES PUBLIC lardataobj::RecoBase
             PRIVATE ROOT::Core)

install_headers()
install_source()
