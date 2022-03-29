add_subdirectory(ProxyBase)

art_make(LIB_LIBRARIES lardataobj_RecoBase ROOT::Core) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_source()
