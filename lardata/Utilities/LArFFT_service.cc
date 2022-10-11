////////////////////////////////////////////////////////////////////////
//
// \file LArFFT_plugin
//
//  This class simplifies implementation of Fourier transforms.
//  Because all data inputs and outputs are purely real,  the
//  transforms implemented in this way get a substantial performance
//  increase ~2x.
//
// \author pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////

#include "lardata/Utilities/LArFFT.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

DEFINE_ART_SERVICE(util::LArFFT)
