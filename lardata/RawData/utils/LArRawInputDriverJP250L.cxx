////////////////////////////////////////////////////////////////////////
/// \file  LArRawInputDriverJP250L.cxx
/// \brief Source to convert JP250L files to LArSoft files
///
/// \author  eito@post.kek.jp, brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "lardata/RawData/utils/LArRawInputDriverJP250L.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "lardataobj/RawData/DAQHeader.h"
#include "lardataobj/RawData/RawDigit.h"

#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "canvas/Persistency/Provenance/FileFormatVersion.h"
#include "canvas/Persistency/Provenance/Timestamp.h"

#include "TFile.h"
#include "TObject.h"
#include "TTree.h"

#include <memory>
#include <vector>

// ======================================================================
// JP250L conversion of raw DAQ file to ART format

namespace lris {
  // ======================================================================
  // class c'tor/d'tor:
  LArRawInputDriverJP250L::LArRawInputDriverJP250L(fhicl::ParameterSet const&, // Not used
                                                   art::ProductRegistryHelper& helper,
                                                   art::SourceHelper const& pm)
    : principalMaker_(pm), m_current(0), m_data(0)
  {
    helper.reconstitutes<raw::DAQHeader, art::InEvent>("daq");
    helper.reconstitutes<std::vector<raw::RawDigit>, art::InEvent>("daq");
    helper.reconstitutes<sumdata::RunData, art::InRun>("daq");
  }

  // ======================================================================
  void LArRawInputDriverJP250L::closeCurrentFile() { delete[] m_data; }

  // ======================================================================
  void LArRawInputDriverJP250L::readFile(std::string const& name, art::FileBlock*& fb)
  {
    TFile m_f(name.c_str());
    TTree* runTree = dynamic_cast<TTree*>(m_f.Get("runTree"));

    m_eventTree = dynamic_cast<TTree*>(m_f.Get("eventTree"));
    m_eventTree->SetDirectory(0);
    m_nEvent = m_eventTree->GetEntries();

    // run information
    runTree->SetBranchAddress("runID", &m_runID);
    runTree->SetBranchAddress("unixtime", &m_unixtime);
    runTree->SetBranchAddress("nChannels", &m_nChannels);
    runTree->SetBranchAddress("nSamples", &m_nSamples);
    runTree->GetEntry(0);

    // have to add 1 to the runID because it can't be zero
    // adding 1 to every runID ensures that all runIDs are bumped in the same way
    m_runID += 1;

    const int nLength = (m_nSamples + 4) * m_nChannels;
    m_data = new unsigned short[nLength];
    m_eventTree->SetBranchAddress("data", m_data);

    // Fill and return a new Fileblock.
    // The string tells you what the version of the LArRawInputDriver is
    fb = new art::FileBlock(art::FileFormatVersion(1, "LArRawInputJP250L 2013_01"), name);

    return;
  }

  // ======================================================================
  bool LArRawInputDriverJP250L::readNext(art::RunPrincipal* const& /* inR */,
                                         art::SubRunPrincipal* const& /* inSR */,
                                         art::RunPrincipal*& outR,
                                         art::SubRunPrincipal*& outSR,
                                         art::EventPrincipal*& outE)
  {

    if (m_current > m_nEvent) return false;

    m_eventTree->GetEntry(m_current);

    raw::DAQHeader daqHeader;
    daqHeader.SetRun(m_runID);
    daqHeader.SetTimeStamp(m_unixtime);
    daqHeader.SetNChannels(m_nChannels);
    // set the event number to start at 1 - ART likes to start numbering
    // things from 1
    daqHeader.SetEvent(m_current + 1);

    // the following DAQHeader fields are not used by JP250L at this time
    //daqHeader.SetStatus(1);
    //daqHeader.SetFixedWord(h1.fixed);
    //daqHeader.SetFileFormat(h1.format);
    //daqHeader.SetSoftwareVersion(h1.software);
    //daqHeader.SetSpareWord(h1.spare);

    // make unique_ptrs for the data products to store in the event
    std::unique_ptr<raw::DAQHeader> daqcol(new raw::DAQHeader(daqHeader));
    std::unique_ptr<std::vector<raw::RawDigit>> rdcol(new std::vector<raw::RawDigit>);

    // loop over the signals and break them into
    // one RawDigit for each channel
    std::vector<short> adcVec(m_nSamples, 0);
    for (unsigned int n = 0; n < m_nChannels; ++n) {
      adcVec.clear();
      for (unsigned int i = 0; i < m_nSamples; ++i) {
        adcVec.push_back(m_data[(m_nSamples + 4) * n + (4 + i)]);
      }
      rdcol->push_back(raw::RawDigit(n, m_nSamples, adcVec));
    }

    art::RunNumber_t rn = daqHeader.GetRun();
    art::SubRunNumber_t sn = 1;
    art::EventNumber_t en = daqHeader.GetEvent();
    art::Timestamp tstamp = daqHeader.GetTimeStamp();

    // Make the Run and SubRun principals
    // this step is done once per run.
    if (m_current < 1) {
      std::unique_ptr<sumdata::RunData> rundata(new sumdata::RunData("jpl250l"));
      outR = principalMaker_.makeRunPrincipal(rn, tstamp);
      outSR = principalMaker_.makeSubRunPrincipal(rn, sn, tstamp);
      art::put_product_in_principal(std::move(rundata), *outR, "daq");
    }

    // make the Event principal
    outE = principalMaker_.makeEventPrincipal(rn, sn, en, tstamp);

    // Put products in the event.
    // first argument places the desired data product in the file
    // second argument is event record to associate it to
    // third argument is the label of the module storing the
    // information.  "daq" is standard in LArSoft
    art::put_product_in_principal(std::move(rdcol), *outE,
                                  "daq"); // Module label
    art::put_product_in_principal(std::move(daqcol), *outE,
                                  "daq"); // Module label

    // increment the entry to get out of the TTree for the next time through the loop
    ++m_current;

    return true;
  }

}
