//////////////////////////////////////////////////////////////////////////////
// \version
//
// \brief Wrappers for accessing MVA results and associated data products
//
// \author robert.sulej@cern.ch
//
//////////////////////////////////////////////////////////////////////////////
#ifndef ANAB_MVAREADER_H
#define ANAB_MVAREADER_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "lardata/ArtDataHelper/MVAWrapperBase.h"

namespace anab {

  /// Helper for reading the reconstructed objects of type T together with associated
  /// N-ellement feature vectors with their metadata (this class is not a data product).
  template <class T, size_t N>
  class FVectorReader : public FVectorWrapperBase {
  public:
    /// Create the helper for feature vectors stored in the event evt with the provided input tag
    /// (the same tag which was used to save vectors with FVectorWriter class).
    /// Returns nullptr if data products not found in the event.
    static std::unique_ptr<FVectorReader> create(const art::Event& evt, const art::InputTag& tag)
    {
      bool success;
      std::unique_ptr<FVectorReader> ptr(new FVectorReader(evt, tag, success));
      if (success) { return ptr; }
      else {
        return nullptr;
      }
    }

    /// Create the wrapper for feature vectors stored in the event evt with the provided input tag
    /// (the same tag which was used to save vectors with FVectorWriter class).
    /// Throws exception if data products not found in the event.
    FVectorReader(const art::Event& evt, const art::InputTag& tag);

    /// Access data product at index "key".
    T const& item(size_t key) const { return (*fDataHandle)[key]; }
    std::vector<T> const& items() const { return *fDataHandle; }

    /// Access the vector of the feature vectors.
    std::vector<FeatureVector<N>> const& vectors() const { return *fVectors; }

    /// Access feature vector data at index "key".
    /// *** WOULD LIKE TO CHANGE TYPE OF FVEC DATA MEMBER TO std::array AND THEN ENABLE THIS FUNCTION ***
    //const std::array<float, N> & getVector(size_t key) const { return (*fVectors)[key].data(); }

    /// Get copy of the feature vector at index "key".
    std::array<float, N> getVector(size_t key) const
    {
      std::array<float, N> vout;
      for (size_t i = 0; i < N; ++i)
        vout[i] = (*fVectors)[key][i];
      return vout;
    }

    /// Get copy of the feature vector idicated with art::Ptr::key().
    std::array<float, N> getVector(art::Ptr<T> const& item) const { return getVector(item.key()); }

    /// Get the number of contained items (no. of data product objects equal to no. of feature vectors).
    size_t size() const { return fVectors->size(); }

    /// Get the length of a single feature vector.
    size_t length() const { return N; }

    /// Get the input tag (string representation) of data product used to calculate feature vectors.
    const std::string& dataTag() const { return fDescription->dataTag(); }

    /// Access the data product handle.
    const art::Handle<std::vector<T>>& dataHandle() const { return fDataHandle; }

    /// Meaning/name of the index'th column in the collection of feature vectors.
    const std::string& columnName(size_t index) const { return fDescription->outputName(index); }

    /// Index of column with given name, or -1 if name not found.
    int getIndex(const std::string& name) const { return fDescription->getIndex(name); }

    friend std::ostream& operator<<(std::ostream& o, FVectorReader const& a)
    {
      o << "FVectorReader:" << std::endl << *(a.fDescription) << std::endl;
      return o;
    }

  protected:
    /// Not-throwing constructor.
    FVectorReader(const art::Event& evt, const art::InputTag& tag, bool& success);

  private:
    FVecDescription<N> const* fDescription;
    std::vector<FeatureVector<N>> const* fVectors;
    art::Handle<std::vector<T>> fDataHandle;
  };

  /// Helper for reading the reconstructed objects of type T together with associated
  /// N-outputs MVA results with their metadata (this class is not a data product).
  template <class T, size_t N>
  class MVAReader : public FVectorReader<T, N>, public MVAWrapperBase {
  public:
    /// Create the wrapper for MVA data stored in the event evt with the provided input tag
    /// (the same tag which was used to save MVA results with MVAWriter class).
    /// Returns nullptr if data products not found in the event.
    static std::unique_ptr<MVAReader> create(const art::Event& evt, const art::InputTag& tag)
    {
      bool success;
      std::unique_ptr<MVAReader> ptr(new MVAReader(evt, tag, success));
      if (success) { return ptr; }
      else {
        return nullptr;
      }
    }

    /// Create the wrapper for MVA data stored in the event evt with the provided input tag
    /// (the same tag which was used to save MVA results with MVAWriter class).
    /// Throws exception if data products not found in the event.
    MVAReader(const art::Event& evt, const art::InputTag& tag) : FVectorReader<T, N>(evt, tag) {}

    /// Access the vector of the feature vectors.
    std::vector<FeatureVector<N>> const& outputs() const { return FVectorReader<T, N>::vectors(); }

    /// Get copy of the MVA output vector at index "key".
    std::array<float, N> getOutput(size_t key) const { return FVectorReader<T, N>::getVector(key); }

    /// Get copy of the MVA output vector idicated with art::Ptr::key().
    std::array<float, N> getOutput(art::Ptr<T> const& item) const
    {
      return FVectorReader<T, N>::getVector(item.key());
    }

    /// Get MVA results accumulated over the vector of items (eg. over hits associated to a cluster).
    std::array<float, N> getOutput(std::vector<art::Ptr<T>> const& items) const
    {
      return pAccumulate(items, FVectorReader<T, N>::vectors());
    }

    /// Get MVA results accumulated with provided weights over the vector of items
    /// (eg. over clusters associated to a track, weighted by the cluster size; or
    /// over hits associated to a cluster, weighted by the hit area).
    std::array<float, N> getOutput(std::vector<art::Ptr<T>> const& items,
                                   std::vector<float> const& weights) const
    {
      return pAccumulate(items, weights, FVectorReader<T, N>::vectors());
    }

    /// Get MVA results accumulated with provided weighting function over the vector
    /// of items (eg. over clusters associated to a track, weighted by the cluster size;
    /// or over hits associated to a cluster, weighted by the hit area).
    std::array<float, N> getOutput(std::vector<art::Ptr<T>> const& items,
                                   std::function<float(T const&)> fweight) const
    {
      return pAccumulate(items, fweight, FVectorReader<T, N>::vectors());
    }

    /// Meaning/name of the index'th column in the collection of MVA output vectors.
    const std::string& outputName(size_t index) const
    {
      return FVectorReader<T, N>::columnName(index);
    }

  private:
    /// Not-throwing constructor.
    MVAReader(const art::Event& evt, const art::InputTag& tag, bool& success)
      : FVectorReader<T, N>(evt, tag, success)
    {}
  };

} // namespace anab

//----------------------------------------------------------------------------
// FVectorReader functions.
//
template <class T, size_t N>
anab::FVectorReader<T, N>::FVectorReader(const art::Event& evt, const art::InputTag& tag)
  : fDescription(0)
{
  if (!N) { throw cet::exception("FVectorReader") << "Vector size should be > 0." << std::endl; }

  auto descriptionHandle = evt.getValidHandle<std::vector<anab::FVecDescription<N>>>(tag);

  // search for FVecDescription<N> produced for the type T, with the instance name from the tag
  std::string outputInstanceName = tag.instance() + getProductName(typeid(T));
  for (auto const& dscr : *descriptionHandle) {
    if (dscr.outputInstance() == outputInstanceName) {
      fDescription = &dscr;
      break;
    }
  }
  if (!fDescription) {
    throw cet::exception("FVectorReader")
      << "Vectors description not found for " << outputInstanceName << std::endl;
  }

  fVectors = &*(evt.getValidHandle<std::vector<FeatureVector<N>>>(
    art::InputTag(tag.label(), fDescription->outputInstance(), tag.process())));

  if (!evt.getByLabel(fDescription->dataTag(), fDataHandle)) {
    throw cet::exception("FVectorReader")
      << "Associated data product handle failed: " << *(fDataHandle.whyFailed()) << std::endl;
  }

  if (fVectors->size() != fDataHandle->size()) {
    throw cet::exception("FVectorReader")
      << "Feature vectors and data products sizes inconsistent: " << fVectors->size()
      << "!=" << fDataHandle->size() << std::endl;
  }
}
//----------------------------------------------------------------------------

template <class T, size_t N>
anab::FVectorReader<T, N>::FVectorReader(const art::Event& evt,
                                         const art::InputTag& tag,
                                         bool& success)
  : fDescription(0)
{
  success = false; // until all is done correctly

  if (!N) {
    std::cout << "FVectorReader: Vector size should be > 0." << std::endl;
    return;
  }

  art::Handle<std::vector<anab::FVecDescription<N>>> descriptionHandle;
  if (!evt.getByLabel(tag, descriptionHandle)) { return; }

  // search for FVecDescription<N> produced for the type T, with the instance name from the tag
  std::string outputInstanceName = tag.instance() + getProductName(typeid(T));
  for (auto const& dscr : *descriptionHandle) {
    if (dscr.outputInstance() == outputInstanceName) {
      fDescription = &dscr;
      break;
    }
  }
  if (!fDescription) {
    std::cout << "FVectorReader: Vectors description not found for " << outputInstanceName
              << std::endl;
    return;
  }

  fVectors = &*(evt.getValidHandle<std::vector<FeatureVector<N>>>(
    art::InputTag(tag.label(), fDescription->outputInstance(), tag.process())));

  if (!evt.getByLabel(fDescription->dataTag(), fDataHandle)) {
    std::cout << "FVectorReader: Associated data product handle failed: "
              << *(fDataHandle.whyFailed()) << std::endl;
    return;
  }

  if (fVectors->size() != fDataHandle->size()) {
    std::cout << "FVectorReader: Feature vectors and data products sizes inconsistent: "
              << fVectors->size() << "!=" << fDataHandle->size() << std::endl;
    return;
  }

  success = true; // ok, all data found in the event
}
//----------------------------------------------------------------------------

#endif //ANAB_MVAREADER
