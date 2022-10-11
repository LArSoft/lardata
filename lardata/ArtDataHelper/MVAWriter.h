//////////////////////////////////////////////////////////////////////////////
// \version
//
// \brief Wrapper for saving MVA results into art::Event
//
// \author robert.sulej@cern.ch
//
//////////////////////////////////////////////////////////////////////////////
#ifndef ANAB_MVAWRITER_H
#define ANAB_MVAWRITER_H

#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

#include "lardata/ArtDataHelper/MVAWrapperBase.h"

namespace anab {

  /// Index to the MVA output / FeatureVector collection, used when result vectors are added or set.
  typedef size_t FVector_ID;
  typedef size_t MVAOutput_ID;

  template <size_t N>
  class FVectorWriter : public FVectorWrapperBase {
  public:
    /// Name provided to the constructor is used as an instance name for FVecDescription<N>
    /// and FeatureVector<N> (for which it is combined with the processed data product names).
    /// The name is used as an instance name for the FVecDescription data product
    /// which lets you to save multiple vector collections from a single art module.
    FVectorWriter(art::ProducesCollector& collector, const char* name = "")
      : fCollector(collector)
      , fInstanceName(name)
      , fIsDescriptionRegistered(false)
      , fDescriptions(nullptr)
    {}

    /// Register the collection of metadata type FVecDescription<N> (once for all data types
    /// for which vectors are saved) and the collection of FeatureVectors<N> (using data type name
    /// added to fInstanceName as instance name of the collection made for the type T).
    template <class T>
    void produces_using();

    /// Initialize container for FeatureVectors and, if not yet done, the container for
    /// metadata, then creates metadata for data products of type T. FeatureVector container
    /// is initialized to hold dataSize vectors (if dataSize > 0): use setOutput() to
    /// store values.
    /// Returns index of collection which should be used when saving actual output values.
    template <class T>
    FVector_ID initOutputs(std::string const& dataTag,
                           size_t dataSize,
                           std::vector<std::string> const& names = std::vector<std::string>(N, ""));

    template <class T>
    FVector_ID initOutputs(art::InputTag const& dataTag,
                           size_t dataSize,
                           std::vector<std::string> const& names = std::vector<std::string>(N, ""))
    {
      return initOutputs<T>(dataTag.encode(), dataSize, names);
    }

    void setVector(FVector_ID id, size_t key, std::array<float, N> const& values)
    {
      (*(fVectors[id]))[key] = values;
    }
    void setVector(FVector_ID id, size_t key, std::array<double, N> const& values)
    {
      (*(fVectors[id]))[key] = values;
    }
    void setVector(FVector_ID id, size_t key, std::vector<float> const& values)
    {
      (*(fVectors[id]))[key] = values;
    }
    void setVector(FVector_ID id, size_t key, std::vector<double> const& values)
    {
      (*(fVectors[id]))[key] = values;
    }

    /// Initialize container for FeatureVectors and, if not yet done, the container for
    /// metadata, then creates metadata for data products of type T. FeatureVector container
    /// is initialized as EMPTY and vectors should be added with addOutput() function.
    /// Returns index of collection which should be used when adding actual output values.
    template <class T>
    FVector_ID initOutputs(art::InputTag const& dataTag,
                           std::vector<std::string> const& names = std::vector<std::string>(N, ""))
    {
      return initOutputs<T>(dataTag.encode(), 0, names);
    }

    template <class T>
    FVector_ID initOutputs(std::vector<std::string> const& names = std::vector<std::string>(N, ""))
    {
      return initOutputs<T>(std::string(""), 0, names);
    }

    void addVector(FVector_ID id, std::array<float, N> const& values)
    {
      fVectors[id]->emplace_back(values);
    }
    void addVector(FVector_ID id, std::array<double, N> const& values)
    {
      fVectors[id]->emplace_back(values);
    }
    void addVector(FVector_ID id, std::vector<float> const& values)
    {
      fVectors[id]->emplace_back(values);
    }
    void addVector(FVector_ID id, std::vector<double> const& values)
    {
      fVectors[id]->emplace_back(values);
    }

    /// Set tag of associated data products in case it was not ready at the initialization time.
    void setDataTag(FVector_ID id, art::InputTag const& dataTag)
    {
      (*fDescriptions)[id].setDataTag(dataTag.encode());
    }

    /// Check consistency and save all the results in the event.
    void saveOutputs(art::Event& evt);

    /// Get the number of contained feature vectors.
    size_t size(FVector_ID id) const { return fVectors[id]->size(); }

    /// Get the length of a single feature vector.
    size_t length() const { return N; }

    /// Get copy of the feature vector for the type T, at index "key".
    template <class T>
    std::array<float, N> getVector(size_t key) const
    {
      std::array<float, N> vout;
      auto const& src = (*(fVectors[getProductID<T>()]))[key];
      for (size_t i = 0; i < N; ++i)
        vout[i] = src[i];
      return vout;
    }

    /// Get copy of the feature vector for the type T, idicated with art::Ptr::key().
    template <class T>
    std::array<float, N> getVector(art::Ptr<T> const& item) const
    {
      std::array<float, N> vout;
      auto const& src = (*(fVectors[getProductID<T>()]))[item.key()];
      for (size_t i = 0; i < N; ++i)
        vout[i] = src[i];
      return vout;
    }

    friend std::ostream& operator<<(std::ostream& o, FVectorWriter const& a)
    {
      o << "FVectorWriter for " << a.fInstanceName << ", " << N << " outputs";
      if (!a.fRegisteredDataTypes.empty()) {
        o << ", ready to write results made for:" << std::endl;
        for (auto const& n : a.fRegisteredDataTypes) {
          o << "\t" << n << std::endl;
        }
      }
      else {
        o << ", nothing registered for writing to the events" << std::endl;
      }
      return o;
    }

  protected:
    // Data collected for each event:
    template <class T>
    FVector_ID getProductID() const;

    std::vector<std::unique_ptr<std::vector<anab::FeatureVector<N>>>> fVectors;

  private:
    // Data initialized for the module life:
    art::ProducesCollector& fCollector;
    std::string fInstanceName;

    std::vector<std::string> fRegisteredDataTypes;
    bool fIsDescriptionRegistered;

    std::unordered_map<size_t, FVector_ID> fTypeHashToID;

    std::unique_ptr<std::vector<anab::FVecDescription<N>>> fDescriptions;
    void clearEventData()
    {
      fTypeHashToID.clear();
      fVectors.clear();
      fDescriptions.reset(nullptr);
    }

    /// Check if the the writer is configured to write results for data product type name.
    bool dataTypeRegistered(const std::string& dname) const;
    /// Check if the containers for results prepared for "tname" data type are ready.
    bool descriptionExists(const std::string& tname) const;
  };

  /// Helper for registering in the art::EDProducer all data products needed for
  /// N-output MVA results: keep MVADescriptions<N> for all types T in one collection
  /// while separate instance names are used for the MVA output value collections for
  /// each type T.
  /// Use one instance of this class per one MVA model, applied to one or more types.
  template <size_t N>
  class MVAWriter : public FVectorWriter<N>, public MVAWrapperBase {
  public:
    /// Name provided to the constructor is used as an instance name for MVADescription<N>
    /// and FeatureVector<N> (for which it is combined with the processed data product names).
    /// Good idea is to use the name as an indication of what MVA model was used on the data
    /// (like eg. "emtrack" for outputs from a model distinguishing EM from track-like hits
    /// and clusters). The name is used as an instance name for the MVADescription data product
    /// which lets you to save multiple MVA results from a single art module.
    MVAWriter(art::ProducesCollector& collector, const char* name = "")
      : FVectorWriter<N>(collector, name)
    {}

    void setOutput(FVector_ID id, size_t key, std::array<float, N> const& values)
    {
      FVectorWriter<N>::setVector(id, key, values);
    }
    void setOutput(FVector_ID id, size_t key, std::array<double, N> const& values)
    {
      FVectorWriter<N>::setVector(id, key, values);
    }
    void setOutput(FVector_ID id, size_t key, std::vector<float> const& values)
    {
      FVectorWriter<N>::setVector(id, key, values);
    }
    void setOutput(FVector_ID id, size_t key, std::vector<double> const& values)
    {
      FVectorWriter<N>::setVector(id, key, values);
    }

    void addOutput(FVector_ID id, std::array<float, N> const& values)
    {
      FVectorWriter<N>::addVector(id, values);
    }
    void addOutput(FVector_ID id, std::array<double, N> const& values)
    {
      FVectorWriter<N>::addVector(id, values);
    }
    void addOutput(FVector_ID id, std::vector<float> const& values)
    {
      FVectorWriter<N>::addVector(id, values);
    }
    void addOutput(FVector_ID id, std::vector<double> const& values)
    {
      FVectorWriter<N>::addVector(id, values);
    }

    /// Get MVA results accumulated over the vector of items (eg. over hits associated to a cluster).
    /// NOTE: MVA outputs for these items has to be added to the MVAWriter first!
    template <class T>
    std::array<float, N> getOutput(std::vector<art::Ptr<T>> const& items) const
    {
      return pAccumulate<T, N>(
        items, *(FVectorWriter<N>::fVectors[FVectorWriter<N>::template getProductID<T>()]));
    }

    /// Get MVA results accumulated with provided weights over the vector of items
    /// (eg. over clusters associated to a track, weighted by the cluster size; or
    /// over hits associated to a cluster, weighted by the hit area).
    /// NOTE: MVA outputs for these items has to be added to the MVAWriter first!
    template <class T>
    std::array<float, N> getOutput(std::vector<art::Ptr<T>> const& items,
                                   std::vector<float> const& weights) const
    {
      return pAccumulate<T, N>(
        items,
        weights,
        *(FVectorWriter<N>::fVectors[FVectorWriter<N>::template getProductID<T>()]));
    }

    /// Get MVA results accumulated with provided weighting function over the vector
    /// of items (eg. over clusters associated to a track, weighted by the cluster size;
    /// or over hits associated to a cluster, weighted by the hit area).
    /// NOTE: MVA outputs for these items has to be added to the MVAWriter first!
    template <class T>
    std::array<float, N> getOutput(std::vector<art::Ptr<T>> const& items,
                                   std::function<float(T const&)> fweight) const
    {
      return pAccumulate<T, N>(
        items,
        fweight,
        *(FVectorWriter<N>::fVectors[FVectorWriter<N>::template getProductID<T>()]));
    }

    template <class T>
    std::array<float, N> getOutput(std::vector<art::Ptr<T>> const& items,
                                   std::function<float(art::Ptr<T> const&)> fweight) const
    {
      return pAccumulate<T, N>(
        items,
        fweight,
        *(FVectorWriter<N>::fVectors[FVectorWriter<N>::template getProductID<T>()]));
    }

    /// Get copy of the MVA output vector for the type T, at index "key".
    template <class T>
    std::array<float, N> getOutput(size_t key) const
    {
      return FVectorWriter<N>::template getVector<T>(key);
    }

    /// Get copy of the MVA output vector for the type T, idicated with art::Ptr::key().
    template <class T>
    std::array<float, N> getOutput(art::Ptr<T> const& item) const
    {
      return FVectorWriter<N>::template getVector<T>(item);
    }
  };

} // namespace anab

//----------------------------------------------------------------------------
// FVectorWriter functions.
//
template <size_t N>
template <class T>
anab::FVector_ID anab::FVectorWriter<N>::getProductID() const
{
  auto const& ti = typeid(T);
  auto search = fTypeHashToID.find(getProductHash(ti));
  if (search != fTypeHashToID.end()) { return search->second; }
  else {
    throw cet::exception("FVectorWriter")
      << "Feature vectors not initialized for product " << getProductName(ti) << std::endl;
  }
}
//----------------------------------------------------------------------------

template <size_t N>
bool anab::FVectorWriter<N>::dataTypeRegistered(const std::string& dname) const
{
  for (auto const& s : fRegisteredDataTypes) {
    if (s == dname) { return true; }
  }
  return false;
}
//----------------------------------------------------------------------------

template <size_t N>
template <class T>
void anab::FVectorWriter<N>::produces_using()
{
  std::string dataName = getProductName(typeid(T));
  if (dataTypeRegistered(dataName)) {
    throw cet::exception("FVectorWriter")
      << "Type " << dataName << "was already registered." << std::endl;
  }

  if (!fIsDescriptionRegistered) {
    fCollector.produces<std::vector<anab::FVecDescription<N>>>(fInstanceName);
    fIsDescriptionRegistered = true;
  }

  fCollector.produces<std::vector<anab::FeatureVector<N>>>(fInstanceName + dataName);
  fRegisteredDataTypes.push_back(dataName);
}
//----------------------------------------------------------------------------

template <size_t N>
bool anab::FVectorWriter<N>::descriptionExists(const std::string& tname) const
{
  if (!fDescriptions) return false;

  std::string n = fInstanceName + tname;
  for (auto const& d : *fDescriptions) {
    if (d.outputInstance() == n) { return true; }
  }
  return false;
}
//----------------------------------------------------------------------------

template <size_t N>
template <class T>
anab::FVector_ID anab::FVectorWriter<N>::initOutputs(std::string const& dataTag,
                                                     size_t dataSize,
                                                     std::vector<std::string> const& names)
{
  size_t dataHash = getProductHash(typeid(T));
  std::string dataName = getProductName(typeid(T));

  if (!dataTypeRegistered(dataName)) {
    throw cet::exception("FVectorWriter")
      << "Type " << dataName << "not registered with produces_using() function." << std::endl;
  }

  if (!fDescriptions) { fDescriptions = std::make_unique<std::vector<anab::FVecDescription<N>>>(); }
  else if (descriptionExists(dataName)) {
    throw cet::exception("FVectorWriter")
      << "FVecDescription<N> already initialized for " << dataName << std::endl;
  }
  fDescriptions->emplace_back(dataTag, fInstanceName + dataName, names);

  fVectors.push_back(std::make_unique<std::vector<anab::FeatureVector<N>>>());
  anab::FVector_ID id = fVectors.size() - 1;
  fTypeHashToID[dataHash] = id;

  if (dataSize) { fVectors[id]->resize(dataSize, anab::FeatureVector<N>(0.0F)); }

  return id;
}
//----------------------------------------------------------------------------

template <size_t N>
void anab::FVectorWriter<N>::saveOutputs(art::Event& evt)
{
  for (auto const& n : fRegisteredDataTypes) {
    if (!descriptionExists(n)) {
      throw cet::exception("FVectorWriter")
        << "No FVecDescription<N> prepared for type " << n << std::endl;
    }
  }

  if (fVectors.size() != fDescriptions->size()) {
    throw cet::exception("FVectorWriter")
      << "FVecDescription<N> vector length not equal to the number of FeatureVector<N> vectors"
      << std::endl;
  }

  for (size_t i = 0; i < fVectors.size(); ++i) {
    auto const& outInstName = (*fDescriptions)[i].outputInstance();
    if ((*fDescriptions)[i].dataTag().empty()) {
      throw cet::exception("FVectorWriter")
        << "FVecDescription<N> reco data tag not set for " << outInstName << std::endl;
    }
    evt.put(std::move(fVectors[i]), outInstName);
  }
  evt.put(std::move(fDescriptions), fInstanceName);
  clearEventData();
}
//----------------------------------------------------------------------------

#endif //ANAB_MVAREADER
