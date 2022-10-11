//////////////////////////////////////////////////////////////////////////////
// \version
//
// \brief Helper functions for MVAReader and MVAWriter wrappers
//
// \author robert.sulej@cern.ch
//
//////////////////////////////////////////////////////////////////////////////
#ifndef ANAB_MVAWRAPPERBASE_H
#define ANAB_MVAWRAPPERBASE_H

#include "canvas/Persistency/Common/Ptr.h"

#include "lardataobj/AnalysisBase/MVAOutput.h"

#include <cmath>
#include <functional>
#include <string>
#include <typeinfo>
#include <unordered_map>

namespace anab {

  /// Helper functions for MVAReader/Writer and FVecReader/Writer wrappers.
  class FVectorWrapperBase {
  public:
  protected:
    std::string getProductName(std::type_info const& ti) const;
    size_t getProductHash(std::type_info const& ti) const { return ti.hash_code(); }
  };

  /// Helper functions for MVAReader and MVAWriter wrappers.
  class MVAWrapperBase {
  public:
  protected:
    // all mva outputs in the feature vecor sum up to p=1

    template <class T, size_t N>
    std::array<float, N> pAccumulate(std::vector<art::Ptr<T>> const& items,
                                     std::vector<FeatureVector<N>> const& outs) const;

    template <class T, size_t N>
    std::array<float, N> pAccumulate(std::vector<art::Ptr<T>> const& items,
                                     std::vector<float> const& weights,
                                     std::vector<FeatureVector<N>> const& outs) const;

    template <class T, size_t N>
    std::array<float, N> pAccumulate(std::vector<art::Ptr<T>> const& items,
                                     std::function<float(T const&)> fweight,
                                     std::vector<FeatureVector<N>> const& outs) const;

    template <class T, size_t N>
    std::array<float, N> pAccumulate(std::vector<art::Ptr<T>> const& items,
                                     std::function<float(art::Ptr<T> const&)> fweight,
                                     std::vector<FeatureVector<N>> const& outs) const;

    // outputs in the feature vecor sum up to p=1 in groups:
    // - members of a group are tagged in the mask with the same non-negative number
    // - entries with negative tag in the mask are ignored

    template <class T, size_t N>
    std::array<float, N> pAccumulate(std::vector<art::Ptr<T>> const& items,
                                     std::vector<FeatureVector<N>> const& outs,
                                     std::array<char, N> const& mask) const;
  };

} // namespace anab

//----------------------------------------------------------------------------
// MVAReader functions.
//
template <class T, size_t N>
std::array<float, N> anab::MVAWrapperBase::pAccumulate(
  std::vector<art::Ptr<T>> const& items,
  std::vector<anab::FeatureVector<N>> const& outs) const
{
  std::array<double, N> acc;
  acc.fill(0);

  float pmin = 1.0e-6, pmax = 1.0 - pmin;
  float log_pmin = std::log(pmin), log_pmax = std::log(pmax);

  for (auto const& ptr : items) {
    auto const& vout = outs[ptr.key()];
    for (size_t i = 0; i < vout.size(); ++i) {
      float v;
      if (vout[i] < pmin)
        v = log_pmin;
      else if (vout[i] > pmax)
        v = log_pmax;
      else
        v = std::log(vout[i]);

      acc[i] += v;
    }
  }

  if (!items.empty()) {
    double totp = 0.0;
    for (size_t i = 0; i < N; ++i) {
      acc[i] = exp(acc[i] / items.size());
      totp += acc[i];
    }
    for (size_t i = 0; i < N; ++i) {
      acc[i] /= totp;
    }
  }
  else
    std::fill(acc.begin(), acc.end(), 1.0 / N);

  std::array<float, N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = acc[i];
  return result;
}
//----------------------------------------------------------------------------

template <class T, size_t N>
std::array<float, N> anab::MVAWrapperBase::pAccumulate(
  std::vector<art::Ptr<T>> const& items,
  std::vector<float> const& weights,
  std::vector<anab::FeatureVector<N>> const& outs) const
{
  std::array<double, N> acc;
  acc.fill(0);

  float pmin = 1.0e-6, pmax = 1.0 - pmin;
  float log_pmin = std::log(pmin), log_pmax = std::log(pmax);
  double totw = 0.0;

  for (size_t k = 0; k < items.size(); ++k) {
    auto const& ptr = items[k];
    float w = weights[k];

    if (w == 0) continue;

    auto const& vout = outs[ptr.key()];
    for (size_t i = 0; i < vout.size(); ++i) {
      float v;
      if (vout[i] < pmin)
        v = log_pmin;
      else if (vout[i] > pmax)
        v = log_pmax;
      else
        v = std::log(vout[i]);

      acc[i] += w * v;
    }
    totw += w;
  }

  if (!items.empty()) {
    double totp = 0.0;
    for (size_t i = 0; i < N; ++i) {
      acc[i] = exp(acc[i] / totw);
      totp += acc[i];
    }
    for (size_t i = 0; i < N; ++i) {
      acc[i] /= totp;
    }
  }
  else
    std::fill(acc.begin(), acc.end(), 1.0 / N);

  std::array<float, N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = acc[i];
  return result;
}
//----------------------------------------------------------------------------

template <class T, size_t N>
std::array<float, N> anab::MVAWrapperBase::pAccumulate(
  std::vector<art::Ptr<T>> const& items,
  std::function<float(T const&)> fweight,
  std::vector<anab::FeatureVector<N>> const& outs) const
{
  std::array<double, N> acc;
  acc.fill(0);

  float pmin = 1.0e-6, pmax = 1.0 - pmin;
  float log_pmin = std::log(pmin), log_pmax = std::log(pmax);
  double totw = 0.0;

  for (size_t k = 0; k < items.size(); ++k) {
    auto const& ptr = items[k];
    float w = fweight(*ptr);

    if (w == 0) continue;

    auto const& vout = outs[ptr.key()];
    for (size_t i = 0; i < vout.size(); ++i) {
      float v;
      if (vout[i] < pmin)
        v = log_pmin;
      else if (vout[i] > pmax)
        v = log_pmax;
      else
        v = std::log(vout[i]);

      acc[i] += w * v;
    }
    totw += w;
  }

  if (!items.empty()) {
    double totp = 0.0;
    for (size_t i = 0; i < N; ++i) {
      acc[i] = exp(acc[i] / totw);
      totp += acc[i];
    }
    for (size_t i = 0; i < N; ++i) {
      acc[i] /= totp;
    }
  }
  else
    std::fill(acc.begin(), acc.end(), 1.0 / N);

  std::array<float, N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = acc[i];
  return result;
}
//----------------------------------------------------------------------------

template <class T, size_t N>
std::array<float, N> anab::MVAWrapperBase::pAccumulate(
  std::vector<art::Ptr<T>> const& items,
  std::function<float(art::Ptr<T> const&)> fweight,
  std::vector<anab::FeatureVector<N>> const& outs) const
{
  std::array<double, N> acc;
  acc.fill(0);

  float pmin = 1.0e-6, pmax = 1.0 - pmin;
  float log_pmin = std::log(pmin), log_pmax = std::log(pmax);
  double totw = 0.0;

  for (size_t k = 0; k < items.size(); ++k) {
    auto const& ptr = items[k];
    float w = fweight(ptr);

    if (w == 0) continue;

    auto const& vout = outs[ptr.key()];
    for (size_t i = 0; i < vout.size(); ++i) {
      float v;
      if (vout[i] < pmin)
        v = log_pmin;
      else if (vout[i] > pmax)
        v = log_pmax;
      else
        v = std::log(vout[i]);

      acc[i] += w * v;
    }
    totw += w;
  }

  if (!items.empty()) {
    double totp = 0.0;
    for (size_t i = 0; i < N; ++i) {
      acc[i] = exp(acc[i] / totw);
      totp += acc[i];
    }
    for (size_t i = 0; i < N; ++i) {
      acc[i] /= totp;
    }
  }
  else
    std::fill(acc.begin(), acc.end(), 1.0 / N);

  std::array<float, N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = acc[i];
  return result;
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// functions with the mask tagging groups og labels
//----------------------------------------------------------------------------

template <class T, size_t N>
std::array<float, N> anab::MVAWrapperBase::pAccumulate(
  std::vector<art::Ptr<T>> const& items,
  std::vector<anab::FeatureVector<N>> const& outs,
  std::array<char, N> const& mask) const
{
  size_t n_groups = 0;
  std::unordered_map<char, size_t> label2group;
  std::vector<size_t> nb_entries;
  std::array<int, N> groupidx;
  for (size_t i = 0; i < N; ++i) {
    int idx = -1;
    if (mask[i] >= 0) {
      auto search = label2group.find(mask[i]);
      if (search == label2group.end()) {
        idx = n_groups;
        label2group[mask[i]] = idx;
        nb_entries.push_back(0);
        ++n_groups;
      }
      else {
        idx = search->second;
        nb_entries[idx]++;
      }
    }
    groupidx[i] = idx;
  }

  std::array<double, N> acc;
  acc.fill(0);

  float pmin = 1.0e-6, pmax = 1.0 - pmin;
  float log_pmin = std::log(pmin), log_pmax = std::log(pmax);

  for (auto const& ptr : items) {
    auto const& vout = outs[ptr.key()];
    for (size_t i = 0; i < vout.size(); ++i) {
      if (groupidx[i] < 0) continue;

      float v;
      if (vout[i] < pmin)
        v = log_pmin;
      else if (vout[i] > pmax)
        v = log_pmax;
      else
        v = std::log(vout[i]);

      acc[i] += v;
    }
  }

  if (!items.empty()) {
    std::vector<double> totp(n_groups, 0.0);
    for (size_t i = 0; i < N; ++i) {
      if (groupidx[i] >= 0) {
        acc[i] = exp(acc[i] / items.size());
        totp[groupidx[i]] += acc[i];
      }
    }
    for (size_t i = 0; i < N; ++i) {
      if (groupidx[i] >= 0) { acc[i] /= totp[groupidx[i]]; }
    }
  }
  else {
    for (size_t i = 0; i < N; ++i) {
      if (groupidx[i] >= 0) { acc[i] = 1 / nb_entries[groupidx[i]]; }
    }
  }

  std::array<float, N> result;
  for (size_t i = 0; i < N; ++i)
    result[i] = acc[i];
  return result;
}
//----------------------------------------------------------------------------

#endif //ANAB_MVAWRAPPERBASE
