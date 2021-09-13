#ifndef lardata_ArtDataHelper_ToElement_h
#define lardata_ArtDataHelper_ToElement_h

// 'to_element' is a small function object that can be presented to a
// transform algorithm, enabling iteration over the underlying
// reference instead of the Ptr (e.g.):
//
//   using lar::to_element;
//   std::vector<art::Ptr<recob::Hit>> hits = from_somewhere();
//   for (recob::Hit const& hit : hits | ranges::views::transform(to_element)) {
//     ...
//   }

#include "canvas/Persistency/Common/Ptr.h"

namespace lar {
  struct to_element_t {
    template <typename T>
    T const& operator()(art::Ptr<T> const& ptr) const {
      return *ptr;
    }
  };

  inline constexpr to_element_t to_element;
}

#endif // lardata_ArtDataHelper_ToElement_h
