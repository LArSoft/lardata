
#ifndef LARDATA_LARDATA_ARTDATAHELPER_REGEXTAGMATCHER_H
#define LARDATA_LARDATA_ARTDATAHELPER_REGEXTAGMATCHER_H

#include "canvas/Utilities/InputTag.h"

namespace lar::util::detail {
  class RegexTagMatcher {
  public:
    explicit RegexTagMatcher(art::InputTag const&);

    bool match(art::InputTag const&) const;

  private:
    art::InputTag fRegexTag;
  };
}

#endif // LARDATA_LARDATA_ARTDATAHELPER_REGEXTAGMATCHER_H
