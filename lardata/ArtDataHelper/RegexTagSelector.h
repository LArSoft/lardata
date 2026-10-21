
#ifndef LARDATA_LARDATA_ARTDATAHELPER_REGEXTAGSELECTOR_H
#define LARDATA_LARDATA_ARTDATAHELPER_REGEXTAGSELECTOR_H

#include "art/Framework/Principal/SelectorBase.h"

#include <regex>
#include <string>

namespace lar::util::detail {
  /// art selector accepting the products whose label, instance and process
  /// names each fully match the corresponding regex pattern.
  class RegexTagSelector : public art::SelectorBase {
  public:
    /// @throws std::regex_error if a pattern is not a valid regex.
    RegexTagSelector(std::string label, std::string instance, std::string process);

  private:
    bool doMatch(art::BranchDescription const&) const override;
    std::string doPrint(std::string const& indent) const override;

    std::string fLabelPattern;
    std::string fInstancePattern;
    std::string fProcessPattern;
    std::regex fLabelRe;
    std::regex fInstanceRe;
    std::regex fProcessRe;
  };
}

#endif // LARDATA_LARDATA_ARTDATAHELPER_REGEXTAGSELECTOR_H
