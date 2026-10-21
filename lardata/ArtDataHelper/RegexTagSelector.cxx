#include "lardata/ArtDataHelper/RegexTagSelector.h"

#include "canvas/Persistency/Provenance/BranchDescription.h"

#include <utility>

namespace lar::util::detail {
  RegexTagSelector::RegexTagSelector(std::string label, std::string instance, std::string process)
    : fLabelPattern{std::move(label)}
    , fInstancePattern{std::move(instance)}
    , fProcessPattern{std::move(process)}
    , fLabelRe{fLabelPattern}
    , fInstanceRe{fInstancePattern}
    , fProcessRe{fProcessPattern}
  {}

  bool RegexTagSelector::doMatch(art::BranchDescription const& p) const
  {
    return std::regex_match(p.moduleLabel(), fLabelRe) &&
           std::regex_match(p.productInstanceName(), fInstanceRe) &&
           std::regex_match(p.processName(), fProcessRe);
  }

  std::string RegexTagSelector::doPrint(std::string const& indent) const
  {
    return indent + "Module label matching regex: \"" + fLabelPattern + "\"\n" + indent +
           "Instance name matching regex: \"" + fInstancePattern + "\"\n" + indent +
           "Process name matching regex: \"" + fProcessPattern + "\"";
  }
}
