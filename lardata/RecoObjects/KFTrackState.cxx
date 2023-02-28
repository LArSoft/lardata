#include "KFTrackState.h"

using namespace trkf;

bool KFTrackState::updateWithHitState(const HitState& hitstate)
{
  // if track and hit not on same plane do not update and return false
  if ((hitstate.plane().position() - fTrackState.plane().position()).Mag2() > 10e-6) return false;
  if ((hitstate.plane().direction() - fTrackState.plane().direction()).Mag2() > 10e-6) return false;
  // Kalman Filter update (simplified case: 1D measurement along the same coordinate as element 0 of the track parameters)
  SMatrixSym55 tmp;
  tmp(0, 0) = 1. / (hitstate.hitMeasErr2() + fTrackState.covariance()(0, 0));
  fTrackState.setParameters(fTrackState.parameters() +
                            fTrackState.covariance() * tmp.Col(0) *
                              (hitstate.hitMeas() - fTrackState.parameters()(0)));
  fTrackState.setCovariance(
    (fTrackState.covariance() - ROOT::Math::Similarity(fTrackState.covariance(), tmp)));
  return true;
}

bool KFTrackState::combineWithTrackState(const TrackState& trackstate)
{
  // if tracks not on same plane do not update and return false
  if ((trackstate.plane().position() - fTrackState.plane().position()).Mag2() > 10e-6) return false;
  if ((trackstate.plane().direction() - fTrackState.plane().direction()).Mag2() > 10e-6)
    return false;
  // compute the weighted average of the two states
  const SVector5& par1 = fTrackState.parameters();
  const SVector5& par2 = trackstate.parameters();
  const SMatrixSym55& cov1 = fTrackState.covariance();
  const SMatrixSym55& cov2 = trackstate.covariance();
  SMatrixSym55&& cov = cov1 + cov2;
  bool success = cov.Invert();
  if (!success) return false;
  SMatrix55 K = cov1 * cov;
  fTrackState.setParameters(par1 + K * (par2 - par1));
  K = K * cov2;
  fTrackState.setCovariance(K.LowerBlock());
  return true;
}
