function calibration = calibrateThreshold(cfg, ref)
% H0 calibration covers selection, phase estimation, retiming and prior gating.
assert(cfg.targetPfa>0 && cfg.targetPfa<1 && cfg.calibrationTrials>=1);
scores = zeros(cfg.calibrationTrials,1);
for trial = 1:cfg.calibrationTrials
    rx = (randn(cfg.length,cfg.nRx)+1i*randn(cfg.length,cfg.nRx))/sqrt(2);
    estimate = railwayPairedRecovery(rx,cfg,ref);
    scores(trial) = estimate.score;
end
ordered = sort(scores);
rank = ceil((numel(scores)+1)*(1-cfg.targetPfa));
assert(rank<=numel(scores),'More H0 calibration trials needed for target Pfa.');
calibration.threshold = ordered(rank);
calibration.scores = scores;
calibration.targetPfa = cfg.targetPfa;
end
