function det = prachDetection5GR(recovery, threshold)
assert(isfinite(threshold) && threshold>=0);
det.detected = recovery.valid && recovery.score>threshold;
det.delaySamples = NaN; det.cfoHz = NaN;
if det.detected
    det.delaySamples = recovery.delaySamples;
    det.cfoHz = recovery.cfoHz;
end
end
