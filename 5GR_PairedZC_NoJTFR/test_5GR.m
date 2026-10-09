function test_5GR()
cfg = defaultConfig();
[s, ref] = pairedZCGenerate(cfg);
assert(abs(mean(abs(s).^2)-1)<1e-12);
assert(max(abs(ref.minus-conj(ref.plus)))<1e-12);
reverse = mod(-(0:cfg.length-1),cfg.length)+1;
assert(max(abs(ref.fdMinus-conj(ref.fdPlus(reverse))))<1e-9);
n = (0:cfg.length-1).';
prior = railwayDoppler(cfg.priorMotion,cfg.carrierHz,cfg.lightSpeed);
cases = 0;
for root = [1 23 129]
    cfg.root = root; [s, ref] = pairedZCGenerate(cfg);
    for delay = [0 1 21 cfg.maxDelaySamples]
        for residual = [-250 -50 0 50 250]
            % Opposite antenna phases test noncoherent combining explicitly.
            r = circshift(s,delay).*exp(1i*2*pi*(prior+residual)*n/cfg.fs);
            rx = r*[1 -1];
            estimate = railwayPairedRecovery(rx,cfg,ref);
            assert(estimate.valid && estimate.delaySamples==delay,'Timing recovery failed.');
            assert(abs(estimate.cfoHz-prior-residual)<1e-6,'CFO sign/scale failed.');
            assert(estimate.score>1-1e-9,'Composite reference match failed.');
            cases = cases+1;
        end
    end
end
% Near-half-integer CFO: global peak pairing alone can fail without noise.
cfg.priorUncertaintyHz = 700;
for root = [1 23 129]
    cfg.root = root; [s, ref] = pairedZCGenerate(cfg);
    for residual = cfg.scsHz*[-0.51 -0.49 0.49 0.51]
        delay = cfg.maxDelaySamples;
        r = circshift(s,delay).*exp(1i*2*pi*(prior+residual)*n/cfg.fs);
        estimate = railwayPairedRecovery(r*[1 -1],cfg,ref);
        assert(estimate.valid && estimate.delaySamples==delay,'Half-integer pairing failed.');
        assert(abs(estimate.cfoHz-prior-residual)<1e-6);
        cases = cases+1;
    end
end
cfg.priorUncertaintyHz = 300;
cfg.root = 1; [s, ref] = pairedZCGenerate(cfg);
estimate = railwayPairedRecovery(complex(zeros(cfg.length,cfg.nRx)),cfg,ref);
det = prachDetection5GR(estimate,0.01);
assert(~det.detected && isnan(det.delaySamples));
% Known-reference phase estimator independently handles both CFO signs.
for residual = [-500 500]
    rx = s.*exp(1i*2*pi*residual*n/cfg.fs);
    phase = phaseCFOEstimator(rx,s,cfg.fs,cfg.phaseLagSamples);
    assert(abs(phase.cfoHz-residual)<1e-6);
end
rng(7,'twister');
cfg.channelType = 'AWGN';
[rx, truth] = railwayChannel(s,cfg,Inf);
estimate = railwayPairedRecovery(rx,cfg,ref);
assert(estimate.valid && estimate.delaySamples==truth.delaySamples);
assert(abs(estimate.cfoHz-truth.cfoHz)<1e-6);
fprintf('Passed %d paired delay/CFO cases, phase-sign, empty-input and railway-channel checks.\n',cases);
end
