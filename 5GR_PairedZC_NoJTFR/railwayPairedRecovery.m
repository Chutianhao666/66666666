function result = railwayPairedRecovery(rx, cfg, ref)
% Bounded modular peak disambiguation + phase CFO + ONE-dimensional retiming.
% No JTFR, no two-dimensional delay/frequency grid, no access to channel truth.
assert(size(rx,1)==cfg.length && size(rx,2)==cfg.nRx);
assert(cfg.priorUncertaintyHz>0 && cfg.phaseIterations>=1);
L = cfg.length;
n = (0:L-1).';
fPrior = railwayDoppler(cfg.priorMotion,cfg.carrierHz,cfg.lightSpeed) ...
    +cfg.priorOscillatorHz;
rPrior = rx.*exp(-1i*2*pi*fPrior*n/cfg.fs);
corrResult = pairedZCCorrelation(rPrior,ref);
[divisor, inverseRoot, ~] = gcd(ref.root,L);
assert(divisor==1);
inverseRoot = mod(inverseRoot,L);
% Enumerate only algebraically admissible integer aliases, not a frequency grid.
% For the default prior bound, aliasLimit=0 and this is a single delay curve.
aliasLimit = floor((cfg.priorUncertaintyHz+cfg.scsHz/2)/cfg.scsHz);
assert(aliasLimit<(L-1)/2,'Prior too wide to resolve integer CFO aliases.');
best = -Inf; d0 = NaN; integerCfo = NaN;
for m = -aliasLimit:aliasLimit
    for d = 0:cfg.maxDelaySamples
        ka = mod(d+inverseRoot*m,L);
        kb = mod(d-inverseRoot*m,L);
        coarseHz = m*cfg.scsHz;
        evidence = sqrt(corrResult.power(ka+1,1)*corrResult.power(kb+1,2));
        evidence = evidence/(1+(coarseHz/cfg.priorUncertaintyHz)^2);
        if evidence>best
            best = evidence; d0 = d; integerCfo = coarseHz;
        end
    end
end
result.valid = isfinite(d0);
result.score = 0; result.delaySamples = NaN; result.cfoHz = NaN;
result.priorHz = fPrior; result.coarseDelaySamples = d0;
result.phaseCoherence = NaN;
result.timingPower = zeros(cfg.maxDelaySamples+1,1);
if ~result.valid, return; end
frequency = fPrior+integerCfo;
delay = d0;
for iteration = 1:cfg.phaseIterations
    corrected = rx.*exp(-1i*2*pi*frequency*n/cfg.fs);
    delayedReference = circshift(ref.pair,delay);
    phase = phaseCFOEstimator(corrected,delayedReference,cfg.fs,cfg.phaseLagSamples);
    if ~phase.valid
        result.valid = false; return;
    end
    frequency = frequency+phase.cfoHz;
    corrected = rx.*exp(-1i*2*pi*frequency*n/cfg.fs);
    matched = ifft(bsxfun(@times,fft(corrected,[],1),conj(fft(ref.pair))),[],1);
    % Normalized noncoherent composite-reference match, bounded by 0..1.
    denom = sum(abs(ref.pair).^2)*sum(abs(rx(:)).^2);
    power = sum(abs(matched).^2,2)/max(denom,realmin);
    power = power(1:cfg.maxDelaySamples+1);
    [~, position] = max(power);
    delay = position-1;
end
% A prior-support rejection is observable behavior, never silently oracle-correct.
result.valid = abs(frequency-fPrior)<=cfg.priorUncertaintyHz;
if ~result.valid, return; end
result.delaySamples = delay;
result.cfoHz = frequency;
result.score = power(position);
result.phaseCoherence = phase.coherence;
result.timingPower = power;
end
