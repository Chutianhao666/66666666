function cfg = defaultConfig()
% 839-point circular sequence-domain research model, not a full NR waveform.
cfg.length = 839;
cfg.root = 1; % original logical index 22 maps to ACTUAL root 1
cfg.scsHz = 1250;
cfg.fs = cfg.length*cfg.scsHz;
cfg.nRx = 2;
cfg.carrierHz = 1.9e9;
cfg.lightSpeed = 299792458;
cfg.maxRangeM = 6000;
cfg.maxDelaySamples = ceil(2*cfg.maxRangeM/cfg.lightSpeed*cfg.fs);
cfg.trueMotion = struct('speedKmh',350,'direction',-1, ...
    'alongTrackM',3000,'trackOffsetM',50);
cfg.priorMotion = struct('speedKmh',330,'direction',-1, ...
    'alongTrackM',3000,'trackOffsetM',50);
cfg.oscillatorHz = 150; % simulation truth ONLY; never read by recovery
cfg.priorOscillatorHz = 0;
cfg.priorUncertaintyHz = 300; % prior error bound assumed, not measured truth
cfg.channelType = 'Rician'; % 'AWGN' or block-flat 'Rician'
cfg.ricianKDb = 10;
cfg.phaseLagSamples = floor(cfg.length/4);
cfg.phaseIterations = 2;
cfg.timingToleranceSamples = 1;
cfg.targetPfa = 0.01;
cfg.calibrationTrials = 500;
cfg.validationNoiseTrials = 500;
cfg.epochs = 200;
cfg.snrDb = [-20 -15 -10 -5 0];
cfg.seed = 20261009;
cfg.makePlots = true;
end
