function results = main_5GR(cfg)
% Complete base-MATLAB entry. Defaults are diagnostic, not paper-level counts.
if nargin==0, cfg = defaultConfig(); end
rng(cfg.seed,'twister');
[s, ref] = pairedZCGenerate(cfg);
calibration = calibrateThreshold(cfg,ref);
falseAlarms = 0;
for trial = 1:cfg.validationNoiseTrials
    rx = (randn(cfg.length,cfg.nRx)+1i*randn(cfg.length,cfg.nRx))/sqrt(2);
    estimate = railwayPairedRecovery(rx,cfg,ref);
    det = prachDetection5GR(estimate,calibration.threshold);
    falseAlarms = falseAlarms+det.detected;
end
results.cfg = cfg;
results.threshold = calibration.threshold;
results.falseAlarms = falseAlarms;
results.noiseTrials = cfg.validationNoiseTrials;
results.pfa = falseAlarms/cfg.validationNoiseTrials;
results.pfaWilson95 = wilson(falseAlarms,cfg.validationNoiseTrials);
results.summary = zeros(numel(cfg.snrDb),7);
results.timeErrors = cell(numel(cfg.snrDb),1);
results.frequencyErrors = cell(numel(cfg.snrDb),1);
for k = 1:numel(cfg.snrDb)
    count = 0; acquired = 0; invalid = 0;
    dt = []; df = [];
    for trial = 1:cfg.epochs
        [rx, truth] = railwayChannel(s,cfg,cfg.snrDb(k));
        estimate = railwayPairedRecovery(rx,cfg,ref);
        invalid = invalid+~estimate.valid;
        det = prachDetection5GR(estimate,calibration.threshold);
        count = count+det.detected;
        if det.detected
            dt(end+1) = det.delaySamples-truth.delaySamples; %#ok<AGROW>
            df(end+1) = det.cfoHz-truth.cfoHz; %#ok<AGROW>
            acquired = acquired+(abs(dt(end))<=cfg.timingToleranceSamples);
        end
    end
    results.timeErrors{k} = dt/cfg.fs*1e6;
    results.frequencyErrors{k} = df;
    results.summary(k,:) = [cfg.snrDb(k),count/cfg.epochs,acquired/cfg.epochs, ...
        sqrt(mean((dt/cfg.fs*1e6).^2)),sqrt(mean(df.^2)),invalid/cfg.epochs,count];
end
fprintf('Threshold %.6f; held-out H0 false alarms %d/%d; Wilson95 [%.4f %.4f]\n', ...
    results.threshold,falseAlarms,cfg.validationNoiseTrials,results.pfaWilson95);
disp('SNR_dB Pd correctTimingRate conditionalTimingRMSE_us conditionalCfoRMSE_Hz invalidRate detectedCount');
disp(results.summary);
if cfg.makePlots
    figure;
    plot(cfg.snrDb,results.summary(:,2),'o-',cfg.snrDb,results.summary(:,3),'s-');
    xlabel('Nominal input per-sample SNR (dB)'); ylabel('Probability'); ylim([0 1]);
    legend('Detection','Detection with correct timing','Location','best'); grid on;
    figure;
    for k = 1:numel(cfg.snrDb)
        errors = sort(abs(results.frequencyErrors{k}));
        if isempty(errors), continue; end
        plot(errors,(1:numel(errors))/numel(errors),'DisplayName',sprintf('%g dB',cfg.snrDb(k))); hold on;
    end
    xlabel('Absolute CFO error (Hz), detected trials only'); ylabel('Empirical CDF'); legend('show'); grid on;
    figure;
    for k = 1:numel(cfg.snrDb)
        errors = sort(abs(results.timeErrors{k}));
        if isempty(errors), continue; end
        plot(errors,(1:numel(errors))/numel(errors),'DisplayName',sprintf('%g dB',cfg.snrDb(k))); hold on;
    end
    xlabel('Absolute timing error (us), detected trials only'); ylabel('Empirical CDF'); legend('show'); grid on;
end
end

function interval = wilson(successes, total)
assert(total>0);
z = 1.95996398454005; p = successes/total;
center = (p+z^2/(2*total))/(1+z^2/total);
half = z*sqrt(p*(1-p)/total+z^2/(4*total^2))/(1+z^2/total);
interval = [max(0,center-half),min(1,center+half)];
end
