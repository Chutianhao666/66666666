clc;
clear;
close all;

% 评估不同阈值因子对PRACH检测性能的影响
%% Parameter Configuration
simConfig = struct('nTx', 1, ...
    'nRx', 2, ...
    'puschBW', 100, ...
    'puschSCS', 30, ...
    'prachFormat', 'F0', ...
    'prachZeroConfigZone', 13, ...
    'prachRootSeqIdx', 22, ...
    'channelType', 'TDL', ...
    'snrdB', -6:-4:-18, ...
    'simEpoch', 1000, ...
    'simAlgo', 0, ...
    'SatelliteAltitude', 600000, ...
    'MobileAltitude', 0, ...
    'CarrierFrequency', 2e9, ...
    'ElevationAngle', 50, ...
    'MobileSpeed', 1, ...
    'TransmissionDelay', 0, ...
    'c', physconst('lightspeed'));

%% 
thrFactor = 0.1:0.05:0.9;

falseAlarm = zeros(length(thrFactor),1); %存储虚警率
missDetection = zeros(length(thrFactor),1); % 存储漏检率



%% Simulation Starts
for thrIdx = 1:length(thrFactor) % 对每个阈值因子，独立统计不同SNR下的FA和MD
    currentFactor = thrFactor(thrIdx);

    false_alarm_snr = zeros(length(simConfig.snrdB),1);
    miss_detection_snr = zeros(length(simConfig.snrdB),1);
    
    for snrIdx = 1:length(simConfig.snrdB) % 针对每个SNR，执行多次独立仿真
        currentSNR = simConfig.snrdB(snrIdx);

        fa = [];
        md = [];

        for epoch = 1:simConfig.simEpoch
            [txSig, prachConfig] = prachSigGen(simConfig);

            txSig = [zeros(simConfig.TransmissionDelay,simConfig.nTx); txSig(1:end-simConfig.TransmissionDelay,:)];
            

            % pass through wireless channel
            if strcmp(simConfig.channelType, 'TDL')
                mobileMaxDoppler = simConfig.MobileSpeed*simConfig.CarrierFrequency/simConfig.c;
                satelliteDopplerShift = dopplerShiftCircularOrbit( ...
                    simConfig.ElevationAngle, simConfig.SatelliteAltitude, ...
                    simConfig.MobileAltitude, simConfig.CarrierFrequency);

                ntnChannel = nrTDLChannel;
                ntnChannel.DelayProfile = 'NTN-TDL-C';
                ntnChannel.DelaySpread = 5.1286e-9;
                ntnChannel.SampleRate = prachConfig.bwpInfo.SamplingRate;
                ntnChannel.MaximumDopplerShift = mobileMaxDoppler;
                ntnChannel.SatelliteDopplerShift = satelliteDopplerShift;

                ntnChannel.MIMOCorrelation = 'Low';
                ntnChannel.NumTransmitAntennas = 1;
                ntnChannel.NumReceiveAntennas = 2;
                % ntnChannel.RandomStream = 'mt19937ar with seed';
                % ntnChannel.Seed = simConfig.seed;


                % rng(simConfig.seed);
                tdlChannelInfo = info(ntnChannel);
                [tdlOut, tdlPathGains, tdlSampleTimes] = ntnChannel(txSig);

                noiseSig = zeros(size(tdlOut));

                rxSig = awgn(tdlOut,currentSNR);
                rxNoiseSig = awgn(noiseSig,currentSNR);
                % rxSig = tdlOut;
                % figure; subplot(211); plot(abs(rxSig(:,1))); subplot(212); plot(abs(rxSig(:,2)));



            elseif strcmp(simConfig.channelType, 'AWGN')
                rxSig = txSig;
                for nrx = 1:simConfig.nRx
                    rxSig(:,nrx) = sum(txSig,2);
                end

                rxSig = awgn(rxSig,currentSNR);
                % figure; subplot(211); plot(abs(rxSig(:,1))); subplot(212); plot(abs(rxSig(:,2)));

            end


            timeData = prachFrontendProcessing(rxSig, prachConfig, simConfig);
            noiseData = prachFrontendProcessing(rxNoiseSig, prachConfig, simConfig);

            delayProfile = prachPowerDelayProfileGen(currentSNR, timeData, prachConfig, simConfig);
            noiseProfile = prachPowerDelayProfileGen(currentSNR, noiseData, prachConfig, simConfig);


            detResults = prachThr4F0(currentFactor, currentSNR, delayProfile, prachConfig);
            falResults = fa_cal(currentFactor, currentSNR, noiseProfile, prachConfig);

            fa = [fa; falResults.false_alarm];
            md = [md; detResults.miss_detection];

            if epoch == simConfig.simEpoch
                false_alarm_snr(snrIdx) = sum(fa)/simConfig.simEpoch;
                miss_detection_snr(snrIdx) = sum(md)/simConfig.simEpoch;
            end

            % if epoch == simConfig.simEpoch && snrIdx == length(simConfig.snrdB)
            %     figure; plot(simConfig.snrdB, peak_snr, 'r-*');
            %     figure; plot(simConfig.snrdB, mean_snr, 'r-*');
            %     figure; plot(simConfig.snrdB, peak2ave_snr, 'r-*');
            % end

        end % end of epoch
    end % end of snr
    fprintf('Finish thrFac = %d\n', thrIdx);

    falseAlarm(thrIdx) = sum(false_alarm_snr)/length(simConfig.snrdB);
    missDetection(thrIdx) = sum(miss_detection_snr)/length(simConfig.snrdB);

    if epoch == simConfig.simEpoch && snrIdx == length(simConfig.snrdB) && thrIdx == length(thrFactor)
        figure; plot(thrFactor, falseAlarm, 'r-*');
        figure; plot(thrFactor, missDetection, 'r-*');
    end
end





