clc;
clear;
close all;

% 使用主序列+互补序列

%% Parameter Configuration
simConfig = struct('nTx', 1, ...
    'nRx', 2, ...
    'puschBW', 100, ...
    'puschSCS', 30, ...
    'prachFormat', 'F0', ...
    'prachZeroConfigZone', 13, ...
    'prachRootSeqIdx', 22, ...
    'channelType', 'TDL', ...
    'snrdB', -5:-1:-20, ...
    'simEpoch', 20, ...
    'simAlgo', 1, ...
    'SatelliteAltitude', 600000, ...
    'MobileAltitude', 0, ...
    'CarrierFrequency', 2e9, ...
    'ElevationAngle', 50, ...
    'MobileSpeed', 1, ...
    'TransmissionDelay', 0, ...
    'c', physconst('lightspeed'));




%% Simulation Starts

detection_snr = zeros(length(simConfig.snrdB),1);

for snrIdx = 1:length(simConfig.snrdB)
    currentSNR = simConfig.snrdB(snrIdx);


    de = [];

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

            rxSig = awgn(tdlOut,currentSNR);
            
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

        delayProfile = prachPowerDelayProfileGen(timeData, prachConfig, simConfig);


        detResults = prachDetection(0.5, currentSNR, delayProfile, prachConfig);


        de = [de; detResults.detection];

        if epoch == simConfig.simEpoch
            detection_snr(snrIdx) = sum(de)/simConfig.simEpoch;
            fprintf('SNR = %2.1f dB, PRACH_Detection_Rate: %1.3f\n',...
                         currentSNR, detection_snr(snrIdx));
        end

        if epoch == simConfig.simEpoch && snrIdx == length(simConfig.snrdB)
            figure; plot(simConfig.snrdB, detection_snr, 'r-*');
        end

    end % end of epoch
end % end of snr




