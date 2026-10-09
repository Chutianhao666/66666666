clc;
clear;
close all;

% 地面铁路信道下，评估不同阈值因子对PRACH检测性能的影响
% 必须从本目录运行，避免调用到 PRACH 中的同名三输入 PDP 函数。
rng(20261010, 'twister');
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
    'CarrierFrequency', 1.9e9, ...       % 实验假设：1.9 GHz，不是强制规范参数
    'MobileSpeedKmh', 350, ...          % 高铁速度，单位 km/h
    'DelayProfile', 'TDL-D', ...        % 地面视距配置；不是完整铁路专用信道
    'DelaySpread', 100e-9, ...          % 实验假设：100 ns RMS 时延扩展
    'NoiseReferencePowerDbw', 0, ...    % 延续原 awgn 的 1 W 名义功率约定
    'TransmissionDelay', 0, ...         % 非负整数，单位为原始波形样点
    'c', 299792458);

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

            validateattributes(simConfig.TransmissionDelay, {'numeric'}, ...
                {'scalar','integer','nonnegative'});
            txSig = [zeros(simConfig.TransmissionDelay,simConfig.nTx); txSig];
            

            % pass through wireless channel
            if strcmp(simConfig.channelType, 'TDL')
                speedMps = simConfig.MobileSpeedKmh/3.6;
                mobileMaxDoppler = speedMps*simConfig.CarrierFrequency/simConfig.c;
                % fmax 表示信道的 Doppler 尺度，不是可直接用作估频真值的公共 CFO。
                railwayChannel = nrTDLChannel;
                railwayChannel.DelayProfile = simConfig.DelayProfile;
                railwayChannel.DelaySpread = simConfig.DelaySpread;
                railwayChannel.SampleRate = prachConfig.bwpInfo.SamplingRate;
                railwayChannel.MaximumDopplerShift = mobileMaxDoppler;
                railwayChannel.MIMOCorrelation = 'Low';
                railwayChannel.NumTransmitAntennas = simConfig.nTx;
                railwayChannel.NumReceiveAntennas = simConfig.nRx;
                % 使用全局随机流；每次 trial 获得不同实现，同时可复现整次实验。
                railwayChannel.RandomStream = 'Global stream';
                cleanRx = railwayChannel(txSig);
                release(railwayChannel);



            elseif strcmp(simConfig.channelType, 'AWGN')
                cleanRx = repmat(sum(txSig,2),1,simConfig.nRx);
            else
                error('Unsupported channelType: %s',simConfig.channelType);
            end

            % H1 和 H0 使用相同名义噪声功率，分别抽取独立噪声。
            % SNR 仍沿用原工程的名义输入约定，不是自动测量的有效突发 SNR。
            rxSig = awgn(cleanRx,currentSNR,simConfig.NoiseReferencePowerDbw);
            rxNoiseSig = awgn(zeros(size(cleanRx)),currentSNR,simConfig.NoiseReferencePowerDbw);


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





