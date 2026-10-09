clc;
clear;
close all;

% 评估PRACH检测算法在不同信噪比下的性能
% 生成PRACH信号并模拟信道传输
% 接收端处理信号并计算功率延迟分布
% 基于阈值检测判断前导码存在性
% 统计检测率并可视化结果
% 仅使用主序列

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
    'simEpoch', 5, ...
    'simAlgo', 0, ...
    'SatelliteAltitude', 600000, ...
    'MobileAltitude', 0, ...
    'CarrierFrequency', 2e9, ...
    'ElevationAngle', 50, ...
    'MobileSpeed', 1, ...
    'TransmissionDelay', 0, ...
    'c', physconst('lightspeed'));




%% Simulation Starts

detection_snr = zeros(length(simConfig.snrdB),1); % 初始化检测率存储数组

for snrIdx = 1:length(simConfig.snrdB)
    currentSNR = simConfig.snrdB(snrIdx);


    de = [];

    for epoch = 1:simConfig.simEpoch
        [txSig, prachConfig] = prachSigGen(simConfig);

        txSig = [zeros(simConfig.TransmissionDelay,simConfig.nTx); txSig(1:end-simConfig.TransmissionDelay,:)]; % 添加传输延迟


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

            ntnChannel.MIMOCorrelation = 'Low'; % 低MIMO相关性
            ntnChannel.NumTransmitAntennas = 1; % 单发
            ntnChannel.NumReceiveAntennas = 2; % 双发
            % ntnChannel.RandomStream = 'mt19937ar with seed';
            % ntnChannel.Seed = simConfig.seed;


            % rng(simConfig.seed);
            tdlChannelInfo = info(ntnChannel);
            [tdlOut, tdlPathGains, tdlSampleTimes] = ntnChannel(txSig);

            rxSig = awgn(tdlOut,currentSNR); % 添加AWGN噪声
            
            % rxSig = tdlOut;
            % figure; subplot(211); plot(abs(rxSig(:,1))); subplot(212); plot(abs(rxSig(:,2)));



        elseif strcmp(simConfig.channelType, 'AWGN')
            rxSig = txSig;
            for nrx = 1:simConfig.nRx 
                rxSig(:,nrx) = sum(txSig,2); % 合并发射天线信号
            end

            rxSig = awgn(rxSig,currentSNR);
            % figure; subplot(211); plot(abs(rxSig(:,1))); subplot(212); plot(abs(rxSig(:,2)));

        end


        timeData = prachFrontendProcessing(rxSig, prachConfig, simConfig);

        delayProfile = prachPowerDelayProfileGen(timeData, prachConfig, simConfig);


        detResults = prachDetection(0.5, currentSNR, delayProfile, prachConfig); % PRACH检测（阈值为0.5）


        de = [de; detResults.detection]; % 记录检测结果（1成功/0失败）

        if epoch == simConfig.simEpoch % 每完成一个SNR仿真，计算检测率
            detection_snr(snrIdx) = sum(de)/simConfig.simEpoch;
            fprintf('SNR = %2.1f dB, PRACH_Detection_Rate: %1.3f\n',...
                         currentSNR, detection_snr(snrIdx));
        end

        if epoch == simConfig.simEpoch && snrIdx == length(simConfig.snrdB) % 绘制检测率曲线
            figure; plot(simConfig.snrdB, detection_snr, 'r-*');
        end

    end % end of epoch
end % end of snr





