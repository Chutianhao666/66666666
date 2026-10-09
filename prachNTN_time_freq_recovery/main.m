clc;
clear;
close all;

% 模拟卫星通信场景下的PRACH信号传输，考虑TDL多径和卫星多普勒效应
% 生成PRACH信号，模拟传输延迟，通过信道模型并添加噪声
% 前端处理提取时域数据，计算功率延迟分布，分布完成时频同步
% 统计时频估计误差，通过CDF和核密度图可视化分析同步精度

%% Parameter Configuration
simConfig = struct('nTx', 1, ...
    'nRx', 2, ...
    'puschBW', 100, ...
    'puschSCS', 30, ...
    'prachFormat', 'F0', ...
    'prachZeroConfigZone', 13, ...
    'prachRootSeqIdx', 22, ...
    'channelType', 'TDL', ...
    'snrdB', -10, ...
    'simEpoch', 1, ...
    'simAlgo', 1, ...
    'SatelliteAltitude', 600000, ...
    'MobileAltitude', 0, ...
    'CarrierFrequency', 2e9, ...
    'ElevationAngle', 50, ...
    'MobileSpeed', 1, ...
    'TransmissionDelay', randi([200,2000]), ...
    'c', physconst('lightspeed'));




%% Simulation Starts



currentSNR = simConfig.snrdB;


fErrors = [];
tErrors = [];

for epoch = 1:simConfig.simEpoch
    [txSig, prachConfig] = prachSigGen(simConfig); %生成时域PRACH信号及配置

    delay_samples = floor(simConfig.TransmissionDelay*1e-9*4096*30e3);
    delay_samples = 200;

    txSig = [zeros(delay_samples,simConfig.nTx); txSig(1:end-delay_samples,:)];



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


        % rng(1);
        tdlChannelInfo = info(ntnChannel);
        [tdlOut, tdlPathGains, tdlSampleTimes] = ntnChannel(txSig);

        rxSig = awgn(tdlOut,currentSNR);
        % rxSig = tdlOut;
        % figure; subplot(211); plot(abs(rxSig(:,1))); subplot(212); plot(abs(rxSig(:,2)));

        simConfig.SatelliteDopplerShift = satelliteDopplerShift;

    elseif strcmp(simConfig.channelType, 'AWGN')
        rxSig = txSig;
        for nrx = 1:simConfig.nRx
            rxSig(:,nrx) = sum(txSig,2);
        end

        rxSig = awgn(rxSig,currentSNR);
       % figure; subplot(211); plot(abs(rxSig(:,1))); subplot(212); plot(abs(rxSig(:,2)));

    end


    timeData = prachFrontendProcessing(rxSig, prachConfig, simConfig); % 前端处理（下采样、去CP、频偏补偿）

    delayProfile = prachPowerDelayProfileGen(timeData, prachConfig, simConfig); % 生成功率延迟分布


    [time_comp, freq_comp] = time_freq_recovery(delayProfile, prachConfig); % 时频粗同步

    % time and frequency compensation
    rxWaveform = rxSig(:,1); % 取第一天线信号
    rxWaveform_comp = [rxWaveform(time_comp+1:end); zeros(time_comp,1)]; % 时间补偿
    doppler_comp = exp(-1i*2*pi*(freq_comp/(prachConfig.bwpInfo.SamplingRate).*(1:size(rxWaveform, 1)))); % 频偏补偿
    rxWaveform_comp = rxWaveform_comp.*(doppler_comp.');

    [time_refine, freq_refine] = time_freq_refinement(prachConfig, rxWaveform_comp); % 进一步优化时频估计
    
    t = time_comp + time_refine; % 总时间偏移
    f = freq_comp + freq_refine; % 总频率偏移

    fErrors = [fErrors; f - satelliteDopplerShift];
    tErrors = [tErrors; t - delay_samples];

    

    % fErrors = [fErrors; recoveryResult.fError];
    % tErrors = [tErrors; recoveryResult.tError];
    

    % if epoch == simConfig.simEpoch 
    %     figure; plot(cdf());
    % end

    if epoch == simConfig.simEpoch
        figure; cdfplot(fErrors);title('fErrors');
        figure; cdfplot(tErrors);title('tErrors');
    end

end % end of epoch

[ffTime,xxTime] = ksdensity(tErrors);
cdfTime = cumsum(ffTime)/sum(ffTime);
figure; plot(xxTime,cdfTime,'lineWidth',1);title('CDFTime');

[ffFreq,xxFreq] = ksdensity(fErrors);
cdfFreq = cumsum(ffFreq)/sum(ffFreq);
figure; plot(xxFreq,cdfFreq,'lineWidth',1);title('CDFFreq');

save('tErrors.mat','tErrors');
save('fErrors.mat','fErrors');

figure; 
cdfplot(tErrors);
hold on;
plot(xxTime,cdfTime);
ylim([0,1])
hold off;

figure;
cdfplot(fErrors);
hold on;
plot(xxFreq,cdfFreq);
ylim([0,1])
hold off









