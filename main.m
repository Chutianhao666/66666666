%%
%                    _ooOoo_
%                   o8888888o
%                   88" . "88
%                   (| -_- |)
%                   O\  =  /O
%                ____/`---'\____
%              .'  \\|     |//  `.
%             /  \\|||  :  |||//  \
%            /  _||||| -:- |||||-  \
%            |   | \\\  -  /// |   |
%            | \_|  ''\---/''  |   |
%            \  .-\__  `-`  ___/-. /
%          ___`. .'  /--.--\  `. . __
%       ."" '<  `.___\_<|>_/___.'  >'"".
%      | | :  `- \`.;`\ _ /`;.`/ - ` : | |
%      \  \ `-.   \_ __\ /__ _/   .-` /  /
% ======`-.____`-.___\_____/___.-`____.-'======
%                    `=---='
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%             佛祖保佑       永无BUG
%%
clc;
clear;
close all;

%% Parameter Configuration
% 配置仿真参数
simConfig = struct('nTx', 1, ...
    'nRx', 2, ...
    'puschBW', 100, ...
    'puschSCS', 30, ...
    'prachFormat', 'F0', ...
    'prachZeroConfigZone', 13, ...
    'prachRootSeqIdx', 22, ...
    'channelType', 'AWGN', ...  % 使用AWGN信道模型进行地面通信仿真
    'snrdB', -20, ...
    'distances', 3000, ... % 设置不同的传输距离
    'simEpoch', 1, ... 
    'simAlgo', 2, ...
    'CarrierFrequency', 3e9, ...  % 设置为3.6 GHz频段
    'MobileSpeed', 3400 , ...  % 移动设备速度为3400 m/s（10 Ma）
    'c', physconst('lightspeed')); % 光速


detection_snr = zeros(length(simConfig.distances),1); % 存储每个SNR和距离的检测结果

%% Simulation Starts

for distIdx = 1:length(simConfig.distances)
    % 设置当前传输距离
    currentDistance = simConfig.distances(distIdx);

    % 设置传输时延
    transmissionDelay = 2 * currentDistance / simConfig.c;
    % transmissionDelay = 2e-9;
    
    de = [];

        for epoch = 1:simConfig.simEpoch
            % 生成PRACH信号
            [txSig, prachConfig] = prachSigGen(simConfig);

           % 传播延迟
           fs = prachConfig.bwpInfo.SamplingRate; 
           numDelaySamples = round(transmissionDelay * fs);  
           
           txSig = [zeros(numDelaySamples, simConfig.nTx); txSig(1:end - numDelaySamples, :)];

           % pathLoss = (4*pi*currentDistance*simConfig.CarrierFrequency / simConfig.c) ^ 2;
           % txSigAtt = txSig / sqrt(pathLoss);
           
           if strcmp(simConfig.channelType, 'TDL')
                mobileMaxDoppler = simConfig.MobileSpeed*simConfig.CarrierFrequency/simConfig.c;

                ntnChannel = nrTDLChannel;
                ntnChannel.DelayProfile = 'TDL-E';
                ntnChannel.DelaySpread = 3e-9;
                ntnChannel.MaximumDopplerShift = mobileMaxDoppler;
                ntnChannel.MIMOCorrelation = 'Low';
                ntnChannel.NumTransmitAntennas = 1;
                ntnChannel.NumReceiveAntennas = 2;

                tdlChannelInfo = info(ntnChannel);
                [tdlOut, tdlPathGains, tdlSampleTimes] = ntnChannel(txSig);

                rxSig = awgn(tdlOut,simConfig.snrdB);
                
   
            % 信道模型：使用AWGN模型（地面通信）
           elseif strcmp(simConfig.channelType, 'AWGN')
                rxSig = txSig;
                for nrx = 1:simConfig.nRx
                    rxSig(:, nrx) = sum(txSig, 2); % 多天线接收处理
                end
                % rxSig = awgn(rxSig, simConfig.snrdB); % 添加噪声
           end

            % rxPrachPart = prachpssDetection(rxSig, prachConfig, simConfig);

            % PRACH信号处理
             timeData = prachFrontendProcessing(rxSig, prachConfig, simConfig);
            % timeData = prachFrontendProcessing( rxPrachPart, prachConfig, simConfig);
            delayProfile = prachPowerDelayProfileGen(timeData, prachConfig, simConfig);
            
            
            % PRACH信号检测
            detResults = prachDetection(0.5, simConfig.snrdB, delayProfile, prachConfig);
            de = [de; detResults.detection];

            if epoch == simConfig.simEpoch
                detection_snr(distIdx) = sum(de) / simConfig.simEpoch;
                fprintf('Distance = %2.1f m,  PRACH_Detection_Rate: %1.3f\n', ...
                         currentDistance,  detection_snr(distIdx));
            end
           
            % if epoch == simConfig.simEpoch && distIdx == length(simConfig.distances)
            % figure; plot(simConfig.distances, detection_snr, 'r-*');
            % end
        end % end of epoch
end % end of distances

