function delayProfile = prachPowerDelayProfileGen(timeData, prachConfig, simConfig)
nRx = simConfig.nRx;
Nrep = prachConfig.Nrep;
NIFFT = prachConfig.bwpInfo.NFFT;
TBHead = prachConfig.pssRootseq;

fdRootSeq = prachConfig.fdRootSeq;               % 频域ZC序列
fdRootSeq_ = prachConfig.fdRootSeq_;             % 频域互补序列
tdRootSeq = ifft(fdRootSeq,NIFFT)*sqrt(NIFFT);   % 时域序列
tdRootSeq_ = ifft(fdRootSeq_,NIFFT)*sqrt(NIFFT); % 时域互补序列
TBHeadseq = ifft(TBHead, NIFFT)*sqrt(NIFFT);

agc = comm.AGC('DesiredOutputPower',1, 'AveragingLength',1); % 自动增益控制
tdRootSeq_agc = agc(tdRootSeq.');
tdRootSeq_sym_agc = agc(tdRootSeq_.');

    if nRx > 1
        timeData_agc = zeros(NIFFT*Nrep,nRx);
        
        for nrx = 1:nRx
            timeData_agc(:,nrx) = agc(timeData(nrx,:).');
        end
    else
        timeData_agc = agc(timeData.');
    end



% figure; subplot(211); plot(abs(timeData_agc(:,1))); subplot(212); plot(abs(timeData_agc(:,2)));
% figure; plot(abs(tdRootSeq_agc));
pdp_temp = zeros(NIFFT*Nrep,1);  % 主序列互相关结果
pdp_temp_ = zeros(NIFFT*Nrep,1); % 互补序列互相关结果
pdp_tempTb = zeros(NIFFT*Nrep,1);

%% 同步头检测
if simConfig.simAlgo == 2
    for nrx = 1:nRx
        corrIn = timeData_agc(:,nrx);
        % figure;plot(abs(TBHeadseq));
        [TBc,lag] = xcorr(corrIn,TBHeadseq);
        maxCorr = max(abs(TBc));
        % figure;plot(abs(TBc))
        TBc = fftshift(TBc);
        startIndex = find(lag == 0);
        pdp_temp = pdp_tempTb + TBc(startIndex:startIndex+(NIFFT*Nrep)-1);    
    end
end
%% 共轭序列检测
if simConfig.simAlgo == 1
    for nrx = 1:nRx
        corrIn = timeData_agc(:,nrx); % 当前天线接收信号    
        [c,lag] = xcorr(corrIn,tdRootSeq_agc.'); % 与主序列互相关
        % figure;plot(abs(c));
        startIndex = find(lag==0); % 找到零延迟位置
        pdp_temp = pdp_temp + c(startIndex:startIndex+(NIFFT*Nrep)-1); % 累加多天线结果
        % figure;plot(abs(pdp_temp));
        [c_,lag_] = xcorr(corrIn,tdRootSeq_sym_agc.'); % 与互补序列互相关
        c_ = fftshift(c_); % 频域中心对齐
        startIndex = find(lag_==0);
        pdp_temp_ = pdp_temp_ + c_(startIndex:startIndex+(NIFFT*Nrep)-1);
        % figure;plot(abs(pdp_temp_));
    end
end

%% 普通检测
if simConfig.simAlgo == 0
    for nrx = 1:nRx
        corrIn = timeData_agc(:,nrx); % 当前天线接收信号    
        [c,lag] = xcorr(corrIn,tdRootSeq_agc.'); % 与主序列互相关
        % figure;plot(abs(c));
        startIndex = find(lag==0); % 找到零延迟位置
        pdp_temp = pdp_temp + c(startIndex:startIndex+(NIFFT*Nrep)-1); % 累加多天线结果
        % figure;plot(abs(pdp_temp));
    end
end
%%




pdp_temp = circshift(pdp_temp,-2); % 循环左移2样本
pdp_temp_ = circshift(pdp_temp_,-2); % 对齐

if strcmp(simConfig.channelType,'AWGN')
    
    if simConfig.simAlgo
        pdp = abs(pdp_temp).^2/59422300; % AWGN主序列功率归一化
        pdp_ = abs(pdp_temp_).^2/60577000; % 互补序列功率归一化
        delayProfile = [pdp pdp_]; % 合并结果
    else
        pdp = abs(pdp_temp).^2/66148500;
        delayProfile = pdp;
    end
    
elseif strcmp(simConfig.channelType,'TDL')
    
    if simConfig.simAlgo
        pdp = abs(pdp_temp).^2/54458400;
        pdp_ = abs(pdp_temp_).^2/55480600;

        delayProfile = [pdp pdp_];
    else
        pdp = abs(pdp_temp).^2/60875200;
 
        delayProfile = pdp;
    end
    
end







end
