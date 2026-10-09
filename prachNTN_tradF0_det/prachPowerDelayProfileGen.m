function delayProfile = prachPowerDelayProfileGen(timeData, prachConfig, simConfig)
nRx = simConfig.nRx;
Nrep = prachConfig.Nrep;
NIFFT = prachConfig.bwpInfo.NFFT;

fdRootSeq = prachConfig.fdRootSeq;
fdRootSeq_ = prachConfig.fdRootSeq_;
tdRootSeq = ifft(fdRootSeq,NIFFT)*sqrt(NIFFT);
tdRootSeq_ = ifft(fdRootSeq_,NIFFT)*sqrt(NIFFT);

agc = comm.AGC('DesiredOutputPower',1, 'AveragingLength',1);
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
pdp_temp = zeros(NIFFT*Nrep,1);
pdp_temp_ = zeros(NIFFT*Nrep,1);




for nrx = 1:nRx
    corrIn = timeData_agc(:,nrx);
    [c,lag] = xcorr(corrIn,tdRootSeq_agc.');
    startIndex = find(lag==0);
    pdp_temp = pdp_temp + c(startIndex:startIndex+(NIFFT*Nrep)-1);
    
    

    if simConfig.simAlgo
        [c_,lag_] = xcorr(corrIn,tdRootSeq_sym_agc.');
        c_ = fftshift(c_);
        startIndex = find(lag_==0);
        pdp_temp_ = pdp_temp_ + c_(startIndex:startIndex+(NIFFT*Nrep)-1);
    end
end





pdp_temp = circshift(pdp_temp,-2);
pdp_temp_ = circshift(pdp_temp_,-2);

if strcmp(simConfig.channelType,'AWGN')
    
    if simConfig.simAlgo
        pdp = abs(pdp_temp).^2/59422300;
        pdp_ = abs(pdp_temp_).^2/60577000;

        delayProfile = [pdp pdp_];
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