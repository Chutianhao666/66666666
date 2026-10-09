% 遍历预设的时频偏移范围，计算接收信号与本地参考信号的互相关峰值，找到最优补偿值

function [time, freq] = time_freq_refinement(prachConfig, rxWaveform)
fdRootSeq = prachConfig.fdRootSeq;
nIFFT = prachConfig.nFFT;
firstSC = prachConfig.firstSC;
L_RA = prachConfig.L_RA;
Nu = prachConfig.Nu; % 有效序列长度
NCP_RA = prachConfig.NCP_RA;

ifftin = zeros(1,nIFFT); % 初始化全零频域信号
ifftin(firstSC + (1:L_RA)) = fdRootSeq; % 将根序列映射到指定子载波
ifftout = ifft(fftshift(ifftin))*sqrt(nIFFT); % IFFT转换并能量归一化

fRange = -500:10:500;
tRange = -15:1:15;
peaks = zeros(length(fRange),length(tRange));

for fIdx = 1:length(fRange)
    f = fRange(fIdx);
    doppler = exp(-1i*2*pi*(f/(prachConfig.bwpInfo.SamplingRate).*(1:size(rxWaveform, 1))));
    rxComp1 = rxWaveform.*(doppler.');
    for tIdx = 1:length(tRange)
        t = tRange(tIdx);
        if t < 0
            rxComp2 = [rxComp1(-t+1:end); zeros(-t,1)];
        else
            rxComp2 = [zeros(t,1); rxComp1(1:end-t)];
        end

        rxComp2 = rxComp2(NCP_RA+1:NCP_RA+Nu);
        corr = xcorr(rxComp2, ifftout.');
        maxCorr = max(abs(corr));
        
        peaks(fIdx,tIdx) = maxCorr;
    end
end


figure; contour(tRange,fRange,peaks)
[fmax,tmax] = find(peaks==max(max(peaks)));

time = tRange(tmax);
freq = fRange(fmax);

hold on
scatter3(time,freq,peaks(fmax,tmax),200,'black','x','LineWidth',2.5);
scatter3(-12,-103,peaks(fmax,tmax),200,'red','.','LineWidth',2.5);


end
