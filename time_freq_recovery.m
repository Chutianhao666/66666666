% 通过主序列和互补序列的功率延迟分布双峰定位，联合估计时频偏移

function [time_comp, freq_comp] = time_freq_recovery(pdp, prachConfig)
NIFFT = prachConfig.bwpInfo.NFFT;
L_RA = prachConfig.L_RA;
N_CS = prachConfig.N_CS;
scs_RA = prachConfig.scs_RA;
scs_PUSCH = prachConfig.bwpInfo.SCS;

overSampling = NIFFT/L_RA;
windowSize = ceil(overSampling*N_CS);

pdp_cut = zeros(windowSize,2); % 初始化截取窗口

pdp_cut(:,1) = pdp(1:windowSize,1); % 主序列窗口
pdp_cut(:,2) = pdp(NIFFT-windowSize+1:end,2); % 互补序列窗口

pdp_all = pdp_cut(:,1) + pdp_cut(:,2);





% figure; subplot(211); plot(pdp_cut(:,1),'LineWidth',1); subplot(212); plot(pdp_cut(:,2),'LineWidth',1);

[~,peak_idx1] = max(pdp_cut(:,1)); % 主序列峰值位置
[~,peak_idx2] = max(pdp_cut(:,2)); % 互补序列峰值位置

% figure; 
% hold on
% 
% plot(pdp_all)
% 
% scatter(peak_idx1,pdp_all(peak_idx1),200,'black','x','LineWidth',2.5)
% scatter(peak_idx2,pdp_all(peak_idx2),200,'black','x','LineWidth',2.5)
% hold off

delt_n1 = peak_idx1; % 主序列峰值相对窗口起始的偏移
delt_n2 = peak_idx2 - windowSize; % 互补序列峰相对窗口末尾的偏移

det_freq = (delt_n1 - delt_n2)/2;
det_time = (delt_n1 + delt_n2)/2;

% det_freq = 116.15;
% det_time = 8.85;

freq_comp = det_freq*scs_RA*1e3/overSampling; % in Hz
time_comp = det_time*scs_PUSCH/scs_RA; % in sample

end