function bwp = cal_bwp_info(simConfig)
bandWidth = simConfig.puschBW;
scs = simConfig.puschSCS;

rbTable = [38 51 65 78 92 106 119 133 162 189 217 245 273];
nrBWConfig = [15,20,25,30,35,40,45,50,60,70,80,90,100];

bwIdx = nrBWConfig == bandWidth;
nRB = rbTable(bwIdx);
nFFT = power(2,ceil(log2(nRB*12/0.85)));
samplingRate = scs*nFFT*1e3;

bwp.SCS = scs;
bwp.NRB = nRB;
bwp.NFFT = nFFT;
bwp.SamplingRate  = samplingRate;
bwp.BandWidth = bandWidth;
bwp.CyclicPrefix = 'Normal';
end