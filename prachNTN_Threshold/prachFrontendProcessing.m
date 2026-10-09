function timeData = prachFrontendProcessing(rxSig, prachConfig, simConfig)

nRx = simConfig.nRx;
NCP_RA = prachConfig.NCP_RA;
scs_RA = prachConfig.scs_RA;
L_RA = prachConfig.L_RA;
Nrep = prachConfig.Nrep;
nFFT = prachConfig.nFFT;
K = prachConfig.K;
k1 = prachConfig.k1;
k_bar = prachConfig.k_bar;

rx_prach_sig = rxSig.';

if nRx >1
    rx_prach_sig = rx_prach_sig(:,NCP_RA+1:end);
else
    rx_prach_sig = rx_prach_sig(NCP_RA+1:end);
end

rx_prach_sig_fs = rx_prach_sig(:,1:end).* exp(-1i*2*pi*(K*k1+k_bar)*((1:length(rx_prach_sig))-1)./nFFT);

for nrx = 1:nRx
    % Down-sampling.
    % For L_RA = 839, decimition 24 is used for 1.25kHz and 6 is used for 5KHz.
    % For L_RA = 139, FlexRAN says down-sampling is not needed.
    if L_RA == 839
        decFilter = cascadeDecimatorGen(simConfig,prachConfig);
        switch prachConfig.mu
            case 0 % 15KHz
                if scs_RA == 1.25
                    prach_sig_fs = rx_prach_sig_fs(nrx, 1:floor(length(rx_prach_sig_fs)/12)*12);
                    prach_sig_dec = decFilter.cascadeDecimator_12(prach_sig_fs.').';
%                     prach_sig_dec = prach_sig_dec(1:nFFT/12*Nrep);
                    prach_sig_dec = prach_sig_dec(decFilter.delay+1:decFilter.delay+nFFT/12*Nrep);
                else % scs_RA = 5
                    prach_sig_fs = rx_prach_sig_fs(nrx, 1:floor(length(rx_prach_sig_fs)/3)*3);
                    prach_sig_dec = decFilter.cascadeDecimator_3(prach_sig_fs.').';
%                     prach_sig_dec = prach_sig_dec(1:nFFT/3*Nrep);
                     prach_sig_dec = prach_sig_dec(decFilter.delay+1:decFilter.delay+nFFT/3*Nrep);
                end
            case 1 % 30KHz
                if scs_RA == 1.25
                    prach_sig_fs = rx_prach_sig_fs(nrx, 1:floor(length(rx_prach_sig_fs)/24)*24);
                    prach_sig_dec = decFilter.cascadeDecimator_24(prach_sig_fs.').';
%                     prach_sig_dec = prach_sig_dec(1:nFFT/24*Nrep);
                    prach_sig_dec = prach_sig_dec(decFilter.delay+1:decFilter.delay+nFFT/24*Nrep);
                else % scs_RA = 5
                    prach_sig_fs = rx_prach_sig_fs(nrx, 1:floor(length(rx_prach_sig_fs)/6)*6);
                    prach_sig_dec = decFilter.cascadeDecimator_6(prach_sig_fs.').';
%                     prach_sig_dec = prach_sig_dec(1:nFFT/6*Nrep);
                    prach_sig_dec = prach_sig_dec(decFilter.delay+1:decFilter.delay+nFFT/6*Nrep);
                end
        end 
    else % L_RA = 139
        prach_sig_dec = rx_prach_sig_fs(nrx, 1:nFFT*Nrep);
    end
    
    % FFT
    prach_sig_dec = reshape(prach_sig_dec,[],Nrep);
    % prach_sig_freq = fft(prach_sig_dec);
    % freqData(nrx,:,:) = prach_sig_freq(1:L_RA,:);
    timeData(nrx,:,:) = prach_sig_dec;

end

% allSig = timeData(:);
% mean_pwr = allSig'*allSig/length(allSig);
% mean_amp = sqrt(mean_pwr);
% timeData = timeData/mean_amp;

end

%% Auxiliary Function
function [decFilter] = cascadeDecimatorGen(simConfig,prachConfig)
%
% Design of Cascade Filter, "FlexRAN L1 Algorithm Description" 5.5.3.
%
% cascadeDecimator_24: decimation 24, 4*3*2
% cascadeDecimator_12: decimation 12, 4*3
% cascadeDecimator_6: decimation 6, 3*2
% cascadeDecimator_3: decimation 3
%
% 1st decimator
% Fs = 122.88;                % Sampling Frequency
% Fpass = 0.54;               % Passband Frequency
% Fstop = 30.18;              % Stopband Frequency
% Dpass = 0.0057563991496;    % Passband Ripple
% Dstop = 0.001;              % Stopband Attenuation
% dens  = 20;                 % Density Factor
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
% coeff1 = firpm(N, Fo, Ao, W, {dens});
% firdecim1 = dsp.FIRDecimator(4,'Numerator',coeff1);
% 
% % 2nd decimator
% Fs = 30.72;                 % Sampling Frequency
% Fpass = 0.54;               % Passband Frequency
% Fstop = 9.7;                % Stopband Frequency
% Dpass = 0.0057563991496;    % Passband Ripple
% Dstop = 0.001;              % Stopband Attenuation
% dens  = 20;                 % Density Factor
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
% coeff2 = firpm(N, Fo, Ao, W, {dens});
% firdecim2 = dsp.FIRDecimator(3,'Numerator',coeff2);
% 
% % 3rd decimator
% Fs = 10.24;                 % Sampling Frequency
% Fpass = 0.54;               % Passband Frequency
% Fstop = 2.56;               % Stopband Frequency
% Dpass = 0.0057563991496;    % Passband Ripple
% Dstop = 0.001;              % Stopband Attenuation
% dens  = 20;                 % Density Factor
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
% coeff3 = firpm(N, Fo, Ao, W, {dens});
% firdecim3 = dsp.FIRDecimator(2,'Numerator',coeff3);
% 
% cascadeDecimator.cascadeDecimator_24 = cascade(firdecim1,firdecim2,firdecim3);
% cascadeDecimator.cascadeDecimator_6 = cascade(firdecim2,firdecim3);
% 
% decFilter.cascadeDecimator_24 = cascade(firdecim1,firdecim2,firdecim3);
% decFilter.cascadeDecimator_12 = cascade(firdecim1,firdecim2);
% decFilter.cascadeDecimator_6 = cascade(firdecim2,firdecim3);
% decFilter.cascadeDecimator_3 = cascade(firdecim2);

samplingrate = prachConfig.bwpInfo.SamplingRate/10^6;
Fpass = 0.54;               % Passband Frequency 839-->864-->864*1.25KHz/2=1.08MHz/2=0.54MHz
Dpass = 0.0057563991496;    % Passband Ripple
Dstop = 0.001;              % Stopband Attenuation
dens  = 20;                 % Density Factor

L_RA = prachConfig.L_RA; 
scs_RA = prachConfig.scs_RA;
nFFT_RA = (L_RA == 139)*1024 + (L_RA == 839)*2048 + (L_RA == 251)*1024; 
switch prachConfig.mu
    case 0  % 15KHz, sampling rate can be [15.360       30.720       61.440]MHz, 12x, 3x downsampling
%         % 1st decimator
%         Fs = 122.88;                % Sampling Frequency
%         Fpass = 0.54;               % Passband Frequency
%         Fstop = Fs/4-0.54;              % Stopband Frequency
%         Dpass = 0.0057563991496;    % Passband Ripple
%         Dstop = 0.001;              % Stopband Attenuation
%         dens  = 20;                 % Density Factor
%         [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
%         coeff1 = firpm(N, Fo, Ao, W, {dens});
%         firdecim1 = dsp.FIRDecimator(4,'Numerator',coeff1);
%         
%         % 2nd decimator
%         Fs = Fs/4;                 % Sampling Frequency
%         Fpass = 0.54;               % Passband Frequency
%         Fstop = Fs/3-0.54;                % Stopband Frequency
%         Dpass = 0.0057563991496;    % Passband Ripple
%         Dstop = 0.001;              % Stopband Attenuation
%         dens  = 20;                 % Density Factor
%         [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
%         coeff2 = firpm(N, Fo, Ao, W, {dens});
%         firdecim2 = dsp.FIRDecimator(3,'Numerator',coeff2);
%         
%         decFilter.cascadeDecimator_12 = cascade(firdecim1,firdecim2);
%         decFilter.cascadeDecimator_3 = cascade(firdecim2);
        if samplingrate == 61.44   % 4x, 3x,  FFT size =4096
            % 1st decimator
            bold =4;
            Fs = samplingrate;          % Sampling Frequency, 61.44
            targetFs = Fs/4;            % target sampling rate, 15.36MHz
            Fstop = targetFs-Fpass;     % Stopband Frequency, 14.82MHz
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff1 = firpm(N, Fo, Ao, W, {dens});
            firdecim1 = dsp.FIRDecimator(4,'Numerator',coeff1);
            
             % 2nd decimator
            bold =3;
            Fs = targetFs;                % Sampling Frequency, 15.36MHz
            targetFs = Fs/3;              % target sampling rate,5.12MHz
            Fstop = targetFs/2;           % Stopband Frequency, 2.56MHz
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff2 = firpm(N, Fo, Ao, W, {dens});
            firdecim2 = dsp.FIRDecimator(3,'Numerator',coeff2);
            
            decFilter.cascadeDecimator_12 = cascade(firdecim1,firdecim2);
            decFilter.cascadeDecimator_3 = cascade(firdecim2); 
        elseif  samplingrate ==30.72   % 4x, 3x, FFT size = 2048
            % 1st decimator
            bold =4;
            Fs = samplingrate;          % Sampling Frequency
            targetFs = Fs/4;            % target sampling rate, 7.68MHz
            Fstop = targetFs-Fpass;     % Stopband Frequency, 7.14MHz
%             Dpass = 0.0057563991496;    % Passband Ripple
%             Dstop = 0.001;              % Stopband Attenuation
%             dens  = 20;                 % Density Factor
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff1 = firpm(N, Fo, Ao, W, {dens});
            firdecim1 = dsp.FIRDecimator(4,'Numerator',coeff1);
           
             % 2nd decimator
            bold =3;
            Fs = targetFs;          % Sampling Frequency, 7.68MHz
            targetFs = Fs/3;        % target sampling rate, 2.56MHz
            Fstop = targetFs/2;     % Stopband Frequency, 1.28MHz
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff2 = firpm(N, Fo, Ao, W, {dens});
            firdecim2 = dsp.FIRDecimator(3,'Numerator',coeff2);
            
            decFilter.cascadeDecimator_12 = cascade(firdecim1,firdecim2);
            decFilter.cascadeDecimator_3 = cascade(firdecim2); 
        else % 15,36MHz
        end
        
        Fs1= samplingrate;
        Fs2= Fs1/4;
        delay_12 = length(coeff1)/2/Fs1+length(coeff2)/2/Fs2; % us
        delay_12 = floor(delay_12*nFFT_RA*scs_RA/10^3);
        
        delay_3 = length(coeff2)/2/Fs2; % us
        delay_3 = floor(delay_3*nFFT_RA*scs_RA/10^3);
        if scs_RA==1.25
            decFilter.delay= delay_12;
        else
            decFilter.delay= delay_3;
        end 
    case 1  % 30KHz
        if samplingrate == 122.88
            % 1st decimator  4x
            Fs = 122.88;                % Sampling Frequency
            targetFs = Fs/4;            % target sampling rate, 30.72MHz
            Fpass = 0.54;               % Passband Frequency
            Fstop = targetFs - Fpass;    % Stopband Frequency,  30.18MHz
%             Dpass = 0.0057563991496;    % Passband Ripple
%             Dstop = 0.001;              % Stopband Attenuation
%             dens  = 20;                 % Density Factor
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff1 = firpm(N, Fo, Ao, W, {dens});
            firdecim1 = dsp.FIRDecimator(4,'Numerator',coeff1);
            
            % 2nd decimator    3x
            Fs = 30.72;                 % Sampling Frequency
            Fpass = 0.54;               % Passband Frequency
            Fstop = 9.7;                % Stopband Frequency
%             Dpass = 0.0057563991496;    % Passband Ripple
%             Dstop = 0.001;              % Stopband Attenuation
%             dens  = 20;                 % Density Factor
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff2 = firpm(N, Fo, Ao, W, {dens});
            firdecim2 = dsp.FIRDecimator(3,'Numerator',coeff2);
            
            % 3rd decimator   2x
            Fs = 10.24;                 % Sampling Frequency
            Fpass = 0.54;               % Passband Frequency
            Fstop = 2.56;               % Stopband Frequency
%             Dpass = 0.0057563991496;    % Passband Ripple
%             Dstop = 0.001;              % Stopband Attenuation
%             dens  = 20;                 % Density Factor
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff3 = firpm(N, Fo, Ao, W, {dens});
            firdecim3 = dsp.FIRDecimator(2,'Numerator',coeff3);
            
            decFilter.cascadeDecimator_24 = cascade(firdecim1,firdecim2,firdecim3);
            decFilter.cascadeDecimator_6 = cascade(firdecim2,firdecim3);
        elseif samplingrate == 61.44
        else  % samplingrate ==30.72

        end
        Fs1= samplingrate;
        Fs2= Fs1/4;
        Fs3= Fs2/3; 
        delay_24 = length(coeff1)/2/Fs1+length(coeff2)/2/Fs2+length(coeff3)/2/Fs3; % us
        delay_24 = floor(delay_24*nFFT_RA*scs_RA/10^3);
        
             
        delay_6 = length(coeff2)/2/Fs2+length(coeff3)/2/Fs3; % us
        delay_6 = floor(delay_6*nFFT_RA*scs_RA/10^3);
        if scs_RA==1.25
            decFilter.delay= delay_24;
        else
            decFilter.delay= delay_6;
        end 
end
end

