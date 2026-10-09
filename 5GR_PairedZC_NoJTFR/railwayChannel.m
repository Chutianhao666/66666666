function [rx, truth] = railwayChannel(s, cfg, snrDb)
% Circular delay is intentional: equivalent useful-block model with ideal CP.
% Rician h is constant over this block. This is not a multipath TDL model.
[fd, rangeM] = railwayDoppler(cfg.trueMotion,cfg.carrierHz,cfg.lightSpeed);
delay = round(2*rangeM/cfg.lightSpeed*cfg.fs);
assert(delay<=cfg.maxDelaySamples,'True delay outside configured receive window.');
if strcmpi(cfg.channelType,'AWGN')
    h = ones(1,cfg.nRx);
elseif strcmpi(cfg.channelType,'Rician')
    K = 10^(cfg.ricianKDb/10);
    h = sqrt(K/(K+1))*exp(1i*2*pi*rand(1,cfg.nRx)) ...
        +sqrt(1/(2*(K+1)))*(randn(1,cfg.nRx)+1i*randn(1,cfg.nRx));
else
    error('Supported channel types: AWGN, Rician.');
end
f = fd+cfg.oscillatorHz;
n = (0:cfg.length-1).';
clean = circshift(s,delay).*exp(1i*2*pi*f*n/cfg.fs);
rx = clean*h;
% Nominal input per-complex-sample SNR, E|h|^2=1. No per-fade renormalization.
sigma2 = 10^(-snrDb/10);
rx = rx+sqrt(sigma2/2)*(randn(size(rx))+1i*randn(size(rx)));
truth.delaySamples = delay;
truth.cfoHz = f;
truth.dopplerHz = fd;
truth.rangeM = rangeM;
end
