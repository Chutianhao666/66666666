function [s, ref] = pairedZCGenerate(cfg)
% Simultaneous paired roots on the same resource, matching legacy simAlgo=1.
L = cfg.length;
assert(mod(L,2)==1 && gcd(cfg.root,L)==1 && cfg.root>0 && cfg.root<L);
assert(cfg.fs==L*cfg.scsHz && cfg.maxDelaySamples<L);
n = (0:L-1).';
ref.plus = exp(-1i*pi*cfg.root*n.*(n+1)/L);
ref.minus = conj(ref.plus);
pair = ref.plus+ref.minus;
ref.scale = sqrt(mean(abs(pair).^2));
s = pair/ref.scale; % exact power normalization, not assumed sqrt(2)
ref.pair = s;
ref.length = L;
ref.root = cfg.root;
ref.fdPlus = fft(ref.plus)/sqrt(L);
ref.fdMinus = fft(ref.minus)/sqrt(L);
end
