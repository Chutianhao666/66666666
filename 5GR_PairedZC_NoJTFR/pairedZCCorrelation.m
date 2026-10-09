function corrResult = pairedZCCorrelation(rx, ref)
% C_u[k]=sum_n r[n] conj(x_u[(n-k) mod L]); k is zero-based.
assert(size(rx,1)==ref.length);
R = fft(rx,[],1);
plus = ifft(bsxfun(@times,R,conj(fft(ref.plus))),[],1);
minus = ifft(bsxfun(@times,R,conj(fft(ref.minus))),[],1);
% Add antenna POWERS, not complex correlations with unknown relative phases.
corrResult.plus = plus;
corrResult.minus = minus;
corrResult.power = [sum(abs(plus).^2,2),sum(abs(minus).^2,2)];
end
