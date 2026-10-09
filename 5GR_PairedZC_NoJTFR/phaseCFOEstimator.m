function estimate = phaseCFOEstimator(rx, reference, fs, lag)
% Weighted lag phase after removing the KNOWN COMPLETE paired reference.
% Never divide by reference samples, which can approach zero for paired ZC.
assert(size(rx,1)==numel(reference) && lag>=1 && lag<size(rx,1));
reference = reference(:);
first = 1:(size(rx,1)-lag);
second = first+lag;
weights = conj(reference(second)).*reference(first);
terms = bsxfun(@times,rx(second,:).*conj(rx(first,:)),weights);
z = sum(terms(:));
estimate.cfoHz = angle(z)*fs/(2*pi*lag);
estimate.unambiguousHz = fs/(2*lag);
estimate.coherence = abs(z)/max(sum(abs(terms(:))),realmin);
estimate.valid = abs(z)>0 && isfinite(estimate.cfoHz);
end
